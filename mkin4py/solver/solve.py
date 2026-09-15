"""Fourth-order Runge–Kutta integration of the Newton correction flow.

This is the original steady-state algorithm, not physical time integration.
"""
from functools import partial
from time import perf_counter
import math
import jax
import jax.numpy as np
from .derivatives import power_law
from .coverage import update
from .linsolver import correction
from .params import convergence_params
from ..bases import mkmodel


@partial(jax.jit, static_argnames=('param', 'method'))
def attempt(cov, ms, msa, surface, xsurface, ndof, h=1., hfun=.995, delta_min=1e-30,
            criteria=1e-8, inner_criteria=1e-9, convtol=100, convtolH=20,
            inner_convtol=300, param=1, method='qmr'):
    """One compiled attempt. All model values are explicit dynamic arguments."""
    msas = msa[xsurface]
    def residual(c):
        return np.max(np.abs((msa @ power_law(c, ms))[surface]))
    def step(c):
        return correction(c, ms, msas, xsurface, ndof, hfun, inner_criteria,
                          inner_convtol, convtolH, param, method)
    def move(c, delta, fraction):
        return update(c, delta, fraction, surface, xsurface, ndof, delta_min)
    def cond(state):
        count, c, res = state
        return (count < convtol) & (res > criteria) & np.isfinite(res) & np.all(np.isfinite(c))
    def body(state):
        count, c, _ = state
        k1 = step(c)
        k2 = step(move(c, k1, .5))
        k3 = step(move(c, k2, .5))
        k4 = step(move(c, k3, 1.))
        c = move(c, (k1+k4)/6. + (k2+k3)/3., h)
        return count+1, c, residual(c)
    cov = move(cov, np.zeros(xsurface.size, dtype=cov.dtype), 0.)
    return jax.lax.while_loop(cond, body, (0, cov, residual(cov)))


def rk4(param=1, *, model=None, **options):
    """Solve the configured model, preserving original result keys.

    ``time`` includes dispatch, compilation and synchronization. ``max_time``
    is checked between compiled attempts; a running kernel cannot be preempted.
    For high-stiffness problems enable JAX 64-bit mode before creating arrays.
    """
    t0 = perf_counter()
    if param not in (1, 2):
        raise ValueError("param must be 1 or 2")
    unknown = options.keys() - convergence_params.keys()
    if unknown:
        raise ValueError(f"Unknown solver options: {sorted(unknown)}")
    p = dict(convergence_params, **options)
    if p['linear_solver'] not in ('qmr', 'dense'):
        raise ValueError("linear_solver must be qmr or dense")
    for key in ('convtol', 'convtolH', 'inner_convtol', 'max_restarts'):
        if isinstance(p[key], bool) or not isinstance(p[key], int) or p[key] < (0 if key == 'max_restarts' else 1):
            raise ValueError(f"{key} must be an integer within its supported range")
    for key in ('h', 'hfun', 'delta_min', 'criteria', 'inner_criteria', 'max_time'):
        if not math.isfinite(p[key]) or p[key] <= 0:
            raise ValueError(f"{key} must be finite and positive")
    if p['h'] > 1 or p['hfun'] > 1:
        raise ValueError("h and hfun must not exceed 1")
    model = mkmodel if model is None else model
    model.update_model()
    total, restarts, status = 0, 0, 'iteration_limit'
    kernel_options = {k: p[k] for k in ('h', 'hfun', 'delta_min', 'criteria',
        'inner_criteria', 'convtol', 'convtolH', 'inner_convtol')}
    for restart in range(p['max_restarts']+1):
        restarts = restart
        m = model.maps
        count, cov, residual = attempt(model.coverage, model.ms, m['msa'], m['surface'],
            m['xsurface'], m['ndof'], **kernel_options, param=param, method=p['linear_solver'])
        jax.block_until_ready((count, cov, residual))
        total += int(count)
        model.coverage = cov
        if bool(np.isfinite(residual)) and float(residual) <= p['criteria']:
            status = 'converged'
            break
        status = 'numerical_failure' if not bool(np.isfinite(residual)) else 'iteration_limit'
        if perf_counter()-t0 >= p['max_time']:
            status = 'time_limit'
            break
        if restart < p['max_restarts']:
            model.init_coverage()
    psi = power_law(model.coverage, model.ms)
    elem = model.kinetic_parameters['k'][:, None] * psi
    rates = model.maps['msa'] @ psi
    jax.block_until_ready((elem, rates))
    success = status == 'converged'
    return {'coverage': model.coverage[:, None], 'rates': rates, 'elem_rate': elem,
            'msg': 'Convergence achieved' if success else 'Convergence NOT achieved.',
            'time': perf_counter()-t0, 'success': success, 'status': status,
            'iterations': total, 'restarts': restarts, 'residual': float(residual)}
