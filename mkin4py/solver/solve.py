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


@partial(jax.jit, static_argnames=('param', 'method'))
def steady_state(cov, ms, k, surface, xsurface, ndof, *, h=1., hfun=.995,
                 delta_min=1e-30, criteria=1e-8, inner_criteria=1e-9,
                 convtol=100, convtolH=20, inner_convtol=300,
                 param=1, method='qmr'):
    """Return steady-state coverage with implicit forward/reverse derivatives.

    Uses the original Newton/RK4 ``attempt`` to find the root. All arrays are
    explicit: ``cov`` is the initial species vector (gas entries are fixed
    partial pressures), ``ms`` is species by reactions, ``k`` is the reaction
    constant vector, and the index arrays are the original model maps.

    Gradients describe a locally isolated steady state, not the iterations or
    initial surface guess. Differentiate rate constants or fixed gas entries;
    hold the stoichiometry and index maps fixed. The reduced surface Jacobian
    must be nonsingular. Failed convergence returns NaN values and derivatives.
    No Python model mutation, random restart or wall-time cutoff occurs here.
    """
    if param not in (1, 2) or method not in ('qmr', 'dense'):
        raise ValueError("param must be 1 or 2 and method must be qmr or dense")
    cov, ms, k = np.asarray(cov, dtype=float), np.asarray(ms), np.asarray(k, dtype=float)
    surface, xsurface, ndof = np.asarray(surface), np.asarray(xsurface), np.asarray(ndof)
    if (cov.ndim != 1 or ms.ndim != 2 or ms.shape[0] != cov.size
            or k.shape != (ms.shape[1],) or surface.ndim != 1
            or surface.size != xsurface.size + 1 or xsurface.ndim != 1
            or ndof.ndim != 0):
        raise ValueError("Provide a species vector, reaction vector and matching model index maps")
    if not all(np.issubdtype(x.dtype, np.integer) for x in (ms, surface, xsurface, ndof)):
        raise TypeError("Stoichiometry and model index maps must have integer dtype")
    msa = ms * k[None, :]

    def expand(x):
        return cov.at[xsurface].set(x).at[ndof].set(1. - np.sum(x))

    def residual(x):
        return (msa[xsurface] @ power_law(expand(x), ms)).ravel()

    def solve_root(_, initial):
        _, solution, _ = attempt(expand(initial), ms, msa, surface, xsurface, ndof,
            h=h, hfun=hfun, delta_min=delta_min, criteria=criteria,
            inner_criteria=inner_criteria, convtol=convtol, convtolH=convtolH,
            inner_convtol=inner_convtol, param=param, method=method)
        return solution[xsurface]

    def tangent_solve(linear, rhs):
        # Small dense surface systems: materialize the linearized residual,
        # equilibrate it, and let JAX transpose the solve for reverse mode.
        matrix = jax.jacfwd(linear)(np.zeros_like(rhs))
        scale = np.maximum(np.max(np.abs(matrix), axis=1), np.finfo(cov.dtype).tiny)
        return np.linalg.solve(matrix / scale[:, None], rhs / scale)

    if xsurface.size:
        root = jax.lax.custom_root(residual, cov[xsurface], solve_root, tangent_solve)
        solution = expand(root)
    else:
        solution = expand(np.empty((0,), dtype=cov.dtype))
    error = np.max(np.abs((msa @ power_law(solution, ms))[surface]))
    ordered = np.sort(surface)
    indices_valid = (np.all((surface >= 0) & (surface < cov.size))
                     & np.all(ordered[1:] > ordered[:-1])
                     & np.all(ordered == np.sort(np.concatenate((xsurface, ndof[None])))))
    valid = (indices_valid & np.all(np.isfinite(solution)) & np.all(np.isfinite(k))
             & np.all(k >= 0) & np.all(solution >= 0) & (error <= criteria))
    # Multiplication by NaN also marks gradients of an unconverged solve as
    # invalid; a constant NaN branch would misleadingly differentiate to zero.
    return solution * jax.lax.stop_gradient(np.where(valid, 1., np.nan))


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
