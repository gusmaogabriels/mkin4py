"""JAX QMR and the original Newton / two-step Hessian corrections.

The QMR recurrence follows SciPy's unpreconditioned real implementation,
expressed as a JAX loop. See NOTICE.md for its BSD license. No SciPy solver
is called at runtime.
"""
from functools import partial
import jax
import jax.numpy as np
from .derivatives import evaluate
from .params import convergence_params
from ..bases import mkmodel


@jax.jit
def qmr(a, b, rtol=1e-9, atol=0., maxiter=300):
    """Real dense QMR; return (solution, info), 0 success, 1 limit, -1 breakdown."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float).reshape(-1)
    zeros = np.zeros_like(b)
    norm = np.linalg.norm(b)
    tolerance = np.maximum(atol, rtol * norm)
    state = dict(i=0, x=zeros, r=b, v=b, w=b, rho=norm, xi=norm,
                 gamma=np.ones_like(norm), eta=-np.ones_like(norm), theta=np.zeros_like(norm),
                 epsilon=np.ones_like(norm), p=zeros, q=zeros, d=zeros, s=zeros, failed=False)
    tiny = np.finfo(b.dtype).tiny
    safe = lambda x: np.where(np.abs(x) > tiny, x, np.ones_like(x))

    def cond(s):
        return (s['i'] < maxiter) & ~s['failed'] & (np.linalg.norm(s['r']) > tolerance)

    def body(s):
        v, w = s['v'] / safe(s['rho']), s['w'] / safe(s['xi'])
        delta = np.dot(w, v)
        p = v - np.where(s['i'] > 0, s['xi'] * delta / safe(s['epsilon']), 0.) * s['p']
        q = w - np.where(s['i'] > 0, s['rho'] * delta / safe(s['epsilon']), 0.) * s['q']
        ap = a @ p
        epsilon = np.dot(q, ap)
        beta = epsilon / safe(delta)
        vn, wn = ap - beta * v, a.T @ q - beta * w
        rho, xi = np.linalg.norm(vn), np.linalg.norm(wn)
        theta = rho / safe(s['gamma'] * np.abs(beta))
        gamma = 1 / np.hypot(1., theta)
        eta = -s['eta'] * (s['rho'] / safe(beta)) * (gamma / safe(s['gamma']))**2
        factor = (s['theta'] * gamma)**2
        d, ds = eta * p + factor * s['d'], eta * ap + factor * s['s']
        x, r = s['x'] + d, s['r'] - ds
        failed = ((np.abs(s['rho']) <= tiny) | (np.abs(s['xi']) <= tiny) |
                  (np.abs(delta) <= tiny) | (np.abs(epsilon) <= tiny) |
                  (np.abs(beta) <= tiny) | ~np.all(np.isfinite(x)))
        return dict(i=s['i']+1, x=np.where(failed, s['x'], x), r=r, v=vn, w=wn,
                    rho=rho, xi=xi, gamma=gamma, eta=eta, theta=theta,
                    epsilon=epsilon, p=p, q=q, d=d, s=ds, failed=failed)

    result = jax.lax.while_loop(cond, body, state)
    # Check the true residual, including after the last permitted iteration.
    success = np.linalg.norm(b - a @ result['x']) <= tolerance
    info = np.where(success, 0, np.where(result['failed'], -1, 1))
    return result['x'], info


def _linear(a, b, tolerance, maxiter, method):
    # Row equilibration preserves the Newton equation and reduces stiffness.
    scale = np.maximum(np.max(np.abs(a), axis=1), np.finfo(a.dtype).tiny)
    matrix, rhs = a / scale[:, None], b / scale
    if method == 'dense':
        return np.linalg.solve(matrix, rhs)
    result, info = qmr(matrix, rhs, tolerance, 0., maxiter)
    # QMR can break down even for nonsingular matrices. A direct factorization
    # is a local JAX fallback; nonfinite singular solutions propagate to status.
    return jax.lax.cond(info == 0, lambda: result, lambda: np.linalg.solve(matrix, rhs))


@partial(jax.jit, static_argnames=('param', 'method'))
def correction(cov, ms, msas, xsurface, ndof, hfun=.995, tolerance=1e-9,
               maxiter=300, hessian_steps=20, param=1, method='qmr'):
    psi, jac, hess = evaluate(cov, ms, xsurface, ndof, param)
    f, matrix = (msas @ psi).ravel(), msas @ jac
    if xsurface.size == 0:
        return np.zeros((0,), dtype=cov.dtype)
    dx = _linear(matrix, -hfun * f, tolerance, maxiter, method)
    if param == 2:
        hessian = np.einsum('rm,ijm->rij', msas, hess)
        def cond(s):
            i, delta, change = s
            return (i < hessian_steps) & (np.max(np.abs(change)) > tolerance) & np.all(np.isfinite(delta))
        def body(s):
            i, delta, _ = s
            hd = np.einsum('rij,j->ri', hessian, delta)
            change = _linear(matrix + hd, -(f + matrix @ delta + .5 * hd @ delta),
                             tolerance, maxiter, method)
            return i+1, delta+change, change
        dx = jax.lax.cond(np.max(np.abs(dx)) < 1.,
            lambda: jax.lax.while_loop(cond, body, (0, dx, dx))[1], lambda: dx)
    return dx


def newton_type(param=1):
    if param not in (1, 2):
        raise ValueError("param must be 1 or 2")
    p, m = convergence_params, mkmodel.maps
    return correction(mkmodel.coverage, mkmodel.ms, m['msas'], m['xsurface'], m['ndof'],
                      p['hfun'], p['inner_criteria'], p['inner_convtol'], p['convtolH'],
                      param, p['linear_solver'])
