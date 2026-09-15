"""Explicit mass-action derivatives with the original site elimination."""
from functools import partial
import jax
import jax.numpy as np
from .coverage import update
from .params import convergence_params
from ..bases import mkmodel


def _product(cov, powers):
    # Substituting 1 for unused bases avoids 0**0 derivative singularities.
    return np.prod(np.where(powers > 0, cov, 1.) ** powers, axis=-1)


@jax.jit
def power_law(cov, ms):
    """Concentration parcel psi, shape (reactions, 1); accepts jit/grad/vmap."""
    return _product(cov, np.maximum(-ms.T, 0))[:, None]


@partial(jax.jit, static_argnames=('param',))
def evaluate(cov, ms, xsurface, ndof, param=1):
    """Pure analytic psi/Jacobian/Hessian; cov must satisfy the site balance.

    Jacobian shape: (reactions, independent surface species).
    Hessian shape: (independent species, independent species, reactions).
    """
    powers = np.maximum(-ms.T, 0)
    psi = _product(cov, powers)[:, None]
    if param == 0:
        return psi
    indices = np.concatenate((xsurface, np.reshape(ndof, (1,))))
    eye = np.eye(cov.size, dtype=ms.dtype)[indices]
    reduced = np.maximum(powers[:, None, :] - eye[None, :, :], 0)
    first = powers[:, indices] * _product(cov, reduced)
    jacobian = first[:, :-1] - first[:, -1:]
    if param == 1:
        return psi, jacobian, []
    exponents = np.maximum(reduced[:, :, None, :] - eye[None, None, :, :], 0)
    coefficients = powers[:, indices, None] * (
        powers[:, None, indices] - np.eye(indices.size, dtype=ms.dtype)[None, :, :])
    second = coefficients * _product(cov, exponents)
    hess = second[:, :-1, :-1] - second[:, :-1, -1:] - second[:, -1:, :-1] + second[:, -1:, -1:]
    return psi, jacobian, np.moveaxis(hess, 0, -1)


def analytical(param=1, cov=None):
    """Original configured-model API. Use ``evaluate`` inside JAX transforms."""
    if param not in (0, 1, 2):
        raise ValueError("param must be 0, 1 or 2")
    if cov is None or (isinstance(cov, (tuple, list)) and len(cov) == 0):
        cov = mkmodel.coverage
    maps = mkmodel.maps
    mkmodel.coverage = update(np.asarray(cov), np.zeros(len(maps['xsurface'])), 0.,
                             maps['surface'], maps['xsurface'], maps['ndof'],
                             convergence_params['delta_min'])
    return evaluate(mkmodel.coverage, mkmodel.ms, maps['xsurface'], maps['ndof'], param)
