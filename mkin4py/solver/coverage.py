"""Functional surface-site updates used by the original Newton/RK4 method."""
import jax
import jax.numpy as np
from .params import convergence_params
from ..bases import mkmodel


@jax.jit
def update(coverage, delta, h, surface, xsurface, ndof, delta_min):
    cov = coverage.at[xsurface].add(h * delta)
    cov = cov.at[ndof].add(-h * np.sum(delta))
    values = np.maximum(cov[surface], delta_min)
    return cov.at[surface].set(values / np.sum(values))


def coverage_update(delta, h=1, cov=None):
    """Update the configured model; ``cov`` optionally supplies the starting point."""
    if cov is None or (isinstance(cov, (tuple, list)) and len(cov) == 0):
        cov = mkmodel.coverage
    maps = mkmodel.maps
    mkmodel.coverage = update(np.asarray(cov), np.asarray(delta), h,
                             maps['surface'], maps['xsurface'], maps['ndof'],
                             convergence_params['delta_min'])
    return mkmodel.coverage
