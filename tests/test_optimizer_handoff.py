"""Fit an adsorption constant by composing mkin4py into released Optinpy."""
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from mkin4py.solver.solve import steady_state


def test_fit_adsorption_constant_with_optinpy():
    optinpy = pytest.importorskip('optinpy')
    pressures = jnp.array([.2, .5, 1., 2.])
    observed = 4*pressures / (1+4*pressures)
    ms = jnp.array([[-1, 1], [-1, 1], [1, -1]])
    def predicted(log_k, pressure):
        cov = jnp.array([pressure, .5, .5])
        return steady_state(cov, ms, jnp.array([jnp.exp(log_k), 1.]),
            jnp.array([1, 2]), jnp.array([2]), jnp.array(1),
            method='dense', criteria=1e-12)[2]
    def loss(x):
        residuals = jax.vmap(predicted, in_axes=(None, 0))(x[0], pressures) - observed
        return jnp.sum(residuals**2)
    solve = optinpy.compile_minimizer(loss, method='bfgs', tol=1e-9, max_iter=100)
    result = jax.block_until_ready(solve(jnp.array([jnp.log(2.)])))
    assert bool(result['success']), result
    np.testing.assert_allclose(jnp.exp(result['x'][0]), 4., rtol=1e-6)
    assert float(loss(result['x'])) < 1e-14
