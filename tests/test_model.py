import jax
import jax.numpy as jnp
import numpy as np
import pytest
from mkin4py.cli import example, load_model
from mkin4py.solver.solve import rk4, attempt
from mkin4py.solver.derivatives import power_law, evaluate
from mkin4py.solver.linsolver import qmr


@pytest.mark.parametrize('method', ['qmr', 'dense'])
@pytest.mark.parametrize('param', [1, 2])
def test_langmuir_equilibrium(method, param):
    model = load_model(example())
    result = rk4(param, model=model, linear_solver=method, max_restarts=0)
    assert result['success'], result
    np.testing.assert_allclose(result['coverage'].ravel(), [1., 1/3, 2/3], atol=1e-8)
    np.testing.assert_allclose(result['rates'], model.ms @ result['elem_rate'], atol=1e-14)
    assert result['time'] > 0


def test_recompile_does_not_capture_old_parameters():
    model = load_model(example())
    a = rk4(model=model)
    model.set_kinetic_params([[8.], [1.]], [[0.], [0.]])
    b = rk4(model=model)
    np.testing.assert_allclose(b['coverage'].ravel(), [1., 1/9, 8/9], atol=1e-8)
    assert not np.allclose(a['coverage'], b['coverage'])
    model.environment.set_pressure(2.)
    c = rk4(model=model)
    np.testing.assert_allclose(c['coverage'].ravel(), [2., 1/17, 16/17], atol=1e-8)


def test_failure_returns_consistent_rates_and_coverage():
    model = load_model(example())
    result = rk4(model=model, convtol=1, max_restarts=1)
    assert not result['success']
    assert result['iterations'] == 2
    np.testing.assert_allclose(result['elem_rate'], model.kinetic_parameters['k'][:, None] *
                               power_law(result['coverage'].ravel(), model.ms))
    np.testing.assert_allclose(result['rates'], model.ms @ result['elem_rate'])


@pytest.mark.parametrize('bad', [(-1, 2, [0]), (3, 2, [0,0]), (3,2,[3]), (3,2,[.5]), (3,2,[0,1,2])])
def test_invalid_dimensions(bad):
    model = load_model(example())
    with pytest.raises((TypeError, ValueError)):
        model.create(*bad)


@pytest.mark.parametrize('coverage', [[0., .4, .6, 0.], [.1, .2, .3, .5]])
def test_analytic_derivatives_match_autodiff(coverage):
    # Repeated reactants, zero gas/surface concentration, and a cubic reaction.
    ms = jnp.array([[-1,0,0],[-2,-1,0],[1,-1,-3],[1,2,3]])
    cov = jnp.array(coverage)
    independent = jnp.array([1,2])
    dof = jnp.array(3)
    def parcel(x):
        c = cov.at[independent].set(x).at[dof].set(1-jnp.sum(x))
        return power_law(c, ms).ravel()
    cov = cov.at[dof].set(1-jnp.sum(cov[independent]))
    psi, jac, hess = evaluate(cov, ms, independent, dof, 2)
    np.testing.assert_allclose(psi.ravel(), parcel(cov[independent]), atol=1e-15)
    np.testing.assert_allclose(jac, jax.jacfwd(parcel)(cov[independent]), atol=1e-14)
    np.testing.assert_allclose(np.moveaxis(hess, -1, 0), jax.jacfwd(jax.jacrev(parcel))(cov[independent]), atol=1e-14)


def test_jit_vmap_rates():
    ms = jnp.array([[-1,1],[-1,1],[1,-1]])
    xs = jnp.array([[1., .3,.7], [2., .4,.6]])
    actual = jax.jit(jax.vmap(power_law, in_axes=(0,None)))(xs, ms)
    np.testing.assert_allclose(actual[:,:,0], [[.3,.7],[.8,.6]])


@pytest.mark.parametrize('size', [1,3,10])
def test_qmr_against_scipy(size):
    from scipy.sparse.linalg import qmr as reference
    rng = np.random.default_rng(size)
    a = rng.normal(size=(size,size)) + size*np.eye(size)
    b = rng.normal(size=size)
    expected, info = reference(a,b,rtol=1e-10,atol=1e-12)
    assert info == 0
    actual, info = qmr(jnp.asarray(a), jnp.asarray(b), rtol=1e-10, atol=1e-12)
    assert info == 0
    np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-10)


def test_qmr_zero_rhs_and_singular():
    x, info = qmr(jnp.eye(2), jnp.zeros(2))
    assert info == 0
    np.testing.assert_array_equal(x, [0,0])
    _, info = qmr(jnp.zeros((2,2)), jnp.ones(2))
    assert info != 0


def test_original_module_api():
    import mkin4py as m
    model = load_model(example())
    original = m.mkmodel
    # Configure the original singleton via its setters, as in README.
    m.environment.set_temperature(500)
    m.environment.set_pressure(1)
    m.environment.set_gas_constant(.00831456)
    original.create(3,2,[0])
    original.set_ms(model.ms)
    original.set_splabels(model.splabels)
    original.set_kinetic_params([[2.],[1.]], [[0.],[0.]])
    original.set_concentrations([1.])
    assert isinstance(m.solver.derivatives.analytical(0), jax.Array)
    assert m.solver.derivatives.analytical(2)[2].shape == (1,1,2)
    assert m.solver.linsolver.newton_type(1).shape == (1,)
    assert m.solver.coverage.coverage_update(jnp.zeros(1)).shape == (3,)
    assert m.solver.solve.rk4()['success']
