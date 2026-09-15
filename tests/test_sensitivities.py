"""Implicit sensitivities checked against independently known equilibria."""
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from mkin4py.solver.solve import steady_state


def adsorption(parameters, method='qmr', param=1, guess=.4, convtol=100):
    ka, kd, pressure = parameters
    return steady_state(jnp.array([pressure, 1-guess, guess]),
                        jnp.array([[-1, 1], [-1, 1], [1, -1]]),
                        jnp.array([ka, kd]), jnp.array([1, 2]), jnp.array([2]),
                        jnp.array(1), criteria=1e-12, convtol=convtol,
                        method=method, param=param)


def langmuir(parameters):
    ka, kd, pressure = parameters
    fraction = ka * pressure / (kd + ka * pressure)
    return jnp.array([pressure, 1-fraction, fraction])


@pytest.mark.parametrize('method', ['qmr', 'dense'])
@pytest.mark.parametrize('param', [1, 2])
def test_implicit_jacobians_and_hessian(method, param):
    parameters = jnp.array([2., 3., .7])
    solve = lambda p: adsorption(p, method, param)
    np.testing.assert_allclose(jax.jit(solve)(parameters), langmuir(parameters), atol=1e-12)
    for derivative in (jax.jacfwd, jax.jacrev):
        np.testing.assert_allclose(jax.jit(derivative(solve))(parameters),
                                   derivative(langmuir)(parameters), atol=2e-11)
    np.testing.assert_allclose(jax.jit(jax.hessian(lambda p: solve(p)[2]))(parameters),
                               jax.hessian(lambda p: langmuir(p)[2])(parameters), atol=2e-10)


def competitive(parameters):
    # A + * <=> A* and B + * <=> B*: competition for one surface balance.
    k, pressures = parameters[:4], parameters[4:]
    ms = jnp.array([[-1, 1, 0, 0], [0, 0, -1, 1],
                    [-1, 1, -1, 1], [1, -1, 0, 0], [0, 0, 1, -1]])
    return steady_state(jnp.concatenate((pressures, jnp.array([.4, .3, .3]))),
                        ms, k, jnp.array([2, 3, 4]), jnp.array([3, 4]),
                        jnp.array(2), criteria=1e-12, method='dense')


def competitive_reference(parameters):
    a = parameters[0] / parameters[1] * parameters[4]
    b = parameters[2] / parameters[3] * parameters[5]
    surface = jnp.array([1., a, b]) / (1 + a + b)
    return jnp.concatenate((parameters[4:], surface))


def test_competitive_adsorption_batch_and_arrhenius_sensitivity():
    ps = jnp.array([[2., 3., 4., 5., .7, .3], [4., 3., 2., 5., .2, .8]])
    np.testing.assert_allclose(jax.jit(jax.vmap(competitive))(ps),
                               jax.vmap(competitive_reference)(ps), atol=1e-12)
    np.testing.assert_allclose(jax.jit(jax.vmap(jax.jacrev(competitive)))(ps),
                               jax.vmap(jax.jacrev(competitive_reference))(ps), atol=2e-11)
    def at_temperature(t, solve):
        rates = jnp.array([20., 30., 40., 50.]) * jnp.exp(-jnp.array([5., 6., 7., 8.]) / (.00831456*t))
        return solve(jnp.concatenate((rates, jnp.array([.7, .3]))))[3]
    derivative = jax.jit(jax.grad(lambda t: at_temperature(t, competitive)))(500.)
    expected = jax.grad(lambda t: at_temperature(t, competitive_reference))(500.)
    np.testing.assert_allclose(derivative, expected, rtol=1e-8, atol=1e-12)
    step = .1
    finite_difference = (at_temperature(500.+step, competitive)-at_temperature(500.-step, competitive))/(2*step)
    np.testing.assert_allclose(derivative, finite_difference, rtol=2e-6, atol=1e-12)


def test_initial_surface_guess_is_not_a_physical_parameter():
    f = lambda guess: adsorption(jnp.array([2., 3., .7]), guess=guess)[2]
    np.testing.assert_allclose(jax.grad(f)(.4), 0., atol=1e-15)


def test_failed_solve_does_not_return_usable_gradients():
    f = lambda p: adsorption(p, convtol=1)[2]
    p = jnp.array([2., 3., .7])
    assert jnp.isnan(f(p))
    assert jnp.all(jnp.isnan(jax.grad(f)(p)))


def test_single_surface_species():
    f = lambda pressure: steady_state(jnp.array([pressure, 1.]),
        jnp.array([[0], [0]]), jnp.array([2.]), jnp.array([1]),
        jnp.array([], dtype=int), jnp.array(1))
    np.testing.assert_allclose(jax.jit(f)(3.), [3., 1.])
    np.testing.assert_allclose(jax.jacrev(f)(3.), [1., 0.])


def test_catalytic_conversion_with_nonzero_throughput():
    # A + * <=> A*, A* <=> B*, B* <=> B + *; gas is held out of equilibrium.
    ms = jnp.array([[-1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, -1],
                    [-1, 1, 0, 0, 1, -1], [1, -1, -1, 1, 0, 0],
                    [0, 0, 1, -1, -1, 1]])
    def flux(k):
        cov = steady_state(jnp.array([.8, .05, .4, .3, .3]), ms, k,
            jnp.array([2, 3, 4]), jnp.array([3, 4]), jnp.array(2), criteria=1e-12)
        return k[4]*cov[4] - k[5]*.05*cov[2]
    def reference(k):
        # Stationary linear balances and normalization, independent of RK4.
        matrix = np.array([[k[0]*.8, -(k[1]+k[2]), k[3]],
                           [k[5]*.05, k[2], -(k[3]+k[4])], [1., 1., 1.]])
        vac, _, b = np.linalg.solve(matrix, [0., 0., 1.])
        return k[4]*b - k[5]*.05*vac
    k = np.array([3., 2., 4., 1., 5., .5])
    assert float(flux(k)) > 0
    np.testing.assert_allclose(flux(k), reference(k), atol=1e-12)
    step = 1e-4
    finite = np.array([(reference(k+step*direction)-reference(k-step*direction))/(2*step)
                       for direction in np.eye(k.size)])
    np.testing.assert_allclose(jax.jit(jax.grad(flux))(jnp.asarray(k)), finite,
                               rtol=2e-8, atol=1e-10)


def test_singular_equilibrium_has_invalid_sensitivity():
    # At zero adsorption and desorption every surface composition is a root.
    f = lambda k: adsorption(jnp.array([k, 0., 1.]))[2]
    assert jnp.isfinite(f(0.))
    assert not jnp.isfinite(jax.grad(f)(0.))


@pytest.mark.parametrize('scale', [1e-20, 1e-12, 1., 1e12, 1e20])
def test_rate_timescale_does_not_change_equilibrium_or_sensitivity(scale):
    f = lambda pressure: adsorption(jnp.array([2.*scale, 3.*scale, pressure]))[2]
    p = .7
    expected = 2*p/(3+2*p)
    derivative = 6/(3+2*p)**2
    np.testing.assert_allclose(jax.jit(f)(p), expected, atol=2e-12)
    np.testing.assert_allclose(jax.jit(jax.grad(f))(p), derivative, atol=2e-11)
    step = 1e-4
    np.testing.assert_allclose((f(p+step)-f(p-step))/(2*step), derivative, rtol=2e-8)


@pytest.mark.parametrize('option,value', [('criteria', jnp.inf), ('criteria', 0.),
    ('inner_criteria', -1.), ('h', jnp.nan), ('h', 1.1), ('hfun', 0.),
    ('delta_min', 2.), ('convtol', 0), ('convtolH', -1), ('inner_convtol', 0)])
def test_invalid_dynamic_controls_are_not_reported_as_roots(option, value):
    def f(k):
        return steady_state(jnp.array([.7, .6, .4]),
            jnp.array([[-1, 1], [-1, 1], [1, -1]]), k,
            jnp.array([1, 2]), jnp.array([2]), jnp.array(1), **{option: value})[2]
    k = jnp.array([2., 3.])
    assert jnp.isnan(jax.jit(f)(k))
    assert jnp.all(jnp.isnan(jax.jit(jax.grad(f))(k)))
