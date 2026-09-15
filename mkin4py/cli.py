"""Local JSON interface to the original microkinetic model and solver."""
import argparse
import json
from pathlib import Path
import sys


def example():
    return {
        "environment": {"temperature": 500., "pressure": 1., "gas_constant": .00831456},
        "stoichiometry": [[-1, 1], [-1, 1], [1, -1]],
        "gas_species": [0], "pre_exponential": [2., 1.],
        "activation_energies": [0., 0.], "concentrations": [1.],
        "labels": ["A", "*", "A*"], "seed": 0,
    }


def load_model(document):
    import jax.numpy as jnp
    from .bases import Environment, MKmodel
    env = Environment()
    values = document['environment']
    env.set_temperature(values['temperature'])
    env.set_pressure(values['pressure'])
    env.set_gas_constant(values['gas_constant'])
    ms = jnp.asarray(document['stoichiometry'])
    if ms.ndim != 2:
        raise ValueError("stoichiometry must be a two-dimensional matrix")
    model = MKmodel(env)
    model.create(*ms.shape, document['gas_species'])
    model.set_ms(ms)
    model.set_kinetic_params(jnp.asarray(document['pre_exponential']).reshape(-1, 1),
                            jnp.asarray(document['activation_energies']).reshape(-1, 1))
    model.set_splabels(document.get('labels', [f'species_{i}' for i in range(ms.shape[0])]))
    model.set_concentrations(document['concentrations'])
    model.reset_model(seed=document.get('seed', 0))
    if 'coverage' in document:
        cov = jnp.asarray(document['coverage'], dtype=float)
        if cov.shape != (ms.shape[0],) or not bool(jnp.all(jnp.isfinite(cov) & (cov >= 0))):
            raise ValueError("coverage must contain one finite nonnegative value per species")
        if not bool(jnp.isclose(jnp.sum(cov[model.maps['surface']]), 1., atol=1e-8, rtol=0)):
            raise ValueError("Surface coverages must sum to one")
        if not bool(jnp.allclose(cov[model.maps['stoichs']], model.concs * env.pressure)):
            raise ValueError("Gas coverages must equal concentration times pressure")
        model.coverage = cov
    return model


def main(argv=None):
    from . import __version__
    parser = argparse.ArgumentParser(description="Solve local microkinetic catalytic systems with JAX.")
    parser.add_argument('--version', action='version', version=__version__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('methods', 'example'):
        sub.add_parser(name).add_argument('--json', action='store_true')
    solve = sub.add_parser('solve', help='Read a JSON model file (or - for stdin) and solve its steady state')
    solve.add_argument('model')
    solve.add_argument('--param', type=int, choices=(1, 2), default=1)
    solve.add_argument('--linear-solver', choices=('qmr', 'dense'), default='qmr')
    solve.add_argument('--max-iterations', type=int, default=100)
    solve.add_argument('--max-restarts', type=int, default=100)
    solve.add_argument('--max-time', type=float, default=60.)
    solve.add_argument('--float32', action='store_true')
    solve.add_argument('--json', action='store_true')
    args = parser.parse_args(argv)
    try:
        if args.command == 'example':
            result = example()
        elif args.command == 'methods':
            result = {"package": "mkin4py", "version": __version__, "execution": "local",
                      "method": "Newton correction flow with RK4",
                      "derivatives": ["analytic Jacobian", "analytic Hessian two-step correction"],
                      "linear_solvers": ["qmr", "dense"],
                      "input": "JSON environment, stoichiometry, gas_species, pre_exponential, activation_energies, concentrations",
                      "output": "coverage, species rates, elementary rates, convergence status and wall time"}
        else:
            import jax
            from .solver.solve import rk4
            jax.config.update('jax_enable_x64', not args.float32)
            source = sys.stdin.read() if args.model == '-' else Path(args.model).read_text()
            model = load_model(json.loads(source))
            result = rk4(args.param, model=model, linear_solver=args.linear_solver,
                         convtol=args.max_iterations, max_restarts=args.max_restarts, max_time=args.max_time)
            result = {k: v.tolist() if hasattr(v, 'tolist') else v for k, v in result.items()}
        print(json.dumps(result, allow_nan=False, indent=2))
        return 0 if result.get('success', True) else 1
    except (OSError, ValueError, TypeError, KeyError) as exc:
        print(json.dumps({'success': False, 'error': str(exc)}), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
