"""Measure compiled equilibrium values and gradients with synchronized calls."""
import argparse
import json
from pathlib import Path
import platform
import statistics
from time import perf_counter

import jax
import jax.numpy as jnp
from mkin4py.solver.solve import steady_state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeats', type=int, default=20)
    parser.add_argument('--output', default='build/sensitivities.json')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('repeats must be positive')
    jax.config.update('jax_enable_x64', True)
    ms = jnp.array([[-1, 1, 0, 0], [0, 0, -1, 1],
                    [-1, 1, -1, 1], [1, -1, 0, 0], [0, 0, 1, -1]])
    cov = jnp.array([.7, .3, .4, .3, .3])
    k = jnp.array([2., 3., 4., 5.])
    surface, independent, dependent = jnp.array([2, 3, 4]), jnp.array([3, 4]), jnp.array(2)
    report = {'jax': jax.__version__, 'python': platform.python_version(),
              'platform': platform.platform(), 'device': str(jax.devices()[0]),
              'dtype': 'float64', 'case': 'competitive_adsorption',
              'note': 'Identical device-resident inputs; each call synchronized. Lowering/compilation reported separately. Excludes interpreter, imports and model setup.',
              'cases': []}
    jax.block_until_ready((cov, ms, k, surface, independent, dependent))
    for method in ('qmr', 'dense'):
        def coverage(constants):
            return steady_state(cov, ms, constants, surface, independent, dependent,
                                criteria=1e-12, method=method)[3]
        for name, function in [('value', coverage), ('value_and_grad', jax.value_and_grad(coverage))]:
            start = perf_counter()
            lowered = jax.jit(function).lower(k)
            lowering = perf_counter()-start
            start = perf_counter()
            compiled = lowered.compile()
            compilation = perf_counter()-start
            start = perf_counter()
            result = jax.block_until_ready(compiled(k))
            first = perf_counter()-start
            assert all(bool(jnp.all(jnp.isfinite(x))) for x in jax.tree_util.tree_leaves(result))
            samples = []
            for _ in range(args.repeats):
                start = perf_counter()
                jax.block_until_ready(compiled(k))
                samples.append(perf_counter()-start)
            report['cases'].append({'method': method, 'operation': name,
                'lowering_seconds': lowering, 'compilation_seconds': compilation,
                'first_execution_seconds': first, 'cold_total_seconds': lowering+compilation+first,
                'warm_median_seconds': statistics.median(samples),
                'warm_min_seconds': min(samples), 'warm_max_seconds': max(samples),
                'warm_samples_seconds': samples})
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
