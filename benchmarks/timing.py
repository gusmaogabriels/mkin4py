"""Synchronized compile and warm timings on identical initial coverages."""
import argparse
import json
from pathlib import Path
import platform
import statistics
import textwrap
from time import perf_counter
import jax
from mkin4py.cli import example, load_model
from mkin4py.solver.solve import attempt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repeats', type=int, default=20)
    parser.add_argument('--output', default='build/timing.json')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('repeats must be positive')
    jax.config.update('jax_enable_x64', True)
    import mkin4py
    simple = load_model(example())
    readme = (Path(__file__).parents[1]/'README.md').read_text()
    code = textwrap.dedent(readme[readme.index('        import mkin4py'):readme.index('  - **Evaluation**:')])
    exec(compile(code,'README.md','exec'), {})
    mkin4py.mkmodel.reset_model(seed=0)
    report = {'jax': jax.__version__, 'python': platform.python_version(),
              'platform': platform.platform(), 'device': str(jax.devices()[0]),
              'dtype': 'float64', 'repeats': args.repeats,
              'note': 'Each warm call starts from the same unsolved coverage. Timings synchronize device work. Import/startup excluded.',
              'cases': []}
    for name, model in [('adsorption', simple), ('ethylene_epoxidation', mkin4py.mkmodel)]:
        m = model.maps
        inputs = (model.coverage, model.ms, m['msa'], m['surface'], m['xsurface'], m['ndof'])
        jax.block_until_ready(inputs)
        for method in ('qmr','dense'):
            start = perf_counter()
            lowered = attempt.lower(*inputs, method=method)
            lowering = perf_counter()-start
            start = perf_counter()
            kernel = lowered.compile()
            compilation = perf_counter()-start
            start = perf_counter()
            result = jax.block_until_ready(kernel(*inputs))
            first = perf_counter()-start
            samples = []
            for _ in range(args.repeats):
                start = perf_counter()
                result = jax.block_until_ready(kernel(*inputs))
                samples.append(perf_counter()-start)
            assert float(result[2]) <= 1e-8, result
            report['cases'].append({'case': name, 'linear_solver': method,
                'lowering_seconds': lowering, 'compilation_seconds': compilation,
                'first_execution_seconds': first, 'cold_total_seconds': lowering+compilation+first,
                'warm_median_seconds': statistics.median(samples),
                'warm_min_seconds': min(samples), 'warm_max_seconds': max(samples),
                'warm_samples_seconds': samples, 'iterations': int(result[0]), 'residual': float(result[2])})
    path = Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
