"""Synchronized compile and warm timings on identical initial coverages."""
import argparse
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import textwrap
from time import perf_counter
import jax
from mkin4py.cli import example, load_model
from mkin4py.solver.solve import attempt, rk4
from mkin4py.solver.params import default


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
              'note': 'Identical initial coverages; synchronized device work. Kernel timings exclude Python setup. API timings include model refresh. CLI timing includes a fresh process, imports, setup, compile and solve.',
              'cases': []}
    controls = {k: default[k] for k in ('h', 'hfun', 'delta_min', 'criteria',
        'inner_criteria', 'convtol', 'convtolH', 'inner_convtol')}
    for name, model in [('adsorption', simple), ('ethylene_epoxidation', mkin4py.mkmodel)]:
        m = model.maps
        inputs = (model.coverage, model.ms, m['msa'], m['surface'], m['xsurface'], m['ndof'])
        jax.block_until_ready(inputs)
        for method in ('qmr','dense'):
            model.coverage = inputs[0]
            start = perf_counter()
            lowered = attempt.lower(*inputs, **controls, param=1, method=method)
            lowering = perf_counter()-start
            start = perf_counter()
            kernel = lowered.compile()
            compilation = perf_counter()-start
            start = perf_counter()
            result = jax.block_until_ready(kernel(*inputs, **controls))
            first = perf_counter()-start
            samples = []
            for _ in range(args.repeats):
                start = perf_counter()
                result = jax.block_until_ready(kernel(*inputs, **controls))
                samples.append(perf_counter()-start)
            assert float(result[2]) <= 1e-8, result
            api_samples = []
            for _ in range(args.repeats+1):
                model.coverage = inputs[0]
                start = perf_counter()
                api_result = rk4(model=model, linear_solver=method, max_restarts=0)
                api_samples.append(perf_counter()-start)
                assert api_result['success'], api_result
            document = {
                'environment': {k: getattr(model.environment, k) for k in
                    ('temperature', 'pressure', 'gas_constant')},
                'stoichiometry': model.ms.tolist(), 'gas_species': m['stoichs'].tolist(),
                'pre_exponential': model.kinetic_parameters['va'].ravel().tolist(),
                'activation_energies': model.kinetic_parameters['vea'].ravel().tolist(),
                'concentrations': model.concs.tolist(), 'seed': 0,
                'coverage': inputs[0].tolist(),
            }
            start = perf_counter()
            process = subprocess.run([sys.executable, '-m', 'mkin4py', 'solve', '-',
                '--linear-solver', method, '--max-restarts', '0'],
                input=json.dumps(document), capture_output=True, text=True, check=True)
            cli_time = perf_counter()-start
            cli_result = json.loads(process.stdout)
            assert cli_result['success'], cli_result
            report['cases'].append({'case': name, 'linear_solver': method,
                'lowering_seconds': lowering, 'compilation_seconds': compilation,
                'first_execution_seconds': first, 'cold_total_seconds': lowering+compilation+first,
                'warm_median_seconds': statistics.median(samples),
                'warm_min_seconds': min(samples), 'warm_max_seconds': max(samples),
                'warm_samples_seconds': samples,
                'api_after_kernel_compile_seconds': api_samples[0],
                'api_warm_median_seconds': statistics.median(api_samples[1:]),
                'cli_fresh_process_seconds': cli_time,
                'cli_reported_solver_seconds': cli_result['time'], 'iterations': int(result[0]), 'residual': float(result[2])})
    path = Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
