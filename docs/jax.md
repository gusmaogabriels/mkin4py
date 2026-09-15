# JAX execution and local CLI

Install this checkout with Python 3.11 or later: `python -m pip install .`.
JAX 0.4.38 or later supplies the numerical runtime. SciPy is used only as a
reference in development tests (JAX may itself install SciPy transitively).
For the packaged JAX alpha, use the [README installation instructions](../README.md#instructions).
PyPI currently provides the legacy 1.0 release.

```sh
mkin4py methods --json
mkin4py example --json > model.json
mkin4py solve model.json --json
mkin4py solve model.json --linear-solver dense --json
```

The CLI runs on the user's machine. It reads JSON and never imports user-supplied
Python, sends telemetry, or submits a job to a server. Exit codes: 0 converged,
1 no convergence, 2 invalid input. Nonfinite numbers in a failed result are serialized as `null`. All numeric model fields and units match the
original Python setup API. `pre_exponential` and `activation_energies` are flat
reaction vectors; `stoichiometry` is species by reactions, `gas_species` lists
fixed-pressure species indices, and `concentrations` follows that same order.
Optional `coverage` contains every species, with surface fractions summing to 1.
`seed` makes initial coverage and the eliminated surface species reproducible.

The existing Environment, MKmodel, solver modules, singleton setup pattern,
analytic Jacobian/Hessian and Newton/RK4 algorithm remain. RK4 integrates the
Newton correction flow to a steady state, not physical time. QMR now has a local
JAX implementation; row equilibration and a JAX dense fallback handle QMR
breakdown. Select `linear_solver='dense'` to use dense factorization directly.

Convergence requires the maximum absolute **surface species rate** to be at most
`criteria`. The old comparison against gas rates prevented convergence at zero
net gas flux and has been removed. Nonfinite or exhausted solves return an
explicit failure, with rates computed from the returned coverage. `max_time` is
checked between compiled attempts, so one running attempt can exceed the budget.
The underlying model retains the original single normalized surface balance;
it does not infer multiple site types or additional conservation constraints.

Enable double precision before creating arrays for stiff kinetics:

```python
import jax
jax.config.update('jax_enable_x64', True)
import mkin4py
```

The CLI enables it by default; `--float32` opts into reduced precision. Importing
the library does not change JAX global precision. CI executes the historical
README setup with 64-bit mode enabled, then checks both linear solvers and both
Newton correction variants against its published coverages and rates.

The configuration API is Python. Use the explicit array kernels for composition:

```python
import jax
from mkin4py.solver.derivatives import power_law

# ms: species x reactions; k: reaction constants; coverage: species vector
rates = jax.jit(lambda coverage, ms, k: ms @ (k[:, None] * power_law(coverage, ms)))
```

`power_law` supports `jit`, `vmap`, `jacfwd` and `jacrev`. `evaluate` provides the
original analytic derivatives after eliminating one surface species. The solver's
compiled `attempt` takes all arrays explicitly; changing parameter values does
not reuse stale globals. Its dynamic convergence loop supports forward-mode
transforms, but reverse-mode differentiation through the iterative solve is not
provided. Different array shapes/dtypes or solver methods can require compilation.

```sh
python -m pip install '.[dev]'
python -m pytest -q
python -m benchmarks.timing --repeats 20
```

The benchmark reports lowering, compilation, first execution and repeated warm
execution separately, alongside the Python API and a fresh CLI process. Every
call uses the same initial coverage and synchronizes JAX results. The CLI process
timing includes imports, model setup, compilation, and the solve. `rk4()['time']` includes any compilation and synchronization incurred
by that call. Timings are informational in CI; numerical accuracy is required.
