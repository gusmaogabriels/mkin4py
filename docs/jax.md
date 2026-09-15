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
transforms. The source checkout also provides `steady_state` for forward- and
reverse-mode sensitivities of converged solutions. Different array shapes/dtypes
or solver methods can require compilation.

## Steady-state gradients

`steady_state` uses the original Newton/RK4 numerical attempt, with
[implicit differentiation](https://docs.jax.dev/en/latest/_autosummary/jax.lax.custom_root.html)
of the reduced surface balance. It returns a one-dimensional coverage vector and
supports `jit`, `grad`, `jacfwd`, `jacrev`, `hessian` and `vmap`. Model configuration
stays outside the transformed function; pass changing physical parameters as arrays:

```python
import jax
import jax.numpy as jnp
from mkin4py.cli import example, load_model
from mkin4py.solver.solve import steady_state

jax.config.update('jax_enable_x64', True)
model = load_model(example())
maps = model.maps

def adsorbed_fraction(k):
    coverage = steady_state(model.coverage, model.ms, k, maps['surface'],
                            maps['xsurface'], maps['ndof'], criteria=1e-12)
    return coverage[2]

value, gradient = jax.jit(jax.value_and_grad(adsorbed_fraction))(jnp.array([2., 1.]))
# value = 2/3; gradient = [1/9, -2/9]
```

Rate constants and fixed gas entries in the initial coverage can be differentiated;
the initial surface guess is not a physical parameter. Hold integer stoichiometry
and index maps fixed. Chain the rate constants through the Arrhenius expression
to differentiate temperature, activation energies or pre-exponential factors.
The callable can be composed into a local Optinpy objective.

Sensitivities require an isolated converged root with a nonsingular reduced
Jacobian. This interface normalizes each residual row by its largest kinetic
coefficient before applying `criteria`; uniformly slow rates cannot make an
incorrect initial coverage pass the root check. The original `rk4` retains its
absolute physical-rate criterion. Invalid dynamic solver controls return NaN.
Sensitivities describe the selected steady-state branch; they do not differentiate
restart decisions or jumps between multiple steady states. Unconverged outputs and
their sensitivities are NaN; a singular root can have finite coverage but invalid
sensitivities. Check finiteness before using a result. The pure function has no
random restarts or Python wall-time limit. The original `rk4` driver and its result
dictionary remain available for restart/status handling.

## Local usage counts

Usage accounting is opt-in and local. Choose a database path to enable it:

```sh
export MKIN4PY_USAGE_DB="$HOME/mkin4py-usage.sqlite3"
mkin4py methods --json
mkin4py usage --json
```

PowerShell: `$env:MKIN4PY_USAGE_DB = "$HOME/mkin4py-usage.sqlite3"`.
The JSON report contains daily counts by package version, command and exit code,
plus command duration excluding interpreter startup and package imports. It stores
no command arguments, model contents, file paths or user/device identifiers and sends nothing
over the network. `usage` does not count itself. Unset the environment variable to
stop recording; remove your database file to clear counts. Storage errors leave
the command's result and exit code intact. These local counts are separate from
public downloads and MCP request attempts; they are not automatically uploaded
to a dashboard. `steady_state` and `usage` are source additions after the existing
2.0.0a1 GitHub artifacts; use the source installation at the top of this page.

## Verification and timing

```sh
python -m pip install '.[dev]'
python -m pytest -q
python -m benchmarks.timing --repeats 20
python -m benchmarks.sensitivities --repeats 20
```

The benchmark reports lowering, compilation, first execution and repeated warm
execution separately, alongside the Python API and a fresh CLI process. Every
call uses the same initial coverage and synchronizes JAX results. The CLI process
timing includes imports, model setup, compilation, and the solve. `rk4()['time']` includes any compilation and synchronization incurred
by that call. Timings are informational in CI; numerical accuracy is required.
