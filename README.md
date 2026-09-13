# TrafficAssignmentApproaches

**Traffic assignment and continuous network design with a shared C++ core and Python API.**

TrafficAssignmentApproaches implements solvers for the **Traffic Assignment
Problem (TAP)** and the **Continuous Network Design Problem (CNDP)**. It combines
user-equilibrium traffic assignment with bilevel optimization of road-network
capacities, supporting both algorithm research and computational experiments.
The project grew out of graduate-level research at Saint Petersburg State
University.

You can use the solvers through the `traffic_assignment` Python package, the
`traffic_assignment::core` C++ library, or the `tap_solver` and `cndp_solver`
command-line applications. The repository also includes DataFrame adapters,
notebooks, configurable optimization pipelines, and experiment runners.

## Problems and Methods

### Traffic assignment (TAP)

Given a directed road network, link travel-time functions, and an
origin–destination (OD) demand matrix, TAP computes traffic flows at user
equilibrium. Under **Wardrop's first principle**, no traveler can improve their
travel time by changing routes alone.

Links use the Bureau of Public Roads (BPR) travel-time function:

```text
travel_time = free_flow_time * (1 + b * (flow / capacity) ** power)
```

Two TAP approaches are available, both for standalone assignment and as the
lower-level solver in CNDP:

| Approach | Python selector | Implementation |
| --- | --- | --- |
| TAPAS (Traffic Assignment by Paired Alternative Segments) | `"tapas"` | Bush-based assignment using paired alternative segments, with Halley flow shifts and Armijo backtracking. Default approach. |
| Route-based assignment | `"routebased"` | Explicit OD paths with Newton-step flow redistribution and configurable parallel route search. |

TAP results include link flows, link travel times, total travel time, the
Beckmann objective, relative gap, and elapsed solve time. Relative gap measures
how close the computed assignment is to user equilibrium.

### Continuous network design (CNDP)

CNDP chooses continuous link capacities on a fixed network topology. The
upper-level optimizer evaluates capacity changes by solving TAP at the lower
level, accounting for how travelers reroute after an investment.

The implemented design model is:

```text
investment = theta * sum(capacity - lower_bound)
minimize:    total_travel_time_at_user_equilibrium + investment
subject to: lower_bound <= capacity <= upper_bound
            investment <= budget
```

Bounds are **absolute capacities**. CNDP starts at the lower bounds; use the
original capacities as lower bounds to model expansion only. Equal bounds fix a
link's capacity. `theta` is a uniform investment-cost multiplier; the optional
per-link `investment_cost_param` field is currently metadata and does not change
the objective. See the [capacity and budget guide](docs/dataframes.md#cndp-constraints).

Optimization steps run sequentially, each continuing from the capacities left
by the previous step. A pipeline can contain one method or combine several:

| Pipeline type | Available methods |
| --- | --- |
| `nlopt` | NLopt algorithms, including COBYLA, BOBYQA, and ISRES. |
| `optimality_condition` | Sensitivity-based capacity updates using network-design optimality conditions. |
| `optimlib` | OptimLib methods, including differential evolution, particle swarm optimization, and Nelder–Mead. |
| `gradient_descent` | Projected gradient descent with finite-difference, SPSA, or sensitivity-based gradient estimates. |

Budget handling depends on the method: supported NLopt algorithms use an
explicit inequality constraint, while other methods use penalties and/or
capacity projection. Inspect the returned budget usage and effective bounds
when evaluating a design; completing a pipeline does not certify a global
optimum.

## Installation

The core uses **C++20**, **CMake 3.24+**, Ninja, Eigen3, NLopt, and Boost.
The Python package requires **Python 3.11+** and NumPy, with optional pandas and
notebook dependencies.

On Ubuntu/Debian, install the native build prerequisites:

```bash
sudo apt install build-essential cmake ninja-build pkg-config git \
                 libeigen3-dev libnlopt-cxx-dev libboost-all-dev \
                 python3-dev python3-venv
```

CMake fetches OptimLib during configuration and also fetches toml++ when
building the standalone applications, so the initial build needs network access.

### Python package

From the repository root, create an environment using Python 3.11 or newer and
install the package:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

This builds the native extension through scikit-build-core and pybind11.
Optional extras provide DataFrame adapters, notebooks, and test dependencies:

```bash
python -m pip install -e '.[dataframes]'  # pandas input adapters
python -m pip install -e '.[examples]'    # pandas, plotting, and Jupyter
python -m pip install -e '.[test]'        # pytest and pandas
```

### Standalone C++ applications

Configure and build from the repository root:

```bash
cmake --preset linux-release
cmake --build --preset build-release
```

The executables are written to `build/linux-release/`. Debug and RelWithDebInfo
configure/build presets are also available in [CMakePresets.json](CMakePresets.json).

## Python Quick Start

This example solves both problems on a small network constructed in memory.
It requires no dataset files or pandas installation.

```python
import traffic_assignment as ta

# Two alternatives from node 0 to node 1: a direct link and a route via node 2.
network = ta.network_from_arrays(
    "TwoRoutes", init_node=[0, 0, 2], term_node=[1, 2, 1],
    capacity=10.0, free_flow_time=[2.0, 1.0, 1.0], b=1.0, power=1.0,
    demand=[[0.0, 10.0], [0.0, 0.0]], n_nodes=3,
)

# TAP: find user-equilibrium flows at the original capacities.
assignment = ta.solve_tap(network, approach="tapas")
print("Link flows:", assignment.flows)  # [5, 5, 5]
print("Total travel time:", assignment.total_travel_time)  # 30
print("Relative gap:", assignment.relative_gap)

# CNDP: allow each link to expand from capacity 10 to at most 20.
constraints = ta.constraints_from_arrays(network, lower=10.0, upper=20.0)
design = ta.solve_cndp(
    network,
    [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=100)],
    constraints=constraints,
    theta=0.1,
    budget=2.0,
)
print("Designed capacities:", design.capacities)
print("Final total travel time:", design.total_travel_time)
print("Investment used:", design.budget)
print("Objective (travel time + investment):", design.objective)
```

`solve_tap` and `solve_cndp` return dataclasses containing metrics and NumPy
arrays. Solving updates the supplied network, and CNDP changes its capacities;
construct a fresh network for independent scenarios. Result arrays are snapshots.
The standalone [two-route example](examples/two_routes.py) also demonstrates
typed TAP options.

For more control, use `TapOptions`, `CndpOptions`, and the exposed solver objects
(`Network`, `TapasApproach`, `RouteBasedApproach`, `BilevelCND`). Python and the
executables share the same native implementation: TAP arithmetic uses
`long double`, while optimizer interfaces and returned NumPy arrays use
double precision. See the [architecture guide](docs/architecture.md).

CNDP file outputs are off by default in Python. Pass a `metrics=` configuration
to enable trace CSV, metadata JSON, and summary CSV output;
`final_diagnostics=True` additionally enables the final diagnostic and solution
CSV files.

## Inputs and Notebooks

Both solvers consume the same native network. Build it from pandas DataFrames,
NumPy arrays, or preprocessed benchmark CSV files. CNDP also requires capacity
constraints, which can be prepared from DataFrames, arrays, or a CSV file.

### DataFrames

With the `dataframes` extra installed and your input tables prepared:

```python
import traffic_assignment as ta

# links_df: per-link columns; demand_df: square matrix labelled by zone ID
network = ta.network_from_dataframes(
    "MyNetwork", links_df, demand_df, node_index_base=1,
)
result = ta.solve_tap(network)

# constraints_df: link_index, lower_bound, upper_bound (one row per link)
constraints = ta.constraints_from_dataframe(network, constraints_df, node_index_base=1)
design = ta.solve_cndp(
    network,
    [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=100)],
    constraints=constraints,
    budget=10000.0,
)
```

The [DataFrame input guide](docs/dataframes.md) documents columns, node and link
indexing, validation, and capacity constraints. The
[solver input guide](docs/inputs.md) covers the shared numerical contract and
CSV formats.

### Self-contained notebooks

Both Sioux Falls notebooks construct their input DataFrames from embedded
benchmark values. They need no dataset files or downloads after installing the
`examples` extra:

- [TAP with DataFrames](examples/siouxfalls_tap.ipynb): construct the network,
  solve equilibrium, inspect flows and congestion, compare approaches, and
  change demand.
- [CNDP with DataFrames](examples/siouxfalls_cndp.ipynb): construct capacity
  constraints, solve the baseline TAP and network design, check the budget and
  bounds, and compare capacity and flow changes.

Start `jupyter lab examples/` and select a kernel with the package installed.

### Benchmark datasets

The benchmark workflows use networks from the
[Transportation Networks for Research](https://github.com/bstabler/TransportationNetworks)
collection. Solver configurations are provided for Sioux Falls, Anaheim,
Barcelona, Chicago Sketch, Sydney, Winnipeg, and other scenarios under
[`configs/`](configs/).

Dataset files are not tracked in this repository. File-based workflows expect
preprocessed `*_net.csv` and `*_trips.csv` files, plus capacity constraints for
CNDP, under `data/TransportationNetworks/<Dataset>/`. The loader reads these
CSV formats, not raw `.tntp` files. The
[dataset-processing notebook](scripts/dataset_processor.ipynb) contains the
TNTP conversion workflow; see [solver inputs](docs/inputs.md#file-loading-and-indexing)
for the expected formats.

In Python, load prepared data explicitly:

```python
import traffic_assignment as ta

network = ta.load_network("SiouxFalls")
constraints = ta.load_constraints(ta.default_constraints_path("SiouxFalls"))
assignment = ta.solve_tap(network)
```

Python dataset discovery walks up from the working directory looking for
`data/TransportationNetworks`. Override it with `data_root=` on the loaders or
the `TRAFFIC_ASSIGNMENT_DATA` environment variable.

## Command-Line Usage and Configuration

After building the applications and preparing the dataset files, run from the
repository root:

```bash
./build/linux-release/tap_solver --config configs/tap.siouxfalls.toml
./build/linux-release/cndp_solver --config configs/cnd.siouxfalls.toml
./build/linux-release/tap_solver --help
./build/linux-release/cndp_solver --help
```

`cndp_solver` resolves settings in order of increasing priority:
**defaults → TOML file → environment variables → CLI arguments**. Select a
config with `--config` or `CND_CONFIG`; environment overrides use the `CND_`
prefix.

A CNDP config specifies a network, solver settings, and at least one
`[[pipeline]]` entry. For example, this pipeline starts with COBYLA and then
applies optimality-condition updates:

```toml
[network]
dataset = "SiouxFalls"

[solver]
approach = "Tapas"
approach_alpha = 1e-14
budget_function_multiplier = 5.0
budget_upper_bound = 100_000.0

[[pipeline]]
type = "nlopt"
algorithm = "LN_COBYLA"
max_iterations = 100
tolerance = 1e-4

[[pipeline]]
type = "optimality_condition"
max_iterations = 20
```

The same sequence can be expressed as a list of `ta.step(...)` objects in
Python. See [the reference CNDP config](configs/cnd.siouxfalls.toml) for an
example with metrics enabled.

CLI arguments override both the config and environment settings:

```bash
CND_ROUTE_THREADS=1 ./build/linux-release/cndp_solver \
    --config configs/cnd.siouxfalls.toml \
    --route-threads 2 \
    --max-standard-iterations 150
```

## Experiments and Analysis

Experiment runners launch the standalone applications and organize configs,
metrics, and plots in run folders under `performance_results/`. They require
the built executables and prepared datasets; the Python `examples` extra
provides the plotting dependencies.

```bash
# Run one TAP or CNDP experiment.
python scripts/run_experiment.py tap --dataset SiouxFalls
python scripts/run_experiment.py cndp --config configs/cnd.siouxfalls.toml

# Compare CNDP scenarios defined in a TOML manifest.
python scripts/run_cndp_comparison.py configs/comparisons/siouxfalls.toml
```

Comparison manifests under [`configs/comparisons/`](configs/comparisons/)
define shared defaults and per-scenario pipelines. The comparison runner
supports scenario selection and repeated runs. Use
[`scripts/plot_cndp_comparison.py`](scripts/plot_cndp_comparison.py) with
`--run-dir` to produce objective/convergence plots, budget-usage plots, and
summary tables for a run folder. Each script exposes its options through
`--help`.

For TAP comparisons, examine relative gap and runtime together. For CNDP,
compare the objective, total travel time, investment, and capacity feasibility
under the same demand, bounds, budget, and TAP tolerance.

## Repository Structure

| Path | Contents |
| --- | --- |
| [`cpp/include/traffic_assignment/`](cpp/include/traffic_assignment/) | Public C++ API for networks, options, solvers, and results. |
| [`cpp/src/`](cpp/src/) | Compiled core and private TAP/CNDP implementations under `detail/`. |
| [`cpp/bindings/`](cpp/bindings/) | pybind11 bindings. |
| [`cpp/apps/`](cpp/apps/) | Standalone solvers and runtime configuration. |
| [`python/traffic_assignment/`](python/traffic_assignment/) | Python API, input adapters, and dataset loaders. |
| [`examples/`](examples/) | Small Python example and Sioux Falls notebooks. |
| [`configs/`](configs/) | TAP/CNDP settings and comparison manifests. |
| [`scripts/`](scripts/) | Dataset preparation, experiment runners, and analysis. |
| [`tests/`](tests/) | Python tests and native C++ API tests. |

The Python extension and standalone applications link the same C++ library,
`traffic_assignment::core`. See the [architecture guide](docs/architecture.md)
for C++ integration, ownership, build options, and current implementation limits.

## Development Checks

Run the Python tests and the small TAP example after installing the package:

```bash
python -m pip install -e '.[test]'
python -m pytest tests/
python examples/two_routes.py
```

Enable and run the native API tests with:

```bash
cmake --preset linux-release -DBUILD_NATIVE_TESTS=ON
cmake --build --preset build-release
ctest --test-dir build/linux-release --output-on-failure
```
