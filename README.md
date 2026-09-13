# TrafficAssignmentApproaches

This repository presents a set of algorithmic implementations designed to solve the classical **Traffic Assignment Problem (TAP)**. The work stems from graduate-level research conducted at Saint Petersburg State University and includes both standard and novel approaches to equilibrium-based flow distribution on transportation networks.

## Problem Overview

The Traffic Assignment Problem (TAP) involves determining an equilibrium distribution of traffic flows across a transportation network, such that no traveler can reduce their travel time by unilaterally changing routes. This condition corresponds to **Wardrop's first principle of user equilibrium**. TAP is fundamental to transportation systems analysis and infrastructure planning, as it captures the interaction between route choices and network congestion.

Mathematically, TAP is often formulated as a convex optimization problem with nonlinear cost functions, representing travel time as a function of link load. Solving this problem with high precision is particularly important in bilevel applications such as network design and congestion pricing.

## Implemented Approaches

The repository includes the following classes of traffic assignment methods:

- **Path-based Methods**
  Manage flow distributions across explicit paths for each origin–destination (OD) pair. This approach supports high-precision computation and direct access to route-level flow information. The implementation follows a route-equilibration strategy grounded in recent academic literature.


- **TAPAS (Traffic Assignment by Paired Alternative Segments)**
  A refinement of bush-based algorithms that employs paired segments with common endpoints but disjoint paths. TAPAS supports both line-search and Newton-step flow updates, allowing for fast and accurate convergence. This method is especially effective in high-precision applications.

- **Demand-based Method (Novel Contribution)**
  A newly developed algorithm introduced by the author, which models dynamic redistribution of OD demand through the network. The approach uses a pump-out mechanism to iteratively reassign unmet demand. While preliminary results indicate lower accuracy and slower convergence, the method presents a promising direction for further investigation.

## Technical Implementation

- **Programming Languages**
  - Core algorithm implementations: **C++20**
  - Data preprocessing and experiment orchestration: **Python**

- **Build System**
  CMake (>= 3.24) with Ninja, using system-installed Eigen3, NLopt, and Boost.

The shared C++ library target is `traffic_assignment::core`. Its public headers
are in `cpp/include/traffic_assignment/`; algorithm templates are private under
`cpp/src/detail/`. The Python extension (`cpp/bindings/`) and standalone
applications (`cpp/apps/`) both use this compiled library. See the
[architecture guide](docs/architecture.md) for the API, ownership model, and build options.

## Build and Run

### Prerequisites (Ubuntu/Debian)

```bash
sudo apt install build-essential cmake ninja-build pkg-config \
                 libeigen3-dev libnlopt-cxx-dev libboost-all-dev
```

(`OptimLib` is fetched automatically by CMake `FetchContent`; `toml++` is fetched
when building the standalone applications.)

### Configure

```bash
cmake --preset linux-release           # or linux-debug, linux-relwithdebinfo
```

### Build

```bash
cmake --build --preset build-release   # or build-debug, build-relwithdebinfo
```

### Run

```bash
./build/linux-release/cndp_solver --config configs/cnd.siouxfalls.toml
./build/linux-release/cndp_solver --help
```

A separate standalone TAP solver is also produced:

```bash
./build/linux-release/tap_solver --config configs/tap.siouxfalls.toml
```

## Python Package

The solvers are also available as a Python package (`traffic_assignment`) with
pybind11 bindings to the same C++ core (same `long double` numerics as the
executables). Install it in editable mode from the repository root:

```bash
pip install -e .            # builds the native extension via scikit-build-core
pip install -e .[test]      # with pytest, then run: pytest tests/
```

Prepare a network, then pass it to TAP or CNDP. DataFrame adapters are the main
workflow for tabular inputs; arrays and explicit dataset loaders are also
supported. See [Providing solver inputs](docs/inputs.md) for their shared contract.

For a benchmark dataset:

```python
import traffic_assignment as ta

# User equilibrium (TAP); load from data/TransportationNetworks
network = ta.load_network("SiouxFalls")
result = ta.solve_tap(network, approach="tapas")
print(result.relative_gap, result.total_travel_time, result.flows)

# Bilevel CNDP with a pipeline of optimization steps
constraints = ta.load_constraints(ta.default_constraints_path("SiouxFalls"))
design = ta.solve_cndp(
    network,
    [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=100)],
    constraints=constraints,
    budget=10000.0,
)
print(design.objective, design.capacities)

# Networks can also be built programmatically from arrays
network = ta.network_from_arrays(
    "TwoRoutes", init_node=[0, 0, 2], term_node=[1, 2, 1],
    capacity=10.0, free_flow_time=[2.0, 1.0, 1.0], b=1.0, power=1.0,
    demand=[[0.0, 10.0], [0.0, 0.0]], n_nodes=3,
)
```

`solve_tap` / `solve_cndp` return result dataclasses with NumPy arrays, using
objective values and metrics computed by the native library. A self-contained
example is available in [`examples/two_routes.py`](examples/two_routes.py).
The typed `TapOptions` and `CndpOptions` classes expose native solver settings; the
underlying objects (`Network`, `TapasApproach`, `RouteBasedApproach`,
`BilevelCND`, ...) are exposed for object-level workflows. CNDP file outputs
(trace CSV, metadata JSON, summary CSV) are off by default and enabled by
passing a `metrics=` config. Dataset discovery walks up from the working
directory looking for `data/TransportationNetworks`; override with the
`data_root=` argument or the `TRAFFIC_ASSIGNMENT_DATA` environment variable.
Dataset-name and constraint-path shortcuts in solver calls remain available
for compatibility. Explicit loaders make the chosen inputs visible before solving.

### DataFrame inputs and notebooks

Install the optional pandas adapters with `pip install -e '.[dataframes]'`, or
install the notebook environment with `pip install -e '.[examples]'`.

```python
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
indexing, validation, and CNDP capacity/budget semantics. Both Sioux Falls
notebooks construct their input DataFrames from embedded benchmark values;
they need no dataset files or downloads after installing the package:

- [TAP with DataFrames](examples/siouxfalls_tap.ipynb): construct the network,
  solve equilibrium, inspect flows and congestion, compare approaches, and
  change demand.
- [CNDP with DataFrames](examples/siouxfalls_cndp.ipynb): construct capacity
  constraints, solve the baseline TAP and network design, check the budget and
  bounds, and compare capacity and flow changes.

Start `jupyter lab examples/` and select a kernel with the package installed.

## Layered Runtime Config (defaults → config → env → CLI)

`cndp_solver` supports layered configuration:

- Defaults in code
- TOML config file (`--config ...` or `CND_CONFIG`)
- Environment variable overrides (`CND_*`)
- CLI overrides (highest priority)

Example config file: [`configs/cnd.siouxfalls.toml`](configs/cnd.siouxfalls.toml).

Run using config only:

```bash
./build/linux-release/cndp_solver --config configs/cnd.siouxfalls.toml
```

Run with environment override:

```bash
CND_ROUTE_THREADS=1 ./build/linux-release/cndp_solver --config configs/cnd.siouxfalls.toml
```

Run with CLI override (wins over env + file):

```bash
CND_ROUTE_THREADS=1 ./build/linux-release/cndp_solver \
    --config configs/cnd.siouxfalls.toml \
    --route-threads 2 \
    --max-standard-iters 150
```

Show all available options:

```bash
./build/linux-release/cndp_solver --help
```

## Benchmark Networks

The following standardized test networks are used to evaluate algorithmic performance:

- **Sioux Falls Network**
  24 nodes, 24 zones, 76 links

- **Anaheim Network**
  416 nodes, 38 zones, 914 links

These networks are sourced from the open-access [Transportation Networks for Research](https://github.com/bstabler/TransportationNetworks) repository and conform to the TNTP format. Each instance includes a demand matrix and BPR-formulated travel time functions.

## Experimental Results (Summary)

The implemented algorithms were evaluated with respect to solution precision (measured by the Relative Gap, RGAP) and computation time. Key findings include:

- The **TAPAS-NewtonStep** method achieved the best performance, reaching machine-precision accuracy (RGAP ≈ 10⁻¹⁵) in the shortest runtime.
- The **TAPAS-LineSearch** variant also performed favorably, with slightly longer convergence times but similarly high precision.
- The **Path-based** implementation achieved acceptable precision (RGAP ≈ 10⁻¹²), though it incurred higher computational cost due to the maintenance of explicit paths.
- The **Demand-based** method, in its current form, did not achieve comparable precision and exhibited slower convergence, suggesting it may benefit from further optimization and theoretical refinement.
