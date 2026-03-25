# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

C++20 project using CMake + Ninja + MinGW + vcpkg. Dependencies: Eigen3, Boost.Multiprecision, NLopt, OptimLib (FetchContent).

```bash
# Configure
cmake --preset mingw-vcpkg-release    # or mingw-vcpkg-debug / mingw-vcpkg-relwithdebinfo

# Build
cmake --build --preset build-release  # or build-debug / build-relwithdebinfo

# Run
./build/mingw-vcpkg-release/cndp_solver.exe --config configs/cnd.siouxfalls.ini
./build/mingw-vcpkg-release/cndp_solver.exe --help
```

Two executables: `cndp_solver.exe` (CNDP solver) and `tap_solver.exe` (standalone TAP solver). Build output goes to `build/<preset-name>/`. No automated tests — validation is done by running experiments and inspecting metrics outputs.

## Architecture

### Domain Layers

```
cndp_solver.cpp                   — Config parsing, wiring, execution
    ↓
cnd/BilevelCND                    — Upper-level optimization over link capacities
    ↓
cnd/OptimizationPipeline          — Sequential [step.N] execution
    ├── NloptOptimizationStep     — NLopt algorithms (COBYLA, BOBYQA, ISRES, ...)
    ├── OptimalityConditionStep   — Optimality condition-based iteration
    └── OptimlibOptimizationStep  — Population-based metaheuristics (DE, PSO, ...)
    ↓
cnd/CndOptimizationContext        — Objective evaluation (TotalTravelTime + BudgetCost)
    ↓
tap/algorithms/                   — Lower-level TAP solvers (user equilibrium)
    ↓
tap/core/Network                  — Network topology, Dijkstra, OD pairs
tap/data/Link                     — BPR delay function: t(x) = t0*(1 + b*(x/c)^p)
```

### Traffic Assignment Approaches (`include/tap/algorithms/`)

Two solver implementations, both implementing `TrafficAssignmentApproach<T>`:

1. **RouteBasedApproach** — Explicit path enumeration; factory-based shift methods (`NewtonStep`, `Krylatov2023`). Supports parallel route search via `route_search_threads`.

2. **TapasApproach** — Bush-based TAPAS algorithm (TAsK-style); fastest, reaches machine precision (~1e-15 RGAP). Uses Halley step with Armijo backtracking, stall detection with bush restart, random origin shuffle.

### Bilevel CNDP (`include/cnd/`)

- **BilevelCND** — Pipeline-based constructor only: takes `std::vector<OptimizationStepConfig>`. Iterates optimization steps ↔ TAP solver to design link capacities minimizing system travel cost under budget constraints.
- **CndOptimizationContext** — Shared evaluation context with crash recovery (SIGSEGV handler). Includes `GetAlgorithmName()` and `SupportsHardBudgetConstraint()`.
- **OptimizationStepFactory** — Creates steps from config: `"nlopt"`, `"optimality_condition"`, `"optimlib"`.
- **CndStatisticsRecorder** — Records trace CSV (quality vs time), metadata JSON, summary CSV (append-only). Summary CSV has a `Scenario` column.
- **DirectedConstraintLoader** — Reads per-link capacity constraints from CSV.

## Configuration

Layered config system with precedence: **defaults → INI file → environment variables → CLI args**.

```bash
./build/mingw-vcpkg-release/cndp_solver.exe \
  --config ./configs/cnd.siouxfalls.ini \
  --route-threads 2 \
  --max-standard-iters 150
```

Environment variable prefix: `CND_` (e.g., `CND_ROUTE_THREADS=1`). CLI metrics args prefixed with `metrics_` (e.g., `--metrics_scenario_name`).

### INI Sections

**`[run]`** — Algorithm selection: `dataset`, `approach` (RouteBased/Tapas), `shift_method` (RouteBased only), `approach_alpha`, `route_search_threads`, `tapas_mu`, `tapas_v`, `budget_function_multiplier`, `budget_upper_bound`, `constraints_file`.

**`[metrics]`** — Output control: `enable_trace`, `output_root`, `write_metadata_json`, `write_summary_csv`, `scenario_name`, `append_dataset_subdir`.

**`[step.N]`** (N = 0, 1, 2, ...) — Pipeline steps executed sequentially:
- `type = nlopt` — `algorithm` (LN_COBYLA, LN_BOBYQA, GN_ISRES, etc.), `max_iterations`, `tolerance`, optional `local_algorithm`
- `type = optimality_condition` — `max_iterations`
- `type = optimlib` — `algorithm` (DE, DE_PRMM, PSO, PSO_DV, NM, GD), `max_iterations`, `population_size` (0 = auto)

At least one `[step.N]` section is required; `cndp_solver.cpp` errors if none found.

## Key File Locations

- `src/cndp_solver.cpp` — CNDP entry point; all config/wiring logic
- `src/tap_solver.cpp` — Standalone TAP solver entry point
- `include/cnd/BilevelCND.h` — Bilevel optimization core
- `include/cnd/CndOptimizationContext.h` — Shared context, objective evaluation, crash recovery
- `include/cnd/OptimizationPipeline.h` — Sequential step execution
- `include/cnd/steps/` — Step implementations (Nlopt, OptimalityCondition, OptimLib)
- `include/tap/algorithms/common/TrafficAssignmentApproach.h` — Base TAP interface
- `include/tap/core/Network.h` — Network data structure + Dijkstra
- `include/tap/data/Link.h` — BPR function implementation
- `include/common/ConfigUtils.h` — Config parsing utilities
- `configs/cnd.siouxfalls.ini` — Reference CNDP config
- `configs/scenarios/` — Auto-generated scenario configs for comparison runs

## Scripts

- `scripts/run_experiment.py` — Unified experiment runner for TAP and CNDP. Creates self-contained run folders with CSV outputs and publication-ready plots. Subcommands: `tap` (standalone TAP), `cndp` (single CNDP scenario).
- `scripts/run_cndp_comparison.py` — Batch runner: generates per-scenario INI configs, launches experiments, organizes outputs into self-contained run folders under `performance_results/{Dataset}/runs/`.
- `scripts/plot_cndp_comparison.py` — Publication-ready analysis: objective vs time plots, convergence, sensitivity, summary tables (CSV + LaTeX). Supports `--run-dir` mode and multi-run averaging.

## Data Format (TNTP)

Network files (`*_net.csv`): `init_node, term_node, capacity, length, free_flow_time, b, power, speed, toll, link_type`. Trip matrices (`*_trips.csv`): OD demand matrix. Datasets live in `data/TransportationNetworks/` (SiouxFalls, Anaheim, Barcelona).

## Code Conventions

- **Headers-only design:** Most implementation is in `.h` files under `include/`.
- **Templates on numeric type:** `template<typename T>` parameterized (typically `double`, `long double`, or Boost high-precision).
- **Factory pattern:** Used for algorithm selection — `RouteBasedShiftMethodFactory`, `OptimizationStepFactory`.
- **Pipeline-based optimization:** All CNDP configs must use `[step.N]` sections; the legacy single-algorithm constructor was removed.
