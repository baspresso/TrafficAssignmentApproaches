# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

C++20 project using CMake + Ninja + MinGW + vcpkg. Dependencies: Eigen3, Boost.Multiprecision, NLopt.

```bash
# Configure
cmake --preset mingw-vcpkg-release    # or mingw-vcpkg-debug / mingw-vcpkg-relwithdebinfo

# Build
cmake --build --preset build-release  # or build-debug / build-relwithdebinfo

# Run (from project root)
c
./build/mingw-vcpkg-release/main.exe --help
```

Build output goes to `build/<preset-name>/`. There are no automated tests — validation is done by running experiments and inspecting metrics outputs.

## Configuration

Layered config system with precedence: **defaults → INI file → environment variables → CLI args**.

```bash
# All three layers combined
./build/mingw-vcpkg-release/main.exe \
  --config ./configs/cnd.siouxfalls.ini \
  --route-threads 2 \
  --max-standard-iters 150
```

Environment variable prefix: `CND_` (e.g., `CND_ROUTE_THREADS=1`).

Key config options: `approach` (RouteBased/Tapas/DemandBased), `shift_method`, `max_standard_iterations`, `max_optimality_condition_iterations`, `nlopt_algorithm` (LN_COBYLA, LN_BOBYQA, etc.), `route_search_threads`, `output_root`, `enable_trace`.

## Architecture

### Domain Layers

```
main.cpp              — Config parsing, wiring, execution
    ↓
cnd/BilevelCND        — Upper-level NLopt optimization over link capacities
    ↓
tap/algorithms/       — Lower-level TAP solvers (user equilibrium)
    ↓
tap/core/Network      — Network topology, Dijkstra, OD pairs
tap/data/Link         — BPR delay function: t(x) = t0*(1 + b*(x/c)^p)
```

### Traffic Assignment Approaches (`include/tap/algorithms/`)

Three solver implementations, all implementing `TrafficAssignmentApproach`:

1. **RouteBasedApproach** — Explicit path enumeration; factory in `RouteBasedShiftMethodFactory`. Shift methods: `NewtonStep`, `Krylatov2023`. Supports parallel route search via `route_search_threads`.

2. **TapasApproach** — Bush-based TAPAS algorithm; fastest, machine precision (~1e-15). Shift methods: `NewtonStep` (best), `LineSearch`, `AdvancedGradientDescent`.

3. **PumpOutDemandBasedApproach** — Experimental demand redistribution approach; lower accuracy.

### Bilevel CND (`include/cnd/`)

- `BilevelCND.h` — Upper-level: iterates NLopt optimizer ↔ TAP solver to design link capacities minimizing system travel cost under budget constraints.
- `DirectedConstraintLoader.h` — Reads per-link capacity constraints from CSV.
- `CndStatisticsRecorder.h` — Records JSON metadata + CSV metrics for each optimization run.

### Data Format (TNTP)

Network files (`*_net.csv`): `init_node, term_node, capacity, length, free_flow_time, b, power, speed, toll, link_type`

Trip matrices (`*_trips.csv`): OD demand matrix. Datasets live in `data/TransportationNetworks/`.

### Statistics Output

TAP convergence: `StatisticsRecorder` (`include/tap/utils/`). CND metrics: `CndStatisticsRecorder`. Both write to `output_root` (default: `performance_results/`).

## Key File Locations

- `src/main.cpp` — Entry point; all config/wiring logic
- `include/cnd/BilevelCND.h` — Bilevel optimization core
- `include/tap/algorithms/common/TrafficAssignmentApproach.h` — Base interface for TAP solvers
- `include/tap/core/Network.h` — Network data structure + Dijkstra
- `include/tap/data/Link.h` — BPR function implementation
- `configs/cnd.siouxfalls.ini` — Reference config example
- `scripts/` — Python notebooks for constraint generation and results plotting

## Code Conventions

- Headers-only design: most implementation is in `.h` files under `include/`
- Factory pattern used for algorithm selection (shift methods, approaches)
- Templates parameterized on numeric type (typically `double`, `long double`, or Boost high-precision)
- Algorithms use `long double` / Boost.Multiprecision for precision-critical calculations
