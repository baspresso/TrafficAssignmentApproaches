# Library architecture

The reusable computation lives in the `traffic_assignment_core` static C++ library.
Both command-line applications and the Python extension link this target, available
to CMake consumers as `traffic_assignment::core`.

```text
Python API                         Standalone applications
python/traffic_assignment/          cpp/apps/
          |                            |
cpp/bindings/                           |
          +-------------+--------------+
                        |
       cpp/include/traffic_assignment/   public C++ API
                        |
       cpp/src/{network,assignment,design}.cpp
                        |
       cpp/src/detail/                  private algorithm templates
```

## Responsibilities

| Layer | Responsibility |
| --- | --- |
| Public C++ headers | Network ownership, input data, options, solvers, result snapshots |
| Compiled implementation | Graph construction, cost functions, TAP, CNDP, metrics and numerical results |
| Private templates | Existing assignment algorithms and optimization pipelines, instantiated with `long double` |
| Bindings | Python/C++ type conversion, NumPy copies, exception translation, GIL release |
| Python package | Dataset discovery, input normalization, convenient calls and result dataclasses |
| Applications | TOML/environment/CLI configuration, console output, CLI result protocol |

The public headers depend only on the C++ standard library. Eigen, Boost, NLopt,
OptimLib, and the algorithm headers are private implementation dependencies.
`toml++` is required only when building the standalone applications. There is no
Python or pybind11 dependency in the core target.

The existing optional C++ statistics recorders remain internal to the library.
Their output is explicitly enabled by callers; the library does not print the
CLI `[RESULT]` protocol. Research runners retain their existing executable names,
output paths, and structured result fields.

## Public C++ API

- `network.hpp`: `NetworkData`, `LinkData`, `Network`, `Link`, CSV loading.
- `options.hpp`: `TapOptions`, `CndpOptions`, `StepConfig`, `MetricsConfig`, `LinkConstraint`.
- `results.hpp`: independent `TapResult` and `CndpResult` value snapshots.
- `solve.hpp`: approach construction, `SolveTap`, and `BilevelCND`.

Network and solver implementations are hidden behind private implementation
pointers. A solver retains its network; a CNDP solver also retains its TAP
approach. A link view retains native network storage. These lifetime guarantees
apply to both C++ and Python callers.

Use the library in a CMake project with:

```cmake
set(BUILD_CPP_APPS OFF CACHE BOOL "" FORCE)
set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
add_subdirectory(path/to/TrafficAssignmentApproaches traffic_assignment_build)
target_link_libraries(my_application PRIVATE traffic_assignment::core)
```

For example:

```cpp
#include <traffic_assignment/solve.hpp>

namespace ta = traffic_assignment;

ta::NetworkData data;
data.name = "TwoRoutes";
data.number_of_nodes = 3;
data.number_of_zones = 2;
data.links = {{0, 1, 10, 0, 2, 1, 1},
              {0, 2, 10, 0, 1, 1, 1},
              {2, 1, 10, 0, 1, 1, 1}};
data.demand = {{0, 10}, {0, 0}};
auto network = std::make_shared<ta::Network>(data);
auto result = ta::SolveTap(network);
```

This first library version supports consumption through `add_subdirectory`.
It does not install a separate C++ SDK or promise a stable binary ABI.

## Python compatibility

The existing top-level functions and classes remain available, including
`solve_tap`, `solve_cndp`, `network_from_arrays`, `TapasApproach`, and `BilevelCND`.
High-level calls still return Python dataclasses with NumPy arrays.
`CndpResult.network` is the live network for compatibility; its array fields are
independent snapshots and do not change when that network is reused.

`alpha=` remains an alias for `relative_gap_tolerance=`. Typed TAP options can be
passed as the `approach` argument:

```python
options = ta.TapOptions()
options.approach = "routebased"
options.relative_gap_tolerance = 1e-10
result = ta.solve_tap(network, approach=options)
```

`CndpOptions` is also available for object-level workflows:

```python
options = ta.CndpOptions()
options.budget_upper_bound = 10000
solver = ta.BilevelCND(
    network, ta.make_approach(network), constraints,
    [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=100)], options,
)
native_result = solver.compute_network_design()
```

The native `compute_network_design()` call now returns a snapshot. It is silent
and has file output disabled by default, including for the original positional
constructor. Use the existing output setters or `CndpOptions` to enable output.
The old positional constructor accepts already prepared constraints; the typed
constructor and high-level function default to native insensitive-link filtering.
CNDP rejects an approach belonging to a different network.

The Python CNDP wrapper no longer filters links or recomputes the budget and
objective. It packages values returned from the native final evaluation,
including the bounds after filtering. Caller-provided constraints are copied
inside the native library.

`io.py` handles array conversion and constraint-file loading. Previous imports
from `datasets.py` continue to work. The `_core.pyi` declarations describe the
native interface for editors and type checkers.

`dataframes.py` adds optional pandas adapters, exported as
`network_from_dataframes` and `constraints_from_dataframe`. They validate numeric
inputs, align demand axes by zone label, and align links/constraints by an explicit
zero-based `link_index` before using the existing array/native constructors.
pandas is imported only when an adapter is called. See the
[DataFrame input guide](dataframes.md) for the schema and the Sioux Falls notebooks
in `examples/` for complete workflows.

## State, precision, and scope

Reuse an approach object to warm-start its internal state. A new TAPAS approach
starts a new assignment. `reset=True` on `solve_tap` resets both the supplied
approach and its network. Serialize access to mutable network/solver state;
releasing the GIL does not establish native thread safety.

TAP arithmetic uses `long double`; optimizer interfaces and CNDP summary scalars
retain their existing `double` representation. Python receives `float64` arrays.
Copies at the Python boundary are intentional so returned arrays remain valid
after native objects are destroyed or reused.

This reorganization preserves the numerical algorithms and CSV data conventions.
A returned result does not certify convergence or feasibility. Unreachable-demand
handling, capacity-derivative caching, TNTP transit-node semantics, and experimental
optimizer crash recovery remain separate correctness work.

## Development checks

```bash
cmake --preset linux-release -DBUILD_NATIVE_TESTS=ON
cmake --build --preset build-release
ctest --test-dir build/linux-release --output-on-failure

python -m pip install -e '.[test]'
python -m pytest tests/
python examples/two_routes.py
```

The C++ test consumes only public headers and the library target. The Python
boundary tests use synthetic networks and check lifetime ownership, result
snapshots, option translation, native CNDP filtering, and silent defaults.
