# Providing solver inputs

DataFrames, NumPy arrays, and preprocessed CSV datasets all build the same native
`Network`. Both TAP and CNDP consume that network; CNDP also needs capacity
constraints. pandas is optional and is only imported by the DataFrame adapters.

Prepare inputs explicitly before solving:

```python
import traffic_assignment as ta

network = ta.network_from_dataframes("MyNetwork", links_df, demand_df)
constraints = ta.constraints_from_dataframe(network, constraints_df)
result = ta.solve_cndp(
    network, [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=100)],
    constraints=constraints, budget=10000,
)
```

See the [DataFrame guide](dataframes.md) for the complete table schemas and
examples. For generated networks, arrays provide the same numerical contract:

```python
network = ta.network_from_arrays(
    "TwoRoutes", init_node=[0, 0, 2], term_node=[1, 2, 1],
    capacity=10, free_flow_time=[2, 1, 1], b=1, power=1,
    demand=[[0, 10], [0, 0]],
)
assignment = ta.solve_tap(network)
constraints = ta.constraints_from_arrays(network, lower=10, upper=20)
```

## Shared numerical contract

- IDs are integers within the native 32-bit range. Python checks IDs before
  casting, so fractional values are rejected rather than truncated.
- Capacities and constraint bounds are finite and strictly positive.
- Demand, free-flow times, BPR parameters, length, speed, toll, and investment
  cost metadata are finite and nonnegative. Zero demand and constant-cost links
  are supported.
- Demand is a nonempty square matrix. Zone IDs correspond to the first
  `n_zones` nodes. Inferred node counts include isolated demand zones.
- Python numerical inputs must have real numeric dtypes. Boolean, complex,
  missing, and string values are rejected. CSV loaders parse complete numeric
  fields and reject malformed numbers.
- Array endpoint inputs are one-dimensional and equally sized. Per-link
  numerical fields and constraint bounds accept either scalars or one entry
  per link. Constraints require `lower_bound <= upper_bound`.

Python shares validation helpers across array and DataFrame adapters. Native
construction independently enforces the domain invariants for direct C++ callers
and CSV loading, before building graph adjacency. CNDP validates every constraint
and its endpoints before filtering insensitive links or modifying capacities.
Public capacity setters also validate all supplied values before applying them.

## File loading and indexing

```python
network = ta.load_network("SiouxFalls", data_root="data/TransportationNetworks")
constraints = ta.load_constraints(
    ta.default_constraints_path("SiouxFalls", data_root="data/TransportationNetworks"),
    node_index_base=1,
)
result = ta.solve_cndp(
    network, [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=100)],
    constraints=constraints,
)
```

`load_network` reads the project's preprocessed `*_net.csv` and `*_trips.csv`
files, including their metadata lines. It does not read raw `.tntp` files.
Network columns are read by name: required and optional fields follow the
DataFrame link schema. Extra metadata columns are ignored. In particular,
Sydney's `critical_speed` and `lanes` fields are not interpreted as `toll` and
`link_type`; the missing optional fields receive defaults of zero and one.
Network endpoints in these files are one-based; demand rows and columns follow
zone order. Declared node/zone/link counts and demand dimensions are checked.

Constraint CSVs require `init_node`, `term_node`, `lower_bound`, and `upper_bound`.
`investment_cost_param` is optional and defaults to one. Comma, semicolon, tab,
and pipe delimiters are supported, along with blank lines and `#` comments.
Malformed data rows raise an error containing the file and line number; they
are never silently discarded.

| Input | Endpoint convention | Link ordering |
| --- | --- | --- |
| Arrays and direct native objects | Zero-based | Array/list position |
| DataFrames | `node_index_base=0` by default | `link_index`, or link row position |
| Network CSV | One-based | Optional `link_index`, otherwise row position |
| Constraint CSV | `node_index_base=1` by default | Optional `link_index`, otherwise row position |

All native endpoints and link indices are zero-based. For a constraint file
whose endpoints already use zero-based IDs, call
`load_constraints(path, node_index_base=0)`. The C++ equivalent is
`LoadConstraints(path, verbose, node_index_base)`.

When present in a CSV, `link_index` must contain every index from zero to the
number of rows minus one exactly once. Rows are sorted by that key. Legacy
constraint CSVs without the key must stay in network link order; CNDP checks
their endpoints against the indexed links. Endpoint checks cannot distinguish
reordered parallel links sharing the same endpoints, so use `link_index` when
reordering those records. Result arrays always follow native link order.

The solver still uses `theta * sum(capacity - lower_bound)` for investment;
per-link investment cost fields remain metadata. Bounds are absolute capacities,
and CNDP initializes capacities to their lower bounds.

## Ownership and compatibility

Input arrays and DataFrames are copied into native storage. Modifying them later
does not update the network. Solving modifies network state, and CNDP modifies
capacities. Build a fresh network for independent scenarios. `Network.reset()`
clears flows and routes; it does not restore initial capacities. Reuse an approach
object to warm-start its solver state. A supplied approach must belong to the
supplied network; TAP now rejects mismatches, as CNDP already does.

Existing top-level functions and compatibility imports from `datasets` remain
available. Dataset names, constraint paths, and automatic constraint-file lookup
in solver calls also remain supported. Explicit loaders are preferred because
they make filesystem access and the chosen constraints visible at the call site.
Omitting CNDP constraints still looks up the conventional file using the network
name; it does not generate bounds for an in-memory network.

Previously accepted invalid numeric inputs now fail validation. Loaded
constraints now expose zero-based endpoints, matching DataFrame and array
constraints; code inspecting those endpoint fields should use native indexing.
