# DataFrame inputs

Install `pip install -e '.[dataframes]'` from the repository root. pandas is an
optional dependency: array inputs, CSV loading, and the native solvers remain
available without it. Install `'.[examples]'` to also get JupyterLab, a Python
kernel, and plotting dependencies.

The two adapters are exported at the package root:

```python
import pandas as pd
import traffic_assignment as ta

network = ta.network_from_dataframes(
    "MyNetwork", links_df, demand_df, node_index_base=1,
)
assignment = ta.solve_tap(network)

constraints = ta.constraints_from_dataframe(
    network, constraints_df, node_index_base=1,
)
design = ta.solve_cndp(
    network,
    [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=100)],
    constraints=constraints, theta=5.0, budget=10000.0,
)
```

The adapters validate and copy values into the existing native types. They do
not modify the input DataFrames. Later DataFrame edits do not update a built
network; construct another network for another scenario. Solving updates the
native network's state, and CNDP changes its capacities. Returned result arrays
are independent snapshots in native link order.

## Links

Each row describes one directed link. Column order and pandas row labels do not
matter. Extra columns may contain application metadata and are ignored.

| Column | Requirement |
| --- | --- |
| `init_node`, `term_node` | Required integer node IDs |
| `capacity` | Required, strictly positive |
| `free_flow_time`, `b`, `power` | Required, nonnegative BPR parameters |
| `length`, `speed`, `toll` | Optional nonnegative values, default zero |
| `link_type` | Optional nonnegative integer, default one |
| `link_index` | Optional unique native link index, described below |

Numeric columns must have real numeric dtypes and finite values; missing values,
booleans, complex values, and numeric strings are rejected. pandas nullable
numeric dtypes are accepted when all values are present. Errors identify the
field and, for invalid numeric values, the original row label.

`node_index_base=0` is the default. Use `node_index_base=1` for the conventional
TNTP labels. Native node IDs are always zero-based, and demand zones must occupy
the first `n_zones` nodes. Arbitrary external node labels require preprocessing.
By default the node count includes all endpoints and demand zones; `n_nodes=`
can explicitly include isolated nodes. IDs must fit the native 32-bit integer
range.

If `link_index` is omitted, input row position becomes the native link index.
If supplied, it must contain every integer from zero to `len(links_df) - 1`
exactly once. The adapter sorts links by this column before construction. Assign
it once, before shuffling rows, and retain it when joining results:

```python
links_df.insert(0, "link_index", range(len(links_df)))
# ... construct and solve the network ...
flows_df = pd.DataFrame({
    "link_index": range(network.number_of_links),
    "flow": assignment.flows,
    "travel_time": assignment.link_costs,
})
links_with_flows = links_df.merge(flows_df, on="link_index", validate="one_to_one")
```

Link indices are always zero-based, regardless of node indexing. Distinct link
indices keep parallel links distinguishable. Filtering out links changes the
network: assign a new contiguous set of link indices and rebuild its constraints.

## Demand

Demand is a nonempty square DataFrame. Rows are origins and columns are
destinations. Both axes must contain every integer zone ID from
`node_index_base` to `node_index_base + n_zones - 1` exactly once. Their orders
may differ; the adapter aligns each axis by label. Values must be finite,
nonnegative numbers, and zero demand is valid.

For a two-zone network using one-based node IDs:

```python
import pandas as pd

demand_df = pd.DataFrame(
    [[0.0, 10.0], [0.0, 0.0]],
    index=pd.Index([1, 2], name="origin"),
    columns=pd.Index([1, 2], name="destination"),
)
```

The initial adapter supports dense square matrices. An OD table with
`origin, destination, demand` columns must be pivoted and reindexed by the
caller; duplicate OD records need an explicit aggregation policy. Native
construction still uses a dense matrix, so memory scales with the square of the
number of zones.

## CNDP constraints

`constraints_from_dataframe(network, constraints_df)` requires `link_index`,
`lower_bound`, and `upper_bound`. It requires exactly one row for every network
link and aligns rows by `link_index`. Missing, duplicate, and unknown link indices
are errors. A pandas index alone does not supply link identity.

Bounds are finite, strictly positive **absolute capacities**, with
`lower_bound <= upper_bound`. Equal bounds fix a link's capacity. CNDP starts at
the lower bounds. For an expansion-only scenario, use the initial capacities as
lower bounds:

```python
constraints_df = links_df[["link_index", "init_node", "term_node"]].copy()
constraints_df["lower_bound"] = links_df["capacity"]
constraints_df["upper_bound"] = 1.2 * links_df["capacity"]
```

Endpoint columns are optional but must be supplied together. When present, the
adapter checks that they match the native link identified by `link_index` after
applying `node_index_base`. This check catches constraints from another network
or a stale link mapping. Zero capacity is not accepted as a road-closure model.

The implemented upper-level objective and budget are:

```text
investment = theta * sum(final_capacity - lower_bound)
objective = final_total_travel_time + investment
investment <= budget
```

An optional finite, nonnegative `investment_cost_param` column is copied into
the native constraints (default one) for compatibility. The current solver does
not use that per-link field in these formulas. Use `theta` for the uniform cost
multiplier. Supporting link-specific investment costs requires a separate solver
change.

The solver may fix insensitive links by setting their effective upper bound to
their lower bound. Inspect `design.lower_bounds` and `design.upper_bounds` for
the effective bounds. Supplying prepared constraints explicitly also avoids the
high-level solver's default lookup of a dataset constraints CSV.

## Examples and validation

[Sioux Falls TAP](../examples/siouxfalls_tap.ipynb) and
[Sioux Falls CNDP](../examples/siouxfalls_cndp.ipynb) each contain the full network
and demand values and construct their own DataFrames. Run all cells in order.
They require neither the other notebook nor local CSV files. Outputs include
metrics, per-link tables, plots, and numerical consistency checks. The CNDP
iteration limit demonstrates the workflow and does not certify a global optimum.

`tests/test_dataframes.py` covers analytic TAP results, zone/link reordering,
parallel-link constraint alignment, invalid inputs, ownership, pandas being
optional, and TAP/CNDP parity with the Sioux Falls CSV path when data is present.
