"""DataFrame input contracts, link/zone alignment, and solver parity."""

import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

import traffic_assignment as ta


@pytest.fixture
def frames():
    links = pd.DataFrame({
        "link_index": [0, 1, 2],
        "init_node": [0, 0, 2], "term_node": [1, 2, 1],
        "capacity": [10.0, 10.0, 10.0], "free_flow_time": [2.0, 1.0, 1.0],
        "b": [1.0, 1.0, 1.0], "power": [1.0, 1.0, 1.0],
    }, index=["direct", "first", "second"])
    demand = pd.DataFrame([[0.0, 10.0], [0.0, 0.0]], index=[0, 1], columns=[0, 1])
    return links, demand


@pytest.mark.parametrize("base", [0, 1])
@pytest.mark.parametrize("approach", ["tapas", "routebased"])
def test_zone_and_link_alignment_and_ownership(frames, base, approach):
    links, demand = frames
    links[["init_node", "term_node"]] += base
    demand.index += base
    demand.columns += base
    # Nullable numeric dtypes are valid as long as no values are missing.
    links["init_node"] = links.init_node.astype("Int64")
    links["capacity"] = links.capacity.astype("Float64")
    links = links.iloc[[2, 0, 1]]
    demand = demand.iloc[[1, 0], [0, 1]]
    original_links, original_demand = links.copy(deep=True), demand.copy(deep=True)
    network = ta.network_from_dataframes("TwoRoutes", links, demand, node_index_base=base)
    result = ta.solve_tap(network, approach=approach)
    np.testing.assert_array_equal(network.link_nodes(), [[0, 1], [0, 2], [2, 1]])
    np.testing.assert_allclose(result.flows, [5, 5, 5], atol=1e-6)
    assert result.total_travel_time == pytest.approx(30)
    pd.testing.assert_frame_equal(links, original_links)
    pd.testing.assert_frame_equal(demand, original_demand)
    links.loc[:, "capacity"] = 20.0
    demand.iloc[:, :] = 0.0
    np.testing.assert_array_equal(network.capacities(), [10, 10, 10])
    assert ta.solve_tap(network).total_travel_time == pytest.approx(30)


def test_row_order_without_link_index_and_optional_columns(frames):
    links, demand = frames
    links = links.drop(columns="link_index").iloc[[2, 0, 1]].assign(
        length=[4.0, 5.0, 6.0], speed=50.0, toll=2.0, link_type=3, street="Main",
    )
    network = ta.network_from_dataframes("RowOrder", links, demand, n_nodes=4)
    np.testing.assert_array_equal(network.link_nodes(), [[2, 1], [0, 1], [0, 2]])
    assert network.number_of_nodes == 4
    first = network.link(0)
    assert (first.length, first.speed, first.toll, first.type) == (4, 50, 2, 3)


@pytest.mark.parametrize(("column", "value", "error"), [
    ("capacity", np.nan, "finite"),
    ("capacity", np.inf, "finite"),
    ("capacity", 0.0, "strictly positive"),
    ("free_flow_time", -1.0, "nonnegative"),
    ("b", -0.1, "nonnegative"),
    ("power", -1.0, "nonnegative"),
    ("init_node", 0.5, "integer"),
    ("term_node", -1, "integers in"),
    ("init_node", 2**40, "integers in"),
    ("link_index", 1, "exactly once"),
])
def test_invalid_link_values(frames, column, value, error):
    links, demand = frames
    links[column] = links[column].astype(float)
    links.loc["direct", column] = value
    with pytest.raises(ValueError, match=error):
        ta.network_from_dataframes("Invalid", links, demand)


@pytest.mark.parametrize("values", [["10", "10", "10"], [True, True, True], [1j, 1j, 1j]])
def test_reject_non_real_dtypes(frames, values):
    links, demand = frames
    links["capacity"] = values
    with pytest.raises(ValueError, match="capacity.*real numeric dtype"):
        ta.network_from_dataframes("Invalid", links, demand)


def test_nullable_missing_value_reports_row(frames):
    links, demand = frames
    links["capacity"] = links.capacity.astype("Float64")
    links.loc["first", "capacity"] = pd.NA
    with pytest.raises(ValueError, match="links.capacity.*row 'first'"):
        ta.network_from_dataframes("Missing", links, demand)


def test_frame_structure_errors(frames):
    links, demand = frames
    with pytest.raises(TypeError, match="links must be a pandas DataFrame"):
        ta.network_from_dataframes("Invalid", links.to_dict(), demand)
    with pytest.raises(ValueError, match="missing required columns: capacity"):
        ta.network_from_dataframes("Invalid", links.drop(columns="capacity"), demand)
    with pytest.raises(ValueError, match="unique column"):
        ta.network_from_dataframes("Invalid", pd.concat([links, links.capacity], axis=1), demand)
    with pytest.raises(ValueError, match="square"):
        ta.network_from_dataframes("Invalid", links, demand.iloc[:1])
    with pytest.raises(ValueError, match="nonempty"):
        ta.network_from_dataframes("Invalid", links, pd.DataFrame())


@pytest.mark.parametrize("axis", ["index", "columns"])
@pytest.mark.parametrize("labels", [[0, 0], [0, 2], [0, 0.5], ["0", "1"]])
def test_invalid_zone_labels(frames, axis, labels):
    links, demand = frames
    setattr(demand, axis, labels)
    with pytest.raises(ValueError, match="demand"):
        ta.network_from_dataframes("Invalid", links, demand)


@pytest.mark.parametrize("value", [-1.0, np.nan, np.inf])
def test_invalid_demand_values(frames, value):
    links, demand = frames
    demand.loc[0, 1] = value
    with pytest.raises(ValueError, match="demand column 1.*row 0"):
        ta.network_from_dataframes("Invalid", links, demand)


@pytest.mark.parametrize("base", [-1, 2, 1.0, True])
def test_invalid_node_index_base(frames, base):
    with pytest.raises(ValueError, match="node_index_base"):
        ta.network_from_dataframes("Invalid", *frames, node_index_base=base)


@pytest.mark.parametrize("nodes", [1, 2, 3.5, True, 2**40])
def test_invalid_node_count(frames, nodes):
    with pytest.raises(ValueError, match="n_nodes|node IDs below"):
        ta.network_from_dataframes("Invalid", *frames, n_nodes=nodes)


def test_inferred_nodes_include_isolated_zones(frames):
    links, demand = frames
    links = links.iloc[:1].drop(columns="link_index")
    demand = demand.reindex(index=range(4), columns=range(4), fill_value=0.0)
    network = ta.network_from_dataframes("Isolated", links, demand)
    assert network.number_of_nodes == network.number_of_zones == 4


def test_shuffled_constraints_and_parallel_links(frames):
    links, demand = frames
    # Parallel links have distinct bounds despite sharing their endpoints.
    links.loc["second", ["init_node", "term_node"]] = [0, 1]
    network = ta.network_from_dataframes("Parallel", links, demand)
    frame = links[["link_index", "init_node", "term_node"]].assign(
        lower_bound=[10.0, 12.0, 14.0], upper_bound=[20.0, 22.0, 24.0],
        investment_cost_param=[1.0, 2.0, 3.0],
    ).iloc[[2, 0, 1]]
    frame[["init_node", "term_node"]] += 1
    original = frame.copy(deep=True)
    constraints = ta.constraints_from_dataframe(network, frame, node_index_base=1)
    assert [c.lower_bound for c in constraints] == [10, 12, 14]
    assert [c.upper_bound for c in constraints] == [20, 22, 24]
    assert [c.investment_cost_param for c in constraints] == [1, 2, 3]
    assert [(c.init_node, c.term_node) for c in constraints] == [(0, 1), (0, 2), (0, 1)]
    pd.testing.assert_frame_equal(frame, original)


@pytest.fixture
def constraint_inputs(frames):
    links, demand = frames
    network = ta.network_from_dataframes("TwoRoutes", links, demand)
    constraints = links[["link_index", "init_node", "term_node"]].assign(
        lower_bound=10.0, upper_bound=20.0,
    )
    return network, constraints


@pytest.mark.parametrize("indices", [[0, 1], [0, 1, 1], [0, 1, 3], [0, 1, 2, 3]])
def test_constraint_coverage(constraint_inputs, indices):
    network, constraints = constraint_inputs
    constraints = constraints.reindex(range(len(indices))).assign(
        link_index=indices, lower_bound=10, upper_bound=20,
    )
    with pytest.raises(ValueError, match="link_index.*exactly once"):
        ta.constraints_from_dataframe(network, constraints)


@pytest.mark.parametrize(("column", "value", "error"), [
    ("lower_bound", 0.0, "strictly positive"),
    ("lower_bound", 21.0, "less than or equal"),
    ("upper_bound", np.inf, "finite"),
    ("upper_bound", np.nan, "finite"),
    ("init_node", 1, "match the network endpoint"),
    ("investment_cost_param", -1, "nonnegative"),
])
def test_invalid_constraints(constraint_inputs, column, value, error):
    network, constraints = constraint_inputs
    constraints[column] = value
    with pytest.raises(ValueError, match=error):
        ta.constraints_from_dataframe(network, constraints)


def test_constraints_require_key_and_paired_endpoints(constraint_inputs):
    network, constraints = constraint_inputs
    with pytest.raises(ValueError, match="missing required columns: link_index"):
        ta.constraints_from_dataframe(network, constraints.drop(columns="link_index"))
    with pytest.raises(ValueError, match="together"):
        ta.constraints_from_dataframe(network, constraints.drop(columns="term_node"))


def test_cndp_from_frames_without_files(frames, constraint_inputs, tmp_path, monkeypatch):
    network, frame = constraint_inputs
    original = frame.copy(deep=True)
    # Fixed links and omitted endpoint columns are supported.
    frame.loc["second", "upper_bound"] = 10.0
    constraints = ta.constraints_from_dataframe(
        network, frame.drop(columns=["init_node", "term_node"]).iloc[::-1],
    )
    monkeypatch.chdir(tmp_path)
    result = ta.solve_cndp(
        network, [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=5)],
        constraints=constraints, budget=100,
    )
    assert result.objective == pytest.approx(result.total_travel_time + result.budget)
    assert result.budget <= 100 + 1e-6
    assert np.all(result.capacities >= result.lower_bounds - 1e-9)
    assert np.all(result.capacities <= result.upper_bounds + 1e-9)
    assert result.capacities[2] == pytest.approx(10)
    assert not list(tmp_path.iterdir())
    # The solver mutates the native network, not the source link DataFrame.
    np.testing.assert_array_equal(frames[0].capacity, [10, 10, 10])
    pd.testing.assert_series_equal(frame.lower_bound, original.lower_bound)


def test_siouxfalls_csv_parity(data_root):
    folder = data_root / "SiouxFalls"
    if not all((folder / f"SiouxFalls_{part}.csv").is_file() for part in ("net", "trips", "constraints")):
        pytest.skip("Sioux Falls CSV files not available")
    links = pd.read_csv(folder / "SiouxFalls_net.csv", comment="#")
    links.insert(0, "link_index", np.arange(len(links)))
    demand = pd.read_csv(folder / "SiouxFalls_trips.csv", comment="#", header=None)
    demand.index = demand.columns = range(1, len(demand) + 1)

    def from_frames():
        return ta.network_from_dataframes(
            "SiouxFalls", links.iloc[::-1], demand.iloc[::-1, ::-1], node_index_base=1,
        )

    csv_tap = ta.solve_tap("SiouxFalls", data_root=data_root)
    frame_tap = ta.solve_tap(from_frames())
    np.testing.assert_allclose(frame_tap.flows, csv_tap.flows, rtol=1e-10, atol=1e-8)
    np.testing.assert_allclose(frame_tap.link_costs, csv_tap.link_costs, rtol=1e-10)
    assert frame_tap.total_travel_time == pytest.approx(csv_tap.total_travel_time, rel=1e-10)

    bounds = pd.read_csv(folder / "SiouxFalls_constraints.csv")
    bounds.insert(0, "link_index", np.arange(len(bounds)))
    network = from_frames()
    constraints = ta.constraints_from_dataframe(network, bounds.iloc[::-1], node_index_base=1)
    pipeline = [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=3)]
    csv_cndp = ta.solve_cndp("SiouxFalls", pipeline, data_root=data_root, budget=10000)
    frame_cndp = ta.solve_cndp(network, pipeline, constraints=constraints, budget=10000)
    np.testing.assert_allclose(frame_cndp.capacities, csv_cndp.capacities, rtol=1e-10, atol=1e-8)
    np.testing.assert_allclose(frame_cndp.flows, csv_cndp.flows, rtol=1e-10, atol=1e-8)
    assert frame_cndp.objective == pytest.approx(csv_cndp.objective, rel=1e-10)


def test_pandas_is_optional():
    script = """
import builtins
original_import = builtins.__import__
def without_pandas(name, *args, **kwargs):
    if name == 'pandas' or name.startswith('pandas.'):
        raise ImportError('pandas disabled for test')
    return original_import(name, *args, **kwargs)
builtins.__import__ = without_pandas
import traffic_assignment as ta
network = ta.network_from_arrays(
    'NoPandas', init_node=[0], term_node=[1], capacity=10,
    free_flow_time=1, b=1, power=1, demand=[[0, 10], [0, 0]],
)
assert ta.solve_tap(network).flows[0] == 10
try:
    ta.network_from_dataframes('MissingPandas', None, None)
except ImportError as exc:
    assert 'traffic-assignment[dataframes]' in str(exc)
else:
    raise AssertionError('expected an actionable pandas import error')
"""
    subprocess.run([sys.executable, "-c", script], check=True, capture_output=True, text=True)
