"""Input parity and validation, using small networks without downloaded data."""

import numpy as np
import pandas as pd
import pytest

import traffic_assignment as ta


@pytest.fixture
def inputs():
    links = pd.DataFrame({
        "init_node": [0, 0, 2], "term_node": [1, 2, 1], "capacity": [10.0] * 3,
        "length": [0.0] * 3, "free_flow_time": [2.0, 1.0, 1.0],
        "b": [1.0] * 3, "power": [1.0] * 3, "speed": [0.0] * 3,
        "toll": [0.0] * 3, "link_type": [1] * 3,
    })
    demand = pd.DataFrame([[0.0, 10.0], [0.0, 0.0]])
    return links, demand


def write_dataset(root, links, demand):
    folder = root / "Tiny"
    folder.mkdir(exist_ok=True)
    net_path = folder / "Tiny_net.csv"
    trips_path = folder / "Tiny_trips.csv"
    links = links.copy()
    links[["init_node", "term_node"]] += 1
    net_path.write_text(
        f"# NUMBER OF ZONES:2,NUMBER OF NODES:3,NUMBER OF LINKS:{len(links)}\n"
        + links.to_csv(index=False)
    )
    trips_path.write_text("# NUMBER OF ZONES:2\n" + demand.to_csv(index=False, header=False))
    return net_path, trips_path


def build(source, inputs, tmp_path):
    links, demand = inputs
    if source == "arrays":
        return ta.network_from_arrays("Tiny", **{c: links[c].to_numpy() for c in links},
                                      demand=demand.to_numpy(), n_nodes=3)
    if source == "dataframes":
        return ta.network_from_dataframes("Tiny", links, demand, n_nodes=3)
    write_dataset(tmp_path, links, demand)
    return ta.load_network("Tiny", data_root=tmp_path)


@pytest.mark.parametrize("source", ["arrays", "dataframes", "csv"])
@pytest.mark.parametrize(("column", "value"), [
    ("capacity", 0), ("capacity", np.inf), ("free_flow_time", -1),
    ("b", -1), ("power", np.nan), ("length", -1), ("speed", -1), ("toll", -1),
    ("init_node", 0.5), ("term_node", 2**40), ("link_type", 1.5),
])
def test_invalid_link_values_across_sources(inputs, tmp_path, source, column, value):
    links, _ = inputs
    links[column] = float(value)
    with pytest.raises((ValueError, RuntimeError), match=column):
        build(source, inputs, tmp_path)


@pytest.mark.parametrize("source", ["arrays", "dataframes", "csv"])
@pytest.mark.parametrize("value", [-1.0, np.nan, np.inf])
def test_invalid_demand_across_sources(inputs, tmp_path, source, value):
    inputs[1].iloc[0, 1] = value
    with pytest.raises((ValueError, RuntimeError), match="demand"):
        build(source, inputs, tmp_path)


@pytest.mark.parametrize("source", ["arrays", "dataframes", "csv"])
def test_tap_cndp_parity_without_datasets(inputs, tmp_path, source):
    network = build(source, inputs, tmp_path)
    result = ta.solve_tap(network)
    np.testing.assert_allclose(result.flows, [5, 5, 5], atol=1e-6)
    assert result.total_travel_time == pytest.approx(30)
    if source == "arrays":
        constraints = ta.constraints_from_arrays(network, 10, 20)
    else:
        frame = inputs[0][["init_node", "term_node"]].assign(
            link_index=range(3), lower_bound=10, upper_bound=20,
        ).iloc[::-1]
        if source == "dataframes":
            constraints = ta.constraints_from_dataframe(network, frame)
        else:
            path = tmp_path / "constraints.csv"
            frame.to_csv(path, index=False)
            constraints = ta.load_constraints(path, node_index_base=0)
    design = ta.solve_cndp(network, [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=3)],
                           constraints=constraints, budget=100)
    reference = build("arrays", inputs, tmp_path)
    expected = ta.solve_cndp(reference, [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=3)],
                             constraints=ta.constraints_from_arrays(reference, 10, 20), budget=100)
    np.testing.assert_allclose(design.capacities, expected.capacities, rtol=1e-10)
    np.testing.assert_allclose(design.flows, expected.flows, rtol=1e-10)
    assert design.objective == pytest.approx(expected.objective, rel=1e-10)


@pytest.mark.parametrize("value", [True, "10", 1j, np.nan, np.inf, -1.0])
def test_array_constraint_validation(inputs, tmp_path, value):
    network = build("arrays", inputs, tmp_path)
    with pytest.raises(ValueError):
        ta.constraints_from_arrays(network, value, 20)


@pytest.mark.parametrize("value", [True, "1", 1j, 0.5, 2**40])
def test_array_ids_checked_before_cast(value):
    with pytest.raises(ValueError, match="init_node"):
        ta.network_from_arrays("Invalid", [value], [1], 10, 1, 1, 1, [[0, 10], [0, 0]])


@pytest.mark.parametrize("field", ["n_nodes", "n_zones"])
@pytest.mark.parametrize("value", [True, 2.5, "2", 0, 2**40])
def test_array_counts_checked_before_cast(field, value):
    with pytest.raises(ValueError, match=field):
        ta.network_from_arrays("Invalid", [0], [1], 10, 1, 1, 1,
                               [[0, 10], [0, 0]], **{field: value})


def test_array_shapes_broadcasting_and_isolated_zones():
    network = ta.network_from_arrays("Isolated", [0], [1], 10, 1, 1, 1, np.zeros((4, 4)))
    assert network.number_of_nodes == 4
    np.testing.assert_array_equal(network.capacities(), [10])
    with pytest.raises(ValueError, match="capacity"):
        ta.network_from_arrays("Invalid", [0], [1], [10, 20], 1, 1, 1, np.zeros((2, 2)))
    with pytest.raises(ValueError, match="nonempty square"):
        ta.network_from_arrays("Invalid", [], [], 10, 1, 1, 1, np.zeros((0, 0)))


@pytest.mark.parametrize(("field", "value"), [
    ("lower_bound", np.nan), ("upper_bound", 0), ("upper_bound", 9),
    ("investment_cost_param", -1), ("init_node", 1),
])
def test_native_constraints_validated_before_mutation(inputs, tmp_path, field, value):
    network = build("arrays", inputs, tmp_path)
    network.set_capacities([15, 15, 15])
    constraints = ta.constraints_from_arrays(network, 10, 20)
    setattr(constraints[0], field, value)
    with pytest.raises(ValueError):
        ta.BilevelCND(network, ta.make_approach(network), constraints,
                      [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=1)], ta.CndpOptions())
    np.testing.assert_array_equal(network.capacities(), [15, 15, 15])


def test_tap_network_identity(inputs, tmp_path):
    first = build("arrays", inputs, tmp_path)
    other = build("arrays", inputs, tmp_path)
    approach = ta.make_approach(first)
    with pytest.raises(ValueError, match="supplied network"):
        ta.solve_tap(other, approach=approach)
    # A matching legacy name reuses the approach without dataset discovery.
    assert ta.solve_tap("Tiny", approach=approach).total_travel_time == pytest.approx(30)


@pytest.mark.parametrize("base", [0, 1])
def test_constraint_csv_parallel_links_and_bases(tmp_path, base):
    path = tmp_path / "constraints.csv"
    path.write_text(
        "# Explicit indices distinguish parallel links.\n"
        "link_index;init_node;term_node;lower_bound;upper_bound\n"
        f"1;{base};{base+1};12;22\n0;{base};{base+1};10;20\n"
    )
    constraints = ta.load_constraints(path, node_index_base=base)
    assert [(c.init_node, c.term_node, c.lower_bound, c.investment_cost_param)
            for c in constraints] == [(0, 1, 10, 1), (0, 1, 12, 1)]


@pytest.mark.parametrize("rows", [
    "0,1,2,10,20\n0,1,3,10,20\n",  # duplicate index
    "0,1,2,10,20\n2,1,3,10,20\n",  # missing index
    "0,1.5,2,10,20\n",              # fractional endpoint
    "0,1,2,10,20junk\n",            # partial number
    "0,1,2,10\n",                   # missing field
    "0,1,2,10,nan\n",               # nonfinite bound
    "0,1,2,10,9\n",                 # reversed bounds
    "0,1,2,0,20\n",                 # nonpositive bound
])
def test_constraint_csv_rejects_bad_rows(tmp_path, rows):
    path = tmp_path / "constraints.csv"
    path.write_text("link_index,init_node,term_node,lower_bound,upper_bound\n" + rows)
    with pytest.raises(RuntimeError, match="line"):
        ta.load_constraints(path)


@pytest.mark.parametrize("header", [
    "init_node,term_node,lower_bound,lower_bound", "init_node,term_node,lower_bound",
])
def test_constraint_csv_header_errors(tmp_path, header):
    path = tmp_path / "constraints.csv"
    path.write_text(header + "\n1,2,10,20\n")
    with pytest.raises(RuntimeError, match="column"):
        ta.load_constraints(path)


@pytest.mark.parametrize("file_kind", ["network", "demand"])
def test_csv_extra_and_missing_rows(inputs, tmp_path, file_kind):
    paths = write_dataset(tmp_path, *inputs)
    path = paths[0 if file_kind == "network" else 1]
    original = path.read_text()
    path.write_text(original + original.splitlines()[-1] + "\n")
    with pytest.raises(RuntimeError, match="row count|one row per zone"):
        ta.load_network("Tiny", data_root=tmp_path)
    path.write_text("\n".join(original.splitlines()[:-1]) + "\n")
    with pytest.raises(RuntimeError, match="row count|one row per zone"):
        ta.load_network("Tiny", data_root=tmp_path)


def test_legacy_csv_order_is_checked_at_solve(inputs, tmp_path):
    network = build("arrays", inputs, tmp_path)
    path = tmp_path / "constraints.csv"
    path.write_text("init_node,term_node,lower_bound,upper_bound\n"
                    "3,2,10,20\n1,3,10,20\n1,2,10,20\n")
    constraints = ta.load_constraints(path)
    with pytest.raises(ValueError, match="endpoints"):
        ta.solve_cndp(network, [ta.step("nlopt", max_iterations=1)], constraints=constraints)


def test_csv_network_uses_column_names_and_optional_defaults(inputs, tmp_path):
    links, demand = inputs
    # Sydney has critical_speed and lanes rather than toll and link_type.
    links = links.drop(columns=["toll", "link_type"]).assign(
        critical_speed=19.2, lanes=2, link_index=range(3),
    ).iloc[::-1, ::-1]
    write_dataset(tmp_path, links, demand)
    network = ta.load_network("Tiny", data_root=tmp_path)
    np.testing.assert_array_equal(network.link_nodes(), [[0, 1], [0, 2], [2, 1]])
    assert network.link(0).toll == 0
    assert network.link(0).type == 1
    assert ta.solve_tap(network).total_travel_time == pytest.approx(30)


@pytest.mark.parametrize("indices", [[0, 0, 2], [0, 1, 3]])
def test_csv_network_rejects_bad_indices(inputs, tmp_path, indices):
    links, demand = inputs
    write_dataset(tmp_path, links.assign(link_index=indices), demand)
    with pytest.raises(RuntimeError, match="link_index.*exactly once"):
        ta.load_network("Tiny", data_root=tmp_path)
