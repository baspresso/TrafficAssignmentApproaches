"""Public API ownership and C++/Python boundary checks, with no external data."""

import gc

import numpy as np
import pytest

import traffic_assignment as ta


@pytest.fixture
def network():
    # Fourth link has no capacity-dependent delay and is excluded from design.
    return ta.network_from_arrays(
        "TwoRoutes", init_node=[0, 0, 2, 2], term_node=[1, 2, 1, 3],
        capacity=10, free_flow_time=[2, 1, 1, 1], b=[1, 1, 1, 0], power=1,
        demand=[[0, 10], [0, 0]], n_nodes=4,
    )


def test_typed_options_and_result_snapshot(network):
    options = ta.TapOptions()
    options.approach = "route_based"
    options.relative_gap_tolerance = 1e-10
    result = ta.solve_tap(network, approach=options, relative_gap_tolerance=1e-12)
    assert options.relative_gap_tolerance == pytest.approx(1e-10)
    np.testing.assert_allclose(result.flows, [5, 5, 5, 0], atol=1e-6)
    network.reset()
    np.testing.assert_allclose(network.flows(), 0)
    np.testing.assert_allclose(result.flows, [5, 5, 5, 0], atol=1e-6)


def test_native_ownership_survives_python_gc(network):
    approach = ta.make_approach(network)
    link = network.link(0)
    del network
    gc.collect()
    approach.compute_traffic_flows()
    assert link.flow == pytest.approx(5)
    del approach
    gc.collect()
    assert link.flow == pytest.approx(5)
    link.capacity = 20
    assert link.capacity == 20


def test_native_cndp_filtering_and_silent_defaults(network, tmp_path, monkeypatch, capfd):
    monkeypatch.chdir(tmp_path)
    constraints = ta.constraints_from_arrays(network, 10, 20)
    result = ta.solve_cndp(
        network, [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=3)],
        constraints=constraints, budget=100,
    )
    assert result.upper_bounds[3] == 10
    assert constraints[3].upper_bound == 20
    assert result.objective == pytest.approx(result.total_travel_time + result.budget)
    assert result.budget <= 100 + 1e-6
    assert not list(tmp_path.iterdir())
    captured = capfd.readouterr()
    assert captured.out == captured.err == ""
    saved_flows = result.flows.copy()
    result.network.reset()
    np.testing.assert_array_equal(result.flows, saved_flows)


def test_native_cndp_options_and_network_identity(network):
    constraints = ta.constraints_from_arrays(network, 10, 20)
    pipeline = [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=3)]
    options = ta.CndpOptions()
    options.budget_upper_bound = 100
    options.exclude_insensitive_links = False
    solver = ta.BilevelCND(network, ta.make_approach(network), constraints, pipeline, options)
    assert solver.active_constraint_count == 4
    result = solver.compute_network_design()
    assert result.upper_bounds[3] == 20

    other = ta.network_from_arrays(
        "Other", init_node=[0], term_node=[1], capacity=10,
        free_flow_time=1, b=1, power=1, demand=[[0, 10], [0, 0]],
    )
    with pytest.raises(ValueError, match="supplied network"):
        ta.BilevelCND(network, ta.make_approach(other), constraints, pipeline, options)


def test_option_errors(network):
    with pytest.raises(TypeError, match="either alpha"):
        ta.solve_tap(network, alpha=1e-8, relative_gap_tolerance=1e-8)
    with pytest.raises(TypeError, match="Unknown TAP option"):
        ta.solve_tap(network, unexpected=True)
    with pytest.raises(ValueError, match="Unsupported approach"):
        ta.solve_tap(network, approach="missing")


def test_io_compatibility_imports():
    from traffic_assignment import datasets, io

    assert datasets.network_from_arrays is io.network_from_arrays
    assert datasets.load_constraints is io.load_constraints
    assert datasets.constraints_from_arrays is io.constraints_from_arrays
