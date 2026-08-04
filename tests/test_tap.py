import numpy as np
import pytest

import traffic_assignment as ta


def test_siouxfalls_tapas(data_root):
    network = ta.load_network("SiouxFalls", data_root=data_root)
    result = ta.solve_tap(network)

    assert result.approach.lower().startswith("tapas") or "tapas" in result.approach.lower()
    assert result.relative_gap < 1e-10
    assert np.isfinite(result.total_travel_time)
    assert result.flows.shape == (network.number_of_links,)
    assert np.all(result.flows >= 0)
    # TSTT must equal sum of flow * cost by definition.
    assert result.total_travel_time == pytest.approx(
        float(np.sum(result.flows * result.link_costs)), rel=1e-9
    )


def test_routebased_matches_tapas(data_root):
    tapas = ta.solve_tap(ta.load_network("SiouxFalls", data_root=data_root))
    route_based = ta.solve_tap(
        ta.load_network("SiouxFalls", data_root=data_root),
        approach="routebased",
        alpha=1e-10,
        max_iterations=100,
    )

    assert route_based.relative_gap < 1e-6
    assert route_based.total_travel_time == pytest.approx(
        tapas.total_travel_time, rel=1e-4
    )


def test_two_route_equilibrium():
    # Demand of 10 travels from node 0 to node 1 over two routes: the direct
    # link, and 0->2->1. Both cost 2*(1 + x/10) in aggregate, so equilibrium
    # splits demand evenly: 5 per route at route cost 3.
    demand = np.zeros((2, 2))
    demand[0, 1] = 10.0
    network = ta.network_from_arrays(
        "TwoRoutes",
        init_node=[0, 0, 2],
        term_node=[1, 2, 1],
        capacity=10.0,
        free_flow_time=[2.0, 1.0, 1.0],
        b=1.0,
        power=1.0,
        demand=demand,
        n_nodes=3,
    )

    result = ta.solve_tap(network)

    assert result.relative_gap < 1e-12
    np.testing.assert_allclose(result.flows, [5.0, 5.0, 5.0], atol=1e-6)
    np.testing.assert_allclose(result.link_costs, [3.0, 1.5, 1.5], atol=1e-6)


def test_network_accessors(data_root):
    network = ta.load_network("SiouxFalls", data_root=data_root)

    assert network.number_of_links == 76
    assert network.number_of_zones == 24
    nodes = network.link_nodes()
    assert nodes.shape == (76, 2)
    first = network.link(0)
    assert (first.init, first.term) == (int(nodes[0, 0]), int(nodes[0, 1]))

    capacities = network.capacities()
    network.set_capacities(capacities * 2.0)
    np.testing.assert_allclose(network.capacities(), capacities * 2.0)
    with pytest.raises(IndexError):
        network.link(network.number_of_links)
