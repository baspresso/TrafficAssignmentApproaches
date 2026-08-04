"""High-level Traffic Assignment (user equilibrium) solving."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from ._core import Network, RouteBasedApproach, TapasApproach, TrafficAssignmentApproach
from .datasets import load_network

_TAPAS_NAMES = {"tapas", "task", "tasktapas"}
_ROUTE_BASED_NAMES = {"routebased"}


@dataclass
class TapResult:
    """Outcome of a user-equilibrium solve."""

    approach: str
    network_name: str
    relative_gap: float
    total_travel_time: float
    beckmann_objective: float
    flows: np.ndarray
    link_costs: np.ndarray
    solve_seconds: float


def make_approach(network: Network, approach: str = "tapas", *, alpha: float = 1e-14,
                  max_iterations: int = 200, mu: float = 0.5, v: float = 0.25,
                  shift_method: str = "NewtonStep", route_search_threads: int = 1,
                  full_iteration_count: int = 3, origin_iteration_count: int = 1,
                  ema_alpha: float = 0.7) -> TrafficAssignmentApproach:
    """Create a TAP solver by name; defaults mirror the ``tap_solver`` executable."""
    key = approach.lower().replace("_", "")
    if key in _TAPAS_NAMES:
        return TapasApproach(network, alpha, max_iterations, mu, v)
    if key in _ROUTE_BASED_NAMES:
        return RouteBasedApproach(network, alpha, shift_method, route_search_threads,
                                  max_iterations, full_iteration_count,
                                  origin_iteration_count, ema_alpha)
    raise ValueError(f"Unsupported approach '{approach}'. Supported: Tapas, RouteBased.")


def solve_tap(network: Network | str,
              approach: str | TrafficAssignmentApproach = "tapas", *,
              data_root=None, reset: bool = False, **approach_options) -> TapResult:
    """Solve user equilibrium and return flows plus convergence metrics.

    ``network`` is a :class:`Network` or a TNTP dataset name. ``approach`` is
    ``"tapas"``, ``"routebased"``, or a pre-built approach object (in which case
    its own network is used and ``approach_options`` must be empty). With
    ``reset=True`` the network is cleared first; the default keeps existing
    flows/routes as a warm start.
    """
    if isinstance(approach, TrafficAssignmentApproach):
        if approach_options:
            raise TypeError("approach_options are only valid when approach is given by name")
        approach_obj = approach
        network = approach_obj.network
        if reset:
            approach_obj.reset()
    else:
        if isinstance(network, str):
            network = load_network(network, data_root=data_root)
        if reset:
            network.reset()
        approach_obj = make_approach(network, approach, **approach_options)

    start = time.perf_counter()
    approach_obj.compute_traffic_flows()
    solve_seconds = time.perf_counter() - start

    return TapResult(
        approach=approach_obj.approach_name,
        network_name=network.name,
        relative_gap=network.relative_gap(),
        total_travel_time=network.total_travel_time(),
        beckmann_objective=network.beckmann_objective(),
        flows=network.flows(),
        link_costs=network.link_costs(),
        solve_seconds=solve_seconds,
    )
