"""High-level Traffic Assignment (user equilibrium) solving."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._core import Network, TapOptions, TrafficAssignmentApproach, _make_approach, _solve_tap
from .datasets import load_network


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


def make_approach(network: Network, approach: str | TapOptions = "tapas", **options) -> TrafficAssignmentApproach:
    """Create a native solver from a name or typed :class:`TapOptions`.

    Legacy ``alpha=`` is accepted as an alias for ``relative_gap_tolerance=``.
    A typed options object is copied before applying keyword overrides.
    """
    config = TapOptions(approach) if isinstance(approach, TapOptions) else TapOptions()
    if not isinstance(approach, TapOptions):
        config.approach = approach
    if "alpha" in options:
        if "relative_gap_tolerance" in options:
            raise TypeError("pass either alpha or relative_gap_tolerance, not both")
        options["relative_gap_tolerance"] = options.pop("alpha")
    for key, value in options.items():
        if not hasattr(config, key):
            raise TypeError(f"Unknown TAP option '{key}'")
        setattr(config, key, value)
    return _make_approach(network, config)


def solve_tap(network: Network | str,
              approach: str | TapOptions | TrafficAssignmentApproach = "tapas", *,
              data_root=None, reset: bool = False, **approach_options) -> TapResult:
    """Solve user equilibrium and return flows plus convergence metrics.

    ``network`` is a :class:`Network` or a TNTP dataset name. ``approach`` is
    ``"tapas"``, ``"routebased"``, or a pre-built approach object (in which case
    its own network is used and ``approach_options`` must be empty). With
    ``reset=True`` the network is cleared first; the default keeps existing
    solver state. Reuse an approach object to warm-start TAPAS; creating a new
    TAPAS approach starts a fresh assignment.
    """
    if isinstance(approach, TrafficAssignmentApproach):
        if approach_options:
            raise TypeError("approach_options are only valid when approach is given by name")
        approach_obj = approach
        network = approach_obj.network
    else:
        if isinstance(network, str):
            network = load_network(network, data_root=data_root)
        approach_obj = make_approach(network, approach, **approach_options)

    native = _solve_tap(approach_obj, reset=reset)
    return TapResult(
        approach=native.approach,
        network_name=native.network_name,
        relative_gap=native.relative_gap,
        total_travel_time=native.total_travel_time,
        beckmann_objective=native.beckmann_objective,
        flows=native.flows,
        link_costs=native.link_costs,
        solve_seconds=native.solve_seconds,
    )
