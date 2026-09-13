"""Input conversion for the native library (arrays and constraint CSV files)."""

from __future__ import annotations

import os
import numpy as np

from . import _core
from ._core import LinkConstraint, Network
from . import _validation as validate

def load_constraints(path: str | os.PathLike, verbose: bool = False, *,
                     node_index_base: int = 1) -> list[LinkConstraint]:
    """Load capacity constraints with zero-based native endpoints.

    CSV endpoints default to one-based benchmark IDs; use ``node_index_base=0``
    for zero-based files. An optional zero-based ``link_index`` column orders
    rows and must contain every index exactly once. Without it, file row order
    must match the network. Missing investment costs default to one.
    """
    return _core._load_constraints(str(path), verbose, validate.index_base(node_index_base))


def network_from_arrays(name: str, init_node, term_node, capacity, free_flow_time,
                        b, power, demand, *, length=0.0, speed=0.0, toll=0.0,
                        link_type=1, n_nodes: int | None = None,
                        n_zones: int | None = None) -> Network:
    """Construct a :class:`Network` from per-link arrays and a demand matrix.

    Node ids are 0-based, and — following the TNTP convention — the demand
    zones are nodes ``0 .. n_zones-1``. ``demand[i, j]`` is the travel demand
    from zone ``i`` to zone ``j``. Scalar values for the BPR parameters and
    ``length``/``speed``/``toll``/``link_type`` broadcast to every link.
    """
    init = validate.integers(init_node, "init_node")
    term = validate.integers(term_node, "term_node")
    if init.ndim != 1 or term.shape != init.shape:
        raise ValueError("init_node/term_node must be 1-D arrays of equal length")
    n_links = init.shape[0]

    def per_link(values, label, *, positive=False, integer=False):
        values = (validate.integers(values, label) if integer else
                  validate.nonnegative(values, label, positive=positive))
        try:
            return np.ascontiguousarray(np.broadcast_to(values, (n_links,)))
        except ValueError as exc:
            raise ValueError(f"{label} must be scalar or have one entry per link") from exc

    demand = validate.nonnegative(demand, "demand")
    if demand.ndim != 2 or demand.shape[0] != demand.shape[1] or not demand.shape[0]:
        raise ValueError("demand must be a nonempty square (n_zones x n_zones) matrix")
    zones = validate.integer_scalar(
        n_zones if n_zones is not None else demand.shape[0], "n_zones", minimum=1,
    )
    if zones != demand.shape[0]:
        raise ValueError("n_zones must match the demand matrix dimension")
    if n_nodes is not None:
        nodes = validate.integer_scalar(n_nodes, "n_nodes", minimum=zones)
    elif n_links > 0:
        nodes = validate.integer_scalar(
            max(zones, int(max(init.max(), term.max())) + 1), "n_nodes", minimum=zones,
        )
    else:
        nodes = zones

    return _core._build_network_from_arrays(
        str(name), nodes, zones, np.ascontiguousarray(init), np.ascontiguousarray(term),
        per_link(capacity, "capacity", positive=True), per_link(length, "length"),
        per_link(free_flow_time, "free_flow_time"), per_link(b, "b"), per_link(power, "power"),
        per_link(speed, "speed"), per_link(toll, "toll"),
        per_link(link_type, "link_type", integer=True), np.ascontiguousarray(demand),
    )


def constraints_from_arrays(network: Network, lower, upper,
                            investment_cost=1.0) -> list[LinkConstraint]:
    """Build constraints in native link order, broadcasting scalar values.

    Bounds are finite, positive absolute capacities with lower <= upper.
    Investment costs are finite and nonnegative, but are currently metadata:
    the solver uses its uniform ``theta`` multiplier.
    """
    if not isinstance(network, Network):
        raise TypeError("network must be a Network")
    n_links = network.number_of_links

    def per_link(values, label, positive=False):
        values = validate.nonnegative(values, label, positive=positive)
        try:
            return np.broadcast_to(values, (n_links,))
        except ValueError as exc:
            raise ValueError(f"{label} must be scalar or have one entry per link") from exc

    lower = per_link(lower, "lower_bound", positive=True)
    upper = per_link(upper, "upper_bound", positive=True)
    cost = per_link(investment_cost, "investment_cost_param")
    validate.reject(lower > upper, "lower_bound", "be less than or equal to upper_bound")
    return [
        LinkConstraint(int(init), int(term), float(lb), float(ub), float(c))
        for (init, term), lb, ub, c in zip(network.link_nodes(), lower, upper, cost)
    ]
