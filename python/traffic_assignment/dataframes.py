"""Optional pandas adapters for the native network and CNDP constraints."""

from __future__ import annotations

from numbers import Integral
from typing import TYPE_CHECKING

import numpy as np

from ._core import LinkConstraint, Network
from .io import network_from_arrays

if TYPE_CHECKING:
    import pandas as pd


def _pandas():
    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            "DataFrame inputs require pandas; install 'traffic-assignment[dataframes]'"
        ) from exc
    return pd


def _check_frame(frame, name: str, required=()):
    if not isinstance(frame, _pandas().DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame")
    if not frame.columns.is_unique:
        raise ValueError(f"{name} must have unique column labels")
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {', '.join(missing)}")


def _reject(series, invalid, label: str, requirement: str):
    positions = np.flatnonzero(invalid)
    if positions.size:
        row = series.index[positions[0]]
        if isinstance(row, np.generic):
            row = row.item()
        raise ValueError(f"{label} must {requirement}; invalid value at row {row!r}")


def _numeric(series, label: str) -> np.ndarray:
    types = _pandas().api.types
    if (not types.is_numeric_dtype(series.dtype)
            or types.is_bool_dtype(series.dtype)
            or types.is_complex_dtype(series.dtype)):
        raise ValueError(f"{label} must have a real numeric dtype")
    values = series.to_numpy(dtype=np.float64, na_value=np.nan)
    _reject(series, ~np.isfinite(values), label, "contain only finite, non-missing values")
    return values


def _integers(series, label: str, minimum: int = 0) -> np.ndarray:
    values = _numeric(series, label)
    _reject(series, values != np.floor(values), label, "contain integer values")
    _reject(series, (values < minimum) | (values > np.iinfo(np.int32).max),
            label, f"contain integers in [{minimum}, {np.iinfo(np.int32).max}]")
    return values.astype(np.int64)


def _nonnegative(series, label: str, *, positive: bool = False) -> np.ndarray:
    values = _numeric(series, label)
    _reject(series, values <= 0 if positive else values < 0, label,
            "be strictly positive" if positive else "be nonnegative")
    return values


def _index_base(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value not in (0, 1):
        raise ValueError("node_index_base must be 0 or 1")
    return int(value)


def _ordered_by_link_index(frame, n_links: int, name: str):
    indices = _integers(frame["link_index"], f"{name}.link_index")
    if len(indices) != n_links or not np.array_equal(np.sort(indices), np.arange(n_links)):
        raise ValueError(f"{name}.link_index must contain every index in [0, {n_links}) exactly once")
    return frame.iloc[np.argsort(indices)]


def network_from_dataframes(name: str, links: pd.DataFrame, demand: pd.DataFrame, *,
                            node_index_base: int = 0,
                            n_nodes: int | None = None) -> Network:
    """Build a native network from links and a labelled square OD matrix.

    ``links`` requires ``init_node``, ``term_node``, ``capacity``,
    ``free_flow_time``, ``b``, and ``power``. Optional ``length``, ``speed``,
    and ``toll`` default to zero; ``link_type`` defaults to one. Other columns
    are metadata and are ignored. The pandas row index is ignored. If a
    ``link_index`` column is present, it must be a permutation of
    ``0 .. len(links)-1`` and links are built in that order; otherwise native
    link indices follow input row order. Result arrays use native link order.

    ``node_index_base`` is zero by default; use one for TNTP node IDs. Demand
    row and column labels must each contain every zone ID from that base to
    ``base + n_zones - 1`` exactly once. The two axes are aligned by label,
    independently of their input order. Zones correspond to the first nodes
    of the network. ``n_nodes`` can include isolated nodes; by default it is
    inferred from endpoints and zones.

    Numeric columns must have real numeric dtypes and finite values. Node IDs
    and link types must be integers, capacities positive, and demand and BPR
    parameters nonnegative. The demand matrix must be nonempty. Zero demand,
    zero free-flow times, and constant-cost links are allowed. Input frames
    are never modified; later edits to them do not change the native network.
    Requires the optional ``dataframes`` extra.
    """
    base = _index_base(node_index_base)
    _check_frame(links, "links", (
        "init_node", "term_node", "capacity", "free_flow_time", "b", "power",
    ))
    _check_frame(demand, "demand")
    if demand.empty or demand.shape[0] != demand.shape[1]:
        raise ValueError("demand must be a nonempty square DataFrame")
    zones = len(demand)
    expected = np.arange(base, base + zones)
    pd = _pandas()
    axes = []
    for axis, labels in (("index", demand.index), ("columns", demand.columns)):
        ids = _integers(pd.Series(labels), f"demand.{axis}", minimum=base)
        if not np.array_equal(np.sort(ids), expected):
            raise ValueError(
                f"demand.{axis} must contain each zone ID in "
                f"[{base}, {base + zones}) exactly once"
            )
        axes.append(np.argsort(ids))
    aligned_demand = demand.iloc[axes[0], axes[1]]
    trips = np.column_stack([
        _nonnegative(aligned_demand[column], f"demand column {column!r}")
        for column in aligned_demand.columns
    ])

    if "link_index" in links.columns:
        links = _ordered_by_link_index(links, len(links), "links")
    init = _integers(links["init_node"], "links.init_node", minimum=base) - base
    term = _integers(links["term_node"], "links.term_node", minimum=base) - base
    if n_nodes is None:
        n_nodes = max(zones, int(max(init.max(), term.max())) + 1) if len(links) else zones
    if (isinstance(n_nodes, bool) or not isinstance(n_nodes, Integral)
            or not zones <= n_nodes <= np.iinfo(np.int32).max):
        raise ValueError("n_nodes must be an integer with n_zones <= n_nodes <= 2147483647")
    for column, values in (("init_node", init), ("term_node", term)):
        _reject(links[column], values >= n_nodes, f"links.{column}",
                f"contain node IDs below {n_nodes + base}")

    fields = {
        column: _nonnegative(links[column], f"links.{column}", positive=column == "capacity")
        for column in ("capacity", "free_flow_time", "b", "power")
    }
    for column in ("length", "speed", "toll"):
        if column in links.columns:
            fields[column] = _nonnegative(links[column], f"links.{column}")
    if "link_type" in links.columns:
        fields["link_type"] = _integers(links["link_type"], "links.link_type")
    return network_from_arrays(name, init, term, demand=trips, n_nodes=int(n_nodes), **fields)


def constraints_from_dataframe(network: Network, constraints: pd.DataFrame, *,
                               node_index_base: int = 0) -> list[LinkConstraint]:
    """Build CNDP constraints, aligned to native links by ``link_index``.

    Required columns are ``link_index``, ``lower_bound``, and ``upper_bound``.
    Every native link index must occur exactly once, including fixed links
    (equal bounds). Bounds are finite, strictly positive absolute capacities,
    with ``lower_bound <= upper_bound``. CNDP starts at the lower bounds.

    Optional ``init_node`` and ``term_node`` must be supplied together and
    match the indexed link, using ``node_index_base`` (zero or one). Link
    indices are always zero-based, even when node IDs are one-based. Parallel
    links are distinguished by link index. Row order and pandas row labels
    do not affect alignment. Inputs are never modified.

    Optional nonnegative ``investment_cost_param`` defaults to one and is
    stored for compatibility. The current solver does not use this per-link
    field: its investment cost is ``theta * sum(capacity - lower_bound)``.
    Requires the optional ``dataframes`` extra.
    """
    base = _index_base(node_index_base)
    if not isinstance(network, Network):
        raise TypeError("network must be a Network")
    _check_frame(constraints, "constraints", ("link_index", "lower_bound", "upper_bound"))
    constraints = _ordered_by_link_index(constraints, network.number_of_links, "constraints")
    nodes = network.link_nodes()
    endpoints = [column in constraints.columns for column in ("init_node", "term_node")]
    if any(endpoints) and not all(endpoints):
        raise ValueError("constraints must supply init_node and term_node together")
    if all(endpoints):
        for i, column in enumerate(("init_node", "term_node")):
            ids = _integers(constraints[column], f"constraints.{column}", minimum=base) - base
            _reject(constraints[column], ids != nodes[:, i], f"constraints.{column}",
                    "match the network endpoint for link_index")
    lower = _nonnegative(constraints["lower_bound"], "constraints.lower_bound", positive=True)
    upper = _nonnegative(constraints["upper_bound"], "constraints.upper_bound", positive=True)
    _reject(constraints["lower_bound"], lower > upper, "constraints.lower_bound",
            "be less than or equal to upper_bound")
    costs = (
        _nonnegative(constraints["investment_cost_param"], "constraints.investment_cost_param")
        if "investment_cost_param" in constraints.columns else np.ones(len(constraints))
    )
    return [
        LinkConstraint(int(init), int(term), float(lb), float(ub), float(cost))
        for (init, term), lb, ub, cost in zip(nodes, lower, upper, costs)
    ]
