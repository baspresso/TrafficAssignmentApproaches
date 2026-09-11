"""Input conversion for the native library (arrays and constraint CSV files)."""

from __future__ import annotations

import os
import numpy as np

from . import _core
from ._core import LinkConstraint, Network

def load_constraints(path: str | os.PathLike, verbose: bool = False) -> list[LinkConstraint]:
    """Load per-link capacity constraints from CSV, in network link order."""
    return _core._load_constraints(str(path), verbose)


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
    init = np.ascontiguousarray(np.asarray(init_node, dtype=np.int64))
    term = np.ascontiguousarray(np.asarray(term_node, dtype=np.int64))
    if init.ndim != 1 or term.shape != init.shape:
        raise ValueError("init_node/term_node must be 1-D arrays of equal length")
    n_links = init.shape[0]

    def per_link(values, dtype):
        return np.ascontiguousarray(np.broadcast_to(np.asarray(values, dtype=dtype), (n_links,)))

    demand = np.ascontiguousarray(np.asarray(demand, dtype=float))
    if demand.ndim != 2 or demand.shape[0] != demand.shape[1]:
        raise ValueError("demand must be a square (n_zones x n_zones) matrix")
    zones = int(n_zones) if n_zones is not None else demand.shape[0]
    if zones != demand.shape[0]:
        raise ValueError("n_zones must match the demand matrix dimension")
    if n_nodes is not None:
        nodes = int(n_nodes)
    elif n_links > 0:
        nodes = int(max(init.max(), term.max())) + 1
    else:
        nodes = zones

    return _core._build_network_from_arrays(
        str(name), nodes, zones, init, term,
        per_link(capacity, float), per_link(length, float),
        per_link(free_flow_time, float), per_link(b, float), per_link(power, float),
        per_link(speed, float), per_link(toll, float), per_link(link_type, np.int64),
        demand,
    )
