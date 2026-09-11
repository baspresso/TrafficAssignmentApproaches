"""Dataset loading and programmatic network construction."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from . import _core
from .io import load_constraints, network_from_arrays
from ._core import LinkConstraint, Network

DATA_ROOT_ENV = "TRAFFIC_ASSIGNMENT_DATA"


def find_data_root(start: str | os.PathLike | None = None) -> Path:
    """Locate the TNTP data root (a ``data/TransportationNetworks`` directory).

    Resolution order: the ``TRAFFIC_ASSIGNMENT_DATA`` environment variable,
    then walking up from ``start`` (default: the current working directory) —
    mirroring the project-root discovery of the CLI executables.
    """
    env = os.environ.get(DATA_ROOT_ENV)
    if env:
        root = Path(env)
        if root.is_dir():
            return root
        raise FileNotFoundError(f"{DATA_ROOT_ENV} points to a missing directory: {root}")
    current = Path(start if start is not None else Path.cwd()).resolve()
    for candidate in (current, *current.parents):
        root = candidate / "data" / "TransportationNetworks"
        if root.is_dir():
            return root
    raise FileNotFoundError(
        f"Could not locate data/TransportationNetworks above {current}; "
        f"pass data_root= or set {DATA_ROOT_ENV}."
    )


def load_network(dataset: str, data_root: str | os.PathLike | None = None) -> Network:
    """Load a TNTP dataset (``{dataset}_net.csv`` + ``{dataset}_trips.csv``)."""
    root = Path(data_root) if data_root is not None else find_data_root()
    return _core._build_network_from_dataset(dataset, str(root))


def default_constraints_path(dataset: str, data_root: str | os.PathLike | None = None) -> Path:
    """Conventional constraints CSV location: ``{root}/{dataset}/{dataset}_constraints.csv``."""
    root = Path(data_root) if data_root is not None else find_data_root()
    return root / dataset / f"{dataset}_constraints.csv"


def constraints_from_arrays(network: Network, lower, upper,
                            investment_cost=1.0) -> list[LinkConstraint]:
    """Build per-link constraints from arrays (or scalars broadcast to all links)."""
    n_links = network.number_of_links
    lower = np.broadcast_to(np.asarray(lower, dtype=float), (n_links,))
    upper = np.broadcast_to(np.asarray(upper, dtype=float), (n_links,))
    cost = np.broadcast_to(np.asarray(investment_cost, dtype=float), (n_links,))
    if np.any(lower > upper):
        raise ValueError("lower bound exceeds upper bound for some links")
    nodes = network.link_nodes()
    return [
        LinkConstraint(int(nodes[i, 0]), int(nodes[i, 1]),
                       float(lower[i]), float(upper[i]), float(cost[i]))
        for i in range(n_links)
    ]
