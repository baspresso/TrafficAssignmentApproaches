"""Dataset discovery and compatibility resolution for solver shortcuts."""

from __future__ import annotations

import os
from pathlib import Path

from . import _core
from .io import constraints_from_arrays, load_constraints, network_from_arrays
from ._core import Network

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


def _resolve_network(network, data_root=None) -> Network:
    if isinstance(network, str):
        network = load_network(network, data_root=data_root)
    if not isinstance(network, Network):
        raise TypeError("network must be a Network or a dataset name")
    return network


def _resolve_constraints(network: Network, constraints, data_root=None):
    if constraints is None:
        constraints = default_constraints_path(network.name, data_root=data_root)
    if isinstance(constraints, (str, os.PathLike)):
        path = Path(constraints)
        if not path.is_file():
            raise FileNotFoundError(f"Constraints file does not exist: {path}")
        constraints = load_constraints(path)
    return constraints
