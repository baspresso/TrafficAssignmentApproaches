"""Pinned TNTP inputs and independent checks for the benchmark study."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import heapq
import math
from pathlib import Path
import re
import subprocess
from urllib.request import urlopen

import numpy as np

DATA_REVISION = "375e0da93858c547230c5cf9ea8a96de4ccff29e"
DATA_REPOSITORY = "https://github.com/bstabler/TransportationNetworks"
SOURCES = {
    "SiouxFalls": {
        "net": ("SiouxFalls/SiouxFalls_net.tntp",
                "ace99b24cec69c273ff0cf3d6d074110177f0cc0ae24b0c7a9f4f4cb5e27635c"),
        "trips": ("SiouxFalls/SiouxFalls_trips.tntp",
                  "56f9566857f3f66730fd5c4232258d7ee3ac2931a476526331afd062f4958de7"),
    },
    "ChicagoSketch": {
        "net": ("Chicago-Sketch/ChicagoSketch_net.tntp",
                "4396bff6101cb5ad3edaf0eb5b9aec051055cb7f0d85be907d24b43f98bf0027"),
        "trips": ("Chicago-Sketch/ChicagoSketch_trips.tntp",
                  "efe68abffc4af09e344cf1e175cfc048c08f4cd8f1f5454f74371b40e8245edc"),
    },
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def prepare_data(cache: Path, datasets: list[str], source_root: Path | None) -> dict:
    """Read committed blobs, never locally edited CSV/TNTP files, or fetch pinned URLs."""
    cache.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name in datasets:
        manifest[name] = {}
        for kind, (upstream_path, checksum) in SOURCES[name].items():
            destination = cache / f"{name}_{kind}.tntp"
            if destination.exists():
                data = destination.read_bytes()
            elif source_root is not None:
                data = subprocess.check_output(
                    ["git", "-C", str(source_root), "show", f"{DATA_REVISION}:{upstream_path}"],
                    stderr=subprocess.PIPE, timeout=30,
                )
            else:
                url = f"https://raw.githubusercontent.com/bstabler/TransportationNetworks/{DATA_REVISION}/{upstream_path}"
                with urlopen(url, timeout=30) as response:
                    data = response.read()
            if sha256(data) != checksum:
                raise ValueError(f"Checksum mismatch for {destination}; remove the cached file and retry")
            if not destination.exists():
                destination.write_bytes(data)
            manifest[name][kind] = {"upstream_path": upstream_path, "sha256": checksum}
    return manifest


def _metadata(text: str) -> tuple[dict[str, str], str]:
    marker = re.search(r"<END OF METADATA>", text, re.IGNORECASE)
    if marker is None:
        raise ValueError("TNTP input is missing <END OF METADATA>")
    metadata = {k.strip().upper(): v.strip() for k, v in
                re.findall(r"<([^>]+)>\s*([^\r\n]*)", text[:marker.start()])}
    body = "\n".join(line.split("~", 1)[0] for line in text[marker.end():].splitlines())
    return metadata, body


@dataclass
class NetworkInput:
    name: str
    n_nodes: int
    links: np.ndarray
    demand: np.ndarray

    def build(self, capacities=None):
        import traffic_assignment as ta
        a = self.links
        return ta.network_from_arrays(
            self.name, a[:, 0].astype(int) - 1, a[:, 1].astype(int) - 1,
            a[:, 2] if capacities is None else capacities, a[:, 4], a[:, 5], a[:, 6],
            self.demand, length=a[:, 3], speed=a[:, 7], toll=a[:, 8],
            link_type=a[:, 9].astype(int), n_nodes=self.n_nodes,
        )

    def dimensions(self) -> dict:
        return {
            "nodes": self.n_nodes, "links": len(self.links), "zones": len(self.demand),
            "positive_od_pairs": int(np.count_nonzero(self.demand)),
            "interzonal_od_pairs": int(np.count_nonzero(self.demand) - np.count_nonzero(np.diag(self.demand))),
            "active_origins": int(np.count_nonzero(self.demand.sum(axis=1))),
            "total_demand": float(self.demand.sum()),
            "intrazonal_demand": float(np.trace(self.demand)),
            "interzonal_demand": float(self.demand.sum() - np.trace(self.demand)),
        }


def parse_tntp(name: str, network_text: str, trips_text: str) -> NetworkInput:
    """Parse the study's travel-time-only networks without rounding demand."""
    meta, body = _metadata(network_text)
    zones, nodes, count = (int(meta[k]) for k in
                           ("NUMBER OF ZONES", "NUMBER OF NODES", "NUMBER OF LINKS"))
    if int(meta["FIRST THRU NODE"]) != 1:
        raise ValueError("This benchmark supports FIRST THRU NODE = 1 only")
    if not 0 < zones <= nodes or count <= 0:
        raise ValueError("Invalid TNTP dimensions")
    rows = [line.replace(";", " ").split() for line in body.splitlines() if line.strip()]
    if len(rows) != count or any(len(row) != 10 for row in rows):
        raise ValueError("TNTP link rows do not match declared dimensions or ten-column schema")
    links = np.asarray(rows, dtype=float)
    if not np.all(np.isfinite(links)) or np.any(links < 0) or np.any(links[:, 2] <= 0):
        raise ValueError("TNTP link values must be finite, nonnegative, with positive capacities")
    for col in (0, 1, 9):
        if np.any(links[:, col] != np.floor(links[:, col])):
            raise ValueError("Node IDs and link types must be integers")
    if np.any(links[:, :2] < 1) or np.any(links[:, :2] > nodes):
        raise ValueError("Link node ID outside declared range")

    trip_meta, body = _metadata(trips_text)
    if int(trip_meta["NUMBER OF ZONES"]) != zones:
        raise ValueError("Network and demand zone counts differ")
    demand = np.zeros((zones, zones))
    tokens = re.compile(r"Origin\s+(\d+)|(\d+)\s*:\s*([+\-\d.eE]+)\s*;", re.IGNORECASE)
    origin = None
    seen = set()
    end = 0
    for token in tokens.finditer(body):
        if body[end:token.start()].strip():
            raise ValueError("Malformed TNTP demand entry")
        end = token.end()
        if token[1] is not None:
            origin = int(token[1]) - 1
            if not 0 <= origin < zones:
                raise ValueError("Demand origin outside declared range")
            continue
        destination, value = int(token[2]) - 1, float(token[3])
        if origin is None or not 0 <= destination < zones:
            raise ValueError("Demand entry has no valid origin/destination")
        if not math.isfinite(value) or value < 0:
            raise ValueError("Demand must be finite and nonnegative")
        if (origin, destination) in seen:
            raise ValueError("Duplicate demand entry")
        seen.add((origin, destination))
        demand[origin, destination] = value
    if body[end:].strip() or not seen:
        raise ValueError("Malformed or empty TNTP demand")
    if "TOTAL OD FLOW" in trip_meta and not math.isclose(
        float(trip_meta["TOTAL OD FLOW"]), float(demand.sum()), rel_tol=1e-8, abs_tol=1e-6,
    ):
        raise ValueError("Demand sum differs from TNTP TOTAL OD FLOW")
    return NetworkInput(name, nodes, links, demand)


def load_input(cache: Path, name: str) -> NetworkInput:
    return parse_tntp(name, (cache / f"{name}_net.tntp").read_text(),
                      (cache / f"{name}_trips.tntp").read_text())


def validate_assignment(data: NetworkInput, flows, costs, capacities, reported_tstt,
                        reported_gap, tolerance: float) -> dict:
    """Check feasibility and recompute the equilibrium gap with Python Dijkstra."""
    flows, costs, capacities = (np.asarray(a, dtype=float) for a in (flows, costs, capacities))
    errors = []
    if any(a.shape != (len(data.links),) for a in (flows, costs, capacities)):
        return {"valid": False, "errors": ["result dimensions"]}
    if not all(np.all(np.isfinite(a)) for a in (flows, costs, capacities)):
        return {"valid": False, "errors": ["nonfinite result"]}
    if np.any(flows < -1e-8) or np.any(costs < 0) or np.any(capacities <= 0):
        return {"valid": False, "errors": ["invalid flow, cost, or capacity"]}
    links = data.links
    bpr = links[:, 4] * (1 + links[:, 5] * (np.maximum(flows, 0) / capacities) ** links[:, 6])
    if not np.allclose(costs, bpr, rtol=1e-10, atol=1e-10):
        errors.append("BPR costs")
    tstt = float(np.dot(flows, bpr))
    if not math.isfinite(reported_tstt) or not math.isclose(
        tstt, reported_tstt, rel_tol=1e-9, abs_tol=1e-7,
    ):
        errors.append("total travel time")
    init, term = (links[:, i].astype(int) - 1 for i in (0, 1))
    balance = np.bincount(init, weights=flows, minlength=data.n_nodes)
    balance -= np.bincount(term, weights=flows, minlength=data.n_nodes)
    expected = np.zeros(data.n_nodes)
    expected[:len(data.demand)] = data.demand.sum(axis=1) - data.demand.sum(axis=0)
    residual = float(np.max(np.abs(balance - expected)))
    if not np.allclose(balance, expected, rtol=1e-9, atol=1e-6):
        errors.append("node flow conservation")
    adjacency = [[] for _ in range(data.n_nodes)]
    for u, v, cost in zip(init, term, bpr):
        adjacency[u].append((int(v), float(cost)))
    shortest_total = 0.0
    unreachable = 0
    for origin in np.flatnonzero(data.demand.sum(axis=1)):
        distances = [math.inf] * data.n_nodes
        distances[origin] = 0.0
        queue = [(0.0, int(origin))]
        while queue:
            distance, node = heapq.heappop(queue)
            if distance > distances[node]:
                continue
            for dest, cost in adjacency[node]:
                candidate = distance + cost
                if candidate < distances[dest]:
                    distances[dest] = candidate
                    heapq.heappush(queue, (candidate, dest))
        for dest in np.flatnonzero(data.demand[origin]):
            if not math.isfinite(distances[dest]):
                unreachable += 1
            else:
                shortest_total += data.demand[origin, dest] * distances[dest]
    if unreachable:
        errors.append("unreachable demand")
    gap = 1 - shortest_total / tstt if tstt > 0 else math.inf
    # Float64 validation has a small absolute roundoff allowance.
    gap_ok = math.isfinite(gap) and abs(gap) <= tolerance + 1e-12
    if not math.isfinite(reported_gap) or abs(reported_gap) > tolerance + 1e-12:
        gap_ok = False
    if not gap_ok:
        errors.append("equilibrium tolerance")
    if not math.isfinite(reported_gap) or not math.isclose(gap, reported_gap, abs_tol=1e-10):
        errors.append("native/independent gap agreement")
    return {
        "valid": not errors, "errors": errors, "independent_relative_gap": float(gap),
        "max_flow_balance_residual": residual, "unreachable_od_pairs": unreachable,
        "recomputed_total_travel_time": tstt,
    }
