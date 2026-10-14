"""Pinned TNTP inputs and independent checks for the benchmark study."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import io
import json
import math
from pathlib import Path
import re
import subprocess
from urllib.request import urlopen
from zipfile import ZipFile

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
    "Anaheim": {
        "net": ("Anaheim/Anaheim_net.tntp",
                "99933b415e9500b13907829c37a43cfa9141714fad5af279081e28e5f9356f9a"),
        "trips": ("Anaheim/Anaheim_trips.tntp",
                  "906893854cd0db4479c0b5f07678ce5616fa8e42e2b997f918c378309c66a94e"),
    },
    "Sydney": {
        "net": ("Sydney/Sydney_net.tntp",
                "cf7162f01eb9518b246e7a761f7bb02b3c154577f9a4514fed419fd957c08923"),
        "trips": ("Sydney/Sydney_trips.tntp.zip",
                  "85acc191d27c7c36c60d86d4f00f5c177caf01b799ca8add355e153b91a1ac33"),
    },
}
SYDNEY_MEMBER = "Sydney_trips.tntp"
SYDNEY_MEMBER_SHA256 = "6f98db94972c535ec0590ac911ccdc816d9db63f81705ac73593bbbe9afda735"
SYDNEY_DEMAND_SHA256 = "eee076405349109bcfea762e24a8478dc8b272c6c6d8fd9ec0e45dd99bca9778"
PREPARATION_VERSION = 1


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def array_sha256(array: np.ndarray) -> str:
    return sha256(np.asarray(array, dtype="<f8", order="C").tobytes())


def prepare_data(cache: Path, datasets: list[str], source_root: Path | None) -> dict:
    """Read committed blobs, never locally edited CSV/TNTP files, or fetch pinned URLs."""
    cache.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name in datasets:
        manifest[name] = {}
        for kind, (upstream_path, checksum) in SOURCES[name].items():
            suffix = ".tntp.zip" if upstream_path.endswith(".zip") else ".tntp"
            destination = cache / f"{name}_{kind}{suffix}"
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
        if name == "Sydney":
            manifest[name]["prepared"] = prepare_sydney(cache, manifest[name])
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
    details: dict = field(default_factory=dict)

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
            "routing": "unrestricted transit",
            **self.details,
        }


def parse_network(name: str, network_text: str) -> tuple[int, int, np.ndarray]:
    meta, body = _metadata(network_text)
    zones, nodes, count = (int(meta[k]) for k in
                           ("NUMBER OF ZONES", "NUMBER OF NODES", "NUMBER OF LINKS"))
    # The study intentionally uses the existing all-nodes-transit solver model.
    # FIRST THRU NODE is source metadata, not a routing constraint here.
    if not 0 < zones <= nodes or count <= 0:
        raise ValueError("Invalid TNTP dimensions")
    rows = [line.replace(";", " ").split() for line in body.splitlines() if line.strip()]
    if len(rows) != count or any(len(row) != 10 for row in rows):
        raise ValueError("TNTP link rows do not match declared dimensions or ten-column schema")
    links = np.asarray(rows, dtype=float)
    if not np.all(np.isfinite(links)) or np.any(links < 0) or np.any(links[:, 2] <= 0):
        raise ValueError("TNTP link values must be finite, nonnegative, with positive capacities")
    if name == "Sydney":
        # Its final columns are critical_speed and lanes, not toll and link_type.
        if not re.search(r"critical_speed\s+lanes", network_text):
            raise ValueError("Unexpected Sydney network columns")
        links[:, 8], links[:, 9] = 0, 1
    for col in (0, 1, 9):
        if np.any(links[:, col] != np.floor(links[:, col])):
            raise ValueError("Node IDs and link types must be integers")
    if np.any(links[:, :2] < 1) or np.any(links[:, :2] > nodes):
        raise ValueError("Link node ID outside declared range")
    return zones, nodes, links


def parse_trips(stream, zones: int, *, truncate: bool = False,
                allow_total_mismatch: bool = False) -> tuple[np.ndarray, dict]:
    """Parse one line at a time; duplicate tracking uses one byte per OD."""
    header = []
    for line in stream:
        header.append(line)
        if "<END OF METADATA>" in line.upper():
            break
    trip_meta, trailing = _metadata("".join(header))
    if trailing.strip():
        raise ValueError("Demand metadata must end before entries")
    if int(trip_meta["NUMBER OF ZONES"]) != zones:
        raise ValueError("Network and demand zone counts differ")
    demand = np.zeros((zones, zones))
    tokens = re.compile(r"Origin\s+(\d+)|(\d+)\s*:\s*([+\-\d.eE]+)\s*;", re.IGNORECASE)
    origin = None
    seen = np.zeros((zones, zones), dtype=bool)
    for line in stream:
        body = line.split("~", 1)[0]
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
            if seen[origin, destination]:
                raise ValueError("Duplicate demand entry")
            seen[origin, destination] = True
            demand[origin, destination] = value
        if body[end:].strip():
            raise ValueError("Malformed TNTP demand entry")
    if not seen.any():
        raise ValueError("Empty TNTP demand")
    actual = float(demand.sum())
    declared = float(trip_meta["TOTAL OD FLOW"]) if "TOTAL OD FLOW" in trip_meta else None
    if declared is not None and (not math.isfinite(declared) or declared < 0):
        raise ValueError("Invalid TNTP TOTAL OD FLOW")
    matches = declared is None or math.isclose(declared, actual, rel_tol=1e-8, abs_tol=1e-6)
    if not matches and not allow_total_mismatch:
        raise ValueError("Demand sum differs from TNTP TOTAL OD FLOW")
    if truncate:
        np.floor(demand, out=demand)
    return demand, {
        "demand_variant": "CSV-compatible truncated demand" if truncate else "fractional source demand",
        "declared_source_demand": declared, "actual_source_demand": actual,
        "source_total_matches_header": matches,
        "demand_transform": "floor each nonnegative OD entry" if truncate else "none",
    }


def parse_tntp(name: str, network_text: str, trips_text: str) -> NetworkInput:
    """Parse unrestricted, travel-time-only inputs, preserving fractional demand."""
    zones, nodes, links = parse_network(name, network_text)
    demand, details = parse_trips(io.StringIO(trips_text), zones)
    return NetworkInput(name, nodes, links, demand, details)


def prepare_sydney(cache: Path, sources: dict) -> dict:
    """Prepare the checksum-pinned legacy workload once, outside worker timing."""
    prepared = cache / "Sydney_prepared.npz"
    metadata = cache / "Sydney_prepared.json"
    if prepared.exists() or metadata.exists():
        if not (prepared.exists() and metadata.exists()):
            raise ValueError("Incomplete Sydney prepared cache; remove Sydney_prepared.* and retry")
        info = json.loads(metadata.read_text())
        if (info["preparation_version"] != PREPARATION_VERSION or info["sources"] != sources
                or sha256(prepared.read_bytes()) != info["sha256"]):
            raise ValueError("Sydney prepared cache mismatch; remove Sydney_prepared.* and retry")
        return info
    zones, nodes, links = parse_network("Sydney", (cache / "Sydney_net.tntp").read_text())
    with ZipFile(cache / "Sydney_trips.tntp.zip") as archive:
        if archive.namelist().count(SYDNEY_MEMBER) != 1:
            raise ValueError("Sydney archive must contain exactly one demand member")
        with archive.open(SYDNEY_MEMBER) as member:
            digest = hashlib.file_digest(member, "sha256").hexdigest()
        if digest != SYDNEY_MEMBER_SHA256:
            raise ValueError("Sydney demand member checksum mismatch")
        # The source checksum above narrowly authorizes this known header mismatch.
        with archive.open(SYDNEY_MEMBER) as member:
            demand, details = parse_trips(io.TextIOWrapper(member), zones,
                                          truncate=True, allow_total_mismatch=True)
    if array_sha256(demand) != SYDNEY_DEMAND_SHA256:
        raise ValueError("Prepared Sydney demand differs from the verified legacy CSV workload")
    details["network_mapping"] = "first eight source columns; toll=0, link_type=1"
    np.savez_compressed(prepared, n_nodes=nodes, links=links, demand=demand)
    info = {
        "preparation_version": PREPARATION_VERSION, "sources": sources.copy(),
        "archive_member": SYDNEY_MEMBER, "member_sha256": digest,
        "sha256": sha256(prepared.read_bytes()), "demand_sha256": array_sha256(demand),
        "links_sha256": array_sha256(links), "details": details,
    }
    metadata.write_text(json.dumps(info, indent=2) + "\n")
    return info


def load_input(cache: Path, name: str) -> NetworkInput:
    if name == "Sydney":
        info = json.loads((cache / "Sydney_prepared.json").read_text())
        prepared = cache / "Sydney_prepared.npz"
        if sha256(prepared.read_bytes()) != info["sha256"]:
            raise ValueError("Sydney prepared cache checksum mismatch")
        with np.load(prepared, allow_pickle=False) as arrays:
            links, demand = arrays["links"], arrays["demand"]
            if (array_sha256(demand) != SYDNEY_DEMAND_SHA256
                    or array_sha256(links) != info["links_sha256"]):
                raise ValueError("Sydney prepared array checksum mismatch")
            return NetworkInput(name, int(arrays["n_nodes"]), links, demand, info["details"])
    return parse_tntp(name, (cache / f"{name}_net.tntp").read_text(),
                      (cache / f"{name}_trips.tntp").read_text())


def shortest_path_total(data: NetworkInput, costs, batch_size: int = 64) -> tuple[float, int]:
    """Independent directed shortest paths, with bounded distance-matrix memory."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import dijkstra

    # CSR construction sums duplicates. Collapse parallel links by MIN instead.
    edges = {}
    for (u, v), cost in zip(data.links[:, :2].astype(int) - 1, costs):
        key = (int(u), int(v))
        edges[key] = min(float(cost), edges.get(key, math.inf))
    rows, columns = zip(*edges) if edges else ([], [])
    graph = csr_matrix((list(edges.values()), (rows, columns)), shape=(data.n_nodes, data.n_nodes))
    # Explicit sparse zeros represent real zero-cost links; do not eliminate them.
    origins = np.flatnonzero(data.demand.sum(axis=1))
    shortest_total, unreachable = 0.0, 0
    for start in range(0, len(origins), batch_size):
        batch = origins[start:start + batch_size]
        distances = dijkstra(graph, directed=True, indices=batch, return_predecessors=False)
        for i, origin in enumerate(batch):
            destinations = np.flatnonzero(data.demand[origin])
            lengths = distances[i, destinations]
            finite = np.isfinite(lengths)
            unreachable += int(np.count_nonzero(~finite))
            shortest_total += float(np.dot(data.demand[origin, destinations[finite]], lengths[finite]))
    return shortest_total, unreachable


def validate_assignment(data: NetworkInput, flows, costs, capacities, reported_tstt,
                        reported_gap, tolerance: float) -> dict:
    """Check feasibility and recompute the equilibrium gap independently with SciPy."""
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
    shortest_total, unreachable = shortest_path_total(data, bpr)
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
