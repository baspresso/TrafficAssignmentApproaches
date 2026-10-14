"""Benchmark evidence must preserve inputs and distinguish failure from speed."""

import heapq
import json
import math
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import benchmark_data as data
from benchmark_report import render, summarize
from run_benchmarks import collect, configurations, run_process, schedule, worker


NETWORK = """<NUMBER OF ZONES> 2
<NUMBER OF NODES> 3
<FIRST THRU NODE> 1
<NUMBER OF LINKS> 3
<END OF METADATA>
~ init term capacity length free_flow b power speed toll type
1 2 10 0 2 1 1 0 0 1;
1 3 10 0 1 1 1 0 0 1;
3 2 10 0 1 1 1 0 0 1;
"""
TRIPS = """<NUMBER OF ZONES> 2
<TOTAL OD FLOW> 10
<END OF METADATA>
Origin 1
2 : 1e1;
"""


def test_tntp_preserves_fractional_demand_and_row_order():
    trips = TRIPS.replace("FLOW> 10", "FLOW> 10.25").replace("1e1", "1.025E+1")
    network = data.parse_tntp("Tiny", NETWORK, trips)
    assert network.demand[0, 1] == 10.25
    np.testing.assert_array_equal(network.links[:, :2], [[1, 2], [1, 3], [3, 2]])
    assert network.dimensions()["positive_od_pairs"] == 1
    assert network.build().number_of_links == 3


def test_intrazonal_demand_is_retained_without_road_flow():
    import traffic_assignment as ta
    trips = TRIPS.replace("FLOW> 10", "FLOW> 12.5").replace("2 : 1e1;", "1 : 2.5; 2 : 1e1;")
    network = data.parse_tntp("Tiny", NETWORK, trips)
    assert network.dimensions()["intrazonal_demand"] == 2.5
    assert network.dimensions()["interzonal_demand"] == 10
    assert network.dimensions()["interzonal_od_pairs"] == 1
    for approach in ("tapas", "routebased"):
        result = ta.solve_tap(network.build(), approach=approach)
        np.testing.assert_allclose(result.flows, [5, 5, 5], atol=1e-7)
        assert data.validate_assignment(
            network, result.flows, result.link_costs, [10, 10, 10],
            result.total_travel_time, result.relative_gap, 1e-8,
        )["valid"]


@pytest.mark.parametrize(("net", "trips"), [
    (NETWORK.replace("LINKS> 3", "LINKS> 4"), TRIPS),
    (NETWORK.replace("1 2 10", "1 2 nan"), TRIPS),
    (NETWORK.replace("1 2 10", "1.5 2 10"), TRIPS),
    (NETWORK, TRIPS.replace("2 : 1e1;", "2 : 1e1; garbage")),
    (NETWORK, TRIPS.replace("2 : 1e1;", "2 : 1e1; 2 : 0;")),
    (NETWORK, TRIPS.replace("2 : 1e1;", "3 : 1e1;")),
    (NETWORK, TRIPS.replace("1e1", "-10")),
    (NETWORK, TRIPS.replace("1e1", "1e999")),
    (NETWORK, TRIPS.replace("FLOW> 10", "FLOW> 11")),
    (NETWORK, TRIPS.replace("ZONES> 2", "ZONES> 3")),
])
def test_tntp_rejects_unsupported_or_invalid_input(net, trips):
    with pytest.raises((ValueError, KeyError)):
        data.parse_tntp("Tiny", net, trips)


def test_zone_nodes_remain_available_for_intermediate_transit():
    import traffic_assignment as ta
    net = """<NUMBER OF ZONES> 3
<NUMBER OF NODES> 4
<FIRST THRU NODE> 4
<NUMBER OF LINKS> 4
<END OF METADATA>
1 2 100 0 1 0 1 0 0 1;
2 3 100 0 1 0 1 0 0 1;
1 4 100 0 3 0 1 0 0 1;
4 3 100 0 3 0 1 0 0 1;
"""
    trips = TRIPS.replace("ZONES> 2", "ZONES> 3").replace("2 : 1e1;", "3 : 1e1;")
    network = data.parse_tntp("Zones", net, trips)
    for approach in ("tapas", "routebased"):
        result = ta.solve_tap(network.build(), approach=approach)
        np.testing.assert_allclose(result.flows, [10, 10, 0, 0])
        assert data.validate_assignment(network, result.flows, result.link_costs,
                                        [100] * 4, result.total_travel_time,
                                        result.relative_gap, 1e-8)["valid"]


def sydney_fixture(cache, monkeypatch):
    net = NETWORK.replace("speed toll type", "speed critical_speed lanes").replace("0 0 1;", "0 19.2 2;")
    trips = TRIPS.replace("FLOW> 10", "FLOW> 999").replace("1e1", "10.25") + "Origin 2\n1 : 0.75;\n"
    (cache / "Sydney_net.tntp").write_text(net)
    archive_path = cache / "Sydney_trips.tntp.zip"
    with ZipFile(archive_path, "w") as archive:
        archive.writestr("Sydney_trips.tntp", trips)
        archive.writestr("__MACOSX/._Sydney_trips.tntp", "ignored resource fork")
    sources = {
        "net": ("Sydney/Sydney_net.tntp", data.sha256(net.encode())),
        "trips": ("Sydney/Sydney_trips.tntp.zip", data.sha256(archive_path.read_bytes())),
    }
    monkeypatch.setattr(data, "SOURCES", {"Sydney": sources})
    monkeypatch.setattr(data, "SYDNEY_MEMBER_SHA256", data.sha256(trips.encode()))
    monkeypatch.setattr(data, "SYDNEY_DEMAND_SHA256", data.array_sha256(np.array([[0., 10.], [0., 0.]])))
    return net, trips


def test_sydney_archive_conversion_and_checked_cache(tmp_path, monkeypatch):
    net, trips = sydney_fixture(tmp_path, monkeypatch)
    # Header discrepancies remain errors outside the verified preparation path.
    with pytest.raises(ValueError, match="TOTAL OD FLOW"):
        data.parse_tntp("Sydney", net, trips)
    manifest = data.prepare_data(tmp_path, ["Sydney"], None)
    network = data.load_input(tmp_path, "Sydney")
    np.testing.assert_array_equal(network.demand, [[0, 10], [0, 0]])
    np.testing.assert_array_equal(network.links[:, 8:], [[0, 1]] * 3)
    assert network.details["declared_source_demand"] == 999
    assert network.details["actual_source_demand"] == 11
    assert network.dimensions()["total_demand"] == 10
    assert not network.details["source_total_matches_header"]
    assert manifest["Sydney"]["prepared"]["demand_sha256"] == data.SYDNEY_DEMAND_SHA256
    assert data.prepare_data(tmp_path, ["Sydney"], None) == manifest
    assert not (tmp_path / "__MACOSX").exists()
    (tmp_path / "Sydney_prepared.npz").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="cache mismatch"):
        data.prepare_data(tmp_path, ["Sydney"], None)
    with pytest.raises(ValueError, match="cache checksum"):
        data.load_input(tmp_path, "Sydney")


@pytest.mark.parametrize("corruption", ["archive", "member", "prepared_demand"])
def test_sydney_preparation_rejects_checksum_failures(tmp_path, monkeypatch, corruption):
    sydney_fixture(tmp_path, monkeypatch)
    if corruption == "archive":
        (tmp_path / "Sydney_trips.tntp.zip").write_bytes(b"bad zip")
    elif corruption == "member":
        monkeypatch.setattr(data, "SYDNEY_MEMBER_SHA256", "wrong")
    else:
        monkeypatch.setattr(data, "SYDNEY_DEMAND_SHA256", "wrong")
    with pytest.raises(ValueError, match="(?i)checksum|verified legacy"):
        data.prepare_data(tmp_path, ["Sydney"], None)


def test_anaheim_keeps_fractional_demand_with_unrestricted_routing():
    net = NETWORK.replace("THRU NODE> 1", "THRU NODE> 3")
    trips = TRIPS.replace("FLOW> 10", "FLOW> 10.25").replace("1e1", "10.25")
    parsed = data.parse_tntp("Anaheim", net, trips)
    assert parsed.demand[0, 1] == 10.25
    assert parsed.dimensions()["routing"] == "unrestricted transit"


@pytest.mark.parametrize("batch_size", [1, 2, 64])
def test_batched_dijkstra_matches_reference_with_parallel_zero_and_unreachable_links(batch_size):
    # Node IDs are one-based in NetworkInput. Node 5 cannot reach node 1.
    links = np.zeros((6, 10))
    links[:, :2] = [[1, 2], [1, 2], [2, 3], [1, 3], [3, 1], [4, 5]]
    costs = [3, 0, 2, 9, 7, 1]
    demand = np.zeros((5, 5))
    demand[0, 2], demand[2, 1], demand[3, 4] = 2, 3, 4
    demand[4, 0], demand[4, 4] = 5, 6
    network = data.NetworkInput("Paths", 5, links, demand)
    expected, unreachable = 0., 0
    adjacency = [[] for _ in range(5)]
    for (u, v), cost in zip(links[:, :2].astype(int) - 1, costs):
        adjacency[u].append((v, cost))
    for origin in np.flatnonzero(demand.sum(axis=1)):
        lengths = [math.inf] * 5
        lengths[origin] = 0.
        queue = [(0., origin)]
        while queue:
            distance, node = heapq.heappop(queue)
            if distance != lengths[node]:
                continue
            for dest, cost in adjacency[node]:
                if distance + cost < lengths[dest]:
                    lengths[dest] = distance + cost
                    heapq.heappush(queue, (lengths[dest], dest))
        for dest in np.flatnonzero(demand[origin]):
            if math.isfinite(lengths[dest]):
                expected += demand[origin, dest] * lengths[dest]
            else:
                unreachable += 1
    assert data.shortest_path_total(network, costs, batch_size) == (expected, unreachable) == (29., 1)


def test_full_and_quick_experiment_matrix():
    full = configurations("full")
    assert len(full) == len({c["id"] for c in full}) == 34
    assert sum(c["repeats"] for c in full) == 168
    assert sum(warmup for _, _, warmup in schedule(full)) == 34
    assert len(list(schedule(full))) == 202
    for name in ("Anaheim", "Sydney"):
        cases = [c for c in full if c["dataset"] == name]
        assert len(cases) == 9
        assert {c["threads"] for c in cases} == {1, 2, 4, 8}
        assert {c["tolerance"] for c in cases} == {1e-6, 1e-8, 1e-10}
        assert all(c["repeats"] == 5 for c in cases)
        assert all(c["timeout"] == (600 if name == "Sydney" else 120) for c in cases)
    assert [c["dataset"] for c in full if c["kind"] == "cndp"] == ["SiouxFalls"]
    quick = configurations("quick")
    assert len(quick) == 3
    assert all(c["dataset"] == "SiouxFalls" and c["repeats"] == 1 for c in quick)


def test_pinned_input_uses_committed_blob_and_checks_hash(tmp_path, monkeypatch):
    payload = b"pinned data"
    monkeypatch.setattr(data, "SOURCES", {"Tiny": {"net": ("upstream/file", data.sha256(payload))}})
    calls = []

    def read_blob(command, **kwargs):
        calls.append(command)
        return payload

    monkeypatch.setattr(subprocess, "check_output", read_blob)
    data.prepare_data(tmp_path, ["Tiny"], Path("/clone"))
    assert calls[0][-1] == f"{data.DATA_REVISION}:upstream/file"
    (tmp_path / "Tiny_net.tntp").write_bytes(b"locally changed")
    with pytest.raises(ValueError, match="Checksum mismatch"):
        data.prepare_data(tmp_path, ["Tiny"], Path("/clone"))


def test_independent_validation_detects_bad_flows_costs_and_gap():
    network = data.parse_tntp("Tiny", NETWORK, TRIPS)
    good = dict(flows=[5, 5, 5], costs=[3, 1.5, 1.5], capacities=[10, 10, 10],
                reported_tstt=30, reported_gap=0, tolerance=1e-8)
    assert data.validate_assignment(network, **good)["valid"]
    for change, error in [
        ({"flows": [6, 5, 5]}, "node flow conservation"),
        ({"costs": [4, 1.5, 1.5]}, "BPR costs"),
        ({"reported_gap": 1e-3}, "equilibrium tolerance"),
        ({"reported_tstt": 40}, "total travel time"),
    ]:
        result = data.validate_assignment(network, **(good | change))
        assert not result["valid"]
        assert error in result["errors"]


def test_independent_validation_detects_unreachable_demand():
    network = data.parse_tntp("Tiny", NETWORK, TRIPS)
    network.links = network.links[1:2]  # Only node 1 -> node 3 remains.
    result = data.validate_assignment(network, [0], [1], [10], 0, 0, 1e-8)
    assert not result["valid"]
    assert result["unreachable_od_pairs"] == 1


def test_worker_tap_and_cndp_validate_real_synthetic_solutions(tmp_path):
    (tmp_path / "Tiny_net.tntp").write_text(NETWORK)
    (tmp_path / "Tiny_trips.tntp").write_text(TRIPS)
    for case in (configurations("quick")[0], configurations("quick")[-1]):
        case = case | {"dataset": "Tiny"}
        if case["kind"] == "cndp":
            case.update(theta=0.1, budget=2)
        result = worker(case, tmp_path)
        assert result["status"] == "ok", result
        assert result["solve_seconds"] >= 0


def test_worker_timeout_and_crash_are_recorded():
    result = run_process([sys.executable, "-c", "import time; time.sleep(10)"], timeout=0.05)
    assert result["status"] == "timeout"
    result = run_process([sys.executable, "-c", "raise RuntimeError('broken')"], timeout=10)
    assert result["status"] == "error"
    assert "broken" in result["error"]
    assert run_process([sys.executable, "-c", "print('invalid')"], timeout=10)["status"] == "error"


def test_collect_keeps_completed_results_when_interrupted(tmp_path):
    calls = 0

    def interrupted(command, timeout):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise KeyboardInterrupt
        return {"status": "ok", "solve_seconds": 1}

    case = configurations("quick")[0]
    with pytest.raises(KeyboardInterrupt):
        collect([case], tmp_path, tmp_path, 10, run=interrupted)
    records = [json.loads(line) for line in (tmp_path / "measurements.jsonl").read_text().splitlines()]
    assert len(records) == 1
    assert records[0]["warmup"]


def test_collect_records_unrun_cases_after_deadline(tmp_path):
    def forbidden(*args):
        pytest.fail("Expired budget must not start a worker")

    records = collect([configurations("quick")[0]], tmp_path, tmp_path, 0, run=forbidden)
    assert len(records) == 2
    assert all(r["status"] == "budget_exhausted" for r in records)


def measured(case, seconds, repeat=1, **extra):
    return dict(case_id=case["id"], repeat=repeat, warmup=False, status="ok",
                solve_seconds=seconds, total_travel_time=30, beckmann_objective=25,
                relative_gap=1e-10, **extra)


def test_summary_excludes_warmup_and_reports_quartiles():
    case = configurations("full")[0]
    records = [measured(case, x, repeat=i) for i, x in enumerate([1, 2, 3, 4, 5], 1)]
    records.append(measured(case, 9999) | {"warmup": True})
    row = summarize([case], records)[0]
    assert row["all_valid"]
    assert (row["solve_seconds"], row["q25_seconds"], row["q75_seconds"]) == (3, 2, 4)


@pytest.mark.parametrize("dataset", ["ChicagoSketch", "Anaheim", "Sydney"])
def test_speedup_requires_complete_and_numerically_matching_runs(dataset):
    cases = [c | {"repeats": 1} for c in configurations("full")
             if c["dataset"] == dataset and c["approach"] == "routebased"
             and c["tolerance"] == 1e-8 and c["threads"] in (1, 2)]
    serial, parallel = cases
    records = [measured(serial, 2), measured(parallel, 1)]
    assert summarize(cases, records)[1]["speedup"] == 2
    records[1]["status"] = "timeout"
    assert summarize(cases, records)[1]["speedup"] is None
    records[1]["status"] = "ok"
    records[1]["total_travel_time"] = 40
    assert summarize(cases, records)[1]["speedup"] is None


def test_scaling_never_uses_another_networks_serial_baseline():
    cases = [c | {"repeats": 1} for c in configurations("full")
             if c["dataset"] in ("Anaheim", "Sydney") and c["approach"] == "routebased"
             and c["tolerance"] == 1e-8 and c["threads"] in (1, 2)]
    times = {("Anaheim", 1): 2, ("Anaheim", 2): 1, ("Sydney", 1): 12, ("Sydney", 2): 4}
    records = [measured(c, times[c["dataset"], c["threads"]]) for c in cases]
    rows = summarize(cases, records)
    assert [r["speedup"] for r in rows if r["threads"] == 2] == [2, 3]
    cases = [c for c in cases if (c["dataset"], c["threads"]) != ("Sydney", 1)]
    assert next(r for r in summarize(cases, records)
                if r["dataset"] == "Sydney")["speedup"] is None


def test_version_one_manifest_remains_renderable(tmp_path, monkeypatch):
    import benchmark_report
    monkeypatch.setattr(benchmark_report, "plots", lambda *args: [])
    cases = configurations("quick")
    manifest = dict(
        schema_version=1, started_utc="2026-01-01T00:00:00", profile="quick",
        time_budget_seconds=3600, data_repository=data.DATA_REPOSITORY,
        data_revision=data.DATA_REVISION, configuration=cases,
        datasets={"SiouxFalls": data.parse_tntp("Tiny", NETWORK, TRIPS).dimensions()},
        cpu="Test CPU", available_cpus=1, platform="Test platform", build={},
        python="3.11", packages={"numpy": "1.24"}, git_commit="test",
        working_tree_dirty=False, source_sha256="test",
    )
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    (tmp_path / "measurements.jsonl").write_text("\n".join(json.dumps(measured(c, 1)) for c in cases))
    render(tmp_path)
    text = (tmp_path / "README.md").read_text()
    assert "Python heap-based Dijkstra" in text
    assert "--time-budget 3600" in text
    assert "SciPy: not used" in text
