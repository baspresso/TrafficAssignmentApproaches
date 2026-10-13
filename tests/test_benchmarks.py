"""Benchmark evidence must preserve inputs and distinguish failure from speed."""

import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import benchmark_data as data
from benchmark_report import summarize
from run_benchmarks import collect, configurations, run_process, worker


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
    (NETWORK.replace("THRU NODE> 1", "THRU NODE> 3"), TRIPS),
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


def test_speedup_requires_complete_and_numerically_matching_runs():
    cases = [c | {"repeats": 1} for c in configurations("full")
             if c["dataset"] == "ChicagoSketch" and c["approach"] == "routebased"
             and c["tolerance"] == 1e-8 and c["threads"] in (1, 2)]
    serial, parallel = cases
    records = [measured(serial, 2), measured(parallel, 1)]
    assert summarize(cases, records)[1]["speedup"] == 2
    records[1]["status"] = "timeout"
    assert summarize(cases, records)[1]["speedup"] is None
    records[1]["status"] = "ok"
    records[1]["total_travel_time"] = 40
    assert summarize(cases, records)[1]["speedup"] is None
