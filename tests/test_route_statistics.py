"""CNDP route history, final snapshots, and scenario comparison tables."""

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

import traffic_assignment as ta


@pytest.fixture
def runner():
    path = Path(__file__).resolve().parents[1] / "scripts" / "run_experiment.py"
    spec = importlib.util.spec_from_file_location("run_experiment", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_snapshot(run_dir, scenario, rows):
    path = run_dir / "solutions" / f"{scenario}_route_counts.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter=";")
        writer.writerow(["od_pair_index", "init_node", "dest_node", "routes_count"])
        writer.writerows(rows)
    return path


def read_table(path):
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.reader(file, delimiter=";"))


def assert_history_matches_trace(run_dir, expected_od_count):
    with next(run_dir.glob("BilevelCND_*_quality_time.csv")).open(newline="") as file:
        trace = list(csv.DictReader(file))
    with next(run_dir.glob("BilevelCND_*_route_counts_history.csv")).open(newline="") as file:
        history = list(csv.DictReader(file, delimiter=";"))
    snapshot = read_table(next(run_dir.glob("BilevelCND_*_route_counts.csv")))[1:]
    assert len(snapshot) == expected_od_count
    assert len(history) == len(trace) * expected_od_count
    for index, point in enumerate(trace):
        group = history[index * expected_od_count:(index + 1) * expected_od_count]
        for row, od in zip(group, snapshot):
            assert row["TracePointIndex"] == str(index)
            for field in ("RunId", "Scenario", "ElapsedTime(s)", "Phase", "Step"):
                assert row[field] == point[field]
            assert [row[f] for f in ("od_pair_index", "init_node", "dest_node")] == od[:3]
            assert int(row["routes_count"]) >= 0
    assert trace[-1]["Phase"] == "final"
    assert [row["routes_count"] for row in history[-expected_od_count:]] == [od[3] for od in snapshot]
    return trace, history


def test_table_aligns_od_pairs_and_preserves_scenario_order(runner, tmp_path):
    names = ["Missing", "Etalon", "COBYLA_run1", "COBYLA_run2", "Method; custom"]
    write_snapshot(tmp_path, "Etalon", [[0, 0, 1, 1], [1, 1, 0, 2]])
    write_snapshot(tmp_path, "COBYLA_run1", [[1, 1, 0, 4], [0, 0, 1, 0]])
    write_snapshot(tmp_path, "COBYLA_run2", [[0, 0, 1, 3], [1, 1, 0, 1]])
    write_snapshot(tmp_path, "Method; custom", [[0, 0, 1, 1], [1, 1, 0, 1]])

    path = runner.write_route_counts_table(tmp_path, names)

    assert path == tmp_path / "tables" / "route_counts.csv"
    assert read_table(path) == [
        ["od_pair_index", "init_node", "dest_node", *[f"routes_count_{n}" for n in names]],
        ["0", "0", "1", "", "1", "0", "3", "1"],
        ["1", "1", "0", "", "2", "4", "1", "1"],
    ]


def test_table_without_snapshots_has_only_header(runner, tmp_path):
    path = runner.write_route_counts_table(tmp_path, ["Failed", "TimedOut"])
    assert read_table(path) == [[
        "od_pair_index", "init_node", "dest_node", "routes_count_Failed", "routes_count_TimedOut",
    ]]


@pytest.mark.parametrize("rows", [[[0, 1, 0, 2]], [[0, 0, 1, 2], [1, 1, 0, 1]], []])
def test_table_rejects_incompatible_od_identities(runner, tmp_path, rows):
    write_snapshot(tmp_path, "A", [[0, 0, 1, 1]])
    write_snapshot(tmp_path, "B", rows)
    with pytest.raises(ValueError, match="Incompatible OD identities"):
        runner.write_route_counts_table(tmp_path, ["A", "B"])
    assert not (tmp_path / "tables" / "route_counts.csv").exists()


def test_table_rejects_duplicate_od_indices(runner, tmp_path):
    write_snapshot(tmp_path, "A", [[0, 0, 1, 1], [0, 1, 0, 1]])
    with pytest.raises(ValueError, match="duplicate OD row"):
        runner.write_route_counts_table(tmp_path, ["A"])


def test_postprocessing_preserves_snapshots_before_cleanup(runner, tmp_path):
    snapshot = write_snapshot(tmp_path, "COBYLA", [[0, 0, 1, 2]])
    raw = tmp_path / "BilevelCND_RouteBasedNewtonStep_123_route_counts.csv"
    snapshot.rename(raw)
    (tmp_path / "traces").mkdir()
    raw_history = tmp_path / "BilevelCND_RouteBasedNewtonStep_123_route_counts_history.csv"
    history_contents = "RunId;Scenario;TracePointIndex;ElapsedTime(s);Phase;Step;od_pair_index;init_node;dest_node;routes_count\n123;COBYLA;0;0.1;final;0;0;0;1;2\n"
    raw_history.write_text(history_contents)
    runner.post_process_cndp_outputs(tmp_path, "COBYLA", tmp_path / "absent.toml")
    assert not raw.exists()
    assert snapshot.exists()
    history = tmp_path / "traces" / "COBYLA_route_counts_history.csv"
    assert history.read_text() == history_contents
    assert not raw_history.exists()
    # Root leftovers are removed, but organized snapshots remain available.
    raw.write_text("unused", encoding="utf-8")
    raw_history.write_text("unused", encoding="utf-8")
    runner.cleanup_cndp_root(tmp_path)
    assert not raw.exists()
    assert not raw_history.exists()
    assert history.read_text() == history_contents
    path = runner.write_route_counts_table(tmp_path, ["COBYLA"])
    assert read_table(path)[1] == ["0", "0", "1", "2"]


@pytest.mark.parametrize("approach, expected_count", [("routebased", 2), ("tapas", 1)])
@pytest.mark.parametrize("diagnostics", [False, True])
@pytest.mark.parametrize("enable_trace", [False, True])
def test_cndp_records_final_counts_without_changing_results(
    tmp_path, monkeypatch, approach, expected_count, diagnostics, enable_trace,
):
    monkeypatch.chdir(tmp_path)

    def solve(metrics=None, final_diagnostics=False):
        network = ta.network_from_arrays(
            "TwoRoutes", init_node=[0, 0, 2], term_node=[1, 2, 1],
            capacity=10, free_flow_time=[2, 1, 1], b=1, power=1,
            demand=[[0, 10], [0, 0]], n_nodes=3,
        )
        return ta.solve_cndp(
            network, [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=1)],
            constraints=ta.constraints_from_arrays(network, 10, 20),
            approach=approach, metrics=metrics, final_diagnostics=final_diagnostics,
        )

    baseline = solve()
    assert not list(tmp_path.iterdir()), "disabled statistics must not produce files"
    recorded = solve({
        "output_root": str(tmp_path), "append_dataset_subdir": False, "run_id": "final",
        "enable_trace": enable_trace, "write_metadata_json": False, "write_summary_csv": False,
    }, diagnostics)

    paths = list(tmp_path.glob("BilevelCND_*_final_route_counts.csv"))
    assert len(paths) == 1
    assert read_table(paths[0]) == [
        ["od_pair_index", "init_node", "dest_node", "routes_count"],
        ["0", "0", "1", str(expected_count)],
    ]
    if enable_trace:
        trace, history = assert_history_matches_trace(tmp_path, 1)
        assert [row["Phase"] for row in trace] == ["standard_eval", "post_nlopt(LN_COBYLA)", "final"]
        assert history[-1]["routes_count"] == str(expected_count)
    else:
        assert not list(tmp_path.glob("*_route_counts_history.csv"))
    if not diagnostics and not enable_trace:
        assert list(tmp_path.iterdir()) == paths
    np.testing.assert_array_equal(recorded.flows, baseline.flows)
    np.testing.assert_array_equal(recorded.capacities, baseline.capacities)
    assert recorded.objective == baseline.objective


@pytest.mark.parametrize("approach", ["routebased", "tapas"])
@pytest.mark.parametrize("flush_period", [0, 2])
def test_route_history_covers_all_pipeline_phases(tmp_path, monkeypatch, approach, flush_period):
    monkeypatch.chdir(tmp_path)
    network = ta.network_from_arrays(
        "TwoODs", init_node=[0, 0, 2, 1], term_node=[1, 2, 1, 0],
        capacity=10, free_flow_time=[2, 1, 1, 2], b=1, power=1,
        demand=[[0, 10], [5, 0]], n_nodes=3,
    )
    pipeline = [
        ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=2),
        ta.step("optimlib", algorithm="NM", max_iterations=1),
        ta.step("gradient_descent", max_iterations=1),
        ta.step("optimality_condition", max_iterations=1),
        ta.step("optimality_condition", max_iterations=1),
    ]
    ta.solve_cndp(
        network, pipeline, constraints=ta.constraints_from_arrays(network, 10, 20),
        approach=approach, metrics={
            "output_root": str(tmp_path), "append_dataset_subdir": False,
            "run_id": 'history;"run"', "scenario_name": 'History; "custom"',
            "enable_trace": True, "enable_relative_gap": False,
            "flush_every_n_points": flush_period,
            "write_metadata_json": False, "write_summary_csv": False,
        },
    )
    trace, history = assert_history_matches_trace(tmp_path, 2)
    phases = [row["Phase"] for row in trace]
    assert {"standard_eval", "gd_eval", "optcond_init", "optcond_iter", "final"} <= set(phases)
    assert sum(phase.startswith("post_") for phase in phases) == len(pipeline)
    assert phases.count("optcond_init") == phases.count("optcond_iter") == 2
    assert all(row["Scenario"] == 'History; "custom"' for row in history)


def test_comparison_manifest_exposes_available_route_history(tmp_path, monkeypatch):
    scripts = Path(__file__).resolve().parents[1] / "scripts"
    monkeypatch.syspath_prepend(str(scripts))
    comparison = importlib.import_module("run_cndp_comparison")
    monkeypatch.setattr(comparison, "get_git_commit", lambda: "test")
    history = tmp_path / "traces" / "Recorded_route_counts_history.csv"
    history.parent.mkdir()
    history.write_text("history")
    results = [{
        "scenario": name, "description": "", "exit_code": 0, "elapsed_seconds": 1.0,
        "final_objective": 20.0, "final_budget": 5.0,
    } for name in ["Recorded", "WithoutTrace"]]
    comparison.write_run_manifest(tmp_path, "TwoODs", {}, results, 2.0)
    scenarios = json.loads((tmp_path / "run_manifest.json").read_text())["scenarios"]
    assert scenarios[0]["route_counts_history"] == "traces/Recorded_route_counts_history.csv"
    assert "route_counts_history" not in scenarios[1]


def test_plot_loader_ignores_route_history(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    scripts = Path(__file__).resolve().parents[1] / "scripts"
    monkeypatch.syspath_prepend(str(scripts))
    plotting = importlib.import_module("plot_cndp_comparison")
    traces = tmp_path / "traces"
    traces.mkdir()
    (traces / "COBYLA_trace.csv").write_text("Scenario,Phase,Step\nCOBYLA,final,0\n")
    (traces / "COBYLA_route_counts_history.csv").write_text("RunId;Scenario;routes_count\n123;COBYLA;2\n")
    loaded = plotting.load_trace_files_from_run_dir(tmp_path)
    assert list(loaded) == ["COBYLA"]
    assert list(loaded["COBYLA"].columns) == ["Scenario", "Phase", "Step"]
