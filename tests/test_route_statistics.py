"""Final CNDP route snapshots and scenario comparison tables."""

import csv
import importlib.util
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
    runner.post_process_cndp_outputs(tmp_path, "COBYLA", tmp_path / "absent.toml")
    assert not raw.exists()
    assert snapshot.exists()
    # Root leftovers are removed, but organized snapshots remain available.
    raw.write_text("unused", encoding="utf-8")
    runner.cleanup_cndp_root(tmp_path)
    assert not raw.exists()
    path = runner.write_route_counts_table(tmp_path, ["COBYLA"])
    assert read_table(path)[1] == ["0", "0", "1", "2"]


@pytest.mark.parametrize("approach, expected_count", [("routebased", 2), ("tapas", 1)])
@pytest.mark.parametrize("diagnostics", [False, True])
def test_cndp_records_final_counts_without_changing_results(
    tmp_path, monkeypatch, approach, expected_count, diagnostics,
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
        "enable_trace": False, "write_metadata_json": False, "write_summary_csv": False,
    }, diagnostics)

    paths = list(tmp_path.glob("BilevelCND_*_final_route_counts.csv"))
    assert len(paths) == 1
    assert read_table(paths[0]) == [
        ["od_pair_index", "init_node", "dest_node", "routes_count"],
        ["0", "0", "1", str(expected_count)],
    ]
    if not diagnostics:
        assert list(tmp_path.iterdir()) == paths
    np.testing.assert_array_equal(recorded.flows, baseline.flows)
    np.testing.assert_array_equal(recorded.capacities, baseline.capacities)
    assert recorded.objective == baseline.objective
