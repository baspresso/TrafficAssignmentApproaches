import numpy as np
import pytest

import traffic_assignment as ta


def test_cndp_siouxfalls_smoke(data_root, tmp_path, monkeypatch):
    # Run from an empty directory to catch any accidental file output.
    monkeypatch.chdir(tmp_path)
    network = ta.load_network("SiouxFalls", data_root=data_root)
    constraints = ta.load_constraints(
        ta.default_constraints_path("SiouxFalls", data_root=data_root)
    )

    result = ta.solve_cndp(
        network,
        [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=5, tolerance=1e-3)],
        constraints=constraints,
        budget=10000.0,
    )

    assert np.isfinite(result.objective)
    assert result.objective == pytest.approx(result.total_travel_time + result.budget)
    assert result.budget <= result.budget_upper_bound + 1e-6
    assert np.all(result.capacities >= result.lower_bounds - 1e-9)
    assert np.all(result.capacities <= result.upper_bounds + 1e-9)
    # Caller-provided constraint objects must not be mutated by the
    # insensitive-link filter.
    assert any(c.upper_bound > c.lower_bound for c in constraints)
    assert not any(tmp_path.iterdir()), "CNDP solve must not write files when metrics are off"


def test_cndp_pipeline_dicts_and_auto_constraints(data_root):
    result = ta.solve_cndp(
        "SiouxFalls",
        [{"type": "nlopt", "algorithm": "LN_BOBYQA", "max_iterations": 3}],
        data_root=data_root,
        budget=5000.0,
    )

    assert np.isfinite(result.objective)
    assert result.network.name == "SiouxFalls"


def test_step_validation():
    with pytest.raises(ValueError):
        ta.step("nlopt", not_an_option=1)
    with pytest.raises(ValueError):
        ta.solve_cndp("SiouxFalls", [])  # empty pipeline rejected before any loading
