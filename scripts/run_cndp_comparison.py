"""
CNDP Algorithm Comparison Runner

Reads a TOML comparison manifest with inline scenarios and runs them.
Imports shared utilities from run_experiment.py to avoid duplication.

Usage:
    python scripts/run_cndp_comparison.py configs/comparisons/siouxfalls.toml
    python scripts/run_cndp_comparison.py configs/comparisons/siouxfalls.toml --skip-etalon
    python scripts/run_cndp_comparison.py configs/comparisons/siouxfalls.toml --scenarios OptCond,COBYLA
    python scripts/run_cndp_comparison.py configs/comparisons/siouxfalls.toml --repeat 3
    python scripts/run_cndp_comparison.py configs/comparisons/siouxfalls.toml --etalon-only
    python scripts/run_cndp_comparison.py configs/comparisons/siouxfalls.toml --exclude DE-only
"""

import argparse
import copy
import csv
import json
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        print("ERROR: Python 3.11+ or 'pip install tomli' required for TOML support")
        sys.exit(1)

try:
    import tomli_w
except ModuleNotFoundError:
    tomli_w = None

# Ensure scripts/ is on sys.path for sibling imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_experiment import (
    stream_subprocess,
    create_run_dir,
    post_process_cndp_outputs,
    cleanup_cndp_root,
    resolve_exe,
    get_git_commit,
    PROJECT_ROOT,
    DEFAULT_CNDP_EXE,
)


# ---------------------------------------------------------------------------
# TOML generation helpers
# ---------------------------------------------------------------------------

def _write_toml(data: dict, path: Path):
    """Write a dict as TOML. Uses tomli_w if available, else manual formatting."""
    if tomli_w is not None:
        path.write_bytes(tomli_w.dumps(data).encode())
        return

    # Manual TOML writer for the subset of types we use
    lines = []
    _write_toml_section(lines, data, [])
    path.write_text("\n".join(lines))


def _write_toml_section(lines: list[str], data: dict, prefix: list[str]):
    """Recursively write TOML sections."""
    # First pass: write scalar values
    for key, value in data.items():
        if isinstance(value, dict) or isinstance(value, list):
            continue
        lines.append(f"{key} = {_toml_value(value)}")

    # Second pass: write sub-tables
    for key, value in data.items():
        if isinstance(value, dict):
            section = prefix + [key]
            lines.append(f"\n[{'.'.join(section)}]")
            _write_toml_section(lines, value, section)

    # Third pass: write array-of-tables
    for key, value in data.items():
        if isinstance(value, list):
            section = prefix + [key]
            for item in value:
                lines.append(f"\n[[{'.'.join(section)}]]")
                _write_toml_section(lines, item, section)


def _toml_value(v) -> str:
    """Format a Python value as a TOML literal."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return str(v)
    if isinstance(v, str):
        return f'"{v}"'
    return str(v)


def _deep_merge(base: dict, override: dict) -> dict:
    """Deep-merge override into a copy of base. Override values win."""
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def generate_scenario_toml(defaults: dict, scenario: dict, output_dir: Path) -> Path:
    """Generate a standalone TOML config file for a scenario.

    Deep-merges [defaults] with scenario-level overrides, writes [[scenarios.pipeline]]
    as [[pipeline]], and returns the path to the generated file.
    """
    # Start from defaults
    config = copy.deepcopy(defaults)

    # Apply per-scenario overrides (e.g., [scenarios.solver])
    for section in ("network", "solver", "metrics", "output"):
        if section in scenario:
            if section in config:
                config[section] = _deep_merge(config[section], scenario[section])
            else:
                config[section] = copy.deepcopy(scenario[section])

    # Pipeline comes from the scenario
    config["pipeline"] = copy.deepcopy(scenario.get("pipeline", []))

    # Write TOML
    safe_name = scenario["name"].replace(" ", "_").replace("/", "_")
    toml_path = output_dir / f"{safe_name}.toml"
    _write_toml(config, toml_path)
    return toml_path


# ---------------------------------------------------------------------------
# TOML loading & validation
# ---------------------------------------------------------------------------

def load_comparison(path: Path) -> dict:
    """Read and validate TOML comparison manifest."""
    if not path.exists():
        print(f"ERROR: Comparison manifest not found: {path}")
        sys.exit(1)

    with open(path, "rb") as f:
        config = tomllib.load(f)

    if "comparison" not in config:
        print(f"ERROR: Missing [comparison] section in {path}")
        sys.exit(1)

    if "scenarios" not in config or not config["scenarios"]:
        print(f"ERROR: No [[scenarios]] entries in {path}")
        sys.exit(1)

    for i, scenario in enumerate(config["scenarios"]):
        if "name" not in scenario:
            print(f"ERROR: Scenario #{i+1} missing 'name' field")
            sys.exit(1)

        # Inline scenarios must have [[scenarios.pipeline]]
        if "pipeline" not in scenario:
            print(f"ERROR: Scenario '{scenario['name']}' missing [[scenarios.pipeline]] entries")
            sys.exit(1)

        scenario.setdefault("description", "")
        scenario.setdefault("etalon", False)
        scenario.setdefault("enabled", True)

    return config


def filter_scenarios(scenarios: list[dict], args) -> list[dict]:
    """Apply CLI filters: --scenarios, --etalon-only, --skip-etalon, --exclude, enabled."""
    result = [s for s in scenarios if s.get("enabled", True)]

    if args.etalon_only:
        result = [s for s in result if s.get("etalon", False)]
    elif args.skip_etalon:
        result = [s for s in result if not s.get("etalon", False)]

    if args.scenarios:
        names = {n.strip() for n in args.scenarios.split(",")}
        result = [s for s in result if s["name"] in names]

    if args.exclude:
        excluded = {n.strip() for n in args.exclude.split(",")}
        result = [s for s in result if s["name"] not in excluded]

    return result


def expand_repeats(scenarios: list[dict], repeat: int) -> list[dict]:
    """For --repeat N: duplicate non-etalon scenarios with _runN suffixes."""
    if repeat <= 1:
        return scenarios

    expanded = []
    for s in scenarios:
        if s.get("etalon", False):
            expanded.append(s)
        else:
            for run_idx in range(1, repeat + 1):
                copy_s = copy.deepcopy(s)
                copy_s["name"] = f"{s['name']}_run{run_idx}"
                copy_s["description"] = f"{s.get('description', '')} (run {run_idx})"
                expanded.append(copy_s)
    return expanded


# ---------------------------------------------------------------------------
# Batch progress rendering
# ---------------------------------------------------------------------------

class BatchProgressRenderer:
    """Two-tier progress: batch header + per-scenario step bars."""

    def start_scenario(self, index: int, total: int, scenario: dict):
        name = scenario["name"]
        desc = scenario.get("description", "")
        etalon_tag = " [ETALON]" if scenario.get("etalon") else ""
        print(f"\n{'='*70}")
        print(f"  [{index}/{total}] {name}{etalon_tag}")
        if desc:
            print(f"  {desc}")
        print(f"{'='*70}")

    def finish_scenario(self, scenario: dict, result: dict):
        name = scenario["name"]
        elapsed = result["elapsed_seconds"]
        obj = result.get("final_objective")
        status = "OK" if result["exit_code"] == 0 else "FAILED"
        obj_str = f"  obj={obj:.2f}" if obj is not None else ""
        print(f"  -> {name}  {status}  {elapsed:.1f}s{obj_str}")


# ---------------------------------------------------------------------------
# Scenario execution
# ---------------------------------------------------------------------------

def run_scenario(
    scenario: dict,
    run_dir: Path,
    exe: Path,
    time_limit: Optional[int] = None,
) -> dict:
    """Run single scenario: launch exe, stream output, return results dict."""
    config_path = scenario["config_path"]

    cmd = [
        str(exe),
        "--config", str(config_path),
        "--metrics_scenario_name", scenario["name"],
        "--metrics_output_root", str(run_dir),
        "--metrics_no_dataset_subdir", "true",
        "--quiet", "true",
    ]

    start_time = time.time()
    exit_code, _, result_dict = stream_subprocess(
        cmd, time_limit=time_limit, render_progress=True,
    )
    elapsed = time.time() - start_time

    final_objective = None
    final_travel_time = None
    final_budget = None
    if result_dict:
        try:
            final_objective = float(result_dict.get("objective_function", ""))
        except (ValueError, TypeError):
            pass
        try:
            final_travel_time = float(result_dict.get("total_travel_time", ""))
        except (ValueError, TypeError):
            pass
        try:
            final_budget = float(result_dict.get("budget_function", ""))
        except (ValueError, TypeError):
            pass

    return {
        "scenario": scenario["name"],
        "description": scenario.get("description", ""),
        "exit_code": exit_code,
        "elapsed_seconds": elapsed,
        "final_objective": final_objective,
        "final_travel_time": final_travel_time,
        "final_budget": final_budget,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_comparison_results(results: list[dict], run_dir: Path):
    """Write comparison results to CSV."""
    output_dir = run_dir / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "summary.csv"

    fieldnames = [
        "scenario", "description", "exit_code",
        "elapsed_seconds", "final_objective", "final_travel_time", "final_budget",
    ]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    print(f"\nComparison results written to: {output_path}")


def write_run_manifest(run_dir: Path, dataset: str, comparison_config: dict,
                       results: list[dict], total_wall_seconds: float):
    """Write machine-readable run_manifest.json."""
    scenarios_out = []
    for r in results:
        entry = {
            "name": r["scenario"],
            "description": r["description"],
            "exit_code": r["exit_code"],
            "elapsed_seconds": round(r["elapsed_seconds"], 2),
            "final_objective": r["final_objective"],
            "final_budget": r["final_budget"],
        }
        trace = run_dir / "traces" / f"{r['scenario']}_trace.csv"
        if trace.exists():
            entry["trace"] = f"traces/{r['scenario']}_trace.csv"
        solution = run_dir / "solutions" / f"{r['scenario']}_solution.csv"
        if solution.exists():
            entry["solution"] = f"solutions/{r['scenario']}_solution.csv"
        scenarios_out.append(entry)

    manifest = {
        "timestamp": datetime.now().isoformat(),
        "git_commit": get_git_commit(),
        "platform": platform.platform(),
        "dataset": dataset,
        "comparison_description": comparison_config.get("comparison", {}).get("description", ""),
        "total_wall_time_seconds": round(total_wall_seconds, 2),
        "scenario_count": len(results),
        "scenarios": scenarios_out,
    }
    manifest_path = run_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"Manifest written to: {manifest_path}")


def write_run_readme(run_dir: Path, dataset: str, comparison_config: dict,
                     results: list[dict]):
    """Write human-readable README.md summarizing the run."""
    description = comparison_config.get("comparison", {}).get("description", "CNDP Comparison")
    lines = [
        f"# {description} - {dataset}",
        "",
        f"- **Dataset**: {dataset}",
        f"- **Timestamp**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- **Git commit**: {get_git_commit()}",
        f"- **Scenarios**: {len(results)}",
        "",
        "## Results",
        "",
        "| Scenario | Status | Time (s) | Best Objective | Budget |",
        "|----------|--------|----------|----------------|--------|",
    ]
    for r in results:
        status = "OK" if r["exit_code"] == 0 else "FAILED"
        obj = f"{r['final_objective']:.4f}" if r["final_objective"] is not None else "N/A"
        budget = f"{r['final_budget']:.4f}" if r["final_budget"] is not None else "N/A"
        lines.append(f"| {r['scenario']} | {status} | {r['elapsed_seconds']:.1f} | {obj} | {budget} |")

    lines += [
        "",
        "## Folder Structure",
        "",
        "```",
        "figures/    - All plots (PNG + PDF)",
        "tables/     - Summary tables (CSV + LaTeX)",
        "traces/     - Per-scenario iteration-level trace data",
        "solutions/  - Per-scenario final link capacities",
        "metadata/   - Per-scenario run metadata JSON",
        "configs/    - Exact config files used (for reproducibility)",
        "```",
    ]
    (run_dir / "README.md").write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run CNDP comparison from a TOML manifest",
    )
    parser.add_argument(
        "comparison", type=Path,
        help="Path to TOML comparison manifest",
    )
    parser.add_argument(
        "--scenarios",
        help="Comma-separated scenario names to run",
    )
    parser.add_argument(
        "--exclude",
        help="Comma-separated scenario names to skip",
    )
    parser.add_argument(
        "--etalon-only", action="store_true",
        help="Run only etalon scenarios",
    )
    parser.add_argument(
        "--skip-etalon", action="store_true",
        help="Skip etalon scenarios",
    )
    parser.add_argument(
        "--repeat", type=int, default=1,
        help="Repeat non-etalon scenarios N times (default: 1)",
    )
    parser.add_argument(
        "--etalon-time-limit", type=int,
        help="Time limit for etalon scenarios in seconds (overrides TOML)",
    )
    parser.add_argument(
        "--scenario-time-limit", type=int,
        help="Time limit for non-etalon scenarios in seconds (overrides TOML)",
    )
    parser.add_argument(
        "--exe", type=str,
        help="Path to cndp_solver.exe (default: build/mingw-vcpkg-release/cndp_solver.exe)",
    )

    args = parser.parse_args()

    # Load and validate TOML
    config = load_comparison(args.comparison)
    comp = config["comparison"]
    defaults = config.get("defaults", {})

    # Resolve time limits: CLI overrides > TOML values > defaults
    etalon_time_limit = args.etalon_time_limit or comp.get("etalon_time_limit", 7200)
    scenario_time_limit = args.scenario_time_limit or comp.get("time_limit", 1800)

    # Resolve executable
    exe = resolve_exe(args.exe, DEFAULT_CNDP_EXE)

    # Filter and expand scenarios
    scenarios = filter_scenarios(config["scenarios"], args)
    scenarios = expand_repeats(scenarios, args.repeat)

    if not scenarios:
        print("No scenarios to run after filtering.")
        sys.exit(0)

    # Dataset comes from defaults.network.dataset
    dataset = defaults.get("network", {}).get("dataset")
    if not dataset:
        print("ERROR: Could not determine dataset. Set [defaults.network] dataset in manifest.")
        sys.exit(1)

    # Create run directory
    run_dir = create_run_dir(dataset, "cndp", "comparison")
    for sub in ["figures", "tables", "traces", "solutions", "metadata", "configs"]:
        (run_dir / sub).mkdir(parents=True, exist_ok=True)

    # Generate per-scenario TOML configs
    configs_dir = run_dir / "configs"
    for scenario in scenarios:
        toml_path = generate_scenario_toml(defaults, scenario, configs_dir)
        scenario["config_path"] = toml_path

    # Print header
    description = comp.get("description", "CNDP Comparison")
    print(f"\n{description}")
    print(f"Run folder: {run_dir}")
    print(f"Dataset: {dataset}")
    print(f"Executable: {exe}")
    print(f"Scenarios: {len(scenarios)}")
    for i, s in enumerate(scenarios):
        tag = " [ETALON]" if s.get("etalon") else ""
        print(f"  {i+1}. {s['name']}{tag} -- {s.get('description', '')}")

    # Run scenarios
    batch = BatchProgressRenderer()
    wall_start = time.time()
    results = []

    for i, scenario in enumerate(scenarios, 1):
        batch.start_scenario(i, len(scenarios), scenario)

        tl = etalon_time_limit if scenario.get("etalon") else scenario_time_limit
        result = run_scenario(scenario, run_dir, exe, tl)
        results.append(result)

        batch.finish_scenario(scenario, result)

        # Post-process outputs into organized subfolders
        post_process_cndp_outputs(run_dir, scenario["name"], scenario["config_path"])

    total_wall = time.time() - wall_start

    # Clean up root-level C++ output files
    cleanup_cndp_root(run_dir)

    # Copy comparison TOML for reproducibility
    shutil.copy2(str(args.comparison.resolve()), str(run_dir / "configs" / args.comparison.name))

    # Write outputs
    write_comparison_results(results, run_dir)
    write_run_manifest(run_dir, dataset, config, results, total_wall)
    write_run_readme(run_dir, dataset, config, results)

    # Generate comparison plots
    plot_script = PROJECT_ROOT / "scripts" / "plot_cndp_comparison.py"
    if plot_script.exists():
        print(f"\nGenerating comparison plots...")
        subprocess.run(
            [sys.executable, str(plot_script), "--run-dir", str(run_dir),
             "--dataset", dataset],
            cwd=str(PROJECT_ROOT),
        )

    # Print summary table
    print(f"\n{'='*90}")
    print(f"SUMMARY: {dataset} -- {description}")
    print(f"{'='*90}")
    print(f"{'Scenario':<30} {'Status':<10} {'Time(s)':<10} {'Objective':<20} {'Budget':<15}")
    print(f"{'-'*30} {'-'*10} {'-'*10} {'-'*20} {'-'*15}")
    for r in results:
        status = "OK" if r["exit_code"] == 0 else "FAILED"
        obj_str = f"{r['final_objective']:.6f}" if r["final_objective"] is not None else "N/A"
        budget_str = f"{r['final_budget']:.4f}" if r["final_budget"] is not None else "N/A"
        print(f"{r['scenario']:<30} {status:<10} {r['elapsed_seconds']:<10.1f} {obj_str:<20} {budget_str:<15}")

    print(f"\nTotal wall time: {total_wall:.1f}s")
    print(f"Self-contained run folder: {run_dir}")


if __name__ == "__main__":
    main()
