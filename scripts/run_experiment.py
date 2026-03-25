"""
Unified Experiment Runner for TAP and CNDP.

Creates a self-contained run folder with CSV outputs and publication-ready plots.

Usage:
    python scripts/run_experiment.py tap --dataset SiouxFalls
    python scripts/run_experiment.py tap --dataset SiouxFalls --approach RouteBased --shift-method NewtonStep
    python scripts/run_experiment.py cndp --config configs/cnd.siouxfalls.ini
    python scripts/run_experiment.py cndp --config configs/cnd.siouxfalls.ini --scenario-name MyTest
"""

import argparse
import json
import platform
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CNDP_EXE = PROJECT_ROOT / "build" / "mingw-vcpkg-release" / "cndp_solver.exe"
DEFAULT_TAP_EXE = PROJECT_ROOT / "build" / "mingw-vcpkg-release" / "tap_solver.exe"


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def get_git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, cwd=str(PROJECT_ROOT), timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def resolve_exe(cli_override: Optional[str], default: Path) -> Path:
    exe = Path(cli_override).resolve() if cli_override else default
    if not exe.exists():
        print(f"ERROR: Executable not found: {exe}")
        print("Build first: cmake --build --preset build-release")
        sys.exit(1)
    return exe


def read_dataset_from_config(config_path: Path) -> Optional[str]:
    """Read dataset name from [run] section of an INI config file."""
    in_run = False
    with open(config_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.lower() == "[run]":
                in_run = True
                continue
            if stripped.startswith("["):
                in_run = False
                continue
            if in_run and stripped.startswith("dataset"):
                _, _, value = stripped.partition("=")
                return value.strip()
    return None


def create_run_dir(dataset: str, mode: str, label: Optional[str] = None) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    suffix = label or mode
    run_dir = PROJECT_ROOT / "performance_results" / dataset / mode / "runs" / f"{timestamp}_{dataset}_{suffix}"
    return run_dir


def stream_subprocess(cmd: list[str], time_limit: Optional[int] = None) -> tuple[int, list[str]]:
    """Run subprocess, streaming stdout live. Returns (exit_code, stdout_lines)."""
    process = subprocess.Popen(
        cmd,
        cwd=str(PROJECT_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    stdout_lines = []
    for line in process.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()
        stdout_lines.append(line.rstrip())

    process.wait(timeout=time_limit)

    stderr_text = process.stderr.read() if process.stderr else ""
    if stderr_text.strip():
        sys.stderr.write(stderr_text)

    return process.returncode, stdout_lines


def write_manifest(run_dir: Path, dataset: str, mode: str, extra: dict):
    manifest = {
        "timestamp": datetime.now().isoformat(),
        "git_commit": get_git_commit(),
        "platform": platform.platform(),
        "dataset": dataset,
        "mode": mode,
        **extra,
    }
    (run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))


def write_readme(run_dir: Path, content: str):
    (run_dir / "README.md").write_text(content)


# ---------------------------------------------------------------------------
# TAP mode
# ---------------------------------------------------------------------------

def generate_tap_config(run_dir: Path, args) -> Path:
    """Generate a temporary INI config from CLI args."""
    lines = ["[run]"]
    lines.append(f"dataset = {args.dataset}")
    lines.append(f"approach = {args.approach}")
    if args.shift_method:
        lines.append(f"shift_method = {args.shift_method}")
    if args.approach_alpha is not None:
        lines.append(f"approach_alpha = {args.approach_alpha}")
    if args.max_iterations is not None:
        lines.append(f"max_standard_iterations = {args.max_iterations}")
    lines.append("")

    config_path = run_dir / "configs" / "generated_tap.ini"
    config_path.write_text("\n".join(lines))
    return config_path


def generate_tap_plot(run_dir: Path, dataset: str):
    """Generate convergence plot from TAP trace CSVs (replicates graphs_plotting.ipynb)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        print("WARNING: matplotlib/pandas not installed, skipping plot generation")
        return

    traces_dir = run_dir / "traces"
    csv_files = sorted(traces_dir.glob("*.csv"))
    if not csv_files:
        print("WARNING: No trace CSVs found, skipping plot")
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    for csv_path in csv_files:
        df = pd.read_csv(csv_path)
        label = csv_path.stem
        ax.plot(df["ElapsedTime(s)"], df["RelativeGap"], label=label)

    ax.set_yscale("log")
    ax.set_yticks([10 ** i for i in range(-18, 1, 3)])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Relative Gap")
    ax.set_title(f"TAP Convergence — {dataset}")
    ax.legend()
    ax.grid(True, which="major", axis="both")

    figures_dir = run_dir / "figures"
    for fmt in ("png", "pdf"):
        fig.savefig(figures_dir / f"convergence_{dataset}.{fmt}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to {figures_dir}/convergence_{dataset}.png")


def run_tap(args):
    exe = resolve_exe(args.exe, DEFAULT_TAP_EXE)

    # Resolve effective dataset: prefer INI config value over CLI default
    effective_dataset = args.dataset
    if args.config:
        config_path = Path(args.config).resolve()
        cfg_dataset = read_dataset_from_config(config_path)
        if cfg_dataset:
            effective_dataset = cfg_dataset

    run_dir = create_run_dir(effective_dataset, "tap", args.label)

    for sub in ["figures", "traces", "configs", "solutions"]:
        (run_dir / sub).mkdir(parents=True, exist_ok=True)

    # Config
    if args.config:
        shutil.copy2(str(config_path), str(run_dir / "configs" / config_path.name))
    else:
        config_path = generate_tap_config(run_dir, args)

    # Build command
    staging_dir = run_dir / "staging"
    staging_dir.mkdir(exist_ok=True)

    cmd = [str(exe)]
    if args.config:
        cmd += ["--config", str(config_path)]
    else:
        cmd += ["--config", str(config_path)]
    cmd += ["--output-root", str(staging_dir)]

    print(f"\nRun folder: {run_dir}")
    print(f"Executable: {exe}")
    print(f"Config: {config_path}")
    print(f"{'='*70}")

    start = time.time()
    exit_code, stdout_lines = stream_subprocess(cmd, args.time_limit)
    elapsed = time.time() - start

    # Move trace CSV from staging into traces/
    # TAP writes to {staging}/{dataset}/{ApproachName}.csv
    staging_dataset = staging_dir / effective_dataset
    if staging_dataset.exists():
        # Move link_flows.csv into solutions/
        link_flows_src = staging_dataset / "link_flows.csv"
        if link_flows_src.exists():
            shutil.move(str(link_flows_src), str(run_dir / "solutions" / "link_flows.csv"))

        # Move remaining CSVs (trace files) into traces/
        for csv_file in staging_dataset.glob("*.csv"):
            shutil.move(str(csv_file), str(run_dir / "traces" / csv_file.name))

    # Clean up staging
    shutil.rmtree(str(staging_dir), ignore_errors=True)

    # Parse final results from stdout
    final_results = {}
    for line in stdout_lines:
        for key in ("RelativeGap", "ObjectiveFunction", "TotalTravelTime"):
            if key + ":" in line:
                match = re.search(rf"{key}:\s+([0-9eE.+\-]+)", line)
                if match:
                    final_results[key] = match.group(1)

    print(f"\n{'='*70}")
    print(f"TAP completed in {elapsed:.1f}s (exit code: {exit_code})")
    for k, v in final_results.items():
        print(f"  {k}: {v}")

    # Generate convergence plot
    generate_tap_plot(run_dir, effective_dataset)

    # Write manifest & readme
    write_manifest(run_dir, effective_dataset, "tap", {
        "exe": str(exe),
        "exit_code": exit_code,
        "elapsed_seconds": round(elapsed, 2),
        "final_results": final_results,
    })

    readme_lines = [
        f"# TAP Experiment — {effective_dataset}",
        "",
        f"- **Dataset**: {effective_dataset}",
        f"- **Approach**: {args.approach}",
        f"- **Timestamp**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- **Git commit**: {get_git_commit()}",
        f"- **Elapsed**: {elapsed:.1f}s",
        f"- **Exit code**: {exit_code}",
        "",
        "## Results",
        "",
    ]
    for k, v in final_results.items():
        readme_lines.append(f"- **{k}**: {v}")
    readme_lines += [
        "",
        "## Folder Structure",
        "",
        "```",
        "figures/    - Convergence plots (PNG + PDF)",
        "traces/     - Per-approach TAP trace CSV",
        "solutions/  - Link flow distribution CSV",
        "configs/    - Config file used",
        "```",
    ]
    write_readme(run_dir, "\n".join(readme_lines))

    print(f"\nSelf-contained run folder: {run_dir}")


# ---------------------------------------------------------------------------
# CNDP mode
# ---------------------------------------------------------------------------

def post_process_cndp_outputs(run_dir: Path, scenario_name: str, config_path: Path):
    """Move/rename C++ output files into organized subfolders."""
    # Move trace CSV
    for f in sorted(run_dir.glob("BilevelCND_*_quality_time.csv")):
        shutil.move(str(f), str(run_dir / "traces" / f"{scenario_name}_trace.csv"))
        break

    # Move solution CSV
    for f in sorted(run_dir.glob("BilevelCND_*_solution.csv")):
        shutil.move(str(f), str(run_dir / "solutions" / f"{scenario_name}_solution.csv"))
        break

    # Move metadata JSON
    for f in sorted(run_dir.glob("BilevelCND_*_metadata.json")):
        shutil.move(str(f), str(run_dir / "metadata" / f"{scenario_name}_metadata.json"))
        break

    # Copy run_summary.csv (append-only)
    summary_src = run_dir / "BilevelCND_run_summary.csv"
    if summary_src.exists():
        shutil.copy2(str(summary_src), str(run_dir / "metadata" / "run_summary.csv"))

    # Copy config
    if config_path.exists():
        shutil.copy2(str(config_path), str(run_dir / "configs" / f"{scenario_name}.ini"))


def cleanup_cndp_root(run_dir: Path):
    """Remove leftover C++ output files from run_dir root."""
    for pattern in ["BilevelCND_*_quality_time.csv", "BilevelCND_*_solution.csv",
                     "BilevelCND_*_metadata.json", "BilevelCND_run_summary.csv"]:
        for f in run_dir.glob(pattern):
            f.unlink(missing_ok=True)


def run_cndp(args):
    exe = resolve_exe(args.exe, DEFAULT_CNDP_EXE)
    config_path = Path(args.config).resolve()

    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)

    # Resolve effective dataset: prefer INI config value over CLI default
    effective_dataset = args.dataset
    cfg_dataset = read_dataset_from_config(config_path)
    if cfg_dataset:
        effective_dataset = cfg_dataset

    # Derive scenario name
    scenario_name = args.scenario_name
    if not scenario_name:
        # e.g. "cnd.siouxfalls.ini" -> "cnd.siouxfalls"
        scenario_name = config_path.stem

    run_dir = create_run_dir(effective_dataset, "cndp", args.label)

    for sub in ["figures", "tables", "traces", "solutions", "metadata", "configs"]:
        (run_dir / sub).mkdir(parents=True, exist_ok=True)

    cmd = [
        str(exe),
        "--config", str(config_path),
        "--metrics_output_root", str(run_dir),
        "--metrics_no_dataset_subdir", "true",
        "--metrics_scenario_name", scenario_name,
    ]

    print(f"\nRun folder: {run_dir}")
    print(f"Executable: {exe}")
    print(f"Config: {config_path}")
    print(f"Scenario: {scenario_name}")
    print(f"{'='*70}")

    start = time.time()
    exit_code, stdout_lines = stream_subprocess(cmd, args.time_limit)
    elapsed = time.time() - start

    # Parse final results from stdout
    final_objective = None
    final_travel_time = None
    final_budget = None
    for line in stdout_lines:
        if "optimization_time" in line:
            parts = line.strip().split()
            for part in parts:
                if part.startswith("objective_function="):
                    final_objective = float(part.split("=")[1])
                elif part.startswith("total_travel_time="):
                    final_travel_time = float(part.split("=")[1])
                elif part.startswith("budget_function="):
                    final_budget = float(part.split("=")[1])

    print(f"\n{'='*70}")
    print(f"CNDP completed in {elapsed:.1f}s (exit code: {exit_code})")
    if final_objective is not None:
        print(f"  Objective: {final_objective:.6f}")
    if final_budget is not None:
        print(f"  Budget:    {final_budget:.2f}")

    # Post-process outputs
    post_process_cndp_outputs(run_dir, scenario_name, config_path)
    cleanup_cndp_root(run_dir)

    # Generate plots via plot_cndp_comparison.py
    plot_script = PROJECT_ROOT / "scripts" / "plot_cndp_comparison.py"
    if plot_script.exists():
        print(f"\nGenerating comparison plots...")
        subprocess.run(
            [sys.executable, str(plot_script), "--run-dir", str(run_dir),
             "--dataset", effective_dataset],
            cwd=str(PROJECT_ROOT),
        )

    # Write manifest & readme
    write_manifest(run_dir, effective_dataset, "cndp", {
        "exe": str(exe),
        "scenario_name": scenario_name,
        "config": str(config_path),
        "exit_code": exit_code,
        "elapsed_seconds": round(elapsed, 2),
        "final_objective": final_objective,
        "final_travel_time": final_travel_time,
        "final_budget": final_budget,
    })

    obj_str = f"{final_objective:.4f}" if final_objective is not None else "N/A"
    budget_str = f"{final_budget:.2f}" if final_budget is not None else "N/A"
    status = "OK" if exit_code == 0 else "FAILED"

    readme_lines = [
        f"# CNDP Experiment — {effective_dataset}",
        "",
        f"- **Dataset**: {effective_dataset}",
        f"- **Scenario**: {scenario_name}",
        f"- **Timestamp**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- **Git commit**: {get_git_commit()}",
        f"- **Elapsed**: {elapsed:.1f}s",
        f"- **Status**: {status}",
        "",
        "## Results",
        "",
        f"| Scenario | Status | Time (s) | Objective | Budget |",
        f"|----------|--------|----------|-----------|--------|",
        f"| {scenario_name} | {status} | {elapsed:.1f} | {obj_str} | {budget_str} |",
        "",
        "## Folder Structure",
        "",
        "```",
        "figures/    - All plots (PNG + PDF)",
        "tables/     - Summary tables (CSV + LaTeX)",
        "traces/     - Iteration-level trace data",
        "solutions/  - Final link capacities",
        "metadata/   - Run metadata JSON",
        "configs/    - Exact config file used",
        "```",
    ]
    write_readme(run_dir, "\n".join(readme_lines))

    print(f"\nSelf-contained run folder: {run_dir}")


# ---------------------------------------------------------------------------
# Main: argparse with subcommands
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Unified experiment runner for TAP and CNDP",
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # --- TAP subcommand ---
    tap_parser = subparsers.add_parser("tap", help="Run standalone Traffic Assignment")
    tap_parser.add_argument("--dataset", default="SiouxFalls",
                            help="Dataset name (default: SiouxFalls)")
    tap_parser.add_argument("--config", default=None,
                            help="INI config file (optional; CLI args used if omitted)")
    tap_parser.add_argument("--approach", default="Tapas",
                            help="RouteBased | Tapas (default: Tapas)")
    tap_parser.add_argument("--shift-method", default=None,
                            help="Shift method for RouteBased (default: NewtonStep)")
    tap_parser.add_argument("--approach-alpha", type=float, default=None,
                            help="Convergence threshold")
    tap_parser.add_argument("--max-iterations", type=int, default=None,
                            help="TAP iteration limit")
    tap_parser.add_argument("--exe", default=None,
                            help="Path to tap_solver.exe")
    tap_parser.add_argument("--time-limit", type=int, default=None,
                            help="Time limit in seconds")
    tap_parser.add_argument("--label", default=None,
                            help="Custom run folder suffix")

    # --- CNDP subcommand ---
    cndp_parser = subparsers.add_parser("cndp", help="Run CNDP optimization")
    cndp_parser.add_argument("--config", required=True,
                             help="INI config file with [step.N] sections")
    cndp_parser.add_argument("--dataset", default="SiouxFalls",
                             help="Dataset name (default: SiouxFalls)")
    cndp_parser.add_argument("--scenario-name", default=None,
                             help="Scenario label (default: derived from config filename)")
    cndp_parser.add_argument("--exe", default=None,
                             help="Path to cndp_solver.exe")
    cndp_parser.add_argument("--time-limit", type=int, default=None,
                             help="Time limit in seconds")
    cndp_parser.add_argument("--label", default=None,
                             help="Custom run folder suffix")

    args = parser.parse_args()

    if args.mode == "tap":
        run_tap(args)
    elif args.mode == "cndp":
        run_cndp(args)


if __name__ == "__main__":
    main()
