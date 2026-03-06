"""
CNDP Algorithm Comparison Runner

Runs multiple optimization scenarios for the Continuous Network Design Problem
and collects trace data for comparison analysis.

Each scenario is defined as a config file with [step.N] sections.
Scenario configs are auto-generated in configs/scenarios/.

Usage:
    python scripts/run_cndp_comparison.py --dataset SiouxFalls
    python scripts/run_cndp_comparison.py --dataset SiouxFalls --etalon-only
    python scripts/run_cndp_comparison.py --dataset SiouxFalls --skip-etalon
    python scripts/run_cndp_comparison.py --dataset SiouxFalls --scenarios 1,2,3
    python scripts/run_cndp_comparison.py --dataset SiouxFalls --tier sensitivity
    python scripts/run_cndp_comparison.py --dataset SiouxFalls --repeat 3
"""

import argparse
import csv
import math
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional


PROJECT_ROOT = Path(__file__).resolve().parent.parent
EXE_PATH = PROJECT_ROOT / "build" / "mingw-vcpkg-release" / "main.exe"
SCENARIO_CONFIG_DIR = PROJECT_ROOT / "configs" / "scenarios"


@dataclass
class StepConfig:
    type: str
    algorithm: str = ""
    max_iterations: int = 100
    tolerance: float = 1e-4
    population_size: int = 0  # optimlib-specific: 0 = auto


@dataclass
class Scenario:
    name: str
    description: str
    steps: list[StepConfig] = field(default_factory=list)
    is_etalon: bool = False
    skip_scaling: bool = False


DATASETS = {
    "SiouxFalls": {
        "base_config": "configs/cnd.siouxfalls.ini",
        "budget_multiplier": 5,
        "compute_scale": 1.0,
    },
    "Anaheim": {
        "base_config": "configs/cnd.anaheim.ini",
        "budget_multiplier": 2,
        "compute_scale": 0.5,
    },
}


def _scale_iterations(steps: list[StepConfig], scale: float) -> list[StepConfig]:
    """Scale iteration counts by a factor, rounding up, minimum 1."""
    if scale == 1.0:
        return steps
    return [
        StepConfig(
            type=s.type,
            algorithm=s.algorithm,
            max_iterations=max(1, math.ceil(s.max_iterations * scale)),
            tolerance=s.tolerance,
            population_size=s.population_size,
        )
        for s in steps
    ]


def build_core_scenarios() -> list[Scenario]:
    """Build the core comparison scenarios (dataset-independent step definitions)."""
    return [
        Scenario(
            name="Etalon",
            description="Reference solution (COBYLA 5000)",
            steps=[
                StepConfig(type="nlopt", algorithm="LN_COBYLA", max_iterations=5000, tolerance=1e-4),
            ],
            is_etalon=True,
            skip_scaling=True,
        ),
        Scenario(
            name="COBYLA-only",
            description="Pure NLopt local optimizer (COBYLA 500 iters)",
            steps=[
                StepConfig(type="nlopt", algorithm="LN_COBYLA", max_iterations=500, tolerance=1e-4),
            ],
        ),
        Scenario(
            name="BOBYQA-only",
            description="Pure NLopt local optimizer (BOBYQA 500 iters)",
            steps=[
                StepConfig(type="nlopt", algorithm="LN_BOBYQA", max_iterations=500, tolerance=1e-4),
            ],
        ),
        Scenario(
            name="OptCond-only",
            description="Pure optimality condition method (35 iters)",
            steps=[
                StepConfig(type="optimality_condition", max_iterations=35),
            ],
        ),
        Scenario(
            name="COBYLA-then-OptCond",
            description="COBYLA 200 + OptCond 10 (default combo)",
            steps=[
                StepConfig(type="nlopt", algorithm="LN_COBYLA", max_iterations=200, tolerance=1e-4),
                StepConfig(type="optimality_condition", max_iterations=10),
            ],
        ),
        Scenario(
            name="ISRES-then-OptCond",
            description="Global optimizer ISRES 200 + OptCond 10",
            steps=[
                StepConfig(type="nlopt", algorithm="GN_ISRES", max_iterations=200, tolerance=1e-4),
                StepConfig(type="optimality_condition", max_iterations=10),
            ],
        ),
        Scenario(
            name="BOBYQA-then-OptCond",
            description="BOBYQA 200 + OptCond 10",
            steps=[
                StepConfig(type="nlopt", algorithm="LN_BOBYQA", max_iterations=200, tolerance=1e-4),
                StepConfig(type="optimality_condition", max_iterations=10),
            ],
        ),
        Scenario(
            name="OptCond-then-COBYLA",
            description="Reversed: OptCond 15 then COBYLA 200",
            steps=[
                StepConfig(type="optimality_condition", max_iterations=15),
                StepConfig(type="nlopt", algorithm="LN_COBYLA", max_iterations=200, tolerance=1e-4),
            ],
        ),
        # OptimLib population-based metaheuristics
        Scenario(
            name="DE-only",
            description="OptimLib Differential Evolution (500 gens)",
            steps=[
                StepConfig(type="optimlib", algorithm="DE", max_iterations=500),
            ],
        ),
        Scenario(
            name="PSO-only",
            description="OptimLib Particle Swarm Optimization (500 gens)",
            steps=[
                StepConfig(type="optimlib", algorithm="PSO", max_iterations=500),
            ],
        ),
        Scenario(
            name="DE-PRMM-only",
            description="OptimLib DE with Parameter Rotation (500 gens)",
            steps=[
                StepConfig(type="optimlib", algorithm="DE_PRMM", max_iterations=500),
            ],
        ),
        Scenario(
            name="NM-only",
            description="OptimLib Nelder-Mead (500 iters)",
            steps=[
                StepConfig(type="optimlib", algorithm="NM", max_iterations=500),
            ],
        ),
        Scenario(
            name="DE-then-OptCond",
            description="DE 200 gens + OptCond 10",
            steps=[
                StepConfig(type="optimlib", algorithm="DE", max_iterations=200),
                StepConfig(type="optimality_condition", max_iterations=10),
            ],
        ),
        Scenario(
            name="PSO-then-OptCond",
            description="PSO 200 gens + OptCond 10",
            steps=[
                StepConfig(type="optimlib", algorithm="PSO", max_iterations=200),
                StepConfig(type="optimality_condition", max_iterations=10),
            ],
        ),
    ]


def build_sensitivity_scenarios() -> list[Scenario]:
    """Build sensitivity sweep scenarios for detailed parameter analysis."""
    scenarios = []

    # OptCond iteration sweep
    for n in [5, 10, 20, 35, 50, 100]:
        scenarios.append(Scenario(
            name=f"OptCond-sweep-{n}",
            description=f"OptCond-only with {n} iterations",
            steps=[StepConfig(type="optimality_condition", max_iterations=n)],
        ))

    # COBYLA iteration sweep
    for n in [50, 100, 200, 500, 1000]:
        scenarios.append(Scenario(
            name=f"COBYLA-sweep-{n}",
            description=f"COBYLA-only with {n} iterations",
            steps=[StepConfig(type="nlopt", algorithm="LN_COBYLA", max_iterations=n, tolerance=1e-4)],
        ))

    # Hybrid ratio sweep: COBYLA-then-OptCond with varying splits
    for nlopt_n, oc_n in [(50, 30), (100, 20), (200, 15), (300, 10), (500, 5)]:
        scenarios.append(Scenario(
            name=f"Hybrid-{nlopt_n}c-{oc_n}oc",
            description=f"COBYLA {nlopt_n} + OptCond {oc_n}",
            steps=[
                StepConfig(type="nlopt", algorithm="LN_COBYLA", max_iterations=nlopt_n, tolerance=1e-4),
                StepConfig(type="optimality_condition", max_iterations=oc_n),
            ],
        ))

    return scenarios


def build_scenarios(tier: str = "core") -> list[Scenario]:
    """Build scenario list based on tier."""
    core = build_core_scenarios()
    if tier == "core":
        return core
    sensitivity = build_sensitivity_scenarios()
    if tier == "sensitivity":
        return sensitivity
    # tier == "all"
    return core + sensitivity


def read_base_config_sections(dataset: str) -> tuple[str, str]:
    """Read [run] and [metrics] sections from the base config, returning them as raw text."""
    base_path = PROJECT_ROOT / DATASETS[dataset]["base_config"]
    run_lines = []
    metrics_lines = []
    current = None

    with open(base_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.lower().startswith("[run]"):
                current = "run"
                continue
            elif stripped.lower().startswith("[metrics]"):
                current = "metrics"
                continue
            elif stripped.startswith("["):
                current = None
                continue

            if current == "run":
                run_lines.append(line.rstrip())
            elif current == "metrics":
                metrics_lines.append(line.rstrip())

    return "\n".join(run_lines), "\n".join(metrics_lines)


def apply_dataset_scaling(scenarios: list[Scenario], dataset: str) -> list[Scenario]:
    """Apply dataset-specific iteration scaling to scenarios."""
    scale = DATASETS[dataset]["compute_scale"]
    if scale == 1.0:
        return scenarios
    result = []
    for s in scenarios:
        if s.skip_scaling:
            result.append(s)
        else:
            result.append(Scenario(
                name=s.name,
                description=s.description,
                steps=_scale_iterations(s.steps, scale),
                is_etalon=s.is_etalon,
                skip_scaling=s.skip_scaling,
            ))
    return result


def generate_scenario_config(
    dataset: str,
    scenario: Scenario,
    run_section: str,
    metrics_section: str,
) -> Path:
    """Generate an INI config file for a scenario."""
    SCENARIO_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    config_path = SCENARIO_CONFIG_DIR / f"{dataset.lower()}.{scenario.name.lower().replace(' ', '_')}.ini"

    lines = []
    lines.append("[run]")
    lines.append(run_section)
    lines.append("")
    lines.append("[metrics]")
    lines.append(metrics_section)
    # Override flush for comparison runs
    lines.append("flush_every_n_points = 10")
    lines.append("")

    for i, step in enumerate(scenario.steps):
        lines.append(f"[step.{i}]")
        lines.append(f"type = {step.type}")
        if step.algorithm:
            lines.append(f"algorithm = {step.algorithm}")
        lines.append(f"max_iterations = {step.max_iterations}")
        if step.type == "nlopt":
            lines.append(f"tolerance = {step.tolerance}")
        if step.population_size > 0:
            lines.append(f"population_size = {step.population_size}")
        lines.append("")

    config_path.write_text("\n".join(lines))
    return config_path


def run_scenario(
    dataset: str,
    scenario: Scenario,
    config_path: Path,
    time_limit: Optional[int] = None,
) -> dict:
    """Run a single scenario and return results."""
    cmd = [
        str(EXE_PATH),
        "--config", str(config_path),
        "--metrics_scenario_name", scenario.name,
    ]

    print(f"\n{'='*70}")
    print(f"Running: {scenario.name}")
    print(f"Description: {scenario.description}")
    print(f"Config: {config_path}")
    print(f"{'='*70}")

    start_time = time.time()

    try:
        # Stream output live so user can see progress bars and metrics
        process = subprocess.Popen(
            cmd,
            cwd=str(PROJECT_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        stdout_lines = []
        # Read stdout line by line, printing live and collecting for parsing
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            stdout_lines.append(line.rstrip())

        process.wait(timeout=time_limit)
        elapsed = time.time() - start_time

        stderr_text = process.stderr.read() if process.stderr else ""
        if stderr_text.strip():
            sys.stderr.write(stderr_text)

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

        return {
            "scenario": scenario.name,
            "dataset": dataset,
            "description": scenario.description,
            "exit_code": process.returncode,
            "elapsed_seconds": elapsed,
            "final_objective": final_objective,
            "final_travel_time": final_travel_time,
            "final_budget": final_budget,
            "stdout_tail": "\n".join(stdout_lines[-10:]),
            "stderr_tail": "\n".join(stderr_text.splitlines()[-5:]) if stderr_text else "",
        }

    except (subprocess.TimeoutExpired, Exception) as e:
        elapsed = time.time() - start_time
        if isinstance(e, subprocess.TimeoutExpired):
            print(f"\n  TIMEOUT after {elapsed:.1f}s")
            process.kill()
        else:
            print(f"\n  ERROR: {e}")
        return {
            "scenario": scenario.name,
            "dataset": dataset,
            "description": scenario.description,
            "exit_code": -1,
            "elapsed_seconds": elapsed,
            "final_objective": None,
            "final_travel_time": None,
            "final_budget": None,
            "stdout_tail": "TIMEOUT" if isinstance(e, subprocess.TimeoutExpired) else str(e),
            "stderr_tail": "",
        }


def write_comparison_results(results: list[dict], dataset: str):
    """Write comparison results to CSV."""
    output_dir = PROJECT_ROOT / "performance_results" / dataset
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = output_dir / f"comparison_runs_{timestamp}.csv"

    fieldnames = [
        "scenario", "dataset", "description", "exit_code",
        "elapsed_seconds", "final_objective", "final_travel_time", "final_budget",
    ]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    print(f"\nComparison results written to: {output_path}")
    return output_path


def main():
    parser = argparse.ArgumentParser(description="Run CNDP comparison scenarios")
    parser.add_argument(
        "--dataset",
        choices=list(DATASETS.keys()),
        default="SiouxFalls",
        help="Dataset to run (default: SiouxFalls)",
    )
    parser.add_argument(
        "--tier",
        choices=["core", "sensitivity", "all"],
        default="core",
        help="Scenario tier: core (default 8 scenarios), sensitivity (sweeps), all",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="Run each scenario N times with suffixed names (default: 1)",
    )
    parser.add_argument(
        "--etalon-only",
        action="store_true",
        help="Run only the etalon (reference) scenario",
    )
    parser.add_argument(
        "--skip-etalon",
        action="store_true",
        help="Skip the etalon scenario, run comparison scenarios only",
    )
    parser.add_argument(
        "--scenarios",
        type=str,
        default=None,
        help="Comma-separated list of scenario indices (1-based) to run, e.g. '1,2,3'",
    )
    parser.add_argument(
        "--etalon-time-limit",
        type=int,
        default=7200,
        help="Time limit for etalon in seconds (default: 7200)",
    )
    parser.add_argument(
        "--scenario-time-limit",
        type=int,
        default=1800,
        help="Time limit for each comparison scenario in seconds (default: 1800)",
    )
    parser.add_argument(
        "--exe",
        type=str,
        default=None,
        help="Path to main.exe (default: build/mingw-vcpkg-release/main.exe)",
    )

    args = parser.parse_args()

    global EXE_PATH
    if args.exe:
        EXE_PATH = Path(args.exe).resolve()

    if not EXE_PATH.exists():
        print(f"ERROR: Executable not found: {EXE_PATH}")
        print("Build first: cmake --build --preset build-release")
        sys.exit(1)

    all_scenarios = build_scenarios(tier=args.tier)

    # Apply dataset-specific iteration scaling
    all_scenarios = apply_dataset_scaling(all_scenarios, args.dataset)

    # Filter scenarios
    scenarios_to_run = []
    if args.etalon_only:
        scenarios_to_run = [s for s in all_scenarios if s.is_etalon]
    elif args.skip_etalon:
        scenarios_to_run = [s for s in all_scenarios if not s.is_etalon]
    elif args.scenarios:
        indices = [int(x.strip()) for x in args.scenarios.split(",")]
        non_etalon = [s for s in all_scenarios if not s.is_etalon]
        for idx in indices:
            if 1 <= idx <= len(non_etalon):
                scenarios_to_run.append(non_etalon[idx - 1])
            else:
                print(f"WARNING: Invalid scenario index {idx}, skipping")
    else:
        scenarios_to_run = all_scenarios

    if not scenarios_to_run:
        print("No scenarios to run.")
        sys.exit(0)

    # Read base config sections
    run_section, metrics_section = read_base_config_sections(args.dataset)

    # Expand repeats: for --repeat N>1, create N copies with suffixed names
    expanded_scenarios = []
    for s in scenarios_to_run:
        if args.repeat > 1 and not s.is_etalon:
            for run_idx in range(1, args.repeat + 1):
                expanded_scenarios.append(Scenario(
                    name=f"{s.name}_run{run_idx}",
                    description=f"{s.description} (run {run_idx})",
                    steps=s.steps,
                    is_etalon=False,
                    skip_scaling=True,  # already scaled
                ))
        else:
            expanded_scenarios.append(s)
    scenarios_to_run = expanded_scenarios

    # Generate scenario configs
    print(f"Dataset: {args.dataset} (compute_scale: {DATASETS[args.dataset]['compute_scale']})")
    print(f"Tier: {args.tier}")
    print(f"Scenarios to run: {len(scenarios_to_run)}")
    scenario_configs = {}
    for i, s in enumerate(scenarios_to_run):
        config_path = generate_scenario_config(args.dataset, s, run_section, metrics_section)
        scenario_configs[s.name] = config_path
        tag = " [ETALON]" if s.is_etalon else ""
        print(f"  {i+1}. {s.name}{tag} - {s.description}")
        print(f"     Config: {config_path}")

    results = []
    for scenario in scenarios_to_run:
        time_limit = args.etalon_time_limit if scenario.is_etalon else args.scenario_time_limit
        config_path = scenario_configs[scenario.name]
        result = run_scenario(args.dataset, scenario, config_path, time_limit=time_limit)
        results.append(result)

        status = "OK" if result["exit_code"] == 0 else f"FAILED (exit={result['exit_code']})"
        obj_str = f"{result['final_objective']:.6f}" if result["final_objective"] else "N/A"
        print(f"\n  Result: {status}")
        print(f"  Elapsed: {result['elapsed_seconds']:.1f}s")
        print(f"  Objective: {obj_str}")
        if result["stderr_tail"]:
            print(f"  Stderr: {result['stderr_tail']}")

    write_comparison_results(results, args.dataset)

    # Print summary table
    print(f"\n{'='*90}")
    print(f"SUMMARY: {args.dataset}")
    print(f"{'='*90}")
    print(f"{'Scenario':<30} {'Status':<10} {'Time(s)':<10} {'Objective':<20} {'Budget':<15}")
    print(f"{'-'*30} {'-'*10} {'-'*10} {'-'*20} {'-'*15}")
    for r in results:
        status = "OK" if r["exit_code"] == 0 else "FAILED"
        obj_str = f"{r['final_objective']:.6f}" if r["final_objective"] else "N/A"
        budget_str = f"{r['final_budget']:.2f}" if r["final_budget"] else "N/A"
        print(f"{r['scenario']:<30} {status:<10} {r['elapsed_seconds']:<10.1f} {obj_str:<20} {budget_str:<15}")


if __name__ == "__main__":
    main()
