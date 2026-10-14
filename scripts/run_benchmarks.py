#!/usr/bin/env python3
"""Run reproducible TAP/scaling/CNDP measurements, or render saved measurements."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

# Set before importing NumPy/SciPy, including in fresh workers. Native route
# search uses its own explicit thread count and is unaffected by these pools.
NUMERICAL_THREAD_ENV = {name: "1" for name in (
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
)}
os.environ.update(NUMERICAL_THREAD_ENV)

import numpy as np

from benchmark_data import (
    DATA_REPOSITORY, DATA_REVISION, load_input, prepare_data, sha256, validate_assignment,
)

ROOT = Path(__file__).resolve().parents[1]
BUILD_SETTINGS = ("CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS",
                  "CMAKE_CXX_FLAGS_RELEASE")


def configurations(profile: str) -> list[dict]:
    datasets = ["SiouxFalls"] if profile == "quick" else ["SiouxFalls", "ChicagoSketch", "Anaheim", "Sydney"]
    tolerances = [1e-8] if profile == "quick" else [1e-6, 1e-8, 1e-10]
    cases = []
    for dataset in datasets:
        for tolerance in tolerances:
            for approach in ("tapas", "routebased"):
                cases.append(dict(
                    id=f"{dataset}_{approach}_{tolerance:.0e}_t1",
                    kind="tap", dataset=dataset, approach=approach, tolerance=tolerance,
                    threads=1, max_iterations=1000, repeats=1 if profile == "quick" else 5,
                    timeout=600 if dataset == "Sydney" else 120,
                ))
    if profile == "full":
        for dataset in ("ChicagoSketch", "Anaheim", "Sydney"):
            for threads in (2, 4, 8):
                cases.append(dict(
                    id=f"{dataset}_routebased_1e-08_t{threads}",
                    kind="tap", dataset=dataset, approach="routebased",
                    tolerance=1e-8, threads=threads, max_iterations=1000, repeats=5,
                    timeout=600 if dataset == "Sydney" else 120,
                ))
    cases.append(dict(
        id="SiouxFalls_cndp_cobyla", kind="cndp", dataset="SiouxFalls", approach="tapas",
        tolerance=1e-8, threads=1, max_iterations=1000,
        max_evaluations=30 if profile == "quick" else 500, optimizer_tolerance=1e-4,
        theta=5.0, budget=100000.0, capacity_multiplier=1.5,
        repeats=1 if profile == "quick" else 3, timeout=300,
    ))
    return cases


def schedule(cases: list[dict]):
    for case in cases:
        yield case, 0, True
    for repeat in range(1, max(case["repeats"] for case in cases) + 1):
        offset = (repeat - 1) % len(cases)
        for case in cases[offset:] + cases[:offset]:
            if repeat <= case["repeats"]:
                yield case, repeat, False


def clean_json(value):
    """Strict JSON represents unavailable/nonfinite measurements as null."""
    if isinstance(value, dict):
        return {k: clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path: Path, value):
    path.write_text(json.dumps(clean_json(value), indent=2, allow_nan=False) + "\n")


def worker(case: dict, cache: Path) -> dict:
    import traffic_assignment as ta

    data = load_input(cache, case["dataset"])
    network = data.build()
    options = dict(alpha=case["tolerance"], max_iterations=case["max_iterations"],
                   route_search_threads=case["threads"])
    if case["kind"] == "tap":
        result = ta.solve_tap(network, approach=case["approach"], **options)
        validation = validate_assignment(
            data, result.flows, result.link_costs, network.capacities(),
            result.total_travel_time, result.relative_gap, case["tolerance"],
        )
        return {
            "status": "ok" if validation["valid"] else "validation_failed",
            "solve_seconds": result.solve_seconds, "relative_gap": result.relative_gap,
            "beckmann_objective": result.beckmann_objective,
            "total_travel_time": result.total_travel_time, "validation": validation,
        }

    baseline = ta.solve_tap(network, approach=case["approach"], **options)
    baseline_validation = validate_assignment(
        data, baseline.flows, baseline.link_costs, network.capacities(),
        baseline.total_travel_time, baseline.relative_gap, case["tolerance"],
    )
    network = data.build()
    lower = network.capacities()
    upper = lower * case["capacity_multiplier"]
    result = ta.solve_cndp(
        network, [ta.step("nlopt", algorithm="LN_COBYLA",
                         max_iterations=case["max_evaluations"],
                         tolerance=case["optimizer_tolerance"])],
        constraints=ta.constraints_from_arrays(network, lower, upper),
        approach=case["approach"], approach_options=options,
        theta=case["theta"], budget=case["budget"],
    )
    final_validation = validate_assignment(
        data, result.flows, network.link_costs(), result.capacities,
        result.total_travel_time, network.relative_gap(), case["tolerance"],
    )
    # Validate the design independently of the optimizer's warm-started network.
    fresh_network = data.build(result.capacities)
    fresh = ta.solve_tap(fresh_network, approach=case["approach"], **options)
    fresh_validation = validate_assignment(
        data, fresh.flows, fresh.link_costs, result.capacities,
        fresh.total_travel_time, fresh.relative_gap, case["tolerance"],
    )
    investment = float(case["theta"] * np.sum(result.capacities - lower))
    bound_violation = float(max(0, np.max(lower - result.capacities),
                                np.max(result.capacities - upper)))
    budget_violation = max(0.0, investment - case["budget"])
    errors = []
    if not baseline_validation["valid"]:
        errors.append("baseline assignment")
    if not final_validation["valid"]:
        errors.append("final assignment")
    if not fresh_validation["valid"]:
        errors.append("fresh assignment")
    if bound_violation > 1e-7 or budget_violation > 1e-6:
        errors.append("design feasibility")
    if not math.isclose(investment, result.budget, rel_tol=1e-9, abs_tol=1e-6):
        errors.append("investment calculation")
    if not math.isclose(result.objective, result.total_travel_time + investment, rel_tol=1e-9):
        errors.append("objective decomposition")
    if not math.isclose(fresh.total_travel_time, result.total_travel_time, rel_tol=1e-5):
        errors.append("fresh/warm-started travel time agreement")
    return {
        "status": "ok" if not errors else "validation_failed",
        "solve_seconds": result.elapsed_seconds, "relative_gap": network.relative_gap(),
        "total_travel_time": result.total_travel_time, "objective": result.objective,
        "baseline_total_travel_time": baseline.total_travel_time,
        "fresh_total_travel_time": fresh.total_travel_time,
        "travel_time_reduction_percent": 100 * (1 - fresh.total_travel_time / baseline.total_travel_time),
        "investment": investment, "budget_violation": budget_violation,
        "bound_violation": bound_violation,
        "validation": {"valid": not errors, "errors": errors, "baseline": baseline_validation,
                       "final": final_validation, "fresh": fresh_validation},
        "solution": [
            dict(link_index=i, init_node=int(data.links[i, 0]), term_node=int(data.links[i, 1]),
                 original_capacity=float(lower[i]), designed_capacity=float(result.capacities[i]),
                 baseline_flow=float(baseline.flows[i]), designed_flow=float(fresh.flows[i]))
            for i in range(len(lower))
        ],
    }


def run_process(command: list[str], timeout: float) -> dict:
    start = time.monotonic()
    try:
        process = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "worker_wall_seconds": time.monotonic() - start}
    if process.returncode:
        return {"status": "error", "returncode": process.returncode,
                "error": process.stderr[-3000:] or process.stdout[-3000:],
                "worker_wall_seconds": time.monotonic() - start}
    try:
        result = json.loads(process.stdout)
        if not isinstance(result, dict) or "status" not in result:
            raise ValueError("Worker response has no status")
    except (ValueError, TypeError) as error:
        return {"status": "error", "error": f"Invalid worker response: {error}",
                "worker_wall_seconds": time.monotonic() - start}
    result["worker_wall_seconds"] = time.monotonic() - start
    return result


def collect(cases: list[dict], output: Path, cache: Path, time_budget: float,
            run=run_process) -> list[dict]:
    """Append and flush after each worker so interruption preserves prior results."""
    records = []
    deadline = time.monotonic() + time_budget
    with (output / "measurements.jsonl").open("x") as stream:
        for case, repeat, warmup in schedule(cases):
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                result = {"status": "budget_exhausted"}
            else:
                command = [sys.executable, str(Path(__file__).resolve()), "--worker",
                           json.dumps(case), "--cache", str(cache)]
                result = run(command, min(case["timeout"], remaining))
            record = dict(case_id=case["id"], repeat=repeat, warmup=warmup, **result)
            record = clean_json(record)
            stream.write(json.dumps(record, allow_nan=False) + "\n")
            stream.flush()
            records.append(record)
            label = "warm-up" if warmup else f"repeat {repeat}"
            seconds = record.get("solve_seconds")
            duration = f" {seconds:.4f}s" if seconds is not None else ""
            print(f"{case['id']} {label}: {record['status']}{duration}", flush=True)
    return records


def build_metadata() -> dict:
    import traffic_assignment._core as core

    cache = ROOT / "build/python/CMakeCache.txt"
    build = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            for key in BUILD_SETTINGS:
                if line.startswith(key + ":"):
                    build[key] = line.split("=", 1)[1]
    compiler = build.get("CMAKE_CXX_COMPILER", "c++")
    try:
        build["compiler_version"] = subprocess.check_output(
            [compiler, "--version"], text=True, timeout=10).splitlines()[0]
    except (OSError, subprocess.SubprocessError):
        build["compiler_version"] = "unavailable"
    build["extension_sha256"] = sha256(Path(core.__file__).read_bytes())
    return build


def provenance() -> dict:
    # Include untracked benchmark code while keeping data/results out of the hash.
    files = []
    for directory in ("cpp", "python", "scripts", "cmake"):
        files.extend(p for p in (ROOT / directory).rglob("*")
                     if p.is_file() and p.suffix in {".py", ".cpp", ".h", ".hpp", ".cmake"}
                     and "__pycache__" not in p.parts)
    files.extend(ROOT / p for p in ("CMakeLists.txt", "cpp/CMakeLists.txt",
                                    "CMakePresets.json", "pyproject.toml"))
    hashes = {str(p.relative_to(ROOT)): sha256(p.read_bytes()) for p in sorted(set(files))}
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        dirty = bool(subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=normal"], cwd=ROOT, text=True))
    except subprocess.CalledProcessError:
        revision, dirty = "unavailable", True
    cpu = platform.processor()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        cpu = next((line.split(":", 1)[1].strip() for line in cpuinfo.read_text().splitlines()
                    if line.startswith("model name")), cpu)
    packages = {}
    for name in ("traffic-assignment", "numpy", "scipy", "pandas", "matplotlib", "pybind11", "scikit-build-core"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "git_commit": revision, "working_tree_dirty": dirty,
        "source_sha256": sha256(json.dumps(hashes, sort_keys=True).encode()),
        "source_files": hashes, "platform": platform.platform(), "cpu": cpu,
        "logical_cpus": os.cpu_count(),
        "available_cpus": len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count(),
        "python": platform.python_version(), "packages": packages,
        "build": build_metadata(), "tapas_seed": 42,
        "numerical_thread_environment": NUMERICAL_THREAD_ENV,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("quick", "full"), default="quick")
    parser.add_argument("--source-root", type=Path,
                        help="TransportationNetworks Git clone containing the pinned revision")
    parser.add_argument("--output", type=Path, help="New run directory (existing runs are preserved)")
    parser.add_argument("--time-budget", type=float,
                        help="Worker budget including warm-ups and validation (full: 14400s; quick: 3600s)")
    parser.add_argument("--render-only", type=Path, metavar="RUN_DIR",
                        help="Regenerate summary, report, and figures without running solvers")
    parser.add_argument("--worker", help=argparse.SUPPRESS)
    parser.add_argument("--cache", type=Path, default=ROOT / "build/benchmark-data", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        result = worker(json.loads(args.worker), args.cache)
        print(json.dumps(clean_json(result), allow_nan=False))
        return 0
    from benchmark_report import render
    if args.render_only:
        render(args.render_only.resolve())
        return 0
    if args.time_budget is None:
        args.time_budget = 14400 if args.profile == "full" else 3600
    if not math.isfinite(args.time_budget) or args.time_budget <= 0:
        parser.error("--time-budget must be finite and positive")
    output = (args.output or ROOT / "performance_results/benchmarks" /
              datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")).resolve()
    if (output / "manifest.json").exists() or (output / "measurements.jsonl").exists():
        parser.error("Output already contains a run; choose a new directory or use --render-only")
    cases = configurations(args.profile)
    datasets = list(dict.fromkeys(c["dataset"] for c in cases))
    cache = args.cache.resolve()
    sources = prepare_data(cache, datasets,
                           args.source_root.resolve() if args.source_root else None)
    info = provenance()
    if info["build"].get("CMAKE_BUILD_TYPE") != "Release":
        parser.error("Benchmarks require a local Release build; install with -Ccmake.build-type=Release")
    output.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": 2, "started_utc": datetime.now(timezone.utc).isoformat(),
        "profile": args.profile, "time_budget_seconds": args.time_budget,
        "data_repository": DATA_REPOSITORY, "data_revision": DATA_REVISION, "data_files": sources,
        "datasets": {name: load_input(cache, name).dimensions() for name in datasets},
        "model": "BPR travel-time-only; no toll or distance cost; unrestricted transit through all nodes",
        "validation_backend": "SciPy directed Dijkstra; 64 origins per batch; minimum parallel-link cost",
        "timing": "Native computation; excludes loading, construction, validation and plotting",
        "configuration": cases, **info,
    }
    write_json(output / "manifest.json", manifest)
    started = time.monotonic()
    records = collect(cases, output, cache, args.time_budget)
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["measurement_wall_seconds"] = time.monotonic() - started
    manifest["complete"] = all(r["status"] != "budget_exhausted" for r in records)
    write_json(output / "manifest.json", manifest)
    render(output)
    print(f"Report: {output / 'README.md'}")
    return 0 if all(r["status"] == "ok" for r in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
