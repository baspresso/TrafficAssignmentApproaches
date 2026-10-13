"""Summarize saved benchmark measurements and regenerate publication figures."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np


def summarize(cases: list[dict], records: list[dict]) -> list[dict]:
    summary = []
    for case in cases:
        rows = [r for r in records if r["case_id"] == case["id"] and not r["warmup"]]
        good = [r for r in rows if r["status"] == "ok"]
        row = {k: case[k] for k in ("id", "kind", "dataset", "approach", "tolerance", "threads")}
        row.update(expected_runs=case["repeats"], recorded_runs=len(rows), successful_runs=len(good),
                   all_valid=len(good) == case["repeats"] == len(rows))
        for field in ("solve_seconds", "total_travel_time", "beckmann_objective", "objective",
                      "baseline_total_travel_time", "fresh_total_travel_time",
                      "travel_time_reduction_percent", "investment"):
            values = [r[field] for r in good if r.get(field) is not None]
            row[field] = float(np.median(values)) if values else None
            if field == "solve_seconds":
                row["q25_seconds"] = float(np.percentile(values, 25)) if values else None
                row["q75_seconds"] = float(np.percentile(values, 75)) if values else None
        row["max_relative_gap"] = max((abs(r["relative_gap"]) for r in good), default=None)
        row["max_budget_violation"] = max((r.get("budget_violation", 0) for r in good), default=None)
        row["max_bound_violation"] = max((r.get("bound_violation", 0) for r in good), default=None)
        summary.append(row)
    scaling = [s for s in summary if s["kind"] == "tap" and s["dataset"] == "ChicagoSketch"
               and s["approach"] == "routebased" and s["tolerance"] == 1e-8]
    reference = next((s for s in scaling if s["threads"] == 1), None)
    for row in scaling:
        consistent = bool(reference and reference["all_valid"] and row["all_valid"])
        if consistent:
            ref_records = [r for r in records if r["case_id"] == reference["id"]
                           and not r["warmup"] and r["status"] == "ok"]
            run_records = [r for r in records if r["case_id"] == row["id"]
                           and not r["warmup"] and r["status"] == "ok"]
            consistent = all(
                math.isclose(r[field], ref[field], rel_tol=1e-5, abs_tol=1e-7)
                for r in run_records for ref in ref_records
                for field in ("total_travel_time", "beckmann_objective")
            )
        row["agrees_with_serial"] = consistent
        row["speedup"] = reference["solve_seconds"] / row["solve_seconds"] if consistent else None
    return summary


def number(value, spec=".4f") -> str:
    return format(value, spec) if value is not None else "—"


def write_csv(path: Path, rows: list[dict]):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def plots(output: Path, summary: list[dict]):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ta-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "traffic-assignment-benchmarks"})
    directory = output / "figures"
    directory.mkdir(exist_ok=True)
    generated = []

    def save(fig, name):
        for extension in ("svg", "png"):
            fig.savefig(directory / f"{name}.{extension}", bbox_inches="tight", dpi=180,
                        metadata={"Date": None} if extension == "svg" else {})
        plt.close(fig)
        generated.append(name)

    tap = [s for s in summary if s["kind"] == "tap" and s["threads"] == 1]
    datasets = list(dict.fromkeys(s["dataset"] for s in tap))
    if datasets:
        fig, axes = plt.subplots(1, len(datasets), figsize=(6 * len(datasets), 4.3), squeeze=False)
        for ax, dataset in zip(axes[0], datasets):
            for approach, label, color in (("tapas", "TAPAS", "#2166ac"),
                                            ("routebased", "Route-based", "#d6604d")):
                rows = sorted((s for s in tap if s["dataset"] == dataset and
                               s["approach"] == approach and s["all_valid"]),
                              key=lambda s: s["tolerance"], reverse=True)
                if not rows:
                    continue
                times = np.array([s["solve_seconds"] for s in rows])
                error = np.array([times - [s["q25_seconds"] for s in rows],
                                  [s["q75_seconds"] for s in rows] - times])
                ax.errorbar([s["tolerance"] for s in rows], times, yerr=error,
                            fmt="o-", capsize=4, color=color, label=label)
            ax.set_xscale("log")
            ax.invert_xaxis()
            ax.set_yscale("log")
            ax.set(xlabel="Requested relative gap", ylabel="Native solve time (s)",
                   title=dataset + ("\nTravel-time-only" if dataset == "ChicagoSketch" else ""))
            ax.grid(alpha=0.2)
            if ax.lines:
                ax.legend()
        fig.suptitle("Assignment runtime at matched accuracy targets\nMedian and interquartile range; validated configurations", fontsize=12)
        fig.tight_layout()
        save(fig, "assignment_runtime")

    scaling = sorted((s for s in summary if "speedup" in s), key=lambda s: s["threads"])
    if scaling:
        fig, ax = plt.subplots(figsize=(7, 4.4))
        valid = [s for s in scaling if s["speedup"] is not None]
        if valid:
            ax.plot([s["threads"] for s in valid], [s["speedup"] for s in valid],
                    "o-", color="#2166ac", label="Measured whole-solve speedup")
            for s in valid:
                ax.annotate(f"{s['speedup']:.2f}×", (s["threads"], s["speedup"]),
                            xytext=(0, 9), textcoords="offset points", ha="center")
        ax.axhline(1, color="#777777", linestyle="--", linewidth=1, label="Serial baseline")
        ax.set(xticks=[s["threads"] for s in scaling], xlabel="Route-search threads",
               ylabel="Serial median / parallel median",
               title="Chicago Sketch: parallel route search\nTravel-time-only; relative-gap target 1e-8")
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.2)
        ax.legend(loc="best")
        fig.tight_layout()
        save(fig, "parallel_scaling")

    cndp = next((s for s in summary if s["kind"] == "cndp" and s["all_valid"]), None)
    if cndp:
        fig, ax = plt.subplots(figsize=(7, 4.4))
        values = [cndp["baseline_total_travel_time"], cndp["fresh_total_travel_time"]]
        bars = ax.bar(["Original capacities", "Designed capacities"],
                      np.array(values) / 1e6, color=["#8c9eaf", "#2166ac"], width=0.55)
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.08,
                    f"{value / 1e6:.3f}", ha="center")
        ax.set(ylabel="Total travel time (million native flow × time units)",
               title=f"Sioux Falls capacity optimization\n{cndp['travel_time_reduction_percent']:.2f}% simulated travel-time reduction")
        ax.set_ylim(0, max(values) / 1e6 * 1.2)
        ax.grid(axis="y", alpha=0.2)
        fig.tight_layout()
        save(fig, "capacity_optimization")
    return generated


def render(output: Path):
    manifest = json.loads((output / "manifest.json").read_text())
    records = [json.loads(line) for line in (output / "measurements.jsonl").read_text().splitlines()
               if line.strip()]
    summary = summarize(manifest["configuration"], records)
    write_csv(output / "summary.csv", summary)
    flat = [{k: v for k, v in r.items() if not isinstance(v, (list, dict))} for r in records]
    write_csv(output / "measurements.csv", flat)
    figures = plots(output, summary)
    lines = [
        "# Traffic assignment benchmark study", "",
        f"Measured {manifest['started_utc'][:10]} · **{manifest['profile']} profile** · "
        f"{sum(r['status'] == 'ok' for r in records if not r['warmup'])}/"
        f"{sum(not r['warmup'] for r in records)} measured runs passed validation.", "",
        "This study compares the repository's implementations under identical inputs and "
        "requested accuracy. Results apply to this machine and build. "
        "External solver comparisons and global CNDP optimality are outside this study.", "",
        "## Inputs and model", "",
        "| Dataset | Nodes | Links | Interzonal OD pairs | Source demand | Road-network demand |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for name, d in manifest["datasets"].items():
        lines.append(f"| {name} | {d['nodes']} | {d['links']} | {d['interzonal_od_pairs']} | "
                     f"{d['total_demand']:,.2f} | {d['interzonal_demand']:,.2f} |")
    lines.extend([
        "", "Source: [Transportation Networks for Research]"
        f"({manifest['data_repository']}/tree/{manifest['data_revision']}). "
        "Exact revision and input checksums are in [manifest.json](manifest.json). "
        "Fractional demand is preserved. Intrazonal trips remain in the input as zero-link "
        "trips; they contribute no road flow or travel time. Road-network demand excludes them.", "",
        "The selected networks use BPR travel-time costs and allow transit through every node. "
        "**Chicago Sketch is a travel-time-only variant**: toll and distance terms from "
        "the published generalized-cost reference are excluded, so its published reference "
        "objective does not apply here. Sioux Falls is a small research/debugging "
        "benchmark. Travel-time totals retain native dataset units.", "",
        "## Assignment accuracy and runtime", "",
        "Times are medians with 25th–75th percentiles, in seconds. Statistics use successful "
        "measured runs; the count exposes incomplete configurations. Warm-ups are excluded.", "",
        "| Dataset | Method | Target gap | Successful runs | Median [Q25, Q75] (s) | Max achieved gap |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ])
    for s in summary:
        if s["kind"] == "tap" and s["threads"] == 1:
            lines.append(
                f"| {s['dataset']} | {s['approach']} | {s['tolerance']:.0e} | "
                f"{s['successful_runs']}/{s['expected_runs']} | {number(s['solve_seconds'])} "
                f"[{number(s['q25_seconds'])}, {number(s['q75_seconds'])}] | "
                f"{number(s['max_relative_gap'], '.2e')} |")
    if "assignment_runtime" in figures:
        lines.extend(["", "![Runtime at matched accuracy](figures/assignment_runtime.svg)"])
    scaling = sorted((s for s in summary if "speedup" in s), key=lambda s: s["threads"])
    if scaling:
        lines.extend(["", "## Parallel route-search scaling", "",
                      "Chicago Sketch, route-based solver, target relative gap 1e-8. "
                      "Speedup measures the entire solve: flow redistribution remains sequential. "
                      "Speedups below 1 indicate a slowdown. Only fully validated configurations "
                      "with matching travel-time and Beckmann objectives receive a speedup.", "",
                      "| Threads | Successful runs | Median [Q25, Q75] (s) | Speedup |",
                      "| ---: | ---: | ---: | ---: |"])
        for s in scaling:
            lines.append(f"| {s['threads']} | {s['successful_runs']}/{s['expected_runs']} | "
                         f"{number(s['solve_seconds'])} [{number(s['q25_seconds'])}, "
                         f"{number(s['q75_seconds'])}] | {number(s['speedup'], '.2f')}× |")
        lines.extend(["", "![Parallel scaling](figures/parallel_scaling.svg)"])
    cndp = next(s for s in summary if s["kind"] == "cndp")
    config = next(s for s in manifest["configuration"] if s["kind"] == "cndp")
    lines.extend([
        "", "## Capacity-optimization case", "",
        f"Sioux Falls with COBYLA, at most {config['max_evaluations']} objective evaluations, "
        "optimizer tolerance 1e-4, and TAPAS at relative gap 1e-8. "
        "Original capacities are lower bounds; upper bounds are 1.5 times those values. "
        "The illustrative investment model uses theta = 5 and budget = 100000. "
        "The objective combines equilibrium travel time and investment.", "",
        f"Validated measured runs: **{cndp['successful_runs']}/{cndp['expected_runs']}**. "
        "The designed travel time below comes from a fresh equilibrium solve at the final capacities.", "",
        "| Measurement | Median |", "| --- | ---: |",
        f"| Baseline total travel time | {number(cndp['baseline_total_travel_time'], ',.3f')} |",
        f"| Designed total travel time | {number(cndp['fresh_total_travel_time'], ',.3f')} |",
        f"| Simulated travel-time reduction | {number(cndp['travel_time_reduction_percent'], '.2f')}% |",
        f"| Investment used | {number(cndp['investment'], ',.3f')} / 100,000 |",
        f"| Upper-level objective | {number(cndp['objective'], ',.3f')} |",
        f"| Optimization time | {number(cndp['solve_seconds'])} s |",
        f"| Maximum budget violation across valid runs | {number(cndp['max_budget_violation'], '.2e')} |",
        f"| Maximum capacity-bound violation across valid runs | {number(cndp['max_bound_violation'], '.2e')} |",
    ])
    if "capacity_optimization" in figures:
        lines.extend(["", "![Capacity optimization](figures/capacity_optimization.svg)"])
    problems = [r for r in records if r["status"] != "ok"]
    lines.extend(["", "## Methodology and validation", "",
                  "- Each configuration receives one unmeasured warm-up. Full-profile TAP "
                  "configurations have five measured repetitions; CNDP has three. Quick-profile "
                  "configurations have one. Configuration order rotates across repetitions.",
                  "- Workers execute sequentially in fresh processes and construct fresh networks. "
                  "Native computation time excludes imports, input loading, solver construction, "
                  "result validation, and plotting. CNDP time includes its lower-level solves "
                  "but excludes the separate baseline and fresh final-design validation.",
                  "- Statistics recording is disabled. TAPAS uses its existing seed of 42. "
                  "Repeated runs measure timing variability, rather than variation across seeds.",
                  "- Independent Python checks cover finite/nonnegative flows, node conservation, "
                  "BPR link costs, total travel time, and a Dijkstra-based equilibrium gap. "
                  "Gap validation allows an absolute 1e-12 roundoff margin.",
                  "- CNDP also checks capacity bounds, investment accounting, and a fresh final "
                  "equilibrium. Budget tolerance is 1e-6 and bound tolerance is 1e-7. "
                  "Fresh and warm-started travel times must agree within 1e-5 relative tolerance.",
                  "- TAP has a 1000-iteration ceiling. Worker limits are 120 seconds for TAP "
                  "and 300 seconds for CNDP, including validation. The full run has a "
                  "3600-second default deadline. Errors, timeouts, and unrun cases remain in raw results.",
                  "- Measurements follow the capacity-dependent derivative-cache correction. "
                  "CNDP results describe simulated outcomes under the specified model.", "",
                  "## Machine and build", "",
                  f"- CPU: {manifest['cpu']}; {manifest['available_cpus']} CPUs available to the process.",
                  f"- Platform: {manifest['platform']}.",
                  f"- Compiler: {manifest['build'].get('compiler_version', 'unavailable')}.",
                  f"- Build: {manifest['build'].get('CMAKE_BUILD_TYPE')}; "
                  f"release flags: {manifest['build'].get('CMAKE_CXX_FLAGS_RELEASE', '')}.",
                  f"- Python: {manifest['python']}; NumPy: {manifest['packages'].get('numpy')}.",
                  f"- Git base: {manifest['git_commit']}; uncommitted changes present: "
                  f"{str(manifest['working_tree_dirty']).lower()}.",
                  f"- Source fingerprint: {manifest['source_sha256']}.",
                  "- The manifest records individual source hashes, the native extension hash, "
                  "dependency versions, and exact configurations. Runtime depends on hardware, "
                  "system load, compiler, and numerical stopping behavior.", "",
                  "## Reproduce", "",
                  "From the repository root, after installing the native prerequisites in the main README:", "",
                  chr(96) * 3 + "bash",
                  "python -m pip install -e '.[examples,test]' -Ccmake.build-type=Release",
                  f"python scripts/run_benchmarks.py --profile {manifest['profile']} \\",
                  "  --time-budget 3600 --output performance_results/benchmark-reproduction",
                  "python scripts/run_benchmarks.py \\",
                  "  --render-only performance_results/benchmark-reproduction",
                  chr(96) * 3, "",
                  "The runner downloads only the pinned TNTP inputs into build/benchmark-data. "
                  "For offline preparation, add --source-root data/TransportationNetworks "
                  "with a Git clone containing the pinned revision; it reads committed blobs "
                  "and ignores local dataset edits. Existing run directories are protected from overwrite.", "",
                  "A short installation check uses --profile quick. Regenerating reports needs "
                  "the plotting dependencies but does not execute solvers.", "",
                  "## Evidence", "",
                  "- [Raw measurements and detailed validation](measurements.jsonl)",
                  "- [Flat measurement table](measurements.csv)",
                  "- [Aggregated results](summary.csv)",
                  "- [Configuration and provenance](manifest.json)", "",
                  f"Recorded unsuccessful workers (including warm-ups): **{len(problems)}**."])
    for row in problems:
        reason = row.get("error") or ", ".join(row.get("validation", {}).get("errors", []))
        reason = str(reason).replace("\n", " ")[:250]
        lines.append(f"- {row['case_id']}, repeat {row['repeat']}: {row['status']} {reason}")
    (output / "README.md").write_text("\n".join(lines) + "\n")
