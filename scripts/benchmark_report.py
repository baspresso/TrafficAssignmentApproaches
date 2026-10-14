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
    def group(row):
        return row["dataset"], row["approach"], row["tolerance"]

    scaling_groups = {group(s) for s in summary if s["kind"] == "tap" and s["threads"] > 1}
    scaling = [s for s in summary if s["kind"] == "tap" and group(s) in scaling_groups]
    for row in scaling:
        reference = next((s for s in scaling if group(s) == group(row) and s["threads"] == 1), None)
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


def dataset_title(name: str) -> str:
    return {
        "ChicagoSketch": "Chicago Sketch\nTravel-time-only",
        "Anaheim": "Anaheim\nUnrestricted transit; fractional demand",
        "Sydney": "Sydney\nUnrestricted transit; CSV-compatible truncated demand",
        "SiouxFalls": "Sioux Falls",
    }.get(name, name)


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
        columns = min(2, len(datasets))
        rows = math.ceil(len(datasets) / columns)
        fig, axes = plt.subplots(rows, columns, figsize=(6 * columns, 4.3 * rows), squeeze=False)
        for ax in axes.flat[len(datasets):]:
            ax.set_visible(False)
        for ax, dataset in zip(axes.flat, datasets):
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
                   title=dataset_title(dataset))
            ax.grid(alpha=0.2)
            if ax.lines:
                ax.legend()
            else:
                ax.text(0.5, 0.5, "No fully validated configurations", transform=ax.transAxes, ha="center")
        fig.suptitle("Assignment runtime at matched accuracy targets\nMedian and interquartile range; validated configurations", fontsize=12)
        fig.tight_layout()
        save(fig, "assignment_runtime")

    scaling = sorted((s for s in summary if "speedup" in s), key=lambda s: (s["dataset"], s["threads"]))
    if scaling:
        datasets = list(dict.fromkeys(s["dataset"] for s in scaling))
        fig, axes = plt.subplots(len(datasets), 1, figsize=(8, 4 * len(datasets)), squeeze=False)
        for ax, dataset in zip(axes.flat, datasets):
            rows = [s for s in scaling if s["dataset"] == dataset]
            valid = [s for s in rows if s["speedup"] is not None]
            if valid:
                ax.plot([s["threads"] for s in valid], [s["speedup"] for s in valid],
                        "o-", color="#2166ac", label="Measured whole-solve speedup")
                for s in valid:
                    ax.annotate(f"{s['speedup']:.2f}×", (s["threads"], s["speedup"]),
                                xytext=(0, 9), textcoords="offset points", ha="center")
            ax.axhline(1, color="#777777", linestyle="--", linewidth=1, label="Serial baseline")
            ax.set(xticks=[s["threads"] for s in rows], xlabel="Route-search threads",
                   ylabel="Serial median / parallel median", title=dataset_title(dataset))
            ax.set_ylim(bottom=0)
            ax.margins(y=0.2)
            ax.grid(alpha=0.2)
            ax.legend(loc="best")
        fig.suptitle("Parallel route search: whole-solve speedup\nRoute-based assignment; relative-gap target 1e-8", fontsize=12)
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
        "| Dataset | Nodes | Links | Interzonal OD pairs | Modeled demand | Road-network demand |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for name, d in manifest["datasets"].items():
        lines.append(f"| {name} | {d['nodes']} | {d['links']} | {d['interzonal_od_pairs']} | "
                     f"{d['total_demand']:,.2f} | {d['interzonal_demand']:,.2f} |")
    lines.extend([
        "", "Source: [Transportation Networks for Research]"
        f"({manifest['data_repository']}/tree/{manifest['data_revision']}). "
        "Exact revision and input checksums are in [manifest.json](manifest.json). "
        "Fractional demand is preserved except for the explicitly identified Sydney variant. "
        "Intrazonal trips remain in the input as zero-link "
        "trips; they contribute no road flow or travel time. Road-network demand excludes them.", "",
        "The selected networks use BPR travel-time costs and allow transit through every node. "
        "**Chicago Sketch is a travel-time-only variant**: toll and distance terms from "
        "the published generalized-cost reference are excluded, so its published reference "
        "objective does not apply here. Sioux Falls is a small research/debugging "
        "benchmark. Travel-time totals retain native dataset units.", "",
    ])
    if "Anaheim" in manifest["datasets"]:
        lines.extend(["**Anaheim uses unrestricted transit and fractional source demand.** "
                      "This project-specific variant allows intermediate zone nodes; published "
                      "reference objectives for restricted routing do not apply.", ""])
    if "Sydney" in manifest["datasets"]:
        d = manifest["datasets"]["Sydney"]
        lines.extend([
            "**Sydney uses unrestricted transit and CSV-compatible truncated demand.** "
            "Every nonnegative source OD value is rounded down, reproducing the existing "
            "CSV-based `scripts/run_experiment.py tap --dataset Sydney` workload. "
            "The study measures native solve time with statistics disabled, so timings "
            "are not directly comparable to that command's elapsed time. "
            "Source critical-speed and lane columns are not interpreted as tolls or link types.", "",
            f"Sydney's source header declares {number(d.get('declared_source_demand'), ',.0f')} "
            f"trips; its entries sum to {number(d.get('actual_source_demand'), ',.6f')}. "
            f"The modeled total after truncation is {d['total_demand']:,.0f}, including "
            f"{d['intrazonal_demand']:,.0f} intrazonal trips. No demand rescaling is applied. "
            "The known header discrepancy is accepted only for the checksum-verified source. "
            "Results describe this modified workload, not the original full-demand reference case.", "",
        ])
    lines.extend([
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
    scaling = sorted((s for s in summary if "speedup" in s), key=lambda s: (s["dataset"], s["threads"]))
    if scaling:
        lines.extend(["", "## Parallel route-search scaling", "",
                      "Route-based solver, target relative gap 1e-8, with a separate serial "
                      "baseline for each network. "
                      "Speedup measures the entire solve: flow redistribution remains sequential. "
                      "Speedups below 1 indicate a slowdown. Only fully validated configurations "
                      "with matching travel-time and Beckmann objectives receive a speedup.", "",
                      "| Dataset | Threads | Successful runs | Median [Q25, Q75] (s) | Speedup |",
                      "| --- | ---: | ---: | ---: | ---: |"])
        for s in scaling:
            speedup = f"{s['speedup']:.2f}×" if s['speedup'] is not None else "—"
            lines.append(f"| {s['dataset']} | {s['threads']} | {s['successful_runs']}/{s['expected_runs']} | "
                         f"{number(s['solve_seconds'])} [{number(s['q25_seconds'])}, "
                         f"{number(s['q75_seconds'])}] | {speedup} |")
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
                  f"Backend: {manifest.get('validation_backend', 'Python heap-based Dijkstra')}. "
                  "Gap validation allows an absolute 1e-12 roundoff margin.",
                  "- CNDP also checks capacity bounds, investment accounting, and a fresh final "
                  "equilibrium. Budget tolerance is 1e-6 and bound tolerance is 1e-7. "
                  "Fresh and warm-started travel times must agree within 1e-5 relative tolerance.",
                  "- TAP has a 1000-iteration ceiling. Worker limits, including validation: "
                  + ", ".join(f"{name} TAP {max(c['timeout'] for c in manifest['configuration'] if c['dataset'] == name and c['kind'] == 'tap')}s"
                              for name in manifest['datasets'])
                  + f"; CNDP {config['timeout']}s. This run's worker budget is "
                  f"{manifest['time_budget_seconds']:g}s. Errors, timeouts, and unrun cases remain in raw results.",
                  "- Measurements follow the capacity-dependent derivative-cache correction. "
                  "CNDP results describe simulated outcomes under the specified model.", "",
                  "## Machine and build", "",
                  f"- CPU: {manifest['cpu']}; {manifest['available_cpus']} CPUs available to the process.",
                  f"- Platform: {manifest['platform']}.",
                  f"- Compiler: {manifest['build'].get('compiler_version', 'unavailable')}.",
                  f"- Build: {manifest['build'].get('CMAKE_BUILD_TYPE')}; "
                  f"release flags: {manifest['build'].get('CMAKE_CXX_FLAGS_RELEASE', '')}.",
                  f"- Python: {manifest['python']}; NumPy: {manifest['packages'].get('numpy')}; "
                  f"SciPy: {manifest['packages'].get('scipy', 'not used')}.",
                  "- Numerical-library thread environment: "
                  + (", ".join(f"{key}={value}" for key, value in manifest.get('numerical_thread_environment', {}).items()) or "not recorded")
                  + ". Native route-search thread counts are configured separately.",
                  f"- Git base: {manifest['git_commit']}; uncommitted changes present: "
                  f"{str(manifest['working_tree_dirty']).lower()}.",
                  f"- Source fingerprint: {manifest['source_sha256']}.",
                  "- The manifest records individual source hashes, the native extension hash, "
                  "dependency versions, and exact configurations. Runtime depends on hardware, "
                  "system load, compiler, and numerical stopping behavior.", "",
                  "## Reproduce", "",
                  "From the repository root, after installing the native prerequisites in the main README:", "",
                  chr(96) * 3 + "bash",
                  "python -m pip install -e '.[benchmarks,test]' -Ccmake.build-type=Release",
                  f"python scripts/run_benchmarks.py --profile {manifest['profile']} \\",
                  f"  --time-budget {manifest['time_budget_seconds']:g} --output performance_results/benchmark-reproduction",
                  "python scripts/run_benchmarks.py \\",
                  "  --render-only performance_results/benchmark-reproduction",
                  chr(96) * 3, "",
                  "The runner downloads only the pinned TNTP inputs into build/benchmark-data. "
                  "For offline preparation, add --source-root data/TransportationNetworks "
                  "with a Git clone containing the pinned revision; it reads committed blobs "
                  "and ignores local dataset edits. Sydney's ZIP demand is streamed and converted once "
                  "into a checked array cache. Existing run directories are protected from overwrite.", "",
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
