"""
CNDP Results Analysis - Derived Metrics and LaTeX Tables

Computes publication-quality metrics from trace CSVs and outputs:
- LaTeX table files for direct \\input{} in papers
- JSON results digest
- Console summary table

Metrics computed:
- Optimality gap: (Z* - Z_etalon) / Z_etalon x 100%
- Time-to-target: time to first reach within X% of etalon (1%, 0.5%, 0.1%)
- TA evaluations to target: trace row count to reach within X% of etalon
- Per-iteration objective reduction: mean delta-Z per OptCond iteration
- Budget feasibility rate: fraction of trace points with BudgetViolation <= 0.1

Usage:
    python scripts/analyze_cndp_results.py --dataset SiouxFalls
    python scripts/analyze_cndp_results.py --dataset Anaheim
    python scripts/analyze_cndp_results.py --datasets SiouxFalls,Anaheim
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent

GAP_TARGETS = [1.0, 0.5, 0.1]  # percent


def load_trace_files(dataset: str) -> dict[str, pd.DataFrame]:
    """Load all trace CSVs for a dataset, grouped by scenario name (latest run)."""
    trace_dir = PROJECT_ROOT / "performance_results" / dataset
    if not trace_dir.exists():
        print(f"ERROR: Trace directory not found: {trace_dir}")
        sys.exit(1)

    traces = {}
    for csv_path in sorted(trace_dir.glob("BilevelCND_*_quality_time.csv")):
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"WARNING: Failed to read {csv_path.name}: {e}")
            continue

        if df.empty:
            continue

        if "Scenario" not in df.columns or df["Scenario"].isna().all():
            run_id = df["RunId"].iloc[0] if "RunId" in df.columns else None
            if run_id:
                meta_path = csv_path.parent / csv_path.name.replace(
                    "_quality_time.csv", "_metadata.json"
                )
                if meta_path.exists():
                    with open(meta_path) as f:
                        meta = json.load(f)
                    scenario_name = meta.get("scenario", "")
                    if scenario_name:
                        df["Scenario"] = scenario_name

        if "Scenario" not in df.columns or df.empty or df["Scenario"].isna().all():
            continue

        scenario_name = df["Scenario"].iloc[0]
        if not scenario_name or pd.isna(scenario_name):
            continue

        if scenario_name in traces:
            existing_run = traces[scenario_name]["RunId"].iloc[0]
            new_run = df["RunId"].iloc[0]
            if str(new_run) > str(existing_run):
                traces[scenario_name] = df
        else:
            traces[scenario_name] = df

    return traces


def filter_budget_feasible(df: pd.DataFrame, tolerance: float = 0.1) -> pd.DataFrame:
    if "BudgetViolation" in df.columns:
        return df[df["BudgetViolation"] <= tolerance].copy()
    return df.copy()


def compute_best_so_far(series: pd.Series) -> pd.Series:
    best = series.copy()
    current_best = float("inf")
    for i in range(len(best)):
        val = best.iloc[i]
        if np.isfinite(val) and val < current_best:
            current_best = val
        best.iloc[i] = current_best if np.isfinite(current_best) else np.nan
    return best


def get_etalon_best_objective(traces: dict[str, pd.DataFrame]) -> float | None:
    if "Etalon" not in traces:
        return None
    df = traces["Etalon"]
    feasible = filter_budget_feasible(df)
    if feasible.empty:
        return None
    if "BestFeasibleObjective" in feasible.columns:
        vals = feasible["BestFeasibleObjective"].dropna()
        if not vals.empty:
            return vals.iloc[-1]
    return compute_best_so_far(feasible["Objective"]).iloc[-1]


def compute_scenario_metrics(
    scenario_name: str, df: pd.DataFrame, etalon_obj: float | None
) -> dict:
    """Compute all derived metrics for a single scenario."""
    feasible = filter_budget_feasible(df)
    total_time = df["ElapsedTime(s)"].max() if not df.empty else 0
    n_evals = len(df)

    # Best feasible objective
    if not feasible.empty:
        best_envelope = compute_best_so_far(feasible["Objective"])
        best_obj = best_envelope.iloc[-1]
    else:
        best_obj = np.nan
        best_envelope = pd.Series(dtype=float)

    # Optimality gap
    if etalon_obj is not None and np.isfinite(etalon_obj) and etalon_obj > 0 and np.isfinite(best_obj):
        gap_pct = (best_obj - etalon_obj) / etalon_obj * 100
    else:
        gap_pct = np.nan

    # Time-to-target and evaluations-to-target
    time_to_target = {}
    evals_to_target = {}
    if etalon_obj is not None and np.isfinite(etalon_obj) and not feasible.empty:
        for target_pct in GAP_TARGETS:
            threshold = etalon_obj * (1 + target_pct / 100)
            reached = best_envelope <= threshold
            if reached.any():
                first_idx = reached.idxmax()
                pos = feasible.index.get_loc(first_idx)
                time_to_target[target_pct] = feasible["ElapsedTime(s)"].iloc[pos]
                evals_to_target[target_pct] = pos + 1
            else:
                time_to_target[target_pct] = np.nan
                evals_to_target[target_pct] = np.nan

    # Per-iteration objective reduction for OptCond phases
    optcond_delta = np.nan
    if "Phase" in df.columns:
        oc_rows = df[df["Phase"].str.startswith("optcond", na=False)]
        if len(oc_rows) >= 2:
            oc_feasible = filter_budget_feasible(oc_rows)
            if len(oc_feasible) >= 2:
                deltas = oc_feasible["Objective"].diff().dropna()
                negative_deltas = deltas[deltas < 0]
                if not negative_deltas.empty:
                    optcond_delta = negative_deltas.mean()

    # Budget feasibility rate
    if "BudgetViolation" in df.columns and not df.empty:
        feasibility_rate = (df["BudgetViolation"] <= 0.1).mean()
    else:
        feasibility_rate = np.nan

    # Average TA compute time
    avg_ta_time = np.nan
    if "TAComputeTime(s)" in df.columns:
        avg_ta_time = df["TAComputeTime(s)"].mean()

    return {
        "Scenario": scenario_name,
        "BestObjective": best_obj,
        "OptimalityGap(%)": gap_pct,
        "ElapsedTime(s)": total_time,
        "Evaluations": n_evals,
        "AvgTATime(s)": avg_ta_time,
        "OptCondDelta": optcond_delta,
        "FeasibilityRate": feasibility_rate,
        "TimeToTarget": time_to_target,
        "EvalsToTarget": evals_to_target,
    }


def analyze_dataset(dataset: str) -> tuple[list[dict], float | None]:
    """Analyze all scenarios for a dataset, return (metrics_list, etalon_obj)."""
    print(f"\nLoading traces for {dataset}...")
    traces = load_trace_files(dataset)

    if not traces:
        print(f"ERROR: No trace files found for {dataset}")
        return [], None

    etalon_obj = get_etalon_best_objective(traces)
    if etalon_obj is not None:
        print(f"Etalon reference objective: {etalon_obj:.4f}")
    else:
        print("WARNING: No Etalon scenario found - gap metrics unavailable")

    all_metrics = []
    for scenario_name in sorted(traces.keys()):
        if scenario_name == "Etalon":
            continue
        metrics = compute_scenario_metrics(scenario_name, traces[scenario_name], etalon_obj)
        all_metrics.append(metrics)

    return all_metrics, etalon_obj


def print_console_summary(metrics_list: list[dict], dataset: str):
    """Print a formatted console summary table."""
    print(f"\n{'='*110}")
    print(f"ANALYSIS RESULTS: {dataset}")
    print(f"{'='*110}")
    header = (
        f"{'Scenario':<25} {'BestObj':>14} {'Gap(%)':>8} {'Time(s)':>8} "
        f"{'Evals':>6} {'FeasRate':>8} {'AvgTA(s)':>9}"
    )
    print(header)
    print("-" * 110)
    for m in metrics_list:
        obj_str = f"{m['BestObjective']:.2f}" if np.isfinite(m["BestObjective"]) else "N/A"
        gap_str = f"{m['OptimalityGap(%)']:.3f}" if np.isfinite(m["OptimalityGap(%)"]) else "N/A"
        feas_str = f"{m['FeasibilityRate']:.2%}" if np.isfinite(m["FeasibilityRate"]) else "N/A"
        ta_str = f"{m['AvgTATime(s)']:.3f}" if np.isfinite(m["AvgTATime(s)"]) else "N/A"
        print(
            f"{m['Scenario']:<25} {obj_str:>14} {gap_str:>8} "
            f"{m['ElapsedTime(s)']:>8.1f} {m['Evaluations']:>6} "
            f"{feas_str:>8} {ta_str:>9}"
        )

    # Time-to-target sub-table
    has_targets = any(m["TimeToTarget"] for m in metrics_list)
    if has_targets:
        print(f"\nTime-to-target (seconds to reach within X% of etalon):")
        print(f"{'Scenario':<25}", end="")
        for pct in GAP_TARGETS:
            print(f" {'<' + str(pct) + '%':>10}", end="")
        print()
        print("-" * (25 + 11 * len(GAP_TARGETS)))
        for m in metrics_list:
            print(f"{m['Scenario']:<25}", end="")
            for pct in GAP_TARGETS:
                val = m["TimeToTarget"].get(pct, np.nan)
                if np.isfinite(val):
                    print(f" {val:>10.1f}", end="")
                else:
                    print(f" {'---':>10}", end="")
            print()

        print(f"\nTA evaluations to reach within X% of etalon:")
        print(f"{'Scenario':<25}", end="")
        for pct in GAP_TARGETS:
            print(f" {'<' + str(pct) + '%':>10}", end="")
        print()
        print("-" * (25 + 11 * len(GAP_TARGETS)))
        for m in metrics_list:
            print(f"{m['Scenario']:<25}", end="")
            for pct in GAP_TARGETS:
                val = m["EvalsToTarget"].get(pct, np.nan)
                if np.isfinite(val):
                    print(f" {int(val):>10}", end="")
                else:
                    print(f" {'---':>10}", end="")
            print()


def write_latex_table(metrics_list: list[dict], dataset: str, etalon_obj: float | None, output_dir: Path):
    """Write a LaTeX table file for the dataset."""
    output_dir.mkdir(parents=True, exist_ok=True)
    tex_path = output_dir / f"results_table_{dataset}.tex"

    lines = []
    lines.append("% Auto-generated by analyze_cndp_results.py")
    if etalon_obj is not None:
        lines.append(f"% Etalon reference: {etalon_obj:.4f}")
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(f"\\caption{{CNDP Results -- {dataset}}}")
    lines.append(f"\\label{{tab:cndp-{dataset.lower()}}}")
    lines.append(r"\begin{tabular}{lrrrrrr}")
    lines.append(r"\toprule")
    lines.append(
        r"Scenario & Best Obj. & Gap (\%) & Time (s) & Evals & Feas. Rate & Avg TA (s) \\"
    )
    lines.append(r"\midrule")

    for m in metrics_list:
        name_tex = m["Scenario"].replace("_", r"\_").replace("-", r"-")
        obj_str = f"{m['BestObjective']:.2f}" if np.isfinite(m["BestObjective"]) else "---"
        gap_str = f"{m['OptimalityGap(%)']:.3f}" if np.isfinite(m["OptimalityGap(%)"]) else "---"
        time_str = f"{m['ElapsedTime(s)']:.1f}"
        evals_str = str(m["Evaluations"])
        feas_str = f"{m['FeasibilityRate']:.1%}" if np.isfinite(m["FeasibilityRate"]) else "---"
        ta_str = f"{m['AvgTATime(s)']:.3f}" if np.isfinite(m["AvgTATime(s)"]) else "---"
        lines.append(
            f"  {name_tex} & {obj_str} & {gap_str} & {time_str} & {evals_str} & {feas_str} & {ta_str} \\\\"
        )

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    # Time-to-target table
    has_targets = any(m["TimeToTarget"] for m in metrics_list)
    if has_targets:
        lines.append("")
        lines.append(r"\begin{table}[htbp]")
        lines.append(r"\centering")
        lines.append(f"\\caption{{Time-to-Target -- {dataset}}}")
        lines.append(f"\\label{{tab:ttt-{dataset.lower()}}}")
        cols = "l" + "r" * len(GAP_TARGETS) * 2
        lines.append(f"\\begin{{tabular}}{{{cols}}}")
        lines.append(r"\toprule")
        header_parts = ["Scenario"]
        for pct in GAP_TARGETS:
            header_parts.append(f"$t_{{<{pct}\\%}}$ (s)")
            header_parts.append(f"$n_{{<{pct}\\%}}$")
        lines.append(" & ".join(header_parts) + r" \\")
        lines.append(r"\midrule")

        for m in metrics_list:
            name_tex = m["Scenario"].replace("_", r"\_").replace("-", r"-")
            parts = [name_tex]
            for pct in GAP_TARGETS:
                t_val = m["TimeToTarget"].get(pct, np.nan)
                e_val = m["EvalsToTarget"].get(pct, np.nan)
                parts.append(f"{t_val:.1f}" if np.isfinite(t_val) else "---")
                parts.append(f"{int(e_val)}" if np.isfinite(e_val) else "---")
            lines.append("  " + " & ".join(parts) + r" \\")

        lines.append(r"\bottomrule")
        lines.append(r"\end{tabular}")
        lines.append(r"\end{table}")

    tex_path.write_text("\n".join(lines))
    print(f"LaTeX table written to: {tex_path}")


def write_json_digest(metrics_list: list[dict], dataset: str, etalon_obj: float | None, output_dir: Path):
    """Write a JSON results digest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"results_digest_{dataset}.json"

    digest = {
        "dataset": dataset,
        "etalon_objective": etalon_obj,
        "scenarios": [],
    }

    for m in metrics_list:
        entry = {
            "scenario": m["Scenario"],
            "best_objective": m["BestObjective"] if np.isfinite(m["BestObjective"]) else None,
            "optimality_gap_pct": m["OptimalityGap(%)"] if np.isfinite(m["OptimalityGap(%)"]) else None,
            "elapsed_time_s": m["ElapsedTime(s)"],
            "evaluations": m["Evaluations"],
            "avg_ta_time_s": m["AvgTATime(s)"] if np.isfinite(m["AvgTATime(s)"]) else None,
            "optcond_delta": m["OptCondDelta"] if np.isfinite(m["OptCondDelta"]) else None,
            "feasibility_rate": m["FeasibilityRate"] if np.isfinite(m["FeasibilityRate"]) else None,
            "time_to_target": {
                str(k): (v if np.isfinite(v) else None)
                for k, v in m["TimeToTarget"].items()
            },
            "evals_to_target": {
                str(k): (int(v) if np.isfinite(v) else None)
                for k, v in m["EvalsToTarget"].items()
            },
        }
        digest["scenarios"].append(entry)

    json_path.write_text(json.dumps(digest, indent=2))
    print(f"JSON digest written to: {json_path}")


def write_cross_dataset_table(
    all_metrics: dict[str, list[dict]],
    all_etalons: dict[str, float | None],
    datasets: list[str],
    output_dir: Path,
):
    """Write a cross-dataset LaTeX comparison table."""
    output_dir.mkdir(parents=True, exist_ok=True)
    tex_path = output_dir / "results_table_cross_dataset.tex"

    # Collect common scenarios
    scenario_sets = [set(m["Scenario"] for m in mlist) for mlist in all_metrics.values()]
    common = sorted(set.intersection(*scenario_sets)) if scenario_sets else []

    lines = []
    lines.append("% Auto-generated by analyze_cndp_results.py (cross-dataset)")
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Cross-Dataset CNDP Comparison}")
    lines.append(r"\label{tab:cndp-cross}")

    # Columns: Scenario + (Gap, Time) per dataset
    col_spec = "l" + "rr" * len(datasets)
    lines.append(f"\\begin{{tabular}}{{{col_spec}}}")
    lines.append(r"\toprule")

    header_parts = ["Scenario"]
    for ds in datasets:
        header_parts.append(f"Gap (\\%) [{ds}]")
        header_parts.append(f"Time (s) [{ds}]")
    lines.append(" & ".join(header_parts) + r" \\")
    lines.append(r"\midrule")

    # Build lookup
    ds_lookup = {}
    for ds, mlist in all_metrics.items():
        ds_lookup[ds] = {m["Scenario"]: m for m in mlist}

    for scenario in common:
        name_tex = scenario.replace("_", r"\_").replace("-", r"-")
        parts = [name_tex]
        for ds in datasets:
            m = ds_lookup.get(ds, {}).get(scenario)
            if m:
                gap = m["OptimalityGap(%)"]
                parts.append(f"{gap:.3f}" if np.isfinite(gap) else "---")
                parts.append(f"{m['ElapsedTime(s)']:.1f}")
            else:
                parts.extend(["---", "---"])
        lines.append("  " + " & ".join(parts) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    tex_path.write_text("\n".join(lines))
    print(f"Cross-dataset LaTeX table written to: {tex_path}")


def main():
    parser = argparse.ArgumentParser(description="Analyze CNDP comparison results")
    parser.add_argument(
        "--dataset",
        choices=["SiouxFalls", "Anaheim"],
        default=None,
        help="Single dataset to analyze",
    )
    parser.add_argument(
        "--datasets",
        type=str,
        default=None,
        help="Comma-separated datasets for cross-dataset analysis (e.g. SiouxFalls,Anaheim)",
    )

    args = parser.parse_args()

    if not args.dataset and not args.datasets:
        args.dataset = "SiouxFalls"

    if args.datasets:
        dataset_list = [d.strip() for d in args.datasets.split(",")]
    elif args.dataset:
        dataset_list = [args.dataset]
    else:
        dataset_list = []

    all_metrics = {}
    all_etalons = {}

    for ds in dataset_list:
        metrics_list, etalon_obj = analyze_dataset(ds)
        if not metrics_list:
            continue

        all_metrics[ds] = metrics_list
        all_etalons[ds] = etalon_obj

        print_console_summary(metrics_list, ds)

        output_dir = PROJECT_ROOT / "performance_results" / ds
        write_latex_table(metrics_list, ds, etalon_obj, output_dir)
        write_json_digest(metrics_list, ds, etalon_obj, output_dir)

    # Cross-dataset table if multiple datasets
    if len(all_metrics) >= 2:
        output_dir = PROJECT_ROOT / "performance_results"
        write_cross_dataset_table(all_metrics, all_etalons, dataset_list, output_dir)


if __name__ == "__main__":
    main()
