"""
CNDP Algorithm Comparison - Analysis and Publication-Ready Plots

Reads trace CSVs produced by run_cndp_comparison.py and generates:
1. Objective vs Time - all scenarios overlaid (etalon as reference line)
2. Best-so-far objective vs Time - monotonic envelope (etalon as reference line)
3. Optimality gap vs TA evaluations (log scale)
4. Budget utilization vs time
5. Per-iteration convergence character (side-by-side OptCond vs COBYLA)
6. Sensitivity sweep plots (when sweep data exists)
7. Summary table (printed + saved as CSV)

Multi-run support: detects _runN suffixed scenarios and plots mean +/- 1 sigma.
Cross-dataset mode: --datasets SiouxFalls,Anaheim for scalability plots.

Usage:
    python scripts/plot_cndp_comparison.py --dataset SiouxFalls
    python scripts/plot_cndp_comparison.py --dataset Anaheim --format pdf
    python scripts/plot_cndp_comparison.py --datasets SiouxFalls,Anaheim
"""

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

PROJECT_ROOT = Path(__file__).resolve().parent.parent

SCENARIO_STYLES = {
    "OptCond-only": {"color": "#d62728", "linestyle": "-", "linewidth": 2.0},
    "COBYLA-only": {"color": "#1f77b4", "linestyle": "-", "linewidth": 1.5},
    "BOBYQA-only": {"color": "#17becf", "linestyle": "-", "linewidth": 1.5},
    "ISRES-only": {"color": "#ff7f0e", "linestyle": "-", "linewidth": 1.5},
    "DE-only": {"color": "#e377c2", "linestyle": "-", "linewidth": 1.5},
    "PSO-only": {"color": "#7f7f7f", "linestyle": "-", "linewidth": 1.5},
    "COBYLA-then-OptCond": {"color": "#2ca02c", "linestyle": "-", "linewidth": 1.5},
    "OptCond-then-COBYLA": {"color": "#8c564b", "linestyle": "-.", "linewidth": 1.5},
    "DE-then-OptCond": {"color": "#9467bd", "linestyle": "--", "linewidth": 1.5},
}

# Marker placed at each step switch point
SWITCH_MARKER = {"marker": "D", "markersize": 7, "zorder": 5}

# Network sizes for cross-dataset plots
NETWORK_LINKS = {"SiouxFalls": 76, "Anaheim": 914}


def load_trace_files(dataset: str) -> dict[str, pd.DataFrame]:
    """Load all trace CSVs for a dataset, grouped by scenario name.

    For each scenario name, keeps the latest run (by RunId).
    """
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
            print(f"  Skipping {csv_path.name} (no scenario label)")
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


def load_trace_files_from_run_dir(run_dir: Path) -> dict[str, pd.DataFrame]:
    """Load trace CSVs from a self-contained run folder's traces/ subfolder."""
    traces_dir = run_dir / "traces"
    if not traces_dir.exists():
        print(f"ERROR: Traces directory not found: {traces_dir}")
        sys.exit(1)

    traces = {}
    for csv_path in sorted(traces_dir.glob("*_trace.csv")):
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"WARNING: Failed to read {csv_path.name}: {e}")
            continue
        if df.empty:
            continue

        # Derive scenario name from filename: {ScenarioName}_trace.csv
        scenario_name = csv_path.stem.replace("_trace", "")

        # Ensure Scenario column is populated
        if "Scenario" not in df.columns or df["Scenario"].isna().all():
            df["Scenario"] = scenario_name

        traces[scenario_name] = df

    return traces


def get_etalon_best_objective(traces: dict[str, pd.DataFrame]) -> float | None:
    """Extract the best feasible objective from the Etalon scenario."""
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


def detect_multi_run_groups(traces: dict[str, pd.DataFrame]) -> dict[str, list[str]]:
    """Detect multi-run scenario groups matching pattern 'ScenarioName_runN'.

    Returns {base_name: [scenario_name_run1, scenario_name_run2, ...]}.
    """
    pattern = re.compile(r"^(.+)_run(\d+)$")
    groups: dict[str, list[str]] = {}
    for name in traces:
        m = pattern.match(name)
        if m:
            base = m.group(1)
            groups.setdefault(base, []).append(name)
    # Sort each group by run number
    for base in groups:
        groups[base].sort(key=lambda n: int(pattern.match(n).group(2)))
    return groups


def interpolate_to_common_grid(
    dfs: list[pd.DataFrame], x_col: str, y_col: str, n_points: int = 500
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate multiple traces to a common time grid, return (x, mean_y, std_y)."""
    max_x = min(df[x_col].max() for df in dfs)
    x_grid = np.linspace(0, max_x, n_points)

    y_interp = []
    for df in dfs:
        y_interp.append(np.interp(x_grid, df[x_col].values, df[y_col].values))
    y_matrix = np.array(y_interp)
    return x_grid, np.mean(y_matrix, axis=0), np.std(y_matrix, axis=0)


def filter_budget_feasible(df: pd.DataFrame, tolerance: float = 0.1) -> pd.DataFrame:
    """Filter to only budget-feasible points."""
    if "BudgetViolation" in df.columns:
        return df[df["BudgetViolation"] <= tolerance].copy()
    return df.copy()


def compute_best_so_far(series: pd.Series) -> pd.Series:
    """Compute the monotonically decreasing best-so-far envelope."""
    best = series.copy()
    current_best = float("inf")
    for i in range(len(best)):
        val = best.iloc[i]
        if np.isfinite(val) and val < current_best:
            current_best = val
        best.iloc[i] = current_best if np.isfinite(current_best) else np.nan
    return best


def find_step_switch_indices(df: pd.DataFrame) -> list[int]:
    """Find row indices where the optimization step switches.

    Detected by Phase column: a post_* row marks the end of a step.
    The switch point is the post_* row itself (last point of the completed step).
    """
    if "Phase" not in df.columns:
        return []

    indices = []
    phases = df["Phase"].values
    for i, phase in enumerate(phases):
        if isinstance(phase, str) and phase.startswith("post_"):
            indices.append(i)
    return indices


def get_style(scenario_name: str) -> dict:
    """Get plot style for a scenario."""
    return SCENARIO_STYLES.get(scenario_name, {
        "color": "gray", "linestyle": "-", "linewidth": 1.0,
    })


def exclude_etalon(traces: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    """Return traces without the Etalon scenario."""
    return {k: v for k, v in traces.items() if k != "Etalon"}


def _add_etalon_reference(ax: plt.Axes, etalon_obj: float | None):
    """Draw a horizontal reference line at the etalon best objective."""
    if etalon_obj is not None and np.isfinite(etalon_obj):
        ax.axhline(
            etalon_obj,
            color="gray",
            linestyle="--",
            linewidth=1.0,
            alpha=0.7,
            label=f"Etalon ({etalon_obj:.0f})",
            zorder=1,
        )


def plot_objective_vs_time(traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path):
    """Plot 1: Objective vs elapsed time for all scenarios (etalon as reference line)."""
    etalon_obj = get_etalon_best_objective(traces)
    plot_traces = exclude_etalon(traces)
    fig, ax = plt.subplots(figsize=(10, 6))

    _add_etalon_reference(ax, etalon_obj)

    for scenario_name, df in sorted(plot_traces.items()):
        feasible = filter_budget_feasible(df)
        if feasible.empty:
            continue
        style = get_style(scenario_name)
        ax.plot(
            feasible["ElapsedTime(s)"],
            feasible["Objective"],
            label=scenario_name,
            **style,
        )

        # Mark step switches
        switch_indices = find_step_switch_indices(df)
        feasible_idx_set = set(feasible.index)
        for idx in switch_indices:
            if idx in feasible_idx_set:
                row = df.loc[idx]
                ax.plot(
                    row["ElapsedTime(s)"],
                    row["Objective"],
                    color=style.get("color", "gray"),
                    **SWITCH_MARKER,
                )

    ax.set_xlabel("Elapsed Time (s)")
    ax.set_ylabel("Objective Function")
    ax.set_title(f"CNDP Objective vs Time - {dataset}")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.savefig(fig_dir / f"objective_vs_time_{dataset}.png")
    fig.savefig(fig_dir / f"objective_vs_time_{dataset}.pdf")
    plt.close(fig)
    print(f"  Saved: objective_vs_time_{dataset}")


def plot_best_so_far(traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path):
    """Plot 2: Best-so-far feasible objective vs time (etalon as reference line)."""
    etalon_obj = get_etalon_best_objective(traces)
    plot_traces = exclude_etalon(traces)
    fig, ax = plt.subplots(figsize=(10, 6))

    _add_etalon_reference(ax, etalon_obj)

    for scenario_name, df in sorted(plot_traces.items()):
        feasible = filter_budget_feasible(df)
        if feasible.empty:
            continue

        best_envelope = compute_best_so_far(feasible["Objective"])
        style = get_style(scenario_name)
        ax.step(
            feasible["ElapsedTime(s)"],
            best_envelope,
            label=scenario_name,
            where="post",
            **style,
        )

        # Mark step switches on the envelope
        switch_indices = find_step_switch_indices(df)
        feasible_idx_set = set(feasible.index)
        for idx in switch_indices:
            if idx in feasible_idx_set:
                pos = feasible.index.get_loc(idx)
                ax.plot(
                    feasible["ElapsedTime(s)"].iloc[pos],
                    best_envelope.iloc[pos],
                    color=style.get("color", "gray"),
                    **SWITCH_MARKER,
                )

    ax.set_xlabel("Elapsed Time (s)")
    ax.set_ylabel("Best Feasible Objective")
    ax.set_title(f"Best-So-Far Objective vs Time - {dataset}")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.savefig(fig_dir / f"best_so_far_{dataset}.png")
    fig.savefig(fig_dir / f"best_so_far_{dataset}.pdf")
    plt.close(fig)
    print(f"  Saved: best_so_far_{dataset}")


def plot_optimality_gap_vs_evaluations(
    traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path
):
    """Plot 3: Optimality gap vs TA evaluations (log scale).

    Y-axis: log10((Z - Z_etalon) / Z_etalon). Decouples algorithm efficiency
    from per-evaluation cost.
    """
    etalon_obj = get_etalon_best_objective(traces)
    if etalon_obj is None or not np.isfinite(etalon_obj) or etalon_obj <= 0:
        print("  Skipping optimality gap plot (no etalon reference)")
        return

    plot_traces = exclude_etalon(traces)
    fig, ax = plt.subplots(figsize=(10, 6))

    for scenario_name, df in sorted(plot_traces.items()):
        feasible = filter_budget_feasible(df)
        if feasible.empty:
            continue

        best_envelope = compute_best_so_far(feasible["Objective"])
        gap = (best_envelope - etalon_obj) / etalon_obj
        # Only plot positive gaps (log scale)
        valid = gap > 0
        if not valid.any():
            continue

        eval_idx = np.arange(len(feasible))
        style = get_style(scenario_name)
        ax.plot(
            eval_idx[valid],
            np.log10(gap[valid]),
            label=scenario_name,
            **style,
        )

        # Mark step switches
        switch_indices = find_step_switch_indices(df)
        feasible_idx_set = set(feasible.index)
        for idx in switch_indices:
            if idx in feasible_idx_set:
                pos = feasible.index.get_loc(idx)
                if valid.iloc[pos]:
                    ax.plot(
                        pos,
                        np.log10(gap.iloc[pos]),
                        color=style.get("color", "gray"),
                        **SWITCH_MARKER,
                    )

    ax.set_xlabel("TA Evaluations (cumulative)")
    ax.set_ylabel(r"$\log_{10}$ Optimality Gap $(Z^* - Z_{etalon}) / Z_{etalon}$")
    ax.set_title(f"Optimality Gap vs TA Evaluations - {dataset}")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.savefig(fig_dir / f"gap_vs_evaluations_{dataset}.png")
    fig.savefig(fig_dir / f"gap_vs_evaluations_{dataset}.pdf")
    plt.close(fig)
    print(f"  Saved: gap_vs_evaluations_{dataset}")


def plot_budget_utilization(
    traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path,
    metadata_dir: Path | None = None,
):
    """Plot 4: Budget utilization vs time.

    Shows how each method explores the feasibility region.
    """
    plot_traces = exclude_etalon(traces)
    fig, ax = plt.subplots(figsize=(10, 6))

    # Try to find budget upper bound from metadata
    budget_upper = None
    search_dirs = []
    if metadata_dir is not None:
        search_dirs.append(metadata_dir)
    search_dirs.append(PROJECT_ROOT / "performance_results" / dataset)
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
        for meta_path in search_dir.glob("*_metadata.json"):
            with open(meta_path) as f:
                meta = json.load(f)
            if "budget_upper_bound" in meta:
                budget_upper = meta["budget_upper_bound"]
                break
        if budget_upper is not None:
            break

    if budget_upper is not None:
        ax.axhline(
            budget_upper,
            color="red",
            linestyle="--",
            linewidth=1.0,
            alpha=0.7,
            label=f"Budget bound ({budget_upper:,.0f})",
            zorder=1,
        )

    for scenario_name, df in sorted(plot_traces.items()):
        if "Budget" not in df.columns:
            continue
        style = get_style(scenario_name)
        ax.plot(
            df["ElapsedTime(s)"],
            df["Budget"],
            label=scenario_name,
            **style,
        )

    ax.set_xlabel("Elapsed Time (s)")
    ax.set_ylabel("Budget Utilization")
    ax.set_title(f"Budget Utilization vs Time - {dataset}")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

    fig.savefig(fig_dir / f"budget_vs_time_{dataset}.png")
    fig.savefig(fig_dir / f"budget_vs_time_{dataset}.pdf")
    plt.close(fig)
    print(f"  Saved: budget_vs_time_{dataset}")


def plot_convergence_character(traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path):
    """Plot 5: Per-iteration convergence character (side-by-side).

    Left: OptCond-only (smooth, monotonic). Right: COBYLA-only (noisy).
    Raw objective per iteration, not best-so-far envelope.
    """
    optcond_df = traces.get("OptCond-only")
    cobyla_df = traces.get("COBYLA-only")

    if optcond_df is None and cobyla_df is None:
        print("  Skipping convergence character plot (need OptCond-only or COBYLA-only)")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left panel: OptCond-only
    if optcond_df is not None:
        feasible = filter_budget_feasible(optcond_df)
        if not feasible.empty:
            axes[0].plot(
                range(len(feasible)),
                feasible["Objective"].values,
                color="#d62728",
                linewidth=1.5,
            )
    axes[0].set_xlabel("Iteration")
    axes[0].set_ylabel("Objective Function")
    axes[0].set_title("OptCond-only (deterministic)")
    axes[0].grid(True, alpha=0.3)

    # Right panel: COBYLA-only
    if cobyla_df is not None:
        feasible = filter_budget_feasible(cobyla_df)
        if not feasible.empty:
            axes[1].plot(
                range(len(feasible)),
                feasible["Objective"].values,
                color="#1f77b4",
                linewidth=1.0,
            )
    axes[1].set_xlabel("Iteration")
    axes[1].set_ylabel("Objective Function")
    axes[1].set_title("COBYLA-only (derivative-free)")
    axes[1].grid(True, alpha=0.3)

    fig.suptitle(f"Convergence Character Comparison - {dataset}", fontsize=13)
    fig.tight_layout()

    fig.savefig(fig_dir / f"convergence_character_{dataset}.png")
    fig.savefig(fig_dir / f"convergence_character_{dataset}.pdf")
    plt.close(fig)
    print(f"  Saved: convergence_character_{dataset}")


def plot_sensitivity_sweeps(traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path):
    """Plot 6: Sensitivity sweep results (when sweep data exists).

    Three subplots: OptCond sweep, COBYLA sweep, Hybrid ratio sweep.
    """
    optcond_sweep = {}
    cobyla_sweep = {}
    hybrid_sweep = {}

    for name, df in traces.items():
        m = re.match(r"OptCond-sweep-(\d+)", name)
        if m:
            optcond_sweep[int(m.group(1))] = df
            continue
        m = re.match(r"COBYLA-sweep-(\d+)", name)
        if m:
            cobyla_sweep[int(m.group(1))] = df
            continue
        m = re.match(r"Hybrid-(\d+)c-(\d+)oc", name)
        if m:
            nlopt_n, oc_n = int(m.group(1)), int(m.group(2))
            total = nlopt_n + oc_n
            frac = nlopt_n / total if total > 0 else 0
            hybrid_sweep[frac] = (name, df)

    if not optcond_sweep and not cobyla_sweep and not hybrid_sweep:
        return  # No sweep data

    n_panels = sum(1 for s in [optcond_sweep, cobyla_sweep, hybrid_sweep] if s)
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    panel_idx = 0

    if optcond_sweep:
        ax = axes[panel_idx]
        iters = sorted(optcond_sweep.keys())
        best_objs = []
        for n in iters:
            df = optcond_sweep[n]
            feasible = filter_budget_feasible(df)
            if not feasible.empty:
                best = compute_best_so_far(feasible["Objective"]).iloc[-1]
                best_objs.append(best)
            else:
                best_objs.append(np.nan)
        ax.plot(iters, best_objs, "o-", color="#d62728", linewidth=1.5, markersize=6)
        ax.set_xlabel("OptCond Iterations")
        ax.set_ylabel("Best Feasible Objective")
        ax.set_title("OptCond Iteration Sweep")
        ax.grid(True, alpha=0.3)
        panel_idx += 1

    if cobyla_sweep:
        ax = axes[panel_idx]
        iters = sorted(cobyla_sweep.keys())
        best_objs = []
        for n in iters:
            df = cobyla_sweep[n]
            feasible = filter_budget_feasible(df)
            if not feasible.empty:
                best = compute_best_so_far(feasible["Objective"]).iloc[-1]
                best_objs.append(best)
            else:
                best_objs.append(np.nan)
        ax.plot(iters, best_objs, "s-", color="#1f77b4", linewidth=1.5, markersize=6)
        ax.set_xlabel("COBYLA Iterations")
        ax.set_ylabel("Best Feasible Objective")
        ax.set_title("COBYLA Iteration Sweep")
        ax.grid(True, alpha=0.3)
        panel_idx += 1

    if hybrid_sweep:
        ax = axes[panel_idx]
        fracs = sorted(hybrid_sweep.keys())
        best_objs = []
        labels = []
        for frac in fracs:
            name, df = hybrid_sweep[frac]
            feasible = filter_budget_feasible(df)
            if not feasible.empty:
                best = compute_best_so_far(feasible["Objective"]).iloc[-1]
                best_objs.append(best)
            else:
                best_objs.append(np.nan)
            labels.append(name.replace("Hybrid-", ""))
        ax.plot(fracs, best_objs, "D-", color="#2ca02c", linewidth=1.5, markersize=6)
        ax.set_xlabel("NLopt Fraction of Total Iterations")
        ax.set_ylabel("Best Feasible Objective")
        ax.set_title("Hybrid Ratio Sweep")
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"Sensitivity Sweeps - {dataset}", fontsize=13)
    fig.tight_layout()

    fig.savefig(fig_dir / f"sensitivity_sweeps_{dataset}.png")
    fig.savefig(fig_dir / f"sensitivity_sweeps_{dataset}.pdf")
    plt.close(fig)
    print(f"  Saved: sensitivity_sweeps_{dataset}")


def plot_multi_run(traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path):
    """Plot multi-run scenarios with mean line and +/-1 sigma shaded band."""
    groups = detect_multi_run_groups(traces)
    if not groups:
        return

    etalon_obj = get_etalon_best_objective(traces)
    fig, ax = plt.subplots(figsize=(10, 6))

    _add_etalon_reference(ax, etalon_obj)

    for base_name, run_names in sorted(groups.items()):
        dfs = []
        for name in run_names:
            feasible = filter_budget_feasible(traces[name])
            if not feasible.empty:
                bsf = compute_best_so_far(feasible["Objective"])
                df_bsf = feasible[["ElapsedTime(s)"]].copy()
                df_bsf["BestObj"] = bsf.values
                dfs.append(df_bsf)

        if len(dfs) < 2:
            continue

        x_grid, mean_y, std_y = interpolate_to_common_grid(
            dfs, "ElapsedTime(s)", "BestObj"
        )
        style = get_style(base_name)
        color = style.get("color", "gray")
        ax.plot(x_grid, mean_y, label=f"{base_name} (n={len(dfs)})", color=color, linewidth=1.5)
        ax.fill_between(x_grid, mean_y - std_y, mean_y + std_y, color=color, alpha=0.2)

    ax.set_xlabel("Elapsed Time (s)")
    ax.set_ylabel("Best Feasible Objective")
    ax.set_title(f"Multi-Run Best-So-Far (mean +/- 1σ) - {dataset}")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.savefig(fig_dir / f"multi_run_{dataset}.png")
    fig.savefig(fig_dir / f"multi_run_{dataset}.pdf")
    plt.close(fig)
    print(f"  Saved: multi_run_{dataset}")


def plot_ta_compute_time(traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path):
    """Boxplot of TA computation time per scenario."""
    plot_traces = exclude_etalon(traces)
    data = []
    labels = []
    for scenario_name, df in sorted(plot_traces.items()):
        if "TAComputeTime(s)" not in df.columns:
            continue
        vals = df["TAComputeTime(s)"].dropna()
        if vals.empty:
            continue
        data.append(vals.values)
        labels.append(scenario_name)

    if not data:
        print("  Skipping TA compute time plot (no data)")
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(data, labels=labels, patch_artist=True, vert=True)
    for i, box in enumerate(bp["boxes"]):
        style = get_style(labels[i])
        box.set_facecolor(style.get("color", "gray"))
        box.set_alpha(0.6)
    ax.set_ylabel("TA Compute Time (s)")
    ax.set_title(f"Traffic Assignment Compute Time per Scenario - {dataset}")
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()

    fig.savefig(fig_dir / "ta_compute_time.png")
    fig.savefig(fig_dir / "ta_compute_time.pdf")
    plt.close(fig)
    print(f"  Saved: ta_compute_time")


def plot_capacity_heatmap(dataset: str, fig_dir: Path, solutions_dir: Path):
    """Heatmap: links (Y) x scenarios (X), color = capacity delta from lower bound."""
    if not solutions_dir.exists():
        return

    solution_files = sorted(solutions_dir.glob("*_solution.csv"))
    if not solution_files:
        print("  Skipping capacity heatmap (no solution files)")
        return

    scenario_names = []
    deltas = []
    for f in solution_files:
        try:
            df = pd.read_csv(f)
        except Exception:
            continue
        if df.empty or "optimized_capacity" not in df.columns or "lower_bound" not in df.columns:
            continue
        name = f.stem.replace("_solution", "")
        scenario_names.append(name)
        deltas.append((df["optimized_capacity"] - df["lower_bound"]).values)

    if len(scenario_names) < 2:
        print("  Skipping capacity heatmap (need >= 2 solutions)")
        return

    # Build matrix: rows=links, cols=scenarios
    matrix = np.array(deltas).T  # shape: (n_links, n_scenarios)

    # Only show links with non-trivial investment in at least one scenario
    max_per_link = matrix.max(axis=1)
    threshold = max_per_link.max() * 0.01
    active_mask = max_per_link > threshold
    if not active_mask.any():
        return

    matrix_filtered = matrix[active_mask]
    link_indices = np.where(active_mask)[0]

    fig, ax = plt.subplots(figsize=(max(8, len(scenario_names) * 1.2), max(6, len(link_indices) * 0.25)))
    im = ax.imshow(matrix_filtered, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax.set_xticks(range(len(scenario_names)))
    ax.set_xticklabels(scenario_names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Link Index")
    ax.set_yticks(range(len(link_indices)))
    ax.set_yticklabels(link_indices, fontsize=7)
    ax.set_title(f"Capacity Investment Heatmap - {dataset}")
    fig.colorbar(im, ax=ax, label="Capacity - Lower Bound")
    fig.tight_layout()

    fig.savefig(fig_dir / "capacity_heatmap.png")
    fig.savefig(fig_dir / "capacity_heatmap.pdf")
    plt.close(fig)
    print(f"  Saved: capacity_heatmap")


def generate_capacity_comparison(tables_dir: Path, solutions_dir: Path):
    """Merge all solution CSVs into wide-format capacity_comparison.csv."""
    if not solutions_dir.exists():
        return

    solution_files = sorted(solutions_dir.glob("*_solution.csv"))
    if not solution_files:
        return

    base_df = None
    scenario_cols = {}
    for f in solution_files:
        try:
            df = pd.read_csv(f)
        except Exception:
            continue
        if df.empty or "optimized_capacity" not in df.columns:
            continue
        name = f.stem.replace("_solution", "")
        if base_df is None:
            base_cols = ["link_id", "init_node", "term_node", "lower_bound", "upper_bound"]
            base_df = df[[c for c in base_cols if c in df.columns]].copy()
        scenario_cols[f"{name}_capacity"] = df["optimized_capacity"].values
        if "flow" in df.columns:
            scenario_cols[f"{name}_flow"] = df["flow"].values

    if base_df is None:
        return

    for col_name, values in scenario_cols.items():
        base_df[col_name] = values

    out_path = tables_dir / "capacity_comparison.csv"
    base_df.to_csv(out_path, index=False, float_format="%.6f")
    print(f"  Saved: capacity_comparison.csv")


def generate_time_to_target(
    traces: dict[str, pd.DataFrame], tables_dir: Path,
    etalon_obj: float | None,
):
    """For each scenario, find time/evals to reach gap thresholds relative to etalon."""
    if etalon_obj is None or not np.isfinite(etalon_obj) or etalon_obj <= 0:
        return

    thresholds = [0.05, 0.02, 0.01, 0.005]
    plot_traces = exclude_etalon(traces)
    rows = []
    for scenario_name, df in sorted(plot_traces.items()):
        feasible = filter_budget_feasible(df)
        if feasible.empty:
            rows.append({"Scenario": scenario_name})
            continue

        best_envelope = compute_best_so_far(feasible["Objective"])
        gap = (best_envelope - etalon_obj) / etalon_obj
        row = {"Scenario": scenario_name}
        for thr in thresholds:
            reached = gap <= thr
            if reached.any():
                first_idx = reached.idxmax()
                row[f"Time_to_{thr*100:.1f}%"] = feasible.loc[first_idx, "ElapsedTime(s)"]
                row[f"Evals_to_{thr*100:.1f}%"] = int(feasible.index.get_loc(first_idx) + 1)
            else:
                row[f"Time_to_{thr*100:.1f}%"] = np.nan
                row[f"Evals_to_{thr*100:.1f}%"] = np.nan
        rows.append(row)

    if rows:
        pd.DataFrame(rows).to_csv(
            tables_dir / "time_to_target.csv", index=False, float_format="%.4f"
        )
        print(f"  Saved: time_to_target.csv")


def generate_latex_summary(traces: dict[str, pd.DataFrame], tables_dir: Path):
    """Write summary.tex with a LaTeX tabular from the summary data."""
    plot_traces = exclude_etalon(traces)
    rows = []
    for scenario_name, df in sorted(plot_traces.items()):
        feasible = filter_budget_feasible(df)
        total_time = df["ElapsedTime(s)"].max() if not df.empty else 0
        n_evals = len(df)
        if not feasible.empty:
            if "BestFeasibleObjective" in feasible.columns:
                best_vals = feasible["BestFeasibleObjective"].dropna()
                best_obj = best_vals.iloc[-1] if not best_vals.empty else np.nan
            else:
                best_obj = compute_best_so_far(feasible["Objective"]).iloc[-1]
        else:
            best_obj = np.nan
        rows.append((scenario_name, best_obj, total_time, n_evals))

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{CNDP Algorithm Comparison Summary}",
        r"\label{tab:cndp-summary}",
        r"\begin{tabular}{lrrr}",
        r"\toprule",
        r"Scenario & Best Objective & Time (s) & Evaluations \\",
        r"\midrule",
    ]
    for name, obj, t, evals in rows:
        obj_str = f"{obj:.2f}" if np.isfinite(obj) else "---"
        # Escape underscores for LaTeX
        latex_name = name.replace("_", r"\_")
        lines.append(f"{latex_name} & {obj_str} & {t:.1f} & {evals} \\\\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
    ]
    (tables_dir / "summary.tex").write_text("\n".join(lines))
    print(f"  Saved: summary.tex")


def generate_run_dir_readme(
    run_dir: Path, dataset: str, traces: dict[str, pd.DataFrame],
):
    """Regenerate README.md in the run folder with updated metrics from plots."""
    etalon_obj = get_etalon_best_objective(traces)
    plot_traces = exclude_etalon(traces)

    lines = [
        f"# CNDP Comparison Results - {dataset}",
        "",
    ]

    # Read manifest for metadata
    manifest_path = run_dir / "run_manifest.json"
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)
        lines += [
            f"- **Timestamp**: {manifest.get('timestamp', 'N/A')}",
            f"- **Git commit**: {manifest.get('git_commit', 'N/A')}",
            f"- **Platform**: {manifest.get('platform', 'N/A')}",
            f"- **Total wall time**: {manifest.get('total_wall_time_seconds', 0):.1f}s",
            "",
        ]

    if etalon_obj is not None and np.isfinite(etalon_obj):
        lines.append(f"**Etalon reference objective**: {etalon_obj:.4f}")
        lines.append("")

    lines += [
        "## Scenario Results",
        "",
        "| Scenario | Best Objective | Time (s) | Evaluations |",
        "|----------|---------------|----------|-------------|",
    ]
    for scenario_name, df in sorted(plot_traces.items()):
        feasible = filter_budget_feasible(df)
        total_time = df["ElapsedTime(s)"].max() if not df.empty else 0
        n_evals = len(df)
        if not feasible.empty:
            if "BestFeasibleObjective" in feasible.columns:
                best_vals = feasible["BestFeasibleObjective"].dropna()
                best_obj = best_vals.iloc[-1] if not best_vals.empty else np.nan
            else:
                best_obj = compute_best_so_far(feasible["Objective"]).iloc[-1]
            obj_str = f"{best_obj:.4f}" if np.isfinite(best_obj) else "N/A"
        else:
            obj_str = "N/A"
        lines.append(f"| {scenario_name} | {obj_str} | {total_time:.1f} | {n_evals} |")

    lines += [
        "",
        "## Files",
        "",
        "- `figures/` - All plots (PNG + PDF)",
        "- `tables/` - Summary tables (CSV + LaTeX)",
        "- `traces/` - Per-scenario iteration trace data",
        "- `solutions/` - Per-scenario final link capacities",
        "- `metadata/` - Per-scenario metadata JSON files",
        "- `configs/` - Exact config files used",
    ]
    (run_dir / "README.md").write_text("\n".join(lines))


def plot_cross_dataset(
    all_traces: dict[str, dict[str, pd.DataFrame]],
    datasets: list[str],
    fig_dir: Path,
):
    """Cross-dataset comparison plots: avg TA time vs network size, gap at fixed time budgets."""
    # 1. Average TA solve time per dataset
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: Avg TA solve time vs network links
    ds_links = []
    ds_ta_times = []
    for ds in datasets:
        traces = all_traces.get(ds, {})
        ta_times = []
        for df in traces.values():
            if "TAComputeTime(s)" in df.columns:
                ta_times.extend(df["TAComputeTime(s)"].dropna().tolist())
        if ta_times:
            ds_links.append(NETWORK_LINKS.get(ds, 0))
            ds_ta_times.append(np.mean(ta_times))

    if ds_links:
        axes[0].bar(
            [str(l) for l in ds_links],
            ds_ta_times,
            color=["#1f77b4", "#d62728"][:len(ds_links)],
            width=0.5,
        )
        for i, (l, t) in enumerate(zip(ds_links, ds_ta_times)):
            axes[0].text(i, t, f"{t:.2f}s", ha="center", va="bottom", fontsize=10)
        axes[0].set_xlabel("Network Links")
        axes[0].set_ylabel("Avg TA Solve Time (s)")
        axes[0].set_title("TA Solve Time vs Network Size")
        axes[0].grid(True, alpha=0.3, axis="y")

    # Panel 2: Gap-to-etalon at fixed time budgets per method per dataset
    time_budgets = [60, 300, 900]
    common_scenarios = set()
    for ds in datasets:
        common_scenarios.update(
            k for k in all_traces.get(ds, {}).keys() if k != "Etalon"
        )
    common_scenarios = sorted(common_scenarios)

    bar_data = {}  # {scenario: {ds: gap_at_best_budget}}
    for ds in datasets:
        traces = all_traces.get(ds, {})
        etalon_obj = get_etalon_best_objective(traces)
        if etalon_obj is None:
            continue
        for scenario_name in common_scenarios:
            if scenario_name not in traces:
                continue
            df = traces[scenario_name]
            feasible = filter_budget_feasible(df)
            if feasible.empty:
                continue
            best_envelope = compute_best_so_far(feasible["Objective"])
            # Find best gap across time budgets
            best_gap = np.inf
            for tb in time_budgets:
                within = feasible["ElapsedTime(s)"] <= tb
                if within.any():
                    val = best_envelope[within].iloc[-1]
                    gap = (val - etalon_obj) / etalon_obj * 100
                    if gap < best_gap:
                        best_gap = gap
            if np.isfinite(best_gap):
                bar_data.setdefault(scenario_name, {})[ds] = best_gap

    if bar_data:
        scenarios_with_data = [s for s in common_scenarios if s in bar_data]
        x = np.arange(len(scenarios_with_data))
        width = 0.35
        ds_colors = {"SiouxFalls": "#1f77b4", "Anaheim": "#d62728"}
        for i, ds in enumerate(datasets):
            vals = [bar_data.get(s, {}).get(ds, np.nan) for s in scenarios_with_data]
            offset = (i - (len(datasets) - 1) / 2) * width
            axes[1].bar(
                x + offset,
                vals,
                width,
                label=ds,
                color=ds_colors.get(ds, "gray"),
            )
        axes[1].set_xlabel("Scenario")
        axes[1].set_ylabel("Gap to Etalon (%)")
        axes[1].set_title(f"Best Gap at t<={max(time_budgets)}s")
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(scenarios_with_data, rotation=45, ha="right", fontsize=8)
        axes[1].legend()
        axes[1].grid(True, alpha=0.3, axis="y")

    fig.suptitle("Cross-Dataset Comparison", fontsize=13)
    fig.tight_layout()

    fig.savefig(fig_dir / "cross_dataset_comparison.png")
    fig.savefig(fig_dir / "cross_dataset_comparison.pdf")
    plt.close(fig)
    print(f"  Saved: cross_dataset_comparison")


def generate_summary_table(
    traces: dict[str, pd.DataFrame], dataset: str, fig_dir: Path
) -> pd.DataFrame:
    """Generate and save summary table (no etalon)."""
    traces = exclude_etalon(traces)

    rows = []
    for scenario_name, df in sorted(traces.items()):
        feasible = filter_budget_feasible(df)
        total_time = df["ElapsedTime(s)"].max() if not df.empty else 0
        n_evals = len(df)

        if not feasible.empty:
            final_obj = feasible["Objective"].iloc[-1]
            if "BestFeasibleObjective" in feasible.columns:
                best_vals = feasible["BestFeasibleObjective"].dropna()
                best_obj = best_vals.iloc[-1] if not best_vals.empty else final_obj
            else:
                best_obj = compute_best_so_far(feasible["Objective"]).iloc[-1]
        else:
            final_obj = np.nan
            best_obj = np.nan

        rows.append({
            "Scenario": scenario_name,
            "FinalObjective": final_obj,
            "BestFeasibleObj": best_obj,
            "ElapsedTime(s)": total_time,
            "Evaluations": n_evals,
        })

    summary_df = pd.DataFrame(rows)
    summary_path = fig_dir / "summary.csv"
    summary_df.to_csv(summary_path, index=False, float_format="%.6f")

    print(f"\n{'='*90}")
    print(f"SUMMARY TABLE: {dataset}")
    print(f"{'='*90}")
    print(f"{'Scenario':<25} {'BestFeasObj':>15} {'Time(s)':>10} {'Evals':>8}")
    print(f"{'-'*25} {'-'*15} {'-'*10} {'-'*8}")
    for _, row in summary_df.iterrows():
        obj_str = f"{row['BestFeasibleObj']:.4f}" if np.isfinite(row['BestFeasibleObj']) else "N/A"
        print(
            f"{row['Scenario']:<25} {obj_str:>15} "
            f"{row['ElapsedTime(s)']:>10.1f} {row['Evaluations']:>8}"
        )

    print(f"\nSummary saved to: {summary_path}")
    return summary_df


def generate_single_dataset_plots(
    dataset: str, traces: dict[str, pd.DataFrame],
    run_dir: Path | None = None,
):
    """Generate all plots for a single dataset.

    When run_dir is provided, outputs go to run_dir/figures/ and run_dir/tables/.
    """
    if run_dir is not None:
        fig_dir = run_dir / "figures"
        tables_dir = run_dir / "tables"
        solutions_dir = run_dir / "solutions"
        metadata_dir = run_dir / "metadata"
    else:
        fig_dir = PROJECT_ROOT / "performance_results" / dataset / "figures"
        tables_dir = fig_dir
        solutions_dir = None
        metadata_dir = None
    fig_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    non_etalon = exclude_etalon(traces)
    if not non_etalon:
        print(f"ERROR: No non-etalon scenario traces found for {dataset}.")
        return

    etalon_tag = " (etalon present)" if "Etalon" in traces else " (no etalon)"
    print(f"\nPlotting {len(non_etalon)} scenarios for {dataset}{etalon_tag}")

    print("\nGenerating plots...")
    plot_objective_vs_time(traces, dataset, fig_dir)
    plot_best_so_far(traces, dataset, fig_dir)
    plot_optimality_gap_vs_evaluations(traces, dataset, fig_dir)
    plot_budget_utilization(traces, dataset, fig_dir, metadata_dir=metadata_dir)
    plot_convergence_character(traces, dataset, fig_dir)
    plot_ta_compute_time(traces, dataset, fig_dir)
    plot_sensitivity_sweeps(traces, dataset, fig_dir)
    plot_multi_run(traces, dataset, fig_dir)

    if solutions_dir is not None and solutions_dir.exists():
        plot_capacity_heatmap(dataset, fig_dir, solutions_dir)
        generate_capacity_comparison(tables_dir, solutions_dir)

    etalon_obj = get_etalon_best_objective(traces)
    generate_summary_table(traces, dataset, tables_dir)
    generate_time_to_target(traces, tables_dir, etalon_obj)
    generate_latex_summary(traces, tables_dir)

    if run_dir is not None:
        generate_run_dir_readme(run_dir, dataset, traces)

    print(f"\nAll figures saved to: {fig_dir}")


def main():
    parser = argparse.ArgumentParser(description="Plot CNDP comparison results")
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
        help="Comma-separated datasets for cross-dataset comparison (e.g. SiouxFalls,Anaheim)",
    )
    parser.add_argument(
        "--format",
        choices=["png", "pdf", "both"],
        default="both",
        help="Output format (default: both)",
    )
    parser.add_argument(
        "--trace-dir",
        type=str,
        default=None,
        help="Override trace directory (default: performance_results/<dataset>)",
    )
    parser.add_argument(
        "--run-dir",
        type=str,
        default=None,
        help="Self-contained run folder (loads from traces/, saves to figures/ and tables/)",
    )

    args = parser.parse_args()

    # Run-dir mode: self-contained folder
    if args.run_dir:
        run_dir = Path(args.run_dir)
        dataset = args.dataset or "SiouxFalls"
        print(f"Loading traces from run folder: {run_dir}")
        traces = load_trace_files_from_run_dir(run_dir)
        if not traces:
            print("ERROR: No trace files found in run folder.")
            sys.exit(1)
        print(f"Found {len(traces)} scenario traces: {', '.join(sorted(traces.keys()))}")
        generate_single_dataset_plots(dataset, traces, run_dir=run_dir)
        return

    # Cross-dataset mode
    if args.datasets:
        dataset_list = [d.strip() for d in args.datasets.split(",")]
        all_traces = {}
        for ds in dataset_list:
            print(f"Loading trace files for {ds}...")
            traces = load_trace_files(ds)
            if traces:
                print(f"  Found {len(traces)} traces: {', '.join(sorted(traces.keys()))}")
                all_traces[ds] = traces
                generate_single_dataset_plots(ds, traces)
            else:
                print(f"  WARNING: No traces found for {ds}")

        if len(all_traces) >= 2:
            fig_dir = PROJECT_ROOT / "performance_results" / "figures"
            fig_dir.mkdir(parents=True, exist_ok=True)
            print("\nGenerating cross-dataset comparison plots...")
            plot_cross_dataset(all_traces, dataset_list, fig_dir)
        return

    # Single dataset mode
    dataset = args.dataset or "SiouxFalls"
    print(f"Loading trace files for {dataset}...")
    traces = load_trace_files(dataset)

    if not traces:
        print("ERROR: No trace files found. Run experiments first:")
        print(f"  python scripts/run_cndp_comparison.py --dataset {dataset}")
        sys.exit(1)

    print(f"Found {len(traces)} scenario traces: {', '.join(sorted(traces.keys()))}")

    generate_single_dataset_plots(dataset, traces)


if __name__ == "__main__":
    main()
