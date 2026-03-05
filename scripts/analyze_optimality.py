"""
Analyze CND optimality results: classify links by bound type, compute lambda,
compare against the reference value, and check optimality conditions.

Reads  <net_name>/<net_name>_optimality.csv
Writes <net_name>/<net_name>_optimality.xlsx

Usage:
    python analyze_optimality.py <net_name> [options]

Example:
    python analyze_optimality.py SiouxFalls
    python analyze_optimality.py SiouxFalls --bound-threshold 0.05 --comparison-threshold 0.05
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd


def classify_bound_type(df: pd.DataFrame, bound_threshold: float) -> pd.DataFrame:
    """Classify each link as Lower / Upper / Middle based on how close its
    optimized capacity is to the constraint bounds."""
    upper_mask = (df["capacity"] - df["upper_bound"]).abs() <= bound_threshold
    lower_mask = (df["capacity"] - df["lower_bound"]).abs() <= bound_threshold

    df["bound_type"] = np.select(
        [upper_mask, lower_mask],
        ["Upper", "Lower"],
        default="Middle",
    )
    return df


def compute_lambda(df: pd.DataFrame) -> pd.DataFrame:
    """Compute lambda for Middle links (condition_result - 1)."""
    middle_mask = df["bound_type"] == "Middle"
    if not middle_mask.any():
        print(
            "Warning: no Middle links found; lambda cannot be computed.",
            file=sys.stderr,
        )
        df["lambda"] = np.nan
        return df

    df["lambda"] = np.nan
    df.loc[middle_mask, "lambda"] = df.loc[middle_mask, "condition_result"] - 1
    return df


def compute_comparison_result(
    df: pd.DataFrame, reference_value: float, comparison_threshold: float
) -> pd.DataFrame:
    """Compare condition_result against the reference value within a threshold."""
    diff = df["condition_result"] - reference_value
    df["comparison_result"] = np.select(
        [diff.abs() <= comparison_threshold, diff > comparison_threshold, diff < -comparison_threshold],
        ["Equal", "More", "Less"],
        default="Equal",
    )
    return df


def check_optimality_condition(row) -> bool:
    """Return True if the (bound_type, comparison_result) pair satisfies the
    KKT-style optimality condition."""
    bound = row["bound_type"]
    comparison = row["comparison_result"]
    return (
        (bound == "Lower" and comparison == "Less")
        or (bound == "Middle" and comparison == "Equal")
        or (bound == "Upper" and comparison == "More")
    )


def analyze_optimality(
    net_name: str,
    data_root: str,
    bound_threshold: float = 0.01,
    comparison_threshold: float = 0.1,
) -> pd.DataFrame:
    net_dir = os.path.join(data_root, net_name)
    in_file = os.path.join(net_dir, f"{net_name}_optimality.csv")

    if not os.path.isfile(in_file):
        print(f"Error: optimality file not found: {in_file}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(in_file, comment="#")

    df = classify_bound_type(df, bound_threshold)
    df = compute_lambda(df)

    middle_lambdas = df.loc[df["bound_type"] == "Middle", "lambda"].dropna()
    if len(middle_lambdas) == 0:
        print("Error: no lambda values for Middle links.", file=sys.stderr)
        sys.exit(1)

    mean_lambda = middle_lambdas.mean()
    reference_value = 1 + mean_lambda
    print(f"Mean lambda (Middle links): {mean_lambda:.6f}")
    print(f"Reference value:            {reference_value:.8f}")
    print(f"Comparison threshold:       ±{comparison_threshold}")

    df["reference_value"] = reference_value
    df = compute_comparison_result(df, reference_value, comparison_threshold)
    df["condition_check"] = df.apply(check_optimality_condition, axis=1)

    # Reorder columns for readability
    preferred_order = [
        "init_node", "term_node", "capacity",
        "lower_bound", "upper_bound",
        "bound_type", "condition_result",
        "reference_value", "lambda",
        "comparison_result", "condition_check",
    ]
    existing = [c for c in preferred_order if c in df.columns]
    remaining = [c for c in df.columns if c not in existing]
    df = df[existing + remaining]

    out_file = os.path.join(net_dir, f"{net_name}_optimality.xlsx")
    df.to_excel(out_file, index=False)
    print(f"Saved: {out_file}")

    n_ok = df["condition_check"].sum()
    print(f"Optimality conditions met: {n_ok}/{len(df)} links")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Analyze CND optimality results and write an annotated Excel file."
    )
    parser.add_argument("net_name", help="Network name (e.g. SiouxFalls)")
    parser.add_argument(
        "--data-root",
        default=os.path.join(os.path.dirname(__file__), "..", "data", "TransportationNetworks"),
        help="Path to the TransportationNetworks data directory (default: ../data/TransportationNetworks relative to this script)",
    )
    parser.add_argument(
        "--bound-threshold",
        type=float,
        default=0.01,
        help="Tolerance for classifying a link as being at its lower/upper bound (default: 0.01)",
    )
    parser.add_argument(
        "--comparison-threshold",
        type=float,
        default=0.1,
        help="Tolerance for declaring condition_result Equal to the reference value (default: 0.1)",
    )

    args = parser.parse_args()
    data_root = os.path.abspath(args.data_root)

    df = analyze_optimality(
        net_name=args.net_name,
        data_root=data_root,
        bound_threshold=args.bound_threshold,
        comparison_threshold=args.comparison_threshold,
    )

    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
