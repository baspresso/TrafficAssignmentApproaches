"""
Generate capacity constraints CSV for a transportation network.

Usage:
    python generate_constraints.py <net_name> [options]

Example:
    python generate_constraints.py SiouxFalls
    python generate_constraints.py SiouxFalls --lower-multiplier 0.8 --upper-multiplier 1.2
    python generate_constraints.py SiouxFalls --investment-cost-param 100
"""

import argparse
import os
import sys
import pandas as pd


def generate_constraints(
    net_name: str,
    data_root: str,
    lower_multiplier: float = 0.9,
    upper_multiplier: float = 1.1,
    investment_cost_param: float = 0.0,
) -> pd.DataFrame:
    net_dir = os.path.join(data_root, net_name)
    net_file = os.path.join(net_dir, f"{net_name}_net.csv")

    if not os.path.isfile(net_file):
        print(f"Error: network file not found: {net_file}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(net_file, comment="#")

    constraints = df[["init_node", "term_node"]].copy()
    constraints["lower_bound"] = df["capacity"] * lower_multiplier
    constraints["upper_bound"] = df["capacity"] * upper_multiplier
    constraints["investment_cost_param"] = investment_cost_param

    out_file = os.path.join(net_dir, f"{net_name}_constraints.csv")
    constraints.to_csv(out_file, index=False)
    print(f"Saved: {out_file}")

    return constraints


def main():
    parser = argparse.ArgumentParser(
        description="Generate capacity constraints CSV for a transportation network."
    )
    parser.add_argument("net_name", help="Network name (e.g. SiouxFalls)")
    parser.add_argument(
        "--data-root",
        default=os.path.join(os.path.dirname(__file__), "..", "data", "TransportationNetworks"),
        help="Path to the TransportationNetworks data directory (default: ../data/TransportationNetworks relative to this script)",
    )
    parser.add_argument(
        "--lower-multiplier",
        type=float,
        default=0.9,
        help="Multiplier applied to capacity for the lower bound (default: 0.9)",
    )
    parser.add_argument(
        "--upper-multiplier",
        type=float,
        default=1.1,
        help="Multiplier applied to capacity for the upper bound (default: 1.1)",
    )
    parser.add_argument(
        "--investment-cost-param",
        type=float,
        default=0.0,
        help="Investment cost parameter assigned to all links (default: 0)",
    )

    args = parser.parse_args()
    data_root = os.path.abspath(args.data_root)

    df = generate_constraints(
        net_name=args.net_name,
        data_root=data_root,
        lower_multiplier=args.lower_multiplier,
        upper_multiplier=args.upper_multiplier,
        investment_cost_param=args.investment_cost_param,
    )

    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
