"""
Create a new transportation network dataset by scaling OD demand.

Reads <data_root>/<net_name>/<net_name>_trips.csv, multiplies every demand entry
by the given multiplier, and writes a new self-contained dataset directory at
<data_root>/<net_name>_demand_<mult>x/ containing the scaled _trips.csv plus
verbatim copies of _net.csv and (if present) _constraints.csv. The new directory
can be selected directly via `dataset = "<net_name>_demand_<mult>x"` in a TOML
config.

Usage:
    python scripts/scale_demand.py <net_name> <multiplier> [--data-root PATH]

Example:
    python scripts/scale_demand.py SiouxFalls 1.5
    python scripts/scale_demand.py SiouxFalls 2
"""

import argparse
import math
import os
import shutil
import sys


def format_multiplier_label(m: float) -> str:
    if m == int(m):
        return f"{int(m)}x"
    return f"{m:g}x"


def format_value(v: float) -> str:
    rounded = round(v)
    if abs(v - rounded) < 1e-9:
        return str(rounded)
    return f"{v:.10g}"


def scale_trips(src_path: str, dst_path: str, multiplier: float):
    """Stream-scale a _trips.csv file. Returns (original_total, scaled_total)."""
    orig_total = 0.0
    scaled_total = 0.0
    with open(src_path, "r") as src, open(dst_path, "w") as dst:
        header = src.readline()
        if not header.startswith("#"):
            raise ValueError(f"Unexpected trips header (expected '# NUMBER OF ZONES:...'): {header!r}")
        dst.write(header)
        for line in src:
            stripped = line.strip()
            if not stripped:
                dst.write(line)
                continue
            values = [float(x) for x in stripped.split(",")]
            scaled = [v * multiplier for v in values]
            orig_total += sum(values)
            scaled_total += sum(scaled)
            dst.write(",".join(format_value(v) for v in scaled) + "\n")
    return orig_total, scaled_total


def scale_demand(net_name: str, multiplier: float, data_root: str) -> str:
    if not math.isfinite(multiplier) or multiplier <= 0:
        print(f"Error: multiplier must be positive and finite, got {multiplier}", file=sys.stderr)
        sys.exit(1)

    src_dir = os.path.join(data_root, net_name)
    src_net = os.path.join(src_dir, f"{net_name}_net.csv")
    src_trips = os.path.join(src_dir, f"{net_name}_trips.csv")
    src_constraints = os.path.join(src_dir, f"{net_name}_constraints.csv")

    if not os.path.isfile(src_net):
        print(f"Error: missing source file: {src_net}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(src_trips):
        print(f"Error: missing source file: {src_trips}", file=sys.stderr)
        sys.exit(1)

    label = format_multiplier_label(multiplier)
    new_name = f"{net_name}_demand_{label}"
    dst_dir = os.path.join(data_root, new_name)

    if os.path.exists(dst_dir):
        print(f"Error: target dataset directory already exists: {dst_dir}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(dst_dir)
    dst_net = os.path.join(dst_dir, f"{new_name}_net.csv")
    dst_trips = os.path.join(dst_dir, f"{new_name}_trips.csv")
    dst_constraints = os.path.join(dst_dir, f"{new_name}_constraints.csv")

    shutil.copyfile(src_net, dst_net)
    orig_total, scaled_total = scale_trips(src_trips, dst_trips, multiplier)
    written = [dst_net, dst_trips]
    if os.path.isfile(src_constraints):
        shutil.copyfile(src_constraints, dst_constraints)
        written.append(dst_constraints)

    print(f"Created dataset: {dst_dir}")
    print(f"  Multiplier:        {multiplier}")
    print(f"  Original OD flow:  {orig_total:g}")
    print(f"  Scaled OD flow:    {scaled_total:g}")
    print(f"  Files written:")
    for f in written:
        print(f"    {f}")
    print(f"  Use in TOML:       dataset = \"{new_name}\"")
    return dst_dir


def main():
    parser = argparse.ArgumentParser(
        description="Create a new transportation network dataset by scaling OD demand by a constant multiplier."
    )
    parser.add_argument("net_name", help="Source network name (e.g. SiouxFalls)")
    parser.add_argument("multiplier", type=float, help="Positive multiplier applied to every OD demand entry")
    parser.add_argument(
        "--data-root",
        default=os.path.join(os.path.dirname(__file__), "..", "data", "TransportationNetworks"),
        help="Path to the TransportationNetworks data directory (default: ../data/TransportationNetworks relative to this script)",
    )
    args = parser.parse_args()
    data_root = os.path.abspath(args.data_root)
    scale_demand(net_name=args.net_name, multiplier=args.multiplier, data_root=data_root)


if __name__ == "__main__":
    main()
