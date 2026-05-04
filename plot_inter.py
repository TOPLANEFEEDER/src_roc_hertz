#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_sphere4_interaction_value.py

Read Rockable conf files (e.g. conf100, conf101, ...) and plot a chosen value
from the Interactions block for contacts involving a target particle.

Supports:
1) direct column index with --icol
2) common derived quantities with --expr, e.g. ftmag

Examples
--------
1) Plot |ft| vs time for Sphere4
   python3 plot_sphere4_interaction_value.py --start 100 --end 200 --expr ftmag

2) Plot Interaction column 22 vs time
   python3 plot_sphere4_interaction_value.py --start 100 --end 200 --icol 22

3) Plot each contact separately
   python3 plot_sphere4_interaction_value.py --start 100 --end 200 --expr ftmag --mode each

4) Use conf index as x-axis
   python3 plot_sphere4_interaction_value.py --start 100 --end 200 --icol 15 --xmode iconf

5) Output CSV
   python3 plot_sphere4_interaction_value.py --start 100 --end 200 --expr ftmag --csv out.csv
"""

import argparse
import math
from pathlib import Path
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot a chosen Interaction value from Rockable conf files."
    )
    p.add_argument("--start", type=int, required=True, help="Start conf index, e.g. 100 for conf100")
    p.add_argument("--end", type=int, required=True, help="End conf index, e.g. 200 for conf200")
    p.add_argument("--step", type=int, default=1, help="Step between conf indices")
    p.add_argument("--directory", type=str, default=".", help="Directory containing conf files")
    p.add_argument("--prefix", type=str, default="conf", help="Filename prefix, default: conf")
    p.add_argument("--particle", type=str, default="Sphere4", help="Particle name to track, default: Sphere4")

    p.add_argument(
        "--xmode",
        type=str,
        choices=["time", "iconf"],
        default="time",
        help="x-axis mode: time (default) or iconf"
    )

    p.add_argument(
        "--mode",
        type=str,
        choices=["sum", "mean", "max", "each"],
        default="sum",
        help="How to reduce multiple contacts involving the particle: sum/mean/max/each"
    )

    group = p.add_mutually_exclusive_group()
    group.add_argument(
        "--icol",
        type=int,
        default=None,
        help="Interaction column index to plot directly (0-based after split)"
    )
    group.add_argument(
        "--expr",
        type=str,
        default="ftmag",
        choices=["ftmag", "fn", "dn", "vd", "vx", "vy", "vz", "nx", "ny", "nz"],
        help="Named quantity to compute from interaction columns"
    )

    p.add_argument("--outfile", type=str, default="interaction_value_vs_time.png", help="Output figure filename")
    p.add_argument("--csv", type=str, default="", help="Optional CSV output filename")
    p.add_argument("--show", action="store_true", help="Show the figure interactively")
    return p.parse_args()


def read_conf_time(conf_path: Path):
    with conf_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s.startswith("t "):
                parts = s.split()
                if len(parts) >= 2:
                    return float(parts[1])
    raise ValueError(f"Could not find time 't' in {conf_path}")


def read_conf_iconf(conf_path: Path):
    with conf_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s.startswith("iconf "):
                parts = s.split()
                if len(parts) >= 2:
                    return int(parts[1])
    raise ValueError(f"Could not find 'iconf' in {conf_path}")


def read_particle_index(conf_path: Path, particle_name: str):
    """
    Read the Particles block and return the 0-based particle index used by Interactions.
    """
    in_particles = False
    idx = -1

    with conf_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()

            if s.startswith("Particles "):
                in_particles = True
                idx = -1
                continue

            if not in_particles:
                continue

            if not s or s.startswith("#"):
                continue

            if s.startswith("Interactions ") or s.startswith("Interfaces "):
                break

            parts = s.split()
            if len(parts) < 1:
                continue

            idx += 1
            name = parts[0]
            if name == particle_name:
                return idx

    raise ValueError(f"Could not find particle '{particle_name}' in {conf_path}")


def compute_interaction_value(parts, icol=None, expr="ftmag"):
    """
    parts = split() result of one interaction line

    Common non-periodic interaction columns assumed:
      0  i
      1  j
      2  type
      3  isub
      4  jsub
      5  nx
      6  ny
      7  nz
      8  dn
      9  pos.x
      10 pos.y
      11 pos.z
      12 vel.x
      13 vel.y
      14 vel.z
      15 fn
      16 ft.x
      17 ft.y
      18 ft.z
      19 mom.x
      20 mom.y
      21 mom.z
      22 vd

    If your file format differs, you can still use --icol directly.
    """
    if icol is not None:
        if icol < 0 or icol >= len(parts):
            raise ValueError(f"Requested column index {icol} out of range for interaction line with {len(parts)} columns")
        return float(parts[icol])

    if expr == "ftmag":
        ft1 = float(parts[16])
        ft2 = float(parts[17])
        ft3 = float(parts[18])
        return math.sqrt(ft1 * ft1 + ft2 * ft2 + ft3 * ft3)

    if expr == "fn":
        return float(parts[15])

    if expr == "dn":
        return float(parts[8])

    if expr == "vd":
        return float(parts[22])

    if expr == "vx":
        return float(parts[12])

    if expr == "vy":
        return float(parts[13])

    if expr == "vz":
        return float(parts[14])

    if expr == "nx":
        return float(parts[5])

    if expr == "ny":
        return float(parts[6])

    if expr == "nz":
        return float(parts[7])

    raise ValueError(f"Unsupported expr: {expr}")


def read_interaction_values(conf_path: Path, particle_index: int, icol=None, expr="ftmag"):
    """
    Read all interactions involving the target particle index and return a list of chosen values.
    """
    in_interactions = False
    values = []

    with conf_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()

            if s.startswith("Interactions "):
                in_interactions = True
                continue

            if not in_interactions:
                continue

            if not s or s.startswith("#"):
                continue

            if s.startswith("Interfaces "):
                break

            parts = s.split()
            if len(parts) < 23:
                continue

            i_idx = int(parts[0])
            j_idx = int(parts[1])

            if i_idx == particle_index or j_idx == particle_index:
                val = compute_interaction_value(parts, icol=icol, expr=expr)
                values.append(val)

    return values


def reduce_values(values, mode):
    if not values:
        return None
    if mode == "sum":
        return sum(values)
    if mode == "mean":
        return sum(values) / len(values)
    if mode == "max":
        return max(values)
    raise ValueError(f"Unsupported mode: {mode}")


def ylabel_from_args(args):
    if args.icol is not None:
        return f"Interaction column[{args.icol}]"
    return args.expr


def title_from_args(args):
    if args.icol is not None:
        return f"column[{args.icol}]"
    return args.expr


def main():
    args = parse_args()

    if args.step <= 0:
        raise ValueError("--step must be > 0")

    directory = Path(args.directory)

    xvals = []
    yvals = []
    conf_ids = []
    each_vals = []

    for i in range(args.start, args.end + 1, args.step):
        conf_path = directory / f"{args.prefix}{i}"

        if not conf_path.exists():
            print(f"[warning] File not found, skipped: {conf_path}")
            continue

        try:
            t = read_conf_time(conf_path)
            iconf = read_conf_iconf(conf_path)
            pidx = read_particle_index(conf_path, args.particle)
            vals = read_interaction_values(conf_path, pidx, icol=args.icol, expr=args.expr)
        except Exception as e:
            print(f"[warning] Failed to parse {conf_path}: {e}")
            continue

        if not vals:
            print(f"[warning] No interaction involving {args.particle} in {conf_path}")
            continue

        x = t if args.xmode == "time" else iconf

        xvals.append(x)
        conf_ids.append(i)

        if args.mode == "each":
            each_vals.append(vals)
        else:
            y = reduce_values(vals, args.mode)
            yvals.append(y)

    if not xvals:
        raise RuntimeError("No valid data found. Please check your conf range and files.")

    # CSV output
    if args.csv:
        with open(args.csv, "w", encoding="utf-8") as f:
            colname = title_from_args(args)

            if args.mode == "each":
                if args.xmode == "time":
                    f.write(f"conf_id,time,contact_id,{colname}\n")
                    for cid, x, vals in zip(conf_ids, xvals, each_vals):
                        for k, v in enumerate(vals):
                            f.write(f"{cid},{x:.16e},{k},{v:.16e}\n")
                else:
                    f.write(f"conf_id,iconf,contact_id,{colname}\n")
                    for cid, x, vals in zip(conf_ids, xvals, each_vals):
                        for k, v in enumerate(vals):
                            f.write(f"{cid},{x},{k},{v:.16e}\n")
            else:
                if args.xmode == "time":
                    f.write(f"conf_id,time,{colname}_{args.mode}\n")
                    for cid, x, y in zip(conf_ids, xvals, yvals):
                        f.write(f"{cid},{x:.16e},{y:.16e}\n")
                else:
                    f.write(f"conf_id,iconf,{colname}_{args.mode}\n")
                    for cid, x, y in zip(conf_ids, xvals, yvals):
                        f.write(f"{cid},{x},{y:.16e}\n")

        print(f"[info] CSV written to: {args.csv}")

    # Plot
    plt.figure(figsize=(8, 5))

    if args.mode == "each":
        max_contacts = max(len(v) for v in each_vals)
        for k in range(max_contacts):
            xs = []
            ys = []
            for x, vals in zip(xvals, each_vals):
                if k < len(vals):
                    xs.append(x)
                    ys.append(vals[k])
            if xs:
                plt.plot(xs, ys, marker="o", linewidth=1, label=f"contact {k}")
        plt.legend()
        ylabel = f"{ylabel_from_args(args)} of each contact involving {args.particle}"
        title_y = f"each {title_from_args(args)}"
    else:
        plt.plot(xvals, yvals, marker="o", linewidth=1)
        ylabel = f"{args.mode} {ylabel_from_args(args)} involving {args.particle}"
        title_y = f"{args.mode} {title_from_args(args)}"

    if args.xmode == "time":
        plt.xlabel("time t")
        title_x = "time"
    else:
        plt.xlabel("iconf")
        title_x = "iconf"

    plt.ylabel(ylabel)
    plt.title(f"{args.particle} {title_y} vs {title_x}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=300)
    print(f"[info] Figure written to: {args.outfile}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()