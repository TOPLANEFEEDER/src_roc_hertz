#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_sphere4_ftmag.py

Read Rockable conf files (e.g. conf100, conf101, ...) and plot tangential-force
magnitude |ft| from the Interactions block.

By default:
- track particle Sphere4
- find all interactions involving that particle
- compute |ft| = sqrt(ft1^2 + ft2^2 + ft3^2)
- sum all |ft| involving that particle in each conf

Usage examples
--------------
1) Default: Sphere4, x-axis = time, y = sum of |ft| over all contacts of Sphere4
   python3 plot_sphere4_ftmag.py --start 100 --end 200

2) Use conf index as x-axis
   python3 plot_sphere4_ftmag.py --start 100 --end 200 --xmode iconf

3) Another particle
   python3 plot_sphere4_ftmag.py --start 100 --end 200 --particle Sphere3

4) Output CSV
   python3 plot_sphere4_ftmag.py --start 100 --end 200 --csv ft_vs_time.csv

5) Plot each contact separately instead of total sum
   python3 plot_sphere4_ftmag.py --start 100 --end 200 --mode each
"""

import argparse
import math
from pathlib import Path
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot |ft| from Rockable conf Interactions block."
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
    p.add_argument("--outfile", type=str, default="sphere4_ftmag_vs_time.png", help="Output figure filename")
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

    Example:
      Particles 4
      Sphere1 ...
      Sphere2 ...
      Sphere3 ...
      Sphere4 ...

    Then:
      Sphere1 -> 0
      Sphere2 -> 1
      Sphere3 -> 2
      Sphere4 -> 3
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


def read_interaction_ftmagnitudes(conf_path: Path, particle_index: int):
    """
    Read all interactions involving the target particle index and return a list of |ft|.

    Interaction columns:
      0  i
      1  j
      2  type
      3  isub
      4  jsub
      5  nx
      6  ny
      7  nz
      8  dn
      9  ix
      10 iy
      11 iz
      12 vx
      13 vy
      14 vz
      15 fn
      16 ft1
      17 ft2
      18 ft3
      19 mom1
      20 mom2
      21 mom3
      22 vd
    """
    in_interactions = False
    ftmags = []

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
                ft1 = float(parts[28])
                ft2 = float(parts[29])
                ft3 = float(parts[30])
                ftmag = math.sqrt(ft1 * ft1 + ft2 * ft2 + ft3 * ft3)
                ftmags.append(ftmag)

    return ftmags


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
            ftmags = read_interaction_ftmagnitudes(conf_path, pidx)
        except Exception as e:
            print(f"[warning] Failed to parse {conf_path}: {e}")
            continue

        if not ftmags:
            print(f"[warning] No interaction involving {args.particle} in {conf_path}")
            continue

        x = t if args.xmode == "time" else iconf

        xvals.append(x)
        conf_ids.append(i)

        if args.mode == "each":
            each_vals.append(ftmags)
        else:
            y = reduce_values(ftmags, args.mode)
            yvals.append(y)

    if not xvals:
        raise RuntimeError("No valid data found. Please check your conf range and files.")

    # Optional CSV output
    if args.csv:
        with open(args.csv, "w", encoding="utf-8") as f:
            if args.mode == "each":
                if args.xmode == "time":
                    f.write("conf_id,time,contact_id,ft_mag\n")
                    for cid, x, vals in zip(conf_ids, xvals, each_vals):
                        for k, v in enumerate(vals):
                            f.write(f"{cid},{x:.16e},{k},{v:.16e}\n")
                else:
                    f.write("conf_id,iconf,contact_id,ft_mag\n")
                    for cid, x, vals in zip(conf_ids, xvals, each_vals):
                        for k, v in enumerate(vals):
                            f.write(f"{cid},{x},{k},{v:.16e}\n")
            else:
                if args.xmode == "time":
                    f.write(f"conf_id,time,ft_mag_{args.mode}\n")
                    for cid, x, y in zip(conf_ids, xvals, yvals):
                        f.write(f"{cid},{x:.16e},{y:.16e}\n")
                else:
                    f.write(f"conf_id,iconf,ft_mag_{args.mode}\n")
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
        ylabel = f"|ft| of each contact involving {args.particle}"
        title_y = "each |ft|"
    else:
        plt.plot(xvals, yvals, marker="o", linewidth=1)
        ylabel = f"{args.mode} |ft| involving {args.particle}"
        title_y = f"{args.mode} |ft|"

    if args.xmode == "time":
        plt.xlabel("time t")
        title_x = "time"
    else:
        plt.xlabel("iconf")
        title_x = "iconf"

    plt.ylabel(ylabel)
    #plt.ylim(0.9,1.1)
    plt.title(f"{args.particle} {title_y} vs {title_x}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=300)
    print(f"[info] Figure written to: {args.outfile}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()