"""
plot_optimum_temperatures.py

Overlay plot of the T(t) trajectory for every archived run_type=="optimum"
entry in raw_data/master_index.csv, each curve independently shifted so
its OWN temperature peak sits at t=0, then windowed to a user-given time
range around that peak.

UNITS: raw_data/run_<id>.csv's "t" and "T" columns are exactly what the
C++ solver wrote after ControlParameters::dimensionalize() -- verified
directly against control_parameters.cpp's dimensionalize(): `t *= t_ref`
with t_ref = 1e-9 [s] (control_parameters.h), `x[2] *= T_ref` for
temperature, no other conversion applied downstream (single_run.py/
archive.py never rescale these). So on disk: t is in SECONDS, T is in
KELVIN. This script converts the shifted time axis to MICROSECONDS
(matching --time-range's unit); T is plotted in Kelvin as-is, unconverted.

COLOR CODING (by energy_demand, read straight from master_index.csv, not
recomputed) -- up to four categories, each mutually exclusive (checked in
this priority order, so a run can't land in more than one):
    - the single run with the global MINIMUM energy_demand -- "optimum",
      red, thick
    - IF --target is given: the single run whose OWN energy_demand is
      numerically closest to that user-specified value -- "target", blue,
      thick (skipped entirely if --target isn't passed; if the closest
      match happens to also be the optimum run, "optimum" wins -- it's
      the same curve either way, just drawn red instead of blue)
    - every remaining run within --within-pct percent of the minimum
      (default 10%, i.e. energy_demand <= min * 1.10) -- "near_optimum",
      orange, thin
    - everything else -- "other", grey, thin and drawn with low alpha so
      a plot with thousands of curves shows density rather than becoming
      a solid grey block
Draw order is other -> near_optimum -> target -> optimum (each later
category on top), so neither highlighted curve is ever buried under the
grey/orange mass. The legend has one entry per category actually plotted
(3 or 4), not one per curve -- with potentially thousands of runs, a
per-curve legend would be both unreadable and slow to render.

Usage:
    python plot_optimum_temperatures.py --time-range 10
    python plot_optimum_temperatures.py --time-range 10 --within-pct 5
    python plot_optimum_temperatures.py --time-range 10 --target 1500
    python plot_optimum_temperatures.py --time-range 10 --no-show
    python plot_optimum_temperatures.py --time-range 10 --master-index raw_data/master_index.csv --output processed_data/optimum_temperatures.png

Saves a PNG to processed_data/optimum_temperatures.png by default (pass
--output to change) and also opens an interactive matplotlib window
unless --no-show is given.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MASTER_INDEX_PATH = Path("raw_data") / "master_index.csv"
OUTPUT_DIR = Path("processed_data")

STYLE = {
    "other":        dict(color="0.6",   linewidth=0.5, alpha=0.25, zorder=1),
    "near_optimum": dict(color="orange", linewidth=0.9, alpha=0.6,  zorder=2),
    "target":       dict(color="blue",   linewidth=2.5, alpha=1.0,  zorder=3),
    "optimum":      dict(color="red",    linewidth=2.5, alpha=1.0,  zorder=4),
}


def load_optimum_rows(master_index_path: Path) -> list[dict]:
    """Reads master_index.csv and returns every row with run_type=='optimum',
    as plain dicts (all values still strings, same as csv.DictReader gives)."""
    if not master_index_path.exists():
        raise SystemExit(f"{master_index_path} does not exist -- run archive.py first.")
    with open(master_index_path, newline="") as f:
        rows = [row for row in csv.DictReader(f) if row.get("run_type") == "optimum"]
    if not rows:
        raise SystemExit(f"No run_type=='optimum' rows found in {master_index_path}.")
    return rows


def read_trajectory(raw_data_path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """Returns (t_seconds, T_kelvin) arrays read from one run's raw CSV, or
    None if the file is missing, empty (e.g. a run whose confirmation run
    hit a fatal preprocessing failure -- archive_single_result() writes an
    empty file in that case, see archive.py), or doesn't carry a usable
    "t"/"T" pair -- all treated the same way: this run can't be plotted,
    the caller counts and reports how many were skipped."""
    if not raw_data_path.exists():
        return None
    try:
        with open(raw_data_path, newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or "t" not in reader.fieldnames or "T" not in reader.fieldnames:
                return None
            t_vals: list[float] = []
            T_vals: list[float] = []
            for r in reader:
                t_vals.append(float(r["t"]))
                T_vals.append(float(r["T"]))
    except (OSError, ValueError):
        return None
    if len(t_vals) < 2:
        return None
    return np.array(t_vals), np.array(T_vals)


def main():
    parser = argparse.ArgumentParser(
        description="Overlay-plot every archived optimum run's T(t) curve, each shifted "
                    "so its own peak sits at t=0, color-coded by energy_demand."
    )
    parser.add_argument("--master-index", type=Path, default=MASTER_INDEX_PATH,
                         help=f"path to master_index.csv (default: {MASTER_INDEX_PATH})")
    parser.add_argument("--time-range", type=float, required=True,
                         help="microseconds shown each side of every curve's own peak, "
                              "e.g. 10 -> window is [-10, +10] us around t=0")
    parser.add_argument("--within-pct", type=float, default=10.0,
                         help="a run within this percent of the global minimum "
                              "energy_demand is drawn orange (default: 10)")
    parser.add_argument("--target", type=float, default=None,
                         help="a target energy_demand [GJ/t] -- the run whose OWN "
                              "energy_demand is numerically closest to this value is drawn "
                              "as an extra thick blue curve. Omit to skip this entirely.")
    parser.add_argument("--output", type=Path, default=None,
                         help=f"PNG output path (default: {OUTPUT_DIR / 'optimum_temperatures.png'})")
    parser.add_argument("--no-show", action="store_true",
                         help="skip the interactive window, just save the PNG")
    args = parser.parse_args()

    rows = load_optimum_rows(args.master_index)

    usable = []
    skipped_bad_energy = 0
    for row in rows:
        try:
            row["_energy_demand"] = float(row["energy_demand"])
        except (KeyError, ValueError):
            skipped_bad_energy += 1
            continue
        usable.append(row)
    if skipped_bad_energy:
        print(f"Skipping {skipped_bad_energy} optimum run(s): missing/invalid energy_demand.")
    if not usable:
        raise SystemExit("No optimum run has a usable energy_demand value -- nothing to plot.")

    min_energy = min(row["_energy_demand"] for row in usable)
    threshold = min_energy * (1.0 + args.within_pct / 100.0)

    # --target: pick the SINGLE closest-matching run up front (by identity,
    # via `is` below) -- not recomputed per-row, so there's exactly one
    # target run even if several rows happen to tie on distance. Picked
    # from `usable` regardless of whether it turns out to be plottable;
    # if it isn't (missing/unreadable/out-of-range trajectory), that's
    # reported explicitly after the main loop rather than silently
    # substituting a different run the user didn't ask about.
    target_row = None
    if args.target is not None:
        target_row = min(usable, key=lambda r: abs(r["_energy_demand"] - args.target))
        print(f"--target {args.target:g} GJ/t -> closest match: run_id="
              f"{target_row.get('run_id', '?')}, energy_demand="
              f"{target_row['_energy_demand']:.6g} GJ/t "
              f"(|diff|={abs(target_row['_energy_demand'] - args.target):.6g}).")

    plotted: dict[str, list[tuple[np.ndarray, np.ndarray, dict]]] = {
        "optimum": [], "target": [], "near_optimum": [], "other": []
    }
    skipped_no_traj = 0
    skipped_out_of_range = 0
    for row in usable:
        traj = read_trajectory(Path(row["raw_data_path"]))
        if traj is None:
            skipped_no_traj += 1
            continue
        t, T = traj
        peak_idx = int(np.argmax(T))
        t_shifted_us = (t - t[peak_idx]) * 1.0e6
        mask = (t_shifted_us >= -args.time_range) & (t_shifted_us <= args.time_range)
        if not mask.any():
            skipped_out_of_range += 1
            continue

        # Checked in this priority order so each run lands in exactly one
        # category -- see COLOR CODING in the module docstring. `row is
        # target_row` (identity, not value equality) is what lets the
        # target run be found even if several rows share the same
        # energy_demand.
        if row["_energy_demand"] == min_energy:
            category = "optimum"
        elif target_row is not None and row is target_row:
            category = "target"
        elif row["_energy_demand"] <= threshold:
            category = "near_optimum"
        else:
            category = "other"
        plotted[category].append((t_shifted_us[mask], T[mask], row))

    if skipped_no_traj:
        print(f"Skipping {skipped_no_traj} optimum run(s): no readable trajectory CSV "
              f"(missing/empty raw_data_path file).")
    if skipped_out_of_range:
        print(f"Skipping {skipped_out_of_range} optimum run(s): peak-shifted trajectory "
              f"never enters the +/-{args.time_range:g} us window.")
    if target_row is not None and not plotted["target"] and target_row["_energy_demand"] != min_energy:
        print(f"WARNING: the closest match to --target {args.target:g} GJ/t (run_id="
              f"{target_row.get('run_id', '?')}) has no plottable trajectory -- "
              f"no blue curve drawn.")

    n_optimum = len(plotted["optimum"])
    n_target = len(plotted["target"])
    n_near = len(plotted["near_optimum"])
    n_other = len(plotted["other"])
    print(f"Plotting {n_optimum} optimum + {n_near} within {args.within_pct:g}% + {n_other} other"
          + (f" + {n_target} target" if target_row is not None else "")
          + f" curve(s) (min energy_demand = {min_energy:.6g} GJ/t).")

    fig, ax = plt.subplots(figsize=(10, 6))
    for category in ("other", "near_optimum", "target", "optimum"):  # draw order: highlights on top
        for t_us, T_vals, _row in plotted[category]:
            ax.plot(t_us, T_vals, **STYLE[category])

    legend_handles = [
        Line2D([0], [0], color=STYLE["optimum"]["color"], linewidth=STYLE["optimum"]["linewidth"],
               label=f"optimum (energy_demand={min_energy:.4g} GJ/t, n={n_optimum})"),
        Line2D([0], [0], color=STYLE["near_optimum"]["color"], linewidth=1.5,
               label=f"within {args.within_pct:g}% of optimum (n={n_near})"),
        Line2D([0], [0], color=STYLE["other"]["color"], linewidth=1.5,
               label=f"other (n={n_other})"),
    ]
    if n_target:
        legend_handles.append(
            Line2D([0], [0], color=STYLE["target"]["color"], linewidth=STYLE["target"]["linewidth"],
                   label=f"closest to target {args.target:g} GJ/t "
                         f"(energy_demand={target_row['_energy_demand']:.4g} GJ/t, n={n_target})")
        )
    ax.legend(handles=legend_handles, loc="best")

    ax.set_xlim(-args.time_range, args.time_range)
    ax.set_xlabel("t - t_peak [us]")
    ax.set_ylabel("T [K]")
    ax.set_title(f"Temperature vs. peak-shifted time "
                 f"({n_optimum + n_near + n_other + n_target} optimum run(s))")
    fig.tight_layout()

    output_path = args.output or (OUTPUT_DIR / "optimum_temperatures.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
