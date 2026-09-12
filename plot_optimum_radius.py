"""
plot_optimum_radius.py

Overlay plot of the R(t) bubble-radius trajectory for every archived
run_type=="optimum" entry in raw_data/master_index.csv, over the RAW
(unshifted) simulation time window t=0 to a user-given --tmax. Unlike
plot_optimum_temperatures.py, curves here are NOT shifted to align any
peak -- t=0 is the simulation's own start, exactly as saved.

UNITS: raw_data/run_<id>.csv's "t" and "R" columns are exactly what the
C++ solver wrote after ControlParameters::dimensionalize() -- verified
directly against control_parameters.cpp's dimensionalize(): `t *= t_ref`
with t_ref = 1e-9 [s] (control_parameters.h), `x[0] *= R_ref` for the
radius state variable, no other conversion applied downstream
(single_run.py/archive.py never rescale these). So on disk: t is in
SECONDS, R is in METERS. By default this script converts t to
MICROSECONDS (matching --tmax's unit) and R to MICROMETERS (matching how
R_E is expressed everywhere else in this project).

--dimless: plots NONDIMENSIONAL axes instead -- R/R_E on the vertical
axis, tau = t/t_ref on the horizontal one, where t_ref = 1/f (f = the
driving frequency in Hz, so t_ref is one driving PERIOD). EACH CURVE uses
its OWN R_E and its OWN f -- read straight from that run's own "R_E" [um]
and "f" [kHz] columns in master_index.csv (already in the units
build_param_columns() in archive.py writes them in for run_type=="optimum"
rows, which are always plain scalars there, never a swept Range) -- NOT
one shared reference across all curves, since different combinations have
different R_E/f. The plotted window is then fixed at tau=0 to tau=5 (five
driving periods) regardless of --tmax, which is ignored in this mode.

COLOR CODING: identical scheme/priority to plot_optimum_temperatures.py,
by energy_demand read straight from master_index.csv (not recomputed):
    - the single run with the global MINIMUM energy_demand -- "optimum",
      red, thick
    - IF --target is given: the single run whose OWN energy_demand is
      numerically closest to that user-specified value -- "target", blue,
      thick (skipped if --target isn't passed; if the closest match is
      also the optimum run, "optimum" wins -- same curve, drawn red)
    - every remaining run within --within-pct percent of the minimum
      (default 10%) -- "near_optimum", orange, thin
    - everything else -- "other", grey, thin, low alpha
Draw order is other -> near_optimum -> target -> optimum (each later
category on top). Legend has one entry per category actually plotted.

Usage:
    python plot_optimum_radius.py --tmax 10
    python plot_optimum_radius.py --tmax 10 --within-pct 5
    python plot_optimum_radius.py --tmax 10 --target 1500
    python plot_optimum_radius.py --dimless
    python plot_optimum_radius.py --dimless --target 1500
    python plot_optimum_radius.py --tmax 10 --no-show
    python plot_optimum_radius.py --tmax 10 --master-index raw_data/master_index.csv --output processed_data/optimum_radius.png

Saves a PNG to processed_data/optimum_radius.png (or
processed_data/optimum_radius_dimless.png in --dimless mode) by default
(pass --output to change) and also opens an interactive matplotlib window
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
    """Returns (t_seconds, R_meters) arrays read from one run's raw CSV, or
    None if the file is missing, empty (e.g. a run whose confirmation run
    hit a fatal preprocessing failure -- archive_single_result() writes an
    empty file in that case, see archive.py), or doesn't carry a usable
    "t"/"R" pair -- all treated the same way: this run can't be plotted,
    the caller counts and reports how many were skipped."""
    if not raw_data_path.exists():
        return None
    try:
        with open(raw_data_path, newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or "t" not in reader.fieldnames or "R" not in reader.fieldnames:
                return None
            t_vals: list[float] = []
            R_vals: list[float] = []
            for r in reader:
                t_vals.append(float(r["t"]))
                R_vals.append(float(r["R"]))
    except (OSError, ValueError):
        return None
    if len(t_vals) < 2:
        return None
    return np.array(t_vals), np.array(R_vals)


def main():
    parser = argparse.ArgumentParser(
        description="Overlay-plot every archived optimum run's R(t) curve from t=0 to "
                    "--tmax, color-coded by energy_demand."
    )
    parser.add_argument("--master-index", type=Path, default=MASTER_INDEX_PATH,
                         help=f"path to master_index.csv (default: {MASTER_INDEX_PATH})")
    parser.add_argument("--tmax", type=float, default=None,
                         help="microseconds -- plotted window is [0, tmax], RAW simulation "
                              "time (no peak-alignment shift, unlike the temperature plot). "
                              "Required unless --dimless is given (--dimless ignores it).")
    parser.add_argument("--dimless", action="store_true",
                         help="plot R/R_E vs. tau=t/t_ref instead (t_ref = 1/f, EACH curve's "
                              "own R_E and f from master_index.csv), window fixed to tau in "
                              "[0, 5] -- see the module docstring. Overrides --tmax.")
    parser.add_argument("--within-pct", type=float, default=10.0,
                         help="a run within this percent of the global minimum "
                              "energy_demand is drawn orange (default: 10)")
    parser.add_argument("--target", type=float, default=None,
                         help="a target energy_demand [GJ/t] -- the run whose OWN "
                              "energy_demand is numerically closest to this value is drawn "
                              "as an extra thick blue curve. Omit to skip this entirely.")
    parser.add_argument("--output", type=Path, default=None,
                         help=f"PNG output path (default: {OUTPUT_DIR / 'optimum_radius.png'})")
    parser.add_argument("--no-show", action="store_true",
                         help="skip the interactive window, just save the PNG")
    args = parser.parse_args()

    if args.dimless:
        if args.tmax is not None:
            print("Note: --tmax is ignored in --dimless mode (window is fixed to tau in [0, 5]).")
    elif args.tmax is None:
        parser.error("--tmax is required unless --dimless is given")

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
    skipped_bad_dimless_params = 0
    for row in usable:
        traj = read_trajectory(Path(row["raw_data_path"]))
        if traj is None:
            skipped_no_traj += 1
            continue
        t, R = traj

        if args.dimless:
            # EACH curve's OWN R_E [um] / f [kHz], straight from this row's
            # master_index.csv columns -- NOT a shared reference across
            # curves, since different (nu_L, f, P_amb) combinations have
            # different R_E/f. See the module docstring for why these
            # columns (rather than reading each run's own JSON) are safe
            # to trust here.
            try:
                R_E_um = float(row["R_E"])
                f_khz = float(row["f"])
                if f_khz == 0.0:
                    raise ValueError("f is zero")
            except (KeyError, ValueError):
                skipped_bad_dimless_params += 1
                continue
            f_Hz = f_khz * 1.0e3
            t_ref = 1.0 / f_Hz                 # [s] -- one driving period
            x_vals = t / t_ref                 # tau = t / t_ref, dimensionless
            y_vals = R / (R_E_um * 1.0e-6)      # R / R_E, dimensionless (both in meters)
            mask = (x_vals >= 0.0) & (x_vals <= 5.0)
        else:
            x_vals = t * 1.0e6                 # s -> us, NO shift -- t=0 is the simulation's own start
            y_vals = R * 1.0e6                 # m -> um
            mask = (x_vals >= 0.0) & (x_vals <= args.tmax)

        if not mask.any():
            skipped_out_of_range += 1
            continue

        # Same priority order as plot_optimum_temperatures.py -- see
        # COLOR CODING in the module docstring.
        if row["_energy_demand"] == min_energy:
            category = "optimum"
        elif target_row is not None and row is target_row:
            category = "target"
        elif row["_energy_demand"] <= threshold:
            category = "near_optimum"
        else:
            category = "other"
        plotted[category].append((x_vals[mask], y_vals[mask], row))

    if skipped_no_traj:
        print(f"Skipping {skipped_no_traj} optimum run(s): no readable trajectory CSV "
              f"(missing/empty raw_data_path file).")
    if skipped_bad_dimless_params:
        print(f"Skipping {skipped_bad_dimless_params} optimum run(s): missing/invalid "
              f"R_E or f in master_index.csv (needed for --dimless).")
    if skipped_out_of_range:
        window_desc = "tau in [0, 5]" if args.dimless else f"[0, {args.tmax:g}] us"
        print(f"Skipping {skipped_out_of_range} optimum run(s): trajectory never enters "
              f"the {window_desc} window.")
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
        for t_vals, R_vals, _row in plotted[category]:
            ax.plot(t_vals, R_vals, **STYLE[category])

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

    n_total = n_optimum + n_near + n_other + n_target
    if args.dimless:
        ax.set_xlim(0.0, 5.0)
        ax.set_xlabel("tau = t / t_ref  (t_ref = 1/f, each curve's own driving period)")
        ax.set_ylabel("R / R_E  (each curve's own R_E)")
        ax.set_title(f"Bubble radius vs. time -- nondimensional ({n_total} optimum run(s))")
    else:
        ax.set_xlim(0.0, args.tmax)
        ax.set_xlabel("t [us]")
        ax.set_ylabel("R [um]")
        ax.set_title(f"Bubble radius vs. time ({n_total} optimum run(s))")
    fig.tight_layout()

    default_output_name = "optimum_radius_dimless.png" if args.dimless else "optimum_radius.png"
    output_path = args.output or (OUTPUT_DIR / default_output_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
