"""
plot_energy_landscape.py

For each distinct nu_L value among archived run_type=="optimum" entries in
raw_data/master_index.csv, draws a heatmap of energy_demand [GJ/t] over
the (P_amb, f) grid actually explored for that nu_L -- a Cartesian batch
(batch_global_optimise.py --combo-mode cartesian, the default) produces
exactly one "optimum" row per (nu_L, f, P_amb) triple, i.e. exactly the
regular grid a heatmap needs, with energy_demand already representing the
OPTIMAL (R_E, p_A) result at that (nu_L, f, P_amb) point.

All panels share ONE log-scale color axis (computed across every plotted
point, not per panel), so energy_demand is directly comparable across
nu_L values -- matching the reference figure. Pass --vmax to cap the top
of that color scale yourself (e.g. to stop one outlier panel from
washing out the rest); anything above --vmax still renders, clipped to
the top color, and the colorbar grows a small triangular arrow to say so.

Each panel marks its own local-minimum cell with a star: RED if that
cell is also the GLOBAL minimum across every nu_L (there's exactly one
red star in the whole figure), ORANGE for every other panel's own best
cell. A small "Min: ... GJ/t" text box sits next to each star, its border
colored to match.

NOT included, per instruction: any boundary_warning / "Boundary Hit"
marking -- every plotted point here is just its own (P_amb, f,
energy_demand) triple.

GRID LAYOUT: panels are placed --max-cols per row (default 5, so a
typical <=5-value nu_L sweep -- like the reference figure -- renders as
one row); once there are more nu_L values than that, panels wrap into
further rows, giving an actual matrix layout for sweeps over many nu_L
values.

AXES: P_amb and f are plotted on an evenly-spaced INDEX grid (one cell
per distinct value actually present for that nu_L, in sorted order), with
tick labels showing the real values -- not a literal linear/log data
axis. This is deliberate: a Cartesian sweep's own P_amb/f values are
rarely evenly (or evenly-log) spaced (see e.g. combos_template_cartesian.csv),
so an index grid is what keeps every explored value visible as its own
full-size cell, exactly like the reference figure's even-looking cells
next to unevenly-spaced tick labels.

A (nu_L, f, P_amb) triple that appears more than once (e.g. an old batch
result never cleaned out, or a rerun) uses whichever duplicate has the
LOWEST energy_demand for that cell -- duplicates are counted and reported,
not silently averaged or arbitrarily picked.

Usage:
    python plot_energy_landscape.py
    python plot_energy_landscape.py --vmax 2000
    python plot_energy_landscape.py --max-cols 4
    python plot_energy_landscape.py --title "My sweep"
    python plot_energy_landscape.py --no-show
    python plot_energy_landscape.py --master-index raw_data/master_index.csv --output processed_data/energy_landscape.png

Saves a PNG to processed_data/energy_landscape.png by default (pass
--output to change) and also opens an interactive matplotlib window
unless --no-show is given.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

MASTER_INDEX_PATH = Path("raw_data") / "master_index.csv"
OUTPUT_DIR = Path("processed_data")

# Rounding used only to GROUP rows into grid cells/panels -- absorbs
# floating-point noise (e.g. nu_L derived from mu_L/rho_L via division)
# without merging genuinely different sweep values. See the module
# docstring in batch_global_optimise.py's own resume feature for the same
# kind of float-noise-vs-real-difference distinction.
GROUP_NDIGITS = 6


def load_grid_rows(master_index_path: Path) -> list[dict]:
    """Reads master_index.csv, keeps run_type=='optimum' rows, and parses
    nu_L/f/P_amb/energy_demand to floats. Rows missing any of these (or
    with a non-numeric value -- e.g. a swept Range string, which shouldn't
    occur for an "optimum" row but is checked defensively) are skipped and
    counted rather than crashing the whole plot."""
    if not master_index_path.exists():
        raise SystemExit(f"{master_index_path} does not exist -- run archive.py first.")
    with open(master_index_path, newline="") as f:
        rows = [row for row in csv.DictReader(f) if row.get("run_type") == "optimum"]
    if not rows:
        raise SystemExit(f"No run_type=='optimum' rows found in {master_index_path}.")

    usable = []
    skipped = 0
    for row in rows:
        try:
            row["_nu_L"] = float(row["nu_L"])
            row["_f"] = float(row["f"])
            row["_P_amb"] = float(row["P_amb"])
            row["_energy_demand"] = float(row["energy_demand"])
        except (KeyError, ValueError):
            skipped += 1
            continue
        usable.append(row)
    if skipped:
        print(f"Skipping {skipped} optimum row(s): missing/non-numeric nu_L, f, P_amb, "
              f"or energy_demand.")
    if not usable:
        raise SystemExit("No optimum row has usable nu_L/f/P_amb/energy_demand values.")
    return usable


def build_grids(rows: list[dict]) -> dict[float, dict[tuple[float, float], float]]:
    """Groups rows by (rounded) nu_L, and within each nu_L by (rounded)
    (P_amb, f) -> energy_demand, keeping the lowest energy_demand on a
    duplicate cell. Returns {nu_L: {(P_amb, f): energy_demand}}."""
    grids: dict[float, dict[tuple[float, float], float]] = {}
    n_duplicates = 0
    for row in rows:
        nu_L = round(row["_nu_L"], GROUP_NDIGITS)
        P_amb = round(row["_P_amb"], GROUP_NDIGITS)
        f = round(row["_f"], GROUP_NDIGITS)
        cell_key = (P_amb, f)
        grid = grids.setdefault(nu_L, {})
        if cell_key in grid:
            n_duplicates += 1
            grid[cell_key] = min(grid[cell_key], row["_energy_demand"])
        else:
            grid[cell_key] = row["_energy_demand"]
    if n_duplicates:
        print(f"Note: {n_duplicates} (nu_L, P_amb, f) triple(s) appeared more than once -- "
              f"using the lowest energy_demand found for each.")
    return grids


def format_tick(value: float) -> str:
    """Compact numeric tick label, e.g. 0.65, 1, 10, 100 -- not 1.0, 10.0."""
    return f"{value:g}"


def main():
    parser = argparse.ArgumentParser(
        description="Heatmap energy_demand over (P_amb, f) for every archived nu_L, one "
                    "panel per nu_L, shared log color scale."
    )
    parser.add_argument("--master-index", type=Path, default=MASTER_INDEX_PATH,
                         help=f"path to master_index.csv (default: {MASTER_INDEX_PATH})")
    parser.add_argument("--max-cols", type=int, default=5,
                         help="panels per row before wrapping into another row (default: 5)")
    parser.add_argument("--title", default="Energy Demand Parameter Landscape across Liquid Viscosities",
                         help="figure suptitle")
    parser.add_argument("--vmax", type=float, default=None,
                         help="cap the color scale's upper bound at this energy_demand "
                              "[GJ/t] (default: the actual max across all plotted points). "
                              "Points above it still plot, clipped to the top color.")
    parser.add_argument("--output", type=Path, default=None,
                         help=f"PNG output path (default: {OUTPUT_DIR / 'energy_landscape.png'})")
    parser.add_argument("--no-show", action="store_true",
                         help="skip the interactive window, just save the PNG")
    args = parser.parse_args()

    rows = load_grid_rows(args.master_index)
    grids = build_grids(rows)
    nu_L_values = sorted(grids.keys())

    all_energy = [ed for grid in grids.values() for ed in grid.values()]
    global_min = min(all_energy)
    data_max = max(all_energy)
    vmax = data_max
    n_clipped = 0
    if args.vmax is not None:
        if args.vmax <= global_min:
            raise SystemExit(f"--vmax ({args.vmax:g}) must be greater than the smallest "
                              f"energy_demand actually present ({global_min:.6g}).")
        vmax = args.vmax
        n_clipped = sum(1 for ed in all_energy if ed > vmax)
        if n_clipped:
            print(f"--vmax {vmax:g}: {n_clipped} of {len(all_energy)} point(s) exceed it and "
                  f"are shown clipped to the top color.")
    norm = LogNorm(vmin=global_min, vmax=vmax)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")  # cells with no data (gaps in the sweep) render blank, not black

    n_panels = len(nu_L_values)
    ncols = min(n_panels, args.max_cols)
    nrows = math.ceil(n_panels / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 4.6 * nrows + 0.6),
                              squeeze=False)

    im = None
    for idx, nu_L in enumerate(nu_L_values):
        ax = axes[idx // ncols][idx % ncols]
        grid = grids[nu_L]

        P_amb_values = sorted({key[0] for key in grid})
        f_values = sorted({key[1] for key in grid})
        P_index = {v: i for i, v in enumerate(P_amb_values)}
        f_index = {v: i for i, v in enumerate(f_values)}

        data = np.full((len(f_values), len(P_amb_values)), np.nan)
        for (P_amb, f), energy_demand in grid.items():
            data[f_index[f], P_index[P_amb]] = energy_demand

        im = ax.imshow(data, origin="lower", aspect="auto", cmap=cmap, norm=norm)

        ax.set_xticks(range(len(P_amb_values)))
        ax.set_xticklabels([format_tick(v) for v in P_amb_values], rotation=90 if len(P_amb_values) > 8 else 0)
        ax.set_yticks(range(len(f_values)))
        ax.set_yticklabels([format_tick(v) for v in f_values])
        ax.set_xlabel("Ambient Pressure $P_{amb}$ [bar]")
        if idx % ncols == 0:
            ax.set_ylabel("Frequency $f$ [kHz]")
        ax.set_title(rf"$\nu_L$ = {format_tick(nu_L)} cSt")

        min_cell = min(grid, key=grid.get)
        min_energy = grid[min_cell]
        # Exactly one panel's local min is ALSO the global min (computed
        # from the same underlying float values, so a direct == is safe
        # here, same reasoning as the other plot scripts' optimum checks)
        # -- that one gets the red star; every other panel's own best
        # cell gets orange instead.
        is_global_best = (min_energy == global_min)
        star_color = "red" if is_global_best else "orange"
        ax.plot(P_index[min_cell[0]], f_index[min_cell[1]], marker="*", markersize=18,
                markerfacecolor=star_color, markeredgecolor="black", markeredgewidth=0.8,
                linestyle="none", zorder=5)
        ax.text(0.03, 0.97, f"Min: {min_energy:.4g} GJ/t", transform=ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85,
                          edgecolor=star_color, linewidth=1.0))

    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle(args.title, fontsize=14, fontweight="bold")
    fig.legend(
        handles=[
            Line2D([0], [0], marker="*", color="none", markerfacecolor="red",
                   markeredgecolor="black", markersize=12, label="global optimum"),
            Line2D([0], [0], marker="*", color="none", markerfacecolor="orange",
                   markeredgecolor="black", markersize=12, label="best in this nu_L slice"),
        ],
        loc="upper right", bbox_to_anchor=(0.99, 0.99), fontsize=9, framealpha=0.9,
    )
    fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.04, pad=0.08,
                 extend="max" if n_clipped else "neither",
                 label="Optimal Energy Demand [GJ/t] (log scale)")

    print(f"Plotted {n_panels} nu_L panel(s), {len(rows)} grid point(s) total, "
          f"energy_demand range {global_min:.4g} to {data_max:.4g} GJ/t "
          f"(color scale capped at {vmax:.4g}).")

    output_path = args.output or (OUTPUT_DIR / "energy_landscape.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
