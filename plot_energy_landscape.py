"""
plot_energy_landscape.py

For each distinct nu_L value among archived run_type=="optimum" entries in
raw_data/master_index.csv, draws a heatmap of energy_demand [GJ/t] over
the (P_amb, f) grid actually explored for that nu_L. Each panel uses its own
independent log-scale color axis capped at panel_optimum * 1.1 (if --vmax-opt-plus is used,
or custom --vmax), with its own individual colorbar.
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

GROUP_NDIGITS = 6


def load_grid_rows(master_index_path: Path) -> list[dict]:
    if not master_index_path.exists():
        raise SystemExit(f"{master_index_path} does not exist -- run archive.py first.")
    with open(master_index_path, newline="") as f:
        rows = [row for row in csv.DictReader(f) if row.get("run_type", "").startswith("optimum")]
    if not rows:
        with open(master_index_path, newline="") as f:
            rows = [row for row in csv.DictReader(f) if "optimum" in row.get("run_type", "")]
    if not rows:
        raise SystemExit(f"No optimum-related rows found in {master_index_path}.")

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
        print(f"Skipping {skipped} optimum row(s): missing/non-numeric nu_L, f, P_amb, or energy_demand.")
    if not usable:
        raise SystemExit("No optimum row has usable nu_L/f/P_amb/energy_demand values.")
    return usable


def build_grids(rows: list[dict]) -> dict[float, dict[tuple[float, float], float]]:
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
        print(f"Note: {n_duplicates} (nu_L, P_amb, f) triple(s) appeared more than once -- using lowest energy_demand.")
    return grids


def format_tick(value: float) -> str:
    return f"{value:g}"


def main():
    parser = argparse.ArgumentParser(
        description="Heatmap energy_demand over (P_amb, f) for every archived nu_L, panel-specific color scales."
    )
    parser.add_argument("--master-index", type=Path, default=MASTER_INDEX_PATH,
                         help=f"path to master_index.csv (default: {MASTER_INDEX_PATH})")
    parser.add_argument("--max-cols", type=int, default=5,
                         help="panels per row before wrapping into another row (default: 5)")
    parser.add_argument("--title", default="Energy Demand Parameter Landscape across Liquid Viscosities",
                         help="figure suptitle")
    parser.add_argument("--vmax-opt-plus", action="store_true",
                         help="Automatically set each panel's vmax to its local optimum + 10% (1.1 * local_min).")
    parser.add_argument("--vmax", type=float, default=None,
                         help="Override panel vmax globally with this value if specified.")
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

    n_panels = len(nu_L_values)
    ncols = min(n_panels, args.max_cols)
    nrows = math.ceil(n_panels / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 5.2 * nrows + 0.6),
                             squeeze=False)

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")

    for idx, nu_L in enumerate(nu_L_values):
        ax = axes[idx // ncols][idx % ncols]
        grid = grids[nu_L]

        P_amb_values = sorted({key[0] for key in grid})
        f_values = sorted({key[1] for key in grid})
        P_index = {v: i for i, v in enumerate(P_amb_values)}
        f_index = {v: i for i, v in enumerate(f_values)}

        data = np.full((len(f_values), len(P_amb_values)), np.nan)
        panel_energies = []
        for (P_amb, f), energy_demand in grid.items():
            data[f_index[f], P_index[P_amb]] = energy_demand
            panel_energies.append(energy_demand)

        panel_min = min(panel_energies)
        panel_max = max(panel_energies)

        if args.vmax is not None:
            vmax_val = args.vmax
        elif args.vmax_opt_plus:
            vmax_val = panel_min * 1.1
        else:
            vmax_val = panel_max

        if vmax_val <= panel_min:
            vmax_val = panel_max if panel_max > panel_min else panel_min * 1.01

        norm = LogNorm(vmin=panel_min, vmax=vmax_val)
        n_clipped = sum(1 for ed in panel_energies if ed > vmax_val)

        im = ax.imshow(data, origin="lower", aspect="auto", cmap=cmap, norm=norm)

        ax.set_xticks(range(len(P_amb_values)))
        ax.set_xticklabels([format_tick(v) for v in P_amb_values], rotation=90 if len(P_amb_values) > 8 else 0)
        ax.set_yticks(range(len(f_values)))
        ax.set_yticklabels([format_tick(v) for v in f_values])
        ax.set_xlabel("Ambient Pressure $P_{amb}$ [bar]")
        if idx % ncols == 0:
            ax.set_ylabel("Frequency $f$ [kHz]")
        ax.set_title(rf"$nu_L$ = {format_tick(nu_L)} cSt")

        min_cell = min(grid, key=grid.get)
        min_energy = grid[min_cell]
        is_global_best = (min_energy == global_min)
        star_color = "red" if is_global_best else "orange"
        ax.plot(P_index[min_cell[0]], f_index[min_cell[1]], marker="*", markersize=18,
                markerfacecolor=star_color, markeredgecolor="black", markeredgewidth=0.8,
                linestyle="none", zorder=5)
        ax.text(0.03, 0.97, f"Min: {min_energy:.4g} GJ/t", transform=ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85,
                          edgecolor=star_color, linewidth=1.0))

        cbar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.046, pad=0.15,
                            extend="max" if n_clipped else "neither")
        cbar.ax.tick_params(labelsize=7)
        cbar.set_label("Energy Demand [GJ/t]", fontsize=7)

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

    output_path = args.output or (OUTPUT_DIR / "energy_landscape.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved panel-specific colorbar energy landscape to: {output_path}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
