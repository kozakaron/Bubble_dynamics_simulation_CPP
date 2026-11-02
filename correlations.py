"""
correlations.py

Generates line plots for all pairwise optimal projections among the 5 core parameters
(P_amb, f, nu_L, p_A, R_E). 

Rules:
- P_amb is ALWAYS on the horizontal (X) axis if present.
- R_E is ALWAYS on the vertical (Y) axis if present.
- Minimized energy_demand over the remaining 3 parameters.
- Plotted as single continuous line plots (plt.plot) arranged in a grid (max 5 per row).
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MASTER_INDEX_PATH = Path("raw_data") / "master_index.csv"
OUTPUT_DIR = Path("processed_data")

GROUP_NDIGITS = 6
CORE_COLUMNS = ["P_amb", "f", "nu_L", "p_A", "R_E"]


def load_raw_data(master_index_path: Path) -> pd.DataFrame:
    """Loads all valid numeric rows from master_index.csv."""
    if not master_index_path.exists():
        raise SystemExit(f"{master_index_path} does not exist -- run archive.py first.")

    with open(master_index_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise SystemExit(f"No rows found in {master_index_path}.")

    records = []
    for row in rows:
        try:
            p_A = float(row.get("p_A", row.get("p_A_opt", np.nan)))
            R_E = float(row.get("R_E", row.get("R_E_opt", np.nan)))
            f_khz = float(row.get("f", np.nan))
            nu_L = float(row.get("nu_L", np.nan))
            P_amb = float(row.get("P_amb", np.nan))
            energy_demand = float(row.get("energy_demand", np.nan))
        except (ValueError, TypeError):
            continue

        if any(math.isnan(v) for v in [p_A, R_E, f_khz, nu_L, P_amb, energy_demand]):
            continue

        records.append({
            "P_amb": P_amb,
            "f": f_khz,
            "nu_L": nu_L,
            "p_A": p_A,
            "R_E": R_E,
            "energy_demand": energy_demand,
        })

    df = pd.DataFrame(records)
    if df.empty:
        raise SystemExit("No valid numeric data found in master index.")
    return df


def get_ordered_pair(a: str, b: str) -> tuple[str, str]:
    """
    Determines X and Y axes based on rules:
    - P_amb must be X if present.
    - R_E must be Y if present.
    """
    s = {a, b}
    if "P_amb" in s:
        x = "P_amb"
        y = list(s - {"P_amb"})[0]
    elif "R_E" in s:
        y = "R_E"
        x = list(s - {"R_E"})[0]
    else:
        lst = [p for p in CORE_COLUMNS if p in s]
        x, y = lst[0], lst[1]
    return x, y


def plot_all_optimal_line_pairs(df: pd.DataFrame, max_cols: int, output_path: Path):
    """Generates a grid of single-line optimal projection plots."""
    raw_pairs = list(itertools.combinations(CORE_COLUMNS, 2))
    unique_pairs = []
    seen = set()
    for a, b in raw_pairs:
        x, y = get_ordered_pair(a, b)
        if (x, y) not in seen:
            seen.add((x, y))
            unique_pairs.append((x, y))

    n_plots = len(unique_pairs)
    ncols = min(n_plots, max_cols)
    nrows = math.ceil(n_plots / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.8 * ncols, 4.2 * nrows + 0.6),
        squeeze=False
    )

    for idx, (col_x, col_y) in enumerate(unique_pairs):
        row_i = idx // ncols
        col_j = idx % ncols
        ax = axes[row_i][col_j]

        # Csoportosítás és kerekítés a zaj szűrésére
        temp_df = df.copy()
        temp_df["x_round"] = temp_df[col_x].round(GROUP_NDIGITS)

        # Optimalizálás: az adott X értékhez tartozó minimális energiaigény kiválasztása
        opt_rows = []
        for _, group in temp_df.groupby("x_round"):
            best_row = group.loc[group["energy_demand"].idxmin()]
            opt_rows.append(best_row)

        opt_df = pd.DataFrame(opt_rows)

        if opt_df.empty:
            ax.axis("off")
            continue

        # Rendezés az X tengely szerint a folytonos vonalhoz
        opt_df = opt_df.sort_values(by=col_x)

        # Vonalkirajzolás pontjelölőkkel
        ax.plot(
            opt_df[col_x], opt_df[col_y],
            color="navy", linestyle="-", linewidth=2.0,
            marker="o", markersize=4, markerfacecolor="red", markeredgecolor="darkred"
        )

        ax.set_xlabel(col_x, fontsize=10, fontweight="bold")
        ax.set_ylabel(col_y, fontsize=10, fontweight="bold")
        ax.set_title(f"{col_y} vs {col_x}", fontsize=11, fontweight="bold")

        # Logaritmikus skála ha R_E szerepel
        if col_x == "R_E":
            ax.set_xscale("log")
        if col_y == "R_E":
            ax.set_yscale("log")

        # Rácsozás a leolvasáshoz
        ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.7)
        ax.tick_params(axis="both", which="major", labelsize=9)

    # Üres panelek eltüntetése a rács végén
    for idx in range(n_plots, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Optimal Parameter Curves (Minimized Energy Demand)", fontsize=14, fontweight="bold", y=0.995)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"All optimal line plots saved to: {output_path}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Generate line plot grid of optimal parameter correlations.")
    parser.add_argument("--master-index", type=Path, default=MASTER_INDEX_PATH)
    parser.add_argument("--max-cols", type=int, default=5, help="Maximum panels per row (default: 5)")
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR / "optimal_line_correlations.png")
    args = parser.parse_args()

    df = load_raw_data(args.master_index)
    plot_all_optimal_line_pairs(df, max_cols=args.max_cols, output_path=args.output)


if __name__ == "__main__":
    main()