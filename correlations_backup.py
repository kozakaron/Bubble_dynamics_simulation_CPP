"""
correlations.py

Loads unique simulation parameters and produced ammonia from master_index.csv,
computes surface tension (sigma) from nu_L, estimates theoretical maximum ammonia
based on bubble volume and internal pressure (P_amb + 2*sigma/R_E), computes conversion rate,
and plots all pairwise correlations in a multi-panel grid layout.
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

import solver_common

MASTER_INDEX_PATH = Path("raw_data") / "master_index.csv"
OUTPUT_DIR = Path("processed_data")

GROUP_NDIGITS = 6
CORE_COLUMNS = ["P_amb", "f", "nu_L", "p_A", "R_E", "energy_demand", "konverziós_ráta"]


def load_and_compute_data(master_index_path: Path) -> pd.DataFrame:
    """Loads master_index.csv, computes surface tension, max ammonia, and conversion rate."""
    if not master_index_path.exists():
        raise SystemExit(f"{master_index_path} does not exist -- run archive.py or update_master_index.py first.")

    with open(master_index_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise SystemExit(f"No rows found in {master_index_path}.")

    usable_records = []
    seen_keys = set()
    skipped_duplicates = 0
    skipped_errors = 0

    for row in rows:
        try:
            p_A = float(row.get("p_A", row.get("p_A_opt", np.nan)))
            R_E = float(row.get("R_E", row.get("R_E_opt", np.nan)))
            f_khz = float(row.get("f", np.nan))
            nu_L = float(row.get("nu_L", np.nan))
            P_amb = float(row.get("P_amb", np.nan))
            energy_demand = float(row.get("energy_demand", np.nan))
            termelt_ammónia = float(row.get("termelt_ammónia", 0.0))
        except (ValueError, TypeError):
            skipped_errors += 1
            continue

        if any(math.isnan(v) for v in [p_A, R_E, f_khz, nu_L, P_amb, energy_demand]):
            skipped_errors += 1
            continue

        # Duplikációk kiszűrése az 5 bemeneti paraméter alapján
        key = (
            round(P_amb, GROUP_NDIGITS),
            round(f_khz, GROUP_NDIGITS),
            round(nu_L, GROUP_NDIGITS),
            round(p_A, GROUP_NDIGITS),
            round(R_E, GROUP_NDIGITS),
        )

        if key in seen_keys:
            skipped_duplicates += 1
            continue
        seen_keys.add(key)

        # 1. Felületi feszültség (sigma) pontos lekérése a nu_L-ből
        liq_props = solver_common.compute_liquid_properties(nu_L)
        sigma = liq_props["surfactant"]["value"]  # [N/m]

        # 2. Kezdeti nyomás a buborékban: P_amb + 2 * sigma / R_E
        # Mértékegységek: P_amb [bar] -> Pa (1 bar = 1e5 Pa), R_E [um] -> m (1e-6 m)
        P_amb_Pa = P_amb * 1.0e5
        R_E_m = R_E * 1.0e-6
        internal_pressure = P_amb_Pa + (2.0 * sigma / R_E_m) if R_E_m > 0 else P_amb_Pa

        # 3. Elméleti maximális ammónia számítás (3:1 H2:N2 arány, ideális gáz törvény)
        bubble_volume = (4.0 / 3.0) * math.pi * (R_E_m ** 3)
        R_gas = 8.314
        T_ref = 298.15  # Referencia hőmérséklet [K]
        
        total_moles_potential = (internal_pressure * bubble_volume) / (R_gas * T_ref)
        max_possible_ammonia = total_moles_potential * 0.75  # A keverék 3/4 része hidrogén

        # 4. Konverziós ráta számítása
        conversion_rate = (termelt_ammónia / max_possible_ammonia) if max_possible_ammonia > 0 else 0.0

        usable_records.append({
            "P_amb": P_amb,
            "f": f_khz,
            "nu_L": nu_L,
            "p_A": p_A,
            "R_E": R_E,
            "energy_demand": energy_demand,
            "konverziós_ráta": conversion_rate,  # Itt már benne van a dictionary-ben!
        })

    if skipped_duplicates > 0:
        print(f"Skipped {skipped_duplicates} duplicate parameter combination(s).")
    if skipped_errors > 0:
        print(f"Skipped {skipped_errors} row(s) due to missing or invalid numeric values.")

    df = pd.DataFrame(usable_records)
    print(f"Successfully loaded and computed {len(df)} unique data point(s).")
    return df


def plot_all_pairwise(df: pd.DataFrame, max_cols: int, output_path: Path, hue_col: str | None = None):
    """Plots all pairwise combinations of CORE_COLUMNS in a subplot grid."""
    pairs = list(itertools.combinations(CORE_COLUMNS, 2))
    n_plots = len(pairs)

    ncols = min(n_plots, max_cols)
    nrows = math.ceil(n_plots / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4.8 * ncols, 4.8 * nrows + 0.6),
        squeeze=False
    )

    for idx, (col_x, col_y) in enumerate(pairs):
        row_i = idx // ncols
        col_j = idx % ncols
        ax = axes[row_i][col_j]

        x_vals = df[col_x]
        y_vals = df[col_y]

        if hue_col and hue_col in df.columns:
            sc = ax.scatter(
                x_vals, y_vals,
                c=df[hue_col], cmap="viridis",
                edgecolor="black", linewidth=0.5, alpha=0.85, s=40
            )
            cbar = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.046, pad=0.18)
            cbar.ax.tick_params(labelsize=8)
            cbar.set_label(hue_col, fontsize=8)
        else:
            ax.scatter(
                x_vals, y_vals,
                color="navy", edgecolor="black", linewidth=0.5, alpha=0.75, s=35
            )

        ax.set_xlabel(col_x, fontsize=10, fontweight="bold")
        ax.set_ylabel(col_y, fontsize=10, fontweight="bold")
        ax.set_title(f"{col_y} vs {col_x}", fontsize=11)

        # Logaritmikus skála alkalmazása energiaigénynél
        if "energy_demand" in col_y:
            ax.set_yscale("log")
        if "energy_demand" in col_x:
            ax.set_xscale("log")

        # Tiszta rácsozás a leolvasás megkönnyítésére
        ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.7)
        ax.tick_params(axis="both", which="major", labelsize=9)

    # Üres panelek eltüntetése a rács végéről
    for idx in range(n_plots, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Pairwise Parameter Correlations Including Conversion Rate", fontsize=15, fontweight="bold", y=0.995)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved correlation panel plot to: {output_path}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Generate pairwise correlation grid including conversion rate.")
    parser.add_argument("--master-index", type=Path, default=MASTER_INDEX_PATH)
    parser.add_argument("--max-cols", type=int, default=5, help="Maximum panels per row (default: 5)")
    parser.add_argument("--hue", type=str, default="nu_L", help="Optional column to color points by (default: nu_L)")
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR / "pairwise_correlations_with_conversion.png")
    args = parser.parse_args()

    df = load_and_compute_data(args.master_index)
    if df.empty:
        raise SystemExit("No usable data points found.")

    plot_all_pairwise(df, max_cols=args.max_cols, output_path=args.output, hue_col=args.hue)


if __name__ == "__main__":
    main()