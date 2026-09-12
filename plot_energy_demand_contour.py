"""
plot_energy_demand_contour.py

Given a run_id (as cataloged by finalize_runs.py in /raw_data/master_index.csv),
loads raw_data/run_<id>.csv + run_<id>.json and draws a filled contour plot of
energy_demand [GJ/t] over the R_E (log, µm) x p_A (linear, bar) plane, with the
best (lowest energy_demand) point marked and annotated.

Usage:
    python plot_energy_demand_contour.py <run_id>
    python plot_energy_demand_contour.py 0 --vmax 10000 --raw-data-dir raw_data --plots-dir plots

Assumes the run being plotted actually swept R_E and p_A (typical for this
project's optimization sweeps) — other axis combinations aren't handled.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_run(run_id: int, raw_data_dir: Path) -> tuple[pd.DataFrame, dict]:
    csv_path = raw_data_dir / f"run_{run_id}.csv"
    json_path = raw_data_dir / f"run_{run_id}.json"
    if not csv_path.exists():
        raise FileNotFoundError(f"{csv_path} not found")
    if not json_path.exists():
        raise FileNotFoundError(f"{json_path} not found")

    df = pd.read_csv(csv_path)
    with open(json_path) as f:
        settings = json.load(f)
    return df, settings


def extract_plot_data(df: pd.DataFrame, vmax: float) -> pd.DataFrame:
    """Build a tidy dataframe with R_E [um], p_A [bar], and a plot-ready
    energy_demand value in GJ/t. Failed simulations are set to NaN so they
    render as blank/white rather than being folded into the color scale."""
    out = pd.DataFrame()
    out["R_E_um"] = df["R_E"].astype(float) * 1.0e6

    # excitation_params is a semicolon-separated string "p_A;freq;" per row
    p_A_pa = df["excitation_params"].str.split(";").str[0].astype(float)
    out["p_A_bar"] = p_A_pa.abs() / 1.0e5

    success = df["success"].astype(str).str.strip().str.lower() == "true"
    energy_demand = pd.to_numeric(df["energy_demand"], errors="coerce")

    plot_value = energy_demand.where(success)  # NaN for failed sims
    out["energy_demand_plot"] = plot_value.clip(upper=vmax)

    # unclipped, successful-only value, used to find the true optimum
    out["energy_demand_true"] = plot_value
    return out


def plot_contour(run_id: int, raw_data_dir: str = "raw_data", plots_dir: str = "plots",
                  vmax: float = 10000.0) -> Path:
    raw_data_dir = Path(raw_data_dir)
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(exist_ok=True)

    df, settings = load_run(run_id, raw_data_dir)
    data = extract_plot_data(df, vmax)

    pivot = data.pivot_table(index="p_A_bar", columns="R_E_um", values="energy_demand_plot")
    R_E_vals = pivot.columns.values
    p_A_vals = pivot.index.values
    Z = np.ma.masked_invalid(pivot.values)  # NaN (failed sims) -> blank/white, not colored

    best_idx = data["energy_demand_true"].idxmin()
    best_R_E = data.loc[best_idx, "R_E_um"]
    best_p_A = data.loc[best_idx, "p_A_bar"]
    best_value = data.loc[best_idx, "energy_demand_true"]

    fig, ax = plt.subplots(figsize=(9, 8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    cmap = plt.get_cmap("jet").copy()
    cmap.set_bad("white")

    levels = np.linspace(0, vmax, 101)
    contour = ax.contourf(R_E_vals, p_A_vals, Z, levels=levels, cmap=cmap, extend="max")
    cbar = fig.colorbar(contour, ax=ax, extend="max")
    cbar.set_label("Energy demand (GJ/t)", fontsize=12, fontweight="bold")

    ax.plot(best_R_E, best_p_A, marker="o", color="lime", markeredgecolor="black",
             markersize=12, markeredgewidth=1.5, linestyle="none")
    ax.text(
        0.03, 0.95,
        f"Optimal solution:\n{best_value:.2f}  GJ/t",
        transform=ax.transAxes, fontsize=13, fontweight="bold",
        va="top", ha="left",
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="black"),
    )

    ax.set_xscale("log")
    ax.set_xlabel(r"R$_E$ ($\mu$m)", fontsize=14, fontweight="bold")
    ax.set_ylabel(r"p$_A$ (bar)", fontsize=14, fontweight="bold")

    out_path = plots_dir / f"run_{run_id}_energy_demand_contour.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Plot energy_demand contour (R_E vs p_A) for a finalized run.")
    parser.add_argument("run_id", type=int)
    parser.add_argument("--vmax", type=float, default=10000.0)
    parser.add_argument("--raw-data-dir", default="raw_data")
    parser.add_argument("--plots-dir", default="plots")
    args = parser.parse_args()

    out_path = plot_contour(args.run_id, args.raw_data_dir, args.plots_dir, args.vmax)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
