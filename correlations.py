from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MASTER_INDEX_PATH = Path("raw_data") / "master_index.csv"
OUTPUT_DIR = Path("processed_data")


def load_all_run_data(master_index_path: Path) -> pd.DataFrame:
    """Loads master_index.csv and extracts postproc / ammonia data from metadata JSONs."""
    if not master_index_path.exists():
        raise SystemExit(f"{master_index_path} does not exist -- run archive.py first.")
    
    df_index = pd.read_csv(master_index_path)
    records = []

    for _, row in df_index.iterrows():
        rec = row.to_dict()
        
        # Alapértelmezett numerikus mezők konverziója
        try:
            rec["p_A_opt"] = float(row.get("p_A", 0))
            rec["R_E_opt"] = float(row.get("R_E", 0))
            rec["f"] = float(row.get("f", 0))
            rec["P_amb"] = float(row.get("P_amb", 0))
            rec["nu_L"] = float(row.get("nu_L", 0))
            rec["energiaigény"] = float(row.get("energy_demand", np.nan))
        except (ValueError, TypeError):
            continue

        # Ammónia adatok kinyerése a metadata JSON-ból (ha létezik)
        meta_path = row.get("metadata_path")
        nh3_moles = np.nan
        nh3_conversion = np.nan

        if meta_path and Path(meta_path).exists():
            try:
                with open(meta_path, "r", encoding="utf-8") as f:
                    meta = json.load(f)
                
                postproc = meta.get("postproc", {})
                species_names = meta.get("mechanism", {}).get("species_names", [])
                
                # Ha a postproc tartalmazza az NH3-at vagy a nyers CSV-ből kiolvasható az utolsó lépés
                # Kereshetjük az NH3 mólszámot a postproc-ban vagy a trajektória utolsó sorában
                raw_csv = row.get("raw_data_path")
                if raw_csv and Path(raw_csv).exists():
                    df_traj = pd.read_csv(raw_csv)
                    # Megkeressük az NH3 oszlopot (lehet 'NH3' vagy hasonló nevű)
                    nh3_col = next((c for c in df_traj.columns if c.upper() in ["NH3", "NH_3"]), None)
                    if nh3_col and not df_traj.empty:
                        nh3_moles = float(df_traj[nh3_col].iloc[-1])
                        
                        # Konverziós ráta számítás: termelt ammónia / potenciális ammónia (3:1 H2:N2 kezdeti keverékből)
                        # A kezdeti 3:1 H2:N2 arányból adódó elméleti maximum vagy referenciamennyiség
                        initial_h2_n2_potential = postproc.get("initial_nh3_potential", 1.0) # vagy becsült referencia
                        if initial_h2_n2_potential > 0:
                            nh3_conversion = nh3_moles / initial_h2_n2_potential
            except Exception:
                pass

        rec["ammónia_mólszám"] = nh3_moles
        rec["ammónia_konverziós_ráta"] = nh3_conversion
        records.append(rec)

    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="Generate X-Y plots for simulation parameters and outputs.")
    parser.add_argument("--master-index", type=Path, default=MASTER_INDEX_PATH)
    parser.add_argument("--x", type=str, default="P_amb", help="X-axis column name")
    parser.add_argument("--y", type=str, default="energiaigény", help="Y-axis column name")
    parser.add_argument("--hue", type=str, default="nu_L", help="Color-coding column name")
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR / "custom_scatter.png")
    args = parser.parse_args()

    df = load_all_run_data(args.master_index)
    if df.empty:
        raise SystemExit("No valid data loaded from master index.")

    if args.x not in df.columns or args.y not in df.columns:
        raise SystemExit(f"Invalid columns. Available columns: {list(df.columns)}")

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(df[args.x], df[args.y], c=df[args.hue] if args.hue in df.columns else None,
                     cmap="viridis", edgecolor="k", alpha=0.8, s=50)
    
    if args.hue in df.columns:
        cbar = plt.colorbar(sc)
        cbar.set_label(args.hue)

    plt.xlabel(args.x)
    plt.ylabel(args.y)
    plt.title(f"{args.y} vs {args.x} (colored by {args.hue})")
    plt.grid(True, linestyle="--", alpha=0.6)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {args.output}")
    plt.show()


if __name__ == "__main__":
    main()