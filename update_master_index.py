"""
update_master_index.py

Scans raw JSON files in raw_data/, extracts postproc.n_target_specie,
and updates/appends it to raw_data/master_index.csv as 'termelt_ammónia'.
"""

from __future__ import annotations
import csv
import json
from pathlib import Path

MASTER_INDEX_PATH = Path("raw_data") / "master_index.csv"
RAW_DATA_DIR = Path("raw_data")


def main():
    if not MASTER_INDEX_PATH.exists():
        raise SystemExit(f"{MASTER_INDEX_PATH} does not exist.")

    # 1. Beolvassuk a master_index.csv-t
    with open(MASTER_INDEX_PATH, mode="r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    if not fieldnames:
        raise SystemExit("Master index has no columns.")

    target_col = "termelt_ammónia"
    if target_col not in fieldnames:
        fieldnames.append(target_col)

    # 2. Összegyűjtjük a JSON fájlokból az n_target_specie értékeket az útvonaluk alapján
    json_ammonia_map = {}
    for jpath in RAW_DATA_DIR.glob("**/*.json"):
        try:
            with open(jpath, "r", encoding="utf-8") as jf:
                content = json.load(jf)
            if "postproc" in content and "n_target_specie" in content["postproc"]:
                n_val = content["postproc"]["n_target_specie"]
                json_ammonia_map[jpath.resolve()] = n_val
        except Exception:
            continue

    updated_count = 0

    # 3. Frissítjük a sorokat
    for row in rows:
        val = None
        meta_path = row.get("metadata_path", "")
        if meta_path:
            p_meta = Path(meta_path)
            if p_meta.resolve() in json_ammonia_map:
                val = json_ammonia_map[p_meta.resolve()]

        if val is None:
            raw_path = row.get("raw_data_path", "")
            if raw_path:
                p_raw = Path(raw_path).with_suffix(".json")
                if p_raw.resolve() in json_ammonia_map:
                    val = json_ammonia_map[p_raw.resolve()]

        if val is not None:
            row[target_col] = str(val)
            updated_count += 1
        else:
            row[target_col] = row.get(target_col, "0.0")

    # 4. Visszaírjuk a frissített master_index.csv-t
    with open(MASTER_INDEX_PATH, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Sikeresen frissítve! {updated_count} sor kapott termelt_ammónia értéket a master_index.csv-ben.")


if __name__ == "__main__":
    main()