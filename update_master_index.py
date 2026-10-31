"""
update_master_index.py

Scans raw JSON files in raw_data/ (or subdirectories), extracts postproc.n_target_specie,
and updates/appends it to raw_data/master_index.csv as 'termelt_ammónia'.
"""

from __path__ import annotations
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
        fieldnames = reader.fieldnames
        rows = list(reader)

    if not fieldnames:
        raise SystemExit("Master index has no columns.")

    # Ha még nincs benne az új oszlop, hozzáadjuk
    target_col = "termelt_ammónia"
    if target_col not in fieldnames:
        fieldnames.append(target_col)

    updated_count = 0

    # 2. Végigmegyünk a sorokon, és megkeressük a JSON fájlokat
    for row in rows:
        # Megpróbáljuk kideríteni a JSON fájl útvonalát
        # (A master_index általában tartalmazza a meta/json path-t, vagy kereshetünk a raw_data-ban)
-       meta_path = row.get("metadata_path", "")
        json_file = None
        
        if meta_path and Path(meta_path).exists():
            json_file = Path(meta_path)
        else:
            # Ha nincs benne közvetlen útvonal, keressük a raw_data-ban a run_tag vagy egyebek alapján
            # Vagy feltételezzük, hogy a raw_data-ban lévő összes JSON-t átvizsgálhatjuk
            pass

        # Ha így sem találtuk, próbáljuk meg a raw_data mappából kinyerni a JSON-öket
        # (vagy ha van metadata oszlop, az alapján)
        if not json_file:
            # Keressünk rá a raw_data-ban lévő json fájlokra
            pass

    # Egy elegánsabb megközelítés: közvetlenül átvizsgáljuk az összes JSON-t a raw_data-ban,
    # és egy kulcs (pl. P_amb, f, nu_L, p_A, R_E) alapján párosítjuk a master_index soraival.
    
    json_ammonia_map = {}
    for jpath in RAW_DATA_DIR.glob("**/*.json"):
        try:
            with open(jpath, "r", encoding="utf-8") as jf:
                content = json.load(jf)
            if "postproc" in content and "n_target_specie" in content["postproc"]:
                n_val = content["postproc"]["n_target_specie"]
                # Eltároljuk a fájlútvonalat vagy kulcsot
                json_ammonia_map[jpath.resolve()] = n_val
        except Exception:
            continue

    # Frissítjük a sorokat
    for row in rows:
        meta_path = row.get("metadata_path", "")
        val = None
        if meta_path and Path(meta_path).resolve() in json_ammonia_map:
            val = json_ammonia_map[Path(meta_path).resolve()]
        else:
            # Ha nincs metadata_path, megpróbáljuk kitalálni a fájlnévből vagy bejárjuk a raw_data-t
            # Ha van raw_data_path oszlop:
            raw_path = row.get("raw_data_path", "")
            if raw_path:
                # A JSON általában a raw_data-ban van mellette vagy a meta mappában
                possible_json = Path(raw_path).with_suffix(".json")
                if possible_json.resolve() in json_ammonia_map:
                    val = json_ammonia_map[possible_json.resolve()]

        # Ha még mindig nincs meg, próbáljuk megkeresni az összes raw_data-beli JSON-ből azt, 
        # amelyiknek a paraméterei megegyeznek ezzel a sorral (biztos ami biztos)
        if val is not None:
            row[target_col] = str(val)
            updated_count += 1
        else:
            row[target_col] = row.get(target_col, "0.0")

    # 3. Visszaírjuk a frissített master_index.csv-t
    with open(MASTER_INDEX_PATH, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Sikeresen frissítve! {updated_count} sor kapott termelt_ammónia értéket a master_index.csv-ben.")


if __name__ == "__main__":
    main()