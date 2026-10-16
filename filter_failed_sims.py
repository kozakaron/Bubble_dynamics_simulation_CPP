from __future__ import annotations

import csv
import re
from pathlib import Path

TMP_DATA_DIR = Path("tmp_data")
OUTPUT_COMBOS_PATH = Path("combos.csv")


def parse_log_files(tmp_dir: Path) -> list[dict]:
    failed_combos = []
    
    if not tmp_dir.exists():
        print(f"Directory {tmp_dir} does not exist.")
        return failed_combos

    # Keresés a batch_*_cpu*.log fájlokban (vagy bármilyen .log fájlban a tmp_data-ban)
    log_files = list(tmp_dir.glob("*cpu*.log")) + list(tmp_dir.glob("*.log"))
    if not log_files:
        print(f"No log files found in {tmp_dir}.")
        return failed_combos

    # Minták a log sorok elemzésére, pl:
    # [batch_low_freq7] nu_L=0.65 f=0.1 P_amb=50.0 -> ... | status=OK | reason=... | boundary_warning=...
    pattern = re.compile(
        r"nu_L\s*=\s*([0-9.]+)\s+f\s*=\s*([0-9.]+)\s+P_amb\s*=\s*([0-9.]+)",
        re.IGNORECASE
    )

    for file_path in log_files:
        print(f"Scanning log file: {file_path}")
        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    # Csekkoljuk, hogy van-e benne hiba, nem konvergálás vagy boundary warning
                    line_lower = line.lower()
                    is_problematic = (
                        "not converged" in line_lower or
                        "maximum number of iterations has been exceeded" in line_lower or
                        #"boundary_warning" in line_lower or
                        "energy_demand=1e+12" in line_lower or
                        "boundary_warning=p_A" in line_lower
                        #"at/near a bound" in line_lower or
                        #"at/near upper bound" in line_lower
                    )

                    if is_problematic:
                        match = pattern.search(line)
                        if match:
                            nu_L_val = float(match.group(1))
                            f_val = float(match.group(2))
                            P_amb_val = float(match.group(3))
                            
                            failed_combos.append({
                                "nu_L": nu_L_val,
                                "P_amb": P_amb_val,
                                "f": f_val
                            })
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    return failed_combos


def save_to_combos_csv(combos: list[dict], output_path: Path):
    if not combos:
        print("No problematic simulations found in log files.")
        return

    # Deduplikálás (nu_L, P_amb, f) alapján
    unique_keys = set()
    unique_combos = []
    for c in combos:
        key = (round(c["nu_L"], 6), round(c["P_amb"], 6), round(c["f"], 6))
        if key not in unique_keys:
            unique_keys.add(key)
            unique_combos.append(c)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["nu_L", "P_amb", "f"])
        writer.writeheader()
        for c in unique_combos:
            writer.writerow(c)
            
    print(f"Sikeresen kiírva {len(unique_combos)} egyedi hibás/határeset kombináció ide: {output_path}.")


def main():
    print(f"Keresés a {TMP_DATA_DIR} logfájljaiban...")
    combos = parse_log_files(TMP_DATA_DIR)
    save_to_combos_csv(combos, OUTPUT_COMBOS_PATH)


if __name__ == "__main__":
    main()