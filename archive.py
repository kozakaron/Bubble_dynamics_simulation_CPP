"""
archive.py

Finalizes raw solver output sitting in /tmp_data into permanent, cataloged
entries in /raw_data + master_index.csv, and provides the management
operations for that archive (cleaning up already-archived tmp_data
folders, wiping raw_data entirely, deleting one specific run).

MASTER_INDEX.CSV HEADER (12 columns):
    run_id, run_type, timestamp, mechanism, p_A, R_E, f, nu_L, P_amb,
    energy_demand, raw_data_path, metadata_path

DETECTION (content-based, not filename-based):
    tmp_data/<folder>/ contains output_*.csv        -> "sweep"  (whole folder = 1 entry)
    tmp_data/<folder>/ contains a *_opt.json file    -> "optimum" (that file only = 1 entry;
                                                         sibling numbered files, if any from
                                                         --keep-history, are exploration history,
                                                         not individually archived)
    tmp_data/<folder>/ otherwise, per .json file      -> "single" (each file = its own entry)

p_A / R_E / f / nu_L / P_amb ENCODING (same for all three run_types):
    - Const or plain scalar (single/optimum always fall in this case)
        -> converted to plot units (p_A: bar, positive magnitude; R_E: um;
           f: kHz; P_amb: bar; nu_L: cSt, derived from mu_L & rho_L)
    - LinearRange/LogRange (sweep only)
        -> "min:max:num_steps:lin" or "...:log", same units as above --
           this string is valid input to single_sweep.py's --pA/--RE or
           global_optimise.py's --pA-bounds/--RE-bounds directly.
    - nu_L specifically, if mu_L/rho_L are themselves swept (not expected
      in this project's normal use, but handled defensively): the string
      "non-constant (mu_L/rho_L swept)", since a derived ratio of two
      independently-varying quantities can't be reduced to one range.

Requires no other project scripts -- reads only the raw JSON/CSV files
solver_common.py-based scripts already produce.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import struct
from datetime import datetime
from pathlib import Path

TMP_DATA_DIR = Path("tmp_data")
RAW_DATA_DIR = Path("raw_data")
MASTER_INDEX_PATH = RAW_DATA_DIR / "master_index.csv"
SETTINGS_FILENAME = "bruteforce_parameter_study_settings.json"
ARCHIVED_MARKER = ".archived"

MASTER_INDEX_HEADER = [
    "run_id", "run_type", "timestamp", "mechanism",
    "p_A", "R_E", "f", "nu_L", "P_amb", "energy_demand",
    "raw_data_path", "metadata_path",
]


# --------------------------------------------------------------------------
# master_index.csv plumbing
# --------------------------------------------------------------------------

def ensure_master_index() -> None:
    """Creates /raw_data and master_index.csv (header only) if either
    doesn't exist yet. Safe to call unconditionally before any read/write."""
    RAW_DATA_DIR.mkdir(exist_ok=True)
    if not MASTER_INDEX_PATH.exists():
        with open(MASTER_INDEX_PATH, "w", newline="") as f:
            csv.writer(f).writerow(MASTER_INDEX_HEADER)
        print(f"Created {MASTER_INDEX_PATH} with header: {', '.join(MASTER_INDEX_HEADER)}")


def get_next_run_id() -> int:
    ensure_master_index()
    max_id = -1
    with open(MASTER_INDEX_PATH, newline="") as f:
        for row in csv.DictReader(f):
            try:
                max_id = max(max_id, int(row["run_id"]))
            except (KeyError, ValueError):
                continue
    return max_id + 1


def append_master_index_row(row: dict) -> None:
    ensure_master_index()
    with open(MASTER_INDEX_PATH, "a", newline="") as f:
        csv.DictWriter(f, fieldnames=MASTER_INDEX_HEADER).writerow(row)


# --------------------------------------------------------------------------
# Unit-converting parameter formatting (used by all three run_types)
# --------------------------------------------------------------------------

def _plain_value(field):
    """None if `field` is itself a swept Range; the raw numeric value
    otherwise -- handles a plain scalar (cpar mode), a Const dict, or a
    Range dict uniformly."""
    if isinstance(field, dict):
        return field["value"] if field.get("type") == "Const" else None
    return field


def format_param(field, converter) -> object:
    """Converts one physical field to its master_index representation:
    a converted plain number (Const/scalar case) or a
    "min:max:num_steps:scale" string (swept Range case)."""
    if isinstance(field, dict) and field.get("type") in ("LinearRange", "LogRange", "GeomRange"):
        scale = {"LinearRange": "lin", "LogRange": "log", "GeomRange": "geom"}[field["type"]]
        return f"{converter(field['start'])}:{converter(field['end'])}:{field['num_steps']}:{scale}"
    return converter(_plain_value(field))


def format_nu_L(mu_L_field, rho_L_field) -> object:
    mu_L = _plain_value(mu_L_field)
    rho_L = _plain_value(rho_L_field)
    if mu_L is None or rho_L is None:
        return "non-constant (mu_L/rho_L swept)"
    return mu_L / (rho_L * 1.0e-6)


def to_bar(pa_value: float) -> float:
    return pa_value / 1.0e5


def to_bar_abs(pa_value: float) -> float:
    return abs(pa_value) / 1.0e5


def to_um(m_value: float) -> float:
    return m_value * 1.0e6


def to_khz(hz_value: float) -> float:
    return hz_value / 1.0e3


def build_param_columns(settings: dict) -> dict:
    """settings: a flat dict of physical fields (works for both the sweep
    settings JSON and a single/optimum result's "cpar" section, since both
    use the same field names -- only whether each field is a scalar/Const
    or a Range differs)."""
    excitation_params = settings.get("excitation_params", [])
    p_A_field = excitation_params[0] if len(excitation_params) > 0 else None
    f_field = excitation_params[1] if len(excitation_params) > 1 else None
    return {
        "mechanism": settings.get("mechanism", ""),
        "p_A": format_param(p_A_field, to_bar_abs) if p_A_field is not None else "",
        "R_E": format_param(settings.get("R_E"), to_um),
        "f": format_param(f_field, to_khz) if f_field is not None else "",
        "nu_L": format_nu_L(settings.get("mu_L"), settings.get("rho_L")),
        "P_amb": format_param(settings.get("P_amb"), to_bar),
    }


# --------------------------------------------------------------------------
# sweep finalization
# --------------------------------------------------------------------------

def merge_output_csvs(output_files: list[Path], dest_path: Path) -> tuple[int, int, float | None]:
    n_total = n_success = 0
    best_energy_demand = None
    with open(dest_path, "w", newline="") as out_f:
        writer = None
        for src_path in output_files:
            with open(src_path, newline="") as in_f:
                reader = csv.DictReader(in_f)
                if writer is None:
                    writer = csv.DictWriter(out_f, fieldnames=reader.fieldnames)
                    writer.writeheader()
                for row in reader:
                    writer.writerow(row)
                    n_total += 1
                    if row.get("success", "").strip().lower() in ("true", "1"):
                        n_success += 1
                        try:
                            ed = float(row["energy_demand"])
                            if best_energy_demand is None or ed < best_energy_demand:
                                best_energy_demand = ed
                        except (KeyError, ValueError):
                            pass
    return n_total, n_success, best_energy_demand


def archive_sweep(folder: Path) -> bool:
    output_files = sorted(folder.glob("output_*.csv"))
    if not output_files:
        return False

    settings_path = folder / SETTINGS_FILENAME
    if not settings_path.exists():
        candidates = list(folder.glob("*.json"))
        if not candidates:
            print(f"  Skipping {folder.name}: no settings JSON found.")
            return False
        settings_path = candidates[0]

    with open(settings_path) as f:
        settings = json.load(f)
    # bruteforce_parameter_study_settings.json is written straight from
    # ParameterCombinator::to_json(), which wraps its fields under a
    # "parameter_study" key (parameter_combinator.cpp) -- unwrap it here the
    # same way single_sweep.py's own load_baseline() unwraps the input
    # template, or every field lookup below silently misses (settings.get(...)
    # returns None instead of KeyError) and format_param() crashes trying to
    # do arithmetic on None.
    settings = settings.get("parameter_study", settings)

    run_id = get_next_run_id()
    raw_csv_path = RAW_DATA_DIR / f"run_{run_id}.csv"
    raw_json_path = RAW_DATA_DIR / f"run_{run_id}.json"

    n_simulations, n_success, best_energy_demand = merge_output_csvs(output_files, raw_csv_path)
    shutil.copy2(settings_path, raw_json_path)
    with open(raw_json_path) as f:
        meta = json.load(f)
    meta["n_simulations"] = n_simulations
    meta["n_success"] = n_success
    with open(raw_json_path, "w") as f:
        json.dump(meta, f, indent=2)

    row = {
        "run_id": run_id, "run_type": "sweep",
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M"),
        **build_param_columns(settings),
        "energy_demand": best_energy_demand if best_energy_demand is not None else "",
        "raw_data_path": str(raw_csv_path), "metadata_path": str(raw_json_path),
    }
    append_master_index_row(row)
    print(f"  [sweep] run_id {run_id}: {n_success}/{n_simulations} successful "
          f"-> {raw_csv_path.name}, {raw_json_path.name}")
    return True


# --------------------------------------------------------------------------
# single / optimum finalization (parses the hybrid JSON+binary format)
# --------------------------------------------------------------------------

def parse_single_result(json_path: Path) -> tuple[dict, list[dict]]:
    """Returns (cpar_settings, trajectory_rows). trajectory_rows is a list
    of one dict per saved timestep: t, R, R_dot, T, <species...>,
    dissipated_energy, p_excitation, p_internal. Works whether the run was
    made with --save (many rows) or without it (2 rows: first/last)."""
    with open(json_path, "rb") as f:
        content = f.read()
    marker = b"<BINARY>"
    idx = content.find(marker)
    header = content[:idx] if idx != -1 else content
    data = json.loads(header.decode("utf-8").strip())

    cpar = data.get("cpar", {})
    if idx == -1:
        return cpar, []  # no trajectory data at all (e.g. a fatal preprocessing failure)

    S = data["sol"]["num_saved_steps"]
    D = data["sol"]["num_dim"]
    species_names = data.get("mechanism", {}).get("species_names", [])

    binary = content[idx + len(marker):]
    offset = 0
    import numpy as np
    t = np.frombuffer(binary[offset:offset + S * 8], dtype=np.float64); offset += S * 8
    x = np.frombuffer(binary[offset:offset + S * D * 8], dtype=np.float64).reshape(S, D); offset += S * D * 8
    p_excitation = np.frombuffer(binary[offset:offset + S * 8], dtype=np.float64); offset += S * 8
    p_internal = np.frombuffer(binary[offset:offset + S * 8], dtype=np.float64)

    rows = []
    for i in range(S):
        row = {"t": t[i], "R": x[i, 0], "R_dot": x[i, 1], "T": x[i, 2]}
        for j, name in enumerate(species_names):
            row[name] = x[i, 3 + j]
        row["dissipated_energy"] = x[i, -1]
        row["p_excitation"] = p_excitation[i]
        row["p_internal"] = p_internal[i]
        rows.append(row)
    return cpar, rows


def archive_single_result(json_path: Path, run_type: str) -> bool:
    try:
        cpar, rows = parse_single_result(json_path)
    except Exception as e:
        print(f"  Skipping {json_path.name}: could not parse ({e}).")
        return False

    with open(json_path, "rb") as f:
        content = f.read()
    marker_idx = content.find(b"<BINARY>")
    header = content[:marker_idx] if marker_idx != -1 else content
    data = json.loads(header.decode("utf-8").strip())
    postproc = data.get("postproc", {})

    run_id = get_next_run_id()
    raw_csv_path = RAW_DATA_DIR / f"run_{run_id}.csv"
    raw_json_path = RAW_DATA_DIR / f"run_{run_id}.json"

    if rows:
        with open(raw_csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    else:
        with open(raw_csv_path, "w") as f:
            f.write("")  # no trajectory available (e.g. fatal preprocessing failure)

    with open(raw_json_path, "w") as f:
        json.dump(data, f, indent=2)

    energy_demand = postproc.get("energy_demand", "")

    row = {
        "run_id": run_id, "run_type": run_type,
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M"),
        **build_param_columns(cpar),
        "energy_demand": energy_demand,
        "raw_data_path": str(raw_csv_path), "metadata_path": str(raw_json_path),
    }
    append_master_index_row(row)
    print(f"  [{run_type}] run_id {run_id}: {len(rows)} timestep(s) "
          f"-> {raw_csv_path.name}, {raw_json_path.name}")
    return True


# --------------------------------------------------------------------------
# top-level dispatch
# --------------------------------------------------------------------------

def archive_folder(folder: Path) -> bool:
    """Returns True if anything in this folder was archived (so it's safe
    to mark .archived), False if the folder was empty/unrecognized."""
    if list(folder.glob("output_*.csv")):
        return archive_sweep(folder)

    opt_files = list(folder.glob("*_opt.json"))
    if opt_files:
        return archive_single_result(opt_files[0], "optimum")

    json_files = sorted(folder.glob("*.json"))
    if not json_files:
        return False
    archived_any = False
    for json_path in json_files:
        if archive_single_result(json_path, "single"):
            archived_any = True
    return archived_any


def archive_all() -> int:
    """Scans /tmp_data for un-archived folders and finalizes each. Returns
    the count of folders successfully archived."""
    if not TMP_DATA_DIR.exists():
        print(f"'{TMP_DATA_DIR}' does not exist -- nothing to do.")
        return 0

    count = 0
    for folder in sorted(p for p in TMP_DATA_DIR.iterdir() if p.is_dir()):
        if (folder / ARCHIVED_MARKER).exists():
            continue
        print(f"Processing {folder.name} ...")
        try:
            if archive_folder(folder):
                (folder / ARCHIVED_MARKER).touch()
                count += 1
            else:
                print(f"  Nothing recognized in {folder.name} -- left unmarked.")
        except Exception as e:
            print(f"  ERROR processing {folder.name}: {e}")

    print(f"\nDone. Archived {count} folder(s)/file(s).")
    return count


# --------------------------------------------------------------------------
# management functions
# --------------------------------------------------------------------------

def clean_archived() -> int:
    """Deletes every tmp_data subfolder carrying the .archived marker, plus
    every loose *.log file sitting directly in tmp_data/ -- the per-worker
    log files batch_global_optimise.py's multi-core mode writes there
    (<run_tag>_cpu<id>.log). Those logs live as siblings of the
    per-combination folders, not inside any of them, so they're not covered
    by the .archived-marker folder scan above; they carry no marker of
    their own, so every *.log file directly in tmp_data/ is deleted
    unconditionally each time this runs.
    Returns the count of folders + log files deleted."""
    if not TMP_DATA_DIR.exists():
        return 0
    count = 0
    for folder in sorted(p for p in TMP_DATA_DIR.iterdir() if p.is_dir()):
        if (folder / ARCHIVED_MARKER).exists():
            shutil.rmtree(folder)
            print(f"Deleted {folder}")
            count += 1
    for log_path in sorted(TMP_DATA_DIR.glob("*.log")):
        log_path.unlink()
        print(f"Deleted {log_path}")
        count += 1
    print(f"Cleaned up {count} archived folder(s)/log file(s).")
    return count


def wipe_raw_data(confirm: bool = False) -> None:
    """Deletes EVERYTHING in /raw_data and resets master_index.csv to just
    its header. Irreversible -- requires confirm=True (CLI: --yes-delete-everything)."""
    if not confirm:
        raise RuntimeError(
            "wipe_raw_data() requires confirm=True (CLI: --yes-delete-everything). "
            "This deletes every archived result permanently."
        )
    if RAW_DATA_DIR.exists():
        shutil.rmtree(RAW_DATA_DIR)
    ensure_master_index()
    print(f"Wiped {RAW_DATA_DIR} and reset master_index.csv.")


def delete_run(run_id: int) -> bool:
    """Deletes one run's files and its master_index.csv row. Leaves a gap
    in run_id numbering deliberately -- see module discussion: renumbering
    would silently invalidate any external reference (e.g. a saved plot
    named after the old run_id)."""
    ensure_master_index()
    rows = []
    target_row = None
    with open(MASTER_INDEX_PATH, newline="") as f:
        for row in csv.DictReader(f):
            if row["run_id"] == str(run_id):
                target_row = row
            else:
                rows.append(row)

    if target_row is None:
        print(f"run_id {run_id} not found in {MASTER_INDEX_PATH}.")
        return False

    for path_key in ("raw_data_path", "metadata_path"):
        p = Path(target_row[path_key])
        if p.exists():
            p.unlink()
            print(f"Deleted {p}")

    with open(MASTER_INDEX_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MASTER_INDEX_HEADER)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Removed run_id {run_id} from {MASTER_INDEX_PATH} (gap left, not renumbered).")
    return True


def main():
    parser = argparse.ArgumentParser(description="Archive tmp_data results into raw_data + master_index.csv.")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("run", help="scan tmp_data and archive everything new (default)")
    sub.add_parser("clean-archived", help="delete tmp_data folders already marked .archived")

    wipe_parser = sub.add_parser("wipe-raw-data", help="delete ALL of raw_data (irreversible)")
    wipe_parser.add_argument("--yes-delete-everything", action="store_true", required=True)

    delete_parser = sub.add_parser("delete-run", help="delete one run_id's files and index row")
    delete_parser.add_argument("run_id", type=int)

    args = parser.parse_args()

    if args.command == "run":
        archive_all()
    elif args.command == "clean-archived":
        clean_archived()
    elif args.command == "wipe-raw-data":
        wipe_raw_data(confirm=args.yes_delete_everything)
    elif args.command == "delete-run":
        delete_run(args.run_id)


if __name__ == "__main__":
    main()
