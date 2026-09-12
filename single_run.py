"""
single_run.py

Runs exactly ONE simulation (the solver's lighter --run mode, not a
degenerate 1-point --parameter_study) for a specific (nu_L, f, P_amb, R_E,
p_A) combination. Intended for use as the objective function inside an
external optimizer (e.g. scipy.optimize.direct) that needs many fast,
individually-addressable evaluations rather than a batch grid.

Requires solver_common.py in the same folder — reuses its liquid-property
formulas, JSON-comment stripping, and executable-path resolution directly,
no duplicated logic. Deliberately does NOT depend on single_sweep.py: the
sweep pathway and the single-run/optimization pathway only ever shared
generic, non-sweep-specific utilities, which now live in solver_common.py
so the two pathways can be used, understood, and modified independently.

Usage:
    python single_run.py --run-tag eval_001 --nuL 0.65 --f 10 --Pamb 1.0 --pA 8.0 --RE 150.0

Units match solver_common.py's convention: nu_L [cSt], f [kHz], P_amb [bar],
p_A [bar] (POSITIVE magnitude, negated internally), R_E [um].

Output: tmp_data/<run_tag>/<run_tag><id>.json — one file per call, id auto-
incrementing (0, 1, 2, ...) within the shared <run_tag> folder, so calling
this repeatedly with the same run_tag (e.g. from an optimizer) keeps every
evaluation instead of overwriting the previous one. Each file contains
both the exact settings used and the result. Without --save, only the
first/last timestep are recorded (still enough for energy_demand, which
the solver computes from the final state regardless of --save — see
read_energy_demand below); with --save, the complete time-series
trajectory is included too, at the cost of a much larger file. Since an
optimizer may call this hundreds of times, --save defaults to OFF; turn it
on only for a specific point you want to inspect afterward (e.g. the final
optimum).

NOTE ON liquid_eos_params: EVERY field is a plain scalar/list in --run
(cpar) mode, including this one — despite this project's example_cpar.json
baseline template showing it Const-wrapped (like the parameter_study
format). That template syntax is itself wrong for --run mode: feeding it
through unchanged produces "Expected floating point number... Using
default value" warnings followed by a fatal "number of liquid EOS
parameters must be equal to 4: params_list={}" error, confirmed directly
against the real solver, not assumed from the template. This script
unwraps every liquid property (including liquid_eos_params) to plain
values before writing the config.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import NamedTuple

import solver_common

CPAR_BASELINE_PATH = Path("config_data") / "example_cpar.json"

# --- condensing the solver's raw diagnostic text down to one short tag -----
#
# The chemistry/ODE side of this project is numerically hard by nature, so
# CVODE warnings (its own, e.g. "[WARNING][rank 0][cvode.c:1457][CVode]
# Internal t = ... t + h = t ...") and this project's own Error/ErrorHandler
# messages (format: "<time>: WARNING (type) <message> in <file>:<line>:
# <function>();", printed to stderr in color by default -- see
# ErrorHandler::log_error/print_when_log in common.cpp) are both routine and
# frequent, not exceptional. Dumping the full captured stdout/stderr for
# every such occurrence during an optimization run (hundreds of evaluations)
# floods the screen; these helpers instead pick out the single most relevant
# line and shrink it to fit on the same line as the eval's normal result.
_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")                       # strip color codes
_SEVERITY_RE = re.compile(r"\b(WARNING|ERROR)\b", re.IGNORECASE)
_LEADING_BRACKETS_RE = re.compile(r"^(?:\[[^\]]*\])+\s*")      # SUNDIALS' own "[WARNING][rank 0][...]" style
_LEADING_TIMESTAMP_RE = re.compile(r"^\d+(?:\.\d+)?:\s*")      # this project's own "12.345: ..." style


def _condense_diagnostics(result: subprocess.CompletedProcess, max_len: int = 70) -> str:
    """Reduces a (possibly large, multi-line) captured stdout+stderr down to
    ONE short 'label: message' tag, e.g.
    'warning: Internal t = 86827.3 and h = 5.10993e-14 are such th...'.

    Picks the LAST line mentioning WARNING/ERROR (most relevant when a
    late fatal error follows earlier benign warnings, which is the normal
    shape of a genuinely failing run), strips ANSI color codes and the
    cosmetic prefixes both diagnostic styles use (bracketed file/line tags
    for CVODE's own messages, a leading timestamp for this project's own
    Error format), then truncates. Falls back to the last line of any kind,
    then to just the exit code, if nothing matched. This is deliberately
    lossy -- it's a screen-width hint, not a replacement for opening the
    run's own file or rerunning with verbose=True for the full picture."""
    text = _ANSI_RE.sub("", (result.stderr or "") + "\n" + (result.stdout or ""))
    lines = [" ".join(line.split()) for line in text.splitlines() if line.strip()]
    tagged = [line for line in lines if _SEVERITY_RE.search(line)]
    chosen = tagged[-1] if tagged else (lines[-1] if lines else None)
    if chosen is None:
        return f"exit code {result.returncode}"

    severity_match = _SEVERITY_RE.search(chosen)
    label = severity_match.group(1).lower() if severity_match else "note"

    rest = _LEADING_BRACKETS_RE.sub("", chosen)
    rest = _LEADING_TIMESTAMP_RE.sub("", rest)
    if severity_match:
        # Drop a redundant leading "WARNING"/"ERROR (type)" now that it's
        # already captured as the label.
        rest = re.sub(rf"^{severity_match.group(1)}\s*(?:\([^)]*\))?\s*", "", rest, flags=re.IGNORECASE)
    rest = rest.strip() or chosen

    if len(rest) > max_len:
        rest = rest[:max_len - 1].rstrip() + "…"
    return f"{label}: {rest}"


class RunResult(NamedTuple):
    """config_path: same as before -- the file to read energy_demand from,
    or None for a dry run.
    condensed_note: None on a clean run (or whenever verbose=True, since
    that path already streams/prints everything live itself); otherwise a
    short 'label: message' string summarizing whatever the solver printed,
    for a caller to fold into its own single-line progress output instead
    of dumping the raw text (see _condense_diagnostics)."""
    config_path: Path | None
    condensed_note: str | None


def next_available_id(run_dir: Path, run_tag: str) -> int:
    """Scan run_dir for existing '<run_tag><N>.json' files and return the
    next unused N (0 if none exist yet). Derived from the filesystem, not
    any counter kept in memory — so this naturally resumes correctly even
    if the calling process (e.g. an optimizer) is restarted partway
    through a run, rather than colliding with or skipping past what's
    already on disk.

    Assumes calls for a given run_tag happen sequentially, not
    concurrently: two processes scanning at the same instant could both
    compute the same "next" id and one would overwrite the other. Fine for
    a single-process optimizer loop; not safe for parallel workers sharing
    one run_tag without extra coordination."""
    pattern = re.compile(rf"^{re.escape(run_tag)}(\d+)\.json$")
    ids = [int(m.group(1)) for p in run_dir.glob(f"{run_tag}*.json")
           if (m := pattern.match(p.name))]
    return max(ids, default=-1) + 1


def load_cpar_baseline() -> dict:
    with open(CPAR_BASELINE_PATH) as f:
        raw_text = f.read()
    return json.loads(solver_common.strip_json_comments(raw_text))["cpar"]


def build_cpar_config(nu_L: float, f_khz: float, P_amb_bar: float, pA_bar: float, RE_um: float) -> dict:
    config = load_cpar_baseline()

    config["R_E"] = RE_um * 1.0e-6                                   # um -> m
    config["P_amb"] = P_amb_bar * 1.0e5                              # bar -> Pa
    config["excitation_params"] = [-abs(pA_bar) * 1.0e5, f_khz * 1.0e3]  # bar,kHz -> Pa,Hz (negative sign applied)
    config["target_specie"] = "NH3"

    props = solver_common.compute_liquid_properties(nu_L)
    config["mu_L"] = props["mu_L"]["value"]
    config["rho_L"] = props["rho_L"]["value"]
    config["c_L"] = props["c_L"]["value"]
    config["surfactant"] = props["surfactant"]["value"]
    # Plain floats, NOT Const-wrapped -- confirmed directly against the
    # solver's own preprocessing error, not just the baseline template's
    # syntax (which turned out to be wrong for --run/cpar mode: it uses
    # Const-wrapped objects here, but the solver expects a plain array and
    # silently drops each entry to "default", leaving zero usable values).
    config["liquid_eos_params"] = [p["value"] for p in props["liquid_eos_params"]]

    return config


def run_single(run_tag: str, nu_L: float, f_khz: float, P_amb_bar: float, pA_bar: float, RE_um: float,
                tmax: float = 1.0, timeout: float = 60.0, save: bool = False,
                dry_run: bool = False, verbose: bool = True) -> RunResult:
    """Returns a RunResult(config_path, condensed_note) -- config_path is
    None only for a dry run. Every call for a given run_tag gets its own
    numbered file inside the same shared folder (tmp_data/<run_tag>/
    <run_tag><id>.json, id = 0, 1, 2, ...) rather than overwriting the
    previous evaluation — see next_available_id.

    verbose=True (the default, matching this script's own standalone CLI)
    prints the "Wrote config"/"Solver command" lines and lets the solver's
    own console output stream straight through live, same as before --
    condensed_note is always None on this path, since the solver's own
    output already told you everything live.
    verbose=False (used by global_optimise.py, which calls this hundreds
    of times per optimization) suppresses those prints entirely AND does
    NOT dump the solver's captured stdout/stderr itself -- chemistry being
    numerically hard means CVODE/solver warnings are frequent and often
    voluminous, and a full dump per evaluation would flood the screen.
    Instead, whenever the run exited non-zero or printed anything
    mentioning WARNING/ERROR, condensed_note carries a single short
    'label: message' summary (see _condense_diagnostics) for the caller to
    fold into its own one-line-per-eval output; on a clean run it's None."""
    run_dir = solver_common.TMP_DATA_DIR / run_tag
    run_dir.mkdir(parents=True, exist_ok=True)

    eval_id = next_available_id(run_dir, run_tag)
    filename = f"{run_tag}{eval_id}.json"

    config = build_cpar_config(nu_L, f_khz, P_amb_bar, pA_bar, RE_um)
    config_path = run_dir / filename
    with open(config_path, "w") as f:
        json.dump({"cpar": config}, f, indent=2)
    if verbose:
        print(f"Wrote config: {config_path}")

    cmd = [
        solver_common.solver_executable(solver_common.SOLVER_DIR),
        "--run", f"../tmp_data/{run_tag}/{filename}",
        "--tmax", str(tmax),
        "--timeout", str(timeout),
    ]
    if save:
        cmd.append("--save")
    if verbose:
        print("Solver command:", " ".join(cmd), f"(cwd={solver_common.SOLVER_DIR})")

    if dry_run:
        print("Dry run — solver not executed.")
        return RunResult(None, None)

    result = subprocess.run(
        cmd, cwd=solver_common.SOLVER_DIR,
        capture_output=not verbose, text=True,
    )

    condensed_note = None
    if verbose:
        # Streamed live already; just flag a non-zero exit, same as before.
        if result.returncode != 0:
            print(f"Note: solver exited with code {result.returncode} — check {config_path} "
                  f"for the postprocessing result (energy_demand is set to a large sentinel "
                  f"value on failure, not omitted — see read_energy_demand).")
    else:
        diagnostic_text = (result.stdout or "") + (result.stderr or "")
        if result.returncode != 0 or _SEVERITY_RE.search(diagnostic_text):
            condensed_note = _condense_diagnostics(result)

    return RunResult(config_path, condensed_note)


def read_energy_demand(result_path: Path) -> float:
    """Read back just the postprocessed energy_demand [MJ/kg = GJ/t] from a
    --run result file — works whether or not --save was used, since
    energy_demand is computed from the final state regardless. Does not
    need to parse the <BINARY> trajectory block at all.

    Two distinct failure shapes are handled, both mapped to float('inf'):
    - The solver's own "ran, but yielded ~zero target species" sentinel is
      DBL_MAX (a huge but finite float) -- returned as-is, no special case
      needed, since it already sorts as "worse than any real result".
    - A FATAL error before the simulation ever ran (e.g. a preprocessing
      error) leaves energy_demand at its C++ default, which is genuine
      infinity -- but JSON has no representation for that, so
      nlohmann::json serializes it as `null`, which Python reads as None.
      Without handling this, callers get a silent None instead of a
      number, surfacing as a confusing crash wherever it's first used
      (e.g. an f-string format), far from the actual cause."""
    with open(result_path, "rb") as f:
        content = f.read()
    marker = b"<BINARY>"
    idx = content.find(marker)
    header = content[:idx] if idx != -1 else content
    data = json.loads(header.decode("utf-8").strip())

    energy_demand = data.get("postproc", {}).get("energy_demand")
    if energy_demand is None:
        return float("inf")
    return energy_demand


def main():
    parser = argparse.ArgumentParser(description="Run a single simulation for one (nu_L, f, P_amb, p_A, R_E) point.")
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--nuL", type=float, required=True, help="[cSt]")
    parser.add_argument("--f", type=float, required=True, help="[kHz]")
    parser.add_argument("--Pamb", type=float, required=True, help="[bar]")
    parser.add_argument("--pA", type=float, required=True, help="positive [bar] magnitude")
    parser.add_argument("--RE", type=float, required=True, help="[um]")
    parser.add_argument("--tmax", type=float, default=1.0)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--save", action="store_true", help="save the full time-series trajectory")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    result = run_single(
        run_tag=args.run_tag, nu_L=args.nuL, f_khz=args.f, P_amb_bar=args.Pamb,
        pA_bar=args.pA, RE_um=args.RE, tmax=args.tmax, timeout=args.timeout,
        save=args.save, dry_run=args.dry_run,
    )
    if result.config_path is not None:
        energy_demand = read_energy_demand(result.config_path)
        print(f"energy_demand = {energy_demand:.4g} GJ/t")


if __name__ == "__main__":
    main()
