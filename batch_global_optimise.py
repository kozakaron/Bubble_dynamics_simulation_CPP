"""
batch_global_optimise.py

Runs global_optimise.py's differential-evolution optimizer once for every
(nu_L, f, P_amb) combination listed in a CSV file, and prints/logs the
found optimum for each one as it finishes.

Requires solver_common.py, single_run.py, and global_optimise.py in the
same folder. No other scripts needed.

INPUT CSV: a single file with exactly the header row "nu_L,f,P_amb"
(units as everywhere else in this project: nu_L [cSt], f [kHz], P_amb
[bar]), given via --params-csv and resolved inside tmp_data/ -- e.g.
--params-csv combos.csv means tmp_data/combos.csv, not a path relative to
wherever this script happens to be run from. Blank lines and lines whose
first non-blank character is '#' are ignored, so a file can carry comments
between groups of rows.

--combo-mode controls how that same 3-column shape is turned into the list
of (nu_L, f, P_amb) triples that actually get optimized:

    --combo-mode cartesian (the default, replaces the old
    --nuL/--f/--Pamb comma-separated lists): each of the three columns is
    its own independent list of values, combined via a full cross
    product -- e.g. 3 nu_L values x 2 f values x 4 P_amb values = 24
    combinations. The three columns don't need the same number of
    entries: since a CSV file is still physically rectangular, a shorter
    column just leaves its remaining cells blank in whichever rows the
    longer column(s) still need. Example (3 nu_L, 2 f, 1 P_amb -> 6
    combinations):
        nu_L,f,P_amb
        0.5,10,1.0
        0.65,20,
        1.0,,

    --combo-mode raw: every row is used exactly as an already-paired
    (nu_L, f, P_amb) triple, no cross product at all -- so all three
    columns must have a value in every row (a blank cell here is an
    error, not "this column is shorter"). This is the mode for a list of
    specific points from anywhere other than a plain grid -- e.g. ones
    picked out by a separate analysis or an AI agent -- where nu_L, f,
    and P_amb must stay matched row-by-row rather than cross-multiplied.
    Example (exactly 3 combinations, not 27):
        nu_L,f,P_amb
        0.65,10,1.0
        1.2,35,2.5
        0.9,22,1.8

Combinations are tagged plainly by their position in the resulting list:
<run_tag>0, <run_tag>1, ... -- matching the project's existing
incrementing convention (single_run.py's own <run_tag><id>.json, the
solver's current_run<N>) rather than encoding the parameter values into
the name. Each log line records the (nu_L, f, P_amb) combination it
corresponds to, so nothing is lost by keeping the tag itself plain.

MULTI-CORE: combinations run concurrently, one per worker OS process (see
--workers, default: all logical cores). This parallelizes ACROSS
combinations only -- each combination's own differential_evolution run
stays exactly as sequential internally as before (workers=1 inside DE
itself is untouched, and unrelated to this), which is what single_run.py's
own auto-incrementing-filename scheme requires (see its docstring) and is
also the natural unit of concurrency here, since each `--run` solver
invocation is already single-threaded internally (checked directly
against the C++ source: --cpu/thread-pool options are read only in the
solver's --parameter_study branch, never in --run mode). Since every
combination already only ever touches its own tmp_data/<combo_tag>/
folder for its actual result files, running several at once needs no
locking or coordination between them there -- confirmed the same way, not
assumed. A worker pool with --workers N processes is created once and
reused for every combination in the batch: with more combinations than
workers, a single worker process runs several combinations one after
another over its lifetime, in whatever order the pool hands them out --
this is exactly the shape LOG FILES below is built around.

SCREEN OUTPUT: exactly one line per combination, printed by the main
process as soon as that combination's worker finishes (or fails) --
nothing is printed while any combination is still running (each worker
calls global_optimise.py with quiet=True, and only ever returns its result
to the main process rather than printing anything itself). Only the main
process ever calls print(), so concurrent workers can never interleave or
garble each other's output. One consequence of real concurrency: lines now
arrive in FINISHING order, not necessarily the [i/N] order combinations
were submitted in -- each line is still fully self-identifying (its own
[i/N] index and tag), so nothing is ambiguous, but don't expect strict
top-to-bottom ordering the way a sequential run would give you.

LOG FILES: one file per WORKER, not per combination --
tmp_data/<run_tag>_cpu<worker_id>.log, worker_id running 0, 1, ...,
--workers-1. Each worker is assigned its id exactly once, when the pool
starts it up (see _init_worker), and keeps that same id -- and the same
open log file -- for as long as the pool lives, appending one line to it
(same content/format this script has always used) every time it finishes
a combination. This mirrors how the C++ solver's own --parameter_study
mode gives each of ITS worker threads its own log file, reused across
every simulation that thread happens to run, rather than one shared file
or one file per simulation. It's safe without any locking for a subtler
reason than "each combination has its own file": a given worker process
only ever runs one combination at a time, so two writes to the same
worker's file can never actually race each other, even though that file
ends up holding lines for many different combinations over the course of
the batch, interleaved in whatever order they landed on that worker (each
line is still self-identifying via its own leading [combo_tag], so this
is never ambiguous to read). There is no single merged tmp_data/<run_tag>.log
for the whole batch, and no longer a log file inside each combo's own
tmp_data/<combo_tag>/ folder either -- use this function's returned list
of per-combo dicts (which carries each row's log_path, i.e. which
worker's file its line landed in) to find any particular combination's
line, or just grep across all tmp_data/<run_tag>_cpu*.log files.

If one combination fails outright (an exception, before any optimum was
found), it's still logged (to its worker's file) and the batch continues
with the remaining combinations -- check the "status" field in the
returned list, or grep the log files, afterward. Each combination's own
tmp_data/<combo_tag>/ folder is still cleaned up to just <combo_tag>_opt.json
by global_optimise.py itself (pass --keep-history to disable that for
every combination) -- this is unaffected by the logging change above,
since it's about the combination's own result files, not its log line.
Per-worker log files are opened in APPEND mode, never truncated -- so
reusing a --run-tag never loses the log history from an earlier
invocation of it, which is what makes RESTART / RESUME below safe.

RESTART / RESUME: safe to interrupt (Ctrl-C, a crash, a power outage) and
simply re-run the exact same command later, as long as you reuse the same
--run-tag AND keep --params-csv/--combo-mode unchanged since the
interrupted run. There's no separate manifest of "what's already done" --
this only ever looks at what's actually on disk for one combination,
checked right before that combination would otherwise run: does
tmp_data/<combo_tag>/<combo_tag>_opt.json already exist? If so, and it
parses cleanly, and its own recorded (nu_L, f, P_amb) -- read back out of
its "cpar" section, same fields/units single_run.py's build_cpar_config()
writes them in -- match what this index computes from the CURRENT CSV
(within a tight rel_tol=1e-6, meant only to absorb floating-point/unit-
conversion/JSON round-trip noise, not real parameter differences), that
combination is skipped entirely: no solver subprocess, no
differential_evolution run, just the already-known result reused for its
row/summary line. A skipped combination does NOT get a new line written
to any worker's log file -- its outcome is already durably recorded in
whichever log captured it the first time. If the file is missing, or
exists but can't be parsed (e.g. left truncated by a crash mid-write),
it's treated the same as "not done yet" and is just (re)computed
normally. Pass --force to ignore all of this and recompute every
combination from scratch regardless of what's already on disk.

If an existing result file parses fine but its parameters DON'T match the
current CSV, that's NOT treated as "not done" -- it raises
ComboMismatchError instead (see that class's docstring), which stops the
whole batch rather than just that one combination: with no manifest, a
changed CSV can silently shift what every later index means, so
continuing past a detected mismatch risks quietly reusing or overwriting
results under the wrong parameters throughout the rest of the run.
Combinations already running in other worker processes at that moment
finish naturally (they can't be interrupted mid-run); nothing new is
submitted after the abort. This is the one hard rule resume relies on:
never edit or rename params-csv for a run_tag you intend to resume --
start a new --run-tag (with its own CSV) for a new sweep instead.

p_A's upper bound defaults to P_amb + 0.5*f (see global_optimise.py) --
left as None here so it's computed FRESH per combination, giving each one
its own physically-appropriate search box rather than one bound shared
across every combination regardless of its own P_amb/f. Pass --pA-bounds
explicitly to override this for all combinations at once instead.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import itertools
import json
import math
import multiprocessing
import os
from pathlib import Path

import global_optimise
import solver_common

REQUIRED_COLUMNS = ("nu_L", "f", "P_amb")

# --- per-worker state (set once by _init_worker when each worker process
# starts, then read back by _run_one_combination for every combination that
# same worker goes on to process) -- see run_batch's own docstring for why
# this is one log file per WORKER rather than one per combination. These
# globals live in the worker process's own memory only: the main process
# never sets or reads them directly, and each worker gets its own separate
# copy (never shared across processes) simply because that's how a fresh
# process's global namespace works.
_worker_id: int | None = None
_worker_log_file = None


def _init_worker(id_counter, id_counter_lock, run_tag: str) -> None:
    """ProcessPoolExecutor's `initializer` -- runs exactly once, right when
    a worker process starts up, before it's given any combination to run.
    id_counter/id_counter_lock are a multiprocessing.Value/Lock shared by
    every worker (picklable specifically so they can cross into each new
    process); using the lock to read-then-increment the counter is what
    guarantees each worker gets its own distinct id (0, 1, 2, ...) even
    though all of them are starting up at roughly the same time. That id
    then picks this worker's own log file, opened once here and kept open
    (module-level, see above) for every combination this worker processes
    for the rest of the pool's life -- write()+flush() per line is enough
    to keep it tailable without needing to reopen it each time."""
    global _worker_id, _worker_log_file
    with id_counter_lock:
        _worker_id = id_counter.value
        id_counter.value += 1
    log_path = solver_common.TMP_DATA_DIR / f"{run_tag}_cpu{_worker_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    # Append mode ("a"), never truncated: reusing a --run-tag (e.g. to
    # resume after a crash) must never discard the log history a previous
    # invocation already wrote -- see RESTART / RESUME in this module's
    # docstring. For a genuinely new --run-tag this file doesn't exist
    # yet, so "a" mode just creates it fresh -- no different from "w" here.
    _worker_log_file = open(log_path, "a")


def _read_csv_rows(csv_path: Path) -> list[tuple[int, dict[str, str]]]:
    """Reads csv_path, skipping blank lines and '#' comment lines, and
    returns (line_number, row_dict) pairs for each DATA row (the header
    itself is consumed, not returned). line_number is the row's actual
    1-based line number in the file -- tracked separately from the
    filtered content -- so an error about row N points at the right place
    in the file even when earlier lines were blank/comments."""
    kept_lines: list[tuple[int, str]] = []
    with open(csv_path, newline="") as f:
        for line_number, line in enumerate(f, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            kept_lines.append((line_number, line))

    if not kept_lines:
        raise ValueError(f"{csv_path}: no header/data found (file is empty, "
                          f"or every line is blank or a '#' comment)")

    header_line_number, header_text = kept_lines[0]
    data_lines = kept_lines[1:]
    reader = csv.DictReader([header_text] + [text for _, text in data_lines])

    missing = [col for col in REQUIRED_COLUMNS if col not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"{csv_path}:{header_line_number}: missing required column(s) "
                          f"{missing} -- found columns {reader.fieldnames}")

    data_line_numbers = [line_number for line_number, _ in data_lines]
    rows = list(reader)
    if len(rows) != len(data_line_numbers):
        # Defensive only -- DictReader should always produce one row per
        # data line here since we fed it exactly that many lines.
        raise ValueError(f"{csv_path}: could not line up parsed rows with source lines")
    return list(zip(data_line_numbers, rows))


def build_cartesian_combinations(csv_path: Path) -> list[tuple[float, float, float]]:
    """--combo-mode cartesian: each column is read as its own independent
    list of values (blank cells simply skipped -- that's how a column
    shorter than the others is represented in a rectangular CSV), then
    combined via a full cross product, same as the old
    --nuL/--f/--Pamb comma-separated lists."""
    columns: dict[str, list[float]] = {col: [] for col in REQUIRED_COLUMNS}
    for line_number, row in _read_csv_rows(csv_path):
        for col in REQUIRED_COLUMNS:
            cell = (row.get(col) or "").strip()
            if not cell:
                continue  # shorter column than the others -- fine in this mode
            try:
                columns[col].append(float(cell))
            except ValueError:
                raise ValueError(f"{csv_path}:{line_number}: column {col!r} has a "
                                  f"non-numeric value: {cell!r}") from None

    empty_columns = [col for col, values in columns.items() if not values]
    if empty_columns:
        raise ValueError(f"{csv_path}: column(s) {empty_columns} have no values at all")

    return list(itertools.product(columns["nu_L"], columns["f"], columns["P_amb"]))


def build_raw_combinations(csv_path: Path) -> list[tuple[float, float, float]]:
    """--combo-mode raw: every row is used exactly as an already-paired
    (nu_L, f, P_amb) triple, no cross product -- so unlike cartesian mode,
    a blank cell here is an error, not "this column is shorter"."""
    combinations: list[tuple[float, float, float]] = []
    seen: set[tuple[float, float, float]] = set()
    for line_number, row in _read_csv_rows(csv_path):
        values = {}
        for col in REQUIRED_COLUMNS:
            cell = (row.get(col) or "").strip()
            if not cell:
                raise ValueError(
                    f"{csv_path}:{line_number}: column {col!r} is blank -- every row needs "
                    f"all three values in --combo-mode raw (use --combo-mode cartesian if "
                    f"the columns are meant to have different lengths)"
                )
            try:
                values[col] = float(cell)
            except ValueError:
                raise ValueError(f"{csv_path}:{line_number}: column {col!r} has a "
                                  f"non-numeric value: {cell!r}") from None

        combo = (values["nu_L"], values["f"], values["P_amb"])
        if combo in seen:
            print(f"Note: {csv_path}:{line_number}: duplicate combination {combo} "
                  f"(already appeared earlier in the file) -- running it again anyway.")
        seen.add(combo)
        combinations.append(combo)

    if not combinations:
        raise ValueError(f"{csv_path}: no data rows found")
    return combinations


def format_log_line(combo_tag: str, nu_L: float, f_khz: float, P_amb_bar: float,
                     ok: bool, R_E_opt=None, p_A_opt=None, energy_demand_opt=None,
                     reason: str = "", boundary_warning: str | None = None) -> str:
    """Deliberately does NOT include opt_path: the combo's tag alone
    ([combo_tag] at the start of the line) already determines the full
    result folder -- tmp_data/<combo_tag>/<combo_tag>_opt.json -- since
    that's exactly how single_run.py/global_optimise.py name it. Spelling
    that path out on every line would just make each one longer without
    adding information."""
    status = "OK" if ok else "FAILED"
    if ok:
        result_part = f"R_E={R_E_opt:.4g} p_A={p_A_opt:.4g} energy_demand={energy_demand_opt:.4g}"
    else:
        result_part = "R_E=- p_A=- energy_demand=-"
    boundary_part = boundary_warning if boundary_warning else "none"
    return (f"[{combo_tag}] nu_L={nu_L} f={f_khz} P_amb={P_amb_bar} -> {result_part} "
            f"| status={status} | reason={reason} | boundary_warning={boundary_part}")


class ComboMismatchError(RuntimeError):
    """Raised by _check_existing_result when an existing
    tmp_data/<combo_tag>/<combo_tag>_opt.json's own recorded (nu_L, f,
    P_amb) doesn't match what this combination's index computes from the
    CURRENT --params-csv -- almost always means the CSV changed since
    whatever earlier invocation produced that file. There is no manifest
    recording what the combinations list looked like back then, so once
    one index is caught disagreeing like this, there's no way to know how
    many OTHER indices are also now silently misaligned (e.g. a value
    inserted into one column shifts the cross product for every index
    after it). Left uncaught by _run_one_combination's own try/except (see
    there) so it propagates all the way out to run_batch's completion
    loop, which stops the batch instead of quietly reusing or overwriting
    a result under the wrong parameters."""


def _read_opt_json(opt_path: Path) -> dict:
    """Parses a <combo_tag>_opt.json the same way archive.py's own
    parse_single_result() does: header JSON, optionally followed by a
    <BINARY> marker and trajectory bytes this function doesn't need."""
    with open(opt_path, "rb") as f:
        content = f.read()
    idx = content.find(b"<BINARY>")
    header = content[:idx] if idx != -1 else content
    return json.loads(header.decode("utf-8").strip())


def _check_existing_result(index: int, combo_dir: Path, combo_tag: str,
                            nu_L: float, f_khz: float, P_amb_bar: float) -> dict | None:
    """Resume support -- see RESTART / RESUME in this module's docstring
    for the full picture. Returns a dict of the already-known result
    (R_E_opt, p_A_opt, energy_demand_opt) if combo_dir/<combo_tag>_opt.json
    already exists, parses cleanly, AND its own recorded parameters agree
    with (nu_L, f_khz, P_amb_bar) -- meaning the caller can skip running
    the optimization entirely. Returns None if no usable result exists yet
    (file missing, or present but unparseable/incomplete -- e.g. left
    truncated by a crash mid-write): both cases mean "not done", safe to
    (re)run fresh. Raises ComboMismatchError if the file parses fine but
    disagrees on the parameters -- see that class's docstring; this is
    deliberately NOT folded into the "just redo it" case."""
    opt_path = combo_dir / f"{combo_tag}_opt.json"
    if not opt_path.exists():
        return None

    try:
        data = _read_opt_json(opt_path)
        cpar = data["cpar"]
        postproc = data["postproc"]
        # Same fields/units single_run.py's build_cpar_config() writes
        # them in (R_E: m, P_amb: Pa, excitation_params: [p_A Pa (negative
        # sign applied), f Hz], mu_L/rho_L: SI) -- converted back to this
        # project's plot units the same way archive.py's build_param_columns
        # does, since we're comparing against nu_L/f_khz/P_amb_bar as given.
        R_E_opt = cpar["R_E"] * 1.0e6                        # m -> um
        p_A_opt = abs(cpar["excitation_params"][0]) / 1.0e5  # Pa -> bar
        f_khz_json = cpar["excitation_params"][1] / 1.0e3    # Hz -> kHz
        P_amb_bar_json = cpar["P_amb"] / 1.0e5                # Pa -> bar
        nu_L_json = cpar["mu_L"] / (cpar["rho_L"] * 1.0e-6)   # -> cSt
        energy_demand_opt = postproc["energy_demand"]
    except (OSError, ValueError, KeyError, IndexError, TypeError, ZeroDivisionError):
        # Corrupt/incomplete -- most likely this combo's confirmation-run
        # file was still being written (e.g. mid-shutil.copy2) when the
        # outage hit. Not trustworthy either way -- treat like "no file".
        return None

    # A SAMENESS check, not a physical-proximity one -- contrast with
    # global_optimise.py's boundary_warning, which uses rel_tol=0.02 to
    # ask "is this near a bound". Here we're asking "is this the same
    # number that merely round-tripped through unit conversions + JSON
    # (de)serialization" -- float64 noise from that is ~1e-15, so 1e-6
    # comfortably absorbs it while still catching any real CSV edit, even
    # a small one (e.g. 100.0 -> 100.001 is a 1e-5 relative change).
    REL_TOL = 1.0e-6
    mismatches = []
    if not math.isclose(nu_L, nu_L_json, rel_tol=REL_TOL):
        mismatches.append(f"nu_L: asked {nu_L} cSt, file has {nu_L_json:.6g} cSt")
    if not math.isclose(f_khz, f_khz_json, rel_tol=REL_TOL):
        mismatches.append(f"f: asked {f_khz} kHz, file has {f_khz_json:.6g} kHz")
    if not math.isclose(P_amb_bar, P_amb_bar_json, rel_tol=REL_TOL):
        mismatches.append(f"P_amb: asked {P_amb_bar} bar, file has {P_amb_bar_json:.6g} bar")
    if mismatches:
        raise ComboMismatchError(
            f"[{combo_tag}] (index {index}) already has a result at {opt_path}, but its "
            f"recorded parameters don't match what --params-csv computes for this index "
            f"now: {'; '.join(mismatches)}. This almost always means params-csv changed "
            f"since the run that produced this file -- with no manifest, that silently "
            f"breaks the index-to-combination mapping for every later index too, not just "
            f"this one. Stopping rather than risk reusing or overwriting a result under "
            f"the wrong parameters. If params-csv really is unchanged, delete or move "
            f"{opt_path} and rerun; if it did change, use a new --run-tag with the new "
            f"CSV instead of resuming this one."
        )

    return {"R_E_opt": R_E_opt, "p_A_opt": p_A_opt, "energy_demand_opt": energy_demand_opt}


def _run_one_combination(index: int, total: int, run_tag: str,
                          nu_L: float, f_khz: float, P_amb_bar: float,
                          pA_bounds: tuple[float, float] | None,
                          RE_bounds: tuple[float, float],
                          tmax: float, timeout: float, maxiter: int, popsize: int,
                          tol: float, atol: float, seed: int | None,
                          polish: bool, cleanup: bool, force: bool) -> dict:
    """The unit of work dispatched to one worker process by run_batch() --
    everything from here down (the solver subprocess calls inside
    run_global_optimization, and this combo's own tmp_data/<combo_tag>/
    result folder) belongs to this ONE combination alone, so nothing here
    needs any coordination with whatever other combinations are running
    concurrently in other worker processes. Every argument is a plain,
    picklable value (no open file handles, no shared objects) since this
    runs in its own separate OS process -- see run_batch. The one thing
    that ISN'T private to this combination is the log line at the very
    end: it goes to _worker_log_file, this worker's own shared, already-
    open file (set up once by _init_worker) -- safe with no locking
    because a given worker only ever runs one combination at a time, so
    there's never a second write racing this one.

    Returns {"row": ..., "summary_line": ...} rather than printing or
    raising: the MAIN process is the only thing that ever calls print()
    (see run_batch), so worker processes can never interleave/garble each
    other's console output; and a combination's own failure is caught and
    reported here rather than propagated, so one bad combination can't take
    down the whole pool or the combinations still running in other workers."""
    combo_tag = f"{run_tag}{index}"
    combo_dir = solver_common.TMP_DATA_DIR / combo_tag
    combo_dir.mkdir(parents=True, exist_ok=True)

    row = {"nu_L": nu_L, "f": f_khz, "P_amb": P_amb_bar, "run_tag": combo_tag,
           "status": "", "R_E_opt": "", "p_A_opt": "", "energy_demand_opt": "",
           "opt_path": "", "log_path": str(_worker_log_file.name), "resumed": False}

    if not force:
        # ComboMismatchError deliberately propagates out of this function
        # uncaught (NOT swallowed by the try/except below) -- see its own
        # docstring and run_batch's handling of it.
        existing = _check_existing_result(index, combo_dir, combo_tag, nu_L, f_khz, P_amb_bar)
        if existing is not None:
            opt_path = combo_dir / f"{combo_tag}_opt.json"
            row.update(status="ok", resumed=True, opt_path=str(opt_path), **existing)
            summary_line = (
                f"[{index+1}/{total}] nu_L={nu_L} f={f_khz} P_amb={P_amb_bar} , "
                f"tag={combo_tag}, R_E={existing['R_E_opt']:.4g} um, "
                f"p_A={existing['p_A_opt']:.4g} bar -> "
                f"energy_demand={existing['energy_demand_opt']:.4g} GJ/t "
                f"| RESUMED (already complete from a previous run)"
            )
            # Per design: a resumed/skipped combo does NOT get a new log
            # line written -- its outcome is already durably recorded in
            # whichever worker's log wrote it the first time this
            # combination actually ran.
            return {"row": row, "summary_line": summary_line}

    try:
        # quiet=True: run_global_optimization() (and every DE evaluation
        # inside it) prints NOTHING -- not even an in-place-updated line --
        # while this combination runs; see its own docstring. This worker
        # doesn't print the resulting summary line itself either -- it's
        # returned to the main process, which is the only thing that prints
        # (see run_batch) -- so nothing here ever touches shared stdout.
        opt_path, (best_RE, best_pA, best_value), diagnostics = global_optimise.run_global_optimization(
            run_tag=combo_tag, nu_L=nu_L, f_khz=f_khz, P_amb_bar=P_amb_bar,
            pA_bounds=pA_bounds, RE_bounds=RE_bounds, tmax=tmax, timeout=timeout,
            maxiter=maxiter, popsize=popsize, tol=tol, atol=atol, seed=seed,
            polish=polish, cleanup=cleanup, quiet=True,
        )
        row.update(status="ok", R_E_opt=best_RE, p_A_opt=best_pA,
                   energy_demand_opt=best_value, opt_path=str(opt_path))
        log_line = format_log_line(
            combo_tag, nu_L, f_khz, P_amb_bar, ok=True,
            R_E_opt=best_RE, p_A_opt=best_pA, energy_demand_opt=best_value,
            reason=diagnostics["message"], boundary_warning=diagnostics["boundary_warning"],
        )

        # DE always returns its best-so-far (R_E, p_A, energy_demand) even
        # when it didn't statistically converge, so those are shown either
        # way; a non-convergent run just gets its reason appended rather
        # than being hidden.
        summary_line = (f"[{index+1}/{total}] nu_L={nu_L} f={f_khz} P_amb={P_amb_bar} , "
                         f"tag={combo_tag}, total eval {diagnostics['nfev']}, "
                         f"R_E={best_RE:.4g} um, p_A={best_pA:.4g} bar "
                         f"-> energy_demand={best_value:.4g} GJ/t")
        if not diagnostics["success"]:
            summary_line += f" | NOT CONVERGED: {diagnostics['message']}"
    except Exception as e:
        # No best (R_E, p_A, energy_demand) exists in this case -- the
        # exception means run_global_optimization() never got that far --
        # so this line/row can only report the combination and the error.
        row["status"] = f"failed: {e}"
        log_line = format_log_line(combo_tag, nu_L, f_khz, P_amb_bar, ok=False, reason=str(e))
        summary_line = (f"[{index+1}/{total}] nu_L={nu_L} f={f_khz} P_amb={P_amb_bar} , "
                         f"tag={combo_tag} -> FAILED: {e}")

    _worker_log_file.write(log_line + "\n")
    _worker_log_file.flush()  # so this worker's file is readable while the batch is still running

    return {"row": row, "summary_line": summary_line}


def run_batch(run_tag: str, combinations: list[tuple[float, float, float]],
              pA_bounds: tuple[float, float] | None = None,
              RE_bounds: tuple[float, float] = global_optimise.DEFAULT_RE_BOUNDS,
              tmax: float = 1.0, timeout: float = 60.0,
              maxiter: int = global_optimise.DEFAULT_MAXITER,
              popsize: int = global_optimise.DEFAULT_POPSIZE,
              tol: float = 0.01, atol: float = 0.0, seed: int | None = None,
              polish: bool = True, cleanup: bool = True,
              workers: int | None = None, force: bool = False) -> list[dict]:
    """combinations: an already-built list of (nu_L, f_khz, P_amb_bar)
    triples -- see build_cartesian_combinations/build_raw_combinations,
    which is how main() builds this from --params-csv/--combo-mode. This
    function itself doesn't care how the list was produced.

    workers: how many combinations run at once, each in its own OS process
    (default: os.cpu_count()). Pass 1 to fall back to one-at-a-time
    execution (still via the same process-pool machinery, just a pool of
    size 1) if you ever want that -- e.g. to watch a single combination's
    solver output live by temporarily setting global_optimise.py's
    quiet=False yourself, which isn't meaningful to do across several
    concurrent workers at once.

    On Windows (and anywhere else using the 'spawn' start method), calling
    this from your own script -- rather than through this file's own
    --run-tag/--params-csv CLI, which already does this -- needs to happen
    inside an `if __name__ == "__main__":` guard, same as any other use of
    multiprocessing/ProcessPoolExecutor; otherwise each worker process
    re-imports and re-runs your top-level script code too.

    Returns a list of one dict per combination (nu_L, f, P_amb, run_tag,
    status, R_E_opt, p_A_opt, energy_demand_opt, opt_path, log_path), in
    ORIGINAL combination order regardless of the order they actually
    finished in -- for programmatic use. The only files this function's
    workers write are the --workers worth of tmp_data/<run_tag>_cpu<id>.log
    files (see _init_worker/_run_one_combination); nothing is written to a
    single shared log or to a CSV."""
    if workers is None:
        workers = os.cpu_count() or 1
    total = len(combinations)
    print(f"Scanning {total} combination(s) of (nu_L, f, P_amb) across {workers} worker process(es)...")

    # Shared (picklable) counter + lock, handed to every worker via
    # _init_worker so each one can atomically claim its own distinct id --
    # see _init_worker's own docstring for why this, rather than e.g.
    # os.getpid(), is what gives clean 0..workers-1 log file names.
    id_counter = multiprocessing.Value("i", 0)
    id_counter_lock = multiprocessing.Lock()

    results: list[dict | None] = [None] * total
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=workers, initializer=_init_worker,
        initargs=(id_counter, id_counter_lock, run_tag),
    ) as executor:
        future_to_index = {
            executor.submit(
                _run_one_combination, i, total, run_tag, nu_L, f_khz, P_amb_bar,
                pA_bounds, RE_bounds, tmax, timeout, maxiter, popsize, tol, atol, seed,
                polish, cleanup, force,
            ): i
            for i, (nu_L, f_khz, P_amb_bar) in enumerate(combinations)
        }
        for future in concurrent.futures.as_completed(future_to_index):
            i = future_to_index[future]
            nu_L, f_khz, P_amb_bar = combinations[i]
            combo_tag = f"{run_tag}{i}"
            try:
                outcome = future.result()
            except ComboMismatchError as e:
                # See ComboMismatchError's own docstring: this means the
                # CSV likely changed since an earlier run of this run_tag,
                # so index-based resume can no longer be trusted for this
                # OR any later index -- stop the whole batch rather than
                # keep going and risk silently reusing/overwriting results
                # under the wrong parameters. Combinations already running
                # in other worker processes at this moment finish naturally
                # (can't be interrupted mid-run); cancel_futures=True just
                # stops anything not yet started from being dispatched.
                print(f"\nFATAL: {e}\n")
                executor.shutdown(cancel_futures=True)
                raise
            except Exception as e:
                # Only reached if the worker PROCESS itself died/couldn't
                # even run _run_one_combination's own try/except (e.g. it
                # crashed or was killed) -- a normal combination-level
                # failure is already caught inside _run_one_combination and
                # comes back as a "failed: ..." row instead of raising here.
                print(f"[{i+1}/{total}] nu_L={nu_L} f={f_khz} P_amb={P_amb_bar} , "
                      f"tag={combo_tag} -> FAILED: worker process error: {e}")
                results[i] = {"nu_L": nu_L, "f": f_khz, "P_amb": P_amb_bar, "run_tag": combo_tag,
                              "status": f"failed: worker process error: {e}", "R_E_opt": "",
                              "p_A_opt": "", "energy_demand_opt": "", "opt_path": "", "log_path": "",
                              "resumed": False}
                continue
            print(outcome["summary_line"])
            results[i] = outcome["row"]

    n_ok = sum(1 for r in results if r["status"] == "ok")
    n_resumed = sum(1 for r in results if r.get("resumed"))
    print(f"\nBatch complete: {n_ok}/{len(results)} succeeded ({n_resumed} reused from a "
          f"previous run). Per-worker logs: tmp_data/{run_tag}_cpu0.log .. "
          f"tmp_data/{run_tag}_cpu{workers - 1}.log")
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Optimize p_A/R_E for each (nu_L, f, P_amb) combination read from a CSV file."
    )
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--params-csv", required=True,
                         help="CSV filename with a 'nu_L,f,P_amb' header, resolved inside "
                              "tmp_data/ (e.g. --params-csv combos.csv means tmp_data/combos.csv)")
    parser.add_argument("--combo-mode", choices=["cartesian", "raw"], default="cartesian",
                         help="'cartesian' (default): each column is its own list of values, "
                              "combined via a full cross product -- columns may have different "
                              "lengths (pad the shorter ones with blank cells). "
                              "'raw': every row is used exactly as an already-paired "
                              "(nu_L, f, P_amb) triple, no cross product -- every row needs all "
                              "three values.")
    parser.add_argument("--pA-bounds", default=None,
                         help="min:max, [bar]. Default: computed per-combination as P_amb + 0.5*f")
    parser.add_argument("--RE-bounds", default="1.0:1000.0", help="min:max, [um]")
    parser.add_argument("--tmax", type=float, default=1.0)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--maxiter", type=int, default=global_optimise.DEFAULT_MAXITER)
    parser.add_argument("--popsize", type=int, default=global_optimise.DEFAULT_POPSIZE)
    parser.add_argument("--tol", type=float, default=0.01)
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=None, help="for reproducible runs")
    parser.add_argument("--no-polish", action="store_true")
    parser.add_argument("--keep-history", action="store_true",
                         help="keep all numbered intermediate evaluation files for every combination")
    parser.add_argument("--workers", type=int, default=None,
                         help="number of combinations to run at once, each in its own process "
                              "(default: os.cpu_count(), i.e. all logical cores)")
    parser.add_argument("--force", action="store_true",
                         help="ignore any existing tmp_data/<combo_tag>/<combo_tag>_opt.json "
                              "files and recompute every combination from scratch, even ones "
                              "that already completed in an earlier invocation of this same "
                              "--run-tag (see RESTART / RESUME above)")
    args = parser.parse_args()

    csv_path = solver_common.TMP_DATA_DIR / args.params_csv
    if args.combo_mode == "cartesian":
        combinations = build_cartesian_combinations(csv_path)
    else:
        combinations = build_raw_combinations(csv_path)

    run_batch(
        run_tag=args.run_tag,
        combinations=combinations,
        pA_bounds=global_optimise.parse_bounds(args.pA_bounds) if args.pA_bounds is not None else None,
        RE_bounds=global_optimise.parse_bounds(args.RE_bounds),
        workers=args.workers,
        tmax=args.tmax, timeout=args.timeout, maxiter=args.maxiter, popsize=args.popsize,
        tol=args.tol, atol=args.atol, seed=args.seed, polish=not args.no_polish,
        cleanup=not args.keep_history, force=args.force,
    )


if __name__ == "__main__":
    main()
