"""
global_optimise.py

Finds the global (R_E, p_A) optimum minimizing energy_demand using
scipy.optimize.differential_evolution.

Renamed from direct_optimise.py / switched from scipy.optimize.direct after
real-data testing showed DIRECT converges very slowly on this project's
actual landscape shape: a thin, curved "ridge" of good performance (visible
directly in the reference contour plots) surrounded by a large, flat,
uniformly bad region. DIRECT explores by recursively splitting the search
box into axis-aligned rectangles — a strategy that struggles badly with a
diagonal/curved feature, since no axis-aligned split ever captures much of
a diagonal line without many rounds of refinement (a well-known weak point
of DIRECT-family methods, the same reason Rosenbrock-shaped valleys are a
standard optimizer benchmark). Differential evolution instead maintains a
population of candidates and moves them via vector differences between
population members, which can travel diagonally through the space directly
— much better suited to curved-valley landscapes like this one.

Requires solver_common.py and single_run.py in the same folder — each DE
evaluation is one single_run.run_single() call (the solver's lighter --run
mode), not a grid sweep, since DE (like DIRECT) probes arbitrary,
non-grid-aligned points.

Usage:
    python global_optimise.py --run-tag nh3_global --nuL 0.65 --f 10 --Pamb 1.0

Every evaluation tried is kept as its own file:
tmp_data/<run_tag>/<run_tag><id>.json (id = 0, 1, 2, ... via
single_run.py's own auto-incrementing). Once DE converges, ONE final
single_run() is made at the exact optimum found (this time with --save,
for the full trajectory) and copied to
tmp_data/<run_tag>/<run_tag>_opt.json.

R_E is searched in log10-space internally (bounds given in real units on
the command line), same reasoning as before: it spans orders of magnitude,
so linear search would waste resolution on the upper end of the range.

Failed simulations: single_run.read_energy_demand() already maps both
solver failure modes to usable floats; anything above FAILURE_PENALTY here
is additionally clipped down to it, since DE's internal arithmetic (like
DIRECT's) can misbehave with values near float64's max.

DE-specific choices:
- polish=True (default here): after DE's population-based search finishes,
  scipy runs one local gradient-based (L-BFGS-B) refinement step from the
  best population member, using finite-difference gradients. This can
  meaningfully sharpen the final answer, but it evaluates the objective a
  handful of extra times (each a real solver call) to estimate those
  gradients, and since our objective is a genuine simulation (not a smooth
  analytic function), any run-to-run numerical noise in the solver could
  make those finite-difference gradient estimates unreliable in principle.
  In practice, for a smoothly-varying physical quantity like energy_demand
  this is usually still a net win -- just worth knowing it's not free, and
  not guaranteed to help if the local landscape is genuinely noisy right at
  the optimum. Disable via polish=False if you'd rather skip it.
- workers=1 (sequential) is required, not just a default: single_run.py's
  auto-incrementing filename scheme assumes calls happen one at a time
  (see its own docstring) — parallel workers sharing one run_tag could
  race and overwrite each other's files.

CONVERGENCE -- exactly when DE stops (from scipy's own source):
1. Statistical convergence (the "normal" stop): every generation, DE checks
       std(population_energies) <= atol + tol * abs(mean(population_energies))
   i.e. once the population's spread of energy_demand values collapses to
   within a band set by --tol (relative) and --atol (absolute) around the
   population mean, it's considered converged and stops -- this is the
   precision control you're looking for. Both default small (tol=0.01,
   atol=0), so a genuinely noisy/flat-plateau-heavy objective (like this
   project's real landscape) can fail to satisfy this for a very long time,
   which is almost certainly why maxiter was hit instead in your run.
2. maxiter reached: if statistical convergence (1) is never satisfied, DE
   just stops after maxiter generations regardless, returning its current
   best. This is what "maxed out" means -- not a sign of failure per se,
   but a sign the precision target in (1) was never reached in the
   generation budget given. result.message / printed output distinguishes
   these two cases explicitly now (see run_global_optimization).
3. (Not currently used here, mentioned for completeness) a user-supplied
   callback returning True can also stop DE early on custom logic.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
from pathlib import Path

from scipy.optimize import differential_evolution

import single_run
import solver_common

FAILURE_PENALTY = 1.0e12  # finite stand-in for any failed/huge energy_demand, safe for DE's internals

DEFAULT_PA_LOWER = 1.0              # [bar], unchanged -- never search below this
DEFAULT_RE_BOUNDS = (1.0, 1000.0)   # [um]
DEFAULT_MAXITER = 40                # generations; actual evals ~= popsize * ndim * generations used
DEFAULT_POPSIZE = 15                # scipy's own default; population = popsize * ndim (ndim=2 here)


def default_pA_upper_bound(P_amb_bar: float, f_khz: float) -> float:
    """Physically-motivated default upper bound: P_amb [bar] + 0.5 * f [kHz],
    used whenever --pA-bounds isn't given explicitly."""
    return P_amb_bar + 0.5 * f_khz


def make_objective(run_tag: str, nu_L: float, f_khz: float, P_amb_bar: float,
                    tmax: float, timeout: float, quiet: bool = False):
    """Returns a function(x) -> energy_demand, x = [p_A_bar, log10(R_E_um)],
    that DIRECT can call directly. Each call is one real single_run().

    quiet=False (default, standalone CLI use of this script): unchanged
    behaviour -- one full scrolling "[eval N] ..." line per eval, exactly
    as before.

    quiet=True (used by batch_global_optimise.py, one DE run per
    Cartesian-product combination): prints NOTHING per eval, not even an
    in-place-updated line -- batch_global_optimise.py prints exactly one
    summary line per combination once run_global_optimization() returns
    instead. Two reasons for going all the way to silent rather than an
    in-place-updated line: chemistry/CVODE warnings are frequent by
    nature, so per-eval detail is mostly noise across a long batch; and a
    terminal has no notion of "my line" vs "someone else's line" -- an
    in-place update from N combinations running concurrently (the
    intended next step, see run_global_optimization) would garble each
    other, whereas printing nothing at all per eval is trivially safe
    under any amount of parallelism."""
    eval_count = 0

    def objective(x):
        nonlocal eval_count
        eval_count += 1
        pA_bar, log10_RE_um = x
        RE_um = 10.0 ** log10_RE_um

        result = single_run.run_single(
            run_tag=run_tag, nu_L=nu_L, f_khz=f_khz, P_amb_bar=P_amb_bar,
            pA_bar=pA_bar, RE_um=RE_um, tmax=tmax, timeout=timeout, save=False,
            verbose=False,  # stays silent itself -- see quiet below for who
                            # (if anyone) prints anything about this eval
        )
        energy_demand = single_run.read_energy_demand(result.config_path)
        if not math.isfinite(energy_demand) or energy_demand > FAILURE_PENALTY:
            energy_demand = FAILURE_PENALTY

        if not quiet:
            result_part = (f"R_E={RE_um:.4g} um, p_A={pA_bar:.4g} bar "
                           f"-> energy_demand={energy_demand:.4g} GJ/t")
            if result.condensed_note:
                result_part += f" | {result.condensed_note}"
            print(f"[eval {eval_count}] {result_part}")
        return energy_demand

    return objective


def cleanup_exploration_files(run_dir: Path, run_tag: str) -> int:
    """Deletes every numbered '<run_tag><id>.json' exploration file in
    run_dir (the ones single_run.py's auto-incrementing produced along the
    way), keeping only '<run_tag>_opt.json'. Matches the same pattern
    single_run.next_available_id() uses, so it can't accidentally catch
    the _opt file itself (that name has a non-digit suffix). Returns the
    count of files deleted."""
    pattern = re.compile(rf"^{re.escape(run_tag)}(\d+)\.json$")
    count = 0
    for p in run_dir.glob(f"{run_tag}*.json"):
        if pattern.match(p.name):
            p.unlink()
            count += 1
    return count


def run_global_optimization(run_tag: str, nu_L: float, f_khz: float, P_amb_bar: float,
                             pA_bounds: tuple[float, float] | None = None,
                             RE_bounds: tuple[float, float] = DEFAULT_RE_BOUNDS,
                             tmax: float = 1.0, timeout: float = 60.0,
                             maxiter: int = DEFAULT_MAXITER, popsize: int = DEFAULT_POPSIZE,
                             tol: float = 0.01, atol: float = 0.0, seed: int | None = None,
                             polish: bool = True, cleanup: bool = True,
                             quiet: bool = False):
    """quiet=False (default, standalone CLI use): unchanged -- prints its
    usual setup/progress/summary lines, plus one scrolling "[eval N] ..."
    line per DE evaluation via make_objective().

    quiet=True (used by batch_global_optimise.py): completely silent --
    no setup line, no per-eval line, no finished/stopping-reason/boundary
    WARNING/save/cleanup lines. batch_global_optimise.py prints its own
    single summary line per combination from the returned diagnostics
    instead, once this call returns. This isn't just "less verbose": it's
    what makes running several combinations concurrently (the intended
    next step) safe to look at, since nothing here writes to the shared
    terminal until the whole combination is done."""
    if pA_bounds is None:
        pA_upper = default_pA_upper_bound(P_amb_bar, f_khz)
        pA_bounds = (DEFAULT_PA_LOWER, pA_upper)
        if not quiet:
            print(f"p_A upper bound not given explicitly -- using P_amb + 0.5*f = "
                  f"{P_amb_bar} + 0.5*{f_khz} = {pA_upper} bar")

    if pA_bounds[1] <= pA_bounds[0]:
        raise ValueError(
            f"p_A upper bound ({pA_bounds[1]} bar) must be greater than the lower bound "
            f"({pA_bounds[0]} bar) -- check P_amb/f if this came from the default formula, "
            f"or pass --pA-bounds explicitly."
        )

    bounds = [pA_bounds, (math.log10(RE_bounds[0]), math.log10(RE_bounds[1]))]
    objective = make_objective(run_tag, nu_L, f_khz, P_amb_bar, tmax, timeout, quiet=quiet)

    if not quiet:
        print(f"Running differential_evolution: p_A in {pA_bounds} bar, R_E in {RE_bounds} um "
              f"(log10: {bounds[1]}), maxiter={maxiter} generations, popsize={popsize}, "
              f"tol={tol}, atol={atol}, polish={polish}")
    result = differential_evolution(
        objective, bounds=bounds, maxiter=maxiter, popsize=popsize, tol=tol, atol=atol,
        seed=seed, polish=polish, workers=1, updating="immediate",
    )

    best_pA, best_log10_RE = result.x
    best_RE = 10.0 ** best_log10_RE
    if not quiet:
        print(f"\ndifferential_evolution finished: R_E={best_RE:.4g} um, p_A={best_pA:.4g} bar, "
              f"energy_demand={result.fun:.4g} GJ/t ({result.nfev} evaluations)")
        print(f"Stopping reason: {result.message}")
        if not result.success:
            print("NOTE: this means maxiter (generations) was exhausted WITHOUT satisfying the "
                  "statistical convergence criterion (tol/atol) -- the result may not be as precise "
                  "as a converged run. Consider raising --maxiter, loosening --tol/--atol, or both.")

    boundary_warnings = []
    if math.isclose(best_pA, pA_bounds[1], rel_tol=0.02):
        warning = (f"p_A ({best_pA:.4g} bar) at/near upper bound ({pA_bounds[1]:.4g} bar)")
        if not quiet:
            print(f"WARNING: optimum {warning} -- the true optimum may lie beyond it. "
                  f"Consider whether P_amb + 0.5*f is actually a hard physical limit here, "
                  f"or pass a wider --pA-bounds explicitly and rerun.")
        boundary_warnings.append(warning)
    if math.isclose(best_RE, RE_bounds[1], rel_tol=0.02) or math.isclose(best_RE, RE_bounds[0], rel_tol=0.02):
        warning = f"R_E ({best_RE:.4g} um) at/near a bound ({RE_bounds})"
        if not quiet:
            print(f"WARNING: optimum {warning} -- consider widening --RE-bounds and rerunning.")
        boundary_warnings.append(warning)

    diagnostics = {
        "success": result.success,
        "message": result.message,
        "boundary_warning": "; ".join(boundary_warnings) if boundary_warnings else None,
        "nfev": result.nfev,
    }

    # Final confirmation run at the exact optimum, this time with the full
    # trajectory saved, copied to a fixed, easy-to-find filename.
    final_result = single_run.run_single(
        run_tag=run_tag, nu_L=nu_L, f_khz=f_khz, P_amb_bar=P_amb_bar,
        pA_bar=best_pA, RE_um=best_RE, tmax=tmax, timeout=timeout, save=True,
        verbose=False,  # "Final optimum run saved" below already announces this
    )
    opt_path = solver_common.TMP_DATA_DIR / run_tag / f"{run_tag}_opt.json"
    shutil.copy2(final_result.config_path, opt_path)
    if not quiet:
        note = f" | {final_result.condensed_note}" if final_result.condensed_note else ""
        print(f"Final optimum run saved: {opt_path}{note}")

    if cleanup:
        run_dir = opt_path.parent
        n_deleted = cleanup_exploration_files(run_dir, run_tag)
        if not quiet:
            print(f"Cleaned up {n_deleted} intermediate evaluation file(s) -- {run_dir} now "
                  f"contains only {opt_path.name}.")

    return opt_path, (best_RE, best_pA, result.fun), diagnostics


def parse_bounds(spec: str) -> tuple[float, float]:
    low, high = spec.split(":")
    return float(low), float(high)


def main():
    parser = argparse.ArgumentParser(description="Find the global (R_E, p_A) optimum via differential_evolution.")
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--nuL", type=float, required=True, help="[cSt]")
    parser.add_argument("--f", type=float, required=True, help="[kHz]")
    parser.add_argument("--Pamb", type=float, required=True, help="[bar]")
    parser.add_argument("--pA-bounds", default=None,
                         help="min:max, [bar]. Default: 1.0 to (P_amb + 0.5*f)")
    parser.add_argument("--RE-bounds", default="1.0:1000.0", help="min:max, [um]")
    parser.add_argument("--tmax", type=float, default=1.0)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--maxiter", type=int, default=DEFAULT_MAXITER,
                         help="max generations (each generation is ~popsize*2 solver evaluations)")
    parser.add_argument("--popsize", type=int, default=DEFAULT_POPSIZE,
                         help="population size multiplier (actual population = popsize * 2)")
    parser.add_argument("--tol", type=float, default=0.01, help="relative convergence tolerance")
    parser.add_argument("--atol", type=float, default=0.0, help="absolute convergence tolerance")
    parser.add_argument("--seed", type=int, default=None, help="for reproducible runs")
    parser.add_argument("--no-polish", action="store_true",
                         help="skip the final local (L-BFGS-B) refinement step")
    parser.add_argument("--keep-history", action="store_true",
                         help="keep all numbered intermediate evaluation files instead of "
                              "deleting them, leaving only <run_tag>_opt.json")
    args = parser.parse_args()

    run_global_optimization(
        run_tag=args.run_tag, nu_L=args.nuL, f_khz=args.f, P_amb_bar=args.Pamb,
        pA_bounds=parse_bounds(args.pA_bounds) if args.pA_bounds is not None else None,
        RE_bounds=parse_bounds(args.RE_bounds),
        tmax=args.tmax, timeout=args.timeout, maxiter=args.maxiter, popsize=args.popsize,
        tol=args.tol, atol=args.atol, seed=args.seed, polish=not args.no_polish,
        cleanup=not args.keep_history,
    )


if __name__ == "__main__":
    main()
