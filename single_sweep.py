"""
single_sweep.py

Builds a parameter_study JSON config from fixed background conditions
(nu_L, f, P_amb) plus swept p_A and R_E ranges, and runs it through the
solver. The baseline template (config_data/example_parameter_study.json)
is loaded read-only; only the fields below are overridden, everything
else (species, fractions, mechanism, ...) comes from the baseline as-is.

Requires solver_common.py in the same folder for shared, non-sweep-
specific utilities (liquid-property physics, JSON-comment stripping,
executable path resolution) -- this file only contains parameter_study/
sweep-specific logic (ranges, scales, resolutions) on top of that.

Usage:
    python single_sweep.py --run-tag=nh3_sweep_01 \
        --nuL=0.65 --f=10 --Pamb=1.0 \
        --pA="1.0:20.0:100:lin" --RE="1.0:1000.0:100:log"

pA is given as a POSITIVE magnitude on the command line — the script applies
the conventional negative sign (rarefaction-first wave) internally when
building the JSON, so you never need to type a leading '-' for it.

Range spec format "min:max:res:scale":
    min, max  - bounds, in the units noted below
    res       - number of steps (integer)
    scale     - "lin" (LinearRange) or "log" (LogRange)

Units taken on the command line (human/plot units, converted internally
to the solver's raw SI units):
    nu_L  [cSt]     f     [kHz]    P_amb [bar]
    p_A   [bar]     R_E   [um]

Optional: --tmax (default 1.0 s), --timeout (default 60.0 s),
--cpu (default: all cores), --dry-run (build the config and print the
solver command without actually running it).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import solver_common
from solver_common import SOLVER_DIR, TMP_DATA_DIR

BASELINE_PATH = Path("config_data") / "example_parameter_study.json"

SCALE_TO_RANGE_TYPE = {"lin": "LinearRange", "log": "LogRange"}


def parse_range_spec(spec: str, require_positive: bool = False) -> tuple[float, float, int, str]:
    """Parse 'min:max:res:scale' -> (min, max, res, scale)."""
    parts = spec.split(":")
    if len(parts) != 4:
        raise ValueError(f"range spec must be 'min:max:res:scale', got: {spec!r}")
    min_val, max_val, res, scale = parts
    scale = scale.strip().lower()
    if scale not in SCALE_TO_RANGE_TYPE:
        raise ValueError(f"scale must be 'lin' or 'log', got: {scale!r}")
    min_val, max_val = float(min_val), float(max_val)
    if require_positive and (min_val < 0 or max_val < 0):
        raise ValueError(f"expected positive magnitudes (sign is applied internally), got: {spec!r}")
    return min_val, max_val, int(res), scale


def build_range_field(min_val: float, max_val: float, res: int, scale: str) -> dict:
    return {"type": SCALE_TO_RANGE_TYPE[scale], "start": min_val, "end": max_val, "num_steps": res}


def load_baseline() -> dict:
    with open(BASELINE_PATH) as f:
        raw_text = f.read()
    return json.loads(solver_common.strip_json_comments(raw_text))["parameter_study"]


def build_config(nu_L: float, f_khz: float, P_amb_bar: float,
                  pA_spec: str, RE_spec: str) -> dict:
    config = load_baseline()

    pA_min, pA_max, pA_res, pA_scale = parse_range_spec(pA_spec, require_positive=True)
    RE_min, RE_max, RE_res, RE_scale = parse_range_spec(RE_spec, require_positive=True)

    config["R_E"] = build_range_field(RE_min * 1.0e-6, RE_max * 1.0e-6, RE_res, RE_scale)  # um -> m
    config["P_amb"] = solver_common.const_field(P_amb_bar * 1.0e5)                          # bar -> Pa
    config["excitation_params"] = [
        # pA entered as a positive magnitude; negated here to match this
        # project's rarefaction-first (negative) excitation convention
        build_range_field(-pA_min * 1.0e5, -pA_max * 1.0e5, pA_res, pA_scale),              # bar -> Pa
        solver_common.const_field(f_khz * 1.0e3),                                           # kHz -> Hz
    ]
    config["target_specie"] = "NH3"
    config.update(solver_common.compute_liquid_properties(nu_L))

    return config


def run_sweep(run_tag: str, nu_L: float, f_khz: float, P_amb_bar: float,
              pA_spec: str, RE_spec: str, tmax: float = 1.0, timeout: float = 60.0,
              cpu: int | None = None, dry_run: bool = False) -> Path | None:
    """Returns the resulting output folder (Path), or None if this was a
    dry run or the solver produced no output at all."""
    TMP_DATA_DIR.mkdir(exist_ok=True)

    config = build_config(nu_L, f_khz, P_amb_bar, pA_spec, RE_spec)
    config_path = TMP_DATA_DIR / f"{run_tag}.json"
    with open(config_path, "w") as f:
        json.dump({"parameter_study": config}, f, indent=2)
    print(f"Wrote config: {config_path}")

    cpu = cpu if cpu is not None else os.cpu_count()
    cmd = [
        solver_common.solver_executable(SOLVER_DIR),
        "--parameter_study", f"../tmp_data/{run_tag}.json",
        "--tmax", str(tmax),
        "--timeout", str(timeout),
        "--directory", f"../tmp_data/{run_tag}",
        "--cpu", str(cpu),
    ]
    print("Solver command:", " ".join(cmd), f"(cwd={SOLVER_DIR})")

    if dry_run:
        print("Dry run — solver not executed.")
        return None

    result = subprocess.run(cmd, cwd=SOLVER_DIR)

    matches = sorted((p for p in TMP_DATA_DIR.glob(f"{run_tag}*") if p.is_dir()), key=os.path.getmtime)
    output_exists = bool(matches) and any(matches[-1].glob("output_*.csv"))

    if output_exists:
        # The solver's own bruteforce_parameter_study_settings.json inside
        # the run folder (written from its own parsed/validated
        # ParameterCombinator) is a superset of this root-level input copy
        # -- same fields, plus an "info" section -- so once we know that
        # file exists, config_path is just clutter accumulating in
        # tmp_data's root. Gated on output_exists rather than deleted
        # unconditionally: a fatal preprocessing error (e.g. malformed
        # JSON) prevents ParameterStudy -- and therefore that settings
        # file -- from ever being created, so this input copy stays behind
        # as the only record of what was attempted in that case. A
        # --dry-run never reaches this line at all, so its config_path is
        # always left in place for inspection.
        config_path.unlink(missing_ok=True)

    if result.returncode != 0:
        if output_exists:
            print(f"Note: solver exited with code {result.returncode} — this happens whenever at "
                  f"least one individual simulation in the sweep failed (see errors.log in the "
                  f"output folder). The sweep itself completed; see the summary above.")
        else:
            print(f"Solver exited with code {result.returncode} and no output was found — "
                  f"this looks like a genuine failure, not just some failed combinations. "
                  f"Check the solver's own console output above for the actual error.")
            return None

    if matches:
        print(f"Output folder: {matches[-1]}")
        return matches[-1]
    return None


def main():
    parser = argparse.ArgumentParser(description="Run a pA vs RE parameter sweep.")
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--nuL", type=float, required=True, help="[cSt]")
    parser.add_argument("--f", type=float, required=True, help="[kHz]")
    parser.add_argument("--Pamb", type=float, required=True, help="[bar]")
    parser.add_argument("--pA", required=True, help="min:max:res:scale, positive [bar] magnitude")
    parser.add_argument("--RE", required=True, help="min:max:res:scale, [um]")
    parser.add_argument("--tmax", type=float, default=1.0)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--cpu", type=int, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    run_sweep(
        run_tag=args.run_tag, nu_L=args.nuL, f_khz=args.f, P_amb_bar=args.Pamb,
        pA_spec=args.pA, RE_spec=args.RE, tmax=args.tmax, timeout=args.timeout,
        cpu=args.cpu, dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
