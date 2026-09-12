"""
solver_common.py

Generic solver-interfacing utilities shared by BOTH the sweep pathway
(single_sweep.py) and the single-run/optimization pathway (single_run.py,
global_optimise.py, batch_global_optimise.py): PDMS liquid-property
physics, JSON-comment stripping, and solver executable/path resolution.

None of this is sweep-specific -- it's needed identically by a single
--run simulation and by a --parameter_study sweep. It used to live inside
single_sweep.py, which meant the single-run/optimization chain had to
import a sweep-specific module just to reach a few generic helpers,
mixing two things that should be independent. This module has no
knowledge of parameter_study/sweep concepts (ranges, scales, resolutions)
at all -- that logic stays in single_sweep.py, which now imports its own
shared pieces from here rather than defining them.
"""

from __future__ import annotations

import platform
from pathlib import Path

TMP_DATA_DIR = Path("tmp_data")
SOLVER_DIR = Path("solver")


def const_field(value: float) -> dict:
    return {"type": "Const", "value": value}


def compute_liquid_properties(nu_L: float) -> dict:
    """PDMS liquid property correlation chain, all derived from nu_L [cSt].
    Returns Const-wrapped fields ({"type": "Const", "value": ...}) -- the
    format the parameter_study/sweep JSON needs directly; callers that
    need plain scalars instead (e.g. --run/cpar mode, see single_run.py)
    unwrap ["value"] themselves rather than this function doing it, since
    which format is needed depends on the JSON mode, not on the physics."""
    rho_L_ref = 972.47 - 157.24 * nu_L ** (-0.6649)          # [kg/m^3]
    mu_L = nu_L * 1.0e-6 * rho_L_ref                          # [Pa*s]
    c_L = 1001.86 - 83.97 * nu_L ** (-0.6281)                 # [m/s]
    surfactant = (21.1750 - 3.9579 * nu_L ** (-0.6163)) / 0.07197 / 1000.0  # [-]
    Gamma_L = 8.4207                                          # [-] constant
    p_L_ref = 100000.0                                        # [Pa] constant
    B_L = rho_L_ref * c_L ** 2 / Gamma_L                      # [Pa]
    return {
        "rho_L": const_field(rho_L_ref),
        "mu_L": const_field(mu_L),
        "c_L": const_field(c_L),
        "surfactant": const_field(surfactant),
        "liquid_eos_params": [
            const_field(Gamma_L), const_field(B_L), const_field(p_L_ref), const_field(rho_L_ref),
        ],
    }


def strip_json_comments(text: str) -> str:
    """Remove '//' line comments from JSON text, respecting string literals
    (so a '//' inside a quoted value is left alone). The baseline template
    files use '//' comments to document the nu_L formulas — valid for the
    solver's own JSON library (nlohmann::json), but not for Python's
    standard json module, which rejects them outright."""
    out = []
    in_string = False
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if in_string:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            if c == '"':
                in_string = False
            i += 1
            continue
        if c == '"':
            in_string = True
            out.append(c)
            i += 1
            continue
        if c == "/" and i + 1 < n and text[i + 1] == "/":
            while i < n and text[i] != "\n":
                i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def solver_executable(solver_dir: Path) -> str:
    """Absolute path to the solver binary. Must be absolute, not relative:
    on Windows, subprocess/CreateProcess resolves a relative executable
    path against the PARENT process's cwd, not the cwd= passed to
    subprocess.run() — so a relative "bin/main.exe" fails to be found
    even though cwd=solver_dir is set correctly for everything else.

    The binary name itself also differs by platform, not just the path
    separator: dev/build_utility.py in the C++ project only appends .exe
    when os.name == 'nt' -- on Linux/Mac the compiled binary is plain
    "main" with no extension. (pathlib's "/" operator parses forward
    slashes correctly on Windows too, so there's no need to special-case
    the separator itself -- only the extension actually differs.)"""
    exe_name = "bin/main.exe" if platform.system() == "Windows" else "bin/main"
    return str((solver_dir / exe_name).resolve())
