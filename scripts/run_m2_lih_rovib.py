#!/usr/bin/env python3
"""Run the OES-M2.1 LiH rovibrational spectrum gate on an OES BO curve."""

import json
import math
from pathlib import Path

import numpy as np

from oes.quantum.diatomic_rovib import (
    fit_rotational_dunham,
    fit_vibrational_dunham,
    solve_diatomic_rovibrational_levels,
)
from oes.quantum.lih_m2 import LIH_REDUCED_NUCLEAR_MASS_ME, run_lih_ground_energy

ROOT = Path(__file__).resolve().parents[1]
TARGET = json.loads((ROOT / "benchmarks" / "lih_m2_nist.json").read_text())
BO_GRID = tuple(2.0 + 0.2 * i for i in range(18))


def main():
    energies = tuple(run_lih_ground_energy(r, basis_name="cc-pVTZ") for r in BO_GRID)
    spectrum = solve_diatomic_rovibrational_levels(
        BO_GRID,
        energies,
        LIH_REDUCED_NUCLEAR_MASS_ME,
        n_vibrational=4,
        j_values=(0, 1, 2),
        radial_grid_points=1800,
    )
    vib = fit_vibrational_dunham(spectrum.term_values_cm_by_j[0])
    rot = fit_rotational_dunham(spectrum.term_values_cm_by_j)

    nist_fundamental = (
        float(TARGET["omega_e_cm-1"])
        - 2.0 * float(TARGET["omega_ex_e_cm-1"])
        + 3.25 * float(TARGET["omega_e_y_e_cm-1"])
    )

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "M2_1_LIH_ROVIBRATIONAL_SPECTRUM_BASELINE",
        "bo_curve": {
            "r_bohr": list(BO_GRID),
            "energy_hartree": list(energies),
        },
        "spectrum": spectrum.as_dict(),
        "vibrational_dunham": vib,
        "rotational_dunham": rot,
        "nist": TARGET,
        "derived_nist_fundamental_cm-1": nist_fundamental,
        "residuals": {
            "fundamental_v0_to_v1_cm-1": vib["fundamental_v0_to_v1_cm-1"] - nist_fundamental,
            "omega_e_cm-1": vib["omega_e_cm-1"] - float(TARGET["omega_e_cm-1"]),
            "omega_ex_e_cm-1": vib["omega_ex_e_cm-1"] - float(TARGET["omega_ex_e_cm-1"]),
            "omega_e_y_e_cm-1": vib["omega_e_y_e_cm-1"] - float(TARGET["omega_e_y_e_cm-1"]),
            "B_e_cm-1": rot["B_e_cm-1"] - float(TARGET["B_e_cm-1"]),
            "alpha_e_cm-1": rot["alpha_e_cm-1"] - float(TARGET["alpha_e_cm-1"]),
            "gamma_e_cm-1": rot["gamma_e_cm-1"] - float(TARGET["gamma_e_cm-1"]),
            "D_v0_cm-1": rot["D_v_cm-1"][0] - float(TARGET["D_e_cm-1"]),
        },
        "benchmark_usage": "OUTPUT_COMPARISON_ONLY_NOT_SOLVER_INPUT",
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if not all(x > 0.0 and math.isfinite(x) for x in spectrum.boundary_margin_cm_by_j.values()):
        raise RuntimeError("M2.1 BO interval confinement gate failed")
    j0 = np.asarray(spectrum.levels_hartree_by_j[0], dtype=float)
    if not np.all(np.diff(j0) > 0.0):
        raise RuntimeError("M2.1 vibrational levels are not strictly ordered")
    if not all(math.isfinite(float(x)) for x in vib.values()):
        raise RuntimeError("M2.1 vibrational Dunham fit is non-finite")
    if vib["fundamental_v0_to_v1_cm-1"] <= 0.0:
        raise RuntimeError("M2.1 vibrational fundamental is non-positive")
    if not all(math.isfinite(float(x)) and float(x) > 0.0 for x in rot["B_v_cm-1"]):
        raise RuntimeError("M2.1 rotational B_v values are invalid")
    if not math.isfinite(float(rot["B_e_cm-1"])) or float(rot["B_e_cm-1"]) <= 0.0:
        raise RuntimeError("M2.1 fitted B_e is invalid")


if __name__ == "__main__":
    main()
