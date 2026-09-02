#!/usr/bin/env python3
"""Run the OES-M2 LiH fixed-20Q heteronuclear molecular baseline."""

import json
import math
from pathlib import Path

from oes.quantum.h2_m1 import ANGSTROM_TO_BOHR, HARTREE_TO_WAVENUMBER_CM
from oes.quantum.lih_m2 import (
    LIH_REDUCED_NUCLEAR_MASS_ME,
    run_lih_curve,
    run_lih_point,
)

ROOT = Path(__file__).resolve().parents[1]
TARGET = json.loads((ROOT / "benchmarks" / "lih_m2_nist.json").read_text())


def rotational_constant_cm(bond_bohr: float) -> float:
    r = float(bond_bohr)
    if not math.isfinite(r) or r <= 0.0:
        raise ValueError("bond_bohr must be finite and positive")
    return HARTREE_TO_WAVENUMBER_CM / (2.0 * LIH_REDUCED_NUCLEAR_MASS_ME * r * r)


def main():
    benchmark_r = float(TARGET["re_angstrom"]) * ANGSTROM_TO_BOHR
    point = run_lih_point(benchmark_r, basis_name="cc-pVTZ")
    curve = run_lih_curve(basis_name="cc-pVTZ")
    predicted_b_e = rotational_constant_cm(curve.fitted_equilibrium_bohr)

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "M2_FIXED_20Q_HETERONUCLEAR_POLARIZED_TWO_CENTRE_BASELINE",
        "equilibrium_point": point.as_dict(),
        "curve": curve.as_dict(),
        "predicted_rotational_constant_cm-1": predicted_b_e,
        "nist": TARGET,
        "residuals": {
            "fitted_re_angstrom": curve.fitted_equilibrium_angstrom - float(TARGET["re_angstrom"]),
            "harmonic_wavenumber_cm": curve.harmonic_wavenumber_cm - float(TARGET["omega_e_cm-1"]),
            "rotational_constant_cm-1": predicted_b_e - float(TARGET["B_e_cm-1"]),
            "static_dipole_debye_vs_v0": point.dipole_debye - float(TARGET["dipole_v0_debye"]),
        },
        "benchmark_usage": "OUTPUT_COMPARISON_ONLY_NOT_SOLVER_INPUT",
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if abs(point.fci_delta_hartree) > 1e-9:
        raise RuntimeError(f"M2 LiH OES/PySCF FCI mismatch: {point.fci_delta_hartree} Ha")
    if point.n_spin_orbitals != 20 or point.n_electrons != 4:
        raise RuntimeError("M2 LiH canonical register/electron-count gate failed")
    if point.full_fixed_particle_dimension != 4_845:
        raise RuntimeError("M2 LiH fixed-N sector dimension gate failed")
    if point.ms_zero_dimension != 2_025:
        raise RuntimeError("M2 LiH M_S=0 sector dimension gate failed")
    if sum(point.selected_group_sizes) != 10:
        raise RuntimeError("M2 LiH complete-degeneracy-block active-space gate failed")
    if not math.isfinite(point.exposure_difference):
        raise RuntimeError("M2 LiH exposure asymmetry is non-finite")
    if not math.isfinite(point.dipole_debye) or point.dipole_debye <= 0.0:
        raise RuntimeError("M2 LiH polarized dipole gate failed")
    if not math.isfinite(curve.fitted_equilibrium_bohr) or curve.fitted_equilibrium_bohr <= 0.0:
        raise RuntimeError("M2 LiH equilibrium geometry gate failed")
    if not math.isfinite(curve.fitted_curvature_hartree_per_bohr2) or curve.fitted_curvature_hartree_per_bohr2 <= 0.0:
        raise RuntimeError("M2 LiH positive-curvature gate failed")
    if not math.isfinite(curve.harmonic_wavenumber_cm) or curve.harmonic_wavenumber_cm <= 0.0:
        raise RuntimeError("M2 LiH harmonic-frequency gate failed")
    if not math.isfinite(predicted_b_e) or predicted_b_e <= 0.0:
        raise RuntimeError("M2 LiH rotational-constant gate failed")


if __name__ == "__main__":
    main()
