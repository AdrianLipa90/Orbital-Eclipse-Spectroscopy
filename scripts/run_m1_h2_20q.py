#!/usr/bin/env python3
"""Run the OES-M1 H2 fixed-20Q two-centre baseline."""

import json
from pathlib import Path

from oes.quantum.h2_m1 import ANGSTROM_TO_BOHR, run_h2_curve, run_h2_point

ROOT = Path(__file__).resolve().parents[1]
TARGET = json.loads((ROOT / "benchmarks" / "h2_m1_nist.json").read_text())


def main():
    benchmark_r = float(TARGET["re_angstrom"]) * ANGSTROM_TO_BOHR
    point = run_h2_point(benchmark_r, basis_name="cc-pVTZ")
    curve = run_h2_curve(basis_name="cc-pVTZ")

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "M1_FIXED_20Q_TWO_CENTRE_SYMMETRY_AND_BONDING_BASELINE",
        "equilibrium_point": point.as_dict(),
        "curve": curve.as_dict(),
        "nist": TARGET,
        "residuals": {
            "fitted_re_angstrom": curve.fitted_equilibrium_angstrom - float(TARGET["re_angstrom"]),
            "harmonic_wavenumber_cm": curve.harmonic_wavenumber_cm - float(TARGET["omega_e_cm-1"]),
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if abs(point.fci_delta_hartree) > 1e-9:
        raise RuntimeError(f"M1 H2 OES/PySCF FCI mismatch: {point.fci_delta_hartree} Ha")
    if point.fixed_particle_dimension != 190:
        raise RuntimeError("M1 H2 fixed-N sector dimension gate failed")
    if abs(point.ground_parity - 1.0) > 1e-8:
        raise RuntimeError(f"M1 H2 ground-state g-parity gate failed: {point.ground_parity}")
    if abs(point.first_bright_parity + 1.0) > 1e-8:
        raise RuntimeError(f"M1 H2 first-bright u-parity gate failed: {point.first_bright_parity}")
    if abs(point.center_exposure_difference) > 1e-8:
        raise RuntimeError(
            f"M1 H2 equivalent-centre exposure gate failed: {point.center_exposure_difference}"
        )
    if abs(point.ground_s2) > 1e-8:
        raise RuntimeError(f"M1 H2 ground singlet gate failed: S2={point.ground_s2}")


if __name__ == "__main__":
    main()
