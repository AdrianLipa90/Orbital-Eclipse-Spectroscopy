#!/usr/bin/env python3
"""Run the blind H35Cl eight-active-electron fixed-20Q spectroscopy gate."""

import json
from pathlib import Path

from oes.quantum.hcl_m4 import run_hcl_curve, run_hcl_point


BENCHMARK_PATH = Path(__file__).resolve().parents[1] / "benchmarks" / "molecules" / "nist_hcl_constants.json"


def relative_error(value: float, reference: float) -> float:
    return (float(value) - float(reference)) / float(reference)


def main():
    # Compute first; benchmark constants are used only after the blind model
    # outputs exist.
    curve = run_hcl_curve()
    point = run_hcl_point(curve.fitted_equilibrium_bohr)
    benchmark = json.loads(BENCHMARK_PATH.read_text())

    residuals = {
        "re_angstrom": curve.fitted_equilibrium_angstrom - benchmark["re_angstrom"],
        "re_relative": relative_error(curve.fitted_equilibrium_angstrom, benchmark["re_angstrom"]),
        "omega_e_cm-1": curve.harmonic_wavenumber_cm - benchmark["omega_e_cm-1"],
        "omega_e_relative": relative_error(curve.harmonic_wavenumber_cm, benchmark["omega_e_cm-1"]),
        "B_e_cm-1": curve.rotational_constant_cm - benchmark["B_e_cm-1"],
        "B_e_relative": relative_error(curve.rotational_constant_cm, benchmark["B_e_cm-1"]),
    }
    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "M4_H35CL_FROZEN_10E_EIGHT_ACTIVE_ELECTRON_20Q_BASELINE",
        "curve": curve.as_dict(),
        "equilibrium_point": point.as_dict(),
        "nist": benchmark,
        "residuals": residuals,
        "predeclared_baseline_relative_tolerance": 0.05,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if point.active_qubits != 20 or point.fixed_ms_dimension != 44100:
        raise RuntimeError("M4 fixed-20Q eight-electron dimension gate failed")
    if abs(point.energy_delta_hartree) > 5e-8:
        raise RuntimeError("M4 OES sparse energy does not reproduce independent active FCI")
    if point.active_rdm1_max_delta > 1e-6:
        raise RuntimeError("M4 active 1-RDM does not reproduce independent active FCI")
    if abs(point.exposure_difference) < 0.1:
        raise RuntimeError("M4 HCl heteronuclear exposure asymmetry was not resolved")
    if point.dipole_debye < 0.1:
        raise RuntimeError("M4 HCl permanent polarization was not resolved")

    # Same predeclared cross-species physical threshold as M3 HF.
    tol = 0.05
    for key in ("re_relative", "omega_e_relative", "B_e_relative"):
        if abs(residuals[key]) > tol:
            raise RuntimeError(f"M4 physical 5% baseline failed for {key}: {residuals[key]}")


if __name__ == "__main__":
    main()
