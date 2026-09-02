#!/usr/bin/env python3
"""Run the first spectroscopic HF baseline in the reduced-active 20Q model."""

import json
from pathlib import Path

from oes.quantum.hf_m3 import run_hf_curve, run_hf_point


BENCHMARK_PATH = Path(__file__).resolve().parents[1] / "benchmarks" / "molecules" / "nist_hf_constants.json"


def relative_error(value: float, reference: float) -> float:
    return (float(value) - float(reference)) / float(reference)


def main():
    benchmark = json.loads(BENCHMARK_PATH.read_text())

    # Blind electronic workflow: the calculation locates its own molecular
    # minimum. NIST data are loaded only after the model outputs exist.
    curve = run_hf_curve()
    point = run_hf_point(curve.fitted_equilibrium_bohr)

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
        "status_semantics": "M3_HF_FROZEN_4E_REDUCED_ACTIVE_20Q_BASELINE",
        "curve": curve.as_dict(),
        "equilibrium_point": point.as_dict(),
        "nist": benchmark,
        "residuals": residuals,
        "predeclared_baseline_relative_tolerance": 0.05,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    # Implementation/reference gates.
    if point.active_qubits != 20 or point.fixed_ms_dimension != 14400:
        raise RuntimeError("M3 fixed-20Q six-active-electron dimension gate failed")
    if abs(point.energy_delta_hartree) > 2e-9:
        raise RuntimeError("M3 OES sparse energy does not reproduce independent active-space FCI")
    if point.active_rdm1_max_delta > 5e-7:
        raise RuntimeError("M3 active 1-RDM does not reproduce independent FCI 1-RDM")
    if abs(point.exposure_difference) < 0.1:
        raise RuntimeError("M3 HF heteronuclear exposure asymmetry was not resolved")
    if point.dipole_debye < 0.1:
        raise RuntimeError("M3 HF permanent polarization was not resolved")

    # Physical baseline gate. Declared before observing the M3 result.
    tol = 0.05
    for key in ("re_relative", "omega_e_relative", "B_e_relative"):
        if abs(residuals[key]) > tol:
            raise RuntimeError(f"M3 physical 5% baseline failed for {key}: {residuals[key]}")


if __name__ == "__main__":
    main()
