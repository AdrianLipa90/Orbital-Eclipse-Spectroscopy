"""Emit a benchmark-blind rigid Ne...Cl2 orientation TDA receipt."""
from __future__ import annotations

import json
import math

import numpy as np

from oes.molecular_tda import run_orientation_tda_scan

R_YY_ANGSTROM = 2.0
X_TO_MIDPOINT_ANGSTROM = 4.0
ANGLES_DEGREES = (0.0, 30.0, 60.0, 90.0)


def geometry(angle_rad: float):
    half = 0.5 * R_YY_ANGSTROM
    d = X_TO_MIDPOINT_ANGSTROM
    return (
        (d * math.sin(angle_rad), 0.0, d * math.cos(angle_rad)),
        (0.0, 0.0, -half),
        (0.0, 0.0, half),
    )


def main() -> None:
    samples = []
    for degrees_value in ANGLES_DEGREES:
        theta = math.radians(degrees_value)
        coords = np.asarray(geometry(theta), dtype=float)
        midpoint = 0.5 * (coords[1] + coords[2])
        if not math.isclose(
            float(np.linalg.norm(coords[2] - coords[1])),
            R_YY_ANGSTROM,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise RuntimeError("rigid NeCl2 control changed Cl-Cl distance")
        if not math.isclose(
            float(np.linalg.norm(coords[0] - midpoint)),
            X_TO_MIDPOINT_ANGSTROM,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise RuntimeError("rigid NeCl2 control changed Ne-midpoint distance")
        samples.append((theta, tuple(tuple(float(x) for x in row) for row in coords)))

    scan = run_orientation_tda_scan(
        species_label="NeCl2-rigid-orientation-control",
        atoms=("Ne", "Cl", "Cl"),
        samples=tuple(samples),
        basis_name="sto-3g",
        nstates=4,
    )

    rows = []
    for degrees_value, point in zip(ANGLES_DEGREES, scan.points):
        rows.append(
            {
                "angle_degrees": degrees_value,
                "angle_rad": point.angle_rad,
                "rhf_energy_hartree": point.result.rhf_energy_hartree,
                "states": [
                    {
                        "root": state.root,
                        "excitation_ev": state.excitation_ev,
                        "transition_dipole_norm_au": state.transition_dipole_norm_au,
                        "oscillator_strength": state.oscillator_strength_oes,
                        "backend_delta": state.oscillator_strength_delta,
                    }
                    for state in point.result.states
                ],
                "sum_oscillator_strength": sum(
                    state.oscillator_strength_oes
                    for state in point.result.states
                ),
            }
        )

    first = np.asarray(
        [state["excitation_ev"] for state in rows[0]["states"]],
        dtype=float,
    )
    last = np.asarray(
        [state["excitation_ev"] for state in rows[-1]["states"]],
        dtype=float,
    )
    payload = {
        "schema": "OES_T3_NECL2_RIGID_ORIENTATION_RECEIPT_V0_1",
        "backend": "PYSCF_RHF_TDA",
        "basis": "sto-3g",
        "geometry_contract": {
            "atoms": ["Ne", "Cl", "Cl"],
            "r_yy_angstrom": R_YY_ANGSTROM,
            "x_to_yy_midpoint_angstrom": X_TO_MIDPOINT_ANGSTROM,
            "angles_degrees": list(ANGLES_DEGREES),
            "varied_coordinate": "relative_orientation_only",
        },
        "state_tracking_status": scan.state_tracking_status,
        "experimental_inputs": [],
        "raw_sorted_line_set_change_0_to_90_ev": float(
            np.max(np.abs(first - last))
        ),
        "rows": rows,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
