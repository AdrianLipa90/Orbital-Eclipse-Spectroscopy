"""Emit T6 spatial transition-density decomposition for rigid Ne...Cl2."""
from __future__ import annotations

import json
import math

from oes.block_response import block_spectral_observables
from oes.molecular_tda import run_closed_shell_tda_runtime
from oes.spatial_transition_density import (
    evaluate_runtime_block_spatial,
)

R_YY_ANGSTROM = 2.0
X_TO_MIDPOINT_ANGSTROM = 4.0
ANGLES_DEGREES = (0.0, 30.0, 60.0, 90.0)
GRID_LEVEL = 3


def geometry(angle_rad: float):
    half = 0.5 * R_YY_ANGSTROM
    d = X_TO_MIDPOINT_ANGSTROM
    return (
        (d * math.sin(angle_rad), 0.0, d * math.cos(angle_rad)),
        (0.0, 0.0, -half),
        (0.0, 0.0, half),
    )


def ratio(numerator: float, denominator: float):
    if denominator == 0.0:
        return None
    return numerator / denominator


def main() -> None:
    rows = []
    for angle_degrees in ANGLES_DEGREES:
        runtime = run_closed_shell_tda_runtime(
            species_label=f"NeCl2-{angle_degrees:g}deg",
            atoms=("Ne", "Cl", "Cl"),
            coordinates_angstrom=geometry(
                math.radians(angle_degrees)
            ),
            basis_name="sto-3g",
            nstates=4,
        )
        bright_spatial = evaluate_runtime_block_spatial(
            runtime,
            (1, 2),
            grid_level=GRID_LEVEL,
        )
        dark_spatial = evaluate_runtime_block_spatial(
            runtime,
            (3, 4),
            grid_level=GRID_LEVEL,
        )
        bright_spectral = block_spectral_observables(
            runtime.result,
            (1, 2),
        )
        dark_spectral = block_spectral_observables(
            runtime.result,
            (3, 4),
        )
        rows.append(
            {
                "angle_degrees": angle_degrees,
                "bright_spatial": bright_spatial.as_dict(),
                "dark_spatial": dark_spatial.as_dict(),
                "bright_oscillator_strength_sum": (
                    bright_spectral.oscillator_strength_sum
                ),
                "dark_oscillator_strength_sum": (
                    dark_spectral.oscillator_strength_sum
                ),
            }
        )

    by_angle = {
        row["angle_degrees"]: row
        for row in rows
    }
    dark0 = by_angle[0.0]["dark_spatial"]
    dark30 = by_angle[30.0]["dark_spatial"]

    payload = {
        "schema": "OES_T6_NECL2_SPATIAL_CANCELLATION_V0_1",
        "backend": "PYSCF_RHF_TDA",
        "basis": "sto-3g",
        "grid_level": GRID_LEVEL,
        "geometry_contract": {
            "atoms": ["Ne", "Cl", "Cl"],
            "r_yy_angstrom": R_YY_ANGSTROM,
            "x_to_yy_midpoint_angstrom": X_TO_MIDPOINT_ANGSTROM,
            "angles_degrees": list(ANGLES_DEGREES),
            "varied_coordinate": "relative_orientation_only",
        },
        "block_contract": {
            "bright_roots": [1, 2],
            "dark_roots": [3, 4],
            "root_basis_invariant_fields": [
                "sum_abs_transition_density_squared",
                "dipole_envelope",
                "block_dipole_frobenius_norm",
                "dipole_survival_ratio",
            ],
        },
        "experimental_inputs": [],
        "rows": rows,
        "dark_30deg_over_0deg": {
            "oscillator_strength_ratio": ratio(
                by_angle[30.0]["dark_oscillator_strength_sum"],
                by_angle[0.0]["dark_oscillator_strength_sum"],
            ),
            "net_block_dipole_norm_ratio": ratio(
                dark30["net_block_dipole_norm_au"],
                dark0["net_block_dipole_norm_au"],
            ),
            "dipole_envelope_ratio": ratio(
                dark30["dipole_envelope_au"],
                dark0["dipole_envelope_au"],
            ),
            "dipole_survival_ratio_factor": ratio(
                dark30["dipole_survival_ratio"],
                dark0["dipole_survival_ratio"],
            ),
            "transition_density_l2_activity_ratio": ratio(
                dark30[
                    "transition_density_l2_activity_bohr_minus3"
                ],
                dark0[
                    "transition_density_l2_activity_bohr_minus3"
                ],
            ),
        },
        "interpretation_status": (
            "CONTROL_ONLY_SPATIAL_MECHANISM_DIAGNOSTIC"
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
