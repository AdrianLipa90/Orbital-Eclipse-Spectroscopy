"""Emit T4 transition-density subspace continuity for rigid Ne...Cl2."""
from __future__ import annotations

import json
import math

from oes.molecular_tda import run_closed_shell_tda_runtime
from oes.tda_tracking import runtime_block_continuity

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
    runtimes = []
    for angle_degrees in ANGLES_DEGREES:
        theta = math.radians(angle_degrees)
        runtimes.append(
            run_closed_shell_tda_runtime(
                species_label=f"NeCl2-{angle_degrees:g}deg",
                atoms=("Ne", "Cl", "Cl"),
                coordinates_angstrom=geometry(theta),
                basis_name="sto-3g",
                nstates=4,
            )
        )

    adjacent = []
    for i in range(len(runtimes) - 1):
        left = runtimes[i]
        right = runtimes[i + 1]
        bright = runtime_block_continuity(
            left,
            right,
            left_roots=(1, 2),
            right_roots=(1, 2),
        )
        dark = runtime_block_continuity(
            left,
            right,
            left_roots=(3, 4),
            right_roots=(3, 4),
        )
        adjacent.append(
            {
                "left_angle_degrees": ANGLES_DEGREES[i],
                "right_angle_degrees": ANGLES_DEGREES[i + 1],
                "bright_block_roots": [1, 2],
                "bright_principal_cosines": list(bright.principal_cosines),
                "bright_minimum_principal_cosine": bright.minimum_principal_cosine,
                "bright_chordal_distance": bright.chordal_distance,
                "dark_block_roots": [3, 4],
                "dark_principal_cosines": list(dark.principal_cosines),
                "dark_minimum_principal_cosine": dark.minimum_principal_cosine,
                "dark_chordal_distance": dark.chordal_distance,
            }
        )

    payload = {
        "schema": "OES_T4_NECL2_TRANSITION_DENSITY_CONTINUITY_V0_1",
        "backend": "PYSCF_RHF_TDA",
        "basis": "sto-3g",
        "geometry_contract": {
            "atoms": ["Ne", "Cl", "Cl"],
            "r_yy_angstrom": R_YY_ANGSTROM,
            "x_to_yy_midpoint_angstrom": X_TO_MIDPOINT_ANGSTROM,
            "angles_degrees": list(ANGLES_DEGREES),
            "varied_coordinate": "relative_orientation_only",
        },
        "tracking_object": "GRAM_WHITENED_AO_TRANSITION_DENSITY_SPAN",
        "root_labels_are_identity": False,
        "experimental_inputs": [],
        "adjacent_pairs": adjacent,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
