"""Emit block-resolved spectral response for rigid Ne...Cl2 orientation."""
from __future__ import annotations

import json
import math

from oes.block_response import (
    block_spectral_observables,
    tracked_block_response,
)
from oes.molecular_tda import run_closed_shell_tda_runtime

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
        runtimes.append(
            run_closed_shell_tda_runtime(
                species_label=f"NeCl2-{angle_degrees:g}deg",
                atoms=("Ne", "Cl", "Cl"),
                coordinates_angstrom=geometry(
                    math.radians(angle_degrees)
                ),
                basis_name="sto-3g",
                nstates=4,
            )
        )

    rows = []
    computed_window_strengths = []
    for angle_degrees, runtime in zip(
        ANGLES_DEGREES,
        runtimes,
    ):
        bright = block_spectral_observables(
            runtime.result,
            (1, 2),
        )
        dark = block_spectral_observables(
            runtime.result,
            (3, 4),
        )
        total = (
            bright.oscillator_strength_sum
            + dark.oscillator_strength_sum
        )
        computed_window_strengths.append(total)
        rows.append(
            {
                "angle_degrees": angle_degrees,
                "bright_block": bright.as_dict(),
                "dark_block": dark.as_dict(),
                "computed_first_four_strength_sum": total,
            }
        )

    adjacent = []
    for i in range(len(runtimes) - 1):
        bright = tracked_block_response(
            runtimes[i],
            runtimes[i + 1],
            left_roots=(1, 2),
            right_roots=(1, 2),
        )
        dark = tracked_block_response(
            runtimes[i],
            runtimes[i + 1],
            left_roots=(3, 4),
            right_roots=(3, 4),
        )
        adjacent.append(
            {
                "left_angle_degrees": ANGLES_DEGREES[i],
                "right_angle_degrees": ANGLES_DEGREES[i + 1],
                "bright": bright.as_dict(),
                "dark": dark.as_dict(),
            }
        )

    dark_strengths = [
        row["dark_block"]["oscillator_strength_sum"]
        for row in rows
    ]
    max_dark_index = max(
        range(len(dark_strengths)),
        key=dark_strengths.__getitem__,
    )
    baseline_dark = dark_strengths[0]
    amplification = None
    if baseline_dark > 0.0:
        amplification = (
            dark_strengths[max_dark_index]
            / baseline_dark
        )

    minimum_total = min(computed_window_strengths)
    maximum_total = max(computed_window_strengths)
    mean_total = (
        sum(computed_window_strengths)
        / len(computed_window_strengths)
    )

    payload = {
        "schema": "OES_T5_NECL2_BLOCK_SPECTRAL_RESPONSE_V0_1",
        "backend": "PYSCF_RHF_TDA",
        "basis": "sto-3g",
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
            "root_labels_are_identity": False,
            "continuity_object": (
                "GRAM_WHITENED_AO_TRANSITION_DENSITY_SPAN"
            ),
        },
        "experimental_inputs": [],
        "rows": rows,
        "adjacent_responses": adjacent,
        "derived_control_readout": {
            "maximum_dark_block_strength_angle_degrees": (
                ANGLES_DEGREES[max_dark_index]
            ),
            "maximum_dark_block_strength": (
                dark_strengths[max_dark_index]
            ),
            "dark_block_amplification_relative_to_0deg": (
                amplification
            ),
            "computed_first_four_strength_min": minimum_total,
            "computed_first_four_strength_max": maximum_total,
            "computed_first_four_strength_relative_span": (
                (maximum_total - minimum_total) / mean_total
                if mean_total > 0.0
                else None
            ),
        },
        "interpretation_status": (
            "CONTROL_ONLY_NO_EXPERIMENTAL_VISIBLE_COLOR_CLAIM"
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
