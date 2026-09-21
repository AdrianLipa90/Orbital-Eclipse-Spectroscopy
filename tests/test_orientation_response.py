import math
import unittest

import numpy as np

from oes.orientation_response import (
    OrientationResponseError,
    orientation_line_point,
    orientation_line_scan,
    oscillator_strength_slope,
)
from oes.transition_contract import (
    ElectronicStateRef,
    NuclearCenter,
    OESTransitionState,
    TransitionProvenance,
)


class OrientationResponseTests(unittest.TestCase):
    def _transition(self, theta: float) -> OESTransitionState:
        t = np.asarray(
            [[0.0, math.cos(theta)], [0.0, 0.0]],
            dtype=complex,
        )
        return OESTransitionState(
            species_label="PP-CONTROL",
            nuclei=(
                NuclearCenter("A", 1, 1.0, (-1.0, 0.0, 0.0)),
                NuclearCenter("B", 1, 1.0, (1.0, 0.0, 0.0)),
            ),
            initial=ElectronicStateRef("i", -1.0),
            final=ElectronicStateRef("f", -0.6),
            transition_rdm=t,
            provenance=TransitionProvenance(
                source_repository="analytic-control",
                source_commit="fixture",
                method="PP_ORIENTATION_CONTROL",
                orbital_basis_id="fixed-two-orbital-control",
                backend_status="ANALYTIC_CONTROL",
            ),
        )

    def _dipole(self):
        operators = np.zeros((3, 2, 2), dtype=complex)
        operators[0, 0, 1] = 1.0
        return operators

    def test_orientation_line_strength_follows_cosine_squared_control(self):
        for theta in (0.0, math.pi / 6.0, math.pi / 3.0):
            point = orientation_line_point(
                theta,
                self._transition(theta),
                self._dipole(),
            )
            expected = (
                (2.0 / 3.0)
                * 0.4
                * math.cos(theta) ** 2
            )
            self.assertAlmostEqual(
                point.oscillator_strength_length_gauge,
                expected,
                places=12,
            )

    def test_quarter_turn_control_is_dark(self):
        point = orientation_line_point(
            math.pi / 2.0,
            self._transition(math.pi / 2.0),
            self._dipole(),
        )
        self.assertLess(
            point.oscillator_strength_length_gauge,
            1.0e-30,
        )

    def test_scan_preserves_angle_order(self):
        angles = (0.7, 0.1, 1.1)
        points = orientation_line_scan(
            (
                (theta, self._transition(theta), self._dipole())
                for theta in angles
            )
        )
        self.assertEqual(
            tuple(point.angle_rad for point in points),
            angles,
        )

    def test_oscillator_strength_slope_decomposes_exactly(self):
        mu = np.asarray([0.8 + 0.1j, -0.2j, 0.1], dtype=complex)
        dmu = np.asarray([-0.2 + 0.03j, 0.05j, -0.04], dtype=complex)
        gap = 0.4
        gap_slope = 0.05

        result = oscillator_strength_slope(
            energy_gap_hartree=gap,
            energy_gap_slope_hartree_per_rad=gap_slope,
            transition_dipole_au=mu,
            transition_dipole_slope_au_per_rad=dmu,
        )

        expected_energy = (
            (2.0 / 3.0)
            * gap_slope
            * float(np.vdot(mu, mu).real)
        )
        expected_moment = (
            (4.0 / 3.0)
            * gap
            * float(np.vdot(mu, dmu).real)
        )
        self.assertAlmostEqual(
            result.energy_gap_term_per_rad,
            expected_energy,
            places=12,
        )
        self.assertAlmostEqual(
            result.transition_moment_term_per_rad,
            expected_moment,
            places=12,
        )
        self.assertAlmostEqual(
            result.total_per_rad,
            expected_energy + expected_moment,
            places=12,
        )

    def test_fail_closed_domains(self):
        with self.assertRaises(OrientationResponseError):
            orientation_line_scan(())
        with self.assertRaises(OrientationResponseError):
            oscillator_strength_slope(
                energy_gap_hartree=0.0,
                energy_gap_slope_hartree_per_rad=0.1,
                transition_dipole_au=(1.0, 0.0, 0.0),
                transition_dipole_slope_au_per_rad=(0.0, 0.0, 0.0),
            )


if __name__ == "__main__":
    unittest.main()
