import unittest

import numpy as np

from oes.quantum.diatomic_rovib import (
    fit_rotational_dunham,
    fit_vibrational_dunham,
    solve_diatomic_rovibrational_levels,
)
from oes.quantum.h2_m1 import HARTREE_TO_WAVENUMBER_CM


class DiatomicRovibrationalTests(unittest.TestCase):
    def test_harmonic_reference_spacing_and_boundary_confinement(self):
        mu = 1000.0
        omega = 0.01
        r = np.linspace(1.5, 4.5, 41)
        energy = 0.5 * mu * omega**2 * (r - 3.0) ** 2
        result = solve_diatomic_rovibrational_levels(
            r,
            energy,
            mu,
            n_vibrational=4,
            j_values=(0, 1, 2),
            radial_grid_points=1200,
        )
        levels = np.asarray(result.levels_hartree_by_j[0])
        self.assertAlmostEqual(levels[1] - levels[0], omega, delta=3e-5)
        self.assertTrue(all(x > 0.0 for x in result.boundary_margin_cm_by_j.values()))

        vib = fit_vibrational_dunham(result.term_values_cm_by_j[0])
        self.assertAlmostEqual(
            vib["fundamental_v0_to_v1_cm-1"],
            omega * HARTREE_TO_WAVENUMBER_CM,
            delta=10.0,
        )
        rot = fit_rotational_dunham(result.term_values_cm_by_j)
        self.assertGreater(rot["B_e_cm-1"], 0.0)
        self.assertTrue(all(x > 0.0 for x in rot["B_v_cm-1"]))

    def test_sampled_minimum_must_be_internal(self):
        r = np.linspace(1.0, 2.0, 9)
        energy = (r - 1.0) ** 2
        with self.assertRaises(RuntimeError):
            solve_diatomic_rovibrational_levels(r, energy, 1000.0)


if __name__ == "__main__":
    unittest.main()
