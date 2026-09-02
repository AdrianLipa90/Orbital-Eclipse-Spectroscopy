import unittest

import numpy as np

from oes.quantum.hcl_m4 import (
    HCL35_REDUCED_NUCLEAR_MASS_ME,
    _adaptive_bracketed_quadratic,
    rotational_constant_cm,
)


class HClM4CoreTests(unittest.TestCase):
    def test_reduced_mass_is_positive(self):
        self.assertTrue(np.isfinite(HCL35_REDUCED_NUCLEAR_MASS_ME))
        self.assertGreater(HCL35_REDUCED_NUCLEAR_MASS_ME, 0.0)

    def test_rotational_constant_scales_as_inverse_r_squared(self):
        b1 = rotational_constant_cm(2.0)
        b2 = rotational_constant_cm(4.0)
        self.assertGreater(b1, b2)
        self.assertAlmostEqual(b1 / b2, 4.0, places=12)

    def test_active_curve_recenters_when_initial_minimum_is_at_edge(self):
        target = 2.42

        def energy(r):
            return 1.5 + (r - target) ** 2

        grid, energies, r_eq, curvature, recenter_count = _adaptive_bracketed_quadratic(
            energy,
            seed_bohr=2.30,
            half_width_bohr=0.12,
            max_recenters=4,
        )

        self.assertGreaterEqual(recenter_count, 1)
        self.assertLess(grid[1], r_eq)
        self.assertGreater(grid[-2], r_eq)
        self.assertAlmostEqual(r_eq, target, places=12)
        self.assertAlmostEqual(curvature, 2.0, places=11)
        self.assertEqual(len(energies), 5)


if __name__ == "__main__":
    unittest.main()
