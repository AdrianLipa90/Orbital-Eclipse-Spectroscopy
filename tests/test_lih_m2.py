import math
import unittest

from oes.quantum.fermions import full_space_dimension, sector_dimension
from oes.quantum.lih_m2 import LIH_REDUCED_NUCLEAR_MASS_ME, lih_ms_zero_basis
from oes.quantum.lithium_a1 import alpha_beta_counts


class LiHM2CoreTests(unittest.TestCase):
    def test_canonical_20q_four_electron_dimensions(self):
        self.assertEqual(full_space_dimension(20), 1_048_576)
        self.assertEqual(sector_dimension(20, 4), math.comb(20, 4))
        self.assertEqual(sector_dimension(20, 4), 4_845)

    def test_ms_zero_sector_dimension_and_spin_counts(self):
        basis = lih_ms_zero_basis(10)
        self.assertEqual(len(basis), math.comb(10, 2) ** 2)
        self.assertEqual(len(basis), 2_025)
        for det in basis:
            self.assertEqual(int(det).bit_count(), 4)
            self.assertEqual(alpha_beta_counts(int(det), 10), (2, 2))

    def test_lih_reduced_nuclear_mass_is_finite_positive(self):
        self.assertTrue(math.isfinite(LIH_REDUCED_NUCLEAR_MASS_ME))
        self.assertGreater(LIH_REDUCED_NUCLEAR_MASS_ME, 0.0)


if __name__ == "__main__":
    unittest.main()
