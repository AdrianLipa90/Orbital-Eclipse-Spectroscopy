import math
import unittest

from oes.quantum.fermions import determinant_basis, full_space_dimension, sector_dimension
from oes.quantum.lithium_a1 import alpha_beta_counts, spin_sector_indices


class LithiumA1CoreTests(unittest.TestCase):
    def test_canonical_20q_three_electron_dimensions(self):
        self.assertEqual(full_space_dimension(20), 1_048_576)
        self.assertEqual(sector_dimension(20, 3), math.comb(20, 3))
        self.assertEqual(sector_dimension(20, 3), 1140)

    def test_ms_half_sector_dimension(self):
        basis = determinant_basis(20, 3)
        indices = spin_sector_indices(basis, 10, 2, 1)
        self.assertEqual(len(indices), math.comb(10, 2) * math.comb(10, 1))
        self.assertEqual(len(indices), 450)
        for i in indices[:25]:
            self.assertEqual(alpha_beta_counts(basis[i], 10), (2, 1))

    def test_ms_sector_and_partner_have_equal_dimensions(self):
        basis = determinant_basis(20, 3)
        plus = spin_sector_indices(basis, 10, 2, 1)
        minus = spin_sector_indices(basis, 10, 1, 2)
        self.assertEqual(len(plus), len(minus))
        self.assertEqual(len(plus), 450)


if __name__ == "__main__":
    unittest.main()
