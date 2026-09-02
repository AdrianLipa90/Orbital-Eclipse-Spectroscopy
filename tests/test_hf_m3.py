import unittest

import numpy as np

from oes.quantum.hf_m3 import (
    HF_REDUCED_NUCLEAR_MASS_ME,
    _total_one_body_expectation,
    rotational_constant_cm,
)


class HFM3CoreTests(unittest.TestCase):
    def test_hf_reduced_mass_is_positive(self):
        self.assertTrue(np.isfinite(HF_REDUCED_NUCLEAR_MASS_ME))
        self.assertGreater(HF_REDUCED_NUCLEAR_MASS_ME, 0.0)

    def test_total_one_body_reconstructs_doubly_occupied_core(self):
        op = np.diag([2.0, 3.0, 5.0, 7.0])
        gamma = np.diag([1.25, 0.75])
        value = _total_one_body_expectation(op, gamma, n_core=2)
        expected = 2 * 2.0 + 2 * 3.0 + 1.25 * 5.0 + 0.75 * 7.0
        self.assertAlmostEqual(value, expected)

    def test_rotational_constant_decreases_with_bond_length_squared(self):
        b1 = rotational_constant_cm(1.5)
        b2 = rotational_constant_cm(3.0)
        self.assertGreater(b1, b2)
        self.assertAlmostEqual(b1 / b2, 4.0, places=12)


if __name__ == "__main__":
    unittest.main()
