import unittest

import numpy as np

from oes.quantum.hcl_m4 import HCL35_REDUCED_NUCLEAR_MASS_ME, rotational_constant_cm


class HClM4CoreTests(unittest.TestCase):
    def test_reduced_mass_is_positive(self):
        self.assertTrue(np.isfinite(HCL35_REDUCED_NUCLEAR_MASS_ME))
        self.assertGreater(HCL35_REDUCED_NUCLEAR_MASS_ME, 0.0)

    def test_rotational_constant_scales_as_inverse_r_squared(self):
        b1 = rotational_constant_cm(2.0)
        b2 = rotational_constant_cm(4.0)
        self.assertGreater(b1, b2)
        self.assertAlmostEqual(b1 / b2, 4.0, places=12)


if __name__ == "__main__":
    unittest.main()
