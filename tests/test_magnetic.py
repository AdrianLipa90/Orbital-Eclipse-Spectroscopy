import unittest

from oes.hydrogen import MU_B_OVER_H
from oes.magnetic import p_mj_half_branches_hz, p_paschen_back_levels_hz, p_spin_orbit_constant_hz
from oes.relativity import fine_structure_split_hz


class MagneticCrossoverTests(unittest.TestCase):
    def test_zero_field_recovers_j_multiplets(self):
        a_so = p_spin_orbit_constant_hz(2)
        levels = p_paschen_back_levels_hz(2, 0.0)
        energies = sorted(round(float(row["energy_hz"]), 3) for row in levels)
        expected = sorted([round(-a_so, 3)] * 2 + [round(0.5 * a_so, 3)] * 4)
        self.assertEqual(energies, expected)
        self.assertAlmostEqual(1.5 * a_so, fine_structure_split_hz(2, 0.5, 1.5), places=3)

    def test_weak_field_recovers_lande_slopes_for_mj_half(self):
        b = 1e-6
        low0, high0 = p_mj_half_branches_hz(2, 0.0, sign=1)
        low1, high1 = p_mj_half_branches_hz(2, b, sign=1)
        low_slope = (low1 - low0) / b / MU_B_OVER_H
        high_slope = (high1 - high0) / b / MU_B_OVER_H
        # j=1/2: g*mj=(2/3)*(1/2)=1/3
        # j=3/2: g*mj=(4/3)*(1/2)=2/3
        self.assertAlmostEqual(low_slope, 1 / 3, places=5)
        self.assertAlmostEqual(high_slope, 2 / 3, places=5)

    def test_high_field_recovers_uncoupled_ml_ms_slopes(self):
        b = 100.0
        db = 1e-3
        low0, high0 = p_mj_half_branches_hz(2, b, sign=1)
        low1, high1 = p_mj_half_branches_hz(2, b + db, sign=1)
        low_slope = (low1 - low0) / db / MU_B_OVER_H
        high_slope = (high1 - high0) / db / MU_B_OVER_H
        # m_j=+1/2 block tends to |m_l=+1,m_s=-1/2> (slope 0)
        # and |m_l=0,m_s=+1/2> (slope 1).
        self.assertLess(abs(low_slope), 1e-3)
        self.assertLess(abs(high_slope - 1.0), 1e-3)


if __name__ == "__main__":
    unittest.main()
