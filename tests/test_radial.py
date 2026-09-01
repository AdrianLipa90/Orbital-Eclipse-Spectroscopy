import math
import unittest

from oes.radial import radial_transition_integrals, radial_wavefunction_dimensionless


class RadialEclipseTests(unittest.TestCase):
    def test_1s_value_at_origin(self):
        self.assertAlmostEqual(radial_wavefunction_dimensionless(1, 0, 0.0), 2.0, places=12)

    def test_2p_1s_exact_radial_dipole_integral(self):
        result = radial_transition_integrals(2, 1, 1, 0, steps=8_000)
        exact = 256.0 / (81.0 * math.sqrt(6.0))
        self.assertAlmostEqual(result["signed_aH"], exact, places=8)
        self.assertAlmostEqual(result["coherence"], 1.0, places=10)

    def test_halpha_channels_have_distinct_radial_flavors(self):
        s_to_p = radial_transition_integrals(3, 0, 2, 1, steps=8_000)
        p_to_s = radial_transition_integrals(3, 1, 2, 0, steps=8_000)
        d_to_p = radial_transition_integrals(3, 2, 2, 1, steps=8_000)
        self.assertAlmostEqual(s_to_p["signed_aH"], 0.93840424, places=6)
        self.assertAlmostEqual(p_to_s["signed_aH"], 3.06481541, places=6)
        self.assertAlmostEqual(d_to_p["signed_aH"], 4.74799161, places=6)
        self.assertLess(s_to_p["coherence"], p_to_s["coherence"])
        self.assertLess(p_to_s["coherence"], d_to_p["coherence"])
        self.assertAlmostEqual(d_to_p["coherence"], 1.0, places=8)

    def test_no_node_high_l_channel_is_coherent(self):
        result = radial_transition_integrals(4, 3, 3, 2, steps=8_000)
        self.assertAlmostEqual(result["signed_aH"], 10.23030262, places=6)
        self.assertAlmostEqual(result["coherence"], 1.0, places=8)


if __name__ == "__main__":
    unittest.main()
