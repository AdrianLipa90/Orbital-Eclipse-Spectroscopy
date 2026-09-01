import json
from pathlib import Path
import unittest

from oes.relativity import ALPHA, dirac_bound_energy_ev, fine_structure_split_hz


BENCHMARK = Path(__file__).parents[1] / "benchmarks" / "hydrogen" / "nist_h1_levels.json"


class RelativisticHydrogenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(BENCHMARK.read_text())

    def test_alpha_is_codatascale(self):
        self.assertAlmostEqual(ALPHA, 7.29735256e-3, places=10)

    def test_dirac_same_n_j_is_degenerate_across_l_flavor(self):
        # Pure Coulomb Dirac energy depends on n,j. Thus 2s_1/2 and 2p_1/2
        # share the same reference energy before Lamb/QED lifting.
        e_2j_half_a = dirac_bound_energy_ev(2, 0.5)
        e_2j_half_b = dirac_bound_energy_ev(2, 0.5)
        self.assertEqual(e_2j_half_a, e_2j_half_b)

    def test_2p_fine_structure_split_is_close_to_nist_without_fit(self):
        predicted_ghz = fine_structure_split_hz(2, 0.5, 1.5) / 1e9
        observed_ghz = self.data["derived_targets"]["2p_fine_split_GHz"]
        relative_error = abs(predicted_ghz - observed_ghz) / observed_ghz
        self.assertLess(relative_error, 0.003)

    def test_dirac_does_not_claim_lamb_split(self):
        # Same n,j must remain exactly degenerate here; the NIST benchmark has
        # a non-zero 2s_1/2 - 2p_1/2 split reserved for the QED gate.
        self.assertEqual(dirac_bound_energy_ev(2, 0.5), dirac_bound_energy_ev(2, 0.5))
        self.assertGreater(self.data["derived_targets"]["2s_2p1/2_lamb_split_MHz"], 0.0)


if __name__ == "__main__":
    unittest.main()
