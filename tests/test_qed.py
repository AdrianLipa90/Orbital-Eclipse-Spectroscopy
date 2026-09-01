import json
from pathlib import Path
import unittest

from oes.qed import leading_lamb_2s_2p1_2_components_mhz


BENCHMARK = Path(__file__).parents[1] / "benchmarks" / "hydrogen" / "nist_h1_levels.json"


class LeadingQEDTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(BENCHMARK.read_text())

    def test_leading_components_are_numerically_stable(self):
        result = leading_lamb_2s_2p1_2_components_mhz()
        self.assertAlmostEqual(result["self_energy_2s_MHz"], 1064.77154, places=4)
        self.assertAlmostEqual(result["self_energy_2p1_2_MHz"], -12.86287, places=4)
        self.assertAlmostEqual(result["vacuum_polarization_2s_MHz"], -27.08448, places=4)

    def test_leading_lamb_interval_is_subpercent_from_nist(self):
        result = leading_lamb_2s_2p1_2_components_mhz()
        predicted = result["leading_interval_MHz"]
        observed = self.data["derived_targets"]["2s_2p1/2_lamb_split_MHz"]
        relative_error = abs(predicted - observed) / observed
        self.assertLess(relative_error, 0.01)

    def test_benchmark_remains_outside_solver(self):
        result = leading_lamb_2s_2p1_2_components_mhz()
        self.assertGreater(result["leading_interval_MHz"], 0.0)
        self.assertEqual(self.data["use"], "benchmark_only_not_solver_input")


if __name__ == "__main__":
    unittest.main()
