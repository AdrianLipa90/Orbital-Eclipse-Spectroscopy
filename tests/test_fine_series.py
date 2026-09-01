import json
from pathlib import Path
import unittest

from oes.relativity import fine_structure_split_hz


BENCHMARK = Path(__file__).parents[1] / "benchmarks" / "hydrogen" / "nist_h1_np_fine_series.json"


class FineStructureSeriesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(BENCHMARK.read_text())

    def test_np_series_without_refit(self):
        for n in (2, 3, 4, 5):
            predicted = fine_structure_split_hz(n, 0.5, 1.5) / 1e9
            observed = self.data["splittings_GHz"][f"{n}p"]
            relative_error = abs(predicted - observed) / observed
            with self.subTest(n=n, predicted=predicted, observed=observed):
                self.assertLess(relative_error, 0.005)

    def test_series_is_consistent_with_leading_n_minus_3_scaling(self):
        values = self.data["splittings_GHz"]
        scaled = [values[f"{n}p"] * n**3 for n in (2, 3, 4, 5)]
        spread = (max(scaled) - min(scaled)) / (sum(scaled) / len(scaled))
        self.assertLess(spread, 0.01)


if __name__ == "__main__":
    unittest.main()
