import unittest

import numpy as np

from oes.holonomy_response import (
    HolonomyResponseError,
    fisher_information,
    identifiability_report,
    identifiable_rank,
    intensity_phase_slope,
    line_frequency_phase_slope,
    transition_frequency,
)


class HolonomyResponseTests(unittest.TestCase):
    def test_transition_frequency_and_feynman_hellmann_slope(self):
        hbar = 2.0
        self.assertAlmostEqual(transition_frequency(7.0, 3.0, hbar), 2.0)
        self.assertAlmostEqual(
            line_frequency_phase_slope(5.0, 1.0, hbar), -2.0
        )

    def test_intensity_phase_slope(self):
        mu = 1.0 + 2.0j
        dmu = -0.5 + 0.25j
        expected = 2.0 * np.real(np.conjugate(mu) * dmu)
        self.assertAlmostEqual(intensity_phase_slope(mu, dmu), expected)

    def test_full_rank_jacobian_is_locally_identifiable(self):
        jac = np.array(
            [
                [1.0, 0.0],
                [0.0, 2.0],
                [1.0, 1.0],
            ]
        )
        sigma = np.diag([1.0, 2.0, 3.0])
        report = identifiability_report(jac, sigma)
        self.assertEqual(report.jacobian_rank, 2)
        self.assertEqual(report.fisher_rank, 2)
        self.assertTrue(report.fully_locally_identifiable)

    def test_collinear_observable_rows_do_not_create_phase_rank(self):
        jac = np.array(
            [
                [1.0, 2.0],
                [2.0, 4.0],
                [-3.0, -6.0],
            ]
        )
        self.assertEqual(identifiable_rank(jac), 1)
        report = identifiability_report(jac, np.eye(3))
        self.assertEqual(report.jacobian_rank, 1)
        self.assertFalse(report.fully_locally_identifiable)

    def test_invertible_phase_basis_change_preserves_rank(self):
        jac = np.array(
            [
                [1.0, 0.0],
                [0.0, 1.0],
                [1.0, -1.0],
            ]
        )
        basis_change = np.array([[2.0, 1.0], [1.0, 1.0]])
        self.assertNotEqual(float(np.linalg.det(basis_change)), 0.0)
        self.assertEqual(
            identifiable_rank(jac), identifiable_rank(jac @ basis_change)
        )

    def test_fisher_matches_direct_inverse_definition(self):
        jac = np.array([[1.0, 2.0], [3.0, 1.0]])
        sigma = np.array([[2.0, 0.3], [0.3, 1.5]])
        expected = jac.T @ np.linalg.inv(sigma) @ jac
        np.testing.assert_allclose(
            fisher_information(jac, sigma), expected, rtol=1e-12, atol=1e-12
        )

    def test_covariance_fail_closed(self):
        covariances = [
            np.array([[1.0, 2.0], [0.0, 1.0]]),
            np.array([[1.0, 0.0], [0.0, 0.0]]),
            np.array([[1.0, 2.0], [2.0, 1.0]]),
        ]
        jac = np.eye(2)
        for covariance in covariances:
            with self.subTest(covariance=covariance):
                with self.assertRaises(HolonomyResponseError):
                    fisher_information(jac, covariance)

    def test_basic_domains_fail_closed(self):
        calls = [
            lambda: transition_frequency(1.0, 0.0, 0.0),
            lambda: line_frequency_phase_slope(1.0, 0.0, -1.0),
            lambda: identifiable_rank([[float("nan")]]),
        ]
        for call in calls:
            with self.subTest(call=call):
                with self.assertRaises(HolonomyResponseError):
                    call()


if __name__ == "__main__":
    unittest.main()
