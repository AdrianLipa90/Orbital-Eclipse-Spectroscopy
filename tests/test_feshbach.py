import unittest

import numpy as np

from oes.quantum.feshbach import effective_hamiltonian, eigenpair_downfolding_residual


class FeshbachTests(unittest.TestCase):
    def test_exact_downfolding_matches_full_eigenpair(self):
        H = np.array(
            [
                [-1.0, 0.08, 0.20],
                [0.08, -0.7, -0.10],
                [0.20, -0.10, 0.4],
            ],
            dtype=float,
        )
        vals, vecs = np.linalg.eigh(H)
        for i, energy in enumerate(vals):
            out = eigenpair_downfolding_residual(H, 2, float(energy), vecs[:, i])
            self.assertLess(out["effective_eigen_residual_hartree"], 1e-11)
            self.assertLess(out["q_reconstruction_error"], 1e-11)

    def test_effective_hamiltonian_is_symmetric(self):
        H = np.array(
            [
                [-1.0, 0.1, 0.2, 0.0],
                [0.1, -0.8, 0.0, 0.15],
                [0.2, 0.0, 0.3, 0.04],
                [0.0, 0.15, 0.04, 0.6],
            ],
            dtype=float,
        )
        heff, diag = effective_hamiltonian(H, 2, -0.9)
        self.assertTrue(np.allclose(heff, heff.T, atol=1e-13))
        self.assertEqual(diag["p_dimension"], 2)
        self.assertEqual(diag["q_dimension"], 2)

    def test_resolvent_fails_closed_at_q_pole(self):
        H = np.array([[-1.0, 0.1], [0.1, 0.5]], dtype=float)
        with self.assertRaisesRegex(RuntimeError, "resolvent gate failed"):
            effective_hamiltonian(H, 1, 0.5, singular_floor=1e-8)


if __name__ == "__main__":
    unittest.main()
