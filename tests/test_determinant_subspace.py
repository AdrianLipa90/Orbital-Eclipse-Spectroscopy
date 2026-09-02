import unittest

import numpy as np

from oes.quantum.determinant_subspace import build_determinant_subspace_hamiltonian
from oes.quantum.fermions import build_sector_hamiltonian, determinant_basis
from oes.quantum.lithium_a1 import spin_sector_indices


class DeterminantSubspaceTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(712)
        n = 4
        h = rng.normal(size=(n, n))
        self.h1 = 0.5 * (h + h.T)
        eri = rng.normal(size=(n, n, n, n))
        # Project a random tensor onto the usual chemists' integral symmetries.
        eri = 0.25 * (
            eri
            + eri.transpose(1, 0, 2, 3)
            + eri.transpose(0, 1, 3, 2)
            + eri.transpose(2, 3, 0, 1)
        )
        self.eri = eri

    def _compare_spin_sector(self, n_electrons, n_alpha, n_beta):
        H_full, basis_full = build_sector_hamiltonian(self.h1, self.eri, n_electrons, ecore=0.17)
        indices = spin_sector_indices(basis_full, 4, n_alpha, n_beta)
        subset = [basis_full[i] for i in indices]
        H_sub, basis_sub = build_determinant_subspace_hamiltonian(
            self.h1, self.eri, subset, ecore=0.17
        )
        expected = H_full[np.ix_(indices, indices)]
        self.assertEqual(tuple(subset), basis_sub)
        self.assertTrue(np.allclose(H_sub, expected, atol=1e-12))

    def test_two_electron_ms_zero_matches_full_sector_slice(self):
        self._compare_spin_sector(2, 1, 1)

    def test_three_electron_ms_half_matches_full_sector_slice(self):
        self._compare_spin_sector(3, 2, 1)

    def test_four_electron_ms_zero_matches_full_sector_slice(self):
        self._compare_spin_sector(4, 2, 2)

    def test_rejects_mixed_particle_number(self):
        with self.assertRaises(ValueError):
            build_determinant_subspace_hamiltonian(self.h1, self.eri, [1, 3])


if __name__ == "__main__":
    unittest.main()
