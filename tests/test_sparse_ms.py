import unittest

import numpy as np

from oes.quantum.determinant_subspace import build_determinant_subspace_hamiltonian
from oes.quantum.sparse_ms import build_sparse_fixed_spin_hamiltonian, fixed_spin_determinant_basis


def symmetric_chemist_eri(seed: int, n: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(n, n, n, n))
    eri = 0.25 * (
        raw
        + raw.transpose(1, 0, 2, 3)
        + raw.transpose(0, 1, 3, 2)
        + raw.transpose(1, 0, 3, 2)
    )
    return 0.5 * (eri + eri.transpose(2, 3, 0, 1))


class SparseFixedSpinTests(unittest.TestCase):
    def _compare(self, n_spatial: int, n_alpha: int, n_beta: int, seed: int):
        rng = np.random.default_rng(seed)
        h = rng.normal(size=(n_spatial, n_spatial))
        h = 0.5 * (h + h.T)
        eri = symmetric_chemist_eri(seed + 1, n_spatial)
        basis = fixed_spin_determinant_basis(n_spatial, n_alpha, n_beta)
        dense, dense_basis = build_determinant_subspace_hamiltonian(h, eri, basis, ecore=0.2718281828)
        sparse, sparse_basis = build_sparse_fixed_spin_hamiltonian(
            h, eri, n_alpha, n_beta, ecore=0.2718281828, zero_tolerance=0.0
        )
        self.assertEqual(dense_basis, sparse_basis)
        self.assertTrue(np.allclose(sparse.toarray(), dense, atol=3e-12, rtol=3e-12))
        self.assertTrue(np.allclose(sparse.toarray(), sparse.toarray().T, atol=1e-13))

    def test_two_electron_ms_zero(self):
        self._compare(4, 1, 1, 210)

    def test_three_electron_doublet_sector(self):
        self._compare(4, 2, 1, 220)

    def test_four_electron_ms_zero(self):
        self._compare(5, 2, 2, 230)

    def test_six_electron_dimension_at_20q(self):
        basis = fixed_spin_determinant_basis(10, 3, 3)
        self.assertEqual(len(basis), 14400)

    def test_eight_electron_dimension_at_20q(self):
        basis = fixed_spin_determinant_basis(10, 4, 4)
        self.assertEqual(len(basis), 44100)


if __name__ == "__main__":
    unittest.main()
