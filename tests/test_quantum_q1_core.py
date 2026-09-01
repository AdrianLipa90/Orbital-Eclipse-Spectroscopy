import unittest

import numpy as np

from oes.quantum.fermions import (
    build_sector_hamiltonian,
    determinant_basis,
    full_space_dimension,
    jw_product,
    one_rdm,
    sector_dimension,
    two_rdm,
)


class TestQuantumQ1Core(unittest.TestCase):
    def test_20_qubit_dimensions(self):
        self.assertEqual(full_space_dimension(20), 1_048_576)
        self.assertEqual(sector_dimension(20, 2), 190)
        self.assertEqual(len(determinant_basis(20, 2)), 190)

    def test_jordan_wigner_number_operator(self):
        expansion = jw_product(2, [("create", 0), ("annihilate", 0)])
        identity = ("I", "I")
        z0 = ("Z", "I")
        self.assertAlmostEqual(expansion[identity].real, 0.5, places=12)
        self.assertAlmostEqual(expansion[z0].real, -0.5, places=12)
        self.assertAlmostEqual(expansion[identity].imag, 0.0, places=12)
        self.assertAlmostEqual(expansion[z0].imag, 0.0, places=12)

    def test_noninteracting_two_electron_sector(self):
        h1 = np.diag([-1.0, 0.5])
        eri = np.zeros((2, 2, 2, 2))
        H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2)
        evals, evecs = np.linalg.eigh(H)
        self.assertEqual(H.shape, (6, 6))
        self.assertAlmostEqual(float(evals[0]), -2.0, places=12)

        gamma1 = one_rdm(evecs[:, 0], basis, 4)
        self.assertAlmostEqual(float(np.trace(gamma1).real), 2.0, places=12)

        gamma2 = two_rdm(evecs[:, 0], basis, 4)
        contraction = sum(gamma2[p, q, p, q] for p in range(4) for q in range(4))
        self.assertAlmostEqual(float(contraction.real), 2.0, places=12)


if __name__ == "__main__":
    unittest.main()
