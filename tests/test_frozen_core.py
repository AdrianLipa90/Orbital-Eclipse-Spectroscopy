import unittest

import numpy as np

from oes.quantum.fermions import build_sector_hamiltonian, determinant_basis
from oes.quantum.frozen_core import frozen_core_effective_hamiltonian


def symmetric_chemist_eri(seed: int, n: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(n, n, n, n))
    eri = 0.25 * (
        raw
        + raw.transpose(1, 0, 2, 3)
        + raw.transpose(0, 1, 3, 2)
        + raw.transpose(1, 0, 3, 2)
    )
    eri = 0.5 * (eri + eri.transpose(2, 3, 0, 1))
    return eri


class FrozenCoreTests(unittest.TestCase):
    def test_effective_hamiltonian_equals_exact_frozen_subspace_projection(self):
        rng = np.random.default_rng(190)
        n_full = 4
        h = rng.normal(size=(n_full, n_full))
        h = 0.5 * (h + h.T)
        eri = symmetric_chemist_eri(191, n_full)
        nuclear = 0.3141592653589793

        full_H, full_basis = build_sector_hamiltonian(h, eri, n_electrons=4, ecore=nuclear)
        full_index = {det: i for i, det in enumerate(full_basis)}

        reduced = frozen_core_effective_hamiltonian(
            h,
            eri,
            core_indices=(0,),
            active_indices=(1, 2, 3),
            nuclear_energy=nuclear,
        )
        active_H, active_basis = build_sector_hamiltonian(
            reduced.h1_active,
            reduced.eri_active,
            n_electrons=2,
            ecore=reduced.ecore,
        )

        # Spatial orbital 0 is the frozen doubly occupied core, i.e. full spin
        # modes 0 and 1. Active spin modes map contiguously to full modes 2..7.
        mapped = [(1 << 0) | (1 << 1) | (int(det) << 2) for det in active_basis]
        rows = [full_index[det] for det in mapped]
        projected = full_H[np.ix_(rows, rows)]

        self.assertEqual(active_H.shape, (15, 15))
        self.assertTrue(np.allclose(active_H, projected, atol=2e-12, rtol=2e-12))

    def test_metadata_and_validation(self):
        h = np.eye(3)
        eri = np.zeros((3, 3, 3, 3))
        out = frozen_core_effective_hamiltonian(h, eri, (0,), (1, 2), nuclear_energy=0.5)
        self.assertEqual(out.metadata()["n_frozen_electrons"], 2)
        self.assertEqual(out.metadata()["n_active_spatial"], 2)
        self.assertAlmostEqual(out.ecore, 2.5)

        with self.assertRaises(ValueError):
            frozen_core_effective_hamiltonian(h, eri, (0,), (0, 1))
        with self.assertRaises(ValueError):
            frozen_core_effective_hamiltonian(h, eri, (), ())


if __name__ == "__main__":
    unittest.main()
