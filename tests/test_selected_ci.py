import unittest

import numpy as np

from oes.quantum.fermions import build_sector_hamiltonian, determinant_basis
from oes.quantum.selected_ci import (
    build_two_electron_subspace_hamiltonian,
    grouped_importance_order,
    grouped_prefix_for_target,
)


def symmetric_eri(seed=7, n=3):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, n, n, n))
    perms = [
        x,
        x.transpose(1, 0, 2, 3),
        x.transpose(0, 1, 3, 2),
        x.transpose(1, 0, 3, 2),
        x.transpose(2, 3, 0, 1),
        x.transpose(3, 2, 0, 1),
        x.transpose(2, 3, 1, 0),
        x.transpose(3, 2, 1, 0),
    ]
    return sum(perms) / len(perms)


class SelectedCITests(unittest.TestCase):
    def test_arbitrary_two_electron_builder_matches_full_sector(self):
        rng = np.random.default_rng(11)
        a = rng.normal(size=(3, 3))
        h1 = 0.5 * (a + a.T)
        eri = symmetric_eri()
        reference, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=0.37)
        selected = build_two_electron_subspace_hamiltonian(h1, eri, basis, ecore=0.37)
        self.assertTrue(np.allclose(selected, reference, atol=1e-12))

    def test_selected_principal_subspace_matches_reference_principal_block(self):
        rng = np.random.default_rng(13)
        a = rng.normal(size=(3, 3))
        h1 = 0.5 * (a + a.T)
        eri = symmetric_eri(seed=17)
        reference, basis = build_sector_hamiltonian(h1, eri, n_electrons=2)
        keep = [0, 2, 5, 7, 11]
        sub_basis = [basis[i] for i in keep]
        selected = build_two_electron_subspace_hamiltonian(h1, eri, sub_basis)
        expected = reference[np.ix_(keep, keep)]
        self.assertTrue(np.allclose(selected, expected, atol=1e-12))

    def test_grouped_order_keeps_tied_scores_together(self):
        scores = np.array([0.9, 0.5, 0.5 * (1 + 1e-9), 0.1, 0.1, 0.01])
        groups = grouped_importance_order(scores, relative_tie_tolerance=1e-7)
        self.assertEqual(set(groups[0]), {0})
        self.assertEqual(set(groups[1]), {1, 2})
        self.assertEqual(set(groups[2]), {3, 4})
        prefix = grouped_prefix_for_target(groups, 2)
        self.assertEqual(len(prefix), 3)
        self.assertEqual(set(prefix), {0, 1, 2})


if __name__ == "__main__":
    unittest.main()
