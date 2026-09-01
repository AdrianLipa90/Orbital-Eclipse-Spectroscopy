import unittest

import numpy as np

from oes.quantum.fermions import determinant_basis
from oes.quantum.h2_m1 import (
    _orbital_inversion_parities,
    _select_complete_energy_blocks,
    determinant_inversion_parity,
    state_inversion_parity,
)


class H2M1CoreTests(unittest.TestCase):
    def test_complete_block_selector_never_cuts_degeneracy(self):
        energies = np.array([
            -1.0,
            -0.5,
            0.1, 0.1,
            0.3, 0.3,
            0.5,
            0.7,
            0.9, 0.9,
            1.2, 1.2, 1.2,
        ])
        selected, sizes, _ = _select_complete_energy_blocks(energies, target_orbitals=10)
        self.assertEqual(len(selected), 10)
        self.assertEqual(sum(sizes), 10)
        for pair in ((2, 3), (4, 5), (8, 9)):
            self.assertEqual(pair[0] in selected, pair[1] in selected)

    def test_orbital_g_u_labels_map_to_parity(self):
        p = _orbital_inversion_parities(["Ag", "B1u", "B2g", "B3u"])
        self.assertTrue(np.array_equal(p, np.array([1, -1, 1, -1])))

    def test_determinant_and_state_parity(self):
        spatial = np.array([1, -1])
        basis = determinant_basis(4, 2)
        values = [determinant_inversion_parity(det, 4, spatial) for det in basis]
        self.assertTrue(set(values).issubset({-1, 1}))
        # Pure determinant occupying spatial orbital 0 with alpha+beta is g.
        target = (1 << 0) | (1 << 1)
        idx = basis.index(target)
        state = np.zeros(len(basis))
        state[idx] = 1.0
        self.assertAlmostEqual(state_inversion_parity(state, basis, spatial), 1.0, places=14)


if __name__ == "__main__":
    unittest.main()
