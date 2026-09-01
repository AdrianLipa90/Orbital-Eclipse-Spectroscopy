import unittest

import numpy as np

from oes.quantum.external_dressing import ExternalCouplingSpace, en2_correction


class ExternalDressingTests(unittest.TestCase):
    def test_en2_matches_closed_form_toy_model(self):
        # One active state |P> coupled to two external determinants.
        # E_P = -1.0, H_aa = {-0.5, 0.0}, V = {0.2, 0.3} Hartree.
        # EN2 = 0.2^2/(-0.5) + 0.3^2/(-1.0) = -0.17 Hartree.
        external = ExternalCouplingSpace(
            external_basis=(1, 2),
            diagonal_hartree=np.array([-0.5, 0.0], dtype=float),
            coupling_qp=np.array([[0.2], [0.3]], dtype=float),
            n_full_spin_orbitals=3,
            n_active_spin_orbitals=1,
        )
        out = en2_correction(-1.0, np.array([1.0]), external)
        self.assertAlmostEqual(out["correction_hartree"], -0.17, places=14)
        self.assertAlmostEqual(out["coupling_norm2_hartree2"], 0.13, places=14)
        self.assertAlmostEqual(out["min_abs_denominator_hartree"], 0.5, places=14)
        self.assertEqual(out["external_determinants"], 2)

    def test_en2_fails_closed_on_intruder_denominator(self):
        external = ExternalCouplingSpace(
            external_basis=(1,),
            diagonal_hartree=np.array([-1.0 + 5e-7], dtype=float),
            coupling_qp=np.array([[0.1]], dtype=float),
            n_full_spin_orbitals=2,
            n_active_spin_orbitals=1,
        )
        with self.assertRaisesRegex(RuntimeError, "intruder gate failed"):
            en2_correction(
                -1.0,
                np.array([1.0]),
                external,
                denominator_floor=1e-5,
            )

    def test_en2_rejects_incompatible_state_dimension(self):
        external = ExternalCouplingSpace(
            external_basis=(1,),
            diagonal_hartree=np.array([0.0], dtype=float),
            coupling_qp=np.array([[0.1, 0.2]], dtype=float),
            n_full_spin_orbitals=3,
            n_active_spin_orbitals=2,
        )
        with self.assertRaises(ValueError):
            en2_correction(-1.0, np.array([1.0]), external)


if __name__ == "__main__":
    unittest.main()
