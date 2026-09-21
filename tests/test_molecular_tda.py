import importlib.util
import math
import unittest

import numpy as np

from oes.molecular_tda import (
    MolecularTDAError,
    run_closed_shell_tda,
    run_orientation_tda_scan,
)

HAS_PYSCF = importlib.util.find_spec("pyscf") is not None


def _water_geometry(angle_rad: float, bond_angstrom: float = 0.96):
    half = 0.5 * angle_rad
    return (
        (0.0, 0.0, 0.0),
        (
            bond_angstrom * math.sin(half),
            0.0,
            bond_angstrom * math.cos(half),
        ),
        (
            -bond_angstrom * math.sin(half),
            0.0,
            bond_angstrom * math.cos(half),
        ),
    )


class MolecularTDAValidationTests(unittest.TestCase):
    def test_fail_closed_before_backend_for_open_shell_request(self):
        with self.assertRaises(MolecularTDAError):
            run_closed_shell_tda(
                species_label="open-shell-rejected",
                atoms=("H",),
                coordinates_angstrom=((0.0, 0.0, 0.0),),
                basis_name="sto-3g",
                spin=1,
                nstates=1,
            )

    def test_fail_closed_on_bad_geometry_shape(self):
        with self.assertRaises(MolecularTDAError):
            run_closed_shell_tda(
                species_label="bad-geometry",
                atoms=("H", "H"),
                coordinates_angstrom=((0.0, 0.0, 0.0),),
                basis_name="sto-3g",
                nstates=1,
            )


@unittest.skipUnless(HAS_PYSCF, "PySCF q1 extra is not installed")
class MolecularTDAPySCFTests(unittest.TestCase):
    def test_global_rotation_leaves_h2_scalar_spectrum_invariant(self):
        z_geometry = (
            (0.0, 0.0, -0.37),
            (0.0, 0.0, 0.37),
        )
        x_geometry = (
            (-0.37, 0.0, 0.0),
            (0.37, 0.0, 0.0),
        )
        z = run_closed_shell_tda(
            species_label="H2-z",
            atoms=("H", "H"),
            coordinates_angstrom=z_geometry,
            basis_name="sto-3g",
            nstates=1,
        )
        x = run_closed_shell_tda(
            species_label="H2-x",
            atoms=("H", "H"),
            coordinates_angstrom=x_geometry,
            basis_name="sto-3g",
            nstates=1,
        )
        self.assertAlmostEqual(
            z.states[0].excitation_hartree,
            x.states[0].excitation_hartree,
            places=10,
        )
        self.assertAlmostEqual(
            z.states[0].transition_dipole_norm_au,
            x.states[0].transition_dipole_norm_au,
            places=10,
        )
        self.assertAlmostEqual(
            z.states[0].oscillator_strength_oes,
            x.states[0].oscillator_strength_oes,
            places=10,
        )
        self.assertLess(z.max_oscillator_strength_delta, 1.0e-10)
        self.assertLess(x.max_oscillator_strength_delta, 1.0e-10)

    def test_internal_water_bend_changes_raw_tda_spectrum(self):
        angles = (math.radians(90.0), math.radians(120.0))
        scan = run_orientation_tda_scan(
            species_label="H2O-rigid-bend-control",
            atoms=("O", "H", "H"),
            samples=tuple(
                (theta, _water_geometry(theta))
                for theta in angles
            ),
            basis_name="sto-3g",
            nstates=2,
        )
        self.assertEqual(
            scan.state_tracking_status,
            "RAW_ROOT_INDEX_UNTRACKED",
        )
        self.assertEqual(
            tuple(point.angle_rad for point in scan.points),
            angles,
        )
        left = np.asarray(
            [
                state.excitation_hartree
                for state in scan.points[0].result.states
            ]
        )
        right = np.asarray(
            [
                state.excitation_hartree
                for state in scan.points[1].result.states
            ]
        )
        self.assertEqual(left.shape, right.shape)
        self.assertGreater(
            float(np.max(np.abs(left - right))),
            1.0e-6,
        )
        for point in scan.points:
            self.assertLess(
                point.result.max_oscillator_strength_delta,
                1.0e-10,
            )


if __name__ == "__main__":
    unittest.main()
