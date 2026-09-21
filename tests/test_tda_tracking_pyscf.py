import importlib.util
import math
import unittest

import numpy as np

from oes.molecular_tda import (
    run_closed_shell_tda,
    run_closed_shell_tda_runtime,
)
from oes.tda_tracking import (
    TDATrackingError,
    runtime_block_continuity,
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


@unittest.skipUnless(HAS_PYSCF, "PySCF q1 extra is not installed")
class TDATrackingPySCFTests(unittest.TestCase):
    def test_runtime_wrapper_preserves_public_result(self):
        kwargs = dict(
            species_label="H2-runtime",
            atoms=("H", "H"),
            coordinates_angstrom=((0.0, 0.0, -0.37), (0.0, 0.0, 0.37)),
            basis_name="sto-3g",
            nstates=1,
        )
        runtime = run_closed_shell_tda_runtime(**kwargs)
        public = run_closed_shell_tda(**kwargs)
        self.assertAlmostEqual(
            runtime.result.rhf_energy_hartree,
            public.rhf_energy_hartree,
            places=12,
        )
        self.assertAlmostEqual(
            runtime.result.states[0].excitation_hartree,
            public.states[0].excitation_hartree,
            places=12,
        )
        self.assertEqual(len(runtime.transition_densities_ao), 1)
        self.assertEqual(
            runtime.transition_densities_ao[0].shape,
            runtime.ao_overlap.shape,
        )
        self.assertTrue(
            np.all(np.isfinite(runtime.transition_densities_ao[0]))
        )

    def test_identical_geometry_has_unit_transition_density_continuity(self):
        kwargs = dict(
            species_label="H2-identical",
            atoms=("H", "H"),
            coordinates_angstrom=((0.0, 0.0, -0.37), (0.0, 0.0, 0.37)),
            basis_name="sto-3g",
            nstates=1,
        )
        left = run_closed_shell_tda_runtime(**kwargs)
        right = run_closed_shell_tda_runtime(**kwargs)
        continuity = runtime_block_continuity(
            left,
            right,
            left_roots=(1,),
            right_roots=(1,),
        )
        self.assertAlmostEqual(
            continuity.minimum_principal_cosine,
            1.0,
            places=10,
        )
        self.assertAlmostEqual(
            continuity.chordal_distance,
            0.0,
            places=10,
        )

    def test_water_bend_two_root_block_has_finite_continuity(self):
        left = run_closed_shell_tda_runtime(
            species_label="H2O-90",
            atoms=("O", "H", "H"),
            coordinates_angstrom=_water_geometry(math.radians(90.0)),
            basis_name="sto-3g",
            nstates=2,
        )
        right = run_closed_shell_tda_runtime(
            species_label="H2O-120",
            atoms=("O", "H", "H"),
            coordinates_angstrom=_water_geometry(math.radians(120.0)),
            basis_name="sto-3g",
            nstates=2,
        )
        continuity = runtime_block_continuity(
            left,
            right,
            left_roots=(1, 2),
            right_roots=(1, 2),
        )
        self.assertEqual(len(continuity.principal_cosines), 2)
        self.assertTrue(
            all(
                0.0 <= value <= 1.0
                for value in continuity.principal_cosines
            )
        )
        self.assertTrue(math.isfinite(continuity.chordal_distance))

    def test_runtime_root_selector_fails_closed(self):
        runtime = run_closed_shell_tda_runtime(
            species_label="H2-selector",
            atoms=("H", "H"),
            coordinates_angstrom=((0.0, 0.0, -0.37), (0.0, 0.0, 0.37)),
            basis_name="sto-3g",
            nstates=1,
        )
        with self.assertRaises(TDATrackingError):
            runtime_block_continuity(
                runtime,
                runtime,
                left_roots=(1,),
                right_roots=(2,),
            )


if __name__ == "__main__":
    unittest.main()
