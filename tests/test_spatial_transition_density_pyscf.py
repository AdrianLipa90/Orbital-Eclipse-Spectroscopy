import importlib.util
import unittest

from oes.molecular_tda import run_closed_shell_tda_runtime
from oes.spatial_transition_density import (
    evaluate_runtime_block_spatial,
)

HAS_PYSCF = importlib.util.find_spec("pyscf") is not None


@unittest.skipUnless(HAS_PYSCF, "PySCF q1 extra is not installed")
class SpatialTransitionPySCFTests(unittest.TestCase):
    def test_runtime_tdm_reconstructs_backend_dipole_exactly(self):
        runtime = run_closed_shell_tda_runtime(
            species_label="H2-tdm-exact",
            atoms=("H", "H"),
            coordinates_angstrom=((0.0, 0.0, -0.37), (0.0, 0.0, 0.37)),
            basis_name="sto-3g",
            nstates=1,
        )
        self.assertEqual(
            len(runtime.transition_dipole_reconstruction_deltas),
            1,
        )
        self.assertLess(
            runtime.transition_dipole_reconstruction_deltas[0],
            1.0e-10,
        )

    def test_h2_real_space_grid_reconstructs_transition_dipole(self):
        runtime = run_closed_shell_tda_runtime(
            species_label="H2-grid",
            atoms=("H", "H"),
            coordinates_angstrom=((0.0, 0.0, -0.37), (0.0, 0.0, 0.37)),
            basis_name="sto-3g",
            nstates=1,
        )
        spatial = evaluate_runtime_block_spatial(
            runtime,
            (1,),
            grid_level=3,
            dipole_grid_atol=5.0e-5,
            dipole_grid_rtol=5.0e-4,
        )
        self.assertGreater(spatial.n_grid_points, 0)
        self.assertLess(spatial.max_grid_dipole_delta_au, 5.0e-5)
        self.assertGreaterEqual(spatial.dipole_survival_ratio, 0.0)
        self.assertLessEqual(spatial.dipole_survival_ratio, 1.0)
        self.assertGreater(
            spatial.transition_density_l2_activity_bohr_minus3,
            0.0,
        )
        self.assertGreater(spatial.dipole_envelope_au, 0.0)


if __name__ == "__main__":
    unittest.main()
