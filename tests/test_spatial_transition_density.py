import math
import unittest

import numpy as np

from oes.spatial_transition_density import (
    SpatialTransitionError,
    block_spatial_invariants,
)


class SpatialTransitionInvariantTests(unittest.TestCase):
    def test_unitary_root_mixing_preserves_block_invariants(self):
        coords = np.asarray(
            [
                [-1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        weights = np.asarray([0.5, 0.75, 0.5, 0.4])
        rho = np.asarray(
            [
                [0.8, -0.2, -0.5, 0.1],
                [0.1j, 0.3, -0.2j, 0.4],
            ],
            dtype=complex,
        )
        dipoles = np.einsum(
            "x,nx,xk->nk",
            weights,
            rho,
            coords,
        )

        theta = 0.41
        phase = np.exp(0.37j)
        unitary = np.asarray(
            [
                [math.cos(theta), phase * math.sin(theta)],
                [-phase.conjugate() * math.sin(theta), math.cos(theta)],
            ],
            dtype=complex,
        )
        mixed_rho = unitary @ rho
        mixed_dipoles = unitary @ dipoles

        reference = block_spatial_invariants(
            rho,
            weights,
            coords,
            dipoles,
        )
        observed = block_spatial_invariants(
            mixed_rho,
            weights,
            coords,
            mixed_dipoles,
        )

        scalar_keys = (
            "transition_density_l2_activity_bohr_minus3",
            "dipole_envelope_au",
            "net_block_dipole_norm_au",
            "dipole_survival_ratio",
            "dipole_cancellation_fraction",
            "activity_rms_radius_bohr",
        )
        for key in scalar_keys:
            self.assertAlmostEqual(
                reference[key],
                observed[key],
                places=12,
            )
        np.testing.assert_allclose(
            reference["activity_centroid_bohr"],
            observed["activity_centroid_bohr"],
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_perfect_two_point_cancellation_has_zero_survival(self):
        coords = np.asarray(
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        )
        weights = np.ones(2)
        rho = np.asarray([[1.0, 1.0]], dtype=complex)
        dipoles = np.einsum(
            "x,nx,xk->nk",
            weights,
            rho,
            coords,
        )
        result = block_spatial_invariants(
            rho,
            weights,
            coords,
            dipoles,
        )
        self.assertAlmostEqual(
            result["net_block_dipole_norm_au"],
            0.0,
            places=15,
        )
        self.assertAlmostEqual(
            result["dipole_survival_ratio"],
            0.0,
            places=15,
        )
        self.assertAlmostEqual(
            result["dipole_cancellation_fraction"],
            1.0,
            places=15,
        )

    def test_aligned_two_point_contribution_has_unit_survival(self):
        coords = np.asarray(
            [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        )
        weights = np.ones(2)
        rho = np.asarray([[1.0, -1.0]], dtype=complex)
        dipoles = np.einsum(
            "x,nx,xk->nk",
            weights,
            rho,
            coords,
        )
        result = block_spatial_invariants(
            rho,
            weights,
            coords,
            dipoles,
        )
        self.assertAlmostEqual(
            result["dipole_survival_ratio"],
            1.0,
            places=15,
        )
        self.assertAlmostEqual(
            result["dipole_cancellation_fraction"],
            0.0,
            places=15,
        )

    def test_negative_weights_fail_closed(self):
        with self.assertRaises(SpatialTransitionError):
            block_spatial_invariants(
                np.ones((1, 2)),
                np.asarray([1.0, -1.0]),
                np.zeros((2, 3)),
                np.zeros((1, 3)),
            )


if __name__ == "__main__":
    unittest.main()
