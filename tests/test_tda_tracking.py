import math
import unittest

import numpy as np

from oes.tda_tracking import (
    TDATrackingError,
    transition_density_block_continuity,
)


class TDATrackingTests(unittest.TestCase):
    def test_same_span_is_invariant_to_internal_mixing(self):
        e00 = np.asarray([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
        e01 = np.asarray([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
        c, s = math.cos(0.4), math.sin(0.4)
        right = (
            c * e00 + s * e01,
            -s * e00 + c * e01,
        )
        result = transition_density_block_continuity(
            (e00, e01),
            right,
            np.eye(2),
            np.eye(2),
            np.eye(2),
        )
        np.testing.assert_allclose(
            result.principal_cosines,
            (1.0, 1.0),
            rtol=0.0,
            atol=1.0e-12,
        )
        self.assertAlmostEqual(
            result.chordal_distance,
            0.0,
            places=12,
        )

    def test_known_principal_angle_is_recovered(self):
        e00 = np.asarray([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
        e01 = np.asarray([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
        e10 = np.asarray([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
        theta = 0.5
        right = (
            e00,
            math.cos(theta) * e01 + math.sin(theta) * e10,
        )
        result = transition_density_block_continuity(
            (e00, e01),
            right,
            np.eye(2),
            np.eye(2),
            np.eye(2),
        )
        np.testing.assert_allclose(
            result.principal_cosines,
            (1.0, math.cos(theta)),
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_gram_whitening_removes_arbitrary_rescaling(self):
        a = np.asarray([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
        b = np.asarray([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
        reference = transition_density_block_continuity(
            (a, b),
            (a, b),
            np.eye(2),
            np.eye(2),
            np.eye(2),
        )
        rescaled = transition_density_block_continuity(
            (2.0 * a, -3.0j * b),
            (0.4 * a, 5.0 * b),
            np.eye(2),
            np.eye(2),
            np.eye(2),
        )
        np.testing.assert_allclose(
            rescaled.principal_cosines,
            reference.principal_cosines,
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_rank_deficient_block_fails_closed(self):
        a = np.asarray([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
        with self.assertRaises(TDATrackingError):
            transition_density_block_continuity(
                (a, 2.0 * a),
                (a, 3.0 * a),
                np.eye(2),
                np.eye(2),
                np.eye(2),
            )


if __name__ == "__main__":
    unittest.main()
