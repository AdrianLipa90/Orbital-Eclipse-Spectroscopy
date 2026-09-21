import math
import unittest

import numpy as np

from oes.subspace_tracking import (
    SubspaceTrackingError,
    align_right_subspace,
    principal_cosines,
)


class SubspaceTrackingTests(unittest.TestCase):
    def test_procrustes_removes_internal_unitary_rotation(self):
        left = np.eye(3, dtype=complex)[:, :2]
        phase = np.diag([np.exp(0.3j), np.exp(-0.7j)])
        c, s = math.cos(0.4), math.sin(0.4)
        rotation = np.asarray([[c, -s], [s, c]], dtype=complex)
        right = left @ phase @ rotation

        result = align_right_subspace(left, right, np.eye(3))
        np.testing.assert_allclose(
            result.principal_cosines,
            np.ones(2),
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            result.aligned_right_coefficients,
            left,
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_principal_cosines_are_invariant_to_subspace_gauge(self):
        left = np.eye(3, dtype=complex)[:, :2]
        right = np.column_stack(
            [
                np.asarray([1.0, 0.0, 0.0]),
                np.asarray([0.0, math.cos(0.5), math.sin(0.5)]),
            ]
        ).astype(complex)
        baseline = principal_cosines(left, right, np.eye(3))

        q_left = np.asarray([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
        q_right = np.diag([np.exp(0.2j), np.exp(-0.1j)])
        transformed = principal_cosines(
            left @ q_left,
            right @ q_right,
            np.eye(3),
        )
        np.testing.assert_allclose(
            transformed,
            baseline,
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_fail_closed_on_invalid_metric_shape(self):
        with self.assertRaises(SubspaceTrackingError):
            principal_cosines(
                np.eye(3)[:, :2],
                np.eye(4)[:, :2],
                np.eye(3),
            )


if __name__ == "__main__":
    unittest.main()
