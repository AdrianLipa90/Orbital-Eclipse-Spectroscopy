"""Gauge-clean subspace tracking for geometry-dependent OES states."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class SubspaceTrackingError(ValueError):
    pass


def _matrix(values, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=complex)
    if (
        array.ndim != 2
        or array.shape[0] < 1
        or array.shape[1] < 1
        or not np.all(np.isfinite(array.real))
        or not np.all(np.isfinite(array.imag))
    ):
        raise SubspaceTrackingError(
            f"{name} must be a finite non-empty 2D matrix"
        )
    return array


def cross_subspace_overlap(
    left_coefficients,
    right_coefficients,
    cross_metric,
) -> np.ndarray:
    """Return C_L^† S_LR C_R for equal-rank subspaces."""
    left = _matrix(left_coefficients, "left_coefficients")
    right = _matrix(right_coefficients, "right_coefficients")
    metric = _matrix(cross_metric, "cross_metric")
    if left.shape[1] != right.shape[1]:
        raise SubspaceTrackingError(
            "left/right subspaces must have equal column count"
        )
    if metric.shape != (left.shape[0], right.shape[0]):
        raise SubspaceTrackingError(
            "cross_metric shape must map right basis into left basis"
        )
    return left.conj().T @ metric @ right


def principal_cosines(
    left_coefficients,
    right_coefficients,
    cross_metric,
    *,
    atol: float = 1.0e-8,
) -> np.ndarray:
    """Return singular values of the cross-subspace overlap matrix."""
    tol = float(atol)
    if not np.isfinite(tol) or tol < 0.0:
        raise SubspaceTrackingError("atol must be finite and non-negative")
    overlap = cross_subspace_overlap(
        left_coefficients,
        right_coefficients,
        cross_metric,
    )
    singular_values = np.linalg.svd(overlap, compute_uv=False)
    if singular_values.size and np.max(singular_values) > 1.0 + tol:
        raise SubspaceTrackingError(
            "principal cosine exceeds one; check metric/orthonormality"
        )
    return np.clip(singular_values, 0.0, 1.0)


@dataclass(frozen=True)
class SubspaceAlignment:
    principal_cosines: np.ndarray
    rotation: np.ndarray
    aligned_right_coefficients: np.ndarray
    overlap_before: np.ndarray
    overlap_after: np.ndarray


def align_right_subspace(
    left_coefficients,
    right_coefficients,
    cross_metric,
    *,
    atol: float = 1.0e-8,
) -> SubspaceAlignment:
    """Align the right subspace to the left by metric orthogonal Procrustes.

    If M = C_L^† S_LR C_R = U Sigma V^†, use Q = V U^† and return C_R Q.
    This changes only the basis inside the right subspace.
    """
    tol = float(atol)
    if not np.isfinite(tol) or tol < 0.0:
        raise SubspaceTrackingError("atol must be finite and non-negative")

    left = _matrix(left_coefficients, "left_coefficients")
    right = _matrix(right_coefficients, "right_coefficients")
    metric = _matrix(cross_metric, "cross_metric")
    overlap = cross_subspace_overlap(left, right, metric)
    u, singular_values, vh = np.linalg.svd(overlap, full_matrices=False)
    if singular_values.size and np.max(singular_values) > 1.0 + tol:
        raise SubspaceTrackingError(
            "principal cosine exceeds one; check metric/orthonormality"
        )

    rotation = vh.conj().T @ u.conj().T
    aligned = right @ rotation
    after = left.conj().T @ metric @ aligned
    return SubspaceAlignment(
        principal_cosines=np.clip(singular_values, 0.0, 1.0),
        rotation=rotation,
        aligned_right_coefficients=aligned,
        overlap_before=overlap,
        overlap_after=after,
    )


__all__ = [
    "SubspaceTrackingError",
    "SubspaceAlignment",
    "cross_subspace_overlap",
    "principal_cosines",
    "align_right_subspace",
]
