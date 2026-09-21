"""Transition-density subspace continuity for geometry-dependent TDA states."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np


class TDATrackingError(ValueError):
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
        raise TDATrackingError(
            f"{name} must be a finite non-empty 2D matrix"
        )
    return array


def transition_density_inner_product(
    left_transition_density,
    right_transition_density,
    cross_metric,
) -> complex:
    """Hilbert-Schmidt overlap with nonorthogonal cross-basis metric.

    For AO transition matrices T_L and T_R,

        <T_L,T_R> = Tr[T_L^† S_LR T_R S_RL].
    """
    left = _matrix(
        left_transition_density,
        "left_transition_density",
    )
    right = _matrix(
        right_transition_density,
        "right_transition_density",
    )
    metric = _matrix(cross_metric, "cross_metric")
    if left.shape[0] != left.shape[1] or right.shape[0] != right.shape[1]:
        raise TDATrackingError("transition densities must be square")
    if metric.shape != (left.shape[0], right.shape[0]):
        raise TDATrackingError("cross_metric shape mismatch")
    return complex(
        np.trace(
            left.conj().T
            @ metric
            @ right
            @ metric.conj().T
        )
    )


def transition_density_gram(
    transition_densities: Iterable[np.ndarray],
    metric,
) -> np.ndarray:
    rows = tuple(
        _matrix(value, "transition_density")
        for value in transition_densities
    )
    if not rows:
        raise TDATrackingError(
            "transition-density block must not be empty"
        )
    dimension = rows[0].shape[0]
    if any(value.shape != (dimension, dimension) for value in rows):
        raise TDATrackingError(
            "transition-density block shapes must match"
        )
    local_metric = _matrix(metric, "metric")
    if local_metric.shape != (dimension, dimension):
        raise TDATrackingError("metric shape mismatch")

    gram = np.empty((len(rows), len(rows)), dtype=complex)
    for i, left in enumerate(rows):
        for j, right in enumerate(rows):
            gram[i, j] = transition_density_inner_product(
                left,
                right,
                local_metric,
            )
    return 0.5 * (gram + gram.conj().T)


def _inverse_sqrt_hermitian(
    gram: np.ndarray,
    rank_tol: float,
) -> np.ndarray:
    values, vectors = np.linalg.eigh(gram)
    if not np.all(np.isfinite(values)):
        raise TDATrackingError("transition-density Gram spectrum is invalid")
    if np.min(values) <= rank_tol:
        raise TDATrackingError(
            "transition-density block is rank deficient"
        )
    return (
        vectors
        * (1.0 / np.sqrt(values))[None, :]
    ) @ vectors.conj().T


@dataclass(frozen=True)
class TDABlockContinuity:
    principal_cosines: tuple[float, ...]
    minimum_principal_cosine: float
    chordal_distance: float
    cross_overlap_whitened: np.ndarray


def _root_indices(
    roots: Iterable[int],
    *,
    nstates: int,
    name: str,
) -> tuple[int, ...]:
    values = tuple(roots)
    if not values:
        raise TDATrackingError(f"{name} must not be empty")
    out = []
    for value in values:
        if isinstance(value, (bool, np.bool_)):
            raise TDATrackingError(f"{name} must contain integer root labels")
        integer = int(value)
        if integer != value or integer < 1 or integer > nstates:
            raise TDATrackingError(
                f"{name} root labels must lie in 1..{nstates}"
            )
        out.append(integer - 1)
    if len(set(out)) != len(out):
        raise TDATrackingError(f"{name} root labels must be unique")
    return tuple(out)


def transition_density_block_continuity(
    left_transition_densities: Iterable[np.ndarray],
    right_transition_densities: Iterable[np.ndarray],
    left_metric,
    right_metric,
    cross_metric,
    *,
    rank_tol: float = 1.0e-12,
    cosine_atol: float = 1.0e-8,
) -> TDABlockContinuity:
    """Compare transition-density spans, not arbitrary root vectors.

    Each block is Gram-whitened in its native AO metric before the cross-block
    singular values are evaluated. The result is invariant to nonsingular
    changes of basis inside either supplied transition-density span.
    """
    rank_tolerance = float(rank_tol)
    cosine_tolerance = float(cosine_atol)
    if (
        not math.isfinite(rank_tolerance)
        or rank_tolerance < 0.0
        or not math.isfinite(cosine_tolerance)
        or cosine_tolerance < 0.0
    ):
        raise TDATrackingError(
            "tracking tolerances must be finite and non-negative"
        )

    left = tuple(
        _matrix(value, "left_transition_density")
        for value in left_transition_densities
    )
    right = tuple(
        _matrix(value, "right_transition_density")
        for value in right_transition_densities
    )
    if not left or not right:
        raise TDATrackingError(
            "transition-density blocks must not be empty"
        )

    gram_left = transition_density_gram(left, left_metric)
    gram_right = transition_density_gram(right, right_metric)
    inverse_left = _inverse_sqrt_hermitian(
        gram_left,
        rank_tolerance,
    )
    inverse_right = _inverse_sqrt_hermitian(
        gram_right,
        rank_tolerance,
    )

    cross = np.empty((len(left), len(right)), dtype=complex)
    for i, left_density in enumerate(left):
        for j, right_density in enumerate(right):
            cross[i, j] = transition_density_inner_product(
                left_density,
                right_density,
                cross_metric,
            )

    whitened = inverse_left @ cross @ inverse_right
    singular_values = np.linalg.svd(
        whitened,
        compute_uv=False,
    )
    if singular_values.size and (
        np.max(singular_values) > 1.0 + cosine_tolerance
    ):
        raise TDATrackingError(
            "principal cosine exceeds one; check transition-density metric"
        )
    singular_values = np.clip(
        singular_values,
        0.0,
        1.0,
    )
    return TDABlockContinuity(
        principal_cosines=tuple(
            float(value) for value in singular_values
        ),
        minimum_principal_cosine=float(
            np.min(singular_values)
        ),
        chordal_distance=float(
            math.sqrt(
                np.sum(1.0 - singular_values * singular_values)
            )
        ),
        cross_overlap_whitened=whitened,
    )



def runtime_block_continuity(
    left_runtime,
    right_runtime,
    *,
    left_roots: Iterable[int],
    right_roots: Iterable[int],
    rank_tol: float = 1.0e-12,
    cosine_atol: float = 1.0e-8,
) -> TDABlockContinuity:
    """Compare selected TDA transition-density blocks across two geometries.

    Root labels are 1-based and are selectors only. Physical identity is the
    resulting transition-density span, not the root numbering.
    """
    from .molecular_tda import (
        MolecularTDARuntime,
        cross_geometry_ao_overlap,
    )

    if not isinstance(left_runtime, MolecularTDARuntime):
        raise TDATrackingError("left_runtime must be MolecularTDARuntime")
    if not isinstance(right_runtime, MolecularTDARuntime):
        raise TDATrackingError("right_runtime must be MolecularTDARuntime")

    left_index = _root_indices(
        left_roots,
        nstates=len(left_runtime.transition_densities_ao),
        name="left_roots",
    )
    right_index = _root_indices(
        right_roots,
        nstates=len(right_runtime.transition_densities_ao),
        name="right_roots",
    )
    if len(left_index) != len(right_index):
        raise TDATrackingError(
            "left/right root blocks must have equal dimension"
        )

    cross_metric = cross_geometry_ao_overlap(
        left_runtime,
        right_runtime,
    )
    return transition_density_block_continuity(
        tuple(
            left_runtime.transition_densities_ao[index]
            for index in left_index
        ),
        tuple(
            right_runtime.transition_densities_ao[index]
            for index in right_index
        ),
        left_runtime.ao_overlap,
        right_runtime.ao_overlap,
        cross_metric,
        rank_tol=rank_tol,
        cosine_atol=cosine_atol,
    )


__all__ = [
    "TDATrackingError",
    "TDABlockContinuity",
    "transition_density_inner_product",
    "transition_density_gram",
    "transition_density_block_continuity",
    "runtime_block_continuity",
]
