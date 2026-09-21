"""Block-resolved orientation spectroscopy for tracked TDA sectors."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

from .molecular_tda import MolecularTDAResult, MolecularTDARuntime
from .tda_tracking import TDABlockContinuity, runtime_block_continuity


class BlockResponseError(ValueError):
    pass


def _root_indices(
    roots: Iterable[int],
    *,
    nstates: int,
) -> tuple[int, ...]:
    values = tuple(roots)
    if not values:
        raise BlockResponseError("root block must not be empty")
    indices = []
    for value in values:
        if isinstance(value, bool):
            raise BlockResponseError("root labels must be integers")
        integer = int(value)
        if integer != value or integer < 1 or integer > nstates:
            raise BlockResponseError(
                f"root labels must lie in 1..{nstates}"
            )
        indices.append(integer - 1)
    if len(set(indices)) != len(indices):
        raise BlockResponseError("root labels must be unique")
    return tuple(indices)


@dataclass(frozen=True)
class BlockSpectralObservables:
    roots: tuple[int, ...]
    excitation_min_ev: float
    excitation_max_ev: float
    excitation_centroid_ev: float
    oscillator_weighted_centroid_ev: float | None
    oscillator_strength_sum: float
    dipole_strength_sum_au2: float
    oscillator_strength_fraction_of_computed_window: float

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class TrackedBlockResponse:
    left: BlockSpectralObservables
    right: BlockSpectralObservables
    continuity: TDABlockContinuity
    oscillator_strength_delta: float
    oscillator_strength_ratio: float | None
    centroid_shift_ev: float
    oscillator_weighted_centroid_shift_ev: float | None

    def as_dict(self) -> dict[str, object]:
        return {
            "left": self.left.as_dict(),
            "right": self.right.as_dict(),
            "continuity": {
                "principal_cosines": list(
                    self.continuity.principal_cosines
                ),
                "minimum_principal_cosine": (
                    self.continuity.minimum_principal_cosine
                ),
                "chordal_distance": (
                    self.continuity.chordal_distance
                ),
            },
            "oscillator_strength_delta": (
                self.oscillator_strength_delta
            ),
            "oscillator_strength_ratio": (
                self.oscillator_strength_ratio
            ),
            "centroid_shift_ev": self.centroid_shift_ev,
            "oscillator_weighted_centroid_shift_ev": (
                self.oscillator_weighted_centroid_shift_ev
            ),
        }


def block_spectral_observables(
    result: MolecularTDAResult,
    roots: Iterable[int],
) -> BlockSpectralObservables:
    indices = _root_indices(
        roots,
        nstates=len(result.states),
    )
    selected = tuple(result.states[index] for index in indices)

    energies = tuple(state.excitation_ev for state in selected)
    strengths = tuple(
        state.oscillator_strength_oes
        for state in selected
    )
    dipole_strengths = tuple(
        state.transition_dipole_norm_au ** 2
        for state in selected
    )
    if any(
        not math.isfinite(value)
        for value in (*energies, *strengths, *dipole_strengths)
    ):
        raise BlockResponseError("block observables must be finite")
    if any(value < 0.0 for value in strengths):
        raise BlockResponseError(
            "oscillator strengths must be non-negative"
        )

    block_f = float(sum(strengths))
    all_f = float(
        sum(
            state.oscillator_strength_oes
            for state in result.states
        )
    )
    if all_f < 0.0 or not math.isfinite(all_f):
        raise BlockResponseError(
            "computed-window oscillator-strength sum is invalid"
        )

    weighted = None
    if block_f > 0.0:
        weighted = float(
            sum(
                energy * strength
                for energy, strength in zip(energies, strengths)
            )
            / block_f
        )

    return BlockSpectralObservables(
        roots=tuple(index + 1 for index in indices),
        excitation_min_ev=float(min(energies)),
        excitation_max_ev=float(max(energies)),
        excitation_centroid_ev=float(
            sum(energies) / len(energies)
        ),
        oscillator_weighted_centroid_ev=weighted,
        oscillator_strength_sum=block_f,
        dipole_strength_sum_au2=float(
            sum(dipole_strengths)
        ),
        oscillator_strength_fraction_of_computed_window=(
            block_f / all_f if all_f > 0.0 else 0.0
        ),
    )


def tracked_block_response(
    left_runtime: MolecularTDARuntime,
    right_runtime: MolecularTDARuntime,
    *,
    left_roots: Iterable[int],
    right_roots: Iterable[int],
) -> TrackedBlockResponse:
    left_root_tuple = tuple(left_roots)
    right_root_tuple = tuple(right_roots)

    continuity = runtime_block_continuity(
        left_runtime,
        right_runtime,
        left_roots=left_root_tuple,
        right_roots=right_root_tuple,
    )
    left = block_spectral_observables(
        left_runtime.result,
        left_root_tuple,
    )
    right = block_spectral_observables(
        right_runtime.result,
        right_root_tuple,
    )

    ratio = None
    if left.oscillator_strength_sum > 0.0:
        ratio = (
            right.oscillator_strength_sum
            / left.oscillator_strength_sum
        )

    weighted_shift = None
    if (
        left.oscillator_weighted_centroid_ev is not None
        and right.oscillator_weighted_centroid_ev is not None
    ):
        weighted_shift = (
            right.oscillator_weighted_centroid_ev
            - left.oscillator_weighted_centroid_ev
        )

    return TrackedBlockResponse(
        left=left,
        right=right,
        continuity=continuity,
        oscillator_strength_delta=(
            right.oscillator_strength_sum
            - left.oscillator_strength_sum
        ),
        oscillator_strength_ratio=ratio,
        centroid_shift_ev=(
            right.excitation_centroid_ev
            - left.excitation_centroid_ev
        ),
        oscillator_weighted_centroid_shift_ev=weighted_shift,
    )


__all__ = [
    "BlockResponseError",
    "BlockSpectralObservables",
    "TrackedBlockResponse",
    "block_spectral_observables",
    "tracked_block_response",
]
