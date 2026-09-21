"""Basis-invariant spatial decomposition of tracked TDA transition blocks."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

import numpy as np

from .molecular_tda import MolecularTDARuntime


class SpatialTransitionError(ValueError):
    pass


def _root_indices(
    roots: Iterable[int],
    *,
    nstates: int,
) -> tuple[int, ...]:
    values = tuple(roots)
    if not values:
        raise SpatialTransitionError("root block must not be empty")
    out = []
    for value in values:
        if isinstance(value, (bool, np.bool_)):
            raise SpatialTransitionError("root labels must be integers")
        integer = int(value)
        if integer != value or integer < 1 or integer > nstates:
            raise SpatialTransitionError(
                f"root labels must lie in 1..{nstates}"
            )
        out.append(integer - 1)
    if len(set(out)) != len(out):
        raise SpatialTransitionError("root labels must be unique")
    return tuple(out)


@dataclass(frozen=True)
class BlockSpatialDecomposition:
    roots: tuple[int, ...]
    grid_level: int
    n_grid_points: int
    max_grid_dipole_delta_au: float
    max_transition_charge_residual: float
    transition_density_l2_activity_bohr_minus3: float
    dipole_envelope_au: float
    net_block_dipole_norm_au: float
    dipole_survival_ratio: float
    dipole_cancellation_fraction: float
    activity_centroid_bohr: tuple[float, float, float]
    activity_rms_radius_bohr: float

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def block_spatial_invariants(
    transition_densities_on_grid,
    weights,
    coordinates_relative_bohr,
    transition_dipoles_au,
    *,
    survival_atol: float = 1.0e-8,
) -> dict[str, object]:
    """Compute root-basis-invariant real-space block diagnostics.

    Under a unitary root mixing U, rho -> U rho and mu -> U mu. The pointwise
    sum of |rho|^2, block dipole Frobenius norm, envelope, centroid, and
    survival ratio are unchanged.
    """
    rho = np.asarray(transition_densities_on_grid, dtype=complex)
    w = np.asarray(weights, dtype=float)
    coords = np.asarray(coordinates_relative_bohr, dtype=float)
    dipoles = np.asarray(transition_dipoles_au, dtype=complex)

    if rho.ndim != 2 or rho.shape[0] < 1 or rho.shape[1] < 1:
        raise SpatialTransitionError(
            "transition_densities_on_grid must have shape (n_roots, n_points)"
        )
    n_roots, n_points = rho.shape
    if w.shape != (n_points,):
        raise SpatialTransitionError("weights shape mismatch")
    if coords.shape != (n_points, 3):
        raise SpatialTransitionError("coordinates shape mismatch")
    if dipoles.shape != (n_roots, 3):
        raise SpatialTransitionError("transition_dipoles_au shape mismatch")
    arrays = (rho.real, rho.imag, w, coords, dipoles.real, dipoles.imag)
    if any(not np.all(np.isfinite(value)) for value in arrays):
        raise SpatialTransitionError("spatial block inputs must be finite")
    if np.any(w < 0.0):
        raise SpatialTransitionError("quadrature weights must be non-negative")

    activity = np.sum(np.abs(rho) ** 2, axis=0)
    l2_activity = float(np.dot(w, activity))
    if l2_activity <= 0.0 or not math.isfinite(l2_activity):
        raise SpatialTransitionError(
            "transition-density block has no finite spatial activity"
        )

    radius = np.linalg.norm(coords, axis=1)
    envelope_density = radius * np.sqrt(activity)
    envelope = float(np.dot(w, envelope_density))
    if envelope < 0.0 or not math.isfinite(envelope):
        raise SpatialTransitionError("dipole envelope is invalid")

    net_dipole = float(
        math.sqrt(
            np.sum(np.abs(dipoles) ** 2)
        )
    )
    if envelope == 0.0:
        if net_dipole > survival_atol:
            raise SpatialTransitionError(
                "non-zero block dipole with zero spatial envelope"
            )
        survival = 0.0
    else:
        survival = net_dipole / envelope
        if survival > 1.0 + survival_atol:
            raise SpatialTransitionError(
                "dipole survival exceeds the triangle-inequality bound"
            )
        survival = min(max(survival, 0.0), 1.0)

    probability = w * activity
    normalization = float(np.sum(probability))
    centroid = np.sum(
        probability[:, None] * coords,
        axis=0,
    ) / normalization
    displacement = coords - centroid[None, :]
    rms_radius = float(
        math.sqrt(
            np.sum(
                probability
                * np.sum(displacement * displacement, axis=1)
            )
            / normalization
        )
    )

    return {
        "transition_density_l2_activity_bohr_minus3": l2_activity,
        "dipole_envelope_au": envelope,
        "net_block_dipole_norm_au": net_dipole,
        "dipole_survival_ratio": survival,
        "dipole_cancellation_fraction": 1.0 - survival,
        "activity_centroid_bohr": tuple(
            float(value) for value in centroid
        ),
        "activity_rms_radius_bohr": rms_radius,
    }


def evaluate_runtime_block_spatial(
    runtime: MolecularTDARuntime,
    roots: Iterable[int],
    *,
    grid_level: int = 3,
    dipole_grid_atol: float = 5.0e-5,
    dipole_grid_rtol: float = 5.0e-4,
) -> BlockSpatialDecomposition:
    """Evaluate a tracked TDA block on a PySCF atom-centered numerical grid."""
    if not isinstance(runtime, MolecularTDARuntime):
        raise SpatialTransitionError(
            "runtime must be MolecularTDARuntime"
        )
    if isinstance(grid_level, bool) or int(grid_level) != grid_level:
        raise SpatialTransitionError("grid_level must be an integer")
    level = int(grid_level)
    if level < 0:
        raise SpatialTransitionError("grid_level must be non-negative")
    atol = float(dipole_grid_atol)
    rtol = float(dipole_grid_rtol)
    if (
        not math.isfinite(atol)
        or atol < 0.0
        or not math.isfinite(rtol)
        or rtol < 0.0
    ):
        raise SpatialTransitionError(
            "grid dipole tolerances must be finite and non-negative"
        )

    indices = _root_indices(
        roots,
        nstates=len(runtime.transition_densities_ao),
    )

    try:
        from pyscf import dft
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "spatial transition density requires the OES q1 extra (PySCF)"
        ) from exc

    grids = dft.gen_grid.Grids(runtime.molecule)
    grids.level = level
    grids.build()
    coordinates = np.asarray(grids.coords, dtype=float)
    weights = np.asarray(grids.weights, dtype=float)
    ao = np.asarray(
        dft.numint.eval_ao(
            runtime.molecule,
            coordinates,
            deriv=0,
        ),
        dtype=complex,
    )
    if ao.shape != (
        coordinates.shape[0],
        runtime.ao_overlap.shape[0],
    ):
        raise RuntimeError("unexpected AO grid-evaluation shape")

    transition_matrices = np.asarray(
        [
            runtime.transition_densities_ao[index]
            for index in indices
        ],
        dtype=complex,
    )
    rho = np.einsum(
        "xp,npq,xq->nx",
        ao.conj(),
        transition_matrices,
        ao,
        optimize=True,
    )

    charges = np.asarray(
        runtime.molecule.atom_charges(),
        dtype=float,
    )
    atom_coordinates = np.asarray(
        runtime.molecule.atom_coords(unit="Bohr"),
        dtype=float,
    )
    charge_sum = float(np.sum(charges))
    if charge_sum <= 0.0:
        raise RuntimeError("molecular nuclear charge must be positive")
    charge_center = np.sum(
        charges[:, None] * atom_coordinates,
        axis=0,
    ) / charge_sum
    relative = coordinates - charge_center[None, :]

    reference_dipoles = np.asarray(
        [
            runtime.result.states[index].transition_dipole_au
            for index in indices
        ],
        dtype=complex,
    )
    grid_dipoles = np.einsum(
        "x,nx,xk->nk",
        weights,
        rho,
        relative,
        optimize=True,
    )
    deltas = np.max(
        np.abs(grid_dipoles - reference_dipoles),
        axis=1,
    )
    scales = np.max(
        np.abs(reference_dipoles),
        axis=1,
    )
    limits = atol + rtol * scales
    if np.any(deltas > limits):
        raise RuntimeError(
            "real-space grid transition dipole does not reproduce "
            "the backend dipole within tolerance"
        )

    charge_residuals = np.einsum(
        "x,nx->n",
        weights,
        rho,
        optimize=True,
    )

    invariants = block_spatial_invariants(
        rho,
        weights,
        relative,
        reference_dipoles,
    )
    return BlockSpatialDecomposition(
        roots=tuple(index + 1 for index in indices),
        grid_level=level,
        n_grid_points=int(coordinates.shape[0]),
        max_grid_dipole_delta_au=float(np.max(deltas)),
        max_transition_charge_residual=float(
            np.max(np.abs(charge_residuals))
        ),
        **invariants,
    )


__all__ = [
    "SpatialTransitionError",
    "BlockSpatialDecomposition",
    "block_spatial_invariants",
    "evaluate_runtime_block_spatial",
]
