"""Rotation-covariant external orbital bath selection for OES-Q1.

Determinant-by-determinant Q selection depends on the arbitrary one-particle
basis used inside a degenerate external subspace. This module instead forms
physical coupling-response states

    |chi_c> = Q H |Psi_c^P>

and compresses their spin-summed external one-body density. Under an orthogonal
rotation of the Q orbitals this density transforms covariantly, so its complete
eigenspaces define the same physical one-particle subspace. Numerically
degenerate occupation groups are always admitted together.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np

from .external_dressing import ExternalCouplingSpace
from .fermions import one_rdm
from .helium_q1 import spatial_transition_rdm


@dataclass(frozen=True)
class ExternalNaturalBath:
    q_rotation: np.ndarray
    occupations: np.ndarray
    occupation_group_sizes: Tuple[int, ...]
    occupation_group_values: Tuple[float, ...]
    selected_group_indices: Tuple[int, ...]
    selected_external_spatial_orbitals: int
    retained_normalized_occupation: float
    class_external_traces: Dict[str, float]


def _degenerate_groups(
    values: np.ndarray,
    relative_tolerance: float = 1e-6,
    absolute_tolerance: float = 1e-10,
) -> List[List[int]]:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1:
        raise ValueError("values must be one-dimensional")
    if len(values) == 0:
        return []
    groups: List[List[int]] = []
    current = [0]
    reference = float(values[0])
    for idx in range(1, len(values)):
        value = float(values[idx])
        tol = max(absolute_tolerance, relative_tolerance * max(abs(reference), abs(value)))
        if abs(value - reference) <= tol:
            current.append(idx)
        else:
            groups.append(current)
            current = [idx]
            reference = value
    groups.append(current)
    return groups


def external_coupling_response_density(
    external: ExternalCouplingSpace,
    class_states: Mapping[str, np.ndarray],
    n_active_spatial: int,
) -> tuple[np.ndarray, Dict[str, float]]:
    """Return equal-class spin-summed Q one-RDM from normalized QH|Psi> states.

    Multi-component classes (for example an exactly degenerate p manifold) are
    averaged over all orthonormal components before the class density is
    trace-normalized. Equal class weights are then used in the final density.
    """
    n_full_spin = external.n_full_spin_orbitals
    n_full_spatial = n_full_spin // 2
    if n_active_spatial < 0 or n_active_spatial >= n_full_spatial:
        raise ValueError("invalid active-spatial partition")
    names = tuple(class_states)
    if not names:
        raise ValueError("at least one state class is required")

    q_slice = slice(n_active_spatial, n_full_spatial)
    rho = np.zeros((n_full_spatial - n_active_spatial,) * 2, dtype=complex)
    traces: Dict[str, float] = {}

    for name in names:
        states = np.asarray(class_states[name], dtype=complex)
        if states.ndim == 1:
            states = states[:, None]
        if states.ndim != 2 or states.shape[0] != external.coupling_qp.shape[1]:
            raise ValueError(f"class {name} has incompatible P-space dimension")
        gram = states.conj().T @ states
        if not np.allclose(gram, np.eye(states.shape[1]), atol=1e-10):
            raise ValueError(f"class {name} states must be orthonormal")

        class_rho = np.zeros_like(rho)
        valid = 0
        for k in range(states.shape[1]):
            chi = external.coupling_qp @ states[:, k]
            norm = float(np.linalg.norm(chi))
            if norm <= 1e-14:
                continue
            chi = chi / norm
            gamma_spin = one_rdm(
                chi,
                external.external_basis,
                external.n_full_spin_orbitals,
            )
            gamma_space = spatial_transition_rdm(gamma_spin)
            class_rho += gamma_space[q_slice, q_slice]
            valid += 1
        if valid == 0:
            raise RuntimeError(f"class {name} has zero QH coupling response")
        class_rho /= valid
        class_rho = 0.5 * (class_rho + class_rho.conj().T)
        trace = float(np.trace(class_rho).real)
        if trace <= 1e-14:
            raise RuntimeError(f"class {name} has zero external one-body trace")
        traces[name] = trace
        rho += (class_rho / trace) / len(names)

    rho = 0.5 * (rho + rho.conj().T)
    if not np.allclose(rho.imag, 0.0, atol=1e-10):
        raise RuntimeError("external response density unexpectedly complex")
    rho = rho.real
    if not np.isclose(float(np.trace(rho)), 1.0, atol=1e-9):
        raise RuntimeError("normalized external response density lost unit trace")
    return rho, traces


def select_external_natural_bath(
    external: ExternalCouplingSpace,
    class_states: Mapping[str, np.ndarray],
    n_active_spatial: int,
    target_external_spatial: int,
    relative_degeneracy_tolerance: float = 1e-6,
    absolute_degeneracy_tolerance: float = 1e-10,
) -> ExternalNaturalBath:
    """Select complete natural-orbital occupation groups up to a target size."""
    if target_external_spatial < 1:
        raise ValueError("target_external_spatial must be positive")
    rho, traces = external_coupling_response_density(
        external,
        class_states,
        n_active_spatial=n_active_spatial,
    )
    occ, vec = np.linalg.eigh(rho)
    order = np.argsort(occ)[::-1]
    occ = np.asarray(occ[order], dtype=float)
    vec = np.asarray(vec[:, order], dtype=float)
    groups = _degenerate_groups(
        occ,
        relative_tolerance=relative_degeneracy_tolerance,
        absolute_tolerance=absolute_degeneracy_tolerance,
    )

    selected_groups: List[int] = []
    selected_indices: List[int] = []
    for gidx, group in enumerate(groups):
        selected_groups.append(gidx)
        selected_indices.extend(group)
        if len(selected_indices) >= target_external_spatial:
            break
    if len(selected_indices) < target_external_spatial:
        raise RuntimeError("target exceeds available external natural-orbital groups")

    retained = float(np.sum(occ[selected_indices]))
    group_values = tuple(float(np.mean(occ[group])) for group in groups)
    return ExternalNaturalBath(
        q_rotation=vec[:, selected_indices],
        occupations=occ,
        occupation_group_sizes=tuple(len(group) for group in groups),
        occupation_group_values=group_values,
        selected_group_indices=tuple(selected_groups),
        selected_external_spatial_orbitals=len(selected_indices),
        retained_normalized_occupation=retained,
        class_external_traces=traces,
    )


def bath_projector_overlap(
    C_a: np.ndarray,
    C_b: np.ndarray,
    overlap_ao: np.ndarray,
) -> Dict[str, float]:
    """Compare two S-orthonormal orbital subspaces by their principal angles."""
    C_a = np.asarray(C_a, dtype=float)
    C_b = np.asarray(C_b, dtype=float)
    S = np.asarray(overlap_ao, dtype=float)
    if C_a.shape[1] != C_b.shape[1]:
        raise ValueError("bath subspaces have different dimensions")
    singular = np.linalg.svd(C_a.T @ S @ C_b, compute_uv=False)
    return {
        "min_principal_cosine": float(np.min(singular)),
        "max_principal_cosine_error": float(np.max(np.abs(1.0 - singular))),
        "projector_frobenius_distance": float(
            np.sqrt(max(0.0, 2.0 * C_a.shape[1] - 2.0 * np.sum(singular**2)))
        ),
    }
