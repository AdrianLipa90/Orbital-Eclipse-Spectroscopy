"""Executable profile for RELATIONAL_INVARIANT_CONTRACT_V0_1.

The primary OES gauge is a unitary change of spatial-orbital basis inside the
same orthonormal one-particle space. Diagonal U(1)^n rephasing is a subgroup.
A separate global U(1) freedom belongs to the initial/final many-electron
state rays.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

INVARIANT_CONTRACT_ID = "RELATIONAL_INVARIANT_CONTRACT_V0_1"
OES_ORBITAL_UNITARY_BASIS_PROFILE_ID = "OES_ORBITAL_UNITARY_BASIS_V0_1"


class InvariantContractError(ValueError):
    pass


def representation_contract_block() -> dict[str, object]:
    return {
        "schema": INVARIANT_CONTRACT_ID,
        "profile": OES_ORBITAL_UNITARY_BASIS_PROFILE_ID,
        "orbital_gauge_group": "U(n)",
        "transition_rdm_covariance": "T'=U^T T U^*",
        "transition_density_status": (
            "INVARIANT_UNDER_MATCHED_ORBITAL_UNITARY_BASIS_CHANGE"
        ),
        "state_ray_phase_status": (
            "GLOBAL_U1_COVARIANT__AMPLITUDE_NODES_RELATIVE_PHASE_INTENSITY_INVARIANT"
        ),
    }


def validate_representation_contract_block(block: Mapping[str, object]) -> None:
    if block.get("schema") != INVARIANT_CONTRACT_ID:
        raise InvariantContractError("unsupported relational invariant contract")
    if block.get("profile") != OES_ORBITAL_UNITARY_BASIS_PROFILE_ID:
        raise InvariantContractError("unsupported orbital-gauge profile")
    expected = representation_contract_block()
    for key in (
        "orbital_gauge_group",
        "transition_rdm_covariance",
        "transition_density_status",
        "state_ray_phase_status",
    ):
        if block.get(key) != expected[key]:
            raise InvariantContractError(f"representation contract mismatch: {key}")


def _finite_matrix(name: str, value) -> np.ndarray:
    out = np.asarray(value, dtype=complex)
    if out.ndim != 2 or out.shape[0] == 0 or out.shape[1] == 0:
        raise InvariantContractError(f"{name} must be a non-empty 2D array")
    if not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise InvariantContractError(f"{name} must be finite")
    return out


def _unitary_matrix(unitary, n: int) -> np.ndarray:
    u = _finite_matrix("unitary", unitary)
    if u.shape != (n, n):
        raise InvariantContractError("unitary must have shape (n_orbitals, n_orbitals)")
    identity = np.eye(n, dtype=complex)
    if not np.allclose(u.conjugate().T @ u, identity, rtol=0.0, atol=1.0e-12):
        raise InvariantContractError("orbital basis transform must be unitary")
    return u


def transform_orbital_values(orbital_values, unitary) -> np.ndarray:
    """Apply phi -> phi U to orbital values stored as (n_points, n_orbitals)."""

    phi = _finite_matrix("orbital_values", orbital_values)
    u = _unitary_matrix(unitary, phi.shape[1])
    return phi @ u


def transform_transition_rdm(transition_rdm, unitary) -> np.ndarray:
    """Apply T -> U^T T U^* for the declared OES transition-RDM convention."""

    t = _finite_matrix("transition_rdm", transition_rdm)
    if t.shape[0] != t.shape[1]:
        raise InvariantContractError("transition_rdm must be square")
    u = _unitary_matrix(unitary, t.shape[0])
    return u.T @ t @ np.conjugate(u)


def transform_one_body_operator(one_body_operator, unitary) -> np.ndarray:
    """Apply d -> U^dagger d U in the same orbital basis change."""

    d = _finite_matrix("one_body_operator", one_body_operator)
    if d.shape[0] != d.shape[1]:
        raise InvariantContractError("one_body_operator must be square")
    u = _unitary_matrix(unitary, d.shape[0])
    return u.conjugate().T @ d @ u


def _phase_unitary(chi, n: int) -> np.ndarray:
    c = np.asarray(chi, dtype=float)
    if c.shape != (n,) or not np.all(np.isfinite(c)):
        raise InvariantContractError("chi must be one finite phase per orbital")
    return np.diag(np.exp(1j * c))


def rephase_orbital_values(orbital_values, chi) -> np.ndarray:
    phi = _finite_matrix("orbital_values", orbital_values)
    return transform_orbital_values(phi, _phase_unitary(chi, phi.shape[1]))


def rephase_transition_rdm(transition_rdm, chi) -> np.ndarray:
    t = _finite_matrix("transition_rdm", transition_rdm)
    if t.shape[0] != t.shape[1]:
        raise InvariantContractError("transition_rdm must be square")
    return transform_transition_rdm(t, _phase_unitary(chi, t.shape[0]))


def relative_transition_phase(
    transition_density,
    *,
    amplitude_floor: float = 1.0e-14,
    reference_index: int | None = None,
) -> tuple[np.ndarray, np.ndarray, int | None]:
    """Return phase differences with the arbitrary global transition phase removed.

    If reference_index is omitted, the maximum-amplitude defined point is used
    deterministically. Undefined points receive zero in the returned phase
    array and are marked false in the mask.
    """

    rho = np.asarray(transition_density, dtype=complex)
    if rho.ndim != 1 or rho.size == 0:
        raise InvariantContractError("transition_density must be a non-empty 1D array")
    if not np.all(np.isfinite(rho.real)) or not np.all(np.isfinite(rho.imag)):
        raise InvariantContractError("transition_density must be finite")
    floor = float(amplitude_floor)
    if not np.isfinite(floor) or floor < 0.0:
        raise InvariantContractError("amplitude_floor must be finite and non-negative")

    amplitude = np.abs(rho)
    defined = amplitude > floor
    if not np.any(defined):
        return np.zeros(rho.shape, dtype=float), defined, None

    if reference_index is None:
        candidates = np.flatnonzero(defined)
        ref = int(candidates[np.argmax(amplitude[candidates])])
    else:
        ref = int(reference_index)
        if ref < 0 or ref >= rho.size:
            raise InvariantContractError("reference_index out of range")
        if not defined[ref]:
            raise InvariantContractError("reference_index must have defined phase")

    reference = rho[ref]
    relative = np.zeros(rho.shape, dtype=float)
    relative[defined] = np.angle(rho[defined] * np.conjugate(reference))
    return relative, defined, ref


__all__ = [
    "INVARIANT_CONTRACT_ID",
    "OES_ORBITAL_UNITARY_BASIS_PROFILE_ID",
    "InvariantContractError",
    "representation_contract_block",
    "validate_representation_contract_block",
    "transform_orbital_values",
    "transform_transition_rdm",
    "transform_one_body_operator",
    "rephase_orbital_values",
    "rephase_transition_rdm",
    "relative_transition_phase",
]
