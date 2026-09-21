"""Executable profile for RELATIONAL_INVARIANT_CONTRACT_V0_1.

The helpers distinguish orbital-basis rephasing invariance from the separate
global U(1) phase freedom of the initial/final many-electron state rays.
"""
from __future__ import annotations

from typing import Mapping

import numpy as np

INVARIANT_CONTRACT_ID = "RELATIONAL_INVARIANT_CONTRACT_V0_1"
OES_ORBITAL_REPHASING_PROFILE_ID = "OES_ORBITAL_U1_REPHASING_V0_1"


class InvariantContractError(ValueError):
    pass


def representation_contract_block() -> dict[str, object]:
    return {
        "schema": INVARIANT_CONTRACT_ID,
        "profile": OES_ORBITAL_REPHASING_PROFILE_ID,
        "orbital_gauge_group": "U(1)^n",
        "transition_rdm_covariance": "T'=D_chi T D_chi^dagger",
        "transition_density_status": "INVARIANT_UNDER_MATCHED_ORBITAL_REPHASING",
        "state_ray_phase_status": (
            "GLOBAL_U1_COVARIANT__AMPLITUDE_NODES_RELATIVE_PHASE_INTENSITY_INVARIANT"
        ),
    }


def validate_representation_contract_block(block: Mapping[str, object]) -> None:
    if block.get("schema") != INVARIANT_CONTRACT_ID:
        raise InvariantContractError("unsupported relational invariant contract")
    if block.get("profile") != OES_ORBITAL_REPHASING_PROFILE_ID:
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


def _phase_vector(chi, n: int) -> np.ndarray:
    c = np.asarray(chi, dtype=float)
    if c.shape != (n,) or not np.all(np.isfinite(c)):
        raise InvariantContractError("chi must be one finite phase per orbital")
    return np.exp(1j * c)


def rephase_orbital_values(orbital_values, chi) -> np.ndarray:
    phi = np.asarray(orbital_values, dtype=complex)
    if phi.ndim != 2 or phi.shape[1] == 0:
        raise InvariantContractError("orbital_values must be a non-empty 2D array")
    if not np.all(np.isfinite(phi.real)) or not np.all(np.isfinite(phi.imag)):
        raise InvariantContractError("orbital_values must be finite")
    phases = _phase_vector(chi, phi.shape[1])
    return phi * phases[None, :]


def rephase_transition_rdm(transition_rdm, chi) -> np.ndarray:
    t = np.asarray(transition_rdm, dtype=complex)
    if t.ndim != 2 or t.shape[0] != t.shape[1] or t.shape[0] == 0:
        raise InvariantContractError("transition_rdm must be a non-empty square matrix")
    if not np.all(np.isfinite(t.real)) or not np.all(np.isfinite(t.imag)):
        raise InvariantContractError("transition_rdm must be finite")
    phases = _phase_vector(chi, t.shape[0])
    return phases[:, None] * np.conjugate(phases[None, :]) * t


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
    "OES_ORBITAL_REPHASING_PROFILE_ID",
    "InvariantContractError",
    "representation_contract_block",
    "validate_representation_contract_block",
    "rephase_orbital_values",
    "rephase_transition_rdm",
    "relative_transition_phase",
]
