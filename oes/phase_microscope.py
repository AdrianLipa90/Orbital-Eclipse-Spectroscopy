"""Phase-microscope fields reconstructed from an OES transition 1-RDM.\n\nMatched orbital rephasing leaves the transition density unchanged. A separate\nglobal phase of the initial/final state rays multiplies the whole transition\ndensity by one U(1) phase, so absolute raw phase is not an observable without\na reference; amplitude, nodes and relative phase are invariant.\n"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .invariant_contract import relative_transition_phase
from .transition_contract import OESTransitionState


@dataclass(frozen=True)
class PhaseMicroscopeField:
    transition_density: np.ndarray
    amplitude: np.ndarray
    phase_rad: np.ndarray
    phase_defined: np.ndarray

    def reconstruct(self) -> np.ndarray:
        return self.amplitude * np.exp(1j * self.phase_rad)

    def relative_phase(
        self,
        *,
        reference_index: int | None = None,
        amplitude_floor: float = 1.0e-14,
    ) -> tuple[np.ndarray, np.ndarray, int | None]:
        """Remove the arbitrary global initial/final state-ray phase.

        Raw ``phase_rad`` is a representation-level diagnostic. Physical phase
        claims require a declared external phase reference or a relative phase
        returned by this method.
        """
        return relative_transition_phase(
            self.transition_density,
            amplitude_floor=amplitude_floor,
            reference_index=reference_index,
        )


def transition_density_on_points(
    transition: OESTransitionState,
    orbital_values: np.ndarray,
) -> np.ndarray:
    """Evaluate rho_FI(x)=sum_pq T_pq phi_p*(x) phi_q(x).

    orbital_values[x,p] must use the same spatial-orbital gauge and ordering
    named by transition.provenance.orbital_basis_id.
    """
    phi = np.asarray(orbital_values, dtype=complex)
    if phi.ndim != 2 or phi.shape[1] != transition.n_spatial_orbitals:
        raise ValueError("orbital_values must have shape (n_points, n_spatial_orbitals)")
    if not np.all(np.isfinite(phi.real)) or not np.all(np.isfinite(phi.imag)):
        raise ValueError("orbital_values must be finite")
    return np.einsum("xp,pq,xq->x", np.conjugate(phi), transition.transition_rdm, phi)


def phase_microscope_field(
    transition: OESTransitionState,
    orbital_values: np.ndarray,
    *,
    amplitude_floor: float = 1.0e-14,
) -> PhaseMicroscopeField:
    if not np.isfinite(amplitude_floor) or amplitude_floor < 0.0:
        raise ValueError("amplitude_floor must be finite and non-negative")
    rho = transition_density_on_points(transition, orbital_values)
    amplitude = np.abs(rho)
    defined = amplitude > float(amplitude_floor)
    phase = np.zeros_like(amplitude)
    phase[defined] = np.angle(rho[defined])
    return PhaseMicroscopeField(
        transition_density=rho,
        amplitude=amplitude,
        phase_rad=phase,
        phase_defined=defined,
    )


__all__ = [
    "PhaseMicroscopeField",
    "transition_density_on_points",
    "phase_microscope_field",
]
