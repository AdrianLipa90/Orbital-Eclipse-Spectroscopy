"""Generic electric-dipole line observables derived from an OES transition state."""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.constants import c

from .transition_contract import OESTransitionState


@dataclass(frozen=True)
class ElectricDipoleLine:
    frequency_hz: float
    wavelength_nm: float
    transition_dipole_au: tuple[complex, complex, complex]
    transition_dipole_norm_au: float
    oscillator_strength_length_gauge: float

    def as_dict(self) -> dict[str, object]:
        return {
            "frequency_hz": self.frequency_hz,
            "wavelength_nm": self.wavelength_nm,
            "transition_dipole_au": [
                {"real": float(value.real), "imag": float(value.imag)}
                for value in self.transition_dipole_au
            ],
            "transition_dipole_norm_au": self.transition_dipole_norm_au,
            "oscillator_strength_length_gauge": self.oscillator_strength_length_gauge,
            "line_model": "ELECTRIC_DIPOLE_LENGTH_GAUGE",
        }


def electric_dipole_absorption_line(
    transition: OESTransitionState,
    dipole_operators_au: np.ndarray,
) -> ElectricDipoleLine:
    """Derive an upward electric-dipole absorption line in atomic units.

    dipole_operators_au[k,p,q] are the x/y/z one-body dipole operator
    coefficients in the same spatial-orbital basis as transition.transition_rdm.
    """
    gap = transition.signed_energy_gap_hartree
    if gap <= 0.0:
        raise ValueError("electric-dipole absorption requires final energy above initial energy")
    operators = np.asarray(dipole_operators_au, dtype=complex)
    expected = (3, transition.n_spatial_orbitals, transition.n_spatial_orbitals)
    if operators.shape != expected:
        raise ValueError(f"dipole_operators_au must have shape {expected}")
    if not np.all(np.isfinite(operators.real)) or not np.all(np.isfinite(operators.imag)):
        raise ValueError("dipole_operators_au must be finite")

    components = tuple(
        transition.transition_amplitude(operators[k])
        for k in range(3)
    )
    mu2 = float(sum(abs(value) ** 2 for value in components))
    mu_norm = math.sqrt(mu2)
    oscillator_strength = (2.0 / 3.0) * gap * mu2
    frequency = transition.frequency_hz
    wavelength_nm = c / frequency * 1.0e9
    return ElectricDipoleLine(
        frequency_hz=frequency,
        wavelength_nm=wavelength_nm,
        transition_dipole_au=components,
        transition_dipole_norm_au=mu_norm,
        oscillator_strength_length_gauge=oscillator_strength,
    )


__all__ = ["ElectricDipoleLine", "electric_dipole_absorption_line"]
