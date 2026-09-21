"""Orientation-dependent transition observables for OES.

The geometry/state generator is upstream. OES owns the transition projection:
line position, dipole strength, and response of those observables to a declared
orientation coordinate.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

import numpy as np

from .line_observables import electric_dipole_absorption_line
from .transition_contract import OESTransitionState


class OrientationResponseError(ValueError):
    pass


def _finite(value: float, name: str) -> float:
    x = float(value)
    if not math.isfinite(x):
        raise OrientationResponseError(f"{name} must be finite")
    return x


@dataclass(frozen=True)
class OrientationLinePoint:
    angle_rad: float
    energy_gap_hartree: float
    frequency_hz: float
    wavelength_nm: float
    transition_dipole_au: tuple[complex, complex, complex]
    transition_dipole_norm_au: float
    oscillator_strength_length_gauge: float

    def as_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["transition_dipole_au"] = [
            {"real": float(value.real), "imag": float(value.imag)}
            for value in self.transition_dipole_au
        ]
        data["schema"] = "OES_ORIENTATION_LINE_POINT_V0_1"
        return data


@dataclass(frozen=True)
class OscillatorStrengthSlope:
    energy_gap_term_per_rad: float
    transition_moment_term_per_rad: float
    total_per_rad: float


def orientation_line_point(
    angle_rad: float,
    transition: OESTransitionState,
    dipole_operators_au: np.ndarray,
) -> OrientationLinePoint:
    theta = _finite(angle_rad, "angle_rad")
    line = electric_dipole_absorption_line(
        transition,
        dipole_operators_au,
    )
    return OrientationLinePoint(
        angle_rad=theta,
        energy_gap_hartree=transition.signed_energy_gap_hartree,
        frequency_hz=line.frequency_hz,
        wavelength_nm=line.wavelength_nm,
        transition_dipole_au=line.transition_dipole_au,
        transition_dipole_norm_au=line.transition_dipole_norm_au,
        oscillator_strength_length_gauge=(
            line.oscillator_strength_length_gauge
        ),
    )


def orientation_line_scan(
    samples: Iterable[
        tuple[float, OESTransitionState, np.ndarray]
    ],
) -> tuple[OrientationLinePoint, ...]:
    rows = tuple(samples)
    if not rows:
        raise OrientationResponseError("orientation scan must not be empty")
    return tuple(
        orientation_line_point(theta, transition, dipole)
        for theta, transition, dipole in rows
    )


def oscillator_strength_slope(
    *,
    energy_gap_hartree: float,
    energy_gap_slope_hartree_per_rad: float,
    transition_dipole_au,
    transition_dipole_slope_au_per_rad,
) -> OscillatorStrengthSlope:
    """Decompose df/dtheta for f=(2/3) DeltaE |mu|^2 in atomic units."""
    gap = _finite(energy_gap_hartree, "energy_gap_hartree")
    gap_slope = _finite(
        energy_gap_slope_hartree_per_rad,
        "energy_gap_slope_hartree_per_rad",
    )
    if gap <= 0.0:
        raise OrientationResponseError(
            "energy_gap_hartree must be positive for absorption"
        )

    mu = np.asarray(transition_dipole_au, dtype=complex)
    dmu = np.asarray(
        transition_dipole_slope_au_per_rad,
        dtype=complex,
    )
    if mu.shape != (3,) or dmu.shape != (3,):
        raise OrientationResponseError(
            "transition dipole and slope must be complex 3-vectors"
        )
    values = (mu.real, mu.imag, dmu.real, dmu.imag)
    if any(not np.all(np.isfinite(value)) for value in values):
        raise OrientationResponseError(
            "transition dipole and slope must be finite"
        )

    mu2 = float(np.vdot(mu, mu).real)
    cross = float(np.vdot(mu, dmu).real)
    energy_term = (2.0 / 3.0) * gap_slope * mu2
    transition_moment_term = (4.0 / 3.0) * gap * cross
    return OscillatorStrengthSlope(
        energy_gap_term_per_rad=energy_term,
        transition_moment_term_per_rad=transition_moment_term,
        total_per_rad=energy_term + transition_moment_term,
    )


__all__ = [
    "OrientationResponseError",
    "OrientationLinePoint",
    "OscillatorStrengthSlope",
    "orientation_line_point",
    "orientation_line_scan",
    "oscillator_strength_slope",
]
