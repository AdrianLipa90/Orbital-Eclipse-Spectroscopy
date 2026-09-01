"""Relativistic hydrogen reference primitives for OES-H0.

The Dirac expression below uses the ordinary-hydrogen reduced mass as a compact
reference approximation.  It is not a full two-body recoil/QED treatment.
"""

from __future__ import annotations

from math import sqrt

from .hydrogen import C, E_CHARGE, PLANCK, REDUCED_MASS_H, COULOMB_G, HBAR

ALPHA = COULOMB_G / (HBAR * C)


def dirac_bound_energy_ev(n: int, j: float) -> float:
    """Hydrogenic Dirac binding energy in eV using reduced-mass substitution."""
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer")
    if j < 0.5 or j > n - 0.5:
        raise ValueError("j outside hydrogenic range")
    k = j + 0.5
    gamma = sqrt(k * k - ALPHA * ALPHA)
    denominator = n - j - 0.5 + gamma
    total_minus_rest_j = REDUCED_MASS_H * C**2 * (
        (1.0 + (ALPHA / denominator) ** 2) ** -0.5 - 1.0
    )
    return total_minus_rest_j / E_CHARGE


def fine_structure_split_hz(n: int, j_a: float, j_b: float) -> float:
    """Absolute Dirac fine-structure separation between two j flavors."""
    delta_ev = abs(dirac_bound_energy_ev(n, j_a) - dirac_bound_energy_ev(n, j_b))
    return delta_ev * E_CHARGE / PLANCK


def dirac_degeneracy_signature(n: int, j: float) -> tuple[int, float]:
    """Energy label showing the pure Coulomb-Dirac n,j degeneracy structure."""
    _ = dirac_bound_energy_ev(n, j)
    return (n, j)
