"""Magnetic-field flavor evolution for hydrogen p states.

Hamiltonian (frequency units):
    H/h = A (L.S / hbar^2) + (mu_B/h) B (m_l + 2 m_s)

For l=1, s=1/2 the spin-orbit constant A is chosen from the OES Dirac
fine-structure splitting: Delta_nu_FS = 3A/2.  No field-dependent spectral
line data are used.
"""

from __future__ import annotations

from math import sqrt

from .hydrogen import MU_B_OVER_H
from .relativity import fine_structure_split_hz


def _eig2(a: float, b: float, d: float) -> tuple[float, float]:
    """Eigenvalues of [[a,b],[b,d]], ascending."""
    center = 0.5 * (a + d)
    radius = sqrt((0.5 * (a - d)) ** 2 + b * b)
    return center - radius, center + radius


def p_spin_orbit_constant_hz(n: int) -> float:
    """A/h in Hz for H_SO/h = A_hz (L.S/hbar^2), l=1."""
    if n < 2:
        raise ValueError("p states require n >= 2")
    return fine_structure_split_hz(n, 0.5, 1.5) / 1.5


def p_paschen_back_levels_hz(n: int, b_tesla: float) -> list[dict[str, float | str]]:
    """Six l=1,s=1/2 eigenlevels across the Zeeman/Paschen-Back crossover.

    Energies are offsets in Hz relative to an arbitrary common n,p gross level.
    The labels identify conserved m_j and an adiabatic branch within mixed
    m_j=+/-1/2 blocks.
    """
    if b_tesla < 0:
        raise ValueError("B magnitude must be non-negative")
    a_so = p_spin_orbit_constant_hz(n)
    z = MU_B_OVER_H * b_tesla

    levels: list[dict[str, float | str]] = []

    # Pure stretched states |m_l=+/-1, m_s=+/-1/2>; L.S = +1/2.
    levels.append({"m_j": 1.5, "branch": "stretched+", "energy_hz": 0.5 * a_so + 2.0 * z})
    levels.append({"m_j": -1.5, "branch": "stretched-", "energy_hz": 0.5 * a_so - 2.0 * z})

    # m_j=+1/2 block in basis |0,+1/2>, |+1,-1/2>.
    # L.S/hbar^2 = [[0, sqrt(2)/2], [sqrt(2)/2, -1/2]].
    e_low, e_high = _eig2(
        z,
        a_so / sqrt(2.0),
        -0.5 * a_so,
    )
    levels.append({"m_j": 0.5, "branch": "lower", "energy_hz": e_low})
    levels.append({"m_j": 0.5, "branch": "upper", "energy_hz": e_high})

    # m_j=-1/2 block in basis |-1,+1/2>, |0,-1/2>.
    e_low, e_high = _eig2(
        -0.5 * a_so,
        a_so / sqrt(2.0),
        -z,
    )
    levels.append({"m_j": -0.5, "branch": "lower", "energy_hz": e_low})
    levels.append({"m_j": -0.5, "branch": "upper", "energy_hz": e_high})

    return sorted(levels, key=lambda row: float(row["energy_hz"]))


def p_mj_half_branches_hz(n: int, b_tesla: float, sign: int = 1) -> tuple[float, float]:
    """Return lower/upper mixed branches for m_j=sign*1/2."""
    if sign not in (-1, 1):
        raise ValueError("sign must be +/-1")
    rows = [
        row
        for row in p_paschen_back_levels_hz(n, b_tesla)
        if row["m_j"] == 0.5 * sign and row["branch"] in ("lower", "upper")
    ]
    values = sorted(float(row["energy_hz"]) for row in rows)
    return values[0], values[1]
