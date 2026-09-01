"""Leading low-Z QED reference model for the hydrogen 2S-2P1/2 Lamb interval.

This module intentionally implements only the declared leading one-loop
self-energy terms A41/A40 plus the leading Uehling vacuum-polarization V40.
Higher-order self-energy, recoil, two-loop, finite-size and related corrections
remain separate OPEN gates.

Reference coefficients:
- ln k0(2S) = 2.811769893
- ln k0(2P) = -0.030016709
- A41(2S) = 4/3
- A40(2S) = 10/9 - (4/3) ln k0(2S)
- A40(2P1/2) = -1/6 - (4/3) ln k0(2P)
- V40(2S) = -4/15, V40(2P) = 0
"""

from __future__ import annotations

from math import log, pi

from .hydrogen import C, E_CHARGE, M_E, PLANCK, REDUCED_MASS_H
from .relativity import ALPHA

BETHE_LOG_2S = 2.811_769_893
BETHE_LOG_2P = -0.030_016_709


def _one_loop_prefactor_j(n: int) -> float:
    """(alpha/pi)(Z alpha)^4/n^3 (m_r/m_e)^3 m_e c^2 for Z=1."""
    if n < 1:
        raise ValueError("n must be positive")
    return (
        (ALPHA / pi)
        * ALPHA**4
        / n**3
        * (REDUCED_MASS_H / M_E) ** 3
        * M_E
        * C**2
    )


def leading_2s_self_energy_hz() -> float:
    pref = _one_loop_prefactor_j(2)
    large_log = log((M_E / REDUCED_MASS_H) * ALPHA**-2)
    a41 = 4.0 / 3.0
    a40 = 10.0 / 9.0 - (4.0 / 3.0) * BETHE_LOG_2S
    return pref * (a41 * large_log + a40) / PLANCK


def leading_2p1_2_self_energy_hz() -> float:
    pref = _one_loop_prefactor_j(2)
    a40 = -1.0 / 6.0 - (4.0 / 3.0) * BETHE_LOG_2P
    return pref * a40 / PLANCK


def leading_2s_vacuum_polarization_hz() -> float:
    pref = _one_loop_prefactor_j(2)
    v40 = -4.0 / 15.0
    return pref * v40 / PLANCK


def leading_lamb_2s_2p1_2_components_mhz() -> dict[str, float]:
    """Return a controlled leading approximation to the H 2S-2P1/2 interval."""
    se_2s = leading_2s_self_energy_hz() / 1e6
    se_2p = leading_2p1_2_self_energy_hz() / 1e6
    vp_2s = leading_2s_vacuum_polarization_hz() / 1e6
    interval = se_2s + vp_2s - se_2p
    return {
        "self_energy_2s_MHz": se_2s,
        "self_energy_2p1_2_MHz": se_2p,
        "vacuum_polarization_2s_MHz": vp_2s,
        "leading_interval_MHz": interval,
    }
