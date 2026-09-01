"""Hydrogen radial geometry for Orbital Eclipse Spectroscopy.

All radial coordinates are expressed in units u=r/a_H.  The returned dipole
radial integral is therefore in units of a_H and the cancellation coherence is
dimensionless.
"""

from __future__ import annotations

from math import comb, exp, factorial, sqrt
from typing import Callable


def generalized_laguerre(k: int, alpha: int, x: float) -> float:
    """Integer-alpha generalized Laguerre polynomial L_k^alpha(x)."""
    if k < 0 or alpha < 0:
        raise ValueError("k and alpha must be non-negative")
    total = 0.0
    for m in range(k + 1):
        total += ((-1) ** m) * comb(k + alpha, k - m) * (x**m) / factorial(m)
    return total


def radial_wavefunction_dimensionless(n: int, l: int, u: float) -> float:
    """Return a_H^(3/2) R_nl(r) at u=r/a_H for hydrogen."""
    if n < 1 or l < 0 or l >= n or u < 0:
        raise ValueError("require n>=1, 0<=l<n, u>=0")
    rho = 2.0 * u / n
    prefactor = sqrt(
        (2.0 / n) ** 3
        * factorial(n - l - 1)
        / (2.0 * n * factorial(n + l))
    )
    return (
        prefactor
        * exp(-rho / 2.0)
        * (rho**l)
        * generalized_laguerre(n - l - 1, 2 * l + 1, rho)
    )


def _composite_simpson(f: Callable[[float], float], a: float, b: float, steps: int) -> float:
    if steps < 2:
        raise ValueError("steps must be >= 2")
    if steps % 2:
        steps += 1
    h = (b - a) / steps
    total = f(a) + f(b)
    for i in range(1, steps):
        total += (4.0 if i % 2 else 2.0) * f(a + i * h)
    return total * h / 3.0


def radial_transition_integrals(
    n_i: int,
    l_i: int,
    n_f: int,
    l_f: int,
    *,
    steps: int = 12_000,
) -> dict[str, float]:
    """Compute signed/absolute E1 radial overlap and cancellation coherence.

    g(u)=R_i(u) R_f(u) u^3.

    signed = integral g(u) du
    absolute = integral |g(u)| du
    coherence = |signed|/absolute

    `signed` and `absolute` are in units of a_H; `coherence` lies in [0,1]
    up to numerical roundoff.  This is the one-electron radial analogue of a
    cancellation diagnostic, not an independent replacement for the E1 matrix
    element.
    """
    for n, l in ((n_i, l_i), (n_f, l_f)):
        if n < 1 or l < 0 or l >= n:
            raise ValueError("require n>=1 and 0<=l<n")
    n_max = max(n_i, n_f)
    u_max = 10.0 * n_max * n_max + 20.0

    def integrand(u: float) -> float:
        return (
            radial_wavefunction_dimensionless(n_i, l_i, u)
            * radial_wavefunction_dimensionless(n_f, l_f, u)
            * u**3
        )

    signed = _composite_simpson(integrand, 0.0, u_max, steps)
    absolute = _composite_simpson(lambda u: abs(integrand(u)), 0.0, u_max, steps)
    coherence = abs(signed) / absolute if absolute else 0.0
    return {
        "signed_aH": signed,
        "absolute_aH": absolute,
        "coherence": coherence,
        "u_max": u_max,
        "steps": float(steps if steps % 2 == 0 else steps + 1),
    }
