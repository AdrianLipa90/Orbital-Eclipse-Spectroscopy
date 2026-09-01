"""Hydrogen reference calculations for OES-H0.

The module intentionally starts from analytic hydrogenic structure rather than
measured spectral wavelengths.  Observed values belong in benchmark fixtures,
not in the solver.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import pi

# SI constants. Exact where SI defines them; mass/epsilon values are CODATA-scale
# numerical inputs and are kept explicit so provenance can later be locked.
PLANCK = 6.626_070_15e-34
HBAR = PLANCK / (2.0 * pi)
C = 299_792_458.0
E_CHARGE = 1.602_176_634e-19
EPSILON_0 = 8.854_187_8128e-12
M_E = 9.109_383_7139e-31
M_P = 1.672_621_925_95e-27
MU_B_OVER_H = 13.996_245_55e9  # Hz/T, conventional Bohr-magneton scale

REDUCED_MASS_H = M_E * M_P / (M_E + M_P)
COULOMB_G = E_CHARGE**2 / (4.0 * pi * EPSILON_0)


@dataclass(frozen=True)
class HydrogenTransition:
    n_i: int
    n_f: int
    delta_e_ev: float
    frequency_hz: float
    wavelength_nm: float


def _validate_n(n: int) -> None:
    if not isinstance(n, int) or n < 1:
        raise ValueError("principal quantum number n must be a positive integer")


def bohr_radius_hydrogen_m() -> float:
    """Reduced-mass Bohr radius for ordinary hydrogen."""
    return HBAR**2 / (REDUCED_MASS_H * COULOMB_G)


def gross_energy_ev(n: int) -> float:
    """Nonrelativistic Coulomb energy E_n in eV, with proton recoil via mu."""
    _validate_n(n)
    energy_j = -(REDUCED_MASS_H * COULOMB_G**2) / (2.0 * HBAR**2 * n**2)
    return energy_j / E_CHARGE


def gross_transition(n_i: int, n_f: int) -> HydrogenTransition:
    """Gross hydrogen emission color from n_i > n_f; no spectral table input."""
    _validate_n(n_i)
    _validate_n(n_f)
    if n_i <= n_f:
        raise ValueError("emission transition requires n_i > n_f")
    delta_e_ev = gross_energy_ev(n_i) - gross_energy_ev(n_f)
    delta_e_j = delta_e_ev * E_CHARGE
    frequency_hz = delta_e_j / PLANCK
    wavelength_nm = (C / frequency_hz) * 1e9
    return HydrogenTransition(n_i, n_f, delta_e_ev, frequency_hz, wavelength_nm)


def exposure_inv_r(n: int) -> float:
    """Hydrogenic central exposure <1/r> in m^-1."""
    _validate_n(n)
    a = bohr_radius_hydrogen_m()
    return 1.0 / (n**2 * a)


def exposure_inv_r3(n: int, l: int) -> float:
    """Hydrogenic <1/r^3> for l>0, in m^-3.

    The l=0 expectation diverges in the point-Coulomb nonrelativistic model and
    is represented by the separate contact exposure instead.
    """
    _validate_n(n)
    if l < 1 or l >= n:
        raise ValueError("exposure_inv_r3 requires 1 <= l < n")
    a = bohr_radius_hydrogen_m()
    return 1.0 / (a**3 * n**3 * l * (l + 0.5) * (l + 1.0))


def contact_density(n: int, l: int) -> float:
    """|psi_nl(0)|^2 in m^-3 for hydrogenic states."""
    _validate_n(n)
    if l < 0 or l >= n:
        raise ValueError("orbital quantum number requires 0 <= l < n")
    if l != 0:
        return 0.0
    a = bohr_radius_hydrogen_m()
    return 1.0 / (pi * a**3 * n**3)


def contact_exposure_dimensionless(n: int, l: int) -> float:
    """pi*a_H^3*|psi(0)|^2 = delta_l0/n^3."""
    return (1.0 / n**3) if l == 0 else 0.0


def radial_node_count(n: int, l: int) -> int:
    _validate_n(n)
    if l < 0 or l >= n:
        raise ValueError("orbital quantum number requires 0 <= l < n")
    return n - l - 1


def dipole_allowed(l_i: int, m_i: int, l_f: int, m_f: int) -> bool:
    """Orbital E1 selection gate, excluding spin/hyperfine structure."""
    if abs(m_i) > l_i or abs(m_f) > l_f:
        raise ValueError("|m| must not exceed l")
    return abs(l_i - l_f) == 1 and abs(m_i - m_f) <= 1


def orbital_flavor_count(n_f: int) -> int:
    """Distinct l_i<->l_f E1 routes for any n_i > n_f."""
    _validate_n(n_f)
    return 2 * n_f - 1


def m_resolved_flavor_count(n_f: int) -> int:
    """m-resolved E1 routes before spin for a gross n_i -> n_f color."""
    _validate_n(n_f)
    return 3 * (n_f**2 + (n_f - 1) ** 2)


def lande_g(l: int, j: float, s: float = 0.5) -> float:
    """LS-coupling Lande factor with g_L=1 and g_S=2."""
    if l < 0 or j <= 0 or s < 0:
        raise ValueError("invalid angular momentum quantum numbers")
    jj = j * (j + 1.0)
    return 1.0 + (jj + s * (s + 1.0) - l * (l + 1.0)) / (2.0 * jj)


def zeeman_transition_shift_hz(
    *,
    b_tesla: float,
    l_i: int,
    j_i: float,
    m_j_i: float,
    l_f: int,
    j_f: float,
    m_j_f: float,
) -> float:
    """Linear-Zeeman frequency displacement of a transition."""
    if b_tesla < 0:
        raise ValueError("B magnitude must be non-negative")
    g_i = lande_g(l_i, j_i)
    g_f = lande_g(l_f, j_f)
    return MU_B_OVER_H * b_tesla * (g_i * m_j_i - g_f * m_j_f)


def hydrogen_fingerprint(
    *, n_i: int, l_i: int, m_i: int, n_f: int, l_f: int, m_f: int
) -> dict[str, float | int | bool | dict[str, float | int]]:
    """Minimal OES color+flavor fingerprint for an orbital E1 candidate."""
    transition = gross_transition(n_i, n_f)
    return {
        "color": asdict(transition),
        "flavor": {
            "l_i": l_i,
            "l_f": l_f,
            "m_i": m_i,
            "m_f": m_f,
            "radial_nodes_i": radial_node_count(n_i, l_i),
            "radial_nodes_f": radial_node_count(n_f, l_f),
            "delta_m": m_i - m_f,
            "e1_allowed": dipole_allowed(l_i, m_i, l_f, m_f),
            "contact_i": contact_exposure_dimensionless(n_i, l_i),
            "contact_f": contact_exposure_dimensionless(n_f, l_f),
        },
    }
