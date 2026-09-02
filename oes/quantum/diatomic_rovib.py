"""One-dimensional diatomic rovibrational solver on a sampled Born-Oppenheimer curve.

The electronic structure layer supplies only internuclear coordinates and total
Born-Oppenheimer energies.  This module interpolates strictly inside that sampled
interval, applies the nuclear kinetic operator in atomic units, and solves the
radial J-resolved eigenproblem with Dirichlet boundaries.  Experimental
spectroscopic constants are not inputs to this solver.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Mapping, Sequence, Tuple

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh_tridiagonal

from .h2_m1 import HARTREE_TO_WAVENUMBER_CM


@dataclass(frozen=True)
class DiatomicRovibrationalResult:
    sample_r_bohr: Tuple[float, ...]
    sample_energy_hartree: Tuple[float, ...]
    reduced_mass_me: float
    radial_grid_points: int
    levels_hartree_by_j: Dict[int, Tuple[float, ...]]
    term_values_cm_by_j: Dict[int, Tuple[float, ...]]
    boundary_margin_cm_by_j: Dict[int, float]

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def solve_diatomic_rovibrational_levels(
    sample_r_bohr: Sequence[float],
    sample_energy_hartree: Sequence[float],
    reduced_mass_me: float,
    *,
    n_vibrational: int = 4,
    j_values: Sequence[int] = (0, 1, 2),
    radial_grid_points: int = 1800,
) -> DiatomicRovibrationalResult:
    r = np.asarray(sample_r_bohr, dtype=float)
    e = np.asarray(sample_energy_hartree, dtype=float)
    mu = float(reduced_mass_me)
    n_v = int(n_vibrational)
    n_grid = int(radial_grid_points)

    if r.ndim != 1 or e.shape != r.shape or r.size < 5:
        raise ValueError("sampled BO curve requires matching one-dimensional arrays with at least five points")
    if not np.all(np.isfinite(r)) or not np.all(np.isfinite(e)):
        raise ValueError("sampled BO curve must be finite")
    if not np.all(np.diff(r) > 0.0):
        raise ValueError("sample_r_bohr must be strictly increasing")
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("reduced_mass_me must be finite and positive")
    if n_v < 1:
        raise ValueError("n_vibrational must be positive")
    if n_grid < 200:
        raise ValueError("radial_grid_points must be at least 200")

    j_tuple = tuple(int(j) for j in j_values)
    if not j_tuple or any(j < 0 for j in j_tuple) or len(set(j_tuple)) != len(j_tuple):
        raise ValueError("j_values must contain unique non-negative integers")
    if 0 not in j_tuple:
        raise ValueError("j_values must include J=0")

    minimum_index = int(np.argmin(e))
    if minimum_index in (0, len(e) - 1):
        raise RuntimeError("sampled BO minimum lies on an interval boundary")

    potential = e - float(np.min(e))
    interpolator = PchipInterpolator(r, potential, extrapolate=False)

    full_grid = np.linspace(float(r[0]), float(r[-1]), n_grid + 2)
    grid = full_grid[1:-1]
    dr = float(full_grid[1] - full_grid[0])
    kinetic_diag = 1.0 / (mu * dr * dr)
    kinetic_offdiag = -1.0 / (2.0 * mu * dr * dr)
    offdiag = np.full(n_grid - 1, kinetic_offdiag, dtype=float)

    levels_h: Dict[int, Tuple[float, ...]] = {}
    terms_cm: Dict[int, Tuple[float, ...]] = {}
    margins_cm: Dict[int, float] = {}

    raw_levels: Dict[int, np.ndarray] = {}
    for j in j_tuple:
        centrifugal = j * (j + 1.0) / (2.0 * mu * grid * grid)
        effective = np.asarray(interpolator(grid), dtype=float) + centrifugal
        diagonal = kinetic_diag + effective
        values = eigh_tridiagonal(
            diagonal,
            offdiag,
            eigvals_only=True,
            select="i",
            select_range=(0, n_v - 1),
            check_finite=True,
        )
        if values.shape != (n_v,) or not np.all(np.isfinite(values)):
            raise RuntimeError(f"non-finite rovibrational spectrum for J={j}")

        left = float(potential[0] + j * (j + 1.0) / (2.0 * mu * r[0] * r[0]))
        right = float(potential[-1] + j * (j + 1.0) / (2.0 * mu * r[-1] * r[-1]))
        margin = min(left, right) - float(values[-1])
        if margin <= 0.0:
            raise RuntimeError(
                f"sampled BO interval does not confine the requested J={j} levels; boundary margin={margin} Ha"
            )
        raw_levels[j] = values
        margins_cm[j] = margin * HARTREE_TO_WAVENUMBER_CM

    origin = float(raw_levels[0][0])
    for j, values in raw_levels.items():
        levels_h[j] = tuple(float(x) for x in values)
        terms_cm[j] = tuple(float((x - origin) * HARTREE_TO_WAVENUMBER_CM) for x in values)

    return DiatomicRovibrationalResult(
        sample_r_bohr=tuple(float(x) for x in r),
        sample_energy_hartree=tuple(float(x) for x in e),
        reduced_mass_me=mu,
        radial_grid_points=n_grid,
        levels_hartree_by_j=levels_h,
        term_values_cm_by_j=terms_cm,
        boundary_margin_cm_by_j=margins_cm,
    )


def fit_vibrational_dunham(term_values_j0_cm: Sequence[float]) -> Dict[str, float]:
    levels = np.asarray(term_values_j0_cm, dtype=float)
    if levels.ndim != 1 or levels.size < 4 or not np.all(np.isfinite(levels)):
        raise ValueError("at least four finite J=0 vibrational term values are required")
    spacings = np.diff(levels[:4])
    matrix = np.asarray(
        [[1.0, -2.0 * (v + 1.0), 3.0 * v * v + 6.0 * v + 3.25] for v in range(3)],
        dtype=float,
    )
    omega_e, omega_ex_e, omega_e_y_e = np.linalg.solve(matrix, spacings)
    return {
        "omega_e_cm-1": float(omega_e),
        "omega_ex_e_cm-1": float(omega_ex_e),
        "omega_e_y_e_cm-1": float(omega_e_y_e),
        "fundamental_v0_to_v1_cm-1": float(spacings[0]),
        "v1_to_v2_cm-1": float(spacings[1]),
        "v2_to_v3_cm-1": float(spacings[2]),
    }


def fit_rotational_dunham(term_values_cm_by_j: Mapping[int, Sequence[float]]) -> Dict[str, object]:
    if not all(j in term_values_cm_by_j for j in (0, 1, 2)):
        raise ValueError("J=0,1,2 term values are required")
    e0 = np.asarray(term_values_cm_by_j[0], dtype=float)
    e1 = np.asarray(term_values_cm_by_j[1], dtype=float)
    e2 = np.asarray(term_values_cm_by_j[2], dtype=float)
    if min(e0.size, e1.size, e2.size) < 3:
        raise ValueError("at least v=0,1,2 are required for rotational Dunham fitting")

    b_v = []
    d_v = []
    for v in range(3):
        delta1 = float(e1[v] - e0[v])
        delta2 = float(e2[v] - e0[v])
        d = (3.0 * delta1 - delta2) / 24.0
        b = (delta1 + 4.0 * d) / 2.0
        b_v.append(b)
        d_v.append(d)

    x = np.asarray([0.5, 1.5, 2.5], dtype=float)
    matrix = np.column_stack([np.ones(3), -x, x * x])
    b_e, alpha_e, gamma_e = np.linalg.solve(matrix, np.asarray(b_v, dtype=float))
    return {
        "B_v_cm-1": tuple(float(z) for z in b_v),
        "D_v_cm-1": tuple(float(z) for z in d_v),
        "B_e_cm-1": float(b_e),
        "alpha_e_cm-1": float(alpha_e),
        "gamma_e_cm-1": float(gamma_e),
    }


__all__ = [
    "DiatomicRovibrationalResult",
    "solve_diatomic_rovibrational_levels",
    "fit_vibrational_dunham",
    "fit_rotational_dunham",
]
