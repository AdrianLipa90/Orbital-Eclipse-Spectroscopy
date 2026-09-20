"""OES holonomy-response and local identifiability reference layer.

The functions consume an upstream phase-dependent model. They do not introduce
or validate a new physical interaction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class HolonomyResponseError(ValueError):
    pass


def _positive(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise HolonomyResponseError(f"{name} must be finite and > 0")
    return value


def transition_frequency(energy_m: float, energy_n: float, hbar: float) -> float:
    h = _positive("hbar", hbar)
    em = float(energy_m)
    en = float(energy_n)
    if not np.isfinite(em) or not np.isfinite(en):
        raise HolonomyResponseError("energies must be finite")
    return (em - en) / h


def line_frequency_phase_slope(
    loop_response_m: float,
    loop_response_n: float,
    hbar: float,
) -> float:
    """Return d omega_mn / d Phi = -(J_m-J_n)/hbar."""

    h = _positive("hbar", hbar)
    jm = float(loop_response_m)
    jn = float(loop_response_n)
    if not np.isfinite(jm) or not np.isfinite(jn):
        raise HolonomyResponseError("loop responses must be finite")
    return -(jm - jn) / h


def intensity_phase_slope(
    transition_moment: complex,
    transition_moment_derivative: complex,
) -> float:
    """Return d|mu|^2/dPhi = 2 Re(mu* dmu/dPhi)."""

    mu = complex(transition_moment)
    dmu = complex(transition_moment_derivative)
    values = (mu.real, mu.imag, dmu.real, dmu.imag)
    if not np.all(np.isfinite(values)):
        raise HolonomyResponseError("transition moments must be finite")
    return float(2.0 * np.real(np.conjugate(mu) * dmu))


def _as_jacobian(jacobian) -> np.ndarray:
    j = np.asarray(jacobian, dtype=float)
    if j.ndim != 2 or j.shape[0] == 0 or j.shape[1] == 0:
        raise HolonomyResponseError("jacobian must be a non-empty 2D array")
    if not np.all(np.isfinite(j)):
        raise HolonomyResponseError("jacobian must be finite")
    return j


def singular_values(jacobian) -> np.ndarray:
    j = _as_jacobian(jacobian)
    return np.linalg.svd(j, compute_uv=False)


def identifiable_rank(jacobian, *, rtol: float = 1e-10) -> int:
    j = _as_jacobian(jacobian)
    tol_scale = _positive("rtol", rtol)
    s = np.linalg.svd(j, compute_uv=False)
    if s.size == 0:
        return 0
    threshold = tol_scale * max(j.shape) * s[0]
    return int(np.count_nonzero(s > threshold))


def fisher_information(jacobian, covariance) -> np.ndarray:
    """Return J^T Sigma^-1 J, failing closed unless Sigma is SPD."""

    j = _as_jacobian(jacobian)
    sigma = np.asarray(covariance, dtype=float)
    if sigma.shape != (j.shape[0], j.shape[0]):
        raise HolonomyResponseError(
            "covariance must be square with dimension equal to observable count"
        )
    if not np.all(np.isfinite(sigma)):
        raise HolonomyResponseError("covariance must be finite")
    if not np.allclose(sigma, sigma.T, rtol=0.0, atol=1e-12):
        raise HolonomyResponseError("covariance must be symmetric")
    try:
        chol = np.linalg.cholesky(sigma)
    except np.linalg.LinAlgError as exc:
        raise HolonomyResponseError("covariance must be positive definite") from exc

    whitened = np.linalg.solve(chol, j)
    return whitened.T @ whitened


@dataclass(frozen=True)
class IdentifiabilityReport:
    observable_count: int
    phase_count: int
    jacobian_rank: int
    fisher_rank: int
    singular_values: tuple[float, ...]
    fully_locally_identifiable: bool


def identifiability_report(
    jacobian,
    covariance,
    *,
    rtol: float = 1e-10,
) -> IdentifiabilityReport:
    j = _as_jacobian(jacobian)
    f = fisher_information(j, covariance)
    j_rank = identifiable_rank(j, rtol=rtol)
    f_rank = identifiable_rank(f, rtol=rtol)
    s = tuple(float(x) for x in singular_values(j))
    return IdentifiabilityReport(
        observable_count=int(j.shape[0]),
        phase_count=int(j.shape[1]),
        jacobian_rank=j_rank,
        fisher_rank=f_rank,
        singular_values=s,
        fully_locally_identifiable=(j_rank == j.shape[1] and f_rank == j.shape[1]),
    )
