"""Exact Feshbach/Schur downfolding into a fixed OES P-space.

For a Hermitian selected-space Hamiltonian

    H = [[H_PP, H_PQ],
         [H_QP, H_QQ]],

the Q amplitudes can be eliminated at energy E to give the exact energy-
dependent P-space operator

    H_eff(E) = H_PP + H_PQ (E I - H_QQ)^(-1) H_QP.

No perturbative truncation is made.  If E is an eigenvalue of the full selected
P+Q problem and is outside the spectrum of H_QQ, its P projection is an exact
eigenvector of H_eff(E).  In Q1 the P-space is the fixed N=2 sector of the 20Q
register (dimension 190); the classical bath is therefore integrated out rather
than appended to the quantum register.
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np


def partition_hamiltonian(H: np.ndarray, p_dim: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return H_PP, H_PQ and H_QQ from a P-first Hermitian matrix."""
    H = np.asarray(H, dtype=float)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError("H must be square")
    if not np.allclose(H, H.T, atol=1e-11):
        raise ValueError("H must be Hermitian/real-symmetric")
    if p_dim < 1 or p_dim >= H.shape[0]:
        raise ValueError("p_dim must leave non-empty P and Q blocks")
    return H[:p_dim, :p_dim], H[:p_dim, p_dim:], H[p_dim:, p_dim:]


def effective_hamiltonian(
    H: np.ndarray,
    p_dim: int,
    energy_hartree: float,
    singular_floor: float = 1e-10,
) -> Tuple[np.ndarray, Dict[str, float]]:
    """Return the exact energy-dependent Feshbach Hamiltonian on P.

    ``singular_floor`` is a fail-closed gate on the distance from E to the QHQ
    spectrum.  No level shift or fitted regularizer is introduced.
    """
    Hpp, Hpq, Hqq = partition_hamiltonian(H, p_dim)
    qevals = np.linalg.eigvalsh(Hqq)
    distance = float(np.min(np.abs(float(energy_hartree) - qevals)))
    if distance < singular_floor:
        raise RuntimeError(
            f"Feshbach resolvent gate failed: distance to QHQ spectrum {distance} Ha "
            f"below {singular_floor} Ha"
        )
    resolvent_rhs = Hpq.T
    solved = np.linalg.solve(float(energy_hartree) * np.eye(Hqq.shape[0]) - Hqq, resolvent_rhs)
    sigma = Hpq @ solved
    sigma = 0.5 * (sigma + sigma.T)
    heff = Hpp + sigma
    return heff, {
        "p_dimension": int(p_dim),
        "q_dimension": int(Hqq.shape[0]),
        "distance_to_qhq_spectrum_hartree": distance,
        "self_energy_frobenius_hartree": float(np.linalg.norm(sigma)),
        "self_energy_spectral_hartree": float(np.linalg.norm(sigma, ord=2)),
    }


def eigenpair_downfolding_residual(
    H: np.ndarray,
    p_dim: int,
    eigenvalue_hartree: float,
    eigenvector: np.ndarray,
    singular_floor: float = 1e-10,
) -> Dict[str, float]:
    """Validate one full-space eigenpair against the exact P-space reduction."""
    H = np.asarray(H, dtype=float)
    vector = np.asarray(eigenvector, dtype=complex)
    if vector.ndim != 1 or vector.shape[0] != H.shape[0]:
        raise ValueError("eigenvector dimension mismatch")
    p = vector[:p_dim]
    q = vector[p_dim:]
    pnorm = float(np.linalg.norm(p))
    if pnorm < 1e-12:
        raise RuntimeError("eigenpair has negligible P-space weight")

    heff, diagnostics = effective_hamiltonian(
        H,
        p_dim,
        eigenvalue_hartree,
        singular_floor=singular_floor,
    )
    residual = (heff - float(eigenvalue_hartree) * np.eye(p_dim)) @ p
    residual_norm = float(np.linalg.norm(residual) / pnorm)

    _, Hpq, Hqq = partition_hamiltonian(H, p_dim)
    q_reconstructed = np.linalg.solve(
        float(eigenvalue_hartree) * np.eye(Hqq.shape[0]) - Hqq,
        Hpq.T @ p,
    )
    q_error = float(np.linalg.norm(q_reconstructed - q))
    full_norm = float(np.linalg.norm(vector))

    return {
        **diagnostics,
        "p_weight": float(np.vdot(p, p).real / (full_norm * full_norm)),
        "q_weight": float(np.vdot(q, q).real / (full_norm * full_norm)),
        "effective_eigen_residual_hartree": residual_norm,
        "q_reconstruction_error": q_error,
    }
