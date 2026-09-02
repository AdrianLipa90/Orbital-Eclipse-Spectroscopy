"""External-space energy dressing for fixed-particle OES active spaces.

The active P-space wavefunction stays inside the fixed qubit register. A larger
one-particle source basis is partitioned into P (the first active orbitals) and
Q (all determinants containing at least one external spin orbital). We build
only the Q<-P Hamiltonian coupling and Q diagonal, never the full dense FCI
matrix.

For an isolated state the initial diagnostic is unshifted Epstein-Nesbet second
order

    dE_I^(2) = sum_a |<a|H|Psi_I>|^2 / (E_I - H_aa),  a in Q.

For an exactly or quasi-degenerate P-space manifold, state-by-state EN2 is not
basis invariant. The corresponding diagnostic therefore uses the full
second-order effective block

    W_ij^(2)(E_ref) = sum_a <Psi_i|H|a><a|H|Psi_j> / (E_ref - H_aa)

and diagonalizes W inside the degenerate manifold. This distinguishes genuine
symmetry breaking by the diagonal-Q approximation from arbitrary rotations of
the active eigenvectors.

No fitted shift or experimental energy enters. Small denominators are reported
and fail closed rather than silently regularized.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .fermions import (
    annihilate,
    antisym_spin_eri,
    create,
    determinant_basis,
    occupied,
)


@dataclass(frozen=True)
class ExternalCouplingSpace:
    external_basis: Tuple[int, ...]
    diagonal_hartree: np.ndarray
    coupling_qp: np.ndarray
    n_full_spin_orbitals: int
    n_active_spin_orbitals: int


def complete_mo_basis_preserving_active(mf, C_active: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Return a complete AO coefficient matrix whose first columns are C_active.

    The canonical RHF MO basis is S-orthonormal. Active vectors are converted
    to that coordinate system and kept exactly; projected canonical unit vectors
    deterministically fill the orthogonal complement.
    """
    C = np.asarray(mf.mo_coeff, dtype=float)
    S = np.asarray(mf.get_ovlp(), dtype=float)
    C_active = np.asarray(C_active, dtype=float)
    U_active = C.T @ S @ C_active
    n_full = C.shape[1]
    n_active = C_active.shape[1]
    if U_active.shape != (n_full, n_active):
        raise ValueError("active coefficient shape incompatible with source MOs")
    if not np.allclose(U_active.T @ U_active, np.eye(n_active), atol=1e-9):
        raise RuntimeError("active MO-coordinate vectors are not orthonormal")

    vectors: List[np.ndarray] = [U_active[:, i].copy() for i in range(n_active)]
    for idx in range(n_full):
        if len(vectors) >= n_full:
            break
        v = np.zeros(n_full, dtype=float)
        v[idx] = 1.0
        for q in vectors:
            v -= q * float(np.dot(q, v))
        for q in vectors:  # second pass
            v -= q * float(np.dot(q, v))
        norm = float(np.linalg.norm(v))
        if norm > tol:
            vectors.append(v / norm)
    if len(vectors) != n_full:
        raise RuntimeError(f"could not complete active basis: {len(vectors)} of {n_full}")
    U = np.column_stack(vectors)
    if not np.allclose(U.T @ U, np.eye(n_full), atol=1e-9):
        raise RuntimeError("completed source basis lost orthonormality")
    if not np.allclose(U[:, :n_active], U_active, atol=1e-9):
        raise RuntimeError("basis completion rotated the active subspace columns")
    return C @ U


def determinant_diagonal(h1: np.ndarray, eri: np.ndarray, det: int, n_spin: int, ecore: float = 0.0) -> float:
    """Slater determinant diagonal <D|H|D>."""
    occ = occupied(det, n_spin)
    value = float(ecore)
    for p in occ:
        P = p // 2
        value += float(h1[P, P])
    for i, p in enumerate(occ):
        for q in occ[i + 1 :]:
            value += antisym_spin_eri(eri, p, q, p, q)
    return float(value)


def build_external_coupling_space(
    h1: np.ndarray,
    eri: np.ndarray,
    active_basis: Sequence[int],
    n_active_spatial: int,
    n_electrons: int = 2,
    ecore: float = 0.0,
) -> ExternalCouplingSpace:
    """Build dense Q<-P coupling for a small active sector and a large source basis."""
    h1 = np.asarray(h1, dtype=float)
    eri = np.asarray(eri, dtype=float)
    n_full_spatial = h1.shape[0]
    if h1.shape != (n_full_spatial, n_full_spatial):
        raise ValueError("h1 must be square")
    if eri.shape != (n_full_spatial,) * 4:
        raise ValueError("eri shape mismatch")
    n_full_spin = 2 * n_full_spatial
    n_active_spin = 2 * n_active_spatial
    if n_electrons != 2:
        raise ValueError("Q1 external dressing primitive is currently validated for two electrons only")

    active_basis = tuple(int(x) for x in active_basis)
    active_set = set(active_basis)
    expected_active = set(determinant_basis(n_active_spin, n_electrons))
    if active_set != expected_active:
        raise ValueError("active basis must be the complete fixed-N determinant basis of the active register")

    full_basis = determinant_basis(n_full_spin, n_electrons)
    external_basis = tuple(det for det in full_basis if det not in active_set)
    qindex = {det: i for i, det in enumerate(external_basis)}
    V = np.zeros((len(external_basis), len(active_basis)), dtype=float)

    # Apply H to each active determinant and retain only Q-space outputs.
    for col, det in enumerate(active_basis):
        occ = occupied(det, n_full_spin)

        # One-body contribution.
        for q in occ:
            q_spatial, q_spin = q // 2, q & 1
            first = annihilate(det, q)
            assert first is not None
            d1, s1 = first
            for p in range(n_full_spin):
                if (p & 1) != q_spin:
                    continue
                p_spatial = p // 2
                hpq = float(h1[p_spatial, q_spatial])
                if abs(hpq) < 1e-15:
                    continue
                second = create(d1, p)
                if second is None:
                    continue
                d2, s2 = second
                row = qindex.get(d2)
                if row is not None:
                    V[row, col] += hpq * s1 * s2

        # Two-body contribution in the same convention as build_sector_hamiltonian.
        for r in occ:
            ar = annihilate(det, r)
            assert ar is not None
            d1, sr_sign = ar
            for s in occ:
                if s == r:
                    continue
                ass = annihilate(d1, s)
                if ass is None:
                    continue
                d2, ss_sign = ass
                for q in range(n_full_spin):
                    cq = create(d2, q)
                    if cq is None:
                        continue
                    d3, cq_sign = cq
                    for p in range(n_full_spin):
                        value = antisym_spin_eri(eri, p, q, r, s)
                        if abs(value) < 1e-15:
                            continue
                        cp = create(d3, p)
                        if cp is None:
                            continue
                        d4, cp_sign = cp
                        row = qindex.get(d4)
                        if row is not None:
                            V[row, col] += 0.25 * value * sr_sign * ss_sign * cq_sign * cp_sign

    diagonal = np.array(
        [determinant_diagonal(h1, eri, det, n_full_spin, ecore=ecore) for det in external_basis],
        dtype=float,
    )
    return ExternalCouplingSpace(
        external_basis=external_basis,
        diagonal_hartree=diagonal,
        coupling_qp=V,
        n_full_spin_orbitals=n_full_spin,
        n_active_spin_orbitals=n_active_spin,
    )


def _checked_denominators(
    energy_hartree: float,
    external: ExternalCouplingSpace,
    denominator_floor: float,
) -> tuple[np.ndarray, float]:
    denominators = float(energy_hartree) - external.diagonal_hartree
    abs_den = np.abs(denominators)
    min_abs = float(np.min(abs_den))
    near = int(np.count_nonzero(abs_den < denominator_floor))
    if near:
        raise RuntimeError(
            f"EN2 intruder gate failed: {near} denominators below {denominator_floor} Ha; min={min_abs} Ha"
        )
    return denominators, min_abs


def en2_correction(
    energy_hartree: float,
    state: np.ndarray,
    external: ExternalCouplingSpace,
    denominator_floor: float = 1e-5,
) -> Dict[str, float]:
    """Return unshifted state-specific Epstein-Nesbet second-order correction."""
    state = np.asarray(state, dtype=complex)
    if state.ndim != 1 or state.shape[0] != external.coupling_qp.shape[1]:
        raise ValueError("state dimension incompatible with active P space")
    v = external.coupling_qp @ state
    denominators, min_abs = _checked_denominators(energy_hartree, external, denominator_floor)
    weights = np.abs(v) ** 2
    correction = float(np.sum(weights / denominators).real)
    return {
        "correction_hartree": correction,
        "coupling_norm2_hartree2": float(np.sum(weights).real),
        "min_abs_denominator_hartree": min_abs,
        "max_abs_term_hartree": float(np.max(np.abs(weights / denominators))),
        "external_determinants": int(len(external.external_basis)),
    }


def en2_degenerate_block(
    reference_energy_hartree: float,
    states: np.ndarray,
    external: ExternalCouplingSpace,
    denominator_floor: float = 1e-5,
) -> Dict[str, object]:
    """Return basis-invariant second-order effective Hamiltonian for a P manifold.

    ``states`` has shape (dim(P), n_states) and its columns must be orthonormal.
    A single shared reference energy is required because this primitive is for an
    exactly/quasi-degenerate manifold. The eigenvalues of the returned W block
    are the second-order corrections in the optimally rotated manifold basis.
    """
    states = np.asarray(states, dtype=complex)
    if states.ndim != 2 or states.shape[0] != external.coupling_qp.shape[1]:
        raise ValueError("state block dimension incompatible with active P space")
    n_states = states.shape[1]
    if n_states < 2:
        raise ValueError("degenerate block requires at least two states")
    gram = states.conj().T @ states
    if not np.allclose(gram, np.eye(n_states), atol=1e-10):
        raise ValueError("degenerate block states must be orthonormal")

    v = external.coupling_qp @ states
    denominators, min_abs = _checked_denominators(
        reference_energy_hartree,
        external,
        denominator_floor,
    )
    W = v.conj().T @ (v / denominators[:, None])
    W = 0.5 * (W + W.conj().T)
    correction_eigenvalues = np.linalg.eigvalsh(W).real
    mean_correction = float(np.trace(W).real / n_states)
    spread = float(np.max(correction_eigenvalues) - np.min(correction_eigenvalues))

    return {
        "correction_matrix_hartree": [[float(x.real) for x in row] for row in W],
        "correction_eigenvalues_hartree": [float(x) for x in correction_eigenvalues],
        "mean_correction_hartree": mean_correction,
        "correction_spread_hartree": spread,
        "coupling_frobenius_norm2_hartree2": float(np.sum(np.abs(v) ** 2).real),
        "min_abs_denominator_hartree": min_abs,
        "external_determinants": int(len(external.external_basis)),
    }
