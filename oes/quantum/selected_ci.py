"""Selected external-determinant dressing for two-electron OES active spaces.

This module keeps the quantum P-space fixed and uses Hamiltonian coupling only
to rank a classical Q-space bath. Selected Q determinants are then diagonalized
together with P, so Q-Q couplings and repeated P<->Q excursions are retained.
No experimental energy enters the selection.
"""

from __future__ import annotations

from typing import Dict, List, Mapping, Sequence

import numpy as np

from .external_dressing import ExternalCouplingSpace
from .fermions import annihilate, antisym_spin_eri, create, occupied


def two_electron_determinant_matrix_element(
    h1: np.ndarray,
    eri: np.ndarray,
    bra: int,
    ket: int,
    n_spin_orbitals: int,
    ecore: float = 0.0,
) -> float:
    """Return <bra|H|ket> for exactly two electrons.

    The implementation uses the same ordered fermionic operators and ERI
    convention as ``build_sector_hamiltonian`` but exploits the fact that both
    determinants contain only two occupied spin orbitals. This makes arbitrary
    selected determinant subspaces cheap to assemble without constructing the
    full C(n_spin, 2) Hamiltonian.
    """
    h1 = np.asarray(h1, dtype=float)
    eri = np.asarray(eri, dtype=float)
    n_spatial = h1.shape[0]
    if h1.shape != (n_spatial, n_spatial):
        raise ValueError("h1 must be square")
    if eri.shape != (n_spatial,) * 4:
        raise ValueError("eri shape mismatch")
    if n_spin_orbitals != 2 * n_spatial:
        raise ValueError("n_spin_orbitals incompatible with spatial integrals")

    occ_bra = occupied(int(bra), n_spin_orbitals)
    occ_ket = occupied(int(ket), n_spin_orbitals)
    if len(occ_bra) != 2 or len(occ_ket) != 2:
        raise ValueError("selected-CI primitive requires exactly two-electron determinants")

    value = float(ecore) if bra == ket else 0.0

    # One-body term sum_pq h_pq a_p^+ a_q. Only p occupied in the target bra
    # can possibly produce the requested determinant.
    for q in occ_ket:
        first = annihilate(ket, q)
        assert first is not None
        d1, s1 = first
        q_spatial, q_spin = q // 2, q & 1
        for p in occ_bra:
            if (p & 1) != q_spin:
                continue
            second = create(d1, p)
            if second is None:
                continue
            d2, s2 = second
            if d2 == bra:
                value += float(h1[p // 2, q_spatial]) * s1 * s2

    # Two-body term 1/4 sum_pqrs <pq||rs> a_p^+ a_q^+ a_s a_r.
    # For N=2, annihilating r,s empties the determinant, and p,q must be the
    # two orbitals occupied by the target bra.
    for r in occ_ket:
        ar = annihilate(ket, r)
        assert ar is not None
        d1, sr_sign = ar
        for s in occ_ket:
            if s == r:
                continue
            ass = annihilate(d1, s)
            if ass is None:
                continue
            d2, ss_sign = ass
            for q in occ_bra:
                cq = create(d2, q)
                if cq is None:
                    continue
                d3, cq_sign = cq
                for p in occ_bra:
                    if p == q:
                        continue
                    cp = create(d3, p)
                    if cp is None:
                        continue
                    d4, cp_sign = cp
                    if d4 != bra:
                        continue
                    g = antisym_spin_eri(eri, p, q, r, s)
                    value += 0.25 * g * sr_sign * ss_sign * cq_sign * cp_sign

    return float(value)


def build_two_electron_subspace_hamiltonian(
    h1: np.ndarray,
    eri: np.ndarray,
    determinants: Sequence[int],
    ecore: float = 0.0,
) -> np.ndarray:
    """Build H on an arbitrary ordered subset of two-electron determinants."""
    h1 = np.asarray(h1, dtype=float)
    n_spin = 2 * h1.shape[0]
    basis = tuple(int(x) for x in determinants)
    if len(set(basis)) != len(basis):
        raise ValueError("determinant subspace contains duplicates")
    H = np.zeros((len(basis), len(basis)), dtype=float)
    for i, bra in enumerate(basis):
        for j in range(i + 1):
            value = two_electron_determinant_matrix_element(
                h1,
                eri,
                bra,
                basis[j],
                n_spin,
                ecore=ecore,
            )
            H[i, j] = value
            H[j, i] = value
    return H


def state_balanced_external_importance(
    external: ExternalCouplingSpace,
    class_states: Mapping[str, np.ndarray],
    class_energies_hartree: Mapping[str, float],
    denominator_floor: float = 1e-5,
) -> Dict[str, object]:
    """Rank external determinants from equal-weight active-state classes.

    Each class contributes a normalized positive importance distribution

        I_a(class) ~ mean_i |<a|H|Psi_i>|^2 / |E_class - H_aa|,

    so a strongly coupled ground state cannot numerically drown a weaker excited
    class. Multi-component degenerate manifolds are averaged before class
    normalization. The final score is the equal-weight mean over classes.
    """
    names = tuple(class_states)
    if not names or set(names) != set(class_energies_hartree):
        raise ValueError("state classes and class energies must have matching non-empty keys")

    total = np.zeros(len(external.external_basis), dtype=float)
    per_class: Dict[str, np.ndarray] = {}
    for name in names:
        states = np.asarray(class_states[name], dtype=complex)
        if states.ndim == 1:
            states = states[:, None]
        if states.ndim != 2 or states.shape[0] != external.coupling_qp.shape[1]:
            raise ValueError(f"class {name} has incompatible P-space dimension")
        gram = states.conj().T @ states
        if not np.allclose(gram, np.eye(states.shape[1]), atol=1e-10):
            raise ValueError(f"class {name} states must be orthonormal")
        denom = np.abs(float(class_energies_hartree[name]) - external.diagonal_hartree)
        if float(np.min(denom)) < denominator_floor:
            raise RuntimeError(f"selected-Q importance hit intruder denominator in class {name}")
        couplings = external.coupling_qp @ states
        raw = np.mean(np.abs(couplings) ** 2, axis=1) / denom
        norm = float(np.sum(raw))
        if norm <= 0.0:
            raise RuntimeError(f"class {name} has zero external importance")
        normalized = raw / norm
        per_class[name] = normalized
        total += normalized / len(names)

    return {
        "scores": total,
        "per_class_scores": per_class,
        "class_names": names,
    }


def grouped_importance_order(
    scores: np.ndarray,
    relative_tie_tolerance: float = 1e-7,
    absolute_tie_tolerance: float = 1e-15,
) -> List[List[int]]:
    """Sort external indices by importance while keeping numerical ties intact."""
    scores = np.asarray(scores, dtype=float)
    if scores.ndim != 1 or np.any(scores < 0) or not np.all(np.isfinite(scores)):
        raise ValueError("scores must be a finite non-negative vector")
    order = np.argsort(scores)[::-1]
    if len(order) == 0:
        return []
    groups: List[List[int]] = []
    current = [int(order[0])]
    reference = float(scores[order[0]])
    for idx_raw in order[1:]:
        idx = int(idx_raw)
        value = float(scores[idx])
        tol = max(
            absolute_tie_tolerance,
            relative_tie_tolerance * max(abs(reference), abs(value)),
        )
        if abs(value - reference) <= tol:
            current.append(idx)
        else:
            groups.append(current)
            current = [idx]
            reference = value
    groups.append(current)
    return groups


def grouped_prefix_for_target(groups: Sequence[Sequence[int]], target: int) -> List[int]:
    """Return complete importance groups until at least ``target`` indices exist."""
    if target < 1:
        raise ValueError("target must be positive")
    out: List[int] = []
    for group in groups:
        out.extend(int(x) for x in group)
        if len(out) >= target:
            break
    if len(out) < target:
        raise ValueError("target exceeds available grouped external space")
    return out
