"""Sparse exact Hamiltonians in fixed (N_alpha, N_beta) sectors.

The register remains one qubit per spin orbital.  This module changes only the
classical reference representation: instead of allocating a dense determinant
matrix, it enumerates the Slater-Condon-connected determinants inside an exact
spin-projection sector and stores the result in CSR form.

Spin-orbital ordering is inherited from ``fermions``:
    (spatial0 alpha, spatial0 beta, spatial1 alpha, spatial1 beta, ...).
"""

from __future__ import annotations

from itertools import combinations
from math import comb
from typing import Iterable, Sequence, Tuple

import numpy as np

from .fermions import annihilate, antisym_spin_eri, create, occupied


def fixed_spin_determinant_basis(
    n_spatial: int,
    n_alpha: int,
    n_beta: int,
) -> Tuple[int, ...]:
    """Return determinants with exact alpha/beta particle counts.

    The ordering is deterministic: alpha spatial combinations are outermost,
    beta combinations innermost.  No full C(2*n, N) sector is materialized.
    """
    if n_spatial < 1:
        raise ValueError("n_spatial must be positive")
    if not (0 <= n_alpha <= n_spatial and 0 <= n_beta <= n_spatial):
        raise ValueError("invalid alpha/beta particle count")
    basis = []
    for alpha_occ in combinations(range(n_spatial), n_alpha):
        alpha_bits = sum(1 << (2 * p) for p in alpha_occ)
        for beta_occ in combinations(range(n_spatial), n_beta):
            beta_bits = sum(1 << (2 * p + 1) for p in beta_occ)
            basis.append(alpha_bits | beta_bits)
    expected = comb(n_spatial, n_alpha) * comb(n_spatial, n_beta)
    if len(basis) != expected:
        raise RuntimeError("fixed-spin determinant enumeration failed")
    return tuple(basis)


def _single_phase(det: int, removed: int, added: int):
    first = annihilate(det, removed)
    if first is None:
        return None
    d1, s1 = first
    second = create(d1, added)
    if second is None:
        return None
    d2, s2 = second
    return d2, s1 * s2


def _double_phase(det: int, removed_a: int, removed_b: int, added_a: int, added_b: int):
    """Apply a_added_a^+ a_added_b^+ a_removed_b a_removed_a.

    This matches the operator order used by the antisymmetrized two-body
    Hamiltonian when p=added_a, q=added_b, r=removed_a, s=removed_b.
    """
    first = annihilate(det, removed_a)
    if first is None:
        return None
    d1, s1 = first
    second = annihilate(d1, removed_b)
    if second is None:
        return None
    d2, s2 = second
    third = create(d2, added_b)
    if third is None:
        return None
    d3, s3 = third
    fourth = create(d3, added_a)
    if fourth is None:
        return None
    d4, s4 = fourth
    return d4, s1 * s2 * s3 * s4


def _same_spin_pairs(modes: Sequence[int]) -> Iterable[Tuple[int, int]]:
    return combinations(modes, 2)


def _mixed_spin_pairs(alpha_modes: Sequence[int], beta_modes: Sequence[int]) -> Iterable[Tuple[int, int]]:
    for a in alpha_modes:
        for b in beta_modes:
            yield (a, b) if a < b else (b, a)


def build_sparse_fixed_spin_hamiltonian(
    h1: np.ndarray,
    eri: np.ndarray,
    n_alpha: int,
    n_beta: int,
    ecore: float = 0.0,
    zero_tolerance: float = 1e-15,
):
    """Build the exact fixed-(N_alpha,N_beta) Hamiltonian as CSR.

    Slater-Condon connectivity is used explicitly:
      diagonal: sum_i h_ii + sum_{i<j} <ij||ij>,
      single i->a: phase * [h_ai + sum_j <aj||ij>],
      double ij->ab: phase * <ab||ij>.

    The returned operator is algebraically the same fermionic Hamiltonian as
    ``build_determinant_subspace_hamiltonian``.  Only its storage and traversal
    differ.
    """
    try:
        from scipy.sparse import coo_matrix
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("sparse fixed-spin Hamiltonians require SciPy") from exc

    h1 = np.asarray(h1, dtype=float)
    eri = np.asarray(eri, dtype=float)
    if h1.ndim != 2 or h1.shape[0] != h1.shape[1]:
        raise ValueError("h1 must be square")
    n_spatial = h1.shape[0]
    if eri.shape != (n_spatial,) * 4:
        raise ValueError("eri shape mismatch")
    if zero_tolerance < 0:
        raise ValueError("zero_tolerance must be non-negative")
    n_spin = 2 * n_spatial
    basis = fixed_spin_determinant_basis(n_spatial, n_alpha, n_beta)
    index = {det: i for i, det in enumerate(basis)}

    rows = []
    cols = []
    data = []

    all_alpha = tuple(2 * p for p in range(n_spatial))
    all_beta = tuple(2 * p + 1 for p in range(n_spatial))

    for col, det in enumerate(basis):
        occ = occupied(det, n_spin)
        occ_set = set(occ)
        occ_alpha = [p for p in occ if (p & 1) == 0]
        occ_beta = [p for p in occ if (p & 1) == 1]
        vir_alpha = [p for p in all_alpha if p not in occ_set]
        vir_beta = [p for p in all_beta if p not in occ_set]

        diagonal = float(ecore)
        for i in occ:
            diagonal += float(h1[i // 2, i // 2])
        for ix, i in enumerate(occ):
            for j in occ[ix + 1 :]:
                diagonal += antisym_spin_eri(eri, i, j, i, j)
        rows.append(col)
        cols.append(col)
        data.append(diagonal)

        # Singles conserve spin.  Only insert one triangle; mirror explicitly.
        for i, virtuals in ((i, vir_alpha if (i & 1) == 0 else vir_beta) for i in occ):
            for a in virtuals:
                phased = _single_phase(det, i, a)
                assert phased is not None
                target, phase = phased
                row = index.get(target)
                if row is None or row <= col:
                    continue
                value = float(h1[a // 2, i // 2])
                for j in occ:
                    if j == i:
                        continue
                    value += antisym_spin_eri(eri, a, j, i, j)
                value *= phase
                if abs(value) <= zero_tolerance:
                    continue
                rows.extend((row, col))
                cols.extend((col, row))
                data.extend((value, value))

        removed_groups = (
            tuple(_same_spin_pairs(occ_alpha)),
            tuple(_same_spin_pairs(occ_beta)),
            tuple(_mixed_spin_pairs(occ_alpha, occ_beta)),
        )
        added_groups = (
            tuple(_same_spin_pairs(vir_alpha)),
            tuple(_same_spin_pairs(vir_beta)),
            tuple(_mixed_spin_pairs(vir_alpha, vir_beta)),
        )
        for removed_pairs, added_pairs in zip(removed_groups, added_groups):
            for i, j in removed_pairs:
                for a, b in added_pairs:
                    phased = _double_phase(det, i, j, a, b)
                    assert phased is not None
                    target, phase = phased
                    row = index.get(target)
                    if row is None or row <= col:
                        continue
                    value = phase * antisym_spin_eri(eri, a, b, i, j)
                    if abs(value) <= zero_tolerance:
                        continue
                    rows.extend((row, col))
                    cols.extend((col, row))
                    data.extend((value, value))

    matrix = coo_matrix((np.asarray(data), (np.asarray(rows), np.asarray(cols))), shape=(len(basis), len(basis))).tocsr()
    matrix.sum_duplicates()
    return matrix, basis
