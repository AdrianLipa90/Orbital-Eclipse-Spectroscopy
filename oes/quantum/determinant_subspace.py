"""Exact Hamiltonian construction on a symmetry-selected determinant subset.

The orbital/register encoding remains unchanged.  This helper only avoids
materializing determinant sectors excluded by an exactly conserved discrete
quantum number such as N_alpha/N_beta (equivalently M_S).

The operator convention is identical to ``fermions.build_sector_hamiltonian``:

    H = sum_pq h_pq a_p^† a_q
      + 1/4 sum_pqrs <pq||rs> a_p^† a_q^† a_s a_r + E_core.

States generated outside the supplied determinant set are simply absent because
the intended use is a Hamiltonian-invariant symmetry sector.  A caller should
therefore construct the subset from a conserved quantum number, not arbitrary
truncation.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

from .fermions import annihilate, antisym_spin_eri, create, occupied


def build_determinant_subspace_hamiltonian(
    h1: np.ndarray,
    eri: np.ndarray,
    determinants: Sequence[int],
    ecore: float = 0.0,
) -> Tuple[np.ndarray, Tuple[int, ...]]:
    h1 = np.asarray(h1, dtype=float)
    eri = np.asarray(eri, dtype=float)
    if h1.ndim != 2 or h1.shape[0] != h1.shape[1]:
        raise ValueError("h1 must be square")
    n_spatial = h1.shape[0]
    if eri.shape != (n_spatial,) * 4:
        raise ValueError("eri shape mismatch")
    n_spin = 2 * n_spatial

    basis = tuple(int(x) for x in determinants)
    if not basis or len(set(basis)) != len(basis):
        raise ValueError("determinant subset must be non-empty and unique")
    electron_counts = {det.bit_count() for det in basis}
    if len(electron_counts) != 1:
        raise ValueError("all determinants must have the same particle number")
    if any(det < 0 or det >= (1 << n_spin) for det in basis):
        raise ValueError("determinant outside spin-orbital register")

    index = {det: i for i, det in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)), dtype=float)

    for col, det in enumerate(basis):
        H[col, col] += float(ecore)
        occ = occupied(det, n_spin)

        for q in occ:
            first = annihilate(det, q)
            assert first is not None
            d1, s1 = first
            q_spatial, q_spin = q // 2, q & 1
            for p in range(n_spin):
                if (p & 1) != q_spin:
                    continue
                hpq = float(h1[p // 2, q_spatial])
                if abs(hpq) < 1e-15:
                    continue
                second = create(d1, p)
                if second is None:
                    continue
                d2, s2 = second
                row = index.get(d2)
                if row is not None:
                    H[row, col] += hpq * s1 * s2

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
                for q in range(n_spin):
                    cq = create(d2, q)
                    if cq is None:
                        continue
                    d3, cq_sign = cq
                    for p in range(n_spin):
                        value = antisym_spin_eri(eri, p, q, r, s)
                        if abs(value) < 1e-15:
                            continue
                        cp = create(d3, p)
                        if cp is None:
                            continue
                        d4, cp_sign = cp
                        row = index.get(d4)
                        if row is not None:
                            H[row, col] += 0.25 * value * sr_sign * ss_sign * cq_sign * cp_sign

    H = 0.5 * (H + H.T)
    return H, basis
