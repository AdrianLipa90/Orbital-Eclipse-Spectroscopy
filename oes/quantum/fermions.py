"""Small, auditable fermionic core for OES-Q1.

The physical qubit register is represented by one qubit per spin orbital.  For
reference simulation we exploit exact particle-number conservation and work in
the fixed-N determinant sector rather than allocating the full 2**n matrix.
"""

from __future__ import annotations

from itertools import combinations
from math import comb
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

PauliString = Tuple[str, ...]
PauliExpansion = Dict[PauliString, complex]


def full_space_dimension(n_spin_orbitals: int) -> int:
    if n_spin_orbitals < 1:
        raise ValueError("n_spin_orbitals must be positive")
    return 1 << n_spin_orbitals


def sector_dimension(n_spin_orbitals: int, n_electrons: int) -> int:
    if not 0 <= n_electrons <= n_spin_orbitals:
        raise ValueError("invalid particle number")
    return comb(n_spin_orbitals, n_electrons)


def determinant_basis(n_spin_orbitals: int, n_electrons: int) -> List[int]:
    """Return bit-string determinants in lexicographic occupied-orbital order."""
    return [sum(1 << p for p in occ) for occ in combinations(range(n_spin_orbitals), n_electrons)]


def occupied(det: int, n_spin_orbitals: int) -> List[int]:
    return [p for p in range(n_spin_orbitals) if (det >> p) & 1]


def _parity_below(det: int, mode: int) -> int:
    return -1 if ((det & ((1 << mode) - 1)).bit_count() & 1) else 1


def annihilate(det: int, mode: int):
    if ((det >> mode) & 1) == 0:
        return None
    return det ^ (1 << mode), _parity_below(det, mode)


def create(det: int, mode: int):
    if ((det >> mode) & 1) == 1:
        return None
    return det | (1 << mode), _parity_below(det, mode)


def _spin_orbital(mode: int) -> Tuple[int, int]:
    """Map spin-orbital mode -> (spatial orbital, spin), spin 0=alpha, 1=beta."""
    return mode // 2, mode & 1


def antisym_spin_eri(eri: np.ndarray, p: int, q: int, r: int, s: int) -> float:
    """Return <pq||rs> from spatial chemists' ERIs (ij|kl).

    Spin-orbital ordering is (spatial0 alpha, spatial0 beta, spatial1 alpha, ...).
    """
    P, sp = _spin_orbital(p)
    Q, sq = _spin_orbital(q)
    R, sr = _spin_orbital(r)
    S, ss = _spin_orbital(s)
    direct = eri[P, R, Q, S] if (sp == sr and sq == ss) else 0.0
    exchange = eri[P, S, Q, R] if (sp == ss and sq == sr) else 0.0
    return float(direct - exchange)


def build_sector_hamiltonian(h1: np.ndarray, eri: np.ndarray, n_electrons: int, ecore: float = 0.0):
    """Build the exact number-conserving Hamiltonian in the determinant sector.

    h1 is a spatial one-electron matrix and eri uses chemists' notation (ij|kl).
    The returned matrix represents

        H = sum_pq h_pq a_p^† a_q
            + 1/4 sum_pqrs <pq||rs> a_p^† a_q^† a_s a_r + ecore.
    """
    h1 = np.asarray(h1, dtype=float)
    eri = np.asarray(eri, dtype=float)
    if h1.ndim != 2 or h1.shape[0] != h1.shape[1]:
        raise ValueError("h1 must be square")
    n_spatial = h1.shape[0]
    if eri.shape != (n_spatial, n_spatial, n_spatial, n_spatial):
        raise ValueError("eri shape must be (n,n,n,n)")
    n_spin = 2 * n_spatial
    basis = determinant_basis(n_spin, n_electrons)
    index = {det: i for i, det in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)), dtype=float)

    for col, det in enumerate(basis):
        H[col, col] += ecore
        occ = occupied(det, n_spin)

        # One-body contribution.
        for q in occ:
            q_spatial, q_spin = _spin_orbital(q)
            first = annihilate(det, q)
            assert first is not None
            d1, s1 = first
            for p in range(n_spin):
                p_spatial, p_spin = _spin_orbital(p)
                if p_spin != q_spin:
                    continue
                hpq = h1[p_spatial, q_spatial]
                if abs(hpq) < 1e-15:
                    continue
                second = create(d1, p)
                if second is None:
                    continue
                d2, s2 = second
                H[index[d2], col] += hpq * s1 * s2

        # Antisymmetrized two-body contribution.  r and s must be occupied.
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
                        H[index[d4], col] += 0.25 * value * sr_sign * ss_sign * cq_sign * cp_sign

    # Numerical construction should be Hermitian; symmetrize only roundoff.
    H = 0.5 * (H + H.T)
    return H, basis


def _transition_one_rdm(final: np.ndarray, initial: np.ndarray, basis: Sequence[int], n_spin: int) -> np.ndarray:
    index = {det: i for i, det in enumerate(basis)}
    gamma = np.zeros((n_spin, n_spin), dtype=complex)
    for col, det in enumerate(basis):
        ci = initial[col]
        if abs(ci) < 1e-15:
            continue
        for q in occupied(det, n_spin):
            a = annihilate(det, q)
            assert a is not None
            d1, s1 = a
            for p in range(n_spin):
                c = create(d1, p)
                if c is None:
                    continue
                d2, s2 = c
                row = index.get(d2)
                if row is None:
                    continue
                gamma[p, q] += np.conjugate(final[row]) * ci * s1 * s2
    return gamma


def transition_one_rdm(final: np.ndarray, initial: np.ndarray, basis: Sequence[int], n_spin_orbitals: int) -> np.ndarray:
    """T_pq = <Psi_f|a_p^† a_q|Psi_i>."""
    return _transition_one_rdm(np.asarray(final), np.asarray(initial), basis, n_spin_orbitals)


def one_rdm(state: np.ndarray, basis: Sequence[int], n_spin_orbitals: int) -> np.ndarray:
    return _transition_one_rdm(np.asarray(state), np.asarray(state), basis, n_spin_orbitals)


def two_rdm(state: np.ndarray, basis: Sequence[int], n_spin_orbitals: int) -> np.ndarray:
    """Gamma_pqrs = <a_p^† a_q^† a_s a_r> in the fixed-N sector."""
    state = np.asarray(state)
    index = {det: i for i, det in enumerate(basis)}
    gamma = np.zeros((n_spin_orbitals,) * 4, dtype=complex)
    for col, det in enumerate(basis):
        ci = state[col]
        if abs(ci) < 1e-15:
            continue
        occ = occupied(det, n_spin_orbitals)
        for r in occ:
            ar = annihilate(det, r)
            assert ar is not None
            d1, s1 = ar
            for s in occ:
                if s == r:
                    continue
                ass = annihilate(d1, s)
                if ass is None:
                    continue
                d2, s2 = ass
                for q in range(n_spin_orbitals):
                    cq = create(d2, q)
                    if cq is None:
                        continue
                    d3, s3 = cq
                    for p in range(n_spin_orbitals):
                        cp = create(d3, p)
                        if cp is None:
                            continue
                        d4, s4 = cp
                        row = index.get(d4)
                        if row is None:
                            continue
                        gamma[p, q, r, s] += np.conjugate(state[row]) * ci * s1 * s2 * s3 * s4
    return gamma


_PAULI_PRODUCT = {
    ("I", "I"): (1, "I"), ("I", "X"): (1, "X"), ("I", "Y"): (1, "Y"), ("I", "Z"): (1, "Z"),
    ("X", "I"): (1, "X"), ("Y", "I"): (1, "Y"), ("Z", "I"): (1, "Z"),
    ("X", "X"): (1, "I"), ("Y", "Y"): (1, "I"), ("Z", "Z"): (1, "I"),
    ("X", "Y"): (1j, "Z"), ("Y", "X"): (-1j, "Z"),
    ("Y", "Z"): (1j, "X"), ("Z", "Y"): (-1j, "X"),
    ("Z", "X"): (1j, "Y"), ("X", "Z"): (-1j, "Y"),
}


def _multiply_pauli_strings(a: PauliString, b: PauliString):
    if len(a) != len(b):
        raise ValueError("Pauli strings must have equal length")
    phase = 1.0 + 0.0j
    out = []
    for x, y in zip(a, b):
        ph, z = _PAULI_PRODUCT[(x, y)]
        phase *= ph
        out.append(z)
    return phase, tuple(out)


def _multiply_expansions(a: PauliExpansion, b: PauliExpansion) -> PauliExpansion:
    out: PauliExpansion = {}
    for pa, ca in a.items():
        for pb, cb in b.items():
            ph, pc = _multiply_pauli_strings(pa, pb)
            out[pc] = out.get(pc, 0.0j) + ca * cb * ph
    return {p: c for p, c in out.items() if abs(c) > 1e-14}


def jw_ladder(n_qubits: int, mode: int, dagger: bool) -> PauliExpansion:
    """Jordan-Wigner image of a_mode or a_mode^†."""
    if not 0 <= mode < n_qubits:
        raise ValueError("mode outside register")
    zprefix = ["Z" if k < mode else "I" for k in range(n_qubits)]
    px = list(zprefix)
    py = list(zprefix)
    px[mode] = "X"
    py[mode] = "Y"
    # a = (X + iY)/2; a^† = (X - iY)/2, with the parity Z prefix.
    return {tuple(px): 0.5, tuple(py): (-0.5j if dagger else 0.5j)}


def jw_product(n_qubits: int, operators: Iterable[Tuple[str, int]]) -> PauliExpansion:
    """Map an ordered fermionic operator product to Pauli strings.

    operators is written left-to-right, e.g. [("create", p), ("annihilate", q)].
    """
    identity = tuple("I" for _ in range(n_qubits))
    out: PauliExpansion = {identity: 1.0 + 0.0j}
    for kind, mode in operators:
        if kind not in {"create", "annihilate"}:
            raise ValueError("operator kind must be create or annihilate")
        out = _multiply_expansions(out, jw_ladder(n_qubits, mode, dagger=(kind == "create")))
    return out
