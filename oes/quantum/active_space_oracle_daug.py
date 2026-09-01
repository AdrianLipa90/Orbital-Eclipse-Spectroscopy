"""Symmetry-aware d-aug oracle capacity bound for the fixed OES-Q1 20Q budget.

This module is deliberately NONPREDICTIVE. It uses full-source FCI states from
four spectral classes — ground 1S, lowest triplet 3S, lowest dark excited 1S,
and the complete three-component first bright 1P manifold — to form a
state-averaged one-body density matrix. The resulting natural-orbital spectrum
is truncated to exactly 10 spatial orbitals only through complete degenerate
occupation blocks.

The question is capacity, not prediction: can a 20-spin-orbital register retain
these source states when given near-ideal information about which one-particle
subspace matters?
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

HARTREE_TO_EV = 27.211_386_245_981


@dataclass(frozen=True)
class DAugOracleReceipt:
    protocol: str
    source_basis: str
    source_spatial_orbitals: int
    target_spatial_orbitals: int
    target_spin_orbitals: int
    compression_backend: str
    source_ground_hartree: float
    source_triplet_excitation_ev: float
    source_dark_singlet_excitation_ev: float
    source_bright_manifold_excitation_ev: float
    source_bright_manifold_oscillator_strength_sum: float
    occupation_group_sizes: Tuple[int, ...]
    occupation_group_values: Tuple[float, ...]
    selected_group_indices: Tuple[int, ...]
    selected_group_sizes: Tuple[int, ...]
    selected_occupation_weight: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _occupation_groups(values: np.ndarray, vectors: np.ndarray, atol: float = 1e-8):
    """Group adjacent natural occupations that are numerically degenerate."""
    groups = []
    start = 0
    n = len(values)
    while start < n:
        ref = float(values[start])
        end = start + 1
        while end < n and abs(float(values[end]) - ref) <= max(atol, atol * abs(ref)):
            end += 1
        groups.append((start, end, float(np.mean(values[start:end])), vectors[:, start:end]))
        start = end
    return groups


def _best_complete_group_subset(groups, target: int) -> Tuple[int, ...]:
    """Exact-size knapsack maximizing retained occupation weight by full groups."""
    # dp[count] = (weight, tuple(group indices))
    dp = {0: (0.0, tuple())}
    for gi, (start, end, occ, _vecs) in enumerate(groups):
        size = end - start
        weight = occ * size
        updates = dict(dp)
        for count, (score, selected) in dp.items():
            new_count = count + size
            if new_count > target:
                continue
            candidate = (score + weight, selected + (gi,))
            if new_count not in updates or candidate[0] > updates[new_count][0]:
                updates[new_count] = candidate
        dp = updates
    if target not in dp:
        sizes = [end - start for start, end, _occ, _vecs in groups[:12]]
        raise RuntimeError(f"no complete-degeneracy group subset sums to {target}; leading sizes={sizes}")
    return dp[target][1]


def build_helium_daug_oracle_20q(
    source_basis: Any,
    source_label: str,
    target_spatial: int = 10,
    singlet_roots: int = 8,
    bright_threshold: float = 1e-5,
):
    if target_spatial != 10:
        raise ValueError("d-aug oracle Q1 capacity gate is fixed at 10 spatial / 20 spin orbitals")
    try:
        from pyscf import ao2mo, fci, gto, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("d-aug oracle capacity gate requires the OES q1 extra (PySCF)") from exc

    mol = gto.M(atom="He 0 0 0", basis=source_basis, unit="Bohr", charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("helium RHF did not converge")

    C = np.asarray(mf.mo_coeff, dtype=float)
    norb = C.shape[1]
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.kernel(mol, C, compact=False).reshape((norb,) * 4)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([C.T @ dip_ao[k] @ C for k in range(3)])

    sing = fci.direct_spin0.FCI()
    sing.conv_tol = 1e-10
    sing.nroots = singlet_roots
    energies, cis = sing.kernel(h1, eri, norb, 2, ecore=float(mol.energy_nuc()))
    energies = np.atleast_1d(np.asarray(energies, dtype=float))
    if not isinstance(cis, (list, tuple)):
        cis = [cis]
    e0 = float(energies[0])
    ci0 = cis[0]
    dm_ground = np.asarray(sing.make_rdm1(ci0, norb, 2), dtype=float)

    rows = []
    for i in range(1, len(cis)):
        tdm = np.asarray(sing.trans_rdm1(cis[i], ci0, norb, 2), dtype=float)
        mu = np.array([np.einsum("pq,qp->", dip[k], tdm) for k in range(3)], dtype=float)
        de = float(energies[i] - e0)
        fosc = (2.0 / 3.0) * de * float(np.dot(mu, mu))
        rows.append((i, de, fosc))

    dark_rows = [row for row in rows if row[2] < 1e-6]
    bright_rows = [row for row in rows if row[2] > bright_threshold]
    if not dark_rows or len(bright_rows) < 3:
        raise RuntimeError("d-aug oracle could not identify dark and complete bright source classes")
    dark_idx, dark_de, _ = min(dark_rows, key=lambda row: row[1])
    first_bright_de = min(row[1] for row in bright_rows)
    bright_manifold = [row for row in bright_rows if abs(row[1] - first_bright_de) < 1e-8]
    if len(bright_manifold) != 3:
        raise RuntimeError(f"d-aug oracle expected three source bright components, got {len(bright_manifold)}")

    dm_dark = np.asarray(sing.make_rdm1(cis[dark_idx], norb, 2), dtype=float)
    dm_bright = np.zeros_like(dm_ground)
    bright_f_sum = 0.0
    for idx, _de, fosc in bright_manifold:
        dm_bright += np.asarray(sing.make_rdm1(cis[idx], norb, 2), dtype=float) / 3.0
        bright_f_sum += float(fosc)

    trip = fci.direct_spin1.FCI()
    trip.conv_tol = 1e-10
    etrip, citrip = trip.kernel(h1, eri, norb, (2, 0), ecore=float(mol.energy_nuc()))
    dm_trip = np.asarray(trip.make_rdm1(citrip, norb, (2, 0)), dtype=float)

    # Equal weight per physical information class. The bright class itself is
    # already averaged over all three Cartesian components, preserving rotation.
    dm_avg = 0.25 * (dm_ground + dm_dark + dm_trip + dm_bright)
    dm_avg = 0.5 * (dm_avg + dm_avg.T)
    occ, U = np.linalg.eigh(dm_avg)
    order = np.argsort(occ)[::-1]
    occ = occ[order]
    U = U[:, order]

    groups = _occupation_groups(occ, U)
    selected_groups = _best_complete_group_subset(groups, target_spatial)
    selected_vectors = []
    selected_sizes = []
    selected_weight = 0.0
    for gi in selected_groups:
        start, end, value, vecs = groups[gi]
        selected_vectors.append(vecs)
        size = end - start
        selected_sizes.append(size)
        selected_weight += value * size
    U_active = np.column_stack(selected_vectors)
    if U_active.shape != (norb, target_spatial):
        raise RuntimeError(f"oracle group selection produced shape {U_active.shape}")
    if not np.allclose(U_active.T @ U_active, np.eye(target_spatial), atol=1e-10):
        raise RuntimeError("oracle selected natural orbitals lost orthonormality")

    receipt = DAugOracleReceipt(
        protocol="D-AUG-ORACLE-SA-FCI-SYMMETRY-GROUP-CAPACITY-20Q",
        source_basis=str(source_label),
        source_spatial_orbitals=norb,
        target_spatial_orbitals=target_spatial,
        target_spin_orbitals=2 * target_spatial,
        compression_backend="NONPREDICTIVE_FULL_SOURCE_FCI_CLASS_AVERAGE_WITH_COMPLETE_OCCUPATION_GROUPS",
        source_ground_hartree=e0,
        source_triplet_excitation_ev=(float(etrip) - e0) * HARTREE_TO_EV,
        source_dark_singlet_excitation_ev=dark_de * HARTREE_TO_EV,
        source_bright_manifold_excitation_ev=first_bright_de * HARTREE_TO_EV,
        source_bright_manifold_oscillator_strength_sum=bright_f_sum,
        occupation_group_sizes=tuple(end - start for start, end, _value, _vecs in groups[:16]),
        occupation_group_values=tuple(value for _start, _end, value, _vecs in groups[:16]),
        selected_group_indices=tuple(int(x) for x in selected_groups),
        selected_group_sizes=tuple(int(x) for x in selected_sizes),
        selected_occupation_weight=float(selected_weight),
    )
    return mol, mf, C @ U_active, receipt
