"""ROAS-v1: relational operator-response active-space selection at fixed 20Q.

No excited-state target energy or target wavefunction enters selection.  A full
source-basis ground-state FCI reference is probed by symmetry-classified one-
body operators:

    scalar classes: 1, r^2, projected r^4
    vector classes: r, symmetrized r^2 r

For each nontrivial class we construct a normalized response state and its
one-body density. Cartesian vector components are averaged as complete classes.
The class-averaged information density is diagonalized and exactly 10 spatial
orbitals are chosen only through complete numerically degenerate occupation
blocks.  This is a REFERENCE_COMPRESSION selector because ground FCI is used;
it does not use measured or full-FCI excited-state data.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Tuple

import numpy as np


@dataclass(frozen=True)
class ROASReceipt:
    protocol: str
    source_basis: str
    source_spatial_orbitals: int
    target_spatial_orbitals: int
    target_spin_orbitals: int
    compression_backend: str
    ground_fci_hartree: float
    class_weights: Dict[str, float]
    response_norms: Dict[str, float]
    occupation_group_sizes: Tuple[int, ...]
    occupation_group_values: Tuple[float, ...]
    selected_group_indices: Tuple[int, ...]
    selected_group_sizes: Tuple[int, ...]
    selected_occupation_weight: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _center_operator(op: np.ndarray, dm1: np.ndarray, nelec: int) -> np.ndarray:
    expectation = float(np.einsum("pq,qp->", op, dm1).real)
    return np.asarray(op, dtype=float) - (expectation / float(nelec)) * np.eye(op.shape[0])


def _response_dm(solver, op: np.ndarray, ci0, norb: int, nelec: int):
    response = np.asarray(solver.contract_1e(np.asarray(op, dtype=float), ci0, norb, nelec))
    norm = float(np.linalg.norm(response))
    if norm <= 1e-12:
        raise RuntimeError("ROAS operator response vanished")
    response = response / norm
    dm = np.asarray(solver.make_rdm1(response, norb, nelec), dtype=float)
    return 0.5 * (dm + dm.T), norm


def _groups(values: np.ndarray, vectors: np.ndarray, atol: float = 1e-8):
    out = []
    start = 0
    while start < len(values):
        ref = float(values[start])
        end = start + 1
        while end < len(values) and abs(float(values[end]) - ref) <= max(atol, atol * abs(ref)):
            end += 1
        out.append((start, end, float(np.mean(values[start:end])), vectors[:, start:end]))
        start = end
    return out


def _best_full_groups(groups, target: int) -> Tuple[int, ...]:
    dp = {0: (0.0, tuple())}
    for gi, (start, end, value, _vecs) in enumerate(groups):
        size = end - start
        weight = value * size
        nxt = dict(dp)
        for count, (score, selected) in dp.items():
            nc = count + size
            if nc > target:
                continue
            cand = (score + weight, selected + (gi,))
            if nc not in nxt or cand[0] > nxt[nc][0]:
                nxt[nc] = cand
        dp = nxt
    if target not in dp:
        raise RuntimeError(
            f"ROAS cannot fill {target} orbitals with complete occupation blocks; "
            f"leading sizes={[end-start for start,end,_v,_u in groups[:16]]}"
        )
    return dp[target][1]


def build_helium_roas_v1_20q(source_basis: Any, source_label: str, target_spatial: int = 10):
    if target_spatial != 10:
        raise ValueError("ROAS-v1 Q1 contract is fixed at 10 spatial / 20 spin orbitals")
    try:
        from pyscf import ao2mo, fci, gto, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("ROAS-v1 requires the OES q1 extra (PySCF)") from exc

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

    solver = fci.direct_spin0.FCI()
    solver.conv_tol = 1e-11
    e0, ci0 = solver.kernel(h1, eri, norb, 2, ecore=float(mol.energy_nuc()))
    dm0 = np.asarray(solver.make_rdm1(ci0, norb, 2), dtype=float)
    dm0 = 0.5 * (dm0 + dm0.T)

    r2_ao = mol.intor("int1e_r2", hermi=1)
    r2 = C.T @ r2_ao @ C
    # Finite-basis projection of r^4. Symmetrization suppresses roundoff.
    r4 = 0.5 * (r2 @ r2 + (r2 @ r2).T)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([C.T @ dip_ao[k] @ C for k in range(3)])
    rdip = np.stack([0.5 * (r2 @ dip[k] + dip[k] @ r2) for k in range(3)])

    response_norms: Dict[str, float] = {}
    dm_r2, response_norms["r2"] = _response_dm(solver, _center_operator(r2, dm0, 2), ci0, norb, 2)
    dm_r4, response_norms["r4"] = _response_dm(solver, _center_operator(r4, dm0, 2), ci0, norb, 2)

    dm_dip = np.zeros_like(dm0)
    dm_rdip = np.zeros_like(dm0)
    for k, axis in enumerate(("x", "y", "z")):
        dm_axis, response_norms[f"r:{axis}"] = _response_dm(
            solver, _center_operator(dip[k], dm0, 2), ci0, norb, 2
        )
        dm_dip += dm_axis / 3.0
        dm_raxis, response_norms[f"r2r:{axis}"] = _response_dm(
            solver, _center_operator(rdip[k], dm0, 2), ci0, norb, 2
        )
        dm_rdip += dm_raxis / 3.0

    # Equal information-class weights, not tuned to NIST or excited-state data.
    weights = {"ground": 0.2, "r2": 0.2, "r4": 0.2, "r": 0.2, "r2r": 0.2}
    dm_info = (
        weights["ground"] * dm0
        + weights["r2"] * dm_r2
        + weights["r4"] * dm_r4
        + weights["r"] * dm_dip
        + weights["r2r"] * dm_rdip
    )
    dm_info = 0.5 * (dm_info + dm_info.T)

    occ, U = np.linalg.eigh(dm_info)
    order = np.argsort(occ)[::-1]
    occ = occ[order]
    U = U[:, order]
    groups = _groups(occ, U)
    selected = _best_full_groups(groups, target_spatial)

    blocks: List[np.ndarray] = []
    sizes: List[int] = []
    retained = 0.0
    for gi in selected:
        start, end, value, vecs = groups[gi]
        blocks.append(vecs)
        size = end - start
        sizes.append(size)
        retained += value * size
    U_active = np.column_stack(blocks)
    if U_active.shape != (norb, target_spatial):
        raise RuntimeError(f"ROAS selected shape {U_active.shape}, expected {(norb, target_spatial)}")
    if not np.allclose(U_active.T @ U_active, np.eye(target_spatial), atol=1e-10):
        raise RuntimeError("ROAS selected orbitals lost orthonormality")

    receipt = ROASReceipt(
        protocol="D-AUG-ROAS-V1-OPERATOR-CLASS-RESPONSE-20Q",
        source_basis=str(source_label),
        source_spatial_orbitals=norb,
        target_spatial_orbitals=target_spatial,
        target_spin_orbitals=2 * target_spatial,
        compression_backend="REFERENCE_GROUND_FCI_OPERATOR_CLASS_RESPONSE_WITH_COMPLETE_OCCUPATION_GROUPS",
        ground_fci_hartree=float(e0),
        class_weights={k: float(v) for k, v in weights.items()},
        response_norms={k: float(v) for k, v in response_norms.items()},
        occupation_group_sizes=tuple(end - start for start, end, _value, _vecs in groups[:16]),
        occupation_group_values=tuple(value for _start, _end, value, _vecs in groups[:16]),
        selected_group_indices=tuple(int(x) for x in selected),
        selected_group_sizes=tuple(int(x) for x in sizes),
        selected_occupation_weight=float(retained),
    )
    return mol, mf, C @ U_active, receipt
