"""ROAS-v2 minimal operator-family selector at fixed 20Q.

Ablation of ROAS-v1: the higher radial-vector probe r^2 r is removed.  The
selector uses four equally weighted information classes only:

    ground density, scalar r^2 response, scalar projected-r^4 response,
    complete Cartesian dipole r response.

The second p-like correlation block, if selected, must therefore emerge from
ground/scalar correlation structure rather than being explicitly driven by a
second vector operator. No excited-state target or measured spectrum enters.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

from .active_space_roas import ROASReceipt, _best_full_groups, _center_operator, _groups, _response_dm


def build_helium_roas_v2_20q(source_basis: Any, source_label: str, target_spatial: int = 10):
    if target_spatial != 10:
        raise ValueError("ROAS-v2 Q1 contract is fixed at 10 spatial / 20 spin orbitals")
    try:
        from pyscf import ao2mo, fci, gto, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("ROAS-v2 requires the OES q1 extra (PySCF)") from exc

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
    r4_raw = r2 @ r2
    r4 = 0.5 * (r4_raw + r4_raw.T)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([C.T @ dip_ao[k] @ C for k in range(3)])

    response_norms: Dict[str, float] = {}
    dm_r2, response_norms["r2"] = _response_dm(solver, _center_operator(r2, dm0, 2), ci0, norb, 2)
    dm_r4, response_norms["r4"] = _response_dm(solver, _center_operator(r4, dm0, 2), ci0, norb, 2)
    dm_dip = np.zeros_like(dm0)
    for k, axis in enumerate(("x", "y", "z")):
        dm_axis, response_norms[f"r:{axis}"] = _response_dm(
            solver, _center_operator(dip[k], dm0, 2), ci0, norb, 2
        )
        dm_dip += dm_axis / 3.0

    weights = {"ground": 0.25, "r2": 0.25, "r4": 0.25, "r": 0.25}
    dm_info = (
        weights["ground"] * dm0
        + weights["r2"] * dm_r2
        + weights["r4"] * dm_r4
        + weights["r"] * dm_dip
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
        raise RuntimeError(f"ROAS-v2 selected shape {U_active.shape}")
    if not np.allclose(U_active.T @ U_active, np.eye(target_spatial), atol=1e-10):
        raise RuntimeError("ROAS-v2 selected orbitals lost orthonormality")

    receipt = ROASReceipt(
        protocol="D-AUG-ROAS-V2-MINIMAL-OPERATOR-FAMILY-20Q",
        source_basis=str(source_label),
        source_spatial_orbitals=norb,
        target_spatial_orbitals=target_spatial,
        target_spin_orbitals=2 * target_spatial,
        compression_backend="REFERENCE_GROUND_FCI_MINIMAL_OPERATOR_FAMILY_RESPONSE_WITH_COMPLETE_OCCUPATION_GROUPS",
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
