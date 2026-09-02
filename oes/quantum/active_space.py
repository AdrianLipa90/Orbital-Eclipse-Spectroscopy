"""Transition-aware active-space compression for the fixed OES-Q1 20Q budget.

TAAS-v1 is a reference-compression protocol.  A larger classical full-basis FCI
is used only to obtain the ground-state 1-RDM.  The 10-orbital target subspace
is then built without using any experimental excited-state energy:

1. dominant ground-state natural orbital,
2. the three Cartesian dipole-response directions generated from it,
3. lowest-energy canonical virtual directions, projected against the existing
   subspace until the 10-orbital budget is filled.

This is deliberately labelled REFERENCE_COMPRESSION: full-basis FCI is not a
scalable ingredient for future large-system quantum advantage.  The same
selection interface can later accept a quantum-measured or approximate 1-RDM.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class ActiveSpaceReceipt:
    protocol: str
    source_basis: str
    source_spatial_orbitals: int
    target_spatial_orbitals: int
    target_spin_orbitals: int
    compression_backend: str
    ground_fci_hartree: float
    dominant_natural_occupation: float
    selected_labels: Tuple[str, ...]

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _append_orthonormal(
    vectors: List[np.ndarray],
    labels: List[str],
    candidate: np.ndarray,
    label: str,
    tol: float = 1e-10,
) -> bool:
    """Modified Gram-Schmidt append in an orthonormal MO coefficient basis."""
    v = np.asarray(candidate, dtype=float).copy()
    for q in vectors:
        v -= q * float(np.dot(q, v))
    # second pass for stability near degenerate atomic subspaces
    for q in vectors:
        v -= q * float(np.dot(q, v))
    norm = float(np.linalg.norm(v))
    if norm <= tol:
        return False
    vectors.append(v / norm)
    labels.append(label)
    return True


def build_helium_taas_v1(source_basis: str = "aug-cc-pVQZ", target_spatial: int = 10):
    """Build the frozen TAAS-v1 He active space and return AO coefficients.

    Returns `(mol, mf, C_active, receipt)`.  `C_active` is AO x target_spatial.
    """
    if target_spatial != 10:
        raise ValueError("TAAS-v1 canonical Q1 contract is fixed at 10 spatial orbitals / 20 qubits")
    try:
        from pyscf import ao2mo, fci, gto, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("TAAS-v1 requires the OES q1 extra (PySCF)") from exc

    mol = gto.M(
        atom="He 0 0 0",
        basis=source_basis,
        unit="Bohr",
        charge=0,
        spin=0,
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("helium RHF did not converge")

    C = np.asarray(mf.mo_coeff, dtype=float)
    n_full = C.shape[1]
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.kernel(mol, C, compact=False).reshape((n_full,) * 4)

    solver = fci.direct_spin0.FCI()
    solver.conv_tol = 1e-11
    e_fci, ci = solver.kernel(h1, eri, n_full, 2, ecore=float(mol.energy_nuc()))
    dm1 = np.asarray(solver.make_rdm1(ci, n_full, 2), dtype=float)
    dm1 = 0.5 * (dm1 + dm1.T)
    occupations, natural_vectors = np.linalg.eigh(dm1)
    order = np.argsort(occupations)[::-1]
    occupations = occupations[order]
    natural_vectors = natural_vectors[:, order]

    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip_mo = np.stack([C.T @ dip_ao[k] @ C for k in range(3)])

    vectors: List[np.ndarray] = []
    labels: List[str] = []
    dominant = natural_vectors[:, 0]
    _append_orthonormal(vectors, labels, dominant, "ground-NO[0]")

    # Operator-reachable flavor subspace.  This is selected from the Hamiltonian
    # representation itself, with no reference to a measured 1P transition energy.
    for axis, name in enumerate(("x", "y", "z")):
        _append_orthonormal(vectors, labels, dip_mo[axis] @ dominant, f"dipole-{name}|ground-NO[0]")

    # Fill remaining budget from low-energy canonical directions.  They carry the
    # diffuse radial response needed for low-lying s/p states while projection
    # prevents duplicate allocation to the dipole subspace.
    canonical_order = np.argsort(np.asarray(mf.mo_energy, dtype=float))
    for idx in canonical_order:
        if len(vectors) >= target_spatial:
            break
        unit = np.zeros(n_full, dtype=float)
        unit[int(idx)] = 1.0
        _append_orthonormal(vectors, labels, unit, f"canonical-MO[{int(idx)}]")

    # If degeneracy/projection prevented filling, fall back to natural orbitals by
    # decreasing ground-state occupation.  This is deterministic and intrinsic.
    for k in range(n_full):
        if len(vectors) >= target_spatial:
            break
        _append_orthonormal(vectors, labels, natural_vectors[:, k], f"ground-NO[{k}]")

    if len(vectors) != target_spatial:
        raise RuntimeError(f"TAAS-v1 constructed only {len(vectors)} independent orbitals")

    U = np.column_stack(vectors)
    gram = U.T @ U
    if not np.allclose(gram, np.eye(target_spatial), atol=1e-10):
        raise RuntimeError("TAAS-v1 MO-space vectors lost orthonormality")
    C_active = C @ U
    receipt = ActiveSpaceReceipt(
        protocol="TAAS-v1",
        source_basis=source_basis,
        source_spatial_orbitals=n_full,
        target_spatial_orbitals=target_spatial,
        target_spin_orbitals=2 * target_spatial,
        compression_backend="REFERENCE_COMPRESSION_FULL_BASIS_FCI_RDM1",
        ground_fci_hartree=float(e_fci),
        dominant_natural_occupation=float(occupations[0]),
        selected_labels=tuple(labels),
    )
    return mol, mf, C_active, receipt


def transform_active_integrals(mol, mf, C_active: np.ndarray):
    """Transform h1, ERI, and Cartesian dipole integrals into an active AO basis."""
    from pyscf import ao2mo

    C_active = np.asarray(C_active, dtype=float)
    n = C_active.shape[1]
    h1 = C_active.T @ mf.get_hcore() @ C_active
    eri = ao2mo.kernel(mol, C_active, compact=False).reshape((n,) * 4)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([C_active.T @ dip_ao[k] @ C_active for k in range(3)])
    return h1, eri, dip
