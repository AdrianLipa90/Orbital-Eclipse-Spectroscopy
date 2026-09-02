"""TAAS-v2: state-balanced operator-response natural orbitals at fixed 20Q.

The source is a full aug-cc-pVQZ helium reference.  No measured excited-state
energy is used to choose the target 10 spatial orbitals.  We state-average
one-body density matrices from three intrinsic information classes:

- 1/3 ground-state correlation density,
- 1/3 centered scalar radial response r^2|Psi0>,
- 1/3 vector dipole response, split equally over x/y/z.

The top 10 natural orbitals of that state-averaged density define the 20-spin-
orbital target.  This remains a REFERENCE_COMPRESSION because the source
response vectors use full-basis classical FCI; later versions can replace the
source RDMs with quantum-measured or approximate RDMs.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Tuple

import numpy as np


@dataclass(frozen=True)
class TAASv2Receipt:
    protocol: str
    source_basis: str
    source_spatial_orbitals: int
    target_spatial_orbitals: int
    target_spin_orbitals: int
    compression_backend: str
    source_ground_fci_hartree: float
    weights: Dict[str, float]
    averaged_natural_occupations: Tuple[float, ...]
    response_norms: Dict[str, float]

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _center_one_body_operator(op: np.ndarray, dm1: np.ndarray, n_electrons: int) -> np.ndarray:
    """Return O-cI so <Psi|sum O_pq E_pq|Psi>=0 for the reference state."""
    expectation = float(np.einsum("pq,qp->", op, dm1).real)
    return np.asarray(op, dtype=float) - (expectation / float(n_electrons)) * np.eye(op.shape[0])


def _normalized_response_density(solver, op: np.ndarray, ci, norb: int, nelec: int):
    response = np.asarray(solver.contract_1e(np.asarray(op, dtype=float), ci, norb, nelec))
    norm = float(np.linalg.norm(response))
    if norm <= 1e-12:
        raise RuntimeError("operator response vanished in full-basis FCI reference")
    response = response / norm
    dm = np.asarray(solver.make_rdm1(response, norb, nelec), dtype=float)
    dm = 0.5 * (dm + dm.T)
    return dm, norm


def build_helium_taas_v2(source_basis: str = "aug-cc-pVQZ", target_spatial: int = 10):
    if target_spatial != 10:
        raise ValueError("TAAS-v2 canonical Q1 target is fixed at 10 spatial / 20 spin orbitals")
    try:
        from pyscf import ao2mo, fci, gto, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("TAAS-v2 requires the OES q1 extra (PySCF)") from exc

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
    e0, ci0 = solver.kernel(h1, eri, n_full, 2, ecore=float(mol.energy_nuc()))
    dm0 = np.asarray(solver.make_rdm1(ci0, n_full, 2), dtype=float)
    dm0 = 0.5 * (dm0 + dm0.T)

    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([C.T @ dip_ao[k] @ C for k in range(3)])
    r2_ao = mol.intor("int1e_r2", hermi=1)
    r2 = C.T @ r2_ao @ C

    response_dms = {}
    response_norms = {}
    for axis, name in enumerate(("x", "y", "z")):
        centered = _center_one_body_operator(dip[axis], dm0, 2)
        response_dms[name], response_norms[name] = _normalized_response_density(solver, centered, ci0, n_full, 2)
    centered_r2 = _center_one_body_operator(r2, dm0, 2)
    response_dms["r2"], response_norms["r2"] = _normalized_response_density(solver, centered_r2, ci0, n_full, 2)

    weights = {
        "ground": 1.0 / 3.0,
        "r2": 1.0 / 3.0,
        "x": 1.0 / 9.0,
        "y": 1.0 / 9.0,
        "z": 1.0 / 9.0,
    }
    dm_avg = weights["ground"] * dm0 + weights["r2"] * response_dms["r2"]
    for axis in ("x", "y", "z"):
        dm_avg += weights[axis] * response_dms[axis]
    dm_avg = 0.5 * (dm_avg + dm_avg.T)

    occ, U = np.linalg.eigh(dm_avg)
    order = np.argsort(occ)[::-1]
    occ = occ[order]
    U = U[:, order]
    U_active = U[:, :target_spatial]
    if not np.allclose(U_active.T @ U_active, np.eye(target_spatial), atol=1e-10):
        raise RuntimeError("TAAS-v2 natural orbitals lost orthonormality")
    C_active = C @ U_active

    receipt = TAASv2Receipt(
        protocol="TAAS-v2-operator-response-NO",
        source_basis=source_basis,
        source_spatial_orbitals=n_full,
        target_spatial_orbitals=target_spatial,
        target_spin_orbitals=2 * target_spatial,
        compression_backend="REFERENCE_COMPRESSION_FULL_BASIS_FCI_OPERATOR_RESPONSE_RDM1",
        source_ground_fci_hartree=float(e0),
        weights={k: float(v) for k, v in weights.items()},
        averaged_natural_occupations=tuple(float(x) for x in occ[:target_spatial]),
        response_norms={k: float(v) for k, v in response_norms.items()},
    )
    return mol, mf, C_active, receipt
