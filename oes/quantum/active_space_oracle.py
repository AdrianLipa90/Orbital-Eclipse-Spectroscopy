"""Oracle state-averaged active space: an upper bound on 20Q capacity.

This module is deliberately NOT a predictive compression method.  It uses
full-basis FCI wavefunctions for the ground state and three target state classes
(dark excited singlet, lowest triplet, first dipole-bright singlet) to form a
state-averaged 1-RDM.  The top 10 natural orbitals are then used as a 20Q active
space.  The result answers a capacity question: can 20 spin-orbital modes retain
these spectral states if the compression is informed by near-ideal source
wavefunctions?
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Tuple

import numpy as np


@dataclass(frozen=True)
class OracleCompressionReceipt:
    protocol: str
    source_basis: str
    source_spatial_orbitals: int
    target_spatial_orbitals: int
    target_spin_orbitals: int
    compression_backend: str
    source_ground_hartree: float
    source_dark_singlet_excitation_ev: float
    source_triplet_excitation_ev: float
    source_bright_singlet_excitation_ev: float
    source_bright_oscillator_strength: float
    averaged_natural_occupations: Tuple[float, ...]

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


HARTREE_TO_EV = 27.211_386_245_981


def build_helium_oracle_capacity_space(
    source_basis: str = "aug-cc-pVQZ",
    target_spatial: int = 10,
    singlet_roots: int = 10,
):
    if target_spatial != 10:
        raise ValueError("oracle Q1 capacity gate is fixed at 10 spatial / 20 spin orbitals")
    try:
        from pyscf import ao2mo, fci, gto, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("oracle Q1 capacity gate requires the q1 extra (PySCF)") from exc

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
    if len(cis) < 3:
        raise RuntimeError("full-basis singlet FCI returned too few roots")
    e0 = float(energies[0])
    ci0 = cis[0]
    dm_ground = np.asarray(sing.make_rdm1(ci0, norb, 2), dtype=float)

    singlet_rows = []
    for i in range(1, len(cis)):
        tdm = np.asarray(sing.trans_rdm1(cis[i], ci0, norb, 2), dtype=float)
        mu = np.array([np.einsum("pq,qp->", dip[k], tdm) for k in range(3)], dtype=float)
        delta_h = float(energies[i] - e0)
        fosc = (2.0 / 3.0) * delta_h * float(np.dot(mu, mu))
        singlet_rows.append((i, delta_h, fosc))

    dark_rows = [row for row in singlet_rows if row[2] < 1e-6]
    bright_rows = [row for row in singlet_rows if row[2] > 1e-4]
    if not dark_rows or not bright_rows:
        raise RuntimeError("could not identify dark and bright source singlet roots")
    dark_idx, dark_delta_h, _ = min(dark_rows, key=lambda row: row[1])
    bright_idx, bright_delta_h, bright_f = min(bright_rows, key=lambda row: row[1])
    dm_dark = np.asarray(sing.make_rdm1(cis[dark_idx], norb, 2), dtype=float)
    dm_bright = np.asarray(sing.make_rdm1(cis[bright_idx], norb, 2), dtype=float)

    trip = fci.direct_spin1.FCI()
    trip.conv_tol = 1e-10
    etrip, citrip = trip.kernel(h1, eri, norb, (2, 0), ecore=float(mol.energy_nuc()))
    dm_trip = np.asarray(trip.make_rdm1(citrip, norb, (2, 0)), dtype=float)

    dm_avg = 0.25 * (dm_ground + dm_dark + dm_trip + dm_bright)
    dm_avg = 0.5 * (dm_avg + dm_avg.T)
    occ, U = np.linalg.eigh(dm_avg)
    order = np.argsort(occ)[::-1]
    occ = occ[order]
    U = U[:, order]
    U_active = U[:, :target_spatial]
    C_active = C @ U_active

    receipt = OracleCompressionReceipt(
        protocol="ORACLE-SA-FCI-NO-CAPACITY-BOUND",
        source_basis=source_basis,
        source_spatial_orbitals=norb,
        target_spatial_orbitals=target_spatial,
        target_spin_orbitals=2 * target_spatial,
        compression_backend="NONPREDICTIVE_FULL_BASIS_FCI_TARGET_STATE_AVERAGE",
        source_ground_hartree=e0,
        source_dark_singlet_excitation_ev=dark_delta_h * HARTREE_TO_EV,
        source_triplet_excitation_ev=(float(etrip) - e0) * HARTREE_TO_EV,
        source_bright_singlet_excitation_ev=bright_delta_h * HARTREE_TO_EV,
        source_bright_oscillator_strength=float(bright_f),
        averaged_natural_occupations=tuple(float(x) for x in occ[:target_spatial]),
    )
    return mol, mf, C_active, receipt
