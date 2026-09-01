"""Full-basis helium spectral diagnostics before active-space compression."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Union

import numpy as np


HARTREE_TO_EV = 27.211_386_245_981


@dataclass(frozen=True)
class SourceSpectrumResult:
    label: str
    n_spatial_orbitals: int
    ground_hartree: float
    dark_singlet_excitation_ev: float
    triplet_excitation_ev: float
    bright_singlet_excitation_ev: float
    bright_oscillator_strength: float
    bright_transition_dipole_norm_au: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def helium_source_spectrum(basis: Union[str, dict], label: str, singlet_roots: int = 6) -> SourceSpectrumResult:
    try:
        from pyscf import ao2mo, fci, gto, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("source spectrum requires the OES q1 extra (PySCF)") from exc

    mol = gto.M(atom="He 0 0 0", basis=basis, unit="Bohr", charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"RHF failed for {label}")

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

    rows = []
    for i in range(1, len(cis)):
        tdm = np.asarray(sing.trans_rdm1(cis[i], ci0, norb, 2), dtype=float)
        mu = np.array([np.einsum("pq,qp->", dip[k], tdm) for k in range(3)], dtype=float)
        mu2 = float(np.dot(mu, mu))
        delta_h = float(energies[i] - e0)
        rows.append((i, delta_h, (2.0 / 3.0) * delta_h * mu2, float(np.sqrt(mu2))))
    dark = [row for row in rows if row[2] < 1e-6]
    bright = [row for row in rows if row[2] > 1e-4]
    if not dark or not bright:
        raise RuntimeError(f"could not identify dark/bright singlet source roots for {label}")
    dark_row = min(dark, key=lambda row: row[1])
    bright_row = min(bright, key=lambda row: row[1])

    trip = fci.direct_spin1.FCI()
    trip.conv_tol = 1e-10
    etrip, _ = trip.kernel(h1, eri, norb, (2, 0), ecore=float(mol.energy_nuc()))

    return SourceSpectrumResult(
        label=label,
        n_spatial_orbitals=norb,
        ground_hartree=e0,
        dark_singlet_excitation_ev=dark_row[1] * HARTREE_TO_EV,
        triplet_excitation_ev=(float(etrip) - e0) * HARTREE_TO_EV,
        bright_singlet_excitation_ev=bright_row[1] * HARTREE_TO_EV,
        bright_oscillator_strength=float(bright_row[2]),
        bright_transition_dipole_norm_au=float(bright_row[3]),
    )
