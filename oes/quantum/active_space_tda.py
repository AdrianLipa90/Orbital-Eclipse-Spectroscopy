"""TAAS-v2: TDA/NTO transition-aware compression at fixed 20 qubits.

The held-out experimental 1s2p 1P energy is not used to choose the active
space.  The bright subspace is obtained from an ab-initio TDA response of the
source Hamiltonian.  Ground correlation is represented by the dominant
full-basis natural orbital, then up to three low-energy bright virtual NTOs are
added, and the remaining budget is filled by low-energy canonical directions.

Like TAAS-v1, this is a REFERENCE_COMPRESSION protocol: the full-basis ground
FCI 1-RDM is used only to provide a clean reference natural orbital for this
small two-electron benchmark.  The TDA response itself is polynomial-cost and
is intended as the scalable selection primitive.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, List, Tuple

import numpy as np

from .active_space import _append_orthonormal


@dataclass(frozen=True)
class TDAActiveSpaceReceipt:
    protocol: str
    source_basis: str
    source_spatial_orbitals: int
    target_spatial_orbitals: int
    target_spin_orbitals: int
    compression_backend: str
    ground_fci_hartree: float
    dominant_natural_occupation: float
    selected_tda_states: Tuple[int, ...]
    selected_tda_excitation_eV: Tuple[float, ...]
    selected_tda_oscillator_strengths: Tuple[float, ...]
    selected_labels: Tuple[str, ...]

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def build_helium_taas_v2(
    source_basis: str = "aug-cc-pVQZ",
    target_spatial: int = 10,
    n_tda_states: int = 12,
    bright_threshold: float = 1e-5,
):
    """Build a 10-spatial-orbital He active space using TDA bright-state NTOs."""
    if target_spatial != 10:
        raise ValueError("TAAS-v2 canonical Q1 contract is fixed at 10 spatial orbitals / 20 qubits")
    try:
        from pyscf import ao2mo, fci, gto, scf, tdscf
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
    overlap = np.asarray(mf.get_ovlp(), dtype=float)
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.kernel(mol, C, compact=False).reshape((n_full,) * 4)

    # Ground-state reference natural orbital.  Kept identical in spirit to v1
    # so the v1/v2 comparison isolates the excited-response selection change.
    fci_solver = fci.direct_spin0.FCI()
    fci_solver.conv_tol = 1e-11
    e_fci, ci = fci_solver.kernel(h1, eri, n_full, 2, ecore=float(mol.energy_nuc()))
    dm1 = np.asarray(fci_solver.make_rdm1(ci, n_full, 2), dtype=float)
    dm1 = 0.5 * (dm1 + dm1.T)
    occupations, natural_vectors = np.linalg.eigh(dm1)
    order = np.argsort(occupations)[::-1]
    occupations = occupations[order]
    natural_vectors = natural_vectors[:, order]

    # Purely theoretical bright-response selector.  No experimental level is
    # consulted here.  TDA X amplitudes are used through PySCF's NTO analysis.
    td = tdscf.TDA(mf)
    td.singlet = True
    td.nstates = n_tda_states
    tda_e, _ = td.kernel()
    tda_e = np.asarray(tda_e, dtype=float)
    osc = np.asarray(td.oscillator_strength(), dtype=float)
    if tda_e.size == 0:
        raise RuntimeError("TAAS-v2 TDA returned no excited states")

    bright = [i for i, f in enumerate(osc) if float(f) > bright_threshold]
    if not bright:
        raise RuntimeError("TAAS-v2 TDA found no bright singlet state")
    # For an isolated atom the first 1P response is triply degenerate.  Taking
    # the first three bright roots captures all Cartesian directions without
    # using the experimental 1P energy.
    bright = bright[:3]

    vectors: List[np.ndarray] = []
    labels: List[str] = []
    dominant = natural_vectors[:, 0]
    _append_orthonormal(vectors, labels, dominant, "ground-NO[0]")

    nocc = int(np.count_nonzero(np.asarray(mf.mo_occ) > 0))
    if nocc < 1:
        raise RuntimeError("TAAS-v2 expected at least one occupied RHF orbital")

    selected_exc_ev: List[float] = []
    selected_osc: List[float] = []
    selected_states: List[int] = []
    HARTREE_TO_EV = 27.211_386_245_981
    for idx in bright:
        weights, nto_ao = td.get_nto(state=int(idx) + 1, threshold=0.0)
        nto_ao = np.asarray(nto_ao, dtype=float)
        if nto_ao.shape[1] <= nocc:
            raise RuntimeError("TAAS-v2 NTO matrix does not contain a virtual NTO")
        # PySCF returns AO-basis NTOs.  Convert the leading virtual NTO to the
        # orthonormal canonical-MO coordinate system used by the selector.
        nto_virtual_ao = nto_ao[:, nocc]
        nto_virtual_mo = C.T @ overlap @ nto_virtual_ao
        if _append_orthonormal(
            vectors,
            labels,
            nto_virtual_mo,
            f"TDA-NTO-virtual[state={int(idx)+1}]",
        ):
            selected_states.append(int(idx) + 1)
            selected_exc_ev.append(float(tda_e[idx] * HARTREE_TO_EV))
            selected_osc.append(float(osc[idx]))

    # Preserve low-energy radial/correlation directions for the dark 2s sector.
    canonical_order = np.argsort(np.asarray(mf.mo_energy, dtype=float))
    for idx in canonical_order:
        if len(vectors) >= target_spatial:
            break
        unit = np.zeros(n_full, dtype=float)
        unit[int(idx)] = 1.0
        _append_orthonormal(vectors, labels, unit, f"canonical-MO[{int(idx)}]")

    for k in range(n_full):
        if len(vectors) >= target_spatial:
            break
        _append_orthonormal(vectors, labels, natural_vectors[:, k], f"ground-NO[{k}]")

    if len(vectors) != target_spatial:
        raise RuntimeError(f"TAAS-v2 constructed only {len(vectors)} independent orbitals")

    U = np.column_stack(vectors)
    if not np.allclose(U.T @ U, np.eye(target_spatial), atol=1e-10):
        raise RuntimeError("TAAS-v2 MO-space vectors lost orthonormality")
    C_active = C @ U

    receipt = TDAActiveSpaceReceipt(
        protocol="TAAS-v2-TDA-NTO",
        source_basis=source_basis,
        source_spatial_orbitals=n_full,
        target_spatial_orbitals=target_spatial,
        target_spin_orbitals=2 * target_spatial,
        compression_backend="REFERENCE_GROUND_FCI_RDM1_PLUS_TDA_NTO",
        ground_fci_hartree=float(e_fci),
        dominant_natural_occupation=float(occupations[0]),
        selected_tda_states=tuple(selected_states),
        selected_tda_excitation_eV=tuple(selected_exc_ev),
        selected_tda_oscillator_strengths=tuple(selected_osc),
        selected_labels=tuple(labels),
    )
    return mol, mf, C_active, receipt
