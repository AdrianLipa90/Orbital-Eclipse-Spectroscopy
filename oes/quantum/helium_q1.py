"""OES-Q1: helium in a 10-spatial / 20-spin-orbital active space.

PySCF is used only to generate standard Gaussian-basis one- and two-electron
integrals and to provide an independent FCI reference.  OES constructs and
diagonalizes its own full fixed-particle spin-orbital Hamiltonian.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import sqrt
from typing import Dict, List, Tuple

import numpy as np

from .fermions import (
    annihilate,
    build_sector_hamiltonian,
    create,
    determinant_basis,
    full_space_dimension,
    occupied,
    sector_dimension,
    transition_one_rdm,
)


@dataclass(frozen=True)
class Q1State:
    index: int
    energy_hartree: float
    excitation_ev: float
    s2: float
    spin_s: float


@dataclass(frozen=True)
class HeliumQ1Result:
    backend: str
    basis_name: str
    n_spatial_orbitals: int
    n_spin_orbitals: int
    n_electrons: int
    full_qubit_dimension: int
    fixed_particle_dimension: int
    rhf_energy_hartree: float
    oes_fci_energy_hartree: float
    pyscf_fci_energy_hartree: float
    fci_delta_hartree: float
    ground_s2: float
    first_triplet_index: int
    first_triplet_energy_hartree: float
    first_singlet_excited_index: int
    first_singlet_excited_energy_hartree: float
    singlet_triplet_spatial_transition_norm: float
    singlet_singlet_spatial_transition_norm: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


HARTREE_TO_EV = 27.211_386_245_981


def _one_body_spin_operator(coeff: np.ndarray, basis: List[int]) -> np.ndarray:
    """Matrix of sum_pq coeff[p,q] a_p^† a_q in a fixed-N basis."""
    coeff = np.asarray(coeff, dtype=complex)
    n_spin = coeff.shape[0]
    index = {det: i for i, det in enumerate(basis)}
    out = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, det in enumerate(basis):
        for q in occupied(det, n_spin):
            a = annihilate(det, q)
            assert a is not None
            d1, s1 = a
            for p in range(n_spin):
                value = coeff[p, q]
                if abs(value) < 1e-15:
                    continue
                c = create(d1, p)
                if c is None:
                    continue
                d2, s2 = c
                row = index.get(d2)
                if row is not None:
                    out[row, col] += value * s1 * s2
    return out


def spin_squared_matrix(n_spatial: int, basis: List[int]) -> np.ndarray:
    """Exact S^2/hbar^2 matrix for alpha/beta interleaved spin orbitals."""
    n_spin = 2 * n_spatial
    sz_coeff = np.zeros((n_spin, n_spin), dtype=complex)
    sp_coeff = np.zeros_like(sz_coeff)
    sm_coeff = np.zeros_like(sz_coeff)
    for p in range(n_spatial):
        a = 2 * p
        b = a + 1
        sz_coeff[a, a] = 0.5
        sz_coeff[b, b] = -0.5
        sp_coeff[a, b] = 1.0
        sm_coeff[b, a] = 1.0
    sz = _one_body_spin_operator(sz_coeff, basis)
    sp = _one_body_spin_operator(sp_coeff, basis)
    sm = _one_body_spin_operator(sm_coeff, basis)
    return sz @ sz + 0.5 * (sp @ sm + sm @ sp)


def _spin_from_s2(s2: float) -> float:
    # Solve S(S+1)=s2 and snap tiny numerical noise.
    s = 0.5 * (-1.0 + sqrt(max(0.0, 1.0 + 4.0 * s2)))
    nearest_half = round(2.0 * s) / 2.0
    return nearest_half if abs(s - nearest_half) < 1e-6 else s


def classify_states(evals: np.ndarray, evecs: np.ndarray, s2_matrix: np.ndarray, limit: int = 20) -> List[Q1State]:
    e0 = float(evals[0])
    states: List[Q1State] = []
    for k in range(min(limit, len(evals))):
        v = evecs[:, k]
        s2 = float(np.real(np.vdot(v, s2_matrix @ v)))
        states.append(
            Q1State(
                index=k,
                energy_hartree=float(evals[k]),
                excitation_ev=(float(evals[k]) - e0) * HARTREE_TO_EV,
                s2=s2,
                spin_s=_spin_from_s2(s2),
            )
        )
    return states


def spatial_transition_rdm(t_spin: np.ndarray) -> np.ndarray:
    """Spin-summed spatial transition 1-RDM from an interleaved spin 1-RDM."""
    n_spatial = t_spin.shape[0] // 2
    out = np.zeros((n_spatial, n_spatial), dtype=complex)
    for p in range(n_spatial):
        for q in range(n_spatial):
            out[p, q] = t_spin[2 * p, 2 * q] + t_spin[2 * p + 1, 2 * q + 1]
    return out


def prepare_helium_integrals(basis_name: str = "cc-pVTZ", n_spatial: int = 10):
    """Generate a canonical-MO active space for isolated helium."""
    try:
        from pyscf import ao2mo, gto, scf
    except ImportError as exc:  # pragma: no cover - exercised in optional Q1 CI
        raise RuntimeError("OES-Q1 helium integral generation requires the 'q1' extra (PySCF)") from exc

    mol = gto.M(
        atom="He 0 0 0",
        basis=basis_name,
        unit="Bohr",
        charge=0,
        spin=0,
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    rhf_energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError("helium RHF did not converge")
    if mf.mo_coeff.shape[1] < n_spatial:
        raise ValueError(f"basis {basis_name} provides fewer than {n_spatial} MOs")

    coeff = mf.mo_coeff[:, :n_spatial]
    hcore_ao = mf.get_hcore()
    h1 = coeff.T @ hcore_ao @ coeff
    eri = ao2mo.kernel(mol, coeff, compact=False).reshape((n_spatial,) * 4)
    return mol, h1, eri, rhf_energy


def pyscf_fci_reference(h1: np.ndarray, eri: np.ndarray, n_electrons: int = 2, ecore: float = 0.0) -> float:
    try:
        from pyscf import fci
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-Q1 FCI reference requires the 'q1' extra (PySCF)") from exc

    n_spatial = h1.shape[0]
    # For helium use M_s=0: one alpha and one beta electron.
    if n_electrons != 2:
        raise ValueError("current Q1 PySCF reference is specialized to two electrons")
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-12
    energy, _ = solver.kernel(h1, eri, n_spatial, (1, 1), ecore=ecore)
    return float(energy)


def run_helium_q1(basis_name: str = "cc-pVTZ", n_spatial: int = 10) -> Tuple[HeliumQ1Result, List[Q1State]]:
    """Execute the 20-spin-orbital helium benchmark and return its receipt data."""
    if n_spatial != 10:
        raise ValueError("OES-Q1 canonical gate is fixed at 10 spatial / 20 spin orbitals")

    mol, h1, eri, rhf_energy = prepare_helium_integrals(basis_name=basis_name, n_spatial=n_spatial)
    n_spin = 2 * n_spatial
    H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    evals, evecs = np.linalg.eigh(H)

    e_oes = float(evals[0])
    e_ref = pyscf_fci_reference(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    s2mat = spin_squared_matrix(n_spatial, basis)
    states = classify_states(evals, evecs, s2mat, limit=min(40, len(evals)))

    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    if not triplets or not singlets:
        raise RuntimeError("could not resolve both singlet and triplet excited sectors")
    triplet = triplets[0]
    singlet = singlets[0]

    ground = evecs[:, 0]
    t_triplet = transition_one_rdm(evecs[:, triplet.index], ground, basis, n_spin)
    t_singlet = transition_one_rdm(evecs[:, singlet.index], ground, basis, n_spin)
    spatial_triplet = spatial_transition_rdm(t_triplet)
    spatial_singlet = spatial_transition_rdm(t_singlet)
    ground_s2 = float(np.real(np.vdot(ground, s2mat @ ground)))

    result = HeliumQ1Result(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        n_spatial_orbitals=n_spatial,
        n_spin_orbitals=n_spin,
        n_electrons=2,
        full_qubit_dimension=full_space_dimension(n_spin),
        fixed_particle_dimension=sector_dimension(n_spin, 2),
        rhf_energy_hartree=rhf_energy,
        oes_fci_energy_hartree=e_oes,
        pyscf_fci_energy_hartree=e_ref,
        fci_delta_hartree=e_oes - e_ref,
        ground_s2=ground_s2,
        first_triplet_index=triplet.index,
        first_triplet_energy_hartree=triplet.energy_hartree,
        first_singlet_excited_index=singlet.index,
        first_singlet_excited_energy_hartree=singlet.energy_hartree,
        singlet_triplet_spatial_transition_norm=float(np.linalg.norm(spatial_triplet)),
        singlet_singlet_spatial_transition_norm=float(np.linalg.norm(spatial_singlet)),
    )
    return result, states
