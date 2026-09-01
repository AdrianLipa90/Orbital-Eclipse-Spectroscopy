"""OES-A1: neutral lithium in a fixed 10-spatial / 20-spin-orbital register.

The canonical A1 register remains 20 Jordan-Wigner modes. Neutral lithium has
three electrons, so the complete fixed-particle sector has C(20,3)=1140
Slater determinants. Because the nonrelativistic Hamiltonian conserves S_z,
the reference diagonalization may be restricted further to M_S=+1/2
(N_alpha=2, N_beta=1), dimension C(10,2) C(10,1)=450, without changing the
underlying 20-mode encoding.

A1 also carries forward the symmetry lesson established by Q1: a finite active
space is assembled from complete degenerate orbital blocks. The canonical
lithium selector keeps four scalar (s-like) one-dimensional energy blocks and
two complete three-dimensional (p-like) blocks, giving 4 + 3 + 3 = 10 spatial
orbitals without cutting an atomic multiplet.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .fermions import (
    build_sector_hamiltonian,
    determinant_basis,
    full_space_dimension,
    sector_dimension,
    transition_one_rdm,
)
from .helium_q1 import spatial_transition_rdm

HARTREE_TO_EV = 27.211_386_245_981


@dataclass(frozen=True)
class LithiumA1Result:
    backend: str
    basis_name: str
    active_protocol: str
    selected_mo_indices: Tuple[int, ...]
    selected_group_sizes: Tuple[int, ...]
    n_spatial_orbitals: int
    n_spin_orbitals: int
    n_electrons: int
    full_qubit_dimension: int
    fixed_particle_dimension: int
    ms_half_dimension: int
    rohf_energy_hartree: float
    oes_fci_energy_hartree: float
    pyscf_fci_energy_hartree: float
    fci_delta_hartree: float
    first_bright_excitation_ev: float
    first_bright_degeneracy: int
    first_bright_spread_ev: float
    first_bright_f_sum: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def alpha_beta_counts(det: int, n_spatial: int) -> Tuple[int, int]:
    """Count alpha/beta occupations for interleaved spin orbitals."""
    na = 0
    nb = 0
    for p in range(n_spatial):
        na += (det >> (2 * p)) & 1
        nb += (det >> (2 * p + 1)) & 1
    return int(na), int(nb)


def spin_sector_indices(
    basis: Sequence[int],
    n_spatial: int,
    n_alpha: int,
    n_beta: int,
) -> List[int]:
    return [
        i
        for i, det in enumerate(basis)
        if alpha_beta_counts(int(det), n_spatial) == (n_alpha, n_beta)
    ]


def energy_degeneracy_groups(
    energies: Sequence[float],
    absolute_tolerance: float = 1e-8,
    relative_tolerance: float = 1e-8,
) -> List[List[int]]:
    """Group adjacent canonical orbitals into numerically degenerate blocks."""
    values = np.asarray(energies, dtype=float)
    if values.ndim != 1 or len(values) == 0 or not np.all(np.isfinite(values)):
        raise ValueError("orbital energies must be a finite non-empty vector")
    groups: List[List[int]] = [[0]]
    reference = float(values[0])
    for idx in range(1, len(values)):
        value = float(values[idx])
        tol = max(absolute_tolerance, relative_tolerance * max(1.0, abs(reference), abs(value)))
        if abs(value - reference) <= tol:
            groups[-1].append(idx)
        else:
            groups.append([idx])
            reference = value
    return groups


def select_atomic_s4p6_indices(energies: Sequence[float]) -> Tuple[List[int], List[int]]:
    """Select four complete scalar blocks and two complete p-like triplets.

    The rule uses only canonical orbital-energy degeneracy, never experimental
    levels. Five-/seven-dimensional groups are skipped rather than truncated.
    """
    groups = energy_degeneracy_groups(energies)
    scalar = [group for group in groups if len(group) == 1]
    p_blocks = [group for group in groups if len(group) == 3]
    if len(scalar) < 4 or len(p_blocks) < 2:
        sizes = [len(group) for group in groups]
        raise RuntimeError(
            f"basis does not expose enough complete 1D/3D atomic blocks for s4+p6; sizes={sizes}"
        )
    chosen_groups = scalar[:4] + p_blocks[:2]
    chosen_groups.sort(key=lambda group: min(group))
    indices = sorted(idx for group in chosen_groups for idx in group)
    if len(indices) != 10:
        raise RuntimeError(f"s4+p6 selector produced {len(indices)} orbitals instead of 10")
    return indices, [len(group) for group in chosen_groups]


def prepare_lithium_integrals(basis_name: str = "cc-pVTZ", n_spatial: int = 10):
    """Generate symmetry-complete ROHF orbital integrals for neutral lithium."""
    try:
        from pyscf import ao2mo, gto, scf
    except ImportError as exc:  # pragma: no cover - optional q1/a1 CI
        raise RuntimeError("OES-A1 lithium integral generation requires the 'q1' extra (PySCF)") from exc

    mol = gto.M(
        atom="Li 0 0 0",
        basis=basis_name,
        unit="Bohr",
        charge=0,
        spin=1,
        verbose=0,
    )
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-12
    rohf_energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError("lithium ROHF did not converge")
    if n_spatial != 10:
        raise ValueError("canonical lithium selector is fixed at ten spatial orbitals")

    selected, group_sizes = select_atomic_s4p6_indices(mf.mo_energy)
    coeff = np.asarray(mf.mo_coeff[:, selected], dtype=float)
    h1 = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.kernel(mol, coeff, compact=False).reshape((n_spatial,) * 4)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([coeff.T @ dip_ao[k] @ coeff for k in range(3)])
    return mol, h1, eri, dip, rohf_energy, selected, group_sizes


def pyscf_lithium_fci_reference(
    h1: np.ndarray,
    eri: np.ndarray,
    ecore: float = 0.0,
) -> float:
    """Independent FCI energy in the N_alpha=2, N_beta=1 sector."""
    try:
        from pyscf import fci
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-A1 FCI reference requires the 'q1' extra (PySCF)") from exc

    norb = h1.shape[0]
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-12
    energy, _ = solver.kernel(h1, eri, norb, (2, 1), ecore=ecore)
    return float(energy)


def first_bright_manifold(
    evals: np.ndarray,
    evecs: np.ndarray,
    basis: Sequence[int],
    dip: np.ndarray,
    n_spin: int = 20,
    f_threshold: float = 1e-6,
    degeneracy_tol_ev: float = 1e-5,
    search_states: int = 80,
):
    """Resolve the first dipole-bright manifold from the M_S=+1/2 FCI sector."""
    ground = evecs[:, 0]
    e0 = float(evals[0])
    rows = []
    for k in range(1, min(search_states, len(evals))):
        t_spin = transition_one_rdm(evecs[:, k], ground, basis, n_spin)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip[a] * t_space) for a in range(3)], dtype=complex)
        de = float(evals[k] - e0)
        f = (2.0 / 3.0) * de * float(np.sum(np.abs(mu) ** 2))
        rows.append((k, de * HARTREE_TO_EV, f))

    bright = [row for row in rows if row[2] > f_threshold]
    if not bright:
        raise RuntimeError("lithium A1 could not resolve a dipole-bright excited state")
    first_ev = min(row[1] for row in bright)
    manifold = [row for row in bright if abs(row[1] - first_ev) < degeneracy_tol_ev]
    return manifold


def run_lithium_a1(basis_name: str = "cc-pVTZ", n_spatial: int = 10) -> LithiumA1Result:
    if n_spatial != 10:
        raise ValueError("OES-A1 canonical gate is fixed at 10 spatial / 20 spin orbitals")

    mol, h1, eri, dip, rohf_energy, selected, group_sizes = prepare_lithium_integrals(
        basis_name=basis_name,
        n_spatial=n_spatial,
    )
    n_spin = 2 * n_spatial

    # Build the complete N=3 Hamiltonian first. The subsequent M_S slice is an
    # exact symmetry reduction of this same 20-mode operator.
    H_full_n3, basis_n3 = build_sector_hamiltonian(
        h1,
        eri,
        n_electrons=3,
        ecore=float(mol.energy_nuc()),
    )
    ms_indices = spin_sector_indices(basis_n3, n_spatial, 2, 1)
    if len(ms_indices) != 450:
        raise RuntimeError(f"unexpected lithium M_S=+1/2 dimension {len(ms_indices)}")
    H = H_full_n3[np.ix_(ms_indices, ms_indices)]
    basis = [basis_n3[i] for i in ms_indices]
    evals, evecs = np.linalg.eigh(H)

    e_oes = float(evals[0])
    e_ref = pyscf_lithium_fci_reference(h1, eri, ecore=float(mol.energy_nuc()))
    manifold = first_bright_manifold(evals, evecs, basis, dip, n_spin=n_spin)
    bright_energies = [row[1] for row in manifold]

    return LithiumA1Result(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        active_protocol="ATOMIC-SYMMETRY-COMPLETE-S4-P3-P3-20Q",
        selected_mo_indices=tuple(selected),
        selected_group_sizes=tuple(group_sizes),
        n_spatial_orbitals=n_spatial,
        n_spin_orbitals=n_spin,
        n_electrons=3,
        full_qubit_dimension=full_space_dimension(n_spin),
        fixed_particle_dimension=sector_dimension(n_spin, 3),
        ms_half_dimension=len(ms_indices),
        rohf_energy_hartree=rohf_energy,
        oes_fci_energy_hartree=e_oes,
        pyscf_fci_energy_hartree=e_ref,
        fci_delta_hartree=e_oes - e_ref,
        first_bright_excitation_ev=float(np.mean(bright_energies)),
        first_bright_degeneracy=len(manifold),
        first_bright_spread_ev=float(max(bright_energies) - min(bright_energies)),
        first_bright_f_sum=float(np.sum([row[2] for row in manifold])),
    )
