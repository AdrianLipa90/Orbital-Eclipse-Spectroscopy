"""OES-M2: LiH as the first heteronuclear inorganic fixed-20Q molecule.

LiH has four electrons in the same 10-spatial / 20-spin-orbital register.  The
complete N=4 sector has C(20,4)=4845 determinants.  Exact spin conservation
allows reference simulation in the M_S=0 subset N_alpha=N_beta=2, dimension
C(10,2)^2=2025, while retaining the same 20-mode Hamiltonian encoding.

The molecular active space is selected as a union of complete canonical energy
degeneracy groups.  No experimental constant enters orbital selection.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import comb
from typing import Dict, Sequence, Tuple

import numpy as np

from .determinant_subspace import build_determinant_subspace_hamiltonian
from .fermions import determinant_basis, one_rdm, sector_dimension
from .h2_m1 import ANGSTROM_TO_BOHR, HARTREE_TO_WAVENUMBER_CM, _select_complete_energy_blocks
from .helium_q1 import HARTREE_TO_EV, spatial_transition_rdm
from .lithium_a1 import spin_sector_indices

DEBYE_PER_E_BOHR = 2.541_746_473
ATOMIC_MASS_UNIT_PER_ELECTRON_MASS = 1822.888_486_209
LI7_ATOMIC_MASS_U = 7.016_003_436_6
H1_ATOMIC_MASS_U = 1.007_825_032_23
LI7_NUCLEAR_MASS_ME = LI7_ATOMIC_MASS_U * ATOMIC_MASS_UNIT_PER_ELECTRON_MASS - 3.0
H1_NUCLEAR_MASS_ME = H1_ATOMIC_MASS_U * ATOMIC_MASS_UNIT_PER_ELECTRON_MASS - 1.0
LIH_REDUCED_NUCLEAR_MASS_ME = (
    LI7_NUCLEAR_MASS_ME * H1_NUCLEAR_MASS_ME / (LI7_NUCLEAR_MASS_ME + H1_NUCLEAR_MASS_ME)
)


@dataclass(frozen=True)
class LiHPointResult:
    backend: str
    basis_name: str
    bond_bohr: float
    active_protocol: str
    selected_mo_indices: Tuple[int, ...]
    selected_group_sizes: Tuple[int, ...]
    selected_irreps: Tuple[str, ...]
    n_spatial_orbitals: int
    n_spin_orbitals: int
    n_electrons: int
    full_fixed_particle_dimension: int
    ms_zero_dimension: int
    rhf_energy_hartree: float
    oes_fci_energy_hartree: float
    pyscf_fci_energy_hartree: float
    fci_delta_hartree: float
    li_rinv_exposure: float
    h_rinv_exposure: float
    exposure_difference: float
    dipole_e_bohr: Tuple[float, float, float]
    dipole_debye: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class LiHCurveResult:
    backend: str
    basis_name: str
    grid_bohr: Tuple[float, ...]
    energies_hartree: Tuple[float, ...]
    fitted_equilibrium_bohr: float
    fitted_equilibrium_angstrom: float
    fitted_curvature_hartree_per_bohr2: float
    harmonic_wavenumber_cm: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def prepare_lih_integrals(bond_bohr: float, basis_name: str = "cc-pVTZ", n_spatial: int = 10):
    try:
        from pyscf import ao2mo, gto, scf, symm
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M2 LiH requires the 'q1' extra (PySCF)") from exc
    if bond_bohr <= 0:
        raise ValueError("bond length must be positive")
    if n_spatial != 10:
        raise ValueError("OES-M2 canonical LiH gate is fixed at ten spatial orbitals")

    z = 0.5 * float(bond_bohr)
    mol = gto.M(
        atom=f"Li 0 0 {-z}; H 0 0 {z}",
        unit="Bohr",
        basis=basis_name,
        charge=0,
        spin=0,
        symmetry=True,
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    rhf_energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"LiH RHF did not converge at R={bond_bohr} bohr")

    labels_all = list(symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mf.mo_coeff))
    selected, group_sizes, all_groups = _select_complete_energy_blocks(mf.mo_energy, target_orbitals=10)
    coeff = np.asarray(mf.mo_coeff[:, selected], dtype=float)
    labels = [str(labels_all[i]) for i in selected]

    h1 = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.kernel(mol, coeff, compact=False).reshape((10,) * 4)
    pos_ao = mol.intor("int1e_r", comp=3, hermi=1)
    pos = np.stack([coeff.T @ pos_ao[k] @ coeff for k in range(3)])

    rinv = []
    for atom_index in (0, 1):
        with mol.with_rinv_origin(mol.atom_coord(atom_index)):
            op_ao = mol.intor("int1e_rinv", hermi=1)
        rinv.append(coeff.T @ op_ao @ coeff)

    return {
        "mol": mol,
        "h1": h1,
        "eri": eri,
        "position": pos,
        "rinv": tuple(rinv),
        "rhf_energy": rhf_energy,
        "selected": selected,
        "group_sizes": group_sizes,
        "all_group_sizes": [len(g) for g in all_groups],
        "irreps": labels,
    }


def lih_ms_zero_basis(n_spatial: int = 10):
    full = determinant_basis(2 * n_spatial, 4)
    indices = spin_sector_indices(full, n_spatial, 2, 2)
    basis = tuple(full[i] for i in indices)
    expected = comb(n_spatial, 2) ** 2
    if len(basis) != expected:
        raise RuntimeError(f"unexpected LiH M_S=0 dimension {len(basis)} != {expected}")
    return basis


def pyscf_lih_fci_reference(h1: np.ndarray, eri: np.ndarray, ecore: float = 0.0) -> float:
    try:
        from pyscf import fci
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M2 FCI reference requires the 'q1' extra (PySCF)") from exc
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-12
    energy, _ = solver.kernel(h1, eri, h1.shape[0], (2, 2), ecore=ecore)
    return float(energy)


def _lowest_eigenpair(H: np.ndarray):
    from scipy.linalg import eigh

    values, vectors = eigh(H, subset_by_index=[0, 0], driver="evr")
    return float(values[0]), np.asarray(vectors[:, 0])


def _lowest_eigenvalue(H: np.ndarray) -> float:
    from scipy.linalg import eigh

    value = eigh(H, subset_by_index=[0, 0], eigvals_only=True, driver="evr")
    return float(value[0])


def run_lih_point(bond_bohr: float, basis_name: str = "cc-pVTZ") -> LiHPointResult:
    data = prepare_lih_integrals(bond_bohr, basis_name=basis_name)
    basis = lih_ms_zero_basis(10)
    H, basis_out = build_determinant_subspace_hamiltonian(
        data["h1"], data["eri"], basis, ecore=float(data["mol"].energy_nuc())
    )
    e_oes, ground = _lowest_eigenpair(H)
    e_ref = pyscf_lih_fci_reference(
        data["h1"], data["eri"], ecore=float(data["mol"].energy_nuc())
    )

    gamma_spin = one_rdm(ground, basis_out, 20)
    gamma = spatial_transition_rdm(gamma_spin)
    exposures = [float(np.sum(op * gamma).real) for op in data["rinv"]]

    electron_position = np.array(
        [float(np.sum(data["position"][k] * gamma).real) for k in range(3)],
        dtype=float,
    )
    nuclear = np.zeros(3, dtype=float)
    for atom_index in range(data["mol"].natm):
        nuclear += float(data["mol"].atom_charge(atom_index)) * np.asarray(data["mol"].atom_coord(atom_index))
    dipole = nuclear - electron_position
    dipole_debye = float(np.linalg.norm(dipole) * DEBYE_PER_E_BOHR)

    return LiHPointResult(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        bond_bohr=float(bond_bohr),
        active_protocol="HETERONUCLEAR-COMPLETE-DEGENERACY-BLOCKS-20Q",
        selected_mo_indices=tuple(int(i) for i in data["selected"]),
        selected_group_sizes=tuple(int(i) for i in data["group_sizes"]),
        selected_irreps=tuple(data["irreps"]),
        n_spatial_orbitals=10,
        n_spin_orbitals=20,
        n_electrons=4,
        full_fixed_particle_dimension=sector_dimension(20, 4),
        ms_zero_dimension=len(basis_out),
        rhf_energy_hartree=float(data["rhf_energy"]),
        oes_fci_energy_hartree=e_oes,
        pyscf_fci_energy_hartree=e_ref,
        fci_delta_hartree=e_oes - e_ref,
        li_rinv_exposure=exposures[0],
        h_rinv_exposure=exposures[1],
        exposure_difference=exposures[0] - exposures[1],
        dipole_e_bohr=tuple(float(x) for x in dipole),
        dipole_debye=dipole_debye,
    )


def run_lih_ground_energy(bond_bohr: float, basis_name: str = "cc-pVTZ") -> float:
    data = prepare_lih_integrals(bond_bohr, basis_name=basis_name)
    H, _ = build_determinant_subspace_hamiltonian(
        data["h1"], data["eri"], lih_ms_zero_basis(10), ecore=float(data["mol"].energy_nuc())
    )
    return _lowest_eigenvalue(H)


def run_lih_curve(
    basis_name: str = "cc-pVTZ",
    grid_bohr: Sequence[float] = (2.82, 2.92, 3.02, 3.12, 3.22),
) -> LiHCurveResult:
    grid = np.asarray(grid_bohr, dtype=float)
    energies = np.asarray([run_lih_ground_energy(float(r), basis_name) for r in grid])
    a, b, c = np.polyfit(grid, energies, 2)
    if a <= 0:
        raise RuntimeError("LiH local potential fit has non-positive curvature")
    r_eq = float(-b / (2.0 * a))
    if not (float(np.min(grid)) <= r_eq <= float(np.max(grid))):
        raise RuntimeError(f"LiH fitted minimum {r_eq} bohr lies outside local scan")
    curvature = float(2.0 * a)
    omega_au = float(np.sqrt(curvature / LIH_REDUCED_NUCLEAR_MASS_ME))
    return LiHCurveResult(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        grid_bohr=tuple(float(x) for x in grid),
        energies_hartree=tuple(float(x) for x in energies),
        fitted_equilibrium_bohr=r_eq,
        fitted_equilibrium_angstrom=r_eq / ANGSTROM_TO_BOHR,
        fitted_curvature_hartree_per_bohr2=curvature,
        harmonic_wavenumber_cm=omega_au * HARTREE_TO_WAVENUMBER_CM,
    )
