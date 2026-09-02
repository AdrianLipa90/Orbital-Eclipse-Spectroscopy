"""OES-M4: H35Cl spectroscopy in an eight-active-electron fixed-20Q space.

HCl is the first cross-species reuse of the M3.1 eight-electron architecture.
Ten closed-shell core electrons are integrated exactly through G1; the eight
remaining valence electrons occupy ten active spatial orbitals / twenty spin
orbitals and are solved in the exact G2 fixed-Ms sparse sector.

Experimental constants are excluded from geometry search, active-space
selection, Hamiltonian construction, density reconstruction and curve fitting.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Sequence, Tuple

import numpy as np

from .fermions import one_rdm
from .frozen_core import frozen_core_effective_hamiltonian
from .h2_m1 import ANGSTROM_TO_BOHR, HARTREE_TO_WAVENUMBER_CM, _select_complete_energy_blocks
from .helium_q1 import spatial_transition_rdm
from .hf_m3 import _total_one_body_expectation
from .lih_m2 import ATOMIC_MASS_UNIT_PER_ELECTRON_MASS, DEBYE_PER_E_BOHR
from .sparse_ms import build_sparse_fixed_spin_hamiltonian

H1_ATOMIC_MASS_U = 1.007_825_032_23
CL35_ATOMIC_MASS_U = 34.968_852_682
H1_NUCLEAR_MASS_ME = H1_ATOMIC_MASS_U * ATOMIC_MASS_UNIT_PER_ELECTRON_MASS - 1.0
CL35_NUCLEAR_MASS_ME = CL35_ATOMIC_MASS_U * ATOMIC_MASS_UNIT_PER_ELECTRON_MASS - 17.0
HCL35_REDUCED_NUCLEAR_MASS_ME = (
    H1_NUCLEAR_MASS_ME * CL35_NUCLEAR_MASS_ME / (H1_NUCLEAR_MASS_ME + CL35_NUCLEAR_MASS_ME)
)


@dataclass(frozen=True)
class HClPointResult:
    backend: str
    basis_name: str
    bond_bohr: float
    active_protocol: str
    core_mo_indices: Tuple[int, ...]
    active_mo_indices: Tuple[int, ...]
    active_group_sizes: Tuple[int, ...]
    n_total_electrons: int
    n_frozen_electrons: int
    n_active_electrons: int
    active_qubits: int
    fixed_ms_dimension: int
    rhf_energy_hartree: float
    oes_sparse_energy_hartree: float
    pyscf_active_fci_energy_hartree: float
    energy_delta_hartree: float
    active_rdm1_max_delta: float
    cl_rinv_exposure: float
    h_rinv_exposure: float
    exposure_difference: float
    dipole_e_bohr: Tuple[float, float, float]
    dipole_debye: float
    csr_nnz: int
    csr_density: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class HClCurveResult:
    backend: str
    basis_name: str
    rhf_scan_bohr: Tuple[float, ...]
    rhf_scan_energies_hartree: Tuple[float, ...]
    rhf_seed_minimum_bohr: float
    active_grid_bohr: Tuple[float, ...]
    active_energies_hartree: Tuple[float, ...]
    fitted_equilibrium_bohr: float
    fitted_equilibrium_angstrom: float
    fitted_curvature_hartree_per_bohr2: float
    harmonic_wavenumber_cm: float
    rotational_constant_cm: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _build_hcl_molecule(bond_bohr: float, basis_name: str):
    try:
        from pyscf import gto
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M4 HCl requires the q1 extra (PySCF)") from exc
    if bond_bohr <= 0:
        raise ValueError("bond length must be positive")
    z = 0.5 * float(bond_bohr)
    return gto.M(
        atom=f"Cl 0 0 {-z}; H 0 0 {z}",
        unit="Bohr",
        basis=basis_name,
        charge=0,
        spin=0,
        symmetry=True,
        verbose=0,
    )


def run_hcl_rhf_energy(bond_bohr: float, basis_name: str = "cc-pVTZ") -> float:
    try:
        from pyscf import scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M4 HCl requires the q1 extra (PySCF)") from exc
    mol = _build_hcl_molecule(bond_bohr, basis_name)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"HCl RHF did not converge at R={bond_bohr} bohr")
    return energy


def _hcl_active_integrals(bond_bohr: float, basis_name: str = "cc-pVTZ"):
    try:
        from pyscf import ao2mo, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M4 HCl requires the q1 extra (PySCF)") from exc

    mol = _build_hcl_molecule(bond_bohr, basis_name)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    rhf_energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"HCl RHF did not converge at R={bond_bohr} bohr")

    # The first five occupied spatial MOs hold the closed Cl 1s/2s/2p core
    # (ten electrons). Eight valence electrons remain active.
    core = tuple(range(5))
    relative, group_sizes, _ = _select_complete_energy_blocks(
        np.asarray(mf.mo_energy[5:], dtype=float), target_orbitals=10
    )
    active = tuple(int(i + 5) for i in relative)
    if len(active) != 10 or set(core) & set(active):
        raise RuntimeError("HCl core/active partition failed")

    selected = core + active
    coeff = np.asarray(mf.mo_coeff[:, selected], dtype=float)
    n_selected = len(selected)
    h = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.kernel(mol, coeff, compact=False).reshape((n_selected,) * 4)
    reduced = frozen_core_effective_hamiltonian(
        h,
        eri,
        core_indices=tuple(range(5)),
        active_indices=tuple(range(5, 15)),
        nuclear_energy=float(mol.energy_nuc()),
    )

    position_ao = mol.intor("int1e_r", comp=3, hermi=1)
    position = np.stack([coeff.T @ position_ao[k] @ coeff for k in range(3)])
    rinv = []
    for atom_index in range(mol.natm):
        with mol.with_rinv_origin(mol.atom_coord(atom_index)):
            op_ao = mol.intor("int1e_rinv", hermi=1)
        rinv.append(coeff.T @ op_ao @ coeff)

    return {
        "mol": mol,
        "mf": mf,
        "rhf_energy": rhf_energy,
        "core": core,
        "active": active,
        "group_sizes": tuple(int(x) for x in group_sizes),
        "reduced": reduced,
        "position": position,
        "rinv": tuple(rinv),
    }


def run_hcl_ground_energy(bond_bohr: float, basis_name: str = "cc-pVTZ") -> float:
    from scipy.sparse.linalg import eigsh

    data = _hcl_active_integrals(bond_bohr, basis_name)
    reduced = data["reduced"]
    H, _ = build_sparse_fixed_spin_hamiltonian(
        reduced.h1_active,
        reduced.eri_active,
        n_alpha=4,
        n_beta=4,
        ecore=reduced.ecore,
    )
    return float(eigsh(H, k=1, which="SA", return_eigenvectors=False, tol=2e-10, maxiter=15000)[0])


def run_hcl_point(bond_bohr: float, basis_name: str = "cc-pVTZ") -> HClPointResult:
    try:
        from pyscf import fci
        from scipy.sparse.linalg import eigsh
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M4 HCl requires q1 and SciPy") from exc

    data = _hcl_active_integrals(bond_bohr, basis_name)
    reduced = data["reduced"]
    H, basis = build_sparse_fixed_spin_hamiltonian(
        reduced.h1_active,
        reduced.eri_active,
        n_alpha=4,
        n_beta=4,
        ecore=reduced.ecore,
    )
    values, vectors = eigsh(H, k=1, which="SA", return_eigenvectors=True, tol=2e-10, maxiter=15000)
    e_oes = float(values[0])
    ground = np.asarray(vectors[:, 0])

    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 2e-10
    e_ref, ci_ref = solver.kernel(
        reduced.h1_active,
        reduced.eri_active,
        10,
        (4, 4),
        ecore=reduced.ecore,
    )
    e_ref = float(e_ref)
    gamma_ref = np.asarray(solver.make_rdm1(ci_ref, 10, (4, 4)), dtype=float)

    gamma_spin = one_rdm(ground, basis, 20)
    gamma_active = np.asarray(spatial_transition_rdm(gamma_spin).real, dtype=float)
    rdm_delta = float(np.max(np.abs(gamma_active - gamma_ref)))

    exposures = [
        _total_one_body_expectation(op, gamma_active, n_core=5)
        for op in data["rinv"]
    ]
    electron_position = np.asarray([
        _total_one_body_expectation(data["position"][k], gamma_active, n_core=5)
        for k in range(3)
    ])
    nuclear = np.zeros(3, dtype=float)
    for atom_index in range(data["mol"].natm):
        nuclear += float(data["mol"].atom_charge(atom_index)) * np.asarray(data["mol"].atom_coord(atom_index))
    dipole = nuclear - electron_position

    return HClPointResult(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        bond_bohr=float(bond_bohr),
        active_protocol="HCL-FROZEN-10E-COMPLETE-DEGENERACY-BLOCKS-20Q",
        core_mo_indices=tuple(int(i) for i in data["core"]),
        active_mo_indices=tuple(int(i) for i in data["active"]),
        active_group_sizes=data["group_sizes"],
        n_total_electrons=18,
        n_frozen_electrons=10,
        n_active_electrons=8,
        active_qubits=20,
        fixed_ms_dimension=len(basis),
        rhf_energy_hartree=float(data["rhf_energy"]),
        oes_sparse_energy_hartree=e_oes,
        pyscf_active_fci_energy_hartree=e_ref,
        energy_delta_hartree=e_oes - e_ref,
        active_rdm1_max_delta=rdm_delta,
        cl_rinv_exposure=float(exposures[0]),
        h_rinv_exposure=float(exposures[1]),
        exposure_difference=float(exposures[0] - exposures[1]),
        dipole_e_bohr=tuple(float(x) for x in dipole),
        dipole_debye=float(np.linalg.norm(dipole) * DEBYE_PER_E_BOHR),
        csr_nnz=int(H.nnz),
        csr_density=float(H.nnz / (len(basis) ** 2)),
    )


def rotational_constant_cm(bond_bohr: float) -> float:
    if bond_bohr <= 0:
        raise ValueError("bond length must be positive")
    return HARTREE_TO_WAVENUMBER_CM / (2.0 * HCL35_REDUCED_NUCLEAR_MASS_ME * bond_bohr**2)


def run_hcl_curve(
    basis_name: str = "cc-pVTZ",
    rhf_scan_bohr: Sequence[float] = (1.90, 2.10, 2.30, 2.50, 2.70, 2.90, 3.10),
    active_half_width_bohr: float = 0.12,
) -> HClCurveResult:
    """Blind HCl curve: RHF locates a seed; 8e/20Q active FCI refines it."""
    rhf_grid = np.asarray(rhf_scan_bohr, dtype=float)
    if len(rhf_grid) < 5 or not np.all(np.diff(rhf_grid) > 0):
        raise ValueError("HCl RHF scan requires an increasing grid with at least five points")
    rhf_energies = np.asarray([run_hcl_rhf_energy(float(r), basis_name) for r in rhf_grid])
    seed_index = int(np.argmin(rhf_energies))
    if seed_index in (0, len(rhf_grid) - 1):
        raise RuntimeError("blind HCl RHF minimum lies at scan boundary")
    seed = float(rhf_grid[seed_index])

    active_grid = seed + np.linspace(-active_half_width_bohr, active_half_width_bohr, 5)
    active_energies = np.asarray([run_hcl_ground_energy(float(r), basis_name) for r in active_grid])
    a, b, c = np.polyfit(active_grid, active_energies, 2)
    if a <= 0:
        raise RuntimeError("HCl active local potential fit has non-positive curvature")
    r_eq = float(-b / (2.0 * a))
    if not (float(active_grid[0]) <= r_eq <= float(active_grid[-1])):
        raise RuntimeError("HCl fitted active minimum lies outside local active scan")
    curvature = float(2.0 * a)
    omega_au = float(np.sqrt(curvature / HCL35_REDUCED_NUCLEAR_MASS_ME))

    return HClCurveResult(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        rhf_scan_bohr=tuple(float(x) for x in rhf_grid),
        rhf_scan_energies_hartree=tuple(float(x) for x in rhf_energies),
        rhf_seed_minimum_bohr=seed,
        active_grid_bohr=tuple(float(x) for x in active_grid),
        active_energies_hartree=tuple(float(x) for x in active_energies),
        fitted_equilibrium_bohr=r_eq,
        fitted_equilibrium_angstrom=r_eq / ANGSTROM_TO_BOHR,
        fitted_curvature_hartree_per_bohr2=curvature,
        harmonic_wavenumber_cm=omega_au * HARTREE_TO_WAVENUMBER_CM,
        rotational_constant_cm=rotational_constant_cm(r_eq),
    )
