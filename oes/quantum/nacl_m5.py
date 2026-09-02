"""OES-M5: 23Na35Cl spectroscopy in an eight-active-electron fixed-20Q space.

NaCl is a large-core ionic stress test of the fixed-width OES molecular
architecture. Twenty closed-shell core electrons are integrated exactly through
the frozen-core Hamiltonian; the eight remaining valence electrons occupy ten
active spatial orbitals / twenty spin orbitals and are solved in the exact
fixed-Ms sparse sector.

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
from .hcl_m4 import _adaptive_bracketed_quadratic
from .helium_q1 import spatial_transition_rdm
from .hf_m3 import _total_one_body_expectation
from .lih_m2 import ATOMIC_MASS_UNIT_PER_ELECTRON_MASS, DEBYE_PER_E_BOHR
from .sparse_ms import build_sparse_fixed_spin_hamiltonian

NA23_ATOMIC_MASS_U = 22.989_769_282_0
CL35_ATOMIC_MASS_U = 34.968_852_682
NA23_NUCLEAR_MASS_ME = NA23_ATOMIC_MASS_U * ATOMIC_MASS_UNIT_PER_ELECTRON_MASS - 11.0
CL35_NUCLEAR_MASS_ME = CL35_ATOMIC_MASS_U * ATOMIC_MASS_UNIT_PER_ELECTRON_MASS - 17.0
NACL35_REDUCED_NUCLEAR_MASS_ME = (
    NA23_NUCLEAR_MASS_ME * CL35_NUCLEAR_MASS_ME / (NA23_NUCLEAR_MASS_ME + CL35_NUCLEAR_MASS_ME)
)


@dataclass(frozen=True)
class NaClPointResult:
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
    na_rinv_exposure: float
    cl_rinv_exposure: float
    exposure_difference_na_minus_cl: float
    dipole_e_bohr: Tuple[float, float, float]
    dipole_debye: float
    csr_nnz: int
    csr_density: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class NaClCurveResult:
    backend: str
    basis_name: str
    rhf_scan_bohr: Tuple[float, ...]
    rhf_scan_energies_hartree: Tuple[float, ...]
    rhf_seed_minimum_bohr: float
    active_grid_bohr: Tuple[float, ...]
    active_energies_hartree: Tuple[float, ...]
    active_recenter_count: int
    fitted_equilibrium_bohr: float
    fitted_equilibrium_angstrom: float
    fitted_curvature_hartree_per_bohr2: float
    harmonic_wavenumber_cm: float
    rotational_constant_cm: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _build_nacl_molecule(bond_bohr: float, basis_name: str):
    try:
        from pyscf import gto
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M5 NaCl requires the q1 extra (PySCF)") from exc
    if bond_bohr <= 0:
        raise ValueError("bond length must be positive")
    z = 0.5 * float(bond_bohr)
    return gto.M(
        atom=f"Na 0 0 {-z}; Cl 0 0 {z}",
        unit="Bohr",
        basis=basis_name,
        charge=0,
        spin=0,
        symmetry=True,
        verbose=0,
    )


def run_nacl_rhf_energy(bond_bohr: float, basis_name: str = "cc-pVTZ") -> float:
    try:
        from pyscf import scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M5 NaCl requires the q1 extra (PySCF)") from exc
    mol = _build_nacl_molecule(bond_bohr, basis_name)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"NaCl RHF did not converge at R={bond_bohr} bohr")
    return energy


def _nacl_active_integrals(bond_bohr: float, basis_name: str = "cc-pVTZ"):
    try:
        from pyscf import ao2mo, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M5 NaCl requires the q1 extra (PySCF)") from exc

    mol = _build_nacl_molecule(bond_bohr, basis_name)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    rhf_energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"NaCl RHF did not converge at R={bond_bohr} bohr")

    # Na 1s/2s/2p plus Cl 1s/2s/2p form ten closed-shell spatial core
    # orbitals (twenty electrons). Four occupied valence orbitals and six
    # additional complete-degeneracy virtual directions remain in the fixed
    # ten-spatial-orbital active register.
    core = tuple(range(10))
    relative, group_sizes, _ = _select_complete_energy_blocks(
        np.asarray(mf.mo_energy[10:], dtype=float), target_orbitals=10
    )
    active = tuple(int(i + 10) for i in relative)
    if len(active) != 10 or set(core) & set(active):
        raise RuntimeError("NaCl core/active partition failed")

    selected = core + active
    coeff = np.asarray(mf.mo_coeff[:, selected], dtype=float)
    n_selected = len(selected)
    h = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.kernel(mol, coeff, compact=False).reshape((n_selected,) * 4)
    reduced = frozen_core_effective_hamiltonian(
        h,
        eri,
        core_indices=tuple(range(10)),
        active_indices=tuple(range(10, 20)),
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


def run_nacl_ground_energy(bond_bohr: float, basis_name: str = "cc-pVTZ") -> float:
    from scipy.sparse.linalg import eigsh

    data = _nacl_active_integrals(bond_bohr, basis_name)
    reduced = data["reduced"]
    H, _ = build_sparse_fixed_spin_hamiltonian(
        reduced.h1_active,
        reduced.eri_active,
        n_alpha=4,
        n_beta=4,
        ecore=reduced.ecore,
    )
    return float(eigsh(H, k=1, which="SA", return_eigenvectors=False, tol=2e-10, maxiter=15000)[0])


def run_nacl_point(bond_bohr: float, basis_name: str = "cc-pVTZ") -> NaClPointResult:
    try:
        from pyscf import fci
        from scipy.sparse.linalg import eigsh
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M5 NaCl requires q1 and SciPy") from exc

    data = _nacl_active_integrals(bond_bohr, basis_name)
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
        _total_one_body_expectation(op, gamma_active, n_core=10)
        for op in data["rinv"]
    ]
    electron_position = np.asarray([
        _total_one_body_expectation(data["position"][k], gamma_active, n_core=10)
        for k in range(3)
    ])
    nuclear = np.zeros(3, dtype=float)
    for atom_index in range(data["mol"].natm):
        nuclear += float(data["mol"].atom_charge(atom_index)) * np.asarray(data["mol"].atom_coord(atom_index))
    dipole = nuclear - electron_position

    return NaClPointResult(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        bond_bohr=float(bond_bohr),
        active_protocol="NACL-FROZEN-20E-COMPLETE-DEGENERACY-BLOCKS-20Q",
        core_mo_indices=tuple(int(i) for i in data["core"]),
        active_mo_indices=tuple(int(i) for i in data["active"]),
        active_group_sizes=data["group_sizes"],
        n_total_electrons=28,
        n_frozen_electrons=20,
        n_active_electrons=8,
        active_qubits=20,
        fixed_ms_dimension=len(basis),
        rhf_energy_hartree=float(data["rhf_energy"]),
        oes_sparse_energy_hartree=e_oes,
        pyscf_active_fci_energy_hartree=e_ref,
        energy_delta_hartree=e_oes - e_ref,
        active_rdm1_max_delta=rdm_delta,
        na_rinv_exposure=float(exposures[0]),
        cl_rinv_exposure=float(exposures[1]),
        exposure_difference_na_minus_cl=float(exposures[0] - exposures[1]),
        dipole_e_bohr=tuple(float(x) for x in dipole),
        dipole_debye=float(np.linalg.norm(dipole) * DEBYE_PER_E_BOHR),
        csr_nnz=int(H.nnz),
        csr_density=float(H.nnz / (len(basis) ** 2)),
    )


def rotational_constant_cm(bond_bohr: float) -> float:
    if bond_bohr <= 0:
        raise ValueError("bond length must be positive")
    return HARTREE_TO_WAVENUMBER_CM / (2.0 * NACL35_REDUCED_NUCLEAR_MASS_ME * bond_bohr**2)


def run_nacl_curve(
    basis_name: str = "cc-pVTZ",
    rhf_scan_bohr: Sequence[float] = (3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2),
    active_half_width_bohr: float = 0.18,
) -> NaClCurveResult:
    """Blind NaCl curve: RHF locates a seed; 8e/20Q active FCI refines it."""
    rhf_grid = np.asarray(rhf_scan_bohr, dtype=float)
    if len(rhf_grid) < 5 or not np.all(np.diff(rhf_grid) > 0):
        raise ValueError("NaCl RHF scan requires an increasing grid with at least five points")
    rhf_energies = np.asarray([run_nacl_rhf_energy(float(r), basis_name) for r in rhf_grid])
    seed_index = int(np.argmin(rhf_energies))
    if seed_index in (0, len(rhf_grid) - 1):
        raise RuntimeError("blind NaCl RHF minimum lies at scan boundary")
    seed = float(rhf_grid[seed_index])

    active_grid, active_energies, r_eq, curvature, recenter_count = _adaptive_bracketed_quadratic(
        lambda r: run_nacl_ground_energy(r, basis_name),
        seed_bohr=seed,
        half_width_bohr=active_half_width_bohr,
    )
    omega_au = float(np.sqrt(curvature / NACL35_REDUCED_NUCLEAR_MASS_ME))

    return NaClCurveResult(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        rhf_scan_bohr=tuple(float(x) for x in rhf_grid),
        rhf_scan_energies_hartree=tuple(float(x) for x in rhf_energies),
        rhf_seed_minimum_bohr=seed,
        active_grid_bohr=tuple(float(x) for x in active_grid),
        active_energies_hartree=tuple(float(x) for x in active_energies),
        active_recenter_count=int(recenter_count),
        fitted_equilibrium_bohr=r_eq,
        fitted_equilibrium_angstrom=r_eq / ANGSTROM_TO_BOHR,
        fitted_curvature_hartree_per_bohr2=curvature,
        harmonic_wavenumber_cm=omega_au * HARTREE_TO_WAVENUMBER_CM,
        rotational_constant_cm=rotational_constant_cm(r_eq),
    )
