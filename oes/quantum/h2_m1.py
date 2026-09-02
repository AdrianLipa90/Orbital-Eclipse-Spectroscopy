"""OES-M1: H2 as the first fixed-20Q two-centre molecular control.

H2 keeps the Q1/A1 register budget (10 spatial / 20 spin orbitals, two
electrons) but replaces atomic spherical flavour blocks with molecular
symmetry classes. The active space is selected only as a union of complete
canonical energy-degeneracy blocks; no experimental constant enters the
selection.

The module validates four new molecular observables:
- inversion parity (g/u) of many-electron states,
- dipole selection between opposite parity sectors,
- equal electron-nuclear exposure of the two equivalent H centres,
- a Born-Oppenheimer potential minimum and harmonic curvature.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Sequence, Tuple

import numpy as np

from .fermions import build_sector_hamiltonian, occupied, one_rdm, sector_dimension, transition_one_rdm
from .helium_q1 import HARTREE_TO_EV, pyscf_fci_reference, spatial_transition_rdm, spin_squared_matrix
from .lithium_a1 import energy_degeneracy_groups

HARTREE_TO_WAVENUMBER_CM = 219_474.631_363_20
PROTON_ELECTRON_MASS_RATIO = 1836.152_673_426
ANGSTROM_TO_BOHR = 1.889_726_125_457_828_1


@dataclass(frozen=True)
class H2PointResult:
    backend: str
    basis_name: str
    bond_bohr: float
    active_protocol: str
    selected_mo_indices: Tuple[int, ...]
    selected_group_sizes: Tuple[int, ...]
    selected_irreps: Tuple[str, ...]
    n_spatial_orbitals: int
    n_spin_orbitals: int
    fixed_particle_dimension: int
    rhf_energy_hartree: float
    oes_fci_energy_hartree: float
    pyscf_fci_energy_hartree: float
    fci_delta_hartree: float
    ground_parity: float
    ground_s2: float
    center_a_rinv_exposure: float
    center_b_rinv_exposure: float
    center_exposure_difference: float
    first_bright_excitation_ev: float
    first_bright_degeneracy: int
    first_bright_parity: float
    first_bright_f_sum: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class H2CurveResult:
    backend: str
    basis_name: str
    grid_bohr: Tuple[float, ...]
    energies_hartree: Tuple[float, ...]
    fitted_equilibrium_bohr: float
    fitted_equilibrium_angstrom: float
    fitted_curvature_hartree_per_bohr2: float
    harmonic_wavenumber_cm: float
    dissociation_probe_bohr: float
    dissociation_probe_energy_hartree: float
    electronic_well_depth_probe_ev: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _select_complete_energy_blocks(energies: Sequence[float], target_orbitals: int = 10):
    """Choose a low-energy union of complete degeneracy groups totalling target."""
    groups = energy_degeneracy_groups(energies, absolute_tolerance=2e-8, relative_tolerance=2e-8)
    dp = {0: (0, ())}
    for gi, group in enumerate(groups):
        size = len(group)
        cost = sum(group)
        for total, (score, chosen) in list(dp.items()):
            new_total = total + size
            if new_total > target_orbitals:
                continue
            candidate = (score + cost, chosen + (gi,))
            if new_total not in dp or candidate[0] < dp[new_total][0]:
                dp[new_total] = candidate
    if target_orbitals not in dp:
        raise RuntimeError(
            f"no complete-degeneracy active space of size {target_orbitals}; group sizes={[len(g) for g in groups]}"
        )
    chosen_groups = [groups[i] for i in dp[target_orbitals][1]]
    indices = sorted(i for group in chosen_groups for i in group)
    if 0 not in indices:
        raise RuntimeError("complete-block selector lost the occupied bonding ground orbital")
    return indices, [len(group) for group in chosen_groups], groups


def _orbital_inversion_parities(irreps: Sequence[str]) -> np.ndarray:
    out = []
    for label in irreps:
        low = str(label).lower()
        if "g" in low:
            out.append(1)
        elif "u" in low:
            out.append(-1)
        else:
            raise ValueError(f"D2h orbital label lacks g/u parity: {label}")
    return np.asarray(out, dtype=int)


def determinant_inversion_parity(det: int, n_spin: int, spatial_parities: Sequence[int]) -> int:
    parity = 1
    for mode in occupied(int(det), n_spin):
        parity *= int(spatial_parities[mode // 2])
    return parity


def state_inversion_parity(state: np.ndarray, basis: Sequence[int], spatial_parities: Sequence[int]) -> float:
    state = np.asarray(state, dtype=complex)
    n_spin = 2 * len(spatial_parities)
    values = np.asarray(
        [determinant_inversion_parity(det, n_spin, spatial_parities) for det in basis],
        dtype=float,
    )
    return float(np.sum(np.abs(state) ** 2 * values).real)


def prepare_h2_integrals(bond_bohr: float, basis_name: str = "cc-pVTZ", n_spatial: int = 10):
    try:
        from pyscf import ao2mo, gto, scf, symm
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("OES-M1 H2 requires the 'q1' extra (PySCF)") from exc
    if bond_bohr <= 0:
        raise ValueError("bond length must be positive")
    if n_spatial != 10:
        raise ValueError("OES-M1 canonical H2 gate is fixed at ten spatial orbitals")

    z = 0.5 * float(bond_bohr)
    mol = gto.M(
        atom=f"H 0 0 {-z}; H 0 0 {z}",
        unit="Bohr",
        basis=basis_name,
        charge=0,
        spin=0,
        symmetry="D2h",
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    rhf_energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"H2 RHF did not converge at R={bond_bohr} bohr")

    labels_all = list(symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mf.mo_coeff))
    selected, group_sizes, groups = _select_complete_energy_blocks(mf.mo_energy, target_orbitals=n_spatial)
    coeff = np.asarray(mf.mo_coeff[:, selected], dtype=float)
    labels = [str(labels_all[i]) for i in selected]
    parities = _orbital_inversion_parities(labels)

    h1 = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.kernel(mol, coeff, compact=False).reshape((n_spatial,) * 4)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([coeff.T @ dip_ao[k] @ coeff for k in range(3)])

    rinv = []
    for atom_index in (0, 1):
        with mol.with_rinv_origin(mol.atom_coord(atom_index)):
            op_ao = mol.intor("int1e_rinv", hermi=1)
        rinv.append(coeff.T @ op_ao @ coeff)

    return {
        "mol": mol,
        "h1": h1,
        "eri": eri,
        "dip": dip,
        "rinv": tuple(rinv),
        "rhf_energy": rhf_energy,
        "selected": selected,
        "group_sizes": group_sizes,
        "all_group_sizes": [len(group) for group in groups],
        "irreps": labels,
        "parities": parities,
    }


def _first_bright_manifold(evals, evecs, basis, dip, parities, threshold=1e-7, degeneracy_tol_ev=2e-6):
    ground = evecs[:, 0]
    e0 = float(evals[0])
    rows = []
    for k in range(1, min(100, len(evals))):
        t_spin = transition_one_rdm(evecs[:, k], ground, basis, 20)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip[a] * t_space) for a in range(3)], dtype=complex)
        de = float(evals[k] - e0)
        f = (2.0 / 3.0) * de * float(np.sum(np.abs(mu) ** 2))
        parity = state_inversion_parity(evecs[:, k], basis, parities)
        rows.append((k, de * HARTREE_TO_EV, f, parity))
    bright = [row for row in rows if row[2] > threshold]
    if not bright:
        raise RuntimeError("H2 M1 found no dipole-bright state in searched active spectrum")
    first_ev = min(row[1] for row in bright)
    return [row for row in bright if abs(row[1] - first_ev) < degeneracy_tol_ev]


def run_h2_point(bond_bohr: float, basis_name: str = "cc-pVTZ") -> H2PointResult:
    data = prepare_h2_integrals(bond_bohr, basis_name=basis_name, n_spatial=10)
    H, basis = build_sector_hamiltonian(
        data["h1"], data["eri"], n_electrons=2, ecore=float(data["mol"].energy_nuc())
    )
    evals, evecs = np.linalg.eigh(H)
    e_oes = float(evals[0])
    e_ref = pyscf_fci_reference(
        data["h1"], data["eri"], n_electrons=2, ecore=float(data["mol"].energy_nuc())
    )

    ground = evecs[:, 0]
    ground_parity = state_inversion_parity(ground, basis, data["parities"])
    s2mat = spin_squared_matrix(10, basis)
    ground_s2 = float(np.real(np.vdot(ground, s2mat @ ground)))
    gamma = spatial_transition_rdm(one_rdm(ground, basis, 20))
    exposures = [float(np.sum(op * gamma).real) for op in data["rinv"]]

    bright = _first_bright_manifold(evals, evecs, basis, data["dip"], data["parities"])
    bright_parities = [row[3] for row in bright]
    if max(bright_parities) - min(bright_parities) > 1e-7:
        raise RuntimeError("first bright H2 manifold mixes inversion parity")

    return H2PointResult(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        bond_bohr=float(bond_bohr),
        active_protocol="MOLECULAR-COMPLETE-DEGENERACY-BLOCKS-20Q",
        selected_mo_indices=tuple(int(i) for i in data["selected"]),
        selected_group_sizes=tuple(int(i) for i in data["group_sizes"]),
        selected_irreps=tuple(data["irreps"]),
        n_spatial_orbitals=10,
        n_spin_orbitals=20,
        fixed_particle_dimension=sector_dimension(20, 2),
        rhf_energy_hartree=float(data["rhf_energy"]),
        oes_fci_energy_hartree=e_oes,
        pyscf_fci_energy_hartree=e_ref,
        fci_delta_hartree=e_oes - e_ref,
        ground_parity=ground_parity,
        ground_s2=ground_s2,
        center_a_rinv_exposure=exposures[0],
        center_b_rinv_exposure=exposures[1],
        center_exposure_difference=exposures[0] - exposures[1],
        first_bright_excitation_ev=float(np.mean([row[1] for row in bright])),
        first_bright_degeneracy=len(bright),
        first_bright_parity=float(np.mean(bright_parities)),
        first_bright_f_sum=float(np.sum([row[2] for row in bright])),
    )


def run_h2_curve(
    basis_name: str = "cc-pVTZ",
    grid_bohr: Sequence[float] = (1.25, 1.325, 1.40, 1.475, 1.55),
    dissociation_probe_bohr: float = 6.0,
) -> H2CurveResult:
    grid = np.asarray(grid_bohr, dtype=float)
    if len(grid) < 3:
        raise ValueError("H2 curve requires at least three points")
    energies = np.asarray([run_h2_point(float(r), basis_name).oes_fci_energy_hartree for r in grid])
    a, b, c = np.polyfit(grid, energies, 2)
    if a <= 0:
        raise RuntimeError("quadratic H2 local potential fit has non-positive curvature")
    r_eq = float(-b / (2.0 * a))
    if not (float(np.min(grid)) <= r_eq <= float(np.max(grid))):
        raise RuntimeError(f"fitted H2 minimum {r_eq} bohr lies outside local scan")
    curvature = float(2.0 * a)
    reduced_nuclear_mass = 0.5 * PROTON_ELECTRON_MASS_RATIO
    harmonic_au = float(np.sqrt(curvature / reduced_nuclear_mass))
    harmonic_cm = harmonic_au * HARTREE_TO_WAVENUMBER_CM
    fit_min = float(a * r_eq * r_eq + b * r_eq + c)

    diss = run_h2_point(float(dissociation_probe_bohr), basis_name).oes_fci_energy_hartree
    depth_ev = float((diss - fit_min) * HARTREE_TO_EV)
    return H2CurveResult(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        grid_bohr=tuple(float(x) for x in grid),
        energies_hartree=tuple(float(x) for x in energies),
        fitted_equilibrium_bohr=r_eq,
        fitted_equilibrium_angstrom=r_eq / ANGSTROM_TO_BOHR,
        fitted_curvature_hartree_per_bohr2=curvature,
        harmonic_wavenumber_cm=harmonic_cm,
        dissociation_probe_bohr=float(dissociation_probe_bohr),
        dissociation_probe_energy_hartree=float(diss),
        electronic_well_depth_probe_ev=depth_ev,
    )
