"""Generic closed-shell molecular TDA spectroscopy for OES.

This is a standard-QM backend adapter over PySCF RHF/TDA. It accepts an
upstream molecular geometry and returns transition energies, length-gauge
transition dipoles, and oscillator strengths. Raw TDA root indices are not
assumed to define state identity across different geometries.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Iterable, Sequence

import numpy as np

HARTREE_TO_EV = 27.211_386_245_981


class MolecularTDAError(ValueError):
    pass


def _finite(value: float, name: str) -> float:
    x = float(value)
    if not math.isfinite(x):
        raise MolecularTDAError(f"{name} must be finite")
    return x


def _integer(value: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise MolecularTDAError(f"{name} must be an integer")
    integer = int(value)
    if integer != value:
        raise MolecularTDAError(f"{name} must be an integer")
    return integer


def _validate_geometry(
    atoms: Sequence[str],
    coordinates_angstrom,
) -> tuple[tuple[str, ...], np.ndarray]:
    symbols = tuple(str(atom).strip() for atom in atoms)
    if not symbols or any(not atom for atom in symbols):
        raise MolecularTDAError("atoms must contain non-empty element symbols")
    coords = np.asarray(coordinates_angstrom, dtype=float)
    if coords.shape != (len(symbols), 3):
        raise MolecularTDAError(
            "coordinates_angstrom must have shape (len(atoms), 3)"
        )
    if not np.all(np.isfinite(coords)):
        raise MolecularTDAError("coordinates_angstrom must be finite")
    return symbols, coords


@dataclass(frozen=True)
class MolecularTDAState:
    root: int
    excitation_hartree: float
    excitation_ev: float
    transition_dipole_au: tuple[float, float, float]
    transition_dipole_norm_au: float
    oscillator_strength_backend: float
    oscillator_strength_oes: float
    oscillator_strength_delta: float

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class MolecularTDAResult:
    backend: str
    species_label: str
    basis_name: str
    charge: int
    spin: int
    atoms: tuple[str, ...]
    coordinates_angstrom: tuple[tuple[float, float, float], ...]
    rhf_energy_hartree: float
    states: tuple[MolecularTDAState, ...]
    max_oscillator_strength_delta: float
    state_tracking_status: str = "SINGLE_GEOMETRY_NO_TRACKING_REQUIRED"

    def as_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["schema"] = "OES_MOLECULAR_RHF_TDA_V0_1"
        return data


@dataclass(frozen=True)
class MolecularTDARuntime:
    """Runtime-only payload for cross-geometry state tracking.

    transition_densities_ao is the spin-summed singlet-TDA transition density
    mapped into the AO basis as 2 C_occ X C_vir^†. OES cross-checks this matrix
    against the backend length-gauge transition dipole. T4 Gram whitening makes
    state-continuity results insensitive to the common normalization.
    """

    result: MolecularTDAResult
    molecule: Any
    ao_overlap: np.ndarray
    transition_densities_ao: tuple[np.ndarray, ...]
    transition_dipole_reconstruction_deltas: tuple[float, ...]


@dataclass(frozen=True)
class OrientationTDAPoint:
    angle_rad: float
    result: MolecularTDAResult

    def as_dict(self) -> dict[str, object]:
        return {
            "angle_rad": self.angle_rad,
            "result": self.result.as_dict(),
        }


@dataclass(frozen=True)
class OrientationTDAScan:
    points: tuple[OrientationTDAPoint, ...]
    state_tracking_status: str = "RAW_ROOT_INDEX_UNTRACKED"

    def as_dict(self) -> dict[str, object]:
        return {
            "schema": "OES_ORIENTATION_TDA_SCAN_V0_1",
            "state_tracking_status": self.state_tracking_status,
            "points": [point.as_dict() for point in self.points],
        }


def run_closed_shell_tda_runtime(
    *,
    species_label: str,
    atoms: Sequence[str],
    coordinates_angstrom,
    basis_name: str = "cc-pVDZ",
    charge: int = 0,
    spin: int = 0,
    nstates: int = 6,
    scf_conv_tol: float = 1.0e-10,
    oscillator_crosscheck_atol: float = 1.0e-10,
    oscillator_crosscheck_rtol: float = 1.0e-10,
) -> MolecularTDAResult:
    """Run a real RHF/TDA line calculation for one fixed geometry."""
    if not species_label:
        raise MolecularTDAError("species_label is required")
    if not isinstance(basis_name, str) or not basis_name.strip():
        raise MolecularTDAError("basis_name must be a non-empty string")
    q = _integer(charge, "charge")
    s = _integer(spin, "spin")
    if s != 0:
        raise MolecularTDAError("v0.1 supports closed-shell spin=0 only")
    roots = _integer(nstates, "nstates")
    if roots < 1:
        raise MolecularTDAError("nstates must be positive")
    conv_tol = _finite(scf_conv_tol, "scf_conv_tol")
    atol = _finite(
        oscillator_crosscheck_atol,
        "oscillator_crosscheck_atol",
    )
    rtol = _finite(
        oscillator_crosscheck_rtol,
        "oscillator_crosscheck_rtol",
    )
    if conv_tol <= 0.0:
        raise MolecularTDAError("scf_conv_tol must be positive")
    if atol < 0.0 or rtol < 0.0:
        raise MolecularTDAError(
            "oscillator cross-check tolerances must be non-negative"
        )

    symbols, coords = _validate_geometry(atoms, coordinates_angstrom)

    try:
        from pyscf import gto, scf, tdscf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "molecular TDA requires the OES q1 extra (PySCF)"
        ) from exc

    atom_spec = [
        (symbol, tuple(float(value) for value in xyz))
        for symbol, xyz in zip(symbols, coords)
    ]
    mol = gto.M(
        atom=atom_spec,
        unit="Angstrom",
        basis=basis_name,
        charge=q,
        spin=s,
        symmetry=False,
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = conv_tol
    mf.max_cycle = 150
    rhf_energy = float(mf.kernel())
    if not mf.converged or not math.isfinite(rhf_energy):
        raise RuntimeError("RHF did not converge for molecular TDA geometry")

    td = tdscf.TDA(mf)
    td.singlet = True
    td.nstates = roots
    excitation, xy = td.kernel()
    excitation = np.atleast_1d(np.asarray(excitation, dtype=float))
    if excitation.size < 1 or not np.all(np.isfinite(excitation)):
        raise RuntimeError("TDA returned no finite excitation energies")

    converged = getattr(td, "converged", None)
    if converged is not None:
        flags = np.atleast_1d(np.asarray(converged, dtype=bool))
        if flags.size and not np.all(flags):
            raise RuntimeError("one or more requested TDA roots did not converge")

    transition_dipoles = np.asarray(
        td.transition_dipole(xy=xy),
        dtype=float,
    )
    oscillator_backend = np.atleast_1d(
        np.asarray(
            td.oscillator_strength(
                e=excitation,
                xy=xy,
                gauge="length",
            ),
            dtype=float,
        )
    )
    if transition_dipoles.ndim == 1:
        transition_dipoles = transition_dipoles.reshape(1, -1)
    expected_shape = (excitation.size, 3)
    if transition_dipoles.shape != expected_shape:
        raise RuntimeError(
            f"unexpected transition-dipole shape {transition_dipoles.shape}"
        )
    if oscillator_backend.shape != (excitation.size,):
        raise RuntimeError("unexpected oscillator-strength shape")
    if (
        not np.all(np.isfinite(transition_dipoles))
        or not np.all(np.isfinite(oscillator_backend))
    ):
        raise RuntimeError("TDA transition observables must be finite")

    states = []
    max_delta = 0.0
    for index, energy in enumerate(excitation):
        de = float(energy)
        if de <= 0.0:
            raise RuntimeError("TDA excitation energy must be positive")
        mu = np.asarray(transition_dipoles[index], dtype=float)
        mu2 = float(np.dot(mu, mu))
        f_oes = (2.0 / 3.0) * de * mu2
        f_backend = float(oscillator_backend[index])
        delta = abs(f_backend - f_oes)
        max_delta = max(max_delta, delta)
        if not math.isclose(
            f_backend,
            f_oes,
            rel_tol=rtol,
            abs_tol=atol,
        ):
            raise RuntimeError(
                "PySCF/OES length-gauge oscillator-strength cross-check failed"
            )
        states.append(
            MolecularTDAState(
                root=index + 1,
                excitation_hartree=de,
                excitation_ev=de * HARTREE_TO_EV,
                transition_dipole_au=tuple(float(x) for x in mu),
                transition_dipole_norm_au=float(math.sqrt(mu2)),
                oscillator_strength_backend=f_backend,
                oscillator_strength_oes=f_oes,
                oscillator_strength_delta=delta,
            )
        )

    result = MolecularTDAResult(
        backend="PYSCF_RHF_TDA",
        species_label=species_label,
        basis_name=basis_name,
        charge=q,
        spin=s,
        atoms=symbols,
        coordinates_angstrom=tuple(
            tuple(float(value) for value in xyz)
            for xyz in coords
        ),
        rhf_energy_hartree=rhf_energy,
        states=tuple(states),
        max_oscillator_strength_delta=max_delta,
    )

    mo_coeff = np.asarray(mf.mo_coeff, dtype=complex)
    mo_occ = np.asarray(mf.mo_occ, dtype=float)
    occupied = np.where(mo_occ > 0.0)[0]
    virtual = np.where(mo_occ <= 0.0)[0]
    if occupied.size < 1 or virtual.size < 1:
        raise RuntimeError(
            "molecular TDA requires occupied and virtual orbital sectors"
        )
    c_occ = mo_coeff[:, occupied]
    c_vir = mo_coeff[:, virtual]

    transition_densities = []
    if len(xy) != excitation.size:
        raise RuntimeError("TDA amplitude count does not match excitation count")
    for root_xy in xy:
        if not isinstance(root_xy, (tuple, list)) or len(root_xy) < 1:
            raise RuntimeError("unexpected TDA amplitude payload")
        x_amplitude = np.asarray(root_xy[0], dtype=complex)
        expected_x_shape = (occupied.size, virtual.size)
        if x_amplitude.shape != expected_x_shape:
            raise RuntimeError(
                f"unexpected TDA X-amplitude shape {x_amplitude.shape}; "
                f"expected {expected_x_shape}"
            )
        if (
            not np.all(np.isfinite(x_amplitude.real))
            or not np.all(np.isfinite(x_amplitude.imag))
        ):
            raise RuntimeError("TDA X amplitudes must be finite")
        transition_densities.append(
            2.0 * c_occ @ x_amplitude @ c_vir.conj().T
        )

    charges = np.asarray(mol.atom_charges(), dtype=float)
    coordinates_bohr = np.asarray(
        mol.atom_coords(unit="Bohr"),
        dtype=float,
    )
    nuclear_charge = float(np.sum(charges))
    if nuclear_charge <= 0.0 or not math.isfinite(nuclear_charge):
        raise RuntimeError("molecular nuclear charge must be positive")
    charge_center_bohr = np.sum(
        charges[:, None] * coordinates_bohr,
        axis=0,
    ) / nuclear_charge
    with mol.with_common_orig(charge_center_bohr):
        dipole_ao = np.asarray(
            mol.intor_symmetric("int1e_r", comp=3),
            dtype=complex,
        )
    dipole_reconstruction_deltas = []
    for index, transition_density in enumerate(transition_densities):
        reconstructed = np.einsum(
            "xpq,pq->x",
            dipole_ao,
            transition_density,
            optimize=True,
        )
        reference = np.asarray(
            transition_dipoles[index],
            dtype=complex,
        )
        delta = float(np.max(np.abs(reconstructed - reference)))
        dipole_reconstruction_deltas.append(delta)
        if not np.allclose(
            reconstructed,
            reference,
            rtol=1.0e-10,
            atol=1.0e-10,
        ):
            raise RuntimeError(
                "AO transition density does not reproduce the PySCF "
                "length-gauge transition dipole"
            )

    ao_overlap = np.asarray(
        mol.intor_symmetric("int1e_ovlp"),
        dtype=float,
    )
    if ao_overlap.shape != (mol.nao_nr(), mol.nao_nr()):
        raise RuntimeError("unexpected AO-overlap shape")
    if not np.all(np.isfinite(ao_overlap)):
        raise RuntimeError("AO overlap must be finite")

    return MolecularTDARuntime(
        result=result,
        molecule=mol,
        ao_overlap=ao_overlap,
        transition_densities_ao=tuple(transition_densities),
        transition_dipole_reconstruction_deltas=tuple(
            dipole_reconstruction_deltas
        ),
    )


def run_closed_shell_tda(
    **kwargs,
) -> MolecularTDAResult:
    """Compatibility wrapper returning the public serializable result only."""
    return run_closed_shell_tda_runtime(**kwargs).result


def cross_geometry_ao_overlap(
    left_runtime: MolecularTDARuntime,
    right_runtime: MolecularTDARuntime,
) -> np.ndarray:
    """Return the AO overlap mapping the right geometry basis into the left."""
    try:
        from pyscf import gto
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "cross-geometry overlap requires the OES q1 extra (PySCF)"
        ) from exc
    overlap = np.asarray(
        gto.intor_cross(
            "int1e_ovlp",
            left_runtime.molecule,
            right_runtime.molecule,
        ),
        dtype=float,
    )
    expected = (
        left_runtime.ao_overlap.shape[0],
        right_runtime.ao_overlap.shape[0],
    )
    if overlap.shape != expected:
        raise RuntimeError(
            f"unexpected cross-geometry AO-overlap shape {overlap.shape}"
        )
    if not np.all(np.isfinite(overlap)):
        raise RuntimeError("cross-geometry AO overlap must be finite")
    return overlap


def run_orientation_tda_scan(
    *,
    species_label: str,
    atoms: Sequence[str],
    samples: Iterable[tuple[float, Sequence[Sequence[float]]]],
    basis_name: str = "cc-pVDZ",
    charge: int = 0,
    spin: int = 0,
    nstates: int = 6,
) -> OrientationTDAScan:
    """Run independent raw-root TDA calculations over an orientation series.

    Root indices are deliberately not treated as state identity across
    geometries. Use OES subspace/state tracking before differentiating or
    following a named transition through the scan.
    """
    rows = tuple(samples)
    if not rows:
        raise MolecularTDAError("orientation TDA scan must not be empty")
    seen: set[float] = set()
    out = []
    for angle_raw, coordinates in rows:
        angle = _finite(angle_raw, "angle_rad")
        if angle in seen:
            raise MolecularTDAError("orientation angles must be unique")
        seen.add(angle)
        out.append(
            OrientationTDAPoint(
                angle_rad=angle,
                result=run_closed_shell_tda(
                    species_label=species_label,
                    atoms=atoms,
                    coordinates_angstrom=coordinates,
                    basis_name=basis_name,
                    charge=charge,
                    spin=spin,
                    nstates=nstates,
                ),
            )
        )
    return OrientationTDAScan(points=tuple(out))


__all__ = [
    "HARTREE_TO_EV",
    "MolecularTDAError",
    "MolecularTDAState",
    "MolecularTDAResult",
    "MolecularTDARuntime",
    "OrientationTDAPoint",
    "OrientationTDAScan",
    "run_closed_shell_tda_runtime",
    "run_closed_shell_tda",
    "cross_geometry_ao_overlap",
    "run_orientation_tda_scan",
]
