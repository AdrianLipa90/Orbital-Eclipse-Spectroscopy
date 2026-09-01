"""Held-out helium spectroscopy gate for the fixed 20Q active-space budget."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, List, Tuple

import numpy as np

from .fermions import build_sector_hamiltonian, transition_one_rdm
from .helium_q1 import HARTREE_TO_EV, classify_states, spatial_transition_rdm, spin_squared_matrix


@dataclass(frozen=True)
class BrightSingletPrediction:
    backend: str
    basis_name: str
    n_spatial_orbitals: int
    n_spin_orbitals: int
    state_index: int
    excitation_ev: float
    energy_hartree: float
    oscillator_strength: float
    transition_dipole_norm_au: float
    s2: float

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def _prepare_with_dipoles(basis_name: str, n_spatial: int):
    try:
        from pyscf import ao2mo, gto, scf
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("helium spectroscopy requires the OES q1 extra (PySCF)") from exc

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
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("helium RHF did not converge")
    if mf.mo_coeff.shape[1] < n_spatial:
        raise ValueError(f"basis {basis_name} provides fewer than {n_spatial} MOs")

    coeff = mf.mo_coeff[:, :n_spatial]
    h1 = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.kernel(mol, coeff, compact=False).reshape((n_spatial,) * 4)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip_mo = np.stack([coeff.T @ dip_ao[k] @ coeff for k in range(3)])
    return mol, h1, eri, dip_mo


def first_bright_singlet_prediction(
    basis_name: str = "aug-cc-pVQZ",
    n_spatial: int = 10,
    bright_threshold: float = 1e-8,
) -> Tuple[BrightSingletPrediction, List[Dict[str, float]]]:
    """Predict the first E1-bright singlet from the ground state.

    The active-space basis was selected using the 2s singlet/triplet sweep.  The
    1s2p 1P experimental energy is not used here and remains held out.
    """
    if n_spatial != 10:
        raise ValueError("canonical Q1 spectroscopy gate is fixed at 10 spatial orbitals / 20 qubits")

    mol, h1, eri, dip_mo = _prepare_with_dipoles(basis_name, n_spatial)
    H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    evals, evecs = np.linalg.eigh(H)
    s2mat = spin_squared_matrix(n_spatial, basis)
    states = classify_states(evals, evecs, s2mat, limit=min(80, len(evals)))
    ground = evecs[:, 0]
    e0 = float(evals[0])
    n_spin = 2 * n_spatial

    candidates: List[Dict[str, float]] = []
    for state in states[1:]:
        if abs(state.s2) > 1e-6:
            continue
        t_spin = transition_one_rdm(evecs[:, state.index], ground, basis, n_spin)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip_mo[k] * t_space) for k in range(3)], dtype=complex)
        mu2 = float(np.sum(np.abs(mu) ** 2))
        delta_e_h = float(evals[state.index] - e0)
        oscillator_strength = (2.0 / 3.0) * delta_e_h * mu2
        candidates.append(
            {
                "index": int(state.index),
                "excitation_ev": delta_e_h * HARTREE_TO_EV,
                "energy_hartree": float(evals[state.index]),
                "oscillator_strength": float(oscillator_strength),
                "transition_dipole_norm_au": float(np.sqrt(mu2)),
                "s2": float(state.s2),
            }
        )

    bright = [row for row in candidates if row["oscillator_strength"] > bright_threshold]
    if not bright:
        raise RuntimeError("no E1-bright singlet found in classified Q1 state window")
    first = min(bright, key=lambda row: row["excitation_ev"])
    result = BrightSingletPrediction(
        backend="SIMULATED_REFERENCE",
        basis_name=basis_name,
        n_spatial_orbitals=n_spatial,
        n_spin_orbitals=2 * n_spatial,
        state_index=int(first["index"]),
        excitation_ev=float(first["excitation_ev"]),
        energy_hartree=float(first["energy_hartree"]),
        oscillator_strength=float(first["oscillator_strength"]),
        transition_dipole_norm_au=float(first["transition_dipole_norm_au"]),
        s2=float(first["s2"]),
    )
    return result, candidates
