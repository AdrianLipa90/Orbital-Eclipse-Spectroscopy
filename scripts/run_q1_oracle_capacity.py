#!/usr/bin/env python3
"""Evaluate the non-predictive oracle state-averaged 20Q capacity bound."""

import json
from pathlib import Path

import numpy as np

from oes.quantum.active_space import transform_active_integrals
from oes.quantum.active_space_oracle import build_helium_oracle_capacity_space
from oes.quantum.fermions import build_sector_hamiltonian, transition_one_rdm
from oes.quantum.helium_q1 import classify_states, pyscf_fci_reference, spatial_transition_rdm, spin_squared_matrix


ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]


def main():
    mol, mf, C_active, compression = build_helium_oracle_capacity_space()
    h1, eri, dip = transform_active_integrals(mol, mf, C_active)
    H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    evals, evecs = np.linalg.eigh(H)
    e0 = float(evals[0])
    e_ref = pyscf_fci_reference(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    if abs(e0 - e_ref) >= 1e-9:
        raise RuntimeError(f"oracle active OES/PySCF FCI mismatch: {e0-e_ref} Ha")

    s2mat = spin_squared_matrix(10, basis)
    states = classify_states(evals, evecs, s2mat, limit=min(100, len(evals)))
    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    triplet = min(triplets, key=lambda s: s.excitation_ev)
    dark_singlet_candidates = []
    bright_rows = []
    ground = evecs[:, 0]
    for state in singlets:
        t_spin = transition_one_rdm(evecs[:, state.index], ground, basis, 20)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip[k] * t_space) for k in range(3)], dtype=complex)
        mu2 = float(np.sum(np.abs(mu) ** 2))
        delta_h = float(evals[state.index] - e0)
        f = (2.0 / 3.0) * delta_h * mu2
        if f < 1e-6:
            dark_singlet_candidates.append((state, f))
        if f > 1e-4:
            bright_rows.append((state, f, float(np.sqrt(mu2))))
    if not dark_singlet_candidates or not bright_rows:
        raise RuntimeError("oracle active space lost dark or bright singlet class")
    dark, dark_f = min(dark_singlet_candidates, key=lambda item: item[0].excitation_ev)
    bright, bright_f, bright_mu = min(bright_rows, key=lambda item: item[0].excitation_ev)

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "NONPREDICTIVE_CAPACITY_UPPER_BOUND",
        "compression": compression.as_dict(),
        "n_spin_orbitals": 20,
        "fixed_particle_dimension": len(basis),
        "active_ground_hartree": e0,
        "active_pyscf_fci_hartree": e_ref,
        "active_fci_delta_hartree": e0 - e_ref,
        "active_ground_loss_hartree": e0 - compression.source_ground_hartree,
        "triplet": {
            "predicted_eV": triplet.excitation_ev,
            "target_eV": TARGETS["1s2s_3S1"],
            "residual_eV": triplet.excitation_ev - TARGETS["1s2s_3S1"],
        },
        "dark_singlet": {
            "predicted_eV": dark.excitation_ev,
            "target_eV": TARGETS["1s2s_1S0"],
            "residual_eV": dark.excitation_ev - TARGETS["1s2s_1S0"],
            "oscillator_strength": dark_f,
        },
        "bright_singlet": {
            "predicted_eV": bright.excitation_ev,
            "target_eV": TARGETS["1s2p_1P1"],
            "residual_eV": bright.excitation_ev - TARGETS["1s2p_1P1"],
            "abs_residual_eV": abs(bright.excitation_ev - TARGETS["1s2p_1P1"]),
            "oscillator_strength": bright_f,
            "transition_dipole_norm_au": bright_mu,
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
