#!/usr/bin/env python3
"""Evaluate TAAS-v1 at the unchanged 10-spatial / 20Q budget."""

import json
from pathlib import Path

import numpy as np

from oes.quantum.active_space import build_helium_taas_v1, transform_active_integrals
from oes.quantum.fermions import build_sector_hamiltonian, transition_one_rdm
from oes.quantum.helium_q1 import (
    HARTREE_TO_EV,
    classify_states,
    pyscf_fci_reference,
    spatial_transition_rdm,
    spin_squared_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]


def main():
    mol, mf, C_active, compression = build_helium_taas_v1()
    h1, eri, dip = transform_active_integrals(mol, mf, C_active)
    H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    evals, evecs = np.linalg.eigh(H)
    e0 = float(evals[0])
    e_ref = pyscf_fci_reference(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    if abs(e0 - e_ref) >= 1e-9:
        raise RuntimeError(f"TAAS active-space OES/PySCF FCI mismatch: {e0-e_ref} Ha")

    s2mat = spin_squared_matrix(10, basis)
    states = classify_states(evals, evecs, s2mat, limit=min(100, len(evals)))
    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    if not triplets or not singlets:
        raise RuntimeError("TAAS-v1 did not resolve both triplet and singlet excited sectors")

    triplet = triplets[0]
    singlet = singlets[0]
    triplet_ev = triplet.excitation_ev
    singlet_ev = singlet.excitation_ev

    ground = evecs[:, 0]
    bright_rows = []
    for state in singlets:
        t_spin = transition_one_rdm(evecs[:, state.index], ground, basis, 20)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip[k] * t_space) for k in range(3)], dtype=complex)
        mu2 = float(np.sum(np.abs(mu) ** 2))
        delta_h = float(evals[state.index] - e0)
        f = (2.0 / 3.0) * delta_h * mu2
        if f > 1e-8:
            bright_rows.append((state, f, float(np.sqrt(mu2))))
    if not bright_rows:
        raise RuntimeError("TAAS-v1 found no E1-bright singlet")
    bright, bright_f, bright_mu = min(bright_rows, key=lambda item: item[0].excitation_ev)

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "compression": compression.as_dict(),
        "n_spin_orbitals": 20,
        "fixed_particle_dimension": len(basis),
        "active_oes_fci_hartree": e0,
        "active_pyscf_fci_hartree": e_ref,
        "active_fci_delta_hartree": e0 - e_ref,
        "source_full_basis_ground_fci_hartree": compression.ground_fci_hartree,
        "active_ground_loss_hartree": e0 - compression.ground_fci_hartree,
        "calibration": {
            "1s2s_3S1_predicted_eV": triplet_ev,
            "1s2s_3S1_target_eV": TARGETS["1s2s_3S1"],
            "1s2s_3S1_residual_eV": triplet_ev - TARGETS["1s2s_3S1"],
            "1s2s_1S0_predicted_eV": singlet_ev,
            "1s2s_1S0_target_eV": TARGETS["1s2s_1S0"],
            "1s2s_1S0_residual_eV": singlet_ev - TARGETS["1s2s_1S0"],
        },
        "held_out_bright": {
            "state_index": bright.index,
            "predicted_eV": bright.excitation_ev,
            "target_eV": TARGETS["1s2p_1P1"],
            "residual_eV": bright.excitation_ev - TARGETS["1s2p_1P1"],
            "abs_residual_eV": abs(bright.excitation_ev - TARGETS["1s2p_1P1"]),
            "oscillator_strength": bright_f,
            "transition_dipole_norm_au": bright_mu,
            "s2": bright.s2,
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
