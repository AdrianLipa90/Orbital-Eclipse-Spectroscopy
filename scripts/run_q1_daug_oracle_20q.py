#!/usr/bin/env python3
"""Evaluate the nonpredictive d-aug symmetry-aware 20Q capacity bound."""

import json
from pathlib import Path

import numpy as np

from oes.quantum.active_space import transform_active_integrals
from oes.quantum.active_space_oracle_daug import build_helium_daug_oracle_20q
from oes.quantum.diffuse_basis import geometric_multi_augment
from oes.quantum.fermions import build_sector_hamiltonian, transition_one_rdm
from oes.quantum.helium_q1 import classify_states, pyscf_fci_reference, spatial_transition_rdm, spin_squared_matrix

ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]
REFERENCE_F = 0.2762


def main():
    d_aug = {"He": geometric_multi_augment("He", "aug-cc-pVQZ", extra_layers=1)}
    mol, mf, C_active, receipt = build_helium_daug_oracle_20q(
        source_basis=d_aug,
        source_label="d-aug-cc-pVQZ/geometric",
    )
    h1, eri, dip = transform_active_integrals(mol, mf, C_active)
    H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    evals, evecs = np.linalg.eigh(H)
    e0 = float(evals[0])
    e_ref = pyscf_fci_reference(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    if abs(e0 - e_ref) >= 1e-9:
        raise RuntimeError(f"d-aug oracle 20Q OES/PySCF FCI mismatch: {e0-e_ref} Ha")

    s2mat = spin_squared_matrix(10, basis)
    states = classify_states(evals, evecs, s2mat, limit=min(180, len(evals)))
    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    if not triplets or not singlets:
        raise RuntimeError("d-aug oracle 20Q did not resolve both spin sectors")
    triplet = min(triplets, key=lambda s: s.excitation_ev)

    ground = evecs[:, 0]
    rows = []
    for state in singlets:
        t_spin = transition_one_rdm(evecs[:, state.index], ground, basis, 20)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip[k] * t_space) for k in range(3)], dtype=complex)
        mu2 = float(np.sum(np.abs(mu) ** 2))
        delta_h = float(evals[state.index] - e0)
        f = (2.0 / 3.0) * delta_h * mu2
        rows.append({"state": state, "f": f, "mu": float(np.sqrt(mu2))})

    dark = [row for row in rows if row["f"] < 1e-6]
    bright = [row for row in rows if row["f"] > 1e-5]
    if not dark or not bright:
        raise RuntimeError("d-aug oracle 20Q lost dark or bright singlet class")
    dark_row = min(dark, key=lambda row: row["state"].excitation_ev)
    first_bright_e = min(row["state"].excitation_ev for row in bright)
    manifold = [row for row in bright if abs(row["state"].excitation_ev - first_bright_e) < 1e-4]
    manifold_f = float(sum(row["f"] for row in manifold))

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "NONPREDICTIVE_CAPACITY_UPPER_BOUND",
        "protocol": receipt.protocol,
        "compression": receipt.as_dict(),
        "n_spin_orbitals": 20,
        "fixed_particle_dimension": len(basis),
        "active_oes_fci_hartree": e0,
        "active_pyscf_fci_hartree": e_ref,
        "active_fci_delta_hartree": e0 - e_ref,
        "active_ground_loss_hartree": e0 - receipt.source_ground_hartree,
        "triplet": {
            "predicted_eV": triplet.excitation_ev,
            "source_eV": receipt.source_triplet_excitation_ev,
            "source_residual_eV": triplet.excitation_ev - receipt.source_triplet_excitation_ev,
            "nist_target_eV": TARGETS["1s2s_3S1"],
            "nist_residual_eV": triplet.excitation_ev - TARGETS["1s2s_3S1"],
        },
        "dark_singlet": {
            "predicted_eV": dark_row["state"].excitation_ev,
            "source_eV": receipt.source_dark_singlet_excitation_ev,
            "source_residual_eV": dark_row["state"].excitation_ev - receipt.source_dark_singlet_excitation_ev,
            "nist_target_eV": TARGETS["1s2s_1S0"],
            "nist_residual_eV": dark_row["state"].excitation_ev - TARGETS["1s2s_1S0"],
            "oscillator_strength": dark_row["f"],
        },
        "first_bright_manifold": {
            "predicted_eV": first_bright_e,
            "source_eV": receipt.source_bright_manifold_excitation_ev,
            "source_residual_eV": first_bright_e - receipt.source_bright_manifold_excitation_ev,
            "nist_target_eV": TARGETS["1s2p_1P1"],
            "nist_residual_eV": first_bright_e - TARGETS["1s2p_1P1"],
            "degeneracy_count": len(manifold),
            "state_indices": [row["state"].index for row in manifold],
            "oscillator_strength_sum": manifold_f,
            "source_f_sum": receipt.source_bright_manifold_oscillator_strength_sum,
            "source_f_retention": manifold_f / receipt.source_bright_manifold_oscillator_strength_sum,
            "reference_f": REFERENCE_F,
            "reference_f_residual": manifold_f - REFERENCE_F,
            "component_f": [row["f"] for row in manifold],
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
