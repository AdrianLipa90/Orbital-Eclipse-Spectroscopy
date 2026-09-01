#!/usr/bin/env python3
"""Benchmark ROAS-v1 operator-class response selection at exactly 20Q."""

import json
from pathlib import Path

import numpy as np

from oes.quantum.active_space import transform_active_integrals
from oes.quantum.active_space_roas import build_helium_roas_v1_20q
from oes.quantum.diffuse_basis import geometric_multi_augment
from oes.quantum.fermions import build_sector_hamiltonian, transition_one_rdm
from oes.quantum.helium_q1 import classify_states, pyscf_fci_reference, spatial_transition_rdm, spin_squared_matrix

ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]
REFERENCE_F = 0.2762
SOURCE = {
    "triplet": 19.80281611421936,
    "dark": 20.61327300402567,
    "bright": 21.29949008913645,
    "bright_f": 0.32933276076136725,
}


def main():
    d_aug = {"He": geometric_multi_augment("He", "aug-cc-pVQZ", extra_layers=1)}
    mol, mf, C_active, receipt = build_helium_roas_v1_20q(
        source_basis=d_aug,
        source_label="d-aug-cc-pVQZ/geometric",
    )
    h1, eri, dip = transform_active_integrals(mol, mf, C_active)
    H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    evals, evecs = np.linalg.eigh(H)
    e0 = float(evals[0])
    e_ref = pyscf_fci_reference(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    if abs(e0 - e_ref) >= 1e-9:
        raise RuntimeError(f"ROAS-v1 20Q OES/PySCF FCI mismatch: {e0-e_ref} Ha")

    s2mat = spin_squared_matrix(10, basis)
    states = classify_states(evals, evecs, s2mat, limit=min(180, len(evals)))
    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    if not triplets or not singlets:
        raise RuntimeError("ROAS-v1 did not resolve both spin sectors")
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
        raise RuntimeError("ROAS-v1 lost dark or bright singlet class")
    dark_row = min(dark, key=lambda row: row["state"].excitation_ev)
    first_bright_e = min(row["state"].excitation_ev for row in bright)
    manifold = [row for row in bright if abs(row["state"].excitation_ev - first_bright_e) < 1e-4]
    manifold_f = float(sum(row["f"] for row in manifold))

    source_residuals = {
        "triplet": triplet.excitation_ev - SOURCE["triplet"],
        "dark": dark_row["state"].excitation_ev - SOURCE["dark"],
        "bright": first_bright_e - SOURCE["bright"],
    }
    source_mean = float(np.mean(list(source_residuals.values())))
    source_centered_rms = float(np.sqrt(np.mean([(x - source_mean) ** 2 for x in source_residuals.values()])))

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "PREDICTIVE_SELECTOR_NO_EXCITED_TARGET_INPUT",
        "protocol": receipt.protocol,
        "compression": receipt.as_dict(),
        "n_spin_orbitals": 20,
        "fixed_particle_dimension": len(basis),
        "active_oes_fci_hartree": e0,
        "active_pyscf_fci_hartree": e_ref,
        "active_fci_delta_hartree": e0 - e_ref,
        "active_ground_loss_hartree": e0 - receipt.ground_fci_hartree,
        "triplet": {
            "predicted_eV": triplet.excitation_ev,
            "source_residual_eV": source_residuals["triplet"],
            "nist_residual_eV": triplet.excitation_ev - TARGETS["1s2s_3S1"],
        },
        "dark_singlet": {
            "predicted_eV": dark_row["state"].excitation_ev,
            "source_residual_eV": source_residuals["dark"],
            "nist_residual_eV": dark_row["state"].excitation_ev - TARGETS["1s2s_1S0"],
            "oscillator_strength": dark_row["f"],
        },
        "first_bright_manifold": {
            "predicted_eV": first_bright_e,
            "source_residual_eV": source_residuals["bright"],
            "nist_residual_eV": first_bright_e - TARGETS["1s2p_1P1"],
            "degeneracy_count": len(manifold),
            "state_indices": [row["state"].index for row in manifold],
            "oscillator_strength_sum": manifold_f,
            "source_f_retention": manifold_f / SOURCE["bright_f"],
            "reference_f": REFERENCE_F,
            "reference_f_residual": manifold_f - REFERENCE_F,
            "component_f": [row["f"] for row in manifold],
        },
        "source_energy_residual_mean_eV": source_mean,
        "source_energy_centered_rms_eV": source_centered_rms,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
