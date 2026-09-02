#!/usr/bin/env python3
"""Converge the fixed 20Q helium core with a state-balanced selected Q bath.

The P-space is always the same predictive d-aug s4+p6 20Q active space. External
determinants are ranked only from Hamiltonian coupling to P-space states. The
selected bath is classical and retains full Q-Q couplings through exact
subspace diagonalization. No NIST value enters selection or Hamiltonian build.
"""

import json
from pathlib import Path

import numpy as np

from oes.quantum.active_space_blocks import build_helium_s4_p6_20q
from oes.quantum.diffuse_basis import geometric_multi_augment
from oes.quantum.external_dressing import build_external_coupling_space, complete_mo_basis_preserving_active
from oes.quantum.fermions import build_sector_hamiltonian, transition_one_rdm
from oes.quantum.helium_q1 import classify_states, spatial_transition_rdm, spin_squared_matrix
from oes.quantum.selected_ci import (
    build_two_electron_subspace_hamiltonian,
    grouped_importance_order,
    grouped_prefix_for_target,
    state_balanced_external_importance,
)

ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]
HARTREE_TO_EV = 27.211_386_245_981
SOURCE = {
    "triplet": 19.80281611421936,
    "dark": 20.61327300402567,
    "bright": 21.29949008913645,
}
TARGET_BATH_COUNTS = (32, 64, 128, 256, 512)


def identify_active_classes(evals, evecs, basis, dip_active):
    e0 = float(evals[0])
    s2mat = spin_squared_matrix(10, basis)
    states = classify_states(evals, evecs, s2mat, limit=min(180, len(evals)))
    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    if not triplets or not singlets:
        raise RuntimeError("selected-Q diagnostic failed to resolve active spin sectors")
    triplet = min(triplets, key=lambda s: s.excitation_ev)

    ground = evecs[:, 0]
    rows = []
    for state in singlets:
        t_spin = transition_one_rdm(evecs[:, state.index], ground, basis, 20)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip_active[k] * t_space) for k in range(3)], dtype=complex)
        mu2 = float(np.sum(np.abs(mu) ** 2))
        de = float(evals[state.index] - e0)
        f = (2.0 / 3.0) * de * mu2
        rows.append({"state": state, "f": f})
    dark_rows = [row for row in rows if row["f"] < 1e-6]
    bright_rows = [row for row in rows if row["f"] > 1e-5]
    if not dark_rows or not bright_rows:
        raise RuntimeError("selected-Q diagnostic lost dark/bright active classes")
    dark = min(dark_rows, key=lambda row: row["state"].excitation_ev)["state"]
    first_bright_e = min(row["state"].excitation_ev for row in bright_rows)
    bright = [row["state"] for row in bright_rows if abs(row["state"].excitation_ev - first_bright_e) < 1e-4]
    if len(bright) != 3:
        raise RuntimeError(f"expected 3 active bright components, got {len(bright)}")
    return triplet, dark, bright


def match_state(target, selected_vectors, p_dim):
    overlaps = np.abs(np.conjugate(target) @ selected_vectors[:p_dim, :]) ** 2
    idx = int(np.argmax(overlaps))
    return idx, float(overlaps[idx])


def main():
    from pyscf import ao2mo

    d_aug = {"He": geometric_multi_augment("He", "aug-cc-pVQZ", extra_layers=1)}
    mol, mf, C_active, receipt = build_helium_s4_p6_20q(
        source_basis=d_aug,
        source_label="d-aug-cc-pVQZ/geometric",
    )
    C_full = complete_mo_basis_preserving_active(mf, C_active)
    nfull = C_full.shape[1]
    nfull_spin = 2 * nfull
    h1_full = C_full.T @ mf.get_hcore() @ C_full
    eri_full = ao2mo.kernel(mol, C_full, compact=False).reshape((nfull,) * 4)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip_active = np.stack([C_active.T @ dip_ao[k] @ C_active for k in range(3)])
    dip_full = np.stack([C_full.T @ dip_ao[k] @ C_full for k in range(3)])

    h1_active = h1_full[:10, :10]
    eri_active = eri_full[:10, :10, :10, :10]
    H_active, p_basis = build_sector_hamiltonian(
        h1_active,
        eri_active,
        n_electrons=2,
        ecore=float(mol.energy_nuc()),
    )
    evals, evecs = np.linalg.eigh(H_active)
    p_dim = len(p_basis)
    if p_dim != 190:
        raise RuntimeError(f"unexpected P dimension {p_dim}")

    triplet, dark, bright = identify_active_classes(evals, evecs, p_basis, dip_active)
    bright_indices = [s.index for s in bright]
    bright_states = np.column_stack([evecs[:, i] for i in bright_indices])
    bright_energy = float(np.mean([evals[i] for i in bright_indices]))

    external = build_external_coupling_space(
        h1_full,
        eri_full,
        active_basis=p_basis,
        n_active_spatial=10,
        n_electrons=2,
        ecore=float(mol.energy_nuc()),
    )
    importance = state_balanced_external_importance(
        external,
        class_states={
            "ground": evecs[:, 0],
            "triplet": evecs[:, triplet.index],
            "dark": evecs[:, dark.index],
            "bright": bright_states,
        },
        class_energies_hartree={
            "ground": float(evals[0]),
            "triplet": float(evals[triplet.index]),
            "dark": float(evals[dark.index]),
            "bright": bright_energy,
        },
    )
    scores = np.asarray(importance["scores"], dtype=float)
    groups = grouped_importance_order(scores, relative_tie_tolerance=1e-7)
    prefixes = {target: grouped_prefix_for_target(groups, target) for target in TARGET_BATH_COUNTS}
    largest = prefixes[max(TARGET_BATH_COUNTS)]
    for target in TARGET_BATH_COUNTS:
        prefix = prefixes[target]
        if prefix != largest[: len(prefix)]:
            raise RuntimeError("importance group prefixes are not nested")

    selected_external_dets = [external.external_basis[i] for i in largest]
    max_basis = list(p_basis) + selected_external_dets
    H_max = build_two_electron_subspace_hamiltonian(
        h1_full,
        eri_full,
        max_basis,
        ecore=float(mol.energy_nuc()),
    )
    p_block_error = float(np.max(np.abs(H_max[:p_dim, :p_dim] - H_active)))
    if p_block_error > 1e-10:
        raise RuntimeError(f"selected-Q P block mismatch: {p_block_error} Ha")

    active_exc = {
        "triplet": float((evals[triplet.index] - evals[0]) * HARTREE_TO_EV),
        "dark": float((evals[dark.index] - evals[0]) * HARTREE_TO_EV),
        "bright": float((bright_energy - evals[0]) * HARTREE_TO_EV),
    }

    sweep = []
    per_class_scores = importance["per_class_scores"]
    for target in TARGET_BATH_COUNTS:
        selected_indices = prefixes[target]
        n_q = len(selected_indices)
        dim = p_dim + n_q
        basis = max_basis[:dim]
        H = H_max[:dim, :dim]
        vals, vecs = np.linalg.eigh(H)

        g_idx, g_overlap = match_state(evecs[:, 0], vecs, p_dim)
        t_idx, t_overlap = match_state(evecs[:, triplet.index], vecs, p_dim)
        d_idx, d_overlap = match_state(evecs[:, dark.index], vecs, p_dim)
        if len({g_idx, t_idx, d_idx}) != 3:
            raise RuntimeError("selected-Q single-state matching collapsed distinct P targets")

        bright_overlap = np.sum(np.abs(np.conjugate(bright_states).T @ vecs[:p_dim, :]) ** 2, axis=0)
        b_candidates = [int(i) for i in np.argsort(bright_overlap)[::-1] if int(i) not in {g_idx, t_idx, d_idx}]
        b_idx = b_candidates[:3]
        if len(b_idx) != 3:
            raise RuntimeError("selected-Q failed to match bright three-state subspace")
        b_idx = sorted(b_idx, key=lambda i: vals[i])

        ground_energy = float(vals[g_idx])
        triplet_ev = float((vals[t_idx] - ground_energy) * HARTREE_TO_EV)
        dark_ev = float((vals[d_idx] - ground_energy) * HARTREE_TO_EV)
        bright_components_ev = np.array([(vals[i] - ground_energy) * HARTREE_TO_EV for i in b_idx])
        bright_ev = float(np.mean(bright_components_ev))
        bright_spread = float(np.max(bright_components_ev) - np.min(bright_components_ev))

        # Observable transition strength in the selected P+Q subspace.
        f_components = []
        for idx in b_idx:
            t_spin = transition_one_rdm(vecs[:, idx], vecs[:, g_idx], basis, nfull_spin)
            t_space = spatial_transition_rdm(t_spin)
            mu = np.array([np.sum(dip_full[k] * t_space) for k in range(3)], dtype=complex)
            mu2 = float(np.sum(np.abs(mu) ** 2))
            de = float(vals[idx] - ground_energy)
            f_components.append((2.0 / 3.0) * de * mu2)

        classes = {
            "triplet": {"eV": triplet_ev, "source": SOURCE["triplet"], "nist": TARGETS["1s2s_3S1"]},
            "dark": {"eV": dark_ev, "source": SOURCE["dark"], "nist": TARGETS["1s2s_1S0"]},
            "bright": {"eV": bright_ev, "source": SOURCE["bright"], "nist": TARGETS["1s2p_1P1"]},
        }
        source_res = np.array([classes[k]["eV"] - classes[k]["source"] for k in ("triplet", "dark", "bright")])
        nist_res = np.array([classes[k]["eV"] - classes[k]["nist"] for k in ("triplet", "dark", "bright")])

        selected_array = np.asarray(selected_indices, dtype=int)
        captured = {
            name: float(np.sum(np.asarray(per_class_scores[name])[selected_array]))
            for name in ("ground", "triplet", "dark", "bright")
        }
        p_weights = {
            "ground": float(np.sum(np.abs(vecs[:p_dim, g_idx]) ** 2)),
            "triplet": float(np.sum(np.abs(vecs[:p_dim, t_idx]) ** 2)),
            "dark": float(np.sum(np.abs(vecs[:p_dim, d_idx]) ** 2)),
            "bright_mean": float(np.mean([np.sum(np.abs(vecs[:p_dim, i]) ** 2) for i in b_idx])),
        }
        sweep.append(
            {
                "target_external_count": target,
                "selected_external_count": n_q,
                "selected_dimension": dim,
                "importance_captured": captured,
                "p_weights": p_weights,
                "match_overlap": {
                    "ground": g_overlap,
                    "triplet": t_overlap,
                    "dark": d_overlap,
                    "bright_components": [float(bright_overlap[i]) for i in b_idx],
                },
                "energies_eV": {"triplet": triplet_ev, "dark": dark_ev, "bright": bright_ev},
                "source_residuals_eV": {k: float(classes[k]["eV"] - classes[k]["source"]) for k in classes},
                "nist_residuals_eV": {k: float(classes[k]["eV"] - classes[k]["nist"]) for k in classes},
                "source_rms_eV": float(np.sqrt(np.mean(source_res**2))),
                "source_centered_rms_eV": float(np.sqrt(np.mean((source_res - np.mean(source_res)) ** 2))),
                "nist_rms_eV": float(np.sqrt(np.mean(nist_res**2))),
                "bright_manifold_spread_eV": bright_spread,
                "bright_f_components": [float(x) for x in f_components],
                "bright_f_sum": float(np.sum(f_components)),
            }
        )

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "STATE_BALANCED_SELECTED_Q_CI_DIAGNOSTIC_NO_EXPERIMENTAL_SELECTION_INPUT",
        "active_protocol": receipt.protocol,
        "n_active_spin_orbitals": 20,
        "active_dimension": p_dim,
        "source_spin_orbitals": nfull_spin,
        "full_external_determinants": len(external.external_basis),
        "p_block_max_error_hartree": p_block_error,
        "active_baseline_eV": active_exc,
        "importance_group_sizes_first20": [len(g) for g in groups[:20]],
        "sweep": sweep,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
