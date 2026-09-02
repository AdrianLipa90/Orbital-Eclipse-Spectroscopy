#!/usr/bin/env python3
"""Test a rotation-covariant external natural-orbital bath around the fixed 20Q core.

Selection uses only the normalized QH response of complete active-state classes:
- ground singlet,
- complete lowest triplet spin manifold,
- lowest dark excited singlet,
- complete first bright singlet p manifold.

No NIST energy enters the bath construction. A deliberate random orthogonal
rotation of the full external one-particle basis provides a hard gauge test:
the selected physical bath projector must remain unchanged.
"""

import json
from pathlib import Path

import numpy as np

from oes.quantum.active_space_blocks import build_helium_s4_p6_20q
from oes.quantum.diffuse_basis import geometric_multi_augment
from oes.quantum.external_dressing import build_external_coupling_space, complete_mo_basis_preserving_active
from oes.quantum.fermions import build_sector_hamiltonian, determinant_basis, transition_one_rdm
from oes.quantum.helium_q1 import classify_states, spatial_transition_rdm, spin_squared_matrix
from oes.quantum.orbital_bath import bath_projector_overlap, select_external_natural_bath

ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]
HARTREE_TO_EV = 27.211_386_245_981
SOURCE = {
    "triplet": 19.80281611421936,
    "dark": 20.61327300402567,
    "bright": 21.29949008913645,
}
TARGET_EXTERNAL_SPATIAL = (4, 7, 10, 13)
MANIFOLD_TOL_EV = 1e-4


def active_classes(evals, evecs, p_basis, dip_active):
    e0 = float(evals[0])
    states = classify_states(evals, evecs, spin_squared_matrix(10, p_basis), limit=min(190, len(evals)))
    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    if not triplets or not singlets:
        raise RuntimeError("could not resolve active singlet/triplet sectors")

    triplet_e = min(s.excitation_ev for s in triplets)
    triplet_manifold = [s for s in triplets if abs(s.excitation_ev - triplet_e) < MANIFOLD_TOL_EV]
    if len(triplet_manifold) != 3:
        raise RuntimeError(f"expected complete three-state lowest triplet manifold, got {len(triplet_manifold)}")

    ground = evecs[:, 0]
    rows = []
    for state in singlets:
        t_spin = transition_one_rdm(evecs[:, state.index], ground, p_basis, 20)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip_active[k] * t_space) for k in range(3)], dtype=complex)
        mu2 = float(np.sum(np.abs(mu) ** 2))
        de = float(evals[state.index] - e0)
        rows.append((state, (2.0 / 3.0) * de * mu2))
    dark_rows = [(s, f) for s, f in rows if f < 1e-6]
    bright_rows = [(s, f) for s, f in rows if f > 1e-5]
    if not dark_rows or not bright_rows:
        raise RuntimeError("could not resolve active dark/bright singlet classes")
    dark = min(dark_rows, key=lambda x: x[0].excitation_ev)[0]
    bright_e = min(s.excitation_ev for s, _ in bright_rows)
    bright_manifold = [s for s, _ in bright_rows if abs(s.excitation_ev - bright_e) < MANIFOLD_TOL_EV]
    if len(bright_manifold) != 3:
        raise RuntimeError(f"expected complete three-state bright manifold, got {len(bright_manifold)}")

    return {
        "ground": evecs[:, 0],
        "triplet": np.column_stack([evecs[:, s.index] for s in triplet_manifold]),
        "dark": evecs[:, dark.index],
        "bright": np.column_stack([evecs[:, s.index] for s in bright_manifold]),
    }, {
        "ground": [0],
        "triplet": [s.index for s in triplet_manifold],
        "dark": [dark.index],
        "bright": [s.index for s in bright_manifold],
    }


def class_overlap_weights(target_states, full_vectors, p_indices):
    target = np.asarray(target_states, dtype=complex)
    if target.ndim == 1:
        target = target[:, None]
    projected = full_vectors[p_indices, :]
    return np.sum(np.abs(target.conj().T @ projected) ** 2, axis=0)


def choose_class_states(target_states, full_vectors, p_indices, count, excluded):
    weights = class_overlap_weights(target_states, full_vectors, p_indices)
    candidates = [int(i) for i in np.argsort(weights)[::-1] if int(i) not in excluded]
    chosen = candidates[:count]
    if len(chosen) != count:
        raise RuntimeError("insufficient matched states for active class")
    return chosen, [float(weights[i]) for i in chosen]


def build_for_source(mol, mf, C_active, C_full, class_states, target_external_spatial):
    from pyscf import ao2mo

    nfull = C_full.shape[1]
    h1_full = C_full.T @ mf.get_hcore() @ C_full
    eri_full = ao2mo.kernel(mol, C_full, compact=False).reshape((nfull,) * 4)
    p_basis = determinant_basis(20, 2)
    external = build_external_coupling_space(
        h1_full,
        eri_full,
        active_basis=p_basis,
        n_active_spatial=10,
        n_electrons=2,
        ecore=float(mol.energy_nuc()),
    )
    bath = select_external_natural_bath(
        external,
        class_states=class_states,
        n_active_spatial=10,
        target_external_spatial=target_external_spatial,
        relative_degeneracy_tolerance=2e-5,
        absolute_degeneracy_tolerance=1e-10,
    )
    C_q = C_full[:, 10:]
    C_bath = C_q @ bath.q_rotation
    return bath, C_bath


def evaluate_combined(mol, mf, C_active, C_bath, class_states):
    from pyscf import ao2mo

    C = np.column_stack([C_active, C_bath])
    nsp = C.shape[1]
    nspin = 2 * nsp
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.kernel(mol, C, compact=False).reshape((nsp,) * 4)
    H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    vals, vecs = np.linalg.eigh(H)
    index = {det: i for i, det in enumerate(basis)}
    p_indices = [index[det] for det in determinant_basis(20, 2)]

    excluded = set()
    g_idx, g_ov = choose_class_states(class_states["ground"], vecs, p_indices, 1, excluded)
    excluded.update(g_idx)
    t_idx, t_ov = choose_class_states(class_states["triplet"], vecs, p_indices, 3, excluded)
    excluded.update(t_idx)
    d_idx, d_ov = choose_class_states(class_states["dark"], vecs, p_indices, 1, excluded)
    excluded.update(d_idx)
    b_idx, b_ov = choose_class_states(class_states["bright"], vecs, p_indices, 3, excluded)

    ground_idx = g_idx[0]
    ground_e = float(vals[ground_idx])
    triplet_components = np.array([(vals[i] - ground_e) * HARTREE_TO_EV for i in t_idx], dtype=float)
    dark_ev = float((vals[d_idx[0]] - ground_e) * HARTREE_TO_EV)
    bright_components = np.array([(vals[i] - ground_e) * HARTREE_TO_EV for i in b_idx], dtype=float)
    triplet_ev = float(np.mean(triplet_components))
    bright_ev = float(np.mean(bright_components))

    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([C.T @ dip_ao[k] @ C for k in range(3)])
    f_components = []
    for idx_b in b_idx:
        t_spin = transition_one_rdm(vecs[:, idx_b], vecs[:, ground_idx], basis, nspin)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip[k] * t_space) for k in range(3)], dtype=complex)
        de = float(vals[idx_b] - ground_e)
        f_components.append((2.0 / 3.0) * de * float(np.sum(np.abs(mu) ** 2)))

    energies = {"triplet": triplet_ev, "dark": dark_ev, "bright": bright_ev}
    source_res = np.array([energies[k] - SOURCE[k] for k in ("triplet", "dark", "bright")])
    nist_targets = {
        "triplet": TARGETS["1s2s_3S1"],
        "dark": TARGETS["1s2s_1S0"],
        "bright": TARGETS["1s2p_1P1"],
    }
    nist_res = np.array([energies[k] - nist_targets[k] for k in ("triplet", "dark", "bright")])
    p_weights = {
        "ground": float(np.sum(np.abs(vecs[p_indices, ground_idx]) ** 2)),
        "triplet_mean": float(np.mean([np.sum(np.abs(vecs[p_indices, i]) ** 2) for i in t_idx])),
        "dark": float(np.sum(np.abs(vecs[p_indices, d_idx[0]]) ** 2)),
        "bright_mean": float(np.mean([np.sum(np.abs(vecs[p_indices, i]) ** 2) for i in b_idx])),
    }
    return {
        "combined_spatial_orbitals": nsp,
        "combined_spin_orbitals": nspin,
        "combined_dimension": len(basis),
        "energies_eV": energies,
        "source_residuals_eV": {k: float(energies[k] - SOURCE[k]) for k in energies},
        "nist_residuals_eV": {k: float(energies[k] - nist_targets[k]) for k in energies},
        "source_rms_eV": float(np.sqrt(np.mean(source_res**2))),
        "source_centered_rms_eV": float(np.sqrt(np.mean((source_res - np.mean(source_res)) ** 2))),
        "nist_rms_eV": float(np.sqrt(np.mean(nist_res**2))),
        "triplet_spread_eV": float(np.max(triplet_components) - np.min(triplet_components)),
        "bright_spread_eV": float(np.max(bright_components) - np.min(bright_components)),
        "bright_f_components": [float(x) for x in f_components],
        "bright_f_sum": float(np.sum(f_components)),
        "class_match_overlaps": {
            "ground": g_ov,
            "triplet": t_ov,
            "dark": d_ov,
            "bright": b_ov,
        },
        "p_weights": p_weights,
    }


def main():
    d_aug = {"He": geometric_multi_augment("He", "aug-cc-pVQZ", extra_layers=1)}
    mol, mf, C_active, receipt = build_helium_s4_p6_20q(
        source_basis=d_aug,
        source_label="d-aug-cc-pVQZ/geometric",
    )
    C_full = complete_mo_basis_preserving_active(mf, C_active)
    S = np.asarray(mf.get_ovlp(), dtype=float)

    # Active classes are defined once in the fixed 20Q core.
    from pyscf import ao2mo
    h1_p = C_active.T @ mf.get_hcore() @ C_active
    eri_p = ao2mo.kernel(mol, C_active, compact=False).reshape((10,) * 4)
    H_p, p_basis = build_sector_hamiltonian(h1_p, eri_p, n_electrons=2, ecore=float(mol.energy_nuc()))
    evals_p, evecs_p = np.linalg.eigh(H_p)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip_p = np.stack([C_active.T @ dip_ao[k] @ C_active for k in range(3)])
    classes, class_indices = active_classes(evals_p, evecs_p, p_basis, dip_p)

    sweep = []
    original_baths = {}
    for target in TARGET_EXTERNAL_SPATIAL:
        bath, C_bath = build_for_source(mol, mf, C_active, C_full, classes, target)
        original_baths[target] = C_bath
        result = evaluate_combined(mol, mf, C_active, C_bath, classes)
        result["target_external_spatial"] = target
        result["selected_external_spatial"] = bath.selected_external_spatial_orbitals
        result["retained_response_occupation"] = bath.retained_normalized_occupation
        result["occupation_group_sizes_first16"] = list(bath.occupation_group_sizes[:16])
        result["occupation_group_values_first16"] = list(bath.occupation_group_values[:16])
        result["selected_group_indices"] = list(bath.selected_group_indices)
        result["class_external_traces"] = bath.class_external_traces
        sweep.append(result)

    # Hard gauge test at target=10: arbitrarily rotate every external spatial
    # orbital, rebuild QH and verify the selected physical natural-orbital
    # subspace is unchanged.
    rng = np.random.default_rng(20260901)
    nq = C_full.shape[1] - 10
    random = rng.normal(size=(nq, nq))
    R, _ = np.linalg.qr(random)
    C_full_rot = np.column_stack([C_active, C_full[:, 10:] @ R])
    bath_rot, C_bath_rot = build_for_source(mol, mf, C_active, C_full_rot, classes, 10)
    C_bath_ref = original_baths[10]
    if C_bath_ref.shape[1] != C_bath_rot.shape[1]:
        raise RuntimeError(
            f"rotation gauge changed selected bath dimension: {C_bath_ref.shape[1]} vs {C_bath_rot.shape[1]}"
        )
    gauge = bath_projector_overlap(C_bath_ref, C_bath_rot, S)
    gauge["reference_selected_external_spatial"] = C_bath_ref.shape[1]
    gauge["rotated_selected_external_spatial"] = C_bath_rot.shape[1]
    gauge["rotated_retained_response_occupation"] = bath_rot.retained_normalized_occupation

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "ROTATION_COVARIANT_COMPLETE_CLASS_EXTERNAL_NATURAL_ORBITAL_BATH",
        "active_protocol": receipt.protocol,
        "active_spin_orbitals": 20,
        "active_dimension": 190,
        "source_spatial_orbitals": C_full.shape[1],
        "source_spin_orbitals": 2 * C_full.shape[1],
        "active_class_indices": class_indices,
        "rotation_gauge_test": gauge,
        "sweep": sweep,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if gauge["min_principal_cosine"] < 1.0 - 1e-7:
        raise RuntimeError(f"external natural-orbital bath failed Q-basis rotation gauge: {gauge}")


if __name__ == "__main__":
    main()
