#!/usr/bin/env python3
"""Validate exact Feshbach elimination of a symmetry-preserving bath into 20Q.

A rotation-covariant natural-orbital bath is constructed from complete active
state classes without experimental input.  The resulting P+Q Hamiltonian is
used only as an independent reference.  Q is then eliminated exactly through
its Schur complement, and every tracked reference eigenpair must satisfy the
energy-dependent 190x190 H_eff(E) acting solely on the fixed N=2 sector of the
20Q P-space.
"""

import json
from pathlib import Path

import numpy as np

from oes.quantum.active_space_blocks import build_helium_s4_p6_20q
from oes.quantum.diffuse_basis import geometric_multi_augment
from oes.quantum.external_dressing import build_external_coupling_space, complete_mo_basis_preserving_active
from oes.quantum.feshbach import eigenpair_downfolding_residual
from oes.quantum.fermions import build_sector_hamiltonian, determinant_basis, transition_one_rdm
from oes.quantum.helium_q1 import classify_states, spatial_transition_rdm, spin_squared_matrix
from oes.quantum.orbital_bath import select_external_natural_bath

ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]
HARTREE_TO_EV = 27.211_386_245_981
TARGET_EXTERNAL_SPATIAL = 13
MANIFOLD_TOL_EV = 1e-4


def active_classes(evals, evecs, p_basis, dip_active):
    e0 = float(evals[0])
    states = classify_states(evals, evecs, spin_squared_matrix(10, p_basis), limit=len(evals))
    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    triplet_e = min(s.excitation_ev for s in triplets)
    triplet = [s for s in triplets if abs(s.excitation_ev - triplet_e) < MANIFOLD_TOL_EV]
    if len(triplet) != 3:
        raise RuntimeError(f"expected 3-state triplet manifold, got {len(triplet)}")

    ground = evecs[:, 0]
    rows = []
    for state in singlets:
        t_spin = transition_one_rdm(evecs[:, state.index], ground, p_basis, 20)
        t_space = spatial_transition_rdm(t_spin)
        mu = np.array([np.sum(dip_active[k] * t_space) for k in range(3)], dtype=complex)
        de = float(evals[state.index] - e0)
        f = (2.0 / 3.0) * de * float(np.sum(np.abs(mu) ** 2))
        rows.append((state, f))
    dark = min((x for x in rows if x[1] < 1e-6), key=lambda x: x[0].excitation_ev)[0]
    bright_e = min(x[0].excitation_ev for x in rows if x[1] > 1e-5)
    bright = [x[0] for x in rows if x[1] > 1e-5 and abs(x[0].excitation_ev - bright_e) < MANIFOLD_TOL_EV]
    if len(bright) != 3:
        raise RuntimeError(f"expected 3-state bright manifold, got {len(bright)}")

    return {
        "ground": ground,
        "triplet": np.column_stack([evecs[:, s.index] for s in triplet]),
        "dark": evecs[:, dark.index],
        "bright": np.column_stack([evecs[:, s.index] for s in bright]),
    }


def class_weights(target, full_vectors, p_dim):
    target = np.asarray(target, dtype=complex)
    if target.ndim == 1:
        target = target[:, None]
    return np.sum(np.abs(target.conj().T @ full_vectors[:p_dim, :]) ** 2, axis=0)


def choose(target, full_vectors, p_dim, count, excluded):
    weights = class_weights(target, full_vectors, p_dim)
    candidates = [int(i) for i in np.argsort(weights)[::-1] if int(i) not in excluded]
    chosen = candidates[:count]
    if len(chosen) != count:
        raise RuntimeError("insufficient full-space states for class match")
    return chosen, [float(weights[i]) for i in chosen]


def main():
    from pyscf import ao2mo

    d_aug = {"He": geometric_multi_augment("He", "aug-cc-pVQZ", extra_layers=1)}
    mol, mf, C_active, receipt = build_helium_s4_p6_20q(
        source_basis=d_aug,
        source_label="d-aug-cc-pVQZ/geometric",
    )
    C_full = complete_mo_basis_preserving_active(mf, C_active)

    h1_p = C_active.T @ mf.get_hcore() @ C_active
    eri_p = ao2mo.kernel(mol, C_active, compact=False).reshape((10,) * 4)
    H_p, p_basis = build_sector_hamiltonian(h1_p, eri_p, n_electrons=2, ecore=float(mol.energy_nuc()))
    vals_p, vecs_p = np.linalg.eigh(H_p)
    p_dim = len(p_basis)
    if p_dim != 190:
        raise RuntimeError(f"unexpected fixed-20Q sector dimension {p_dim}")
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip_p = np.stack([C_active.T @ dip_ao[k] @ C_active for k in range(3)])
    classes = active_classes(vals_p, vecs_p, p_basis, dip_p)

    # Build the complete source response and choose full degenerate orbital
    # classes.  This selector is rotation-covariant and has an independent gauge
    # gate in run_q1_daug_natural_orbital_bath.py.
    nsrc = C_full.shape[1]
    h1_src = C_full.T @ mf.get_hcore() @ C_full
    eri_src = ao2mo.kernel(mol, C_full, compact=False).reshape((nsrc,) * 4)
    external = build_external_coupling_space(
        h1_src,
        eri_src,
        active_basis=p_basis,
        n_active_spatial=10,
        n_electrons=2,
        ecore=float(mol.energy_nuc()),
    )
    bath = select_external_natural_bath(
        external,
        class_states=classes,
        n_active_spatial=10,
        target_external_spatial=TARGET_EXTERNAL_SPATIAL,
        relative_degeneracy_tolerance=2e-5,
        absolute_degeneracy_tolerance=1e-10,
    )
    C_bath = C_full[:, 10:] @ bath.q_rotation
    C = np.column_stack([C_active, C_bath])
    nsp = C.shape[1]
    nspin = 2 * nsp

    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.kernel(mol, C, compact=False).reshape((nsp,) * 4)
    H_raw, basis_raw = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))

    # Reorder determinants so the exact fixed-20Q P sector occupies the first
    # 190 rows/columns.  All remaining determinants are the classical Q bath.
    p_dets = tuple(determinant_basis(20, 2))
    index = {det: i for i, det in enumerate(basis_raw)}
    p_indices = [index[d] for d in p_dets]
    p_set = set(p_indices)
    q_indices = [i for i in range(len(basis_raw)) if i not in p_set]
    order = p_indices + q_indices
    H = H_raw[np.ix_(order, order)]
    p_block_error = float(np.max(np.abs(H[:p_dim, :p_dim] - H_p)))
    if p_block_error > 1e-10:
        raise RuntimeError(f"Feshbach P block mismatch: {p_block_error} Ha")

    vals, vecs = np.linalg.eigh(H)
    excluded = set()
    g, g_ov = choose(classes["ground"], vecs, p_dim, 1, excluded)
    excluded.update(g)
    t, t_ov = choose(classes["triplet"], vecs, p_dim, 3, excluded)
    excluded.update(t)
    d, d_ov = choose(classes["dark"], vecs, p_dim, 1, excluded)
    excluded.update(d)
    b, b_ov = choose(classes["bright"], vecs, p_dim, 3, excluded)

    tracked = {
        "ground": g,
        "triplet": t,
        "dark": d,
        "bright": b,
    }
    residuals = {}
    for name, indices in tracked.items():
        residuals[name] = [
            eigenpair_downfolding_residual(H, p_dim, float(vals[i]), vecs[:, i], singular_floor=1e-9)
            for i in indices
        ]

    ground_e = float(vals[g[0]])
    triplet_components = np.array([(vals[i] - ground_e) * HARTREE_TO_EV for i in t])
    bright_components = np.array([(vals[i] - ground_e) * HARTREE_TO_EV for i in b])
    energies = {
        "triplet": float(np.mean(triplet_components)),
        "dark": float((vals[d[0]] - ground_e) * HARTREE_TO_EV),
        "bright": float(np.mean(bright_components)),
    }
    nist = {
        "triplet": TARGETS["1s2s_3S1"],
        "dark": TARGETS["1s2s_1S0"],
        "bright": TARGETS["1s2p_1P1"],
    }
    nist_res = np.array([energies[k] - nist[k] for k in ("triplet", "dark", "bright")])

    flat = [item for group in residuals.values() for item in group]
    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "EXACT_FESHBACH_SCHUR_DOWNFOLDING_TO_FIXED_20Q_SECTOR",
        "active_protocol": receipt.protocol,
        "active_spin_orbitals": 20,
        "effective_dimension": p_dim,
        "selected_external_spatial_orbitals": bath.selected_external_spatial_orbitals,
        "combined_spatial_orbitals_reference": nsp,
        "combined_spin_orbitals_reference": nspin,
        "combined_dimension_reference": len(basis_raw),
        "integrated_q_dimension": len(q_indices),
        "p_block_error_hartree": p_block_error,
        "class_match_overlaps": {"ground": g_ov, "triplet": t_ov, "dark": d_ov, "bright": b_ov},
        "energies_eV": energies,
        "nist_residuals_eV": {k: float(energies[k] - nist[k]) for k in energies},
        "nist_rms_eV": float(np.sqrt(np.mean(nist_res**2))),
        "triplet_spread_eV": float(np.max(triplet_components) - np.min(triplet_components)),
        "bright_spread_eV": float(np.max(bright_components) - np.min(bright_components)),
        "downfolding": residuals,
        "max_effective_eigen_residual_hartree": max(x["effective_eigen_residual_hartree"] for x in flat),
        "max_q_reconstruction_error": max(x["q_reconstruction_error"] for x in flat),
        "min_distance_to_qhq_spectrum_hartree": min(x["distance_to_qhq_spectrum_hartree"] for x in flat),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if payload["max_effective_eigen_residual_hartree"] > 1e-9:
        raise RuntimeError(f"Feshbach effective-eigenpair gate failed: {payload['max_effective_eigen_residual_hartree']} Ha")
    if payload["max_q_reconstruction_error"] > 1e-8:
        raise RuntimeError(f"Feshbach Q reconstruction gate failed: {payload['max_q_reconstruction_error']}")
    if payload["bright_spread_eV"] > 1e-6 or payload["triplet_spread_eV"] > 1e-6:
        raise RuntimeError("Feshbach reference bath lost rotational/spin manifold degeneracy")


if __name__ == "__main__":
    main()
