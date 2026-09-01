#!/usr/bin/env python3
"""Dress the predictive d-aug s4+p6 20Q states with parameter-free EN2 Q-space energy corrections."""

import json
from pathlib import Path

import numpy as np

from oes.quantum.active_space_blocks import build_helium_s4_p6_20q
from oes.quantum.diffuse_basis import geometric_multi_augment
from oes.quantum.external_dressing import (
    build_external_coupling_space,
    complete_mo_basis_preserving_active,
    en2_correction,
)
from oes.quantum.fermions import build_sector_hamiltonian, transition_one_rdm
from oes.quantum.helium_q1 import classify_states, spatial_transition_rdm, spin_squared_matrix

ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]
HARTREE_TO_EV = 27.211_386_245_981
SOURCE = {
    "triplet": 19.80281611421936,
    "dark": 20.61327300402567,
    "bright": 21.29949008913645,
}


def main():
    from pyscf import ao2mo

    d_aug = {"He": geometric_multi_augment("He", "aug-cc-pVQZ", extra_layers=1)}
    mol, mf, C_active, receipt = build_helium_s4_p6_20q(
        source_basis=d_aug,
        source_label="d-aug-cc-pVQZ/geometric",
    )

    # Complete the one-particle source basis without rotating the ten active
    # columns. First ten spatial modes are exactly the 20Q P-space.
    C_full = complete_mo_basis_preserving_active(mf, C_active)
    nfull = C_full.shape[1]
    h1_full = C_full.T @ mf.get_hcore() @ C_full
    eri_full = ao2mo.kernel(mol, C_full, compact=False).reshape((nfull,) * 4)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip_active = np.stack([C_active.T @ dip_ao[k] @ C_active for k in range(3)])

    h1 = h1_full[:10, :10]
    eri = eri_full[:10, :10, :10, :10]
    H, basis = build_sector_hamiltonian(h1, eri, n_electrons=2, ecore=float(mol.energy_nuc()))
    evals, evecs = np.linalg.eigh(H)
    e0 = float(evals[0])

    s2mat = spin_squared_matrix(10, basis)
    states = classify_states(evals, evecs, s2mat, limit=min(180, len(evals)))
    triplets = [s for s in states[1:] if abs(s.s2 - 2.0) < 1e-6]
    singlets = [s for s in states[1:] if abs(s.s2) < 1e-6]
    if not triplets or not singlets:
        raise RuntimeError("EN2 diagnostic failed to resolve active spin sectors")
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
        raise RuntimeError("EN2 diagnostic lost dark/bright active classes")
    dark = min(dark_rows, key=lambda row: row["state"].excitation_ev)["state"]
    first_bright_e = min(row["state"].excitation_ev for row in bright_rows)
    bright_manifold = [row["state"] for row in bright_rows if abs(row["state"].excitation_ev - first_bright_e) < 1e-4]
    if len(bright_manifold) != 3:
        raise RuntimeError(f"EN2 symmetry gate expected 3 bright components, got {len(bright_manifold)}")

    external = build_external_coupling_space(
        h1_full,
        eri_full,
        active_basis=basis,
        n_active_spatial=10,
        n_electrons=2,
        ecore=float(mol.energy_nuc()),
    )

    selected = {
        "ground": 0,
        "triplet": triplet.index,
        "dark": dark.index,
        "bright_0": bright_manifold[0].index,
        "bright_1": bright_manifold[1].index,
        "bright_2": bright_manifold[2].index,
    }
    dressed = {}
    for name, idx in selected.items():
        diag = en2_correction(float(evals[idx]), evecs[:, idx], external, denominator_floor=1e-5)
        dressed[name] = {
            **diag,
            "active_energy_hartree": float(evals[idx]),
            "dressed_energy_hartree": float(evals[idx] + diag["correction_hartree"]),
        }

    dressed_ground = dressed["ground"]["dressed_energy_hartree"]
    active_ground = float(evals[0])
    def excitation(name):
        return (dressed[name]["dressed_energy_hartree"] - dressed_ground) * HARTREE_TO_EV
    def active_excitation(name):
        return (dressed[name]["active_energy_hartree"] - active_ground) * HARTREE_TO_EV

    bright_dressed = [excitation(f"bright_{i}") for i in range(3)]
    bright_active = [active_excitation(f"bright_{i}") for i in range(3)]
    if max(bright_dressed) - min(bright_dressed) > 1e-4:
        raise RuntimeError(f"EN2 dressing broke bright degeneracy: spread={max(bright_dressed)-min(bright_dressed)} eV")

    classes = {
        "triplet": {
            "active_eV": active_excitation("triplet"),
            "dressed_eV": excitation("triplet"),
            "source_eV": SOURCE["triplet"],
            "nist_eV": TARGETS["1s2s_3S1"],
        },
        "dark": {
            "active_eV": active_excitation("dark"),
            "dressed_eV": excitation("dark"),
            "source_eV": SOURCE["dark"],
            "nist_eV": TARGETS["1s2s_1S0"],
        },
        "bright": {
            "active_eV": float(np.mean(bright_active)),
            "dressed_eV": float(np.mean(bright_dressed)),
            "source_eV": SOURCE["bright"],
            "nist_eV": TARGETS["1s2p_1P1"],
            "dressed_manifold_spread_eV": float(max(bright_dressed) - min(bright_dressed)),
        },
    }
    for item in classes.values():
        item["active_source_residual_eV"] = item["active_eV"] - item["source_eV"]
        item["dressed_source_residual_eV"] = item["dressed_eV"] - item["source_eV"]
        item["active_nist_residual_eV"] = item["active_eV"] - item["nist_eV"]
        item["dressed_nist_residual_eV"] = item["dressed_eV"] - item["nist_eV"]

    active_source = np.array([classes[k]["active_source_residual_eV"] for k in ("triplet", "dark", "bright")])
    dressed_source = np.array([classes[k]["dressed_source_residual_eV"] for k in ("triplet", "dark", "bright")])
    active_nist = np.array([classes[k]["active_nist_residual_eV"] for k in ("triplet", "dark", "bright")])
    dressed_nist = np.array([classes[k]["dressed_nist_residual_eV"] for k in ("triplet", "dark", "bright")])

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "PARAMETER_FREE_EXTERNAL_SPACE_EN2_DIAGNOSTIC",
        "active_protocol": receipt.protocol,
        "n_active_spin_orbitals": 20,
        "n_source_spin_orbitals": 2 * nfull,
        "external_determinants": len(external.external_basis),
        "classes": classes,
        "state_corrections": dressed,
        "metrics": {
            "active_source_rms_eV": float(np.sqrt(np.mean(active_source**2))),
            "dressed_source_rms_eV": float(np.sqrt(np.mean(dressed_source**2))),
            "active_source_centered_rms_eV": float(np.sqrt(np.mean((active_source - np.mean(active_source))**2))),
            "dressed_source_centered_rms_eV": float(np.sqrt(np.mean((dressed_source - np.mean(dressed_source))**2))),
            "active_nist_rms_eV": float(np.sqrt(np.mean(active_nist**2))),
            "dressed_nist_rms_eV": float(np.sqrt(np.mean(dressed_nist**2))),
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
