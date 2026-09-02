#!/usr/bin/env python3
"""Benchmark Li -> Li+ ionization and Li D-multiplet strength in the same 20Q active space."""

import json
from pathlib import Path

import numpy as np

from oes.quantum.fermions import build_sector_hamiltonian
from oes.quantum.helium_q1 import pyscf_fci_reference
from oes.quantum.lithium_a1 import HARTREE_TO_EV, prepare_lithium_integrals, run_lithium_a1

ROOT = Path(__file__).resolve().parents[1]
TARGET = json.loads((ROOT / "benchmarks" / "lithium_a1_nist.json").read_text())["neutral_lithium"]


def main():
    neutral = run_lithium_a1(basis_name="cc-pVTZ", n_spatial=10)
    mol, h1, eri, _dip, _rohf, selected, group_sizes = prepare_lithium_integrals(
        basis_name="cc-pVTZ",
        n_spatial=10,
    )

    # Li+ has two electrons but uses exactly the same selected one-particle
    # active space as neutral Li. This isolates the N=3 -> N=2 color difference.
    H2, basis2 = build_sector_hamiltonian(
        h1,
        eri,
        n_electrons=2,
        ecore=float(mol.energy_nuc()),
    )
    e_li_plus = float(np.linalg.eigvalsh(H2)[0])
    e_li_plus_ref = pyscf_fci_reference(
        h1,
        eri,
        n_electrons=2,
        ecore=float(mol.energy_nuc()),
    )
    ionization_ev = (e_li_plus - neutral.oes_fci_energy_hartree) * HARTREE_TO_EV
    f_ref = float(TARGET["d_multiplet_absorption_f_sum_approx"])

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "A1_1_FIXED_20Q_CORE_SCREENING_AND_IONIZATION_DIAGNOSTIC",
        "active_protocol": neutral.active_protocol,
        "selected_mo_indices": list(selected),
        "selected_group_sizes": list(group_sizes),
        "neutral_fci_hartree": neutral.oes_fci_energy_hartree,
        "li_plus_fci_hartree": e_li_plus,
        "li_plus_pyscf_fci_hartree": e_li_plus_ref,
        "li_plus_fci_delta_hartree": e_li_plus - e_li_plus_ref,
        "li_plus_fixed_particle_dimension": len(basis2),
        "ionization_eV": ionization_ev,
        "nist_ionization_eV": float(TARGET["ionization_energy_eV"]),
        "ionization_nist_residual_eV": ionization_ev - float(TARGET["ionization_energy_eV"]),
        "bright_excitation_eV": neutral.first_bright_excitation_ev,
        "bright_nist_residual_eV": neutral.first_bright_excitation_ev - float(TARGET["first_2s_2p_excitation_eV"]),
        "bright_f_sum": neutral.first_bright_f_sum,
        "bright_f_reference_approx": f_ref,
        "bright_f_relative_residual": neutral.first_bright_f_sum / f_ref - 1.0,
        "bright_degeneracy": neutral.first_bright_degeneracy,
        "bright_spread_eV": neutral.first_bright_spread_ev,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if abs(payload["li_plus_fci_delta_hartree"]) > 1e-9:
        raise RuntimeError(f"A1 Li+ OES/PySCF FCI mismatch: {payload['li_plus_fci_delta_hartree']} Ha")
    if payload["li_plus_fixed_particle_dimension"] != 190:
        raise RuntimeError("A1 Li+ fixed-N sector dimension gate failed")
    if payload["bright_degeneracy"] != 3 or payload["bright_spread_eV"] > 1e-8:
        raise RuntimeError("A1 Li bright-manifold symmetry gate failed")


if __name__ == "__main__":
    main()
