#!/usr/bin/env python3
"""Compare several 10-spatial-orbital helium spaces without changing 20Q budget."""

import json
from pathlib import Path

from oes.quantum.helium_q1 import HARTREE_TO_EV, run_helium_q1


TARGET_PATH = Path(__file__).resolve().parents[1] / "benchmarks" / "helium_q1_nist.json"
CANDIDATES = ["cc-pVTZ", "cc-pVQZ", "aug-cc-pVTZ", "aug-cc-pVQZ"]


def main():
    targets = json.loads(TARGET_PATH.read_text())["targets_eV"]
    triplet_target = targets["1s2s_3S1"]
    singlet_target = targets["1s2s_1S0"]
    rows = []
    for basis_name in CANDIDATES:
        result, _ = run_helium_q1(basis_name=basis_name)
        e0 = result.oes_fci_energy_hartree
        triplet_ev = (result.first_triplet_energy_hartree - e0) * HARTREE_TO_EV
        singlet_ev = (result.first_singlet_excited_energy_hartree - e0) * HARTREE_TO_EV
        dt = triplet_ev - triplet_target
        ds = singlet_ev - singlet_target
        rms = ((dt * dt + ds * ds) / 2.0) ** 0.5
        rows.append({
            "basis": basis_name,
            "n_spin_orbitals": 20,
            "ground_hartree": e0,
            "triplet_excitation_eV": triplet_ev,
            "singlet_excitation_eV": singlet_ev,
            "triplet_residual_eV": dt,
            "singlet_residual_eV": ds,
            "rms_two_level_residual_eV": rms,
        })
    rows.sort(key=lambda row: row["rms_two_level_residual_eV"])
    print(json.dumps({"targets_eV": targets, "candidates": rows, "best": rows[0]}, indent=2))


if __name__ == "__main__":
    main()
