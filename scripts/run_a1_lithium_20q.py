#!/usr/bin/env python3
"""Run OES-A1 neutral-lithium fixed-20Q implementation and first bright-line gate."""

import json
from pathlib import Path

from oes.quantum.lithium_a1 import run_lithium_a1

ROOT = Path(__file__).resolve().parents[1]
TARGET = json.loads((ROOT / "benchmarks" / "lithium_a1_nist.json").read_text())["neutral_lithium"]


def main():
    result = run_lithium_a1(basis_name="cc-pVTZ", n_spatial=10)
    payload = result.as_dict()
    payload["first_bright_nist_eV"] = TARGET["first_2s_2p_excitation_eV"]
    payload["first_bright_nist_residual_eV"] = (
        payload["first_bright_excitation_ev"] - TARGET["first_2s_2p_excitation_eV"]
    )
    payload["status_semantics"] = "A1_0_FIXED_20Q_THREE_ELECTRON_IMPLEMENTATION_GATE"
    print(json.dumps(payload, indent=2, sort_keys=True))

    if abs(payload["fci_delta_hartree"]) > 1e-9:
        raise RuntimeError(f"A1 OES/PySCF FCI mismatch: {payload['fci_delta_hartree']} Ha")
    if payload["fixed_particle_dimension"] != 1140 or payload["ms_half_dimension"] != 450:
        raise RuntimeError("A1 canonical sector dimension gate failed")
    if payload["first_bright_degeneracy"] != 3:
        raise RuntimeError(
            f"A1 expected threefold nonrelativistic first bright orbital manifold, got {payload['first_bright_degeneracy']}"
        )


if __name__ == "__main__":
    main()
