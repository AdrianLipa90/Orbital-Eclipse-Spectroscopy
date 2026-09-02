#!/usr/bin/env python3
"""Compare aug and geometric d-aug source spectra before 20Q compression."""

import json
from pathlib import Path

from oes.quantum.diffuse_basis import geometric_multi_augment
from oes.quantum.source_spectrum import helium_source_spectrum


ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]
REFERENCE_F = 0.2762


def decorate(result):
    data = result.as_dict()
    data["residuals_eV"] = {
        "1s2s_3S1": result.triplet_excitation_ev - TARGETS["1s2s_3S1"],
        "1s2s_1S0": result.dark_singlet_excitation_ev - TARGETS["1s2s_1S0"],
        "1s2p_1P1": result.bright_singlet_excitation_ev - TARGETS["1s2p_1P1"],
    }
    data["bright_f_reference"] = REFERENCE_F
    data["bright_f_residual"] = result.bright_oscillator_strength - REFERENCE_F
    return data


def main():
    aug = helium_source_spectrum("aug-cc-pVQZ", "aug-cc-pVQZ/full-source")
    daug_basis = {"He": geometric_multi_augment("He", "aug-cc-pVQZ", extra_layers=1)}
    daug = helium_source_spectrum(daug_basis, "d-aug-cc-pVQZ/geometric/full-source")
    payload = {
        "status_semantics": "SOURCE_REPRESENTATION_DIAGNOSTIC_NOT_20Q",
        "targets_eV": TARGETS,
        "cases": [decorate(aug), decorate(daug)],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
