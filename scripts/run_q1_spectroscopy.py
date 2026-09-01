#!/usr/bin/env python3
"""Evaluate the held-out He I 1s2p 1P bright-line prediction."""

import json
from pathlib import Path

from oes.quantum.helium_spectroscopy import first_bright_singlet_prediction


ROOT = Path(__file__).resolve().parents[1]
TARGETS = json.loads((ROOT / "benchmarks" / "helium_q1_nist.json").read_text())["targets_eV"]


if __name__ == "__main__":
    prediction, candidates = first_bright_singlet_prediction()
    target = TARGETS["1s2p_1P1"]
    payload = prediction.as_dict()
    payload["held_out_target_eV"] = target
    payload["residual_eV"] = prediction.excitation_ev - target
    payload["abs_residual_eV"] = abs(prediction.excitation_ev - target)
    payload["candidate_count"] = len(candidates)
    print(json.dumps(payload, indent=2, sort_keys=True))
