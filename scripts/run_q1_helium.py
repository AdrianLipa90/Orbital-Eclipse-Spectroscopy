#!/usr/bin/env python3
"""Run OES-Q1 helium 20-spin-orbital benchmark and print a JSON receipt."""

import json

from oes.quantum.helium_q1 import run_helium_q1


if __name__ == "__main__":
    result, states = run_helium_q1()
    payload = result.as_dict()
    payload["states"] = [state.__dict__ for state in states[:12]]
    print(json.dumps(payload, indent=2, sort_keys=True))
