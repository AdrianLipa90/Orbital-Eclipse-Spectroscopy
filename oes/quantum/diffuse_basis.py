"""Geometric diffuse augmentation for helium Rydberg reference spaces.

The construction follows the standard multiple-augmentation idea: for each
angular momentum present in an aug-cc-pVXZ basis, append one uncontracted shell
whose exponent continues the geometric progression of the two most diffuse
existing primitive exponents.  One application produces d-aug; two produce
t-aug.  This changes the source orbital representation, not the final Q1 qubit
budget.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Dict, List


def geometric_multi_augment(element: str, base_basis: str, extra_layers: int = 1):
    """Return a PySCF basis list with `extra_layers` additional diffuse shells."""
    if extra_layers < 1:
        raise ValueError("extra_layers must be >= 1")
    try:
        from pyscf import gto
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("diffuse-basis construction requires the OES q1 extra (PySCF)") from exc

    basis = deepcopy(gto.basis.load(base_basis, element))
    by_l: Dict[int, List[float]] = {}
    for shell in basis:
        l = int(shell[0])
        for primitive in shell[1:]:
            exponent = float(primitive[0])
            by_l.setdefault(l, []).append(exponent)

    additions = []
    for l, values in sorted(by_l.items()):
        exps = sorted(set(values))
        if len(exps) < 2:
            raise RuntimeError(f"need at least two primitive exponents for l={l}")
        e0, e1 = exps[0], exps[1]
        generated = []
        for _ in range(extra_layers):
            new_exp = e0 * e0 / e1
            if not (0.0 < new_exp < e0):
                raise RuntimeError(f"invalid geometric diffuse exponent for l={l}: {new_exp}")
            generated.append(new_exp)
            e1, e0 = e0, new_exp
        for exponent in generated:
            additions.append([l, [float(exponent), 1.0]])

    return basis + additions
