"""Frozen-core reduction for fixed-register OES active-space Hamiltonians.

A set of doubly occupied inactive spatial orbitals is integrated into a scalar
core energy and an effective one-electron operator on the selected active
spatial orbitals.  The active two-electron integrals are retained unchanged.

For real spatial orbitals and chemists' ERI notation ``(pq|rs)``:

    E_core = E_nuc
           + 2 sum_i h_ii
           + sum_ij [2 (ii|jj) - (ij|ji)]

    h_eff[p,q] = h[p,q]
               + sum_i [2 (pq|ii) - (pi|iq)]

where i,j run over the frozen doubly occupied core and p,q over the active
space.  No experimental observable enters this reduction.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class FrozenCoreHamiltonian:
    h1_active: np.ndarray
    eri_active: np.ndarray
    ecore: float
    core_indices: Tuple[int, ...]
    active_indices: Tuple[int, ...]

    def metadata(self) -> Dict[str, object]:
        return {
            "ecore": float(self.ecore),
            "core_indices": list(self.core_indices),
            "active_indices": list(self.active_indices),
            "n_core_spatial": len(self.core_indices),
            "n_active_spatial": len(self.active_indices),
            "n_frozen_electrons": 2 * len(self.core_indices),
        }


def _validated_indices(indices: Sequence[int], n_spatial: int, name: str) -> Tuple[int, ...]:
    out = tuple(int(i) for i in indices)
    if len(set(out)) != len(out):
        raise ValueError(f"{name} contains duplicate orbital indices")
    if any(i < 0 or i >= n_spatial for i in out):
        raise ValueError(f"{name} contains an orbital outside the supplied integral space")
    return out


def frozen_core_effective_hamiltonian(
    h1: np.ndarray,
    eri: np.ndarray,
    core_indices: Sequence[int],
    active_indices: Sequence[int],
    nuclear_energy: float = 0.0,
) -> FrozenCoreHamiltonian:
    """Integrate a closed-shell frozen core into an active-space Hamiltonian.

    ``h1`` and ``eri`` are spatial-orbital integrals in one common orthonormal
    basis.  ``eri`` must use chemists' notation ``(pq|rs)``.  Core orbitals are
    assumed doubly occupied.  Core and active sets must be disjoint; orbitals in
    neither set are external and do not enter the returned Hamiltonian.
    """
    h1 = np.asarray(h1, dtype=float)
    eri = np.asarray(eri, dtype=float)
    if h1.ndim != 2 or h1.shape[0] != h1.shape[1]:
        raise ValueError("h1 must be square")
    n_spatial = h1.shape[0]
    if eri.shape != (n_spatial,) * 4:
        raise ValueError("eri shape must be (n,n,n,n)")
    if not np.all(np.isfinite(h1)) or not np.all(np.isfinite(eri)):
        raise ValueError("integrals must be finite")
    if not np.isfinite(float(nuclear_energy)):
        raise ValueError("nuclear_energy must be finite")

    core = _validated_indices(core_indices, n_spatial, "core_indices")
    active = _validated_indices(active_indices, n_spatial, "active_indices")
    if not active:
        raise ValueError("active_indices must be non-empty")
    if set(core) & set(active):
        raise ValueError("core and active orbital sets must be disjoint")

    ecore = float(nuclear_energy)
    for i in core:
        ecore += 2.0 * float(h1[i, i])
    for i in core:
        for j in core:
            ecore += 2.0 * float(eri[i, i, j, j]) - float(eri[i, j, j, i])

    h_eff = np.asarray(h1[np.ix_(active, active)], dtype=float).copy()
    for a, p in enumerate(active):
        for b, q in enumerate(active):
            correction = 0.0
            for i in core:
                correction += 2.0 * float(eri[p, q, i, i]) - float(eri[p, i, i, q])
            h_eff[a, b] += correction

    eri_active = np.asarray(eri[np.ix_(active, active, active, active)], dtype=float).copy()
    h_eff = 0.5 * (h_eff + h_eff.T)

    return FrozenCoreHamiltonian(
        h1_active=h_eff,
        eri_active=eri_active,
        ecore=float(ecore),
        core_indices=core,
        active_indices=active,
    )
