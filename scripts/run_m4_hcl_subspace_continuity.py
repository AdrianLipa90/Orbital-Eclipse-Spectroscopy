"""OES-M4 diagnostic: continuity of the selected HCl 20Q active subspace.

This diagnostic is benchmark-blind.  It asks whether the independently selected
10-spatial-orbital active spaces at neighboring H-Cl distances represent the
same smooth physical subspace.  Cross-geometry AO overlaps are used to form
MO-subspace overlap matrices; their singular values are the principal cosines.

No experimental HCl constant is read or used by this script.
"""

from __future__ import annotations

import json

import numpy as np

from oes.quantum.h2_m1 import _select_complete_energy_blocks
from oes.quantum.hcl_m4 import _build_hcl_molecule


def _selected_active_space(bond_bohr: float, basis_name: str = "cc-pVTZ"):
    from pyscf import scf

    mol = _build_hcl_molecule(bond_bohr, basis_name)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"HCl RHF did not converge at R={bond_bohr} bohr")

    relative, group_sizes, _ = _select_complete_energy_blocks(
        np.asarray(mf.mo_energy[5:], dtype=float), target_orbitals=10
    )
    active = tuple(int(i + 5) for i in relative)
    if len(active) != 10:
        raise RuntimeError("HCl continuity diagnostic did not obtain a 10-orbital active space")

    coeff = np.asarray(mf.mo_coeff[:, active], dtype=float)
    return {
        "mol": mol,
        "coeff": coeff,
        "active_indices": active,
        "group_sizes": tuple(int(x) for x in group_sizes),
        "rhf_energy_hartree": energy,
    }


def _principal_cosines(left, right) -> np.ndarray:
    from pyscf import gto

    cross_ao = np.asarray(gto.intor_cross("int1e_ovlp", left["mol"], right["mol"]), dtype=float)
    overlap = left["coeff"].T @ cross_ao @ right["coeff"]
    singular_values = np.linalg.svd(overlap, compute_uv=False)
    if singular_values.shape != (10,) or not np.all(np.isfinite(singular_values)):
        raise RuntimeError("invalid principal-cosine spectrum in HCl active-space continuity audit")
    if np.min(singular_values) < -1e-10 or np.max(singular_values) > 1.0 + 1e-7:
        raise RuntimeError(f"principal cosine outside physical numerical range: {singular_values}")
    return np.clip(singular_values, 0.0, 1.0)


def main() -> None:
    # This is the geometry window reached by the benchmark-blind adaptive
    # bracketing procedure after the original 2.30-bohr seed placed the
    # discrete minimum at the right edge of its first local grid.
    grid = (2.30, 2.36, 2.42, 2.48, 2.54)
    spaces = [_selected_active_space(r) for r in grid]

    pairs = []
    global_min = 1.0
    for i in range(len(grid) - 1):
        s = _principal_cosines(spaces[i], spaces[i + 1])
        global_min = min(global_min, float(np.min(s)))
        pairs.append(
            {
                "left_bohr": grid[i],
                "right_bohr": grid[i + 1],
                "left_active_indices": list(spaces[i]["active_indices"]),
                "right_active_indices": list(spaces[i + 1]["active_indices"]),
                "left_group_sizes": list(spaces[i]["group_sizes"]),
                "right_group_sizes": list(spaces[i + 1]["group_sizes"]),
                "principal_cosines": [float(x) for x in s],
                "minimum_principal_cosine": float(np.min(s)),
                "maximum_one_minus_cosine": float(np.max(1.0 - s)),
                "chordal_distance": float(np.sqrt(np.sum(1.0 - s * s))),
            }
        )

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "basis": "cc-pVTZ",
        "active_qubits": 20,
        "active_spatial_orbitals": 10,
        "geometry_grid_bohr": list(grid),
        "experimental_inputs": [],
        "status_semantics": "M4_HCL_ACTIVE_SUBSPACE_CONTINUITY_DIAGNOSTIC",
        "minimum_principal_cosine_over_adjacent_geometries": global_min,
        "pairs": pairs,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
