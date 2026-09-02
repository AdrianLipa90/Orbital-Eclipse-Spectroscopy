#!/usr/bin/env python3
"""Scale the HF fixed-20Q reference to eight active electrons (44,100 states)."""

import json
import time

import numpy as np

from oes.quantum.frozen_core import frozen_core_effective_hamiltonian
from oes.quantum.h2_m1 import _select_complete_energy_blocks
from oes.quantum.hf_m3 import _build_hf_molecule
from oes.quantum.sparse_ms import build_sparse_fixed_spin_hamiltonian, fixed_spin_max_connectivity


def main():
    try:
        from pyscf import ao2mo, fci, scf
        from scipy.sparse.linalg import eigsh
    except ImportError as exc:
        raise RuntimeError("M3 HF 8e stress gate requires the q1 extra") from exc

    # Deliberately not an experimental equilibrium geometry.  This is a pure
    # active-electron convergence/scaling gate.
    bond_bohr = 1.8
    mol = _build_hf_molecule(bond_bohr, "cc-pVTZ")
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    rhf = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError("HF 8e stress RHF did not converge")

    # Freeze only the deepest F-dominated core spatial orbital (2 electrons).
    core = (0,)
    relative, group_sizes, _ = _select_complete_energy_blocks(
        np.asarray(mf.mo_energy[1:], dtype=float), target_orbitals=10
    )
    active = tuple(int(i + 1) for i in relative)
    selected = core + active
    if len(active) != 10 or set(core) & set(active):
        raise RuntimeError("HF 8e complete-block active partition failed")

    coeff = np.asarray(mf.mo_coeff[:, selected], dtype=float)
    h = coeff.T @ mf.get_hcore() @ coeff
    eri = ao2mo.kernel(mol, coeff, compact=False).reshape((11,) * 4)
    reduced = frozen_core_effective_hamiltonian(
        h,
        eri,
        core_indices=(0,),
        active_indices=tuple(range(1, 11)),
        nuclear_energy=float(mol.energy_nuc()),
    )

    t0 = time.perf_counter()
    H, basis = build_sparse_fixed_spin_hamiltonian(
        reduced.h1_active,
        reduced.eri_active,
        n_alpha=4,
        n_beta=4,
        ecore=reduced.ecore,
    )
    build_seconds = time.perf_counter() - t0

    t1 = time.perf_counter()
    oes_energy = float(eigsh(H, k=1, which="SA", return_eigenvectors=False, tol=2e-10, maxiter=15000)[0])
    diagonalize_seconds = time.perf_counter() - t1

    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 2e-10
    ref_energy, _ = solver.kernel(
        reduced.h1_active,
        reduced.eri_active,
        10,
        (4, 4),
        ecore=reduced.ecore,
    )
    ref_energy = float(ref_energy)

    csr_bytes = int(H.data.nbytes + H.indices.nbytes + H.indptr.nbytes)
    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "M3_HF_FROZEN_2E_EIGHT_ACTIVE_ELECTRON_20Q_SCALING_GATE",
        "geometry_bohr": bond_bohr,
        "geometry_role": "NON_BENCHMARK_SCALING_GATE",
        "basis": "cc-pVTZ",
        "rhf_energy_hartree": rhf,
        "core_mo_indices": list(core),
        "active_mo_indices": list(active),
        "active_group_sizes": list(group_sizes),
        "frozen_electrons": 2,
        "active_electrons": 8,
        "active_qubits": 20,
        "n_alpha": 4,
        "n_beta": 4,
        "fixed_ms_dimension": len(basis),
        "maximum_connectivity_per_row": fixed_spin_max_connectivity(10, 4, 4),
        "csr_nnz": int(H.nnz),
        "csr_density": float(H.nnz / (len(basis) ** 2)),
        "csr_storage_bytes": csr_bytes,
        "build_seconds": float(build_seconds),
        "diagonalize_seconds": float(diagonalize_seconds),
        "oes_sparse_energy_hartree": oes_energy,
        "pyscf_active_fci_energy_hartree": ref_energy,
        "energy_delta_hartree": float(oes_energy - ref_energy),
        "experimental_inputs": [],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if len(basis) != 44100:
        raise RuntimeError("HF 8e fixed-Ms 20Q dimension gate failed")
    if abs(payload["energy_delta_hartree"]) > 5e-8:
        raise RuntimeError("HF 8e OES sparse energy does not reproduce independent active FCI")
    if payload["csr_density"] >= 0.1:
        raise RuntimeError("HF 8e operator unexpectedly dense")


if __name__ == "__main__":
    main()
