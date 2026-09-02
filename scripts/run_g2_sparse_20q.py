#!/usr/bin/env python3
"""Stress the sparse exact 20Q reference sector at six active electrons."""

import json
import time

import numpy as np

from oes.quantum.frozen_core import frozen_core_effective_hamiltonian
from oes.quantum.sparse_ms import (
    build_sparse_fixed_spin_hamiltonian,
    fixed_spin_max_connectivity,
)


def main():
    try:
        from pyscf import ao2mo, gto, mcscf, scf
        from scipy.sparse.linalg import eigsh
    except ImportError as exc:
        raise RuntimeError("G2 gate requires the q1 extra") from exc

    # Non-benchmark geometry by design: this is a computational/algebraic gate.
    mol = gto.M(
        atom="H 0 0 -0.9; F 0 0 0.9",
        unit="Bohr",
        basis="cc-pVTZ",
        charge=0,
        spin=0,
        symmetry=False,
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    rhf = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError("G2 HF RHF did not converge")

    ncore = 2
    ncas = 10
    nelecas = 6
    coeff = np.asarray(mf.mo_coeff[:, : ncore + ncas], dtype=float)
    h_mo = coeff.T @ mf.get_hcore() @ coeff
    eri_mo = ao2mo.kernel(mol, coeff, compact=False).reshape((ncore + ncas,) * 4)
    active = frozen_core_effective_hamiltonian(
        h_mo,
        eri_mo,
        core_indices=(0, 1),
        active_indices=tuple(range(2, 12)),
        nuclear_energy=float(mol.energy_nuc()),
    )

    t0 = time.perf_counter()
    H, basis = build_sparse_fixed_spin_hamiltonian(
        active.h1_active,
        active.eri_active,
        n_alpha=3,
        n_beta=3,
        ecore=active.ecore,
    )
    build_seconds = time.perf_counter() - t0

    t1 = time.perf_counter()
    evals = eigsh(H, k=1, which="SA", return_eigenvectors=False, tol=1e-11, maxiter=10000)
    oes_energy = float(evals[0])
    diagonalize_seconds = time.perf_counter() - t1

    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.fcisolver.conv_tol = 1e-12
    if mc.ncore != ncore:
        raise RuntimeError(f"unexpected PySCF core count {mc.ncore}")
    pyscf_energy = float(mc.kernel()[0])

    csr_bytes = int(H.data.nbytes + H.indices.nbytes + H.indptr.nbytes)
    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "G2_SPARSE_FIXED_MS_20Q_SIX_ELECTRON_GATE",
        "molecule": "HF",
        "geometry_bohr": 1.8,
        "geometry_role": "NON_BENCHMARK_ALGEBRAIC_STRESS_GATE",
        "basis": "cc-pVTZ",
        "rhf_energy_hartree": rhf,
        "active_qubits": 20,
        "frozen_core_electrons": 4,
        "active_electrons": 6,
        "n_alpha": 3,
        "n_beta": 3,
        "fixed_ms_dimension": len(basis),
        "maximum_connectivity_per_row": fixed_spin_max_connectivity(10, 3, 3),
        "csr_nnz": int(H.nnz),
        "csr_density": float(H.nnz / (len(basis) ** 2)),
        "csr_storage_bytes": csr_bytes,
        "build_seconds": float(build_seconds),
        "diagonalize_seconds": float(diagonalize_seconds),
        "oes_sparse_energy_hartree": oes_energy,
        "pyscf_casci_energy_hartree": pyscf_energy,
        "energy_delta_hartree": float(oes_energy - pyscf_energy),
        "experimental_inputs": [],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if payload["fixed_ms_dimension"] != 14400:
        raise RuntimeError("G2 six-electron 20Q sector dimension gate failed")
    if abs(payload["energy_delta_hartree"]) > 2e-9:
        raise RuntimeError("G2 sparse OES energy does not reproduce PySCF CASCI")
    if payload["csr_density"] >= 0.1:
        raise RuntimeError("G2 sparse operator unexpectedly dense")


if __name__ == "__main__":
    main()
