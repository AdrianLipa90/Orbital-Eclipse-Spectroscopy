#!/usr/bin/env python3
"""Validate the generic frozen-core reduction on a 20Q active register."""

import json

import numpy as np

from oes.quantum.fermions import build_sector_hamiltonian, sector_dimension
from oes.quantum.frozen_core import frozen_core_effective_hamiltonian


def main():
    try:
        from pyscf import ao2mo, gto, mcscf, scf
        from scipy.linalg import eigh
    except ImportError as exc:
        raise RuntimeError("G1 frozen-core gate requires the q1 extra") from exc

    # Deliberately use a non-benchmark geometry. This gate tests Hamiltonian
    # algebra, not agreement with an experimental LiH observable.
    mol = gto.M(
        atom="Li 0 0 -1.5; H 0 0 1.5",
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
        raise RuntimeError("G1 LiH RHF did not converge")

    ncore = 1
    ncas = 10
    nelecas = 2
    if mf.mo_coeff.shape[1] < ncore + ncas:
        raise RuntimeError("basis does not provide the requested core+20Q active orbitals")

    coeff = np.asarray(mf.mo_coeff[:, : ncore + ncas], dtype=float)
    h_mo = coeff.T @ mf.get_hcore() @ coeff
    eri_mo = ao2mo.kernel(mol, coeff, compact=False).reshape((ncore + ncas,) * 4)

    ours = frozen_core_effective_hamiltonian(
        h_mo,
        eri_mo,
        core_indices=(0,),
        active_indices=tuple(range(1, 11)),
        nuclear_energy=float(mol.energy_nuc()),
    )

    H, _ = build_sector_hamiltonian(
        ours.h1_active,
        ours.eri_active,
        n_electrons=nelecas,
        ecore=ours.ecore,
    )
    oes_energy = float(eigh(H, subset_by_index=[0, 0], eigvals_only=True, driver="evr")[0])

    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.fcisolver.conv_tol = 1e-12
    if mc.ncore != ncore:
        raise RuntimeError(f"unexpected PySCF CASCI core count {mc.ncore}")
    pyscf_h1, pyscf_ecore = mc.get_h1eff()
    pyscf_eri = ao2mo.restore(1, mc.get_h2eff(), ncas)
    pyscf_energy = float(mc.kernel()[0])

    payload = {
        "backend": "SIMULATED_REFERENCE",
        "status_semantics": "G1_FROZEN_CORE_FIXED_20Q_HAMILTONIAN_GATE",
        "molecule": "LiH",
        "geometry_bohr": 3.0,
        "basis": "cc-pVTZ",
        "rhf_energy_hartree": rhf,
        "n_total_electrons": int(mol.nelectron),
        "n_frozen_core_electrons": 2,
        "n_active_electrons": nelecas,
        "n_active_spatial_orbitals": ncas,
        "n_active_spin_orbitals": 2 * ncas,
        "active_qubits": 2 * ncas,
        "active_fixed_particle_dimension": sector_dimension(2 * ncas, nelecas),
        "max_h1_delta_hartree": float(np.max(np.abs(ours.h1_active - pyscf_h1))),
        "max_eri_delta_hartree": float(np.max(np.abs(ours.eri_active - pyscf_eri))),
        "ecore_delta_hartree": float(ours.ecore - float(pyscf_ecore)),
        "oes_active_fci_energy_hartree": oes_energy,
        "pyscf_casci_energy_hartree": pyscf_energy,
        "energy_delta_hartree": float(oes_energy - pyscf_energy),
        "experimental_inputs": [],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

    if payload["active_qubits"] != 20 or payload["active_fixed_particle_dimension"] != 190:
        raise RuntimeError("G1 fixed-20Q active-register gate failed")
    if payload["max_h1_delta_hartree"] > 2e-10:
        raise RuntimeError("G1 frozen-core h1 mismatch against PySCF CASCI")
    if payload["max_eri_delta_hartree"] > 2e-10:
        raise RuntimeError("G1 active ERI mismatch against PySCF CASCI")
    if abs(payload["ecore_delta_hartree"]) > 2e-10:
        raise RuntimeError("G1 frozen-core scalar energy mismatch against PySCF CASCI")
    if abs(payload["energy_delta_hartree"]) > 2e-10:
        raise RuntimeError("G1 OES/PySCF CASCI energy mismatch")


if __name__ == "__main__":
    main()
