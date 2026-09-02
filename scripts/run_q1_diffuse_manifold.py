#!/usr/bin/env python3
"""Resolve the d-aug He I singlet bright manifold and oscillator-strength sum.

This is a source-representation diagnostic, not a 20Q result.  It prints all
low singlet roots with transition oscillator strengths and groups roots that
are nearly degenerate in excitation energy.  The purpose is to distinguish a
true loss of E1 strength from a redistribution across a degenerate 1P manifold.
"""

import json

import numpy as np

from oes.quantum.diffuse_basis import geometric_multi_augment

HARTREE_TO_EV = 27.211_386_245_981
REFERENCE_F = 0.2762


def main():
    from pyscf import ao2mo, fci, gto, scf

    basis = {"He": geometric_multi_augment("He", "aug-cc-pVQZ", extra_layers=1)}
    mol = gto.M(atom="He 0 0 0", basis=basis, unit="Bohr", charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("d-aug helium RHF did not converge")

    C = np.asarray(mf.mo_coeff, dtype=float)
    norb = C.shape[1]
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.kernel(mol, C, compact=False).reshape((norb,) * 4)
    dip_ao = mol.intor("int1e_r", comp=3, hermi=1)
    dip = np.stack([C.T @ dip_ao[k] @ C for k in range(3)])

    solver = fci.direct_spin0.FCI()
    solver.conv_tol = 1e-10
    solver.nroots = 14
    energies, cis = solver.kernel(h1, eri, norb, 2, ecore=float(mol.energy_nuc()))
    energies = np.atleast_1d(np.asarray(energies, dtype=float))
    if not isinstance(cis, (list, tuple)):
        cis = [cis]
    e0 = float(energies[0])
    ci0 = cis[0]

    rows = []
    for i in range(1, len(cis)):
        tdm = np.asarray(solver.trans_rdm1(cis[i], ci0, norb, 2), dtype=float)
        mu = np.array([np.einsum("pq,qp->", dip[k], tdm) for k in range(3)], dtype=float)
        delta_h = float(energies[i] - e0)
        mu2 = float(np.dot(mu, mu))
        f = (2.0 / 3.0) * delta_h * mu2
        rows.append({
            "root": i,
            "excitation_eV": delta_h * HARTREE_TO_EV,
            "oscillator_strength": f,
            "mu_xyz_au": [float(x) for x in mu],
            "mu_norm_au": float(np.sqrt(mu2)),
        })

    bright = [row for row in rows if row["oscillator_strength"] > 1e-5]
    if not bright:
        raise RuntimeError("no bright d-aug singlet roots found")
    first_e = min(row["excitation_eV"] for row in bright)
    # Atomic p components should be essentially degenerate.  Use a deliberately
    # tight numerical window that cannot accidentally combine distinct Rydberg n.
    first_manifold = [row for row in bright if abs(row["excitation_eV"] - first_e) < 1e-4]
    manifold_f = float(sum(row["oscillator_strength"] for row in first_manifold))

    payload = {
        "status_semantics": "SOURCE_MANIFOLD_DIAGNOSTIC_NOT_20Q",
        "basis": "d-aug-cc-pVQZ/geometric",
        "n_spatial_orbitals": norb,
        "ground_hartree": e0,
        "roots": rows,
        "first_bright_manifold": {
            "excitation_eV": first_e,
            "degeneracy_count": len(first_manifold),
            "roots": [row["root"] for row in first_manifold],
            "oscillator_strength_sum": manifold_f,
            "reference_f": REFERENCE_F,
            "f_residual": manifold_f - REFERENCE_F,
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
