# OES-M1 H2 Validation Receipt

Date: 2026-09-02
Branch: `oes-m1-h2-20q`
Parent: `oes-a1-lithium-20q`
Backend: `SIMULATED_REFERENCE`
Status: `PARTIAL_PASS__H2_TWO_CENTRE_20Q_BASELINE`

## Canonical register

H2 uses the same ten-spatial / twenty-spin-orbital register as Q1 and A1.

- electrons: 2
- full qubit dimension: 1,048,576
- exact fixed-N sector: C(20,2) = 190 determinants

## Molecular active-space rule

Protocol: `MOLECULAR-COMPLETE-DEGENERACY-BLOCKS-20Q`.

The selector uses only canonical orbital-energy degeneracies and admits whole groups until exactly ten spatial orbitals are present. No experimental bond length, vibrational frequency, line energy or oscillator strength enters Hamiltonian construction or orbital selection.

At the NIST equilibrium-distance benchmark point the selected D2h irreps are:

`Ag, B1u, Ag, B1u, B2u, B3u, Ag, B2g, B3g, B1u`

with complete degeneracy-group sizes

`1 + 1 + 1 + 1 + 2 + 1 + 2 + 1 = 10`.

## Exact implementation and symmetry gates

Hosted GitHub Actions gate `OES M1 H2 20Q`, run 33571405432, head `3556f8c96dad8b3b0d3fa249094a5c891d9e3873` — PASS.

At R = 1.40111853846 bohr (0.74144 Å benchmark geometry):

- OES FCI: -1.1566093705724203 Ha
- independent PySCF FCI: -1.1566093705724185 Ha
- residual: -1.78e-15 Ha — PASS
- ground inversion parity: +1 within 2e-15 — PASS
- first dipole-bright state inversion parity: -1 within 1e-15 — PASS
- ground S^2: 8.67e-34 — PASS
- centre-A `1/r_A` exposure: 1.799879987426903
- centre-B `1/r_B` exposure: 1.7998799874269027
- equivalent-centre exposure residual: 2.22e-16 — PASS

The first active-space bright electronic transition is 13.31983388 eV with summed oscillator strength 0.49670353. It is used here as a parity/transition-density gate; no experimental electronic-line promotion is made from this value.

## Born-Oppenheimer local curve

Fixed-20Q energies were computed at R/bohr:

`1.25, 1.325, 1.40, 1.475, 1.55`.

A local quadratic fit gives:

- fitted equilibrium distance: 1.447389224 bohr = 0.765925392 Å
- NIST `r_e`: 0.74144 Å
- residual: +0.024485392 Å (~3.30%)
- fitted curvature: 0.369506152 Ha/bohr^2
- harmonic wavenumber: 4403.071805 cm^-1
- NIST `omega_e`: 4401.213 cm^-1
- residual: +1.858805 cm^-1 (~0.0422%)

The 6-bohr dissociation probe gives -0.9989465891 Ha and a finite-grid electronic well-depth diagnostic of 4.3014625 eV. The dissociation limit is not promoted as a precision molecular constant at this milestone.

## Verdict

`M1_IMPLEMENTATION_PASS` — the same 20Q two-electron fermionic core reproduces an independent molecular FCI implementation at machine precision.

`M1_TWO_CENTRE_FLAVOR_PASS` — inversion parity, opposite-parity dipole visibility and equal equivalent-centre exposures are recovered directly from the many-electron state.

`M1_VIBRATIONAL_CURVATURE_PASS_AT_BASELINE` — the local harmonic curvature reproduces the cited NIST `omega_e` within ~0.05% on this fixed active space.

`M1_GEOMETRY_PARTIAL` — the fitted equilibrium distance is high by ~0.0245 Å and remains an active-space/source-representation refinement target.

Broader H2 spectroscopy, dissociation convergence, rotational constants, external-bath refinement and physical-QPU execution remain OPEN.
