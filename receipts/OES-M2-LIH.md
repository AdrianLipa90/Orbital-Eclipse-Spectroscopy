# OES-M2 LiH Validation Receipt

Date: 2026-09-02
Branch: `oes-m2-lih-20q`
Parent: `oes-m1-h2-20q`
Backend: `SIMULATED_REFERENCE`
Status: `PARTIAL_PASS__LIH_HETERONUCLEAR_FIXED_20Q_BASELINE`

## Canonical register

LiH uses the same ten-spatial / twenty-spin-orbital register established by Q1, A1 and M1.

- electrons: 4
- full qubit dimension: 1,048,576
- exact fixed-N sector: C(20,4) = 4,845 determinants
- exact M_S=0 sector: C(10,2)^2 = 2,025 determinants

The hosted core suite validates the four-electron particle sector and every determinant in the M_S=0 reference basis carries N_alpha=N_beta=2.

## Molecular active-space rule

Protocol: `HETERONUCLEAR-COMPLETE-DEGENERACY-BLOCKS-20Q`.

The selector uses canonical orbital-energy degeneracy groups and admits whole groups until exactly ten spatial orbitals are present. Experimental molecular constants remain benchmark outputs and do not enter active-space selection or Hamiltonian construction.

At the NIST equilibrium-distance benchmark geometry the selected group sizes are

`1 + 1 + 1 + 2 + 1 + 2 + 1 + 1 = 10`,

with selected irreps

`A1, A1, A1, E1x, E1y, A1, E1x, E1y, A1, A1`.

## Hosted implementation gate

Hosted GitHub Actions gate `OES M2 LiH 20Q`, run `33576331903`, exact head `3c0a3db4fd9d764ee6df3eabbf968eedac7ae8fd` — PASS.

Full regression on that head: 57 tests — PASS.

At R = 3.015435978393056 bohr (1.5957 Å benchmark geometry):

- OES FCI: -7.988668545192502 Ha
- independent PySCF FCI: -7.988668545192500 Ha
- residual: -2.6645352591e-15 Ha — PASS
- Li-centred `1/r_Li` exposure: 6.075373284338921
- H-centred `1/r_H` exposure: 2.2064867027841704
- exposure difference: +3.868886581554751
- static electronic+nuclear dipole magnitude: 5.786636423819385 D
- NIST v=0 dipole benchmark: 5.8820 D
- diagnostic residual: -0.095363576180615 D (~-1.62%)

The nonzero centre-exposure difference and permanent dipole provide the heteronuclear polarization gate directly from the many-electron state. The dipole comparison is diagnostic because the stored experimental value is a v=0 molecular constant while the solver value is static at the benchmark geometry.

## Born-Oppenheimer local curve

Fixed-20Q energies were computed at R/bohr:

`2.82, 2.92, 3.02, 3.12, 3.22`.

The energies are:

`-7.986635998462766, -7.988072094140859, -7.988680318453642, -7.988598863570611, -7.987945025936497` Ha.

A local quadratic fit gives:

- fitted equilibrium distance: 3.065207030497773 bohr = 1.622037706524885 Å
- NIST `r_e`: 1.5957 Å
- residual: +0.026337706524885 Å (~+1.65%)
- fitted curvature: 0.069564940288866 Ha/bohr^2
- harmonic wavenumber: 1444.649913569344 cm^-1
- NIST `omega_e`: 1405.65 cm^-1
- residual: +38.999913569344 cm^-1 (~+2.77%)

Using the independently computed 7LiH reduced nuclear mass and the predicted equilibrium distance,

B_e = 1 / (2 mu R_e^2)

in atomic units gives:

- predicted rotational constant: 7.274461760325991 cm^-1
- NIST `B_e`: 7.5131 cm^-1
- residual: -0.238638239674009 cm^-1 (~-3.18%)

No NIST value is used to fit the active-space Hamiltonian or the local potential curve.

## Verdict

`M2_IMPLEMENTATION_PASS` — the OES four-electron M_S=0 determinant-subspace Hamiltonian reproduces independent PySCF FCI at machine precision in the identical ten-orbital active space.

`M2_HETERONUCLEAR_POLARIZATION_PASS` — unequal Li/H inverse-radius exposures and a finite permanent dipole emerge from the many-electron state on the fixed 20Q register.

`M2_SPECTROSCOPIC_BASELINE_PARTIAL_PASS` — the same fixed active space predicts equilibrium geometry, harmonic vibration and rotational constant with residuals of approximately 1.65%, 2.77% and 3.18% against the stored 7LiH benchmarks; the static dipole diagnostic is within approximately 1.62% of the stored v=0 value.

Broader LiH rovibrational spectroscopy, excited electronic states, dissociation convergence, external-bath refinement and physical-QPU execution remain OPEN.
