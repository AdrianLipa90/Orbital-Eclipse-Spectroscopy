# OES-M2 LiH Validation Receipt

Date: 2026-09-02
Branch: `oes-m2-lih-20q`
Parent: `oes-m1-h2-20q`
Backend: `SIMULATED_REFERENCE`
Status: `PARTIAL_PASS__LIH_ROVIBRATIONAL_FIXED_20Q_BASELINE`

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

## M2.1 hosted rovibrational spectrum gate

Validated implementation head: `afd87f0b23ffe6e370d84e4cf71df7a07cb53e31`.
Hosted workflow: `OES M2 LiH 20Q`, run `33576868980`, job `100082656678` — PASS.
Full exact-head regression: 59/59 PASS.
The pre-existing M2 fixed-20Q LiH runner also passed unchanged on the same head.

The electronic Born-Oppenheimer curve was sampled at 18 geometries from 2.0 to 5.4 bohr using the same `cc-pVTZ`, four-electron, fixed-20Q active-space protocol. The nuclear solver then used only this sampled OES potential and the independently computed 7LiH reduced nuclear mass. It applied a PCHIP interpolation inside the sampled interval and a finite-difference radial Hamiltonian including the J(J+1)/(2 mu R^2) centrifugal term for J=0,1,2. Experimental constants were not solver inputs.

The requested v=0..3 states remain strongly confined inside the sampled potential interval; the smallest reported endpoint margin exceeds 7078 cm^-1.

### Predicted term values

For J=0:

- v=0: 0.000000 cm^-1
- v=1: 1335.003803 cm^-1
- v=2: 2636.911126 cm^-1
- v=3: 3888.990176 cm^-1

For J=1:

- v=0: 14.419103 cm^-1
- v=1: 1349.029100 cm^-1
- v=2: 2650.499637 cm^-1
- v=3: 3902.215904 cm^-1

For J=2:

- v=0: 43.237747 cm^-1
- v=1: 1377.061116 cm^-1
- v=2: 2677.657587 cm^-1
- v=3: 3928.648907 cm^-1

### Vibrational constants recovered from OES levels

- fundamental v=0→1: 1335.003803 cm^-1
- independent NIST-constant-derived benchmark: 1359.779750 cm^-1
- residual: -24.775947 cm^-1 (~-1.82%)
- omega_e: 1352.065647 cm^-1 vs NIST 1405.65 cm^-1; residual -53.584353 cm^-1 (~-3.81%)
- omega_e x_e: 3.999395 cm^-1 vs NIST 23.20 cm^-1; residual -19.200605 cm^-1 (~-82.76%)
- omega_e y_e: -2.788632 cm^-1 vs NIST +0.163 cm^-1; higher-order anharmonicity is not reproduced by the present fixed-20Q potential.

### Rotational and vibration-rotation constants recovered from OES levels

- B_v(v=0,1,2) = 7.211182, 7.014197, 6.795845 cm^-1
- B_e = 7.301662 cm^-1 vs NIST 7.5131 cm^-1; residual -0.211438 cm^-1 (~-2.81%)
- alpha_e = 0.175619 cm^-1 vs NIST 0.2132 cm^-1; residual -0.037581 cm^-1 (~-17.63%)
- gamma_e = -0.010683 cm^-1 vs NIST +0.00075 cm^-1; higher-order vibration-rotation curvature is not reproduced by the present baseline
- D_v(v=0) = 0.000815104 cm^-1 vs NIST D_e = 0.0008617 cm^-1; residual -0.000046596 cm^-1 (~-5.41%)

The M2.1 result therefore separates two regimes cleanly: the low-order rovibrational observables (fundamental, B_e and D_v0) are already close to the independent experimental benchmark, while the higher-order anharmonic and vibration-rotation curvature terms expose the remaining fixed-active-space / electronic-potential error budget.

## Verdict

`M2_IMPLEMENTATION_PASS` — the OES four-electron M_S=0 determinant-subspace Hamiltonian reproduces independent PySCF FCI at machine precision in the identical ten-orbital active space.

`M2_HETERONUCLEAR_POLARIZATION_PASS` — unequal Li/H inverse-radius exposures and a finite permanent dipole emerge from the many-electron state on the fixed 20Q register.

`M2_SPECTROSCOPIC_BASELINE_PARTIAL_PASS` — equilibrium geometry, local harmonic vibration and rotational scale are recovered at the few-percent level from the same fixed active space.

`M2_1_ROVIBRATIONAL_SOLVER_PASS` — the sampled OES Born-Oppenheimer potential supports stable, confined numerical v,J eigenstates and independently recovered term values without experimental fitting.

`M2_1_LOW_ORDER_SPECTROSCOPY_PARTIAL_PASS` — the fundamental, B_e and D_v0 are respectively within approximately 1.82%, 2.81% and 5.41% of the stored independent 7LiH benchmarks.

`M2_1_HIGHER_ORDER_REFINEMENT_OPEN` — omega_e x_e, omega_e y_e, alpha_e and gamma_e expose unresolved higher-order potential-shape and vibration-rotation curvature errors. These remain targets for active-space/external-bath/basis refinement rather than post-hoc fitting.

Excited electronic states, transition-intensity/dipole-curve spectroscopy, dissociation convergence, external-bath refinement and physical-QPU execution remain OPEN.
