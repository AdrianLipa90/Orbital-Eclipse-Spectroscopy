# OES-M4 H35Cl Fixed-20Q Validation Receipt

Date: 2026-09-02
Branch: `oes-m4-hcl-20q`
Status: `PASS__H35CL_FROZEN_10E_EIGHT_ACTIVE_ELECTRON_FIXED_20Q_BASELINE`
Backend: `SIMULATED_REFERENCE`

## Scope

M4 tests whether the fixed-width OES molecular architecture that passed HF can be reused on a chemically distinct, heavier polar diatomic without changing the active-register width.

- molecule: `H35Cl`, ground electronic state `X 1Sigma+`
- total electrons: 18
- frozen closed-shell core electrons: 10
- active electrons: 8
- active spatial orbitals: 10
- spin orbitals / Jordan-Wigner register width: 20
- exact fixed `Ms=0` determinant dimension: `C(10,4)^2 = 44,100`
- basis: `cc-pVTZ`

The first five occupied spatial MOs form the frozen core. The active space is admitted as complete numerical-degeneracy blocks under the fixed ten-spatial-orbital budget.

Final active indices at the fitted minimum:

`[5, 6, 7, 8, 9, 10, 11, 12, 14, 15]`

Complete group sizes:

`[1, 1, 2, 1, 1, 2, 2]`

## Experimental benchmark separation

NIST Chemistry WebBook / Huber-Herzberg values are benchmark outputs only. They do not enter RHF seed selection, active-space selection, frozen-core construction, Hamiltonian construction, active-curve bracketing or quadratic-fit parameters.

H35Cl benchmark values used after prediction:

- `r_e = 1.27455 Angstrom`
- `omega_e = 2990.946 cm^-1`
- `B_e = 10.59341 cm^-1`
- `omega_e x_e = 52.8186 cm^-1`

Source:
https://webbook.nist.gov/cgi/cbook.cgi?ID=C7647010&Mask=1000

The baseline PASS/FAIL tolerance was predeclared before the final result:

`abs(relative residual) <= 0.05`

for `r_e`, `omega_e` and `B_e`.

## Exact implementation checks

At the final fitted equilibrium geometry:

- OES sparse energy: `-460.1583282170549 Eh`
- independent PySCF active FCI: `-460.15832821704805 Eh`
- energy residual: `-6.8212102633e-12 Eh`
- active 1-RDM maximum residual: `8.4506571e-7`
- sparse matrix nnz: `8,782,452`
- CSR density: `0.00451584062`

Verdict: active Hamiltonian / sparse solver / active 1-RDM implementation PASS for this 8e/20Q sector.

## Blind geometry and spectroscopy result

The RHF-only broad scan selected the seed `R = 2.30 bohr`. The final active 8e/20Q curve was internally bracketed on

`[2.30, 2.36, 2.42, 2.48, 2.54] bohr`

with one benchmark-blind recentering step.

Active energies (Eh):

- 2.30: `-460.1550040868266`
- 2.36: `-460.1573740162031`
- 2.42: `-460.1582998800350`
- 2.48: `-460.15800669747864`
- 2.54: `-460.1566888105412`

Quadratic local fit:

- fitted curvature: `0.34105869112 Eh/bohr^2`
- `r_e = 2.43955738024 bohr = 1.29095816975 Angstrom`
- `omega_e = 3033.97951395 cm^-1`
- `B_e = 10.3314140751 cm^-1`

Residuals against the held-out NIST outputs:

| observable | OES | NIST | residual | relative residual | gate |
|---|---:|---:|---:|---:|---|
| `r_e` (Angstrom) | 1.29095817 | 1.27455 | +0.01640817 | +1.28737% | PASS |
| `omega_e` (cm^-1) | 3033.97951 | 2990.946 | +43.03351 | +1.43879% | PASS |
| `B_e` (cm^-1) | 10.3314141 | 10.59341 | -0.261996 | -2.47320% | PASS |

All three predeclared 5% baseline gates PASS.

## Density-derived diagnostics

At the fitted minimum:

- Cl exposure `<1/r_Cl> = 64.7926869959`
- H exposure `<1/r_H> = 7.89802778451`
- exposure difference: `56.8946592114`
- permanent dipole magnitude: `1.179253545 D`

The dipole is retained as a diagnostic observable in this receipt; it was not part of the predeclared M4 PASS gate.

## Preserved initial FAIL and correction provenance

The first blind M4 run used a five-point active grid

`[2.18, 2.24, 2.30, 2.36, 2.42] bohr`.

Its fitted minimum was `2.41654 bohr`, essentially on the upper sampling boundary. That run produced:

- `r_e` residual about `+0.33%` — PASS
- `B_e` residual about `-0.61%` — PASS
- `omega_e` residual about `+19.96%` — FAIL

This FAIL was not removed or reclassified. Inspection showed that the local minimum was not internally bracketed, so the fitted second derivative was being inferred without sampled active-space energies on both sides of the minimum.

The numerical procedure was therefore corrected without changing the physics model, active-register width, benchmark values or 5% threshold. The replacement routine:

1. evaluates a five-point active grid around the benchmark-blind RHF seed;
2. if the discrete minimum is an endpoint, recenters on that computed point;
3. fits only after the discrete minimum is interior;
4. accepts the quadratic minimum only if it lies between the second and fourth sampled points;
5. otherwise recenters from computed energies and repeats;
6. fails closed if internal bracketing cannot be obtained.

A synthetic parabola unit test explicitly requires recentering and exact recovery of the known minimum / curvature.

The corrected hosted run changed the inferred curvature from approximately

`0.47698 -> 0.34106 Eh/bohr^2`

and the harmonic result from approximately

`3587.96 -> 3033.98 cm^-1`.

Verdict: the original harmonic FAIL was a numerical bracketing failure of the local curvature estimator, not evidence for a failure of the 8e/20Q active Hamiltonian.

## Independent active-subspace continuity diagnostic

Supplemental branch: `oes-m4-hcl-subspace-continuity`
Hosted diagnostic head: `bc78558e0386e7e49ff1bf695be8a01999ac1c9c`

Across adjacent geometries in the final active window `2.30-2.54 bohr`, the selected active indices and complete group sizes remained unchanged. Cross-geometry AO-overlap principal-cosine analysis gave

`minimum adjacent principal cosine = 0.9990927083672345`.

Thus no orbital replacement / active-subspace jump was detected over the final curvature window. This is a supplemental implementation diagnostic, not an experimental benchmark.

## Current M4 verdict

- frozen 10e + active 8e + fixed 20Q construction — PASS
- 44,100-dimensional exact fixed-Ms sparse sector — PASS
- OES vs independent active FCI energy — PASS
- active 1-RDM cross-check — PASS
- blind `r_e` 5% gate — PASS
- blind `omega_e` 5% gate — PASS
- blind `B_e` 5% gate — PASS
- active-subspace continuity diagnostic — PASS
- broader HCl spectroscopy / excited electronic states / anharmonic constants — OPEN
- physical quantum backend execution — OPEN

The promoted simulated-reference M4 baseline is:

`H35Cl : 18e -> frozen 10e + active 8e in one fixed 20Q register`,

with benchmark-blind internally bracketed active-space curvature and common OES sparse / density validation machinery.
