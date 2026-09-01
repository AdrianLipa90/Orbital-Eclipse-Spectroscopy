# OES-H0 Validation Receipt

Date: 2026-09-01
Branch: `oes-h0-hydrogen-closure-20260901`
Status: PARTIAL_PASS__H0_OPEN

## Hosted validation already observed

- Head `dcbe7536d67a481a28fde2ffca53fd9e61776e1d`: `OES reference suite` SUCCESS, run `33527502380`.
- Head `f1cfc1b808e03857dc500c10562a8591790fdb2e`: `OES reference suite` SUCCESS, run `33527782318`.
- Current radial-extension head requires its own exact-head CI before promotion.

## Implemented checks

1. Reduced-mass Bohr radius.
2. Hydrogen ground-state gross energy from constants.
3. Lyman-alpha gross color from `E_2-E_1` without spectral wavelength input.
4. Balmer-alpha gross color.
5. Independent gross-color and E1-flavor gates (`2p->1s` open, `2s->1s` E1 closed).
6. Contact exposure `delta_l0/n^3`.
7. Radial-node code `N_r=n-l-1`.
8. Orbital flavor count `2*n_f-1`.
9. Orbital m-resolved flavor count `3[n_f^2+(n_f-1)^2]`.
10. Lande factors for `s_1/2`, `p_1/2`, `p_3/2` in the LS approximation.
11. Linear Zeeman flavor-to-color shift.
12. Reduced-mass Dirac bound-energy reference solver.
13. NIST H I 2p fine-structure held-out benchmark; no benchmark value is a solver input.
14. Continuous p-state Zeeman -> Paschen-Back diagonalization.
15. Hydrogenic normalized radial functions.
16. Signed E1 radial overlap, absolute overlap and cancellation coherence.
17. Exact `2p->1s` radial dipole integral gate.
18. Distinct H-alpha radial-flavor gates for `3s->2p`, `3p->2s`, `3d->2p`.

## Gates

- H0.1 Gross spectrum: PASS.
- H0.2 Orbital flavor geometry: IMPLEMENTED; current exact-head CI pending after radial extension.
- H0.3 Fine structure + Lamb/contact: DIRAC REFERENCE + CONTACT PASS; FULL LAMB/QED OPEN.
- H0.4 Zeeman -> Paschen-Back: IMPLEMENTED; current exact-head CI pending after radial extension.
- Empirical benchmark: PARTIAL PASS (NIST 2p fine split); broad held-out level/line suite OPEN.

## Fine-structure benchmark

Reduced-mass Dirac prediction for the H I `2p_1/2`/`2p_3/2` split is approximately `10.943683 GHz`. The locked NIST ASD target is `10.969051 GHz`, a relative discrepancy of about `0.23%`. The residual is not silently absorbed: recoil/QED refinements remain outside this reference Dirac gate.

## Promotion rule

Full OES-H0 closure remains OPEN until the exact current head passes hosted CI and the declared full Lamb/QED gate plus broader empirical benchmark are separately resolved. Passing analytic/Dirac/radial gates does not promote the OPEN QED claim.
