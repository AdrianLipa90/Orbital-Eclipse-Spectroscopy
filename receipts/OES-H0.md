# OES-H0 Validation Receipt

Date: 2026-09-01
Branch: `oes-h0-hydrogen-closure-20260901`
Status: OPEN_PENDING_CI

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

## Gates

- H0.1 Gross spectrum: IMPLEMENTED / CI PENDING
- H0.2 Orbital flavor geometry: PARTIAL / CI PENDING
- H0.3 Fine structure + Lamb/contact: CONTACT IMPLEMENTED; FULL RELATIVISTIC/QED OPEN
- H0.4 Zeeman -> Paschen-Back: LINEAR ZEEMAN IMPLEMENTED; PASCHEN-BACK OPEN
- Empirical held-out benchmark: OPEN

## Promotion rule

This receipt becomes PASS only after repository-hosted CI executes the exact branch head successfully. Full OES-H0 closure remains OPEN until the declared OPEN gates are implemented and separately validated.
