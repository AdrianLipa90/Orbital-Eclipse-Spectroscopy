# OES-T6 — Spatial transition-density cancellation

Status: STANDARD_TDA_REAL_SPACE_CONTROL / BLOCK_BASIS_INVARIANT / EXPERIMENTAL_COLOR_BINDING_OPEN

OES-T6 asks where orientation-dependent brightening and darkening arise in
real space without choosing an arbitrary basis inside a nearly degenerate TDA
root block.

## 1. Exact runtime transition density

For restricted singlet TDA, PySCF contracts a one-body operator with the
occupied-virtual amplitude using a factor of two. OES-T6 therefore binds the
runtime AO transition density as

T_n^AO = 2 C_occ X_n C_vir^dagger.

Before any real-space analysis, OES verifies for every root that

Tr[r_alpha T_n^AO] = mu_(n,alpha)

using the same nuclear-charge-center origin as the backend length-gauge
transition dipole.

## 2. Real-space field

On an atom-centered numerical quadrature grid,

rho_n(r) = sum_(mu,nu) phi_mu*(r) T_(n,mu nu)^AO phi_nu(r).

For a tracked block B, define the root-basis-invariant activity field

I_B(r) = sum_(n in B) |rho_n(r)|^2.

Under any unitary mixing of roots inside B, I_B(r) is unchanged.

## 3. Dipole envelope and cancellation

Let r be measured from the molecular nuclear-charge center. Define

A_B = integral |r| sqrt(I_B(r)) d^3r

and the exact block dipole norm

M_B = sqrt(sum_(n in B) |mu_n|^2).

Both are invariant under unitary root mixing.

By the triangle inequality,

0 <= M_B / A_B <= 1

up to numerical quadrature error.

OES-T6 defines

C_B = M_B / A_B

as the dipole survival ratio and

1 - C_B

as the dipole cancellation fraction.

A brightening event can therefore be separated into:

1. increased local transition-density activity or envelope A_B;
2. increased survival C_B, meaning less spatial/vector cancellation;
3. a combination of both.

## 4. Additional spatial descriptors

OES-T6 reports the integral of I_B(r), its activity centroid relative to the
charge center, and its RMS radius.

These are control diagnostics. They are not atomic charges or a population
analysis.

## 5. Numerical firewall

The numerical-grid dipole reconstructed from rho_n(r) must reproduce the
backend transition dipole within declared absolute and relative tolerances.

The AO-matrix contraction is separately required to match the backend dipole
at near machine precision.

## 6. Ne...Cl2 control

The first T6 receipt uses the same rigid Ne...Cl2 scan as T3-T5 with fixed

r_ClCl = 2.0 Angstrom

and

D_Ne-mid = 4.0 Angstrom.

The principal diagnostic compares the weak roots 3-4 sector at 0 degrees and
30 degrees, where T5 found a roughly 550-fold oscillator-strength increase.

T6 determines whether that increase is associated primarily with a change in
the real-space transition-density envelope, the dipole survival ratio, or both.

No experimental spectrum and no visible-color assignment enter this gate.
