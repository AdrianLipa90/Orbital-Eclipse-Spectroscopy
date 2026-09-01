# OES-H0 — Hydrogen Orbital-Eclipse Closure

## Status

Candidate formalism under executable validation. No empirical spectral wavelength is used as an input to the gross hydrogen solver.

## 1. Color coordinate

For ordinary hydrogen with reduced mass `mu`, Coulomb coupling `g=e^2/(4*pi*epsilon_0)`, and principal quantum number `n`:

\[
E_n=-\frac{\mu g^2}{2\hbar^2 n^2}.
\]

For emission `n_i > n_f`:

\[
\omega_{if}=\frac{E_i-E_f}{\hbar},\qquad
\nu_{if}=\frac{E_i-E_f}{h},\qquad
\lambda_{if}=\frac{c}{\nu_{if}}.
\]

This is the OES gross-color coordinate.

## 2. Transition/eclipsing field

The relational transition field is

\[
\tau_{fi}(\mathbf r,t)=\psi_f^*(\mathbf r)\psi_i(\mathbf r)e^{-i\omega_{if}t}.
\]

Its temporal phase carries the gross transition frequency while its spatial structure carries orbital channel information.

## 3. Flavor coordinates

The minimal H0 flavor state records:

- `(l_i,l_f)` orbital route,
- `(m_i,m_f)` and `Delta m`,
- radial node counts `N_r=n-l-1`,
- E1 selection gate,
- central/contact exposure,
- later: radial signed-overlap coherence and full multipole decomposition.

For a gross color terminating on shell `n_f`, the number of distinct orbital E1 routes is

\[
N_{\mathrm{flavor}}=2n_f-1,
\]

and the number of orbital `m`-resolved E1 routes before spin is

\[
N_m=3\left[n_f^2+(n_f-1)^2\right].
\]

## 4. Exposure kernels

OES treats central sensitivity as a family of radial probes:

\[
\mathcal X_K[\rho]=\int \rho(\mathbf r)K(r)d^3r.
\]

Current H0 kernels are

\[
\langle r^{-1}\rangle=\frac{1}{n^2a_H},
\]

\[
\langle r^{-3}\rangle=\frac{1}{a_H^3n^3l(l+1/2)(l+1)}\quad(l>0),
\]

and the contact exposure

\[
\pi a_H^3|\psi_{nl}(0)|^2=\frac{\delta_{l0}}{n^3}.
\]

These coordinates are kept separate rather than collapsed into a single scalar.

## 5. Flavor-to-color conversion under a magnetic field

In the weak-field linear Zeeman regime,

\[
\Delta\nu_Z=\frac{\mu_B}{h}B\left(g_um_u-g_lm_l\right).
\]

Thus an `m_j` flavor coordinate that is degenerate at `B=0` is mapped into a measurable frequency displacement when the degeneracy is lifted.

## 6. Validation boundary

H0 currently separates:

- analytic gross-color reconstruction — implemented,
- orbital E1 flavor gates — implemented,
- central exposure kernels — implemented,
- linear Zeeman flavor-to-color mapping — implemented,
- full relativistic fine-structure solver — OPEN,
- full Lamb/QED solver — OPEN,
- signed radial cancellation/coherence solver — OPEN,
- empirical NIST held-out benchmark ledger — OPEN,
- Paschen–Back continuation — OPEN.

No OPEN item is promoted by the presence of the analytic hydrogen checks.
