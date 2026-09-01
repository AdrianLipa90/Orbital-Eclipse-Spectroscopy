# OES-H0 — Hydrogen Orbital-Eclipse Closure

## Status

Candidate formalism under executable validation. The gross hydrogen solver takes physical constants and quantum numbers as inputs; empirical spectral data are reserved for benchmark fixtures.

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
- signed radial overlap and cancellation coherence.

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

These coordinates remain separately available in the fingerprint.

## 5. Radial eclipse coherence

With `u=r/a_H` and normalized dimensionless radial functions `R_tilde_nl(u)`, define

\[
g_{if}(u)=\widetilde R_i(u)\widetilde R_f(u)u^3.
\]

The signed and unsigned E1 radial overlaps are

\[
R_{if}=\int_0^\infty g_{if}(u)du,
\qquad
A_{if}=\int_0^\infty |g_{if}(u)|du,
\]

and the radial cancellation coherence is

\[
\mathcal C_{if}=\frac{|R_{if}|}{A_{if}}.
\]

`C=1` marks complete radial sign coherence; decreasing values quantify stronger cancellation among radial regions.

## 6. Relativistic flavor splitting

The reference relativistic layer evaluates the hydrogenic Dirac energy by `(n,j)` using the ordinary-hydrogen reduced mass as the compact recoil approximation. The `2p_1/2` / `2p_3/2` split is benchmarked against a locked NIST ASD target.

The pure Coulomb-Dirac `(n,j)` degeneracy keeps `2s_1/2` and `2p_1/2` on the same reference energy. Their measured separation defines the Lamb/QED benchmark gate.

## 7. Leading QED Lamb layer

For the `n=2` Lamb interval, OES-H0 currently evaluates a controlled leading low-Z approximation with

\[
\Delta E_{SE}=\frac{\alpha}{\pi}\frac{\alpha^4}{n^3}
\left(\frac{m_r}{m_e}\right)^3m_ec^2
\left[A_{41}L+A_{40}\right],
\]

where

\[
L=\ln\left[\frac{m_e}{m_r}\alpha^{-2}\right],
\]

and the state dependence includes the Bethe logarithms

\[
\ln k_0(2S)=2.811769893,
\qquad
\ln k_0(2P)=-0.030016709.
\]

The leading Uehling vacuum-polarization term uses

\[
V_{40}=-\frac{4}{15}\delta_{l0}.
\]

For ordinary hydrogen this layer gives a `2S_1/2-2P_1/2` interval of about `1050.55 MHz`, compared with the locked NIST target near `1057.85 MHz`. The residual is retained for higher-order self-energy, recoil, two-loop, finite-size and related QED contributions.

This gate shows that central/contact exposure is one ingredient of the Lamb structure while the Bethe logarithm contributes a state-wide spectral component.

## 8. Flavor-to-color conversion under a magnetic field

In the weak-field linear Zeeman regime,

\[
\Delta\nu_Z=\frac{\mu_B}{h}B\left(g_um_u-g_lm_l\right).
\]

For `p` states the continuous crossover uses

\[
\frac{H}{h}=A\frac{\mathbf L\cdot\mathbf S}{\hbar^2}
+\frac{\mu_B}{h}B(L_z+2S_z)/\hbar.
\]

At weak field the eigenbranches recover the `|j,m_j>` Landé slopes. At strong field they approach the uncoupled `|m_l,m_s>` Paschen–Back slopes.

## 9. Validation boundary

H0 currently separates:

- analytic gross-color reconstruction — PASS on hosted CI,
- orbital E1 flavor gates — PASS on hosted CI,
- central exposure kernels — PASS on hosted CI,
- reduced-mass Dirac fine-structure reference — PASS against the declared NIST tolerance,
- linear Zeeman and p-state Paschen–Back continuation — PASS on hosted CI,
- signed radial cancellation/coherence — PASS on hosted CI,
- leading one-loop Lamb/QED reference — implemented; exact-head CI pending,
- complete QED closure — OPEN,
- broad empirical NIST held-out benchmark ledger — OPEN.

OPEN items retain their independent validation status until their own gates are resolved.
