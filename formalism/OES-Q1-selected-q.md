# OES-Q1 Fixed-20Q Core–Bath Closure

## Scope

This note records the tested many-electron architecture emerging from the helium Q1 benchmark.

The active quantum register remains fixed at ten spatial / twenty spin orbitals:

\[
P_{20Q}:\qquad N_{\mathrm{spin}}=20,\qquad N_e=2,\qquad \dim \mathcal H_{P,N=2}=\binom{20}{2}=190.
\]

A larger source representation is partitioned as

\[
\mathcal H=P\oplus Q,
\]

where `P` is the fixed symmetry/flavor core and `Q` is the external correlation/color bath.

For the geometric d-aug helium source used in Q1,

\[
N_{\mathrm{source,spin}}=124,
\qquad
\dim \mathcal H_{N=2}=\binom{124}{2}=7626,
\]

with

\[
\dim P=190,
\qquad
\dim Q=7436
\]

at the determinant level.

## Flavor-preserving active-space rule

The helium capacity audit supports allocation by complete angular blocks:

\[
N_{\mathrm{spatial}}=\sum_l n_l(2l+1).
\]

The ten-orbital helium capacity oracle selected the group pattern

\[
1+1+3+1+1+3=10,
\]

corresponding to four scalar blocks and two complete three-dimensional vector blocks. The predictive `s4+p6` active space implements the same dimensional pattern as four scalar radial modes plus two complete p-like triplets.

This block rule is the current candidate for extension to larger atoms and symmetry-adapted molecular irreducible representations.

## Color and flavor diagnostics

For a transition class `I`, define the source-energy residual of the active core

\[
r_I^{(P)}=\Delta E_I^{(P)}-\Delta E_I^{(\mathrm{source})}.
\]

The predictive `s4+p6` helium core produced approximately

\[
(r_{^3S}^{(P)},r_{^1S}^{(P)},r_{^1P}^{(P)})
=(-0.67800,-0.64181,-0.66672)\ \mathrm{eV},
\]

with source RMS 0.66235 eV and centered RMS 0.01512 eV. The complete bright manifold and its summed oscillator strength were retained.

The benchmark therefore tracks independent diagnostics

\[
F_{\mathrm{color}}\sim\{\Delta E_I\},
\]

and

\[
F_{\mathrm{flavor}}\sim
\{\text{multiplet dimension},\ f_{\mathrm{manifold}},\ \Delta E_I-\Delta E_J,\ \rho^{(1)}_{FI}\}.
\]

## Complete active-state classes

Bath construction uses complete degenerate classes from the fixed 20Q core:

\[
\mathcal C_P=
\{\,^1S_0^{\mathrm{ground}}\,\}
\oplus
\{\,^3S_1,M_S=-1,0,+1\,\}
\oplus
\{\,^1S_0^{2s}\,\}
\oplus
\{\,^1P_1,x,y,z\,\}.
\]

Using full classes removes arbitrary representative-state choices inside degenerate P-space manifolds.

## External coupling-response states

For each normalized active state

\[
|\Psi_{c,i}^{(P)}\rangle,
\]

define the external response

\[
|\chi_{c,i}\rangle=QH|\Psi_{c,i}^{(P)}\rangle.
\]

Each nonzero response is normalized and mapped to its spin-summed one-body density. Restricting that density to external spatial orbitals gives

\[
\gamma_{c,i}^{(Q)}.
\]

Within a degenerate class of dimension `d_c`, form

\[
\bar\gamma_c^{(Q)}
=
\frac1{d_c}
\sum_{i=1}^{d_c}
\gamma_{c,i}^{(Q)}.
\]

Each class density is trace-normalized, then complete classes receive equal weight:

\[
\Gamma_Q
=
\frac1{|\mathcal C_P|}
\sum_{c\in\mathcal C_P}
\frac{\bar\gamma_c^{(Q)}}{\operatorname{Tr}\bar\gamma_c^{(Q)}}.
\]

The canonical external natural orbitals are eigenvectors of

\[
\boxed{\Gamma_Q u_k=n_k u_k}.
\]

## Rotation covariance

For an arbitrary orthogonal change of the external one-particle basis

\[
Q\to QR,
\qquad R^TR=I,
\]

the response density transforms as

\[
\Gamma_Q\to R^T\Gamma_QR.
\]

Therefore its eigenspaces define the same physical subspaces. Numerically degenerate occupation eigenvalues are retained as complete blocks rather than selecting individual eigenvectors.

The helium hosted gauge test deliberately applied a random orthogonal rotation to all 52 external spatial source orbitals. At the ten-orbital bath prefix:

\[
\cos\theta_{\min}=0.9999999999999991,
\]

\[
\|P_{\mathrm{bath}}-P'_{\mathrm{bath}}\|_F=0,
\]

with maximum principal-cosine error

\[
7.99\times10^{-15}.
\]

Thus the selected bath is invariant at numerical precision under the tested Q-basis rotation.

## Bath block spectrum

The leading external response occupation groups have dimensions

\[
1,3,1,5,3,5,3,7,\ldots
\]

and are admitted whole. This gives a general core–bath allocation rule

\[
\boxed{
\mathcal H_{\mathrm{working}}
=
P_{20Q}
\oplus
\bigoplus_{g\in\mathcal G_K}Q_g
}
\]

where each `Q_g` is a complete natural-orbital occupation eigenspace.

## Hosted natural-bath convergence

The rotation-covariant bath produced:

| external natural orbitals | total spatial orbitals | source RMS (eV) | NIST RMS (eV) | bright spread (eV) | triplet spread (eV) |
|---:|---:|---:|---:|---:|---:|
| 4  | 14 | 0.24474 | 0.22948 | 2.64e-11 | 4.97e-14 |
| 10 | 20 | 0.21509 | 0.20109 | 7.60e-11 | 2.17e-13 |
| 13 | 23 | 0.02862 | 0.05566 | 7.62e-11 | 1.81e-13 |

At the 13-external-orbital closure:

\[
^3S:19.76843\ \mathrm{eV},
\qquad
^1S:20.57865\ \mathrm{eV},
\qquad
^1P:21.29079\ \mathrm{eV}.
\]

Residuals against the full d-aug source are

\[
(-0.03439,-0.03462,-0.00870)\ \mathrm{eV},
\]

with source RMS

\[
0.02862\ \mathrm{eV}.
\]

The complete bright and triplet manifolds remain closed to numerical precision. The bright oscillator-strength sum is 0.31576, compared with approximately 0.32933 in the full d-aug source.

## Determinant-selected diagnostic

A determinant-ranked selected-Q construction independently demonstrated rapid recovery of the missing color energy and reached approximately 0.049-0.050 eV NIST RMS near 512 selected external determinants. Repeated hosted executions exposed sensitivity of the exact selected determinant list to rotations inside degenerate external orbital subspaces. It remains a useful convergence diagnostic and motivates selection at the physical orbital-subspace level.

## Current architecture

The tested helium architecture is

\[
\boxed{
\text{fixed 20Q symmetry/flavor core}
+
\text{complete-class rotation-covariant external natural-orbital bath}
}.
\]

At the tested 13-orbital bath closure it reaches the accuracy scale of the full d-aug source while preserving triplet and bright multiplet structure to numerical precision. Broader He I spectroscopy and physical-backend execution remain the next validation layers.
