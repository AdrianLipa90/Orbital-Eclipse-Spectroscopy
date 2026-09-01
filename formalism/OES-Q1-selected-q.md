# OES-Q1 Selected-Q Core–Bath Closure

## Scope

This note records the tested many-electron architecture emerging from the helium Q1 benchmark.

The active quantum register remains fixed at ten spatial / twenty spin orbitals:

\[
P_{20Q}:\qquad N_{\mathrm{spin}}=20,\qquad N_e=2,\qquad \dim \mathcal H_{P,N=2}=\binom{20}{2}=190.
\]

A larger source representation is partitioned as

\[
\mathcal H = P\oplus Q,
\]

where `P` is the fixed 20Q flavor-preserving core and `Q` contains determinants with at least one external spin orbital.

For the geometric d-aug helium source used in Q1,

\[
N_{\mathrm{source,spin}}=124,
\qquad
\dim \mathcal H_{N=2}=\binom{124}{2}=7626,
\]

and therefore

\[
\dim Q=7626-190=7436.
\]

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

with source RMS 0.66235 eV but centered RMS only 0.01512 eV. The complete bright manifold and its summed oscillator strength were retained.

The benchmark therefore tracks two independent diagnostics:

\[
F_{\mathrm{color}}\sim \{\Delta E_I\},
\]

and

\[
F_{\mathrm{flavor}}\sim
\{\text{multiplet dimension},\ f_{\mathrm{manifold}},\ \Delta E_I-\Delta E_J,\ \rho^{(1)}_{FI}\}.
\]

## Selected external bath

External determinants are ranked without experimental energies. For each active-state class `c`, the raw determinant importance is

\[
I_a^{(c)}=
\frac{1}{d_c}
\sum_{i\in c}
\frac{|\langle a|H|\Psi_i^{(P)}\rangle|^2}
{|E_c^{(P)}-H_{aa}|},
\]

where `d_c` is the number of active states in the class. Each class distribution is normalized independently and the final score is the equal-weight class mean over

\[
\{\mathrm{ground},\ ^3S,\ ^1S,\ ^1P\ \mathrm{manifold}\}.
\]

Numerically tied importance groups are admitted together.

For a selected determinant set `Q_K`, the dressed problem is solved as an exact finite subspace eigenproblem

\[
H_{K}=
\begin{pmatrix}
H_{PP} & H_{PQ_K}\\
H_{Q_KP} & H_{Q_KQ_K}
\end{pmatrix}.
\]

The quantum-core dimension remains 190 throughout; increasing `K` enlarges only the external classical bath used by this simulated-reference diagnostic.

## Hosted convergence result

The state-balanced selected-Q sweep gave:

| selected Q | source RMS (eV) | NIST RMS (eV) | bright f sum | bright spread (eV) |
|---:|---:|---:|---:|---:|
| 32  | 0.43182 | 0.41392 | 0.32616 | 1.36e-9 |
| 64  | 0.36507 | 0.34826 | 0.32426 | 1.36e-9 |
| 128 | 0.24695 | 0.23304 | 0.32335 | 0.00206 |
| 256 | 0.12734 | 0.12168 | 0.32307 | 0.00495 |
| 512 | 0.02483 | 0.04993 | 0.31931 | 0.00503 |

At `K=512`, the three excitation energies were

\[
^3S: 19.78022\ \mathrm{eV},
\qquad
^1S: 20.57933\ \mathrm{eV},
\qquad
^1P: 21.28583\ \mathrm{eV}.
\]

The corresponding residuals against the full d-aug source were

\[
(-0.02259,-0.03394,-0.01366)\ \mathrm{eV}.
\]

Thus the selected bath reduced source RMS from 0.66235 eV for the bare predictive 20Q core to 0.02483 eV while preserving high P-space weights.

## Current closure state

The tested architecture is

\[
\boxed{\text{fixed 20Q flavor core}+\text{Hamiltonian-selected external color bath}}.
\]

The 512-determinant bath reaches the NIST RMS scale of the full d-aug source representation. Exact angular closure of the truncated selected bath remains OPEN: its first bright manifold has a residual 5.03 meV spread. The next symmetry gate is therefore closure of `Q_K` under complete angular/molecular symmetry blocks before extension to additional atoms and molecules.
