# OES-T4 — Transition-density subspace tracking

Status: STANDARD_TDA_STATE_CONTINUITY_CONTROL / BLOCK_IDENTITY_ONLY / INDIVIDUAL_ROOT_IDENTITY_OPEN

OES-T4 replaces raw TDA root-number tracking with a gauge-clean comparison of
transition-density subspaces across neighboring molecular geometries.

## 1. Runtime transition amplitude

For restricted TDA, let X^(n)_ia be the occupied-to-virtual excitation
amplitude for root n. In the AO representation OES-T4 forms

T_n^AO = C_occ X^(n) C_vir^†.

This runtime matrix is used for state continuity only. Its absolute scalar
normalization is not used as a spectroscopic observable.

Oscillator strength remains independently evaluated and cross-checked in T3.

## 2. Cross-geometry metric

For geometries L and R, define the AO cross overlap

S_LR = <chi^(L) | chi^(R)>.

For two AO transition matrices,

<T_L,T_R> = Tr[T_L^† S_LR T_R S_RL].

At one geometry this reduces to the native AO metric with S_LL.

## 3. Block Gram whitening

A nearly degenerate physical sector may be represented by arbitrary mixtures of
individual numerical roots. For transition-density blocks {T_i^(L)} and
{T_j^(R)}, form

(G_L)_ij = <T_i^(L),T_j^(L)>,
(G_R)_ij = <T_i^(R),T_j^(R)>,

and cross matrix

C_ij = <T_i^(L),T_j^(R)>.

The whitened overlap is

C_tilde = G_L^(-1/2) C G_R^(-1/2).

Its singular values are the principal cosines between the two
transition-density spans.

The result is invariant under nonsingular changes of basis inside either
supplied block, including phase changes, rotations, and rescalings.

## 4. Interpretation

For principal cosines sigma_k,

0 <= sigma_k <= 1.

OES-T4 reports the minimum principal cosine and chordal distance

d_chord = sqrt(sum_k(1-sigma_k^2)).

A root permutation or arbitrary mixing inside a stable block does not by
itself imply a physical discontinuity.

## 5. Identity firewall

Root labels such as root 1 or root 2 are selectors only.

OES-T4 currently establishes block continuity, not a unique one-to-one
identity map for individual states inside a degenerate or nearly degenerate
block.

Consequently, derivatives of one named line across a crossing remain OPEN
until the block is further resolved by additional invariant observables.

## 6. Ne...Cl2 control

The first molecular T4 receipt uses the same rigid Ne...Cl2 orientation control
as T3:

r_ClCl = 2.0 Angstrom,
D_Ne-mid = 4.0 Angstrom,

with angles 0, 30, 60, and 90 degrees.

The bright roots 1-2 and dark roots 3-4 are compared as two-dimensional
transition-density spans between adjacent angles.

No experimental spectrum is used.
