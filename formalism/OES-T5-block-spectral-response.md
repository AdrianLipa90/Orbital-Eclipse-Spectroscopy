# OES-T5 — Block-resolved orientation spectroscopy

Status: STANDARD_TDA_BLOCK_RESPONSE_CONTROL / TRACKED_SECTOR_RESPONSE / EXPERIMENTAL_COLOR_BINDING_OPEN

OES-T5 combines T3 spectral observables with T4 transition-density subspace
continuity. Its purpose is to follow a physically continuous sector through an
orientation scan without pretending that numerical root labels are state
identity.

## 1. Block spectral observables

For a selected tracked block B of TDA roots, OES-T5 reports

F_B = sum_(n in B) f_n,

D_B = sum_(n in B) |mu_n|^2,

the energy interval [E_min,E_max], and the arithmetic centroid

E_bar_B = (1/|B|) sum_(n in B) E_n.

When F_B is non-zero, it also reports the oscillator-strength weighted centroid

E_bar_B^(f) = [sum_(n in B) f_n E_n] / F_B.

The fraction of oscillator strength is explicitly relative only to the
currently computed TDA window. It is not a Thomas-Reiche-Kuhn completeness
claim.

## 2. Geometry-to-intensity response

For neighboring orientation samples L and R, the T4 continuity object supplies
the principal cosines of the transition-density spans.

T5 reports, for the same tracked block,

Delta F_B = F_B(R) - F_B(L),

the ratio F_B(R)/F_B(L) when the denominator is non-zero, and centroid shifts.

Thus a large intensity change can be distinguished from a discontinuous change
of the underlying numerical state space.

## 3. Near-degenerate sectors

Inside an exactly or nearly degenerate block, individual roots may rotate,
exchange order, or acquire arbitrary phase. Root labels are therefore
selectors only.

The block-level response is intentionally the primary observable until extra
invariants justify a unique state-by-state continuation.

## 4. Ne...Cl2 control

The first T5 receipt uses the rigid T3/T4 Ne...Cl2 geometry control at fixed

r_ClCl = 2.0 Angstrom

and

D_Ne-mid = 4.0 Angstrom,

with angles 0, 30, 60, and 90 degrees.

The first two TDA roots are treated as one bright two-dimensional sector and
roots 3-4 as a second weak/dark sector. T4 continuity is evaluated between
adjacent angles before interpreting changes in block intensity.

This is a numerical control in a minimal STO-3G basis. No experimental spectrum
is used and no visible-color claim is admitted by this gate.
