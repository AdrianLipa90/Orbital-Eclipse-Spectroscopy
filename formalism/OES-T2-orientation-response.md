# OES-T2 — Relative-orientation response and subspace tracking

Status: `STANDARD_QM_RESPONSE_CONTROL / MOLECULAR_AB_INITIO_SCAN_BINDING_OPEN`

OES-T2 consumes a declared internal molecular orientation coordinate from an
upstream chemistry/state generator. It does not interpret rigid global rotation
of an isolated molecule as a spectral effect.

## 1. Geometry parameter

Let (	heta) denote a relative torsion, bend, or registry coordinate. The
upstream chemistry layer provides

[
H(	heta), qquad E_n(	heta), qquad |Psi_n(	heta)angle.
]

OES owns the transition projection between tracked states.

## 2. Subspace continuity

At neighboring geometries, orbital/state bases can rotate inside the same
physical subspace. For coefficient matrices (C_L,C_R) and cross-geometry
metric (S_{LR}),

[
M=C_L^dagger S_{LR}C_R.
]

The singular values of (M) are the principal cosines and are invariant under
independent unitary basis rotations inside the two subspaces.

For

[
M=USigma V^dagger,
]

OES-T2 uses the metric orthogonal-Procrustes transport

[
Q=VU^dagger,
qquad
C_R^{m aligned}=C_RQ.
]

This changes only the representation inside the right subspace. It does not
change the subspace itself.

The existing HCl continuity diagnostic is recovered as the principal-cosine
part of this more general transport contract.

## 3. Transition observables at fixed orientation

For each tracked state pair,

[
T^{FI}_{pq}(	heta)
=
langlePsi_F(	heta)|a_p^dagger a_q|Psi_I(	heta)angle.
]

The electric-dipole line remains

[
f_{FI}(	heta)
=
rac{2}{3}
Delta E_{FI}(	heta)
|mu_{FI}(	heta)|^2
]

in atomic units, with

[
mu_{FI}(	heta)
=
langle F(	heta)|hatmu|I(	heta)angle.
]

OES-T2 records the angle, energy gap, frequency, wavelength, transition dipole,
and oscillator strength. Resonant-Chemistry does not own these projections.

## 4. Why a line changes with orientation

For a differentiable tracked branch,

[
rac{df}{d	heta}
=
rac{2}{3}
left[
rac{dDelta E}{d	heta}|mu|^2
+
2Delta E,
operatorname{Re}
left(
mu^daggerrac{dmu}{d	heta}
ight)
ight].
]

OES-T2 exposes the two terms separately:

1. energy-gap contribution;
2. transition-moment/geometry contribution.

Therefore a line can strongly brighten or darken even when its central energy
changes only weakly.

## 5. Control before molecular claim

The first exact control uses a synthetic transition moment proportional to
(cos	heta), giving

[
f(	heta)proptocos^2	heta.
]

This is a regression control, not a universal molecular law.

The next physical gate is a rigid molecular orientation scan using a standard
electronic-structure backend. State identity must be tracked by overlap and
subspace continuity rather than root number alone.

## 6. Phase microscope

Once the molecular backend supplies a common mapped orbital representation,
the same tracked transition may be projected as

[
ho_{FI}(mathbf r;	heta)
=
sum_{pq}T^{FI}_{pq}(	heta)
phi_p^*(mathbf r;	heta)phi_q(mathbf r;	heta).
]

That phase-microscope extension remains downstream of the current T2 control.
