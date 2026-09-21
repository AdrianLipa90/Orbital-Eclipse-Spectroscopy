# OES-T3 — Rigid molecular TDA orientation backend

Status: `STANDARD_QM_PYSCF_RHF_TDA / RAW_ROOT_SCAN / CROSS_GEOMETRY_STATE_TRACKING_REQUIRED`

OES-T3 is the first non-analytic molecular orientation backend in the
RC-to-OES chain. It consumes Cartesian molecular geometry and computes
closed-shell RHF/TDA excited-state observables using the pinned PySCF q1
backend.

## 1. Scope

For one fixed geometry (R),

[
R
ightarrow
H_{m RHF}(R)
ightarrow
{E_n^{m TDA}(R),mu_{0n}(R),f_{0n}(R)}.
]

The implementation supports v0.1 closed-shell (S=0) geometries only.

No experimental transition energy or oscillator strength enters the solver.

## 2. Length-gauge cross-check

For every TDA root, OES independently reconstructs

[
f_{0n}^{m OES}
=
rac{2}{3}
Delta E_{0n}
|mu_{0n}|^2
]

in atomic units and compares it with the length-gauge oscillator strength
returned by PySCF.

The gate fails closed if the two disagree beyond the declared numerical
tolerance.

This makes the backend an implementation cross-check rather than a blind
pass-through of one external scalar.

## 3. Global versus internal rotation

A rigid global rotation of an isolated field-free molecule is a coordinate
change. Scalar observables such as excitation energies, transition-dipole norm,
and oscillator strength must remain unchanged.

The T3 control validates this with the same H2 geometry oriented along two
different laboratory axes.

An internal geometry change is different. The rigid H2O bend control keeps both
O-H bond lengths fixed while changing only the H-O-H angle. The resulting raw
TDA spectrum is required to remain finite and is expected to differ between the
two internal geometries.

## 4. Root identity firewall

For an orientation series,

[
	heta_k
mapsto
{E_n(	heta_k),mu_n(	heta_k),f_n(	heta_k)},
]

T3 returns raw roots independently at each geometry.

It explicitly labels the scan

[
	exttt{RAW\_ROOT\_INDEX\_UNTRACKED}.
]

Root number is not accepted as physical state identity across geometry.

Before computing (dE/d	heta), (df/d	heta), or following a named line,
the states must be transported using overlap/subspace continuity from OES-T2.

## 5. RC/OES boundary

Resonant-Chemistry owns the upstream molecular geometry/orientation coordinate.

OES-T3 owns the standard-QM spectral projection.

The intended chain is

[
	ext{RC rigid geometry}
ightarrow
	ext{OES RHF/TDA}
ightarrow
	ext{OES-T2 state tracking}
ightarrow
	ext{tracked line response}.
]

Phase-microscope reconstruction remains downstream of a tracked transition
representation.
