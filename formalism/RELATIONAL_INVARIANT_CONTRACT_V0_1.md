# Relational Invariant Contract v0.1

Status: `EXACT_REPRESENTATION_INVARIANCE_CONTRACT / PHYSICAL_EXTENSION_FIREWALL`

Contract ID: `RELATIONAL_INVARIANT_CONTRACT_V0_1`

## 1. Core rule

Let a representation space (mathcal R) carry an action of a redundancy/gauge group (G). A physical observable (O) is admissible only if it factors through the quotient:

[
pi:mathcal R	omathcal R/G,
qquad
O=ar Ocircpi,
]

equivalently

[
oxed{O(g!cdot!r)=O(r)quadorall gin G.}
]

Intermediate objects may be covariant rather than invariant, but any claimed observable must remove the arbitrary representation dependence or carry an explicit external phase/frame reference.

This contract is shared by Resonant Chemistry (RC) and Orbital Eclipse Spectroscopy (OES).

## 2. RC SCF occupied-subspace profile

For an AO overlap metric (S>0) and occupied coefficient matrix (C), the physically relevant occupied subspace is unchanged by

[
Clongrightarrow C U,
]

for any unitary rotation (U) inside the occupied subspace.

Define the (S)-metric projector

[
oxed{
P_S(C)=C(C^dagger S C)^{-1}C^dagger S.
}
]

Then

[
oxed{P_S(CU)=P_S(C).}
]

Therefore column identity, ordering and sign/phase are representation choices. SCF continuation must compare occupied subspaces, not raw columns.

The (S)-metric orthogonal Procrustes step is a gauge-fixing/alignment operation:

[
U_*=argmin_{U^dagger U=I}|C_{m cand}U-C_{m old}|_S.
]

It is not a new physical interaction and must not change subspace observables.

## 3. OES orbital-gauge profile

For a diagonal orbital rephasing

[
D_chi=operatorname{diag}(e^{ichi_1},ldots,e^{ichi_n}),
qquad
philongrightarrowphi D_chi,
]

the OES transition convention

[
T^{FI}_{pq}=langlePsi_F|a_p^dagger a_q|Psi_Iangle
]

uses the matched covariance rule

[
oxed{T^{FI}longrightarrow D_chi T^{FI}D_chi^dagger.}
]

Hence

[
ho_{FI}(mathbf r)
=
sum_{pq}T^{FI}_{pq}phi_p^*(mathbf r)phi_q(mathbf r)
]

is exactly invariant under matched orbital rephasing.

This is the executable OES profile
`OES_ORBITAL_U1_REPHASING_V0_1`.

## 4. State-ray phase firewall

The phases of the many-electron states are also arbitrary:

[
|Psi_Iangle	o e^{ialpha}|Psi_Iangle,
qquad
|Psi_Fangle	o e^{ieta}|Psi_Fangle.
]

Therefore

[
T^{FI}	o e^{i(alpha-eta)}T^{FI},
qquad
ho_{FI}(mathbf r)	o e^{i(alpha-eta)}ho_{FI}(mathbf r).
]

Consequences:

- (|ho_{FI}(mathbf r)|) is invariant;
- nodes/zero sets are invariant;
- phase differences (argho(mathbf r)-argho(mathbf r_0)) are invariant where defined;
- (|mu_{FI}|^2) and oscillator strengths are invariant;
- the absolute raw phase of (ho_{FI}) is not an observable unless an external phase reference fixes the global (U(1)).

Accordingly, Phase Microscope may expose raw phase as a representation-level diagnostic, but physical phase claims must use a declared reference or relative phase.

## 5. RC graph-holonomy profile

For an oriented graph edge phase (	heta_{ij}),

[
	heta_{ij}	o	heta_{ij}+chi_j-chi_i.
]

For a closed cycle (C),

[
oxed{
W_C=exp!left(isum_{ein C}s_e	heta_eight)
}
]

is invariant.

With the RC Hamiltonian convention used in code,

[
H_{ij}=t_{ij}e^{i	heta_{ij}},
qquad
U_chi=operatorname{diag}(e^{ichi_i}),
]

the transformed Hamiltonian obeys

[
oxed{
H(	heta+B^Tchi)=U_chi^dagger H(	heta)U_chi.
}
]

Therefore eigenvalues and line positions are invariant under vertex rephasing; only cycle-gauge data can affect the spectrum inside the declared graph Hamiltonian.

This fixes the unitary-conjugation orientation for the current incidence convention.

## 6. Cross-repository boundary

RC owns state generation. OES owns transition and spectroscopy inference.

Every RC-to-OES transition packet that opts into this contract must declare:

- `schema = RELATIONAL_INVARIANT_CONTRACT_V0_1`;
- `profile = OES_ORBITAL_UNITARY_BASIS_V0_1`;
- the exact orbital basis identifier;
- the transition-RDM convention;
- provenance for source repository, commit, method and backend.

OES must fail closed on a declared but unsupported contract/profile.

Legacy packets without this additive block remain readable under the existing `OES_TRANSITION_STATE_V0_1` compatibility path.

## 7. Validation obligations

A representation-sensitive algorithm is admissible only when its tests include the corresponding gauge transformation and verify the invariant object:

1. SCF: rotate an occupied subspace and verify (P_S) is unchanged.
2. OES orbital gauge: rephase orbitals and (T^{FI}) coherently and verify (ho_{FI}) is unchanged.
3. OES state-ray gauge: apply a global transition phase and verify amplitude, nodes, relative phase and intensity are unchanged.
4. Graph holonomy: vertex-rephase edge phases and verify Wilson loops and Hamiltonian spectra are unchanged.

Gauge fixing may improve numerical continuity, but a PASS cannot depend on a particular arbitrary gauge.

## 8. Physical-status firewall

This contract formalizes representation invariance already required by standard quantum mechanics and numerical linear algebra.

It does **not** prove:

- (W_{m chem}=W_{m sem});
- a new chemical interaction beyond standard QM;
- that a raw transition-density phase is absolutely observable;
- that numerical gauge fixing creates physical dynamics.

Those remain separate physical hypotheses and validation gates.
