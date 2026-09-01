# OES-Q1 — Quantum Helium 20Q

## Scope

OES-Q1 extends the hydrogen color/flavor construction to the first correlated two-electron atom. The canonical active space is fixed to 10 spatial orbitals, represented as 20 spin orbitals and therefore 20 Jordan–Wigner qubits.

The full register has

\[
2^{20}=1,048,576
\]

computational basis states. Exact reference simulation exploits particle-number conservation and diagonalizes only the two-electron sector

\[
\binom{20}{2}=190.
\]

This symmetry reduction does not alter the 20-mode fermionic encoding.

## Many-body Hamiltonian

OES-Q1 evaluates

\[
H=\sum_{pq}h_{pq}a_p^\dagger a_q
+\frac14\sum_{pqrs}\langle pq\Vert rs\rangle
 a_p^\dagger a_q^\dagger a_s a_r+E_{nuc}.
\]

The one- and two-electron integrals are generated in a declared Gaussian basis. The canonical benchmark uses helium/cc-pVTZ canonical molecular orbitals and retains the lowest 10 spatial orbitals.

## Qubit map

Each spin orbital maps to one qubit through Jordan–Wigner. For mode `p`,

\[
a_p=\frac12\left(\prod_{k<p}Z_k\right)(X_p+iY_p),
\qquad
 a_p^\dagger=\frac12\left(\prod_{k<p}Z_k\right)(X_p-iY_p).
\]

The current repository implements this map directly and tests the number operator identity

\[
a_p^\dagger a_p=\frac{I-Z_p}{2}.
\]

## Many-body flavor coordinates

The one-body state geometry is represented by

\[
\gamma_{pq}=\langle\Psi|a_p^\dagger a_q|\Psi\rangle,
\]

while pair correlation is represented by

\[
\Gamma_{pqrs}=\langle\Psi|a_p^\dagger a_q^\dagger a_s a_r|\Psi\rangle.
\]

For a transition `I -> F`, OES-Q1 evaluates the transition one-body density

\[
T^{FI}_{pq}=\langle\Psi_F|a_p^\dagger a_q|\Psi_I\rangle.
\]

Its spin-summed spatial projection is the many-electron continuation of the hydrogen transition-density eclipse field.

## Spin gate

The full 20-spin-orbital sector retains all `M_S` sectors. OES constructs `S^2` explicitly and classifies eigenstates by

\[
\langle S^2\rangle=S(S+1).
\]

The helium ground state is required to resolve as a singlet. The first triplet and the first excited singlet are separately identified. A spin-independent spatial transition density must close the singlet-to-triplet channel while remaining open for at least one singlet-to-singlet channel.

## Independent validation

PySCF is used for two bounded roles only:

1. generation of standard one- and two-electron integrals in the declared basis;
2. an independent FCI energy on the same active space.

OES constructs its own 20-spin-orbital fixed-particle Hamiltonian and diagonalizes it independently. The primary implementation gate is

\[
|E_{FCI}^{OES}-E_{FCI}^{PySCF}|<10^{-9}\ E_h.
\]

No experimental helium energy is an input to this gate.

## Backend status

The repository benchmark backend is declared

`SIMULATED_REFERENCE`.

A physical-qubit backend may be attached later without changing the Hamiltonian, qubit map, active-space contract, or benchmark targets. Backend provenance must remain explicit in every receipt.
