# OES-G1 — Frozen-Core Fixed-20Q Active Hamiltonian Receipt

Status: `PASS__FROZEN_CORE_ACTIVE_HAMILTONIAN`
Backend: `SIMULATED_REFERENCE`

## Scope

OES-G1 integrates a chosen closed-shell doubly occupied core into a scalar core energy and an effective one-electron Hamiltonian acting on a fixed active register. The canonical gate keeps 10 active spatial orbitals = 20 spin orbitals = 20Q.

For chemists' ERIs `(pq|rs)`:

\[
E_{\rm core}=E_{\rm nuc}+2\sum_i h_{ii}+\sum_{ij}[2(ii|jj)-(ij|ji)],
\]

\[
h^{\rm eff}_{pq}=h_{pq}+\sum_i[2(pq|ii)-(pi|iq)].
\]

The active two-electron integrals are retained unchanged.

## Independent gates

1. Exact random-Hamiltonian projection: the frozen-core effective Hamiltonian equals the principal full-FCI block in which the chosen core orbital is doubly occupied, within `2e-12` matrix tolerance.
2. LiH CASCI integration gate at a deliberately non-benchmark geometry `R=3.0 bohr`, cc-pVTZ, one frozen doubly occupied spatial orbital, 10 active spatial orbitals, 2 active electrons.
3. No experimental observable enters the Hamiltonian construction or the validation geometry.

## Hosted-CI result

Canonical G1 workflow run: `33589975192`.

- active qubits: `20`
- active fixed-particle dimension: `190`
- max `|Delta h_eff|`: `1.1934897514720433e-15 Eh`
- max `|Delta ERI|`: `1.231653667943533e-15 Eh`
- `Delta E_core`: `-8.881784197001252e-16 Eh`
- OES active-FCI energy: `-7.987394067407229 Eh`
- PySCF CASCI energy: `-7.987394067407222 Eh`
- energy residual: `-7.105427357601002e-15 Eh`

## Verdict

`PASS`: the implemented frozen-core reduction reproduces the independently generated CASCI active Hamiltonian and ground-state energy at numerical precision for the validated gate.

This receipt establishes the active-space reduction layer used by later fixed-20Q molecular benchmarks. Broader choices of frozen cores, active-orbital optimization, spectroscopy-specific error budgets, and physical-QPU execution remain separate validation milestones.
