# OES-Q1 Validation Receipt

Date: 2026-09-01
Branch: `oes-q1-quantum-helium-20q`
Status: OPEN_PENDING_HOSTED_CI
Backend: `SIMULATED_REFERENCE`

## Canonical register

- spatial orbitals: 10
- spin orbitals / Jordan-Wigner qubits: 20
- electrons: 2
- full qubit dimension: 1,048,576
- exact fixed-N reference dimension: 190
- orbital basis: helium `cc-pVTZ`, lowest 10 canonical RHF MOs

## Gates

1. 20Q register and N=2 sector dimensions — IMPLEMENTED / CI PENDING
2. Jordan-Wigner number-operator identity — IMPLEMENTED / CI PENDING
3. fixed-sector fermionic Hamiltonian — IMPLEMENTED / CI PENDING
4. OES FCI vs independent PySCF FCI on identical integrals — IMPLEMENTED / CI PENDING
5. S^2 singlet/triplet classification — IMPLEMENTED / CI PENDING
6. transition 1-RDM — IMPLEMENTED / CI PENDING
7. spin-summed singlet→triplet closure — IMPLEMENTED / CI PENDING
8. singlet→singlet transition-density visibility — IMPLEMENTED / CI PENDING
9. two-body RDM primitive — IMPLEMENTED / CI PENDING
10. physical quantum backend execution — OPEN
11. experimental helium spectrum benchmark — OPEN

## Promotion rule

No implemented gate becomes PASS until hosted CI runs the exact branch head successfully. A successful simulated-reference gate does not promote the backend to PHYSICAL.
