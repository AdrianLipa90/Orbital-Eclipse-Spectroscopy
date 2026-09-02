# OES-G2 — Sparse Exact Fixed-Ms 20Q Receipt

Status: `PASS__SPARSE_FIXED_MS_20Q_SIX_ELECTRON_GATE`
Backend: `SIMULATED_REFERENCE`

## Scope

OES-G2 changes only the classical reference representation of the same spin-orbital Hamiltonian. Determinants are enumerated directly in an exact `(N_alpha,N_beta)` sector and Slater-Condon connectivity is stored in CSR form.

For 10 active spatial orbitals and six active electrons:

\[
N_\alpha=N_\beta=3,
\qquad
D=\binom{10}{3}^2=14400.
\]

The operator uses exact diagonal, single-excitation and double-excitation Slater-Condon matrix elements with the same fermionic sign convention as the existing dense OES builder.

## Algebra gates

Sparse and dense matrices are compared directly for random Hamiltonians in:
- two-electron `N_alpha=N_beta=1`,
- three-electron `N_alpha=2,N_beta=1`,
- four-electron `N_alpha=N_beta=2` sectors.

All matrix comparisons pass at `3e-12` tolerance. The 20Q dimensions `14400` for six electrons and `44100` for eight electrons are independently gated.

## Hosted HF stress gate

Workflow run: `33590359956`, job `100122930333`.
Geometry `R=1.8 bohr` is deliberately non-benchmark. No experimental observable is an input.

HF/cc-pVTZ setup:
- total electrons: `10`
- frozen-core electrons: `4`
- active electrons: `6`
- active spatial orbitals: `10`
- active qubits: `20`
- fixed-Ms dimension: `14400`
- maximum connectivity per row: `610`

Sparse result:
- CSR nnz: `4,105,264`
- CSR density: `0.019797762345679013`
- CSR storage: `49,320,772 bytes`
- build time on hosted runner: `17.581992196 s`
- lowest-eigenvalue solve: `0.321180681 s`

Energy cross-check:
- OES sparse energy: `-100.11335481617195 Eh`
- PySCF CASCI energy: `-100.11335481617205 Eh`
- residual: `9.947598300641403e-14 Eh`

## Verdict

`PASS`: the sparse fixed-Ms implementation reproduces the independently generated CASCI ground-state energy at numerical precision while reducing the 14,400-state reference Hamiltonian to a 1.98% sparse operator.

This gate establishes computational access to six-active-electron fixed-20Q molecular sectors. The eight-active-electron `44,100`-determinant sector, spectroscopy observables, broader active-space convergence and physical-QPU execution remain separate milestones.
