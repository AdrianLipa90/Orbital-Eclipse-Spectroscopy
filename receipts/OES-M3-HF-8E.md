# OES-M3.1 — HF Eight-Active-Electron Fixed-20Q Scaling Receipt

Status: `PASS__HF_FROZEN_2E_EIGHT_ACTIVE_ELECTRON_20Q_SCALING_GATE`
Backend: `SIMULATED_REFERENCE`

## Scope

OES-M3.1 tests whether the same fixed 20Q register can retain eight active electrons rather than six. The geometry is deliberately non-benchmark (`R = 1.8 bohr`) so this receipt measures algebraic correctness and computational scaling, not agreement with an experimental HF constant.

The validated partition is:
- HF / cc-pVTZ,
- total electrons: `10`,
- one deepest closed-shell spatial orbital frozen through the exact G1 reduction (`2` frozen electrons),
- `8` active electrons in `10` spatial orbitals = `20` spin orbitals = `20Q`,
- exact `M_S=0` sector

\[
D=\binom{10}{4}^2=44100.
\]

No experimental observable enters the gate.

## Hosted-CI result

Workflow run: `33590901509`.
Job: `100124515216`.
Validated head: `5dce4d2ad3bf0c27bfbd333d8f34cdcd2f5fedd1`.

Algebra regressions for G1/G2/M3 primitives PASS before the scaling calculation.

Eight-electron result:
- active MO indices: `[1,2,3,4,5,6,7,8,9,12]`,
- complete active group sizes: `[1,1,2,1,1,1,2,1]`,
- fixed-Ms dimension: `44,100`,
- maximum Slater-Condon connectivity per row: `805`,
- CSR nnz: `10,021,172`,
- CSR density: `0.005152776877946946`,
- CSR storage: `120,430,468 bytes`,
- sparse build time on hosted runner: `72.447747195 s`,
- lowest-eigenvalue solve: `0.8950386320000234 s`.

Energy cross-check:
- OES sparse energy: `-100.10916226283878 Eh`,
- independent PySCF active-FCI energy: `-100.10916226283166 Eh`,
- residual: `-7.119638212316204e-12 Eh`.

## Verdict

`PASS`: the exact fixed-20Q architecture supports an eight-active-electron `44,100`-determinant molecular reference sector on the hosted CI runner while preserving agreement with independent active-space FCI far inside the declared `5e-8 Eh` gate.

This receipt establishes computational feasibility of eight-active-electron fixed-20Q molecular sectors. Spectroscopic convergence relative to the six-electron HF baseline, broader active-orbital convergence, larger multicentre molecules, and physical-QPU execution remain separate milestones.
