# OES-A1 Lithium Validation Receipt

Date: 2026-09-02
Branch: `oes-a1-lithium-20q`
Parent: `oes-q1-quantum-helium-20q`
Backend: `SIMULATED_REFERENCE`
Status: `PARTIAL_PASS__LITHIUM_FIXED_20Q_BASELINE`

## Canonical register

Neutral lithium is represented in the same ten-spatial / twenty-spin-orbital register used by Q1.

- electrons: 3
- full qubit dimension: 1,048,576
- complete fixed-N sector: C(20,3) = 1,140 determinants
- exact M_S=+1/2 reference sector: C(10,2) C(10,1) = 450 determinants

The 450-state calculation is an exact S_z symmetry reduction of the same 20-mode Hamiltonian; the complete 1,140-state N=3 Hamiltonian is constructed before the symmetry slice.

## Active-space rule

The first `lowest-10-MO` attempt reproduced the FCI energy but broke the expected nonrelativistic 2p orbital manifold: the first bright gate returned only one component.

The canonical A1 selector therefore admits only complete atomic degeneracy blocks. For `cc-pVTZ` it selects

`1 + 1 + 3 + 3 + 1 + 1 = 10`

corresponding to four scalar radial blocks and two complete p-like triplets. Protocol:

`ATOMIC-SYMMETRY-COMPLETE-S4-P3-P3-20Q`.

No NIST energy or oscillator strength enters this selection.

## Implementation gates

Hosted GitHub Actions on commit `72d4f5c4f40caa6cc015d761228d056802e2f2a4`:

- core regression suite: PASS
- N=3 fixed-sector dimension 1,140: PASS
- M_S=+1/2 dimension 450: PASS
- neutral Li OES FCI vs PySCF FCI on identical active integrals: PASS
- Li+ OES FCI vs PySCF FCI on identical active integrals: PASS
- first bright nonrelativistic orbital manifold degeneracy 3: PASS
- bright manifold spread: approximately 2.2e-13 eV

Neutral active-space FCI:

- OES: -7.445015397678315 Ha
- PySCF: -7.445015397678314 Ha

Li+ active-space FCI:

- OES: -7.248879475094965 Ha
- PySCF: -7.248879475094965 Ha

## Spectroscopic diagnostics

Benchmark values are outputs only.

### First 2s -> 2p bright manifold

- prediction: 1.8791112946 eV
- NIST benchmark: approximately 1.8478 eV
- residual: +0.0313113 eV
- degeneracy: 3
- manifold spread: ~1e-13 eV
- predicted oscillator-strength sum: 0.7999762
- approximate NIST D-multiplet f sum derived from the two tabulated A-values: 0.7467372
- relative f residual: +7.13%

### Li -> Li+ ionization

- prediction: 5.3371303461 eV
- NIST: 5.391714996 eV
- residual: -0.05458465 eV

## Verdict

`A1_IMPLEMENTATION_PASS` — the 20Q fermionic architecture generalizes from two to three electrons, with machine-precision agreement against an independent FCI implementation on identical integrals.

`A1_FLAVOR_BLOCK_PASS` — complete atomic degeneracy blocks restore the three-component p manifold after the deliberately retained lowest-MO failure demonstrated the danger of cutting a symmetry block.

`A1_SPECTROSCOPY_PARTIAL` — the first bright color is within 0.032 eV and ionization within 0.055 eV of the cited NIST targets in the undressed ten-orbital baseline; the D-multiplet strength remains about 7% high.

External-bath/Feshbach refinement for lithium and broader Li I spectroscopy remain OPEN. Physical-QPU execution remains OPEN.
