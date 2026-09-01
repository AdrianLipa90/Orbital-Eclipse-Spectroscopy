# OES-Q1 Validation Receipt

Date: 2026-09-01
Branch: `oes-q1-quantum-helium-20q`
Status: PARTIAL_PASS
Backend: `SIMULATED_REFERENCE`

## Canonical register

- spatial active orbitals: 10
- spin orbitals / Jordan-Wigner qubits: 20
- electrons: 2
- full qubit dimension: 1,048,576
- exact fixed-N active reference dimension: 190

## Exact implementation gates

1. 20Q register and N=2 sector dimensions — PASS
2. Jordan-Wigner number-operator identity — PASS
3. fixed-sector fermionic Hamiltonian — PASS
4. OES FCI vs independent PySCF FCI on identical active-space integrals — PASS; residual at machine precision (~1e-15 Ha)
5. S^2 singlet/triplet classification — PASS
6. transition 1-RDM — PASS
7. spin-summed singlet→triplet closure — PASS
8. singlet→singlet transition-density visibility — PASS
9. two-body RDM primitive — PASS
10. arbitrary selected two-electron determinant Hamiltonian vs canonical sector Hamiltonian — PASS
11. state-specific EN2 toy model and intruder fail-closed gate — PASS
12. degenerate EN2 block invariance under active-basis rotation — PASS
13. H0 regression suite on Q1 head — PASS

These gates validate the simulated-reference implementation. Physical-backend execution remains a separate OPEN gate.

## Helium spectroscopy targets

NIST comparison targets are benchmark outputs only and do not enter active-space or bath selection:

- `1s2s 3S1`: 19.81961436 eV
- `1s2s 1S0`: 20.61577465 eV
- `1s2p 1P1`: 21.21802253 eV
- bright-manifold oscillator-strength reference: approximately 0.2762

## Source representation

Full `aug-cc-pVQZ` FCI source residuals:

- triplet: +0.05333 eV
- dark singlet: +0.25034 eV
- bright singlet: +2.82316 eV

Geometrically generated `d-aug-cc-pVQZ` FCI source residuals:

- triplet: -0.01680 eV
- dark singlet: -0.00250 eV
- bright singlet: +0.08147 eV
- NIST RMS for these three levels: approximately 0.04805 eV
- first bright manifold: three components
- full-source bright oscillator-strength sum: approximately 0.32933

Verdict: d-aug resolves the dominant source-representation error of the singly augmented bright Rydberg state and defines the source-quality ceiling used by the fixed-20Q compression diagnostics.

## Fixed-20Q capacity upper bound

A nonpredictive, symmetry-aware state-averaged FCI oracle was used only to measure whether ten spatial orbitals can contain the target information. It selected complete occupation groups

`1 + 1 + 3 + 1 + 1 + 3 = 10`

and retained 1.99953194 of two electrons in the state-averaged one-body occupation measure (99.9766%).

Oracle 20Q results:

- triplet: 19.64739 eV; NIST residual -0.17223 eV
- dark singlet: 20.46866 eV; NIST residual -0.14712 eV
- bright singlet: 21.17287 eV; NIST residual -0.04515 eV
- bright degeneracy: 3
- bright oscillator-strength sum: 0.33179; 100.75% of the full d-aug source value

Verdict: `20Q_CAPACITY_PASS`. The ten-spatial-orbital budget is sufficient for the tested ground/2s/2p information when the subspace is chosen near-optimally. The oracle is a capacity upper bound and is not promoted as a predictive selector.

## Predictive symmetry-block active space

The strongest structure-preserving predictive candidate is

`D-AUG-SYMMETRY-BLOCK-S4-P3-P3R2-20Q`

with four scalar radial modes plus two complete three-dimensional p blocks. It uses exactly ten spatial / twenty spin orbitals.

Active-space results:

- triplet: 19.12482 eV; source residual -0.67800 eV
- dark singlet: 19.97146 eV; source residual -0.64181 eV
- bright singlet: 20.63277 eV; source residual -0.66672 eV
- source RMS: 0.66235 eV
- source centered RMS after removing the common offset: 0.01512 eV
- bright degeneracy: 3
- bright oscillator-strength sum: 0.33060 vs source 0.32933

Verdict: the predictive 20Q core preserves relative flavor structure and transition strength very well while carrying an approximately common color offset from the excluded external space.

## Operator-response selector ablations

`ROAS-v1` used equal-weight response classes `{1, r^2, r^4, r, r^2 r}`. It retained complete block dimensions but produced a bright-state residual of +1.6166 eV and bright oscillator-strength sum 0.74885. Spectral gate: FAIL.

`ROAS-v2` used equal-weight `{1, r^2, r^4, r}`. It produced an accurate triplet (-0.01505 eV vs NIST) and dark residual +0.14767 eV, while the bright state moved to 28.34205 eV with oscillator-strength sum 1.39829. Bright spectral gate: FAIL.

Verdict: top natural directions from mixed ground-state operator-response densities do not yet recover the correct bright-block orientation. The negative ablations are retained as diagnostics; their weights were not fitted to the benchmark spectrum.

## External-space EN2 diagnostic

The fixed 20Q `s4+p6` P-space was coupled to the remaining 7,436 two-electron determinants of the 124-spin-orbital d-aug source. Unshifted Epstein-Nesbet second order used no fitted level shift. Minimum tested denominators remained well above the 1e-5 Ha intruder gate.

A quasi-degenerate 3x3 second-order block was used for the bright manifold. Symmetry-projected color diagnostics gave:

- triplet: 19.58249 eV
- dark singlet: 20.40297 eV
- bright singlet: 21.07801 eV
- NIST RMS: 0.64302 -> 0.20093 eV
- source RMS: 0.66235 -> 0.21743 eV
- source centered RMS: 0.01512 -> 0.00502 eV

Raw quasi-degenerate EN2 bright-block spread: 0.002102 eV.

Verdict: `PARTIAL`. Diagonal-Q dressing recovers a large fraction of the common color offset and improves relative source agreement, while the diagonal-Q approximation introduces a 2.10 meV rotational-symmetry splitting and therefore does not pass the bright-manifold symmetry gate.

## State-balanced selected-Q CI

The same fixed 20Q `s4+p6` P-space was retained. External determinants were ranked only from equal-weight Hamiltonian-coupling importance for four P-space classes: ground, triplet, dark singlet and the complete three-component bright manifold. NIST values did not enter selection. Selected P+Q subspaces retained full Q-Q Hamiltonian couplings and were diagonalized exactly inside each selected space.

Hosted selected-Q sweep:

| selected Q determinants | NIST RMS (eV) | source RMS (eV) | source centered RMS (eV) | bright spread (eV) | bright f sum |
|---:|---:|---:|---:|---:|---:|
| 32  | 0.41392 | 0.43182 | 0.00871 | 1.36e-9 | 0.32616 |
| 64  | 0.34826 | 0.36507 | 0.01291 | 1.36e-9 | 0.32426 |
| 128 | 0.23304 | 0.24695 | 0.01855 | 0.00206 | 0.32335 |
| 256 | 0.12168 | 0.12734 | 0.01755 | 0.00495 | 0.32307 |
| 512 | 0.04993 | 0.02483 | 0.00830 | 0.00503 | 0.31931 |

At 512 selected external determinants:

- triplet: 19.78022 eV; NIST residual -0.03939 eV; source residual -0.02259 eV
- dark singlet: 20.57933 eV; NIST residual -0.03644 eV; source residual -0.03394 eV
- bright singlet: 21.28583 eV; NIST residual +0.06781 eV; source residual -0.01366 eV
- source RMS reduction relative to the bare 20Q core: 0.66235 -> 0.02483 eV (about 26.7x)
- NIST RMS reduction relative to the bare 20Q core: 0.64302 -> 0.04993 eV (about 12.9x)
- full d-aug source NIST RMS: approximately 0.04805 eV
- P-space weights remain above 0.995 for ground/dark/bright and above 0.998 for the matched triplet eigenstate weight; active-target overlap tracking remains high

Verdict: `SELECTED_Q_COLOR_RECOVERY_PASS` for convergence to the tested d-aug source-quality scale. The 20Q core plus a coupling-selected classical bath reproduces the three source excitation energies to 0.02483 eV RMS and reaches approximately the NIST RMS of the full source representation. Exact angular closure of the selected bath remains OPEN because the 512-determinant truncation produces a 5.03 meV bright-manifold spread.

## Current gates

- exact 20Q implementation — PASS
- 20Q information-capacity upper bound for tested He sectors — PASS
- predictive symmetry-block flavor preservation — PASS at the tested active-space diagnostics
- selected-Q recovery of source-scale color accuracy — PASS
- exact angular closure of selected-Q bath — OPEN
- broader He I line/intensity benchmark — OPEN
- physical quantum backend execution — OPEN

## Promotion rule

Q1 remains `PARTIAL_PASS` while angular bath closure, broader helium spectroscopy and physical-backend execution remain open. Hosted simulated-reference PASS gates support the implemented 20Q/core-plus-bath architecture at the tested helium benchmark only.
