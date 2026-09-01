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
13. complete-class external natural-orbital bath rotation gauge — PASS
14. H0 regression suite on Q1 head — PASS

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

Verdict: `20Q_CAPACITY_PASS`. The ten-spatial-orbital budget is sufficient for the tested ground/2s/2p information when the subspace is chosen near-optimally. The oracle remains a capacity upper bound; predictive selection is handled separately.

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

Verdict: top natural directions from these mixed ground-state operator-response densities do not recover the correct bright-block orientation at the tested construction. The negative ablations are retained as diagnostics; their weights were not fitted to the benchmark spectrum.

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

Verdict: `PARTIAL`. Diagonal-Q dressing recovers a large fraction of the common color offset and improves relative source agreement. Its 2.10 meV raw bright-block spread fails the exact bright-manifold symmetry gate. The result remains a preserved perturbative diagnostic.

## Determinant-selected Q diagnostic

A state-balanced determinant selector ranked external determinants by Hamiltonian coupling to active classes and retained full Q-Q couplings inside the selected determinant subspace. It strongly recovered color: independent hosted runs near the nominal 512-determinant point reached NIST RMS around 0.049-0.050 eV and source RMS around 0.018-0.025 eV.

The rerun audit exposed one representation-sensitive degree of freedom. Degenerate external orbitals can rotate without changing the physical source subspace; determinant identities and numerical importance ties then change. Two hosted executions selected 512 and 513 determinants respectively and produced slightly different energies and bright spreads (approximately 5.03 and 3.32 meV), while retaining the same qualitative source-scale convergence.

Verdict: `COLOR_RECOVERY_PASS__GAUGE_OPEN`. This path remains a historical convergence diagnostic. The canonical external-bath selector is the rotation-covariant natural-orbital construction below.

## Rotation-covariant complete-class natural-orbital bath

The canonical bath is built from physical coupling-response states

`|chi_c> = Q H |Psi_c^P>`

for complete active classes:

- ground singlet: 1 state
- lowest triplet: complete three-state spin manifold
- lowest dark excited singlet: 1 state
- first bright singlet: complete three-state p manifold

For every class, normalized Q-coupling states produce a spin-summed external one-body density. Class densities are trace-normalized and equally weighted. Diagonalization yields external natural orbitals; numerically degenerate occupation eigenspaces are admitted as complete blocks.

Hard rotation gauge: all 52 external spatial orbitals were deliberately mixed by a random orthogonal transformation and the bath was rebuilt from transformed integrals.

Gauge result at the 10-external-orbital bath:

- minimum principal cosine: 0.9999999999999991
- maximum principal-cosine error: 7.99e-15
- projector Frobenius distance: 0.0
- selected external spatial dimension: 10 before and after rotation
- retained normalized response occupation: 0.82476801865 in both representations

Verdict: `Q_BASIS_ROTATION_GAUGE_PASS`.

The external response occupation groups begin as

`1, 3, 1, 5, 3, 5, 3, 7, ...`

and are retained as whole eigenspaces. Hosted convergence:

| external natural orbitals | total spatial orbitals | NIST RMS (eV) | source RMS (eV) | bright spread (eV) | triplet spread (eV) | bright f sum |
|---:|---:|---:|---:|---:|---:|---:|
| 4  | 14 | 0.22948 | 0.24474 | 2.64e-11 | 4.97e-14 | 0.32008 |
| 10 | 20 | 0.20109 | 0.21509 | 7.60e-11 | 2.17e-13 | 0.32056 |
| 13 | 23 | 0.05566 | 0.02862 | 7.62e-11 | 1.81e-13 | 0.31576 |

The requested 7- and 10-orbital prefixes both close at the same 10-dimensional complete eigenspace because the fourth occupation group has dimension five.

At 13 external natural orbitals:

- triplet: 19.76843 eV; NIST residual -0.05118 eV; source residual -0.03439 eV
- dark singlet: 20.57865 eV; NIST residual -0.03712 eV; source residual -0.03462 eV
- bright singlet: 21.29079 eV; NIST residual +0.07277 eV; source residual -0.00870 eV
- source RMS: 0.02862 eV
- source centered RMS: 0.01217 eV
- NIST RMS: 0.05566 eV versus approximately 0.04805 eV for the full 62-orbital d-aug source
- bright manifold spread: 7.62e-11 eV
- triplet manifold spread: 1.81e-13 eV
- bright oscillator-strength sum: 0.31576
- active P-space weights: ground 0.99571, triplet mean 0.99875, dark 0.99782, bright mean 0.99914

Verdict: `NATURAL_BATH_COLOR_AND_SYMMETRY_PASS`. The fixed 20Q flavor core plus a 13-spatial-orbital rotation-covariant classical bath reaches the tested full-source accuracy scale while preserving complete triplet and bright manifolds to numerical precision.

## Current gates

- exact 20Q implementation — PASS
- 20Q information-capacity upper bound for tested He sectors — PASS
- predictive symmetry-block flavor preservation — PASS
- external bath Q-basis rotation gauge — PASS
- complete triplet/bright manifold closure in the canonical bath — PASS
- natural-bath recovery of source-scale color accuracy — PASS
- broader He I line/intensity benchmark — OPEN
- physical quantum backend execution — OPEN

## Promotion rule

Q1 remains `PARTIAL_PASS` while broader helium spectroscopy and physical-backend execution remain open. The tested simulated-reference architecture promoted by this receipt is

`fixed 20Q symmetry/flavor core + complete-class rotation-covariant external natural-orbital bath`.
