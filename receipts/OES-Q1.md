# OES-Q1 Validation Receipt

Date: 2026-09-01
Branch: `oes-q1-quantum-helium-20q`
Status: PARTIAL_PASS
Backend: `SIMULATED_REFERENCE`

## Canonical register

- spatial orbitals: 10
- spin orbitals / Jordan-Wigner qubits: 20
- electrons: 2
- full qubit dimension: 1,048,576
- exact fixed-N reference dimension: 190

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
10. H0 regression suite on Q1 head — PASS

All PASS labels above are supported by hosted GitHub Actions on the exact Q1 branch head. They validate the simulated-reference implementation only.

## Helium spectroscopy diagnostics

NIST comparison targets used only for benchmark residuals:

- `1s2s 3S1`: 19.81961436 eV
- `1s2s 1S0`: 20.61577465 eV
- `1s2p 1P1`: 21.21802253 eV

### Source representation

Full `aug-cc-pVQZ` FCI source:

- triplet residual: +0.05333 eV
- dark singlet residual: +0.25034 eV
- bright singlet residual: +2.82316 eV

Geometrically generated `d-aug-cc-pVQZ` FCI source:

- triplet residual: -0.01680 eV
- dark singlet residual: -0.00250 eV
- bright singlet residual: +0.08147 eV

Verdict: the dominant bright-state color error in the singly augmented source is a source-representation error, not a 20Q-capacity error.

### d-aug -> fixed 20Q compression

Predictive TDA/NTO selector; NIST energies are not used to choose orbitals:

- active OES FCI vs active PySCF FCI: PASS, residual ~4.4e-16 Ha
- triplet: 20.49826 eV; residual +0.67865 eV
- dark singlet: 20.84918 eV; residual +0.23341 eV
- first bright `1P` manifold: 21.57940 eV; residual +0.36138 eV
- first bright manifold degeneracy: 3
- summed oscillator strength: 0.32223; benchmark reference 0.2762; residual +0.04603
- active ground-state loss relative to full d-aug source: 0.03759 Ha

Verdict: the fixed 20Q space retains the bright manifold to sub-eV accuracy and recovers the expected three-component flavor degeneracy, but the current selector sacrifices too much radial/ground correlation and the low-lying 2s sectors. This is a compression-selection limitation, not an implementation failure.

## Current OPEN gates

1. state-balanced fixed-20Q compression retaining ground, 2s dark/triplet and 2p bright sectors simultaneously — OPEN
2. d-aug source bright-manifold oscillator-strength accounting — diagnostic running
3. broader He I line/intensity benchmark — OPEN
4. physical quantum backend execution — OPEN

## Promotion rule

No physical-backend claim follows from simulated-reference PASS gates. New spectral/compression gates become PASS only after hosted CI runs their exact branch head successfully.
