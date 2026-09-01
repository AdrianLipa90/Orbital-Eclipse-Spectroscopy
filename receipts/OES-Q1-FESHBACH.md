# OES-Q1 Feshbach / Schur Downfolding Receipt

Date: 2026-09-01
Branch: `oes-q1-quantum-helium-20q`
Backend: `SIMULATED_REFERENCE`
Status: `PASS__EXACT_SELECTED_BATH_DOWNFOLDING_TO_FIXED_20Q_SECTOR`

## Contract

The canonical quantum core remains fixed at 10 spatial / 20 spin orbitals. For helium with two electrons its exact particle-number sector has

\[
\binom{20}{2}=190
\]

states.

A rotation-covariant complete-class natural-orbital bath was first constructed without experimental energies. The strongest tested bath contains 13 external spatial orbitals. The independent P+Q reference therefore contains 23 spatial / 46 spin orbitals and 1,035 two-electron determinants.

The 845 Q-space determinants are eliminated exactly through

\[
H_{\mathrm{eff}}(E)=H_{PP}+H_{PQ}(EI-H_{QQ})^{-1}H_{QP}.
\]

The resulting energy-dependent effective operator acts only on the 190-dimensional fixed-N sector of the original 20Q P-space.

## Exact downfolding gates

Hosted GitHub Actions gate `helium-feshbach-20q` on commit `33cf851faa66dbb3c9abeb13c93329db1df9706e` — PASS.

- P-space dimension: 190
- integrated Q-space dimension: 845
- P-block reconstruction error: 1.67e-16 Ha
- maximum effective-eigenpair residual: 1.42e-13 Ha
- maximum reconstructed-Q amplitude error: 1.02e-14
- minimum distance of tracked eigenvalues from the QHQ spectrum: 1.58498 Ha
- bright manifold spread: 1.33e-10 eV
- triplet manifold spread: 2.66e-13 eV

The exact Schur reduction therefore reproduces the tracked P+Q eigenpairs at numerical precision while preserving the complete bright and triplet manifolds.

## Spectroscopic state of the downfolded model

The exact roots inherited from the symmetry-preserving bath are:

- `1s2s 3S`: 19.76843055 eV; NIST residual -0.05118381 eV
- `1s2s 1S`: 20.57865320 eV; NIST residual -0.03712145 eV
- `1s2p 1P`: 21.29079267 eV; NIST residual +0.07277014 eV
- three-level NIST RMS: 0.05565752 eV

Tracked P-space weights remain high:

- ground: 0.99571
- triplet components: 0.99875
- dark singlet: 0.99782
- bright components: 0.99914

Thus the bath contribution is small in wavefunction weight but materially important for absolute color.

## Interpretation status

Canonical implemented architecture:

`fixed 20Q symmetry/flavor P-space + rotation-covariant natural-orbital Q bath + exact Feshbach/Schur elimination`.

The energy-dependent self-energy is currently evaluated classically in the simulated-reference workflow. A physical-backend implementation requires a compiled or iterative representation of `H_eff(E)` on the 20Q register; that execution gate remains OPEN.

Broader He I line/intensity validation also remains OPEN.
