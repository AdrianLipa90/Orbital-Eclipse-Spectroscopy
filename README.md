# Orbital-Eclipse-Spectroscopy

Orbital Eclipse Spectroscopy (OES) is a calculation-first repository for representing an atomic spectral transition as a **gross frequency/color coordinate plus a resolved orbital-flavor geometry**.

Current milestone: **OES-H0 — Hydrogen Orbital-Eclipse Closure**.

The executable reference layer currently includes:

- reduced-mass Coulomb gross energies and colors,
- orbital E1 channel geometry,
- radial-node structure,
- central and contact exposure kernels,
- normalized hydrogen radial functions,
- signed/unsigned E1 radial overlaps and cancellation coherence,
- reduced-mass Dirac fine-structure reference energies,
- locked NIST H I `np` fine-structure benchmarks for `n=2..5`,
- weak-field Zeeman flavor-to-color conversion,
- continuous p-state Zeeman → Paschen–Back crossover,
- leading low-Z `2S-2P1/2` QED reference from A41/A40 self-energy plus leading Uehling vacuum polarization.

Solver inputs are physical constants, quantum numbers and declared theoretical coefficients. Empirical spectral data live in benchmark fixtures and remain outside solver input paths.

The executable code is in `oes/`, formal definitions are in `formalism/`, held-out targets are in `benchmarks/`, and validation state is recorded in `receipts/`.

## Current status

Hosted CI has validated the gross, flavor, radial, Dirac, magnetic, leading-QED and `np` fine-series reference gates on exact branch heads.

Complete higher-order QED closure and a broad hydrogen line/intensity benchmark remain explicit OPEN gates.
