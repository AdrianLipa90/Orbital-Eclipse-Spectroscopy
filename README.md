# Orbital-Eclipse-Spectroscopy

Orbital Eclipse Spectroscopy (OES) is a calculation-first repository for representing an atomic spectral transition as a **gross frequency/color coordinate plus a resolved orbital-flavor geometry**.

Current integrated frontier: **OES-M3.1 — fixed-20Q molecular reference scaling**, with the earlier **OES-H0 hydrogen orbital-eclipse closure** retained as the atomic spectroscopy baseline.

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
- leading low-Z `2S-2P1/2` QED reference from A41/A40 self-energy plus leading Uehling vacuum polarization,
- H2 fixed-20Q molecular reference baseline (OES-M1),
- LiH fixed-20Q and rovibrational reference baseline (OES-M2/M2.1),
- exact frozen-core active-Hamiltonian reduction (OES-G1),
- sparse exact fixed-Ms 20Q solver (OES-G2),
- blind HF reduced-active fixed-20Q spectroscopy baseline (OES-M3),
- eight-active-electron fixed-20Q scaling gate (OES-M3.1).

Solver inputs are physical constants, quantum numbers and declared theoretical coefficients. Empirical spectral data live in benchmark fixtures and remain outside solver input paths.

The executable code is in `oes/`, formal definitions are in `formalism/`, held-out targets are in `benchmarks/`, and validation state is recorded in `receipts/`.

## Current status

Hosted CI has validated the gross, flavor, radial, Dirac, magnetic, leading-QED and `np` fine-series reference gates on exact branch heads. The integrated M3.1 implementation head `9d732a895f2627adcd40037ccb56de1391db96cc` also completed the hosted OES reference suite and the H2, LiH, frozen-core, sparse fixed-Ms, HF spectroscopy and eight-active-electron fixed-20Q workflow gates successfully.

The molecular/QPU layers described above remain **SIMULATED_REFERENCE** unless a receipt explicitly states otherwise. Physical-QPU execution is not claimed.

Complete higher-order QED closure, broader hydrogen line/intensity benchmarking, systematic molecular active-space/basis convergence and physical-QPU execution remain explicit OPEN gates.
