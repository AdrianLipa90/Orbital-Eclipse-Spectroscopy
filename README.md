# Orbital-Eclipse-Spectroscopy

Orbital Eclipse Spectroscopy (OES) is a calculation-first repository for representing an atomic spectral transition as a **gross frequency/color coordinate plus a resolved orbital-flavor geometry**.

Current milestone: **OES-H0 — Hydrogen Orbital-Eclipse Closure**.

The executable reference layer currently includes:

- reduced-mass Coulomb gross energies and colors,
- orbital E1 channel geometry,
- radial-node structure,
- central and contact exposure kernels,
- reduced-mass Dirac fine-structure reference energies,
- locked NIST H I fine-structure benchmark targets,
- weak-field Zeeman flavor-to-color conversion,
- continuous p-state Zeeman → Paschen–Back crossover,
- normalized hydrogen radial functions,
- signed/unsigned E1 radial overlaps and cancellation coherence.

The gross solver takes physical constants and quantum numbers as inputs. Empirical spectral data live in benchmark fixtures and remain outside solver input paths.

The executable code is in `oes/`, formal definitions are in `formalism/`, held-out targets are in `benchmarks/`, and validation state is recorded in `receipts/`.

## Current status

Hosted CI has already validated the gross, flavor, Dirac and magnetic reference layers on earlier exact branch heads. The current radial-extension head is undergoing its own exact-head validation.

Full Lamb/QED closure and a broader held-out hydrogen level/line benchmark remain explicit OPEN gates.
