# Orbital-Eclipse-Spectroscopy

Orbital Eclipse Spectroscopy (OES) is a calculation-first repository for representing an atomic spectral transition as a **gross frequency/color coordinate plus a resolved orbital-flavor geometry**.

Current milestone: **OES-H0 — Hydrogen Orbital-Eclipse Closure**.

The first implementation separates:

- gross Coulomb color from `E_i-E_f`,
- orbital E1 channel geometry,
- radial-node structure,
- central exposure kernels,
- contact exposure,
- weak-field Zeeman flavor-to-color conversion.

The executable reference code is in `oes/`, the formal definitions are in `formalism/`, and validation state is recorded in `receipts/`.

## Current status

`OES-H0` is under validation. Analytic hydrogen primitives and linear Zeeman mapping are implemented. Full relativistic fine structure, full Lamb/QED treatment, signed radial cancellation/coherence, empirical held-out spectroscopy benchmarking, and the Paschen–Back continuation remain explicit OPEN gates.

No empirical wavelength is used as an input to the gross hydrogen solver.
