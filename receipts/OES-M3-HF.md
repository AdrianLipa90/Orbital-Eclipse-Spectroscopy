# OES-M3 — HF Reduced-Active Fixed-20Q Spectroscopy Receipt

Status: `PASS__HF_FROZEN_4E_REDUCED_ACTIVE_20Q_BASELINE`
Backend: `SIMULATED_REFERENCE`

## Scope

OES-M3 is the first blind spectroscopic benchmark on a ten-electron heteronuclear molecule in the generalized fixed-20Q architecture.

The validated baseline uses:
- HF / cc-pVTZ,
- 10 total electrons,
- two lowest closed-shell spatial orbitals integrated through the exact G1 frozen-core reduction (4 frozen electrons),
- 6 active electrons in 10 active spatial orbitals = 20 spin orbitals = 20Q,
- exact `M_S=0` reference sector

\[
D=\binom{10}{3}^2=14400,
\]

- the exact G2 sparse Slater-Condon Hamiltonian.

Experimental constants do not enter the RHF geometry search, active-space selection, Hamiltonian, FCI solution, density reconstruction, or local potential fit.

## Predeclared physical baseline gate

Before observing the M3 result, the relative tolerance was fixed at 5% independently for:
- equilibrium bond length `r_e`,
- harmonic wavenumber `omega_e`,
- equilibrium rotational constant `B_e`.

The threshold was not modified after the result.

Output benchmarks from the NIST Chemistry WebBook fixture:
- `r_e = 0.916808 Å`,
- `omega_e = 4138.32 cm^-1`,
- `B_e = 20.9557 cm^-1`.

## Blind geometry and fixed-20Q result

Hosted workflow run: `33590666397`, job `100123825543`.
Head used for the calculation: `fc26c14d5d3da2d21c3941255ef84ff4794df6d5`.

The broad RHF scan used only the predetermined grid

`[1.30, 1.45, 1.60, 1.75, 1.90, 2.05, 2.20] bohr`

and found its internal seed minimum at `1.75 bohr`. Five active-space FCI points were then evaluated at

`[1.63, 1.69, 1.75, 1.81, 1.87] bohr`.

The fitted results are:
- `r_e = 0.9266821958995203 Å`, residual `+0.009874195899520322 Å`, relative `+1.0770189504803974%`,
- `omega_e = 4071.1474557680503 cm^-1`, residual `-67.17254423194936 cm^-1`, relative `-1.6231839063182492%`,
- `B_e = 20.52239633675788 cm^-1`, residual `-0.4333036632421212 cm^-1`, relative `-2.067712666444553%`.

All three predeclared physical 5% gates PASS.

## Independent Hamiltonian and density gates

At the fitted equilibrium geometry:
- active qubits: `20`,
- fixed-Ms dimension: `14400`,
- CSR nnz: `2,181,888`,
- CSR density: `0.010522222222222223`,
- OES sparse energy: `-100.11423738242877 Eh`,
- independent PySCF active-FCI energy: `-100.11423738242884 Eh`,
- energy residual: `+7.105427357601002e-14 Eh`,
- maximum active 1-RDM residual against independent PySCF FCI: `3.350412457886476e-08`.

The active 1-RDM is combined with the analytically occupied frozen core for all one-body observables.

## Heteronuclear relational observables

From the same reconstructed many-electron density:
- fluorine-centered `1/r` exposure: `27.10282468710246`,
- hydrogen-centered `1/r` exposure: `6.068337794421346`,
- exposure difference: `+21.034486892681112`,
- permanent dipole magnitude: `1.7341885278732643 D`.

The NIST CCCBDB experimental molecular-beam `mu_0` value `1.827 D` is retained as an independent diagnostic benchmark. The M3 value is approximately `-5.08%` relative to that diagnostic. Dipole agreement was not part of the predeclared M3 physical PASS gate and is not retroactively promoted to one.

## Verdict

- exact sparse Hamiltonian cross-check: `PASS`,
- active 1-RDM cross-check: `PASS`,
- heteronuclear exposure/polarization resolution: `PASS`,
- blind `r_e` 5% gate: `PASS`,
- blind `omega_e` 5% gate: `PASS`,
- blind `B_e` 5% gate: `PASS`,
- independent dipole diagnostic: recorded, not a gate,
- active-electron convergence beyond the 6e baseline: `OPEN`,
- anharmonic/rovibrational HF spectrum: `OPEN`,
- physical-QPU execution: `OPEN`.

The 8-active-electron fixed-20Q convergence branch is tracked separately so that this 6e baseline remains an immutable validation point.
