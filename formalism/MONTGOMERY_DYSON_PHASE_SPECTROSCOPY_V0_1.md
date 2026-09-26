# Montgomery–Dyson Phase Spectroscopy v0.1 — OES adapter

Status split:
- phase-correlation core: DERIVED_IN_FRAMEWORK / FORCED_PREDICTION;
- OES atomic/molecular binding: CANDIDATE / NOT PROMOTED.

Dependency theorem: AdrianLipa90/Infinities, branch formalize/hardy-car-sine-kernel-forced-prediction-v0.1-20260926, head 7008a9121d878674a7190d0966faa727f9eacb55.

OES owns the transition layer. This adapter treats a resolved line or transition coordinate as a scalar spectrum, unfolds it, and maps it to an unwrapped phase:

u_j = (x_j - x_0) / <Delta x>,    Phi_j = 2*pi*u_j.

The Montgomery–Dyson/GUE reference correlation is represented in phase coordinates as

R_2(Delta Phi) = 1 - [sin(Delta Phi/2)/(Delta Phi/2)]^2.

The independent arithmetic reference probe uses prime-power channels

omega_(p,m) = m log p = log(p^m),
w_(p,m) = log(p)/p^(m/2),

with the finite smoothed signal

P_sigma(t) = -(1/pi) sum w_(p,m) exp[-sigma^2 omega_(p,m)^2/2]
             cos(t omega_(p,m) + varphi_(p,m)).

## OES binding

1. Physical line positions and intensities still come from OES atomic/molecular solvers and declared controls.
2. transition_phase_coordinates is a downstream diagnostic on line coordinates; it does not modify energies, intensities, selection rules or Hamiltonians.
3. GUE and Poisson are comparison/null families, not labels forced onto a spectrum.
4. The prime-power signal is a mathematical reference/control basis. No claim is made that atoms or molecules are physically driven by primes.
5. Reverse control preserves frequencies and weights and randomizes channel phases.
6. Domain-specific smooth unfolding must replace the affine control before any serious density-varying spectral claim.

This adapter connects OES transition geometry to the shared spectroscopy chain while preserving the calculation-first evidence hierarchy.


## Forced-prediction provenance

The sinc-square/GUE functional form is no longer introduced merely as a comparison curve. In the declared Hardy–CAR sector it is derived from the consecutive unilateral-shift projector plus the filled CAR/Slater determinant:

g2_N(s) = 1 - [sin(pi s)/(N sin(pi s/N))]^2

and

g2(s) = 1 - [sin(pi s)/(pi s)]^2.

The OES-specific open question is whether a given atomic or molecular transition ensemble instantiates that projector sector. That binding remains candidate and must be tested independently.
