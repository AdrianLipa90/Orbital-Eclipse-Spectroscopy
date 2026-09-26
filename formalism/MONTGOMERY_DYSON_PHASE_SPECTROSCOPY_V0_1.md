# Montgomery–Dyson Phase Spectroscopy v0.1 — OES adapter

Status: CANDIDATE / IMPLEMENTED, not promoted physical mechanism.

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
