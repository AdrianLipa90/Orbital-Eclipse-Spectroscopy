"""Self-consistent Feshbach root solving on the fixed OES P-space.

The exact downfolded operator is energy dependent,

    H_eff(E) = H_PP + H_PQ (E I - H_QQ)^(-1) H_QP.

This module diagonalizes H_QQ once and then solves

    lambda_c(H_eff(E)) - E = 0

for a P-space state class c. The branch is tracked only by overlap with a seed
subspace supplied from the bare fixed-20Q problem. Full P+Q eigenvalues are not
used by the solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Tuple

import numpy as np

from .feshbach import partition_hamiltonian


@dataclass(frozen=True)
class FeshbachKernel:
    hpp: np.ndarray
    q_eigenvalues: np.ndarray
    q_eigen_coupling: np.ndarray
    p_dimension: int
    q_dimension: int

    @classmethod
    def from_hamiltonian(cls, H: np.ndarray, p_dim: int) -> "FeshbachKernel":
        hpp, hpq, hqq = partition_hamiltonian(H, p_dim)
        qevals, qvecs = np.linalg.eigh(hqq)
        # B = U_Q^T H_QP, shape (q,p).
        coupling = qvecs.T @ hpq.T
        return cls(
            hpp=np.asarray(hpp, dtype=float),
            q_eigenvalues=np.asarray(qevals, dtype=float),
            q_eigen_coupling=np.asarray(coupling, dtype=float),
            p_dimension=int(p_dim),
            q_dimension=int(hqq.shape[0]),
        )

    def effective_hamiltonian(
        self,
        energy_hartree: float,
        singular_floor: float = 1e-9,
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        denom = float(energy_hartree) - self.q_eigenvalues
        distance = float(np.min(np.abs(denom)))
        if distance < singular_floor:
            raise RuntimeError(
                f"Feshbach root resolvent gate failed: distance {distance} Ha below {singular_floor} Ha"
            )
        weighted = self.q_eigen_coupling / denom[:, None]
        sigma = self.q_eigen_coupling.T @ weighted
        sigma = 0.5 * (sigma + sigma.T)
        return self.hpp + sigma, {
            "distance_to_qhq_spectrum_hartree": distance,
            "self_energy_frobenius_hartree": float(np.linalg.norm(sigma)),
        }


def _orthonormal_seed(seed_subspace: np.ndarray, p_dim: int) -> np.ndarray:
    seed = np.asarray(seed_subspace, dtype=complex)
    if seed.ndim == 1:
        seed = seed[:, None]
    if seed.ndim != 2 or seed.shape[0] != p_dim or seed.shape[1] < 1:
        raise ValueError("seed_subspace must have shape (p_dim, d), d>=1")
    gram = seed.conj().T @ seed
    if not np.allclose(gram, np.eye(seed.shape[1]), atol=1e-10):
        raise ValueError("seed_subspace columns must be orthonormal")
    return seed


def tracked_effective_branch(
    kernel: FeshbachKernel,
    energy_hartree: float,
    seed_subspace: np.ndarray,
    singular_floor: float = 1e-9,
) -> Dict[str, object]:
    """Evaluate the H_eff branch with maximal overlap with a bare P seed class."""
    seed = _orthonormal_seed(seed_subspace, kernel.p_dimension)
    heff, diagnostics = kernel.effective_hamiltonian(
        energy_hartree,
        singular_floor=singular_floor,
    )
    vals, vecs = np.linalg.eigh(heff)
    weights = np.sum(np.abs(seed.conj().T @ vecs) ** 2, axis=0)
    d = seed.shape[1]
    chosen = np.argsort(weights)[::-1][:d]
    chosen_vals = np.sort(vals[chosen])
    chosen_weights = weights[chosen]
    mean_value = float(np.mean(chosen_vals))
    return {
        **diagnostics,
        "trial_energy_hartree": float(energy_hartree),
        "branch_energy_hartree": mean_value,
        "root_function_hartree": mean_value - float(energy_hartree),
        "branch_spread_hartree": float(np.max(chosen_vals) - np.min(chosen_vals)),
        "branch_overlap_mean": float(np.mean(chosen_weights)),
        "branch_overlap_min": float(np.min(chosen_weights)),
        "branch_dimension": int(d),
    }


def _find_bracket(
    kernel: FeshbachKernel,
    seed_subspace: np.ndarray,
    initial_energy_hartree: float,
    scan_halfwidths: Iterable[float],
    scan_points: int,
    singular_floor: float,
    overlap_floor: float,
):
    if scan_points < 5:
        raise ValueError("scan_points must be >=5")
    for halfwidth in scan_halfwidths:
        if halfwidth <= 0:
            continue
        grid = np.linspace(
            float(initial_energy_hartree) - float(halfwidth),
            float(initial_energy_hartree) + float(halfwidth),
            int(scan_points),
        )
        samples = []
        for energy in grid:
            try:
                out = tracked_effective_branch(
                    kernel,
                    float(energy),
                    seed_subspace,
                    singular_floor=singular_floor,
                )
            except RuntimeError:
                continue
            if out["branch_overlap_min"] < overlap_floor:
                continue
            samples.append((float(energy), float(out["root_function_hartree"]), out))
        candidates = []
        for left, right in zip(samples, samples[1:]):
            e0, f0, _ = left
            e1, f1, _ = right
            if f0 == 0.0:
                return (e0, e0, left[2], left[2])
            if f0 * f1 <= 0.0:
                center = 0.5 * (e0 + e1)
                candidates.append((abs(center - initial_energy_hartree), left, right))
        if candidates:
            _, left, right = min(candidates, key=lambda x: x[0])
            return (left[0], right[0], left[2], right[2])
    raise RuntimeError("Feshbach root solver could not bracket the tracked P-space branch")


def solve_tracked_feshbach_root(
    kernel: FeshbachKernel,
    seed_subspace: np.ndarray,
    initial_energy_hartree: float,
    *,
    scan_halfwidths: Tuple[float, ...] = (0.05, 0.10, 0.20, 0.40, 0.80),
    scan_points: int = 41,
    singular_floor: float = 1e-9,
    overlap_floor: float = 0.25,
    energy_tolerance_hartree: float = 1e-11,
    function_tolerance_hartree: float = 1e-11,
    max_iterations: int = 100,
) -> Dict[str, object]:
    """Solve lambda_c(H_eff(E))=E using only a bare P-space seed class.

    A deterministic scan around the bare P energy locates the nearest continuous
    high-overlap sign-changing branch. Bisection then solves the scalar nonlinear
    equation. No full P+Q eigenvalue is accepted as an input.
    """
    seed = _orthonormal_seed(seed_subspace, kernel.p_dimension)
    left, right, out_left, out_right = _find_bracket(
        kernel,
        seed,
        float(initial_energy_hartree),
        scan_halfwidths,
        scan_points,
        singular_floor,
        overlap_floor,
    )
    if left == right:
        result = out_left
        iterations = 0
    else:
        f_left = float(out_left["root_function_hartree"])
        f_right = float(out_right["root_function_hartree"])
        result = None
        iterations = 0
        for iterations in range(1, max_iterations + 1):
            middle = 0.5 * (left + right)
            out_mid = tracked_effective_branch(
                kernel,
                middle,
                seed,
                singular_floor=singular_floor,
            )
            if out_mid["branch_overlap_min"] < overlap_floor:
                raise RuntimeError("tracked Feshbach branch lost seed overlap during bisection")
            f_mid = float(out_mid["root_function_hartree"])
            result = out_mid
            if abs(f_mid) <= function_tolerance_hartree or abs(right - left) <= energy_tolerance_hartree:
                break
            if f_left * f_mid <= 0.0:
                right = middle
                f_right = f_mid
            else:
                left = middle
                f_left = f_mid
        else:
            raise RuntimeError("Feshbach root solver exceeded max_iterations")
        assert result is not None

    root_energy = float(result["trial_energy_hartree"])
    final = tracked_effective_branch(
        kernel,
        root_energy,
        seed,
        singular_floor=singular_floor,
    )
    return {
        **final,
        "root_energy_hartree": root_energy,
        "initial_energy_hartree": float(initial_energy_hartree),
        "energy_shift_hartree": root_energy - float(initial_energy_hartree),
        "iterations": int(iterations),
        "solver": "BRACKET_SCAN_PLUS_BISECTION",
    }
