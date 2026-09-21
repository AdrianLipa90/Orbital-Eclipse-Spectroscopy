"""Typed cross-repository transition contract for Orbital Eclipse Spectroscopy.

Resonant-Chemistry owns state generation. OES owns transitions. This module
defines the boundary without importing any chemistry solver.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.constants import h, physical_constants

from .invariant_contract import (
    representation_contract_block,
    validate_representation_contract_block,
)

SCHEMA_ID = "OES_TRANSITION_STATE_V0_1"
TRANSITION_RDM_CONVENTION = "T_pq=<Psi_f|a_p^dagger a_q|Psi_i>"
HARTREE_J = physical_constants["Hartree energy"][0]


def _finite_float(value: float, field: str) -> float:
    x = float(value)
    if not math.isfinite(x):
        raise ValueError(f"{field} must be finite")
    return x


@dataclass(frozen=True)
class NuclearCenter:
    label: str
    atomic_number: int
    mass_u: float
    position_bohr: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not self.label:
            raise ValueError("nuclear center label is required")
        if isinstance(self.atomic_number, bool) or int(self.atomic_number) != self.atomic_number or self.atomic_number <= 0:
            raise ValueError("atomic_number must be a positive integer")
        if _finite_float(self.mass_u, "mass_u") <= 0.0:
            raise ValueError("mass_u must be positive")
        if len(self.position_bohr) != 3:
            raise ValueError("position_bohr must have exactly three coordinates")
        for i, value in enumerate(self.position_bohr):
            _finite_float(value, f"position_bohr[{i}]")

    def as_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "atomic_number": int(self.atomic_number),
            "mass_u": float(self.mass_u),
            "position_bohr": [float(x) for x in self.position_bohr],
        }


@dataclass(frozen=True)
class ElectronicStateRef:
    label: str
    energy_hartree: float
    multiplicity: int | None = None
    symmetry: str | None = None

    def __post_init__(self) -> None:
        if not self.label:
            raise ValueError("electronic state label is required")
        _finite_float(self.energy_hartree, "energy_hartree")
        if self.multiplicity is not None:
            if isinstance(self.multiplicity, bool) or int(self.multiplicity) != self.multiplicity or self.multiplicity < 1:
                raise ValueError("multiplicity must be a positive integer when supplied")

    def as_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "energy_hartree": float(self.energy_hartree),
            "multiplicity": self.multiplicity,
            "symmetry": self.symmetry,
        }


@dataclass(frozen=True)
class TransitionProvenance:
    source_repository: str
    source_commit: str
    method: str
    orbital_basis_id: str
    backend_status: str

    def __post_init__(self) -> None:
        for field, value in (
            ("source_repository", self.source_repository),
            ("source_commit", self.source_commit),
            ("method", self.method),
            ("orbital_basis_id", self.orbital_basis_id),
            ("backend_status", self.backend_status),
        ):
            if not value:
                raise ValueError(f"{field} is required")

    def as_dict(self) -> dict[str, str]:
        return {
            "source_repository": self.source_repository,
            "source_commit": self.source_commit,
            "method": self.method,
            "orbital_basis_id": self.orbital_basis_id,
            "backend_status": self.backend_status,
        }


@dataclass(frozen=True)
class OESTransitionState:
    species_label: str
    nuclei: tuple[NuclearCenter, ...]
    initial: ElectronicStateRef
    final: ElectronicStateRef
    transition_rdm: np.ndarray
    provenance: TransitionProvenance
    transition_kind: str = "electronic"

    def __post_init__(self) -> None:
        if not self.species_label:
            raise ValueError("species_label is required")
        if not self.nuclei:
            raise ValueError("at least one nuclear center is required")
        if not self.transition_kind:
            raise ValueError("transition_kind is required")
        t = np.asarray(self.transition_rdm, dtype=complex)
        if t.ndim != 2 or t.shape[0] != t.shape[1] or t.shape[0] < 1:
            raise ValueError("transition_rdm must be a non-empty square matrix")
        if not np.all(np.isfinite(t.real)) or not np.all(np.isfinite(t.imag)):
            raise ValueError("transition_rdm must be finite")
        object.__setattr__(self, "transition_rdm", t.copy())

    @property
    def n_spatial_orbitals(self) -> int:
        return int(self.transition_rdm.shape[0])

    @property
    def signed_energy_gap_hartree(self) -> float:
        return float(self.final.energy_hartree - self.initial.energy_hartree)

    @property
    def frequency_hz(self) -> float:
        return abs(self.signed_energy_gap_hartree) * HARTREE_J / h

    def transition_amplitude(self, one_body_operator: np.ndarray) -> complex:
        operator = np.asarray(one_body_operator, dtype=complex)
        if operator.shape != self.transition_rdm.shape:
            raise ValueError("one_body_operator shape must match transition_rdm")
        if not np.all(np.isfinite(operator.real)) or not np.all(np.isfinite(operator.imag)):
            raise ValueError("one_body_operator must be finite")
        return complex(np.einsum("pq,pq->", operator, self.transition_rdm))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema": SCHEMA_ID,
            "species_label": self.species_label,
            "transition_kind": self.transition_kind,
            "nuclear_framework": [center.as_dict() for center in self.nuclei],
            "initial_state": self.initial.as_dict(),
            "final_state": self.final.as_dict(),
            "transition_rdm": {
                "convention": TRANSITION_RDM_CONVENTION,
                "real": self.transition_rdm.real.tolist(),
                "imag": self.transition_rdm.imag.tolist(),
            },
            "representation_contract": representation_contract_block(),
            "derived": {
                "signed_energy_gap_hartree": self.signed_energy_gap_hartree,
                "frequency_hz": self.frequency_hz,
            },
            "provenance": self.provenance.as_dict(),
            "epistemic_status": "STANDARD_QM_TRANSITION_CONTRACT",
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "OESTransitionState":
        if payload.get("schema") != SCHEMA_ID:
            raise ValueError(f"unsupported transition schema: {payload.get('schema')!r}")
        rdm = payload.get("transition_rdm")
        if not isinstance(rdm, Mapping):
            raise ValueError("transition_rdm object is required")
        if rdm.get("convention") != TRANSITION_RDM_CONVENTION:
            raise ValueError("transition_rdm convention mismatch")
        real = np.asarray(rdm.get("real"), dtype=float)
        imag = np.asarray(rdm.get("imag"), dtype=float)
        if real.shape != imag.shape:
            raise ValueError("transition_rdm real/imag shapes must match")
        initial = payload.get("initial_state")
        final = payload.get("final_state")
        provenance = payload.get("provenance")
        representation_contract = payload.get("representation_contract")
        if representation_contract is not None:
            if not isinstance(representation_contract, Mapping):
                raise ValueError("representation_contract must be an object")
            validate_representation_contract_block(representation_contract)
        nuclei = payload.get("nuclear_framework")
        if not isinstance(initial, Mapping) or not isinstance(final, Mapping):
            raise ValueError("initial_state and final_state objects are required")
        if not isinstance(provenance, Mapping):
            raise ValueError("provenance object is required")
        if not isinstance(nuclei, Sequence) or isinstance(nuclei, (str, bytes)):
            raise ValueError("nuclear_framework array is required")
        return cls(
            species_label=str(payload.get("species_label", "")),
            transition_kind=str(payload.get("transition_kind", "")),
            nuclei=tuple(
                NuclearCenter(
                    label=str(row["label"]),
                    atomic_number=int(row["atomic_number"]),
                    mass_u=float(row["mass_u"]),
                    position_bohr=tuple(float(x) for x in row["position_bohr"]),
                )
                for row in nuclei
            ),
            initial=ElectronicStateRef(
                label=str(initial["label"]),
                energy_hartree=float(initial["energy_hartree"]),
                multiplicity=None if initial.get("multiplicity") is None else int(initial["multiplicity"]),
                symmetry=None if initial.get("symmetry") is None else str(initial["symmetry"]),
            ),
            final=ElectronicStateRef(
                label=str(final["label"]),
                energy_hartree=float(final["energy_hartree"]),
                multiplicity=None if final.get("multiplicity") is None else int(final["multiplicity"]),
                symmetry=None if final.get("symmetry") is None else str(final["symmetry"]),
            ),
            transition_rdm=real + 1j * imag,
            provenance=TransitionProvenance(
                source_repository=str(provenance["source_repository"]),
                source_commit=str(provenance["source_commit"]),
                method=str(provenance["method"]),
                orbital_basis_id=str(provenance["orbital_basis_id"]),
                backend_status=str(provenance["backend_status"]),
            ),
        )


__all__ = [
    "SCHEMA_ID",
    "TRANSITION_RDM_CONVENTION",
    "NuclearCenter",
    "ElectronicStateRef",
    "TransitionProvenance",
    "OESTransitionState",
]
