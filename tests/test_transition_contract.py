import unittest

import numpy as np

from oes.phase_microscope import phase_microscope_field, transition_density_on_points
from oes.transition_contract import (
    ElectronicStateRef,
    NuclearCenter,
    OESTransitionState,
    TransitionProvenance,
)


class OESTransitionContractTests(unittest.TestCase):
    def _fixture(self, t=None):
        matrix = np.asarray(
            [[0.0, 0.25j], [-0.25j, 0.0]] if t is None else t,
            dtype=complex,
        )
        return OESTransitionState(
            species_label="H2",
            nuclei=(
                NuclearCenter("H1", 1, 1.007825, (-0.7, 0.0, 0.0)),
                NuclearCenter("H2", 1, 1.007825, (0.7, 0.0, 0.0)),
            ),
            initial=ElectronicStateRef("X", -1.2, multiplicity=1, symmetry="Sigma_g"),
            final=ElectronicStateRef("A", -0.8, multiplicity=1, symmetry="Sigma_u"),
            transition_rdm=matrix,
            provenance=TransitionProvenance(
                source_repository="AdrianLipa90/Resonant-Chemistry",
                source_commit="deadbeef",
                method="unit-fixture",
                orbital_basis_id="fixture-mo-v1",
                backend_status="SIMULATED_REFERENCE",
            ),
        )

    def test_contract_round_trip_preserves_complex_transition_rdm(self):
        state = self._fixture()
        restored = OESTransitionState.from_dict(state.as_dict())
        self.assertEqual(restored.species_label, state.species_label)
        self.assertEqual(restored.n_spatial_orbitals, 2)
        self.assertTrue(np.allclose(restored.transition_rdm, state.transition_rdm))
        self.assertGreater(restored.frequency_hz, 0.0)

    def test_one_body_transition_amplitude_uses_declared_pq_convention(self):
        state = self._fixture()
        operator = np.asarray([[0.0, 2.0], [3.0, 0.0]], dtype=complex)
        expected = 2.0 * 0.25j + 3.0 * (-0.25j)
        observed = state.transition_amplitude(operator)
        self.assertLess(abs(observed - expected), 1.0e-12)

    def test_contract_fails_closed_on_convention_mismatch(self):
        payload = self._fixture().as_dict()
        payload["transition_rdm"]["convention"] = "ambiguous"
        with self.assertRaisesRegex(ValueError, "convention mismatch"):
            OESTransitionState.from_dict(payload)

    def test_phase_field_reconstructs_transition_density_where_defined(self):
        transition = self._fixture([[0.0, 1.0j], [0.5, 0.0]])
        orbitals = np.asarray([[1.0, 2.0], [1.0j, 0.5]], dtype=complex)
        field = phase_microscope_field(transition, orbitals, amplitude_floor=0.0)
        self.assertTrue(np.all(field.phase_defined))
        self.assertTrue(np.allclose(field.reconstruct(), field.transition_density))

    def test_transition_density_is_invariant_under_matched_orbital_rephasing(self):
        t = np.asarray([[0.2, 0.4 + 0.1j], [-0.3j, 0.7]], dtype=complex)
        transition = self._fixture(t)
        orbitals = np.asarray(
            [[1.0, 0.3j], [0.5 - 0.2j, 1.2], [-0.4j, 0.8 + 0.1j]],
            dtype=complex,
        )
        baseline = transition_density_on_points(transition, orbitals)

        theta = np.asarray([0.37, -0.91])
        phases = np.exp(1j * theta)
        rotated_orbitals = orbitals * phases[None, :]
        rotated_t = phases[:, None] * np.conjugate(phases[None, :]) * t
        rotated_transition = self._fixture(rotated_t)

        observed = transition_density_on_points(rotated_transition, rotated_orbitals)
        self.assertTrue(np.allclose(observed, baseline, atol=1e-12))


if __name__ == "__main__":
    unittest.main()
