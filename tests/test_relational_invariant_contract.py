import unittest

import numpy as np

from oes.invariant_contract import (
    INVARIANT_CONTRACT_ID,
    OES_ORBITAL_REPHASING_PROFILE_ID,
    rephase_orbital_values,
    rephase_transition_rdm,
)
from oes.phase_microscope import phase_microscope_field, transition_density_on_points
from oes.transition_contract import (
    ElectronicStateRef,
    NuclearCenter,
    OESTransitionState,
    TransitionProvenance,
)


class OESRelationalInvariantContractTests(unittest.TestCase):
    def _state(self, transition_rdm):
        return OESTransitionState(
            species_label="H2",
            nuclei=(
                NuclearCenter("H1", 1, 1.007825, (-0.7, 0.0, 0.0)),
                NuclearCenter("H2", 1, 1.007825, (0.7, 0.0, 0.0)),
            ),
            initial=ElectronicStateRef("X", -1.2, multiplicity=1, symmetry="Sigma_g"),
            final=ElectronicStateRef("A", -0.8, multiplicity=1, symmetry="Sigma_u"),
            transition_rdm=np.asarray(transition_rdm, dtype=complex),
            provenance=TransitionProvenance(
                source_repository="AdrianLipa90/Resonant-Chemistry",
                source_commit="fixture",
                method="unit-test",
                orbital_basis_id="fixture-mo-v1",
                backend_status="SIMULATED_REFERENCE",
            ),
        )

    def test_packet_emits_shared_contract_profile(self):
        payload = self._state([[0.0, 0.2j], [-0.2j, 0.0]]).as_dict()
        contract = payload["representation_contract"]
        self.assertEqual(contract["schema"], INVARIANT_CONTRACT_ID)
        self.assertEqual(contract["profile"], OES_ORBITAL_REPHASING_PROFILE_ID)

    def test_declared_unknown_contract_fails_closed(self):
        payload = self._state([[0.0, 0.2j], [-0.2j, 0.0]]).as_dict()
        payload["representation_contract"]["schema"] = "UNKNOWN"
        with self.assertRaisesRegex(ValueError, "unsupported relational invariant contract"):
            OESTransitionState.from_dict(payload)

    def test_legacy_packet_without_additive_contract_remains_readable(self):
        payload = self._state([[0.0, 0.2j], [-0.2j, 0.0]]).as_dict()
        payload.pop("representation_contract")
        restored = OESTransitionState.from_dict(payload)
        self.assertEqual(restored.n_spatial_orbitals, 2)

    def test_matched_orbital_rephasing_leaves_transition_density_invariant(self):
        t = np.asarray([[0.2, 0.4 + 0.1j], [-0.3j, 0.7]], dtype=complex)
        orbitals = np.asarray(
            [[1.0, 0.3j], [0.5 - 0.2j, 1.2], [-0.4j, 0.8 + 0.1j]],
            dtype=complex,
        )
        chi = np.asarray([0.37, -0.91])

        baseline = transition_density_on_points(self._state(t), orbitals)
        transformed = transition_density_on_points(
            self._state(rephase_transition_rdm(t, chi)),
            rephase_orbital_values(orbitals, chi),
        )
        np.testing.assert_allclose(transformed, baseline, atol=1.0e-12)

    def test_state_ray_global_phase_changes_raw_phase_not_relative_phase(self):
        t = np.asarray([[0.0, 1.0 + 0.2j], [0.3 - 0.1j, 0.0]], dtype=complex)
        orbitals = np.asarray(
            [[1.0, 0.4], [0.3j, 1.1], [0.7 - 0.2j, -0.5j]],
            dtype=complex,
        )
        field = phase_microscope_field(self._state(t), orbitals, amplitude_floor=0.0)

        global_phase = np.exp(1j * 0.73)
        shifted = phase_microscope_field(
            self._state(global_phase * t), orbitals, amplitude_floor=0.0
        )

        np.testing.assert_allclose(shifted.amplitude, field.amplitude, atol=1.0e-12)
        relative_a, defined_a, ref = field.relative_phase(amplitude_floor=0.0)
        relative_b, defined_b, ref_b = shifted.relative_phase(
            amplitude_floor=0.0, reference_index=ref
        )
        np.testing.assert_array_equal(defined_b, defined_a)
        self.assertEqual(ref_b, ref)
        np.testing.assert_allclose(
            np.exp(1j * relative_b), np.exp(1j * relative_a), atol=1.0e-12
        )

        # Transition amplitudes are globally phase-covariant; intensities are invariant.
        operator = np.asarray([[0.0, 2.0], [3.0, 0.0]], dtype=complex)
        mu = self._state(t).transition_amplitude(operator)
        mu_shifted = self._state(global_phase * t).transition_amplitude(operator)
        self.assertAlmostEqual(abs(mu_shifted) ** 2, abs(mu) ** 2, places=12)


if __name__ == "__main__":
    unittest.main()
