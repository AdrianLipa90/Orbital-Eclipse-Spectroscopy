import unittest

from oes.block_response import (
    BlockResponseError,
    block_spectral_observables,
)
from oes.molecular_tda import MolecularTDAResult, MolecularTDAState


def _state(root, energy_ev, oscillator, mu_norm):
    return MolecularTDAState(
        root=root,
        excitation_hartree=energy_ev / 27.211386245981,
        excitation_ev=energy_ev,
        transition_dipole_au=(mu_norm, 0.0, 0.0),
        transition_dipole_norm_au=mu_norm,
        oscillator_strength_backend=oscillator,
        oscillator_strength_oes=oscillator,
        oscillator_strength_delta=0.0,
    )


def _result(states):
    return MolecularTDAResult(
        backend="FIXTURE",
        species_label="fixture",
        basis_name="fixture",
        charge=0,
        spin=0,
        atoms=("H", "H"),
        coordinates_angstrom=((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        rhf_energy_hartree=-1.0,
        states=tuple(states),
        max_oscillator_strength_delta=0.0,
    )


class BlockResponseTests(unittest.TestCase):
    def test_block_observables_are_permutation_invariant(self):
        result = _result(
            (
                _state(1, 4.0, 0.2, 0.4),
                _state(2, 5.0, 0.3, 0.5),
                _state(3, 8.0, 0.01, 0.1),
            )
        )
        a = block_spectral_observables(result, (1, 2))
        b = block_spectral_observables(result, (2, 1))
        self.assertAlmostEqual(
            a.oscillator_strength_sum,
            b.oscillator_strength_sum,
            places=15,
        )
        self.assertAlmostEqual(
            a.excitation_centroid_ev,
            b.excitation_centroid_ev,
            places=15,
        )
        self.assertAlmostEqual(
            a.oscillator_weighted_centroid_ev,
            b.oscillator_weighted_centroid_ev,
            places=15,
        )

    def test_weighted_centroid_matches_declared_strengths(self):
        result = _result(
            (
                _state(1, 4.0, 1.0, 0.4),
                _state(2, 8.0, 3.0, 0.5),
            )
        )
        block = block_spectral_observables(result, (1, 2))
        self.assertAlmostEqual(
            block.oscillator_weighted_centroid_ev,
            7.0,
            places=15,
        )
        self.assertAlmostEqual(
            block.oscillator_strength_sum,
            4.0,
            places=15,
        )
        self.assertAlmostEqual(
            block.oscillator_strength_fraction_of_computed_window,
            1.0,
            places=15,
        )

    def test_dark_block_has_no_weighted_centroid(self):
        result = _result(
            (
                _state(1, 4.0, 0.0, 0.0),
                _state(2, 5.0, 0.0, 0.0),
            )
        )
        block = block_spectral_observables(result, (1, 2))
        self.assertIsNone(block.oscillator_weighted_centroid_ev)
        self.assertEqual(block.oscillator_strength_sum, 0.0)

    def test_bad_root_selector_fails_closed(self):
        result = _result((_state(1, 4.0, 0.2, 0.4),))
        with self.assertRaises(BlockResponseError):
            block_spectral_observables(result, (2,))


if __name__ == "__main__":
    unittest.main()
