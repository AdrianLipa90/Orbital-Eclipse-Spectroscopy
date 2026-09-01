import unittest

from oes.quantum.helium_q1 import run_helium_q1


class TestQuantumHeliumQ1(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result, cls.states = run_helium_q1()

    def test_canonical_register(self):
        self.assertEqual(self.result.n_spatial_orbitals, 10)
        self.assertEqual(self.result.n_spin_orbitals, 20)
        self.assertEqual(self.result.full_qubit_dimension, 1_048_576)
        self.assertEqual(self.result.fixed_particle_dimension, 190)
        self.assertEqual(self.result.backend, "SIMULATED_REFERENCE")

    def test_oes_sector_fci_matches_pyscf_fci(self):
        self.assertLess(abs(self.result.fci_delta_hartree), 1e-9)
        self.assertLess(self.result.oes_fci_energy_hartree, self.result.rhf_energy_hartree)

    def test_spin_sectors_resolved(self):
        self.assertLess(abs(self.result.ground_s2), 1e-8)
        triplet = self.states[self.result.first_triplet_index]
        singlet = self.states[self.result.first_singlet_excited_index]
        self.assertLess(abs(triplet.s2 - 2.0), 1e-7)
        self.assertLess(abs(singlet.s2), 1e-7)

    def test_spin_summed_transition_density_selection(self):
        # A spin-independent one-body transition density cannot connect the
        # singlet ground state to a triplet, while a singlet-singlet channel
        # should remain visible in the active space.
        self.assertLess(self.result.singlet_triplet_spatial_transition_norm, 1e-8)
        self.assertGreater(self.result.singlet_singlet_spatial_transition_norm, 1e-5)


if __name__ == "__main__":
    unittest.main()
