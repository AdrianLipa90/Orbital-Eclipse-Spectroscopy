import math
import unittest

from oes.hydrogen import (
    bohr_radius_hydrogen_m,
    contact_exposure_dimensionless,
    dipole_allowed,
    gross_energy_ev,
    gross_transition,
    lande_g,
    m_resolved_flavor_count,
    orbital_flavor_count,
    radial_node_count,
    zeeman_transition_shift_hz,
)


class HydrogenClosureTests(unittest.TestCase):
    def test_reduced_mass_bohr_radius(self):
        self.assertAlmostEqual(bohr_radius_hydrogen_m() * 1e12, 52.94654091, places=6)

    def test_ground_energy_no_spectral_input(self):
        self.assertAlmostEqual(gross_energy_ev(1), -13.59828728, places=7)

    def test_lyman_alpha_gross_color(self):
        line = gross_transition(2, 1)
        self.assertAlmostEqual(line.delta_e_ev, 10.19871546, places=7)
        self.assertAlmostEqual(line.frequency_hz / 1e15, 2.466038427, places=8)
        self.assertAlmostEqual(line.wavelength_nm, 121.56844545, places=7)

    def test_balmer_alpha_gross_color(self):
        line = gross_transition(3, 2)
        self.assertAlmostEqual(line.wavelength_nm, 656.469606, places=5)

    def test_color_and_e1_flavor_are_independent(self):
        # 2p -> 1s is E1-open; 2s -> 1s has the same gross n=2 -> 1 color
        # but its orbital E1 gate is closed.
        self.assertTrue(dipole_allowed(1, 0, 0, 0))
        self.assertFalse(dipole_allowed(0, 0, 0, 0))

    def test_contact_exposure(self):
        self.assertEqual(contact_exposure_dimensionless(2, 1), 0.0)
        self.assertAlmostEqual(contact_exposure_dimensionless(2, 0), 1 / 8)
        self.assertAlmostEqual(contact_exposure_dimensionless(3, 0), 1 / 27)

    def test_radial_node_code(self):
        self.assertEqual(radial_node_count(4, 0), 3)
        self.assertEqual(radial_node_count(4, 1), 2)
        self.assertEqual(radial_node_count(4, 2), 1)
        self.assertEqual(radial_node_count(4, 3), 0)

    def test_orbital_flavor_count(self):
        self.assertEqual([orbital_flavor_count(n) for n in (1, 2, 3, 4)], [1, 3, 5, 7])

    def test_m_resolved_flavor_count(self):
        self.assertEqual([m_resolved_flavor_count(n) for n in (1, 2, 3)], [3, 15, 39])

    def test_lande_factors(self):
        self.assertAlmostEqual(lande_g(0, 0.5), 2.0)
        self.assertAlmostEqual(lande_g(1, 0.5), 2 / 3)
        self.assertAlmostEqual(lande_g(1, 1.5), 4 / 3)

    def test_linear_zeeman_flavor_to_color(self):
        # 2p_1/2(m=+1/2) -> 1s_1/2(m=-1/2) at B=1 T.
        # Shift = (2/3*1/2 - 2*(-1/2))*mu_B/h = 4/3*mu_B/h.
        shift = zeeman_transition_shift_hz(
            b_tesla=1.0,
            l_i=1,
            j_i=0.5,
            m_j_i=0.5,
            l_f=0,
            j_f=0.5,
            m_j_f=-0.5,
        )
        self.assertAlmostEqual(shift / 1e9, (4 / 3) * 13.99624555, places=7)


if __name__ == "__main__":
    unittest.main()
