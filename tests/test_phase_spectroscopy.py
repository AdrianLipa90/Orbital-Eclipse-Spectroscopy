import math
import numpy as np
from oes import phase_spectroscopy as ps


def test_unfolding_has_unit_mean_gap():
    u = ps.mean_spacing_unfold([10.0, 11.0, 14.0, 18.0])
    assert np.isclose(np.mean(np.diff(u)), 1.0)


def test_forced_phase_repulsion_and_legacy_alias():
    x = [0.0, 2.0 * math.pi]
    expected = [0.0, 1.0]
    assert np.allclose(ps.forced_pair_correlation_phase(x), expected, atol=1e-14)
    assert np.allclose(ps.gue_pair_correlation_phase(x), expected, atol=1e-14)


def test_form_factor_global_phase_invariance_and_normalization():
    p = np.array([0.2, 1.3, 2.4, 5.1])
    m = np.array([0.25, 0.5, 1.0, 2.0])
    raw = ps.spectral_form_factor(p, m)
    shifted = ps.spectral_form_factor(p + 0.73, m)
    assert np.allclose(raw, shifted)
    assert np.allclose(ps.normalized_spectral_form_factor(p, m), p.size * raw)


def test_forced_form_factor_target_is_ramp_plateau():
    m = np.array([-1.5, -1.0, -0.5, 0.0, 0.25, 1.0, 1.5])
    assert np.allclose(
        ps.forced_form_factor_target(m),
        [1.0, 1.0, 0.5, 0.0, 0.25, 1.0, 1.0],
    )


def test_prime_power_channels():
    rows = {r[2]: r for r in ps.prime_power_carriers(8)}
    assert {2, 3, 4, 5, 7, 8}.issubset(rows)
    assert math.isclose(rows[8][3], math.log(8.0), abs_tol=1e-14)


def test_phase_scramble_reverse_control():
    t = np.linspace(14.0, 40.0, 64)
    a = ps.prime_power_signal(t, cutoff=100, sigma=0.2)
    b = ps.phase_scrambled_prime_power_signal(t, cutoff=100, sigma=0.2, seed=7)
    assert not np.allclose(a, b)


def test_riemann_phase_coordinates_monotone():
    theta = ps.riemann_zero_phase_coordinates([14.134725141734693, 21.022039638771555, 25.01085758014569])
    assert np.all(np.diff(theta) > 0.0)
