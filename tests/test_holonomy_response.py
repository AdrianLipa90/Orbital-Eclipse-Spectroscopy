import numpy as np
import pytest

from oes.holonomy_response import (
    HolonomyResponseError,
    fisher_information,
    identifiability_report,
    identifiable_rank,
    intensity_phase_slope,
    line_frequency_phase_slope,
    transition_frequency,
)


def test_transition_frequency_and_feynman_hellmann_slope():
    hbar = 2.0
    assert transition_frequency(7.0, 3.0, hbar) == pytest.approx(2.0)
    assert line_frequency_phase_slope(5.0, 1.0, hbar) == pytest.approx(-2.0)


def test_intensity_phase_slope():
    mu = 1.0 + 2.0j
    dmu = -0.5 + 0.25j
    expected = 2.0 * np.real(np.conjugate(mu) * dmu)
    assert intensity_phase_slope(mu, dmu) == pytest.approx(expected)


def test_full_rank_jacobian_is_locally_identifiable():
    j = np.array(
        [
            [1.0, 0.0],
            [0.0, 2.0],
            [1.0, 1.0],
        ]
    )
    sigma = np.diag([1.0, 2.0, 3.0])
    report = identifiability_report(j, sigma)
    assert report.jacobian_rank == 2
    assert report.fisher_rank == 2
    assert report.fully_locally_identifiable


def test_collinear_observable_rows_do_not_create_phase_rank():
    j = np.array(
        [
            [1.0, 2.0],
            [2.0, 4.0],
            [-3.0, -6.0],
        ]
    )
    assert identifiable_rank(j) == 1
    report = identifiability_report(j, np.eye(3))
    assert report.jacobian_rank == 1
    assert not report.fully_locally_identifiable


def test_invertible_phase_basis_change_preserves_rank():
    j = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, -1.0],
        ]
    )
    basis_change = np.array([[2.0, 1.0], [1.0, 1.0]])
    assert np.linalg.det(basis_change) != 0.0
    assert identifiable_rank(j) == identifiable_rank(j @ basis_change)


def test_fisher_matches_direct_inverse_definition():
    j = np.array([[1.0, 2.0], [3.0, 1.0]])
    sigma = np.array([[2.0, 0.3], [0.3, 1.5]])
    expected = j.T @ np.linalg.inv(sigma) @ j
    assert fisher_information(j, sigma) == pytest.approx(expected)


@pytest.mark.parametrize(
    "covariance",
    [
        np.array([[1.0, 2.0], [0.0, 1.0]]),
        np.array([[1.0, 0.0], [0.0, 0.0]]),
        np.array([[1.0, 2.0], [2.0, 1.0]]),
    ],
)
def test_covariance_fail_closed(covariance):
    j = np.eye(2)
    with pytest.raises(HolonomyResponseError):
        fisher_information(j, covariance)


@pytest.mark.parametrize(
    "call",
    [
        lambda: transition_frequency(1.0, 0.0, 0.0),
        lambda: line_frequency_phase_slope(1.0, 0.0, -1.0),
        lambda: identifiable_rank([[float("nan")]]),
    ],
)
def test_basic_domains_fail_closed(call):
    with pytest.raises(HolonomyResponseError):
        call()
