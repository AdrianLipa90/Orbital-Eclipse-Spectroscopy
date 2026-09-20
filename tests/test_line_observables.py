import numpy as np
import pytest

from oes.line_observables import electric_dipole_absorption_line
from oes.transition_contract import (
    ElectronicStateRef,
    NuclearCenter,
    OESTransitionState,
    TransitionProvenance,
)


def _transition(final_energy=-0.5, t=None):
    return OESTransitionState(
        species_label="X",
        nuclei=(NuclearCenter("X", 1, 1.0, (0.0, 0.0, 0.0)),),
        initial=ElectronicStateRef("i", -1.0),
        final=ElectronicStateRef("f", final_energy),
        transition_rdm=np.asarray([[0.0, 0.5], [0.0, 0.0]] if t is None else t, dtype=complex),
        provenance=TransitionProvenance("fixture", "fixture", "fixture", "basis-v1", "SIMULATED_REFERENCE"),
    )


def test_electric_dipole_line_derives_frequency_wavelength_and_strength():
    transition = _transition()
    ops = np.zeros((3, 2, 2), dtype=complex)
    ops[0, 0, 1] = 2.0
    line = electric_dipole_absorption_line(transition, ops)
    assert line.transition_dipole_au[0] == pytest.approx(1.0)
    assert line.transition_dipole_norm_au == pytest.approx(1.0)
    assert line.oscillator_strength_length_gauge == pytest.approx((2.0 / 3.0) * 0.5)
    assert line.frequency_hz > 0.0
    assert line.wavelength_nm > 0.0


def test_zero_transition_density_is_dark():
    transition = _transition(t=np.zeros((2, 2), dtype=complex))
    line = electric_dipole_absorption_line(transition, np.ones((3, 2, 2)))
    assert line.transition_dipole_norm_au == 0.0
    assert line.oscillator_strength_length_gauge == 0.0


def test_absorption_line_fails_closed_for_downward_state_order():
    transition = _transition(final_energy=-1.5)
    with pytest.raises(ValueError, match="final energy above initial"):
        electric_dipole_absorption_line(transition, np.ones((3, 2, 2)))
