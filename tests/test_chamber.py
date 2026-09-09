"""Independent verification of the quasi-steady chamber model.

Expected values are formed from the defining relation ``c* = p_c A_t / m_dot``
rearranged by hand, from hand-computed literals, and from scaling arguments that
do not reuse the production expression.
"""

from __future__ import annotations

import math

import pytest

from hybrid_rocket_motor.chamber import (
    ILLUSTRATIVE_N2O_HTPB_COMBUSTION,
    UNIVERSAL_GAS_CONSTANT_J_MOL_K,
    ChamberStatus,
    CombustionProperties,
    ideal_characteristic_velocity,
    solve_chamber,
)
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.regression import Provenance

COMBUSTION = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
A_T = REFERENCE_NOZZLE.throat_area_m2

M_DOT_OX = 0.100
M_DOT_F = 0.039776644


# --- 1. total mass flow ----------------------------------------------------


def test_total_mass_flow_is_the_sum_of_the_parts():
    state = solve_chamber(COMBUSTION, A_T, M_DOT_OX, M_DOT_F)
    assert state.total_mass_flow_kg_s == pytest.approx(M_DOT_OX + M_DOT_F, rel=1e-15)
    assert state.total_mass_flow_kg_s == pytest.approx(0.139776644, rel=1e-15)


def test_mixture_ratio_matches_the_milestone_1_value():
    state = solve_chamber(COMBUSTION, A_T, M_DOT_OX, M_DOT_F)
    assert state.mixture_ratio == pytest.approx(M_DOT_OX / M_DOT_F, rel=1e-15)
    assert state.mixture_ratio == pytest.approx(2.514038, rel=1e-6)


# --- 2, 3. chamber pressure and the exact c* identity ---------------------


def test_chamber_pressure_matches_the_hand_rearranged_definition():
    """p_c = m_dot c* / A_t, formed independently here."""
    state = solve_chamber(COMBUSTION, A_T, M_DOT_OX, M_DOT_F)
    expected = (M_DOT_OX + M_DOT_F) * COMBUSTION.c_star_m_s / A_T
    assert state.chamber_pressure_pa == pytest.approx(expected, rel=1e-15)


def test_c_star_identity_holds_exactly():
    """c* is *defined* by p_c A_t = m_dot c*; this must close to machine zero."""
    for m_ox, m_f in ((0.05, 0.02), (M_DOT_OX, M_DOT_F), (0.2, 0.06)):
        state = solve_chamber(COMBUSTION, A_T, m_ox, m_f)
        residual = (
            state.chamber_pressure_pa * state.throat_area_m2
            - state.total_mass_flow_kg_s * state.c_star_m_s
        )
        assert abs(residual) <= 1.0e-9
        assert state.chamber_pressure_pa * state.throat_area_m2 / state.total_mass_flow_kg_s == (
            pytest.approx(COMBUSTION.c_star_m_s, rel=1e-14)
        )


def test_reference_chamber_pressure_hardcoded_literal():
    """0.139776644 * 1467.3469189547673 / 7.853981633974483e-5 = 2611424.85 Pa."""
    state = solve_chamber(COMBUSTION, A_T, M_DOT_OX, M_DOT_F)
    assert state.chamber_pressure_pa == pytest.approx(2611424.85, rel=1e-8)
    assert state.chamber_pressure_pa / 1e5 == pytest.approx(26.114, rel=1e-4)


# --- 4, 5, 6. scaling behaviour -------------------------------------------


def test_pressure_scales_linearly_with_total_mass_flow():
    single = solve_chamber(COMBUSTION, A_T, 0.05, 0.02).chamber_pressure_pa
    double = solve_chamber(COMBUSTION, A_T, 0.10, 0.04).chamber_pressure_pa
    assert double == pytest.approx(2.0 * single, rel=1e-14)


def test_pressure_scales_linearly_with_c_star():
    doubled = COMBUSTION.replace(c_star_m_s=2.0 * COMBUSTION.c_star_m_s)
    base = solve_chamber(COMBUSTION, A_T, M_DOT_OX, M_DOT_F).chamber_pressure_pa
    scaled = solve_chamber(doubled, A_T, M_DOT_OX, M_DOT_F).chamber_pressure_pa
    assert scaled == pytest.approx(2.0 * base, rel=1e-14)


def test_pressure_scales_inversely_with_throat_area():
    base = solve_chamber(COMBUSTION, A_T, M_DOT_OX, M_DOT_F).chamber_pressure_pa
    wide = solve_chamber(COMBUSTION, 2.0 * A_T, M_DOT_OX, M_DOT_F).chamber_pressure_pa
    assert wide == pytest.approx(base / 2.0, rel=1e-14)


def test_pressure_is_unchanged_when_flow_is_shifted_between_oxidizer_and_fuel():
    """Only the total matters to the chamber closure."""
    a = solve_chamber(COMBUSTION, A_T, 0.10, 0.04).chamber_pressure_pa
    b = solve_chamber(COMBUSTION, A_T, 0.04, 0.10).chamber_pressure_pa
    assert a == pytest.approx(b, rel=1e-15)


# --- 7. zero-flow idle handling -------------------------------------------


def test_zero_total_flow_gives_an_idle_chamber():
    state = solve_chamber(COMBUSTION, A_T, 0.0, 0.0)
    assert state.status is ChamberStatus.IDLE_NO_FLOW
    assert not state.is_firing
    assert state.chamber_pressure_pa == 0.0
    assert state.total_mass_flow_kg_s == 0.0
    assert math.isnan(state.mixture_ratio)


def test_oxidizer_only_flow_still_fires_but_has_undefined_mixture_ratio():
    state = solve_chamber(COMBUSTION, A_T, 0.05, 0.0)
    assert state.status is ChamberStatus.FIRING
    assert state.chamber_pressure_pa > 0.0
    assert math.isnan(state.mixture_ratio)


def test_tiny_flow_still_fires():
    state = solve_chamber(COMBUSTION, A_T, 1.0e-9, 0.0)
    assert state.status is ChamberStatus.FIRING
    assert state.chamber_pressure_pa > 0.0


# --- 8. invalid input rejection -------------------------------------------


@pytest.mark.parametrize("bad", [-1.0e-9, -0.1, math.nan, math.inf])
def test_negative_or_non_finite_flows_are_rejected(bad):
    with pytest.raises(ValueError):
        solve_chamber(COMBUSTION, A_T, bad, M_DOT_F)
    with pytest.raises(ValueError):
        solve_chamber(COMBUSTION, A_T, M_DOT_OX, bad)


@pytest.mark.parametrize("bad", [0.0, -1.0e-5, math.nan, math.inf])
def test_invalid_throat_area_is_rejected(bad):
    with pytest.raises(ValueError):
        solve_chamber(COMBUSTION, bad, M_DOT_OX, M_DOT_F)


@pytest.mark.parametrize("bad", [0.0, -1500.0, math.nan, math.inf])
def test_invalid_c_star_is_rejected(bad):
    with pytest.raises(ValueError):
        CombustionProperties(
            c_star_m_s=bad,
            gamma=1.2,
            chamber_temperature_k=2600.0,
            molar_mass_kg_mol=0.022,
            label="bad",
            reference="invalid-input test",
        )


@pytest.mark.parametrize("bad_gamma", [1.0, 0.8, 0.0, -1.2, math.nan, math.inf])
def test_invalid_gamma_is_rejected(bad_gamma):
    with pytest.raises(ValueError):
        CombustionProperties(
            c_star_m_s=1500.0,
            gamma=bad_gamma,
            chamber_temperature_k=2600.0,
            molar_mass_kg_mol=0.022,
            label="bad",
            reference="invalid-input test",
        )


@pytest.mark.parametrize("bad", [0.0, -2600.0, math.nan, math.inf])
def test_invalid_temperature_or_molar_mass_is_rejected(bad):
    with pytest.raises(ValueError):
        CombustionProperties(
            c_star_m_s=1500.0,
            gamma=1.2,
            chamber_temperature_k=bad,
            molar_mass_kg_mol=0.022,
            label="bad",
            reference="invalid-input test",
        )
    with pytest.raises(ValueError):
        CombustionProperties(
            c_star_m_s=1500.0,
            gamma=1.2,
            chamber_temperature_k=2600.0,
            molar_mass_kg_mol=bad,
            label="bad",
            reference="invalid-input test",
        )


@pytest.mark.parametrize("bad_eta", [0.0, -0.5, 1.5, math.nan, math.inf])
def test_invalid_c_star_efficiency_is_rejected(bad_eta):
    with pytest.raises(ValueError):
        CombustionProperties.from_ideal_with_efficiency(
            1.2, 2600.0, 0.022, bad_eta, "bad", "invalid-input test"
        )


def test_replace_rejects_unknown_fields():
    with pytest.raises(ValueError):
        COMBUSTION.replace(not_a_field=1.0)


def test_combustion_properties_are_immutable():
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        COMBUSTION.gamma = 1.3


# --- prescribed-property bookkeeping --------------------------------------


def test_gas_constant_is_the_universal_constant_over_the_molar_mass():
    assert COMBUSTION.gas_constant_j_kg_k == pytest.approx(
        UNIVERSAL_GAS_CONSTANT_J_MOL_K / 0.022, rel=1e-15
    )
    assert COMBUSTION.gas_constant_j_kg_k == pytest.approx(377.930119, rel=1e-12)


def test_ideal_c_star_matches_the_hand_written_relation():
    """c*_ideal = sqrt(R T_c / g) [(g+1)/2]^((g+1)/(2(g-1)))."""
    for gamma, r_gas, t_c in ((1.20, 377.93, 2600.0), (1.25, 300.0, 3000.0)):
        exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
        expected = math.sqrt(r_gas * t_c / gamma) * ((gamma + 1.0) / 2.0) ** exponent
        assert ideal_characteristic_velocity(gamma, r_gas, t_c) == pytest.approx(
            expected, rel=1e-14
        )


def test_ideal_c_star_is_recovered_from_the_choked_mass_flow_relation():
    """Independent route: build m_dot from the NASA Glenn choked relation.

    m_dot = (A p_t / sqrt(T_t)) sqrt(g/R) [(g+1)/2]^(-(g+1)/(2(g-1)))
    so c* = p_t A / m_dot must reproduce the closed form.
    """
    gamma, r_gas, t_c, p_c, area = 1.20, 377.93, 2600.0, 2.5e6, 7.853981633974483e-5
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    m_dot = (area * p_c / math.sqrt(t_c)) * math.sqrt(gamma / r_gas) * (
        ((gamma + 1.0) / 2.0) ** (-exponent)
    )
    assert p_c * area / m_dot == pytest.approx(
        ideal_characteristic_velocity(gamma, r_gas, t_c), rel=1e-13
    )


def test_reference_combustion_properties_are_self_consistent():
    assert COMBUSTION.c_star_efficiency == pytest.approx(0.96, rel=1e-12)
    assert COMBUSTION.c_star_efficiency < 1.0
    assert COMBUSTION.c_star_m_s < COMBUSTION.ideal_c_star_m_s


def test_reference_combustion_property_literals():
    """gamma 1.20, T_c 2600 K, M 0.022 kg/mol, eta 0.96 -> c* = 1467.35 m/s."""
    assert COMBUSTION.gamma == 1.20
    assert COMBUSTION.chamber_temperature_k == 2600.0
    assert COMBUSTION.molar_mass_kg_mol == 0.022
    assert COMBUSTION.ideal_c_star_m_s == pytest.approx(1528.486, rel=1e-6)
    assert COMBUSTION.c_star_m_s == pytest.approx(1467.347, rel=1e-6)


def test_prescribed_c_star_lies_inside_the_range_measured_by_the_source_study():
    """Rezaei et al. (2018) measured c* = 1403-1587 m/s for HTPB/N2O.

    A weak consistency check on a PRESCRIBED value, not a validation: the four
    primitive assumptions were fixed first and c* follows from them.
    """
    assert 1403.0 <= COMBUSTION.c_star_m_s <= 1587.0


def test_reference_combustion_properties_are_marked_illustrative():
    assert COMBUSTION.provenance is Provenance.ILLUSTRATIVE
    assert COMBUSTION.is_illustrative
    assert "ILLUSTRATIVE" in COMBUSTION.label
    assert "PRESCRIBED / ILLUSTRATIVE" in COMBUSTION.reference
    assert "NOT experimentally validated" in COMBUSTION.reference


def test_from_ideal_with_efficiency_reproduces_the_stated_efficiency():
    for eta in (0.90, 0.94, 0.96, 0.98, 1.0):
        props = CombustionProperties.from_ideal_with_efficiency(
            1.2, 2600.0, 0.022, eta, "check", "efficiency round-trip"
        )
        assert props.c_star_efficiency == pytest.approx(eta, rel=1e-13)


def test_replace_preserves_untouched_fields():
    hotter = COMBUSTION.replace(chamber_temperature_k=2860.0)
    assert hotter.gamma == COMBUSTION.gamma
    assert hotter.molar_mass_kg_mol == COMBUSTION.molar_mass_kg_mol
    assert hotter.c_star_m_s == COMBUSTION.c_star_m_s
    assert hotter.chamber_temperature_k == 2860.0
    # Raising T_c at fixed c* lowers the implied efficiency.
    assert hotter.c_star_efficiency < COMBUSTION.c_star_efficiency


def test_reference_chamber_pressure_is_inside_the_range_measured_by_the_source_study():
    """Rezaei et al. (2018) measured chamber pressures of 19.9-31.0 bar."""
    p_c_bar = solve_chamber(COMBUSTION, A_T, M_DOT_OX, M_DOT_F).chamber_pressure_pa / 1e5
    assert 19.9 <= p_c_bar <= 31.0
