"""Independent verification of the injector flow models.

Every expected value is rebuilt here from the equations exactly as printed by
Waxman, Zimmerman, Cantwell & Zilliac (NASA NTRS 20190001326), using raw ``math``
and direct CoolProp calls rather than this project's helpers.
"""

from __future__ import annotations

import itertools
import math

import CoolProp.CoolProp as CP
import pytest

from hybrid_rocket_motor.injector import (
    DYER_REFERENCE_DISCHARGE_COEFFICIENT,
    Injector,
    InjectorFlowStatus,
    InjectorModel,
)
from hybrid_rocket_motor.nitrous_properties import saturated_state
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
)

AREA = 3.5e-6
CD = DYER_REFERENCE_DISCHARGE_COEFFICIENT


@pytest.fixture(scope="module")
def tank_state():
    return REFERENCE_TANK.initial_state(
        REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION
    )


def reference_spi(cd_area: float, rho: float, drop: float) -> float:
    """Waxman Eq. (2): m_dot = C_d A sqrt(2 rho dP)."""
    return cd_area * math.sqrt(2.0 * rho * drop)


def reference_hem(cd_area: float, s1: float, h1: float, p2: float) -> float:
    """Waxman Eq. (3): m_dot = C_d A rho_2 sqrt(2 (h_1 - h_2)) at s_2 = s_1."""
    h2 = CP.PropsSI("H", "P", p2, "S", s1, "NitrousOxide")
    rho2 = CP.PropsSI("D", "P", p2, "S", s1, "NitrousOxide")
    return cd_area * rho2 * math.sqrt(2.0 * (h1 - h2))


def reference_kappa(p1: float, p2: float, pv: float) -> float:
    """Waxman Eq. (8): kappa = sqrt((P1 - P2)/(Pv - P2))."""
    return math.sqrt((p1 - p2) / (pv - p2))


# --- 19, 20, 21. SPI behaviour --------------------------------------------


def test_spi_matches_the_printed_equation():
    injector = Injector(AREA, CD, InjectorModel.SPI)
    for rho, drop in ((785.0, 2.4e6), (900.0, 1.0e6), (1000.0, 5.0e5)):
        assert injector.spi_mass_flow_kg_s(rho, drop) == pytest.approx(
            reference_spi(CD * AREA, rho, drop), rel=1e-14
        )


def test_spi_scales_as_the_square_root_of_pressure_drop():
    injector = Injector(AREA, CD, InjectorModel.SPI)
    base = injector.spi_mass_flow_kg_s(785.0, 1.0e6)
    assert injector.spi_mass_flow_kg_s(785.0, 4.0e6) == pytest.approx(2.0 * base, rel=1e-14)
    assert injector.spi_mass_flow_kg_s(785.0, 9.0e6) == pytest.approx(3.0 * base, rel=1e-14)


def test_spi_scales_linearly_with_effective_cd_area():
    small = Injector(AREA, CD, InjectorModel.SPI)
    large = Injector(2.0 * AREA, CD, InjectorModel.SPI)
    assert large.spi_mass_flow_kg_s(785.0, 2.4e6) == pytest.approx(
        2.0 * small.spi_mass_flow_kg_s(785.0, 2.4e6), rel=1e-14
    )
    half_cd = Injector(AREA, CD / 2.0, InjectorModel.SPI)
    assert half_cd.spi_mass_flow_kg_s(785.0, 2.4e6) == pytest.approx(
        0.5 * small.spi_mass_flow_kg_s(785.0, 2.4e6), rel=1e-14
    )


def test_spi_scales_as_the_square_root_of_density():
    injector = Injector(AREA, CD, InjectorModel.SPI)
    base = injector.spi_mass_flow_kg_s(250.0, 2.0e6)
    assert injector.spi_mass_flow_kg_s(1000.0, 2.0e6) == pytest.approx(2.0 * base, rel=1e-14)


def test_effective_cd_area_is_the_product():
    assert Injector(AREA, CD).effective_cd_area_m2 == pytest.approx(CD * AREA, rel=1e-15)


# --- 17, 18. no-flow conditions -------------------------------------------


@pytest.mark.parametrize("drop", [0.0, -1.0, -1.0e6])
def test_zero_or_negative_pressure_drop_gives_zero_spi_flow(drop):
    assert Injector(AREA, CD, InjectorModel.SPI).spi_mass_flow_kg_s(785.0, drop) == 0.0


def test_equal_pressures_give_a_no_flow_status(tank_state):
    result = Injector(AREA, CD).mass_flow(tank_state, tank_state.pressure_pa)
    assert result.status is InjectorFlowStatus.NO_FLOW_NO_PRESSURE_DROP
    assert result.mass_flow_kg_s == 0.0
    assert not result.is_flowing


def test_chamber_above_tank_pressure_gives_no_flow(tank_state):
    result = Injector(AREA, CD).mass_flow(tank_state, tank_state.pressure_pa * 1.2)
    assert result.status is InjectorFlowStatus.NO_FLOW_NO_PRESSURE_DROP
    assert result.mass_flow_kg_s == 0.0


# --- 22. HEM and Dyer -----------------------------------------------------


def test_hem_matches_the_printed_equation(tank_state):
    injector = Injector(AREA, CD, InjectorModel.HEM)
    for p2 in (1.0e6, 2.0e6, 2.6e6, 4.0e6):
        assert injector.hem_mass_flow_kg_s(
            tank_state.liquid_entropy_j_kg_k, tank_state.liquid_enthalpy_j_kg, p2
        ) == pytest.approx(
            reference_hem(
                CD * AREA,
                tank_state.liquid_entropy_j_kg_k,
                tank_state.liquid_enthalpy_j_kg,
                p2,
            ),
            rel=1e-12,
        )


def test_hem_is_finite_and_positive_across_the_operating_range(tank_state):
    injector = Injector(AREA, CD, InjectorModel.HEM)
    for p2 in (1.0e6, 2.0e6, 3.0e6, 4.5e6):
        flow = injector.mass_flow(tank_state, p2).mass_flow_kg_s
        assert math.isfinite(flow)
        assert flow > 0.0


def test_hem_predicts_less_flow_than_spi_for_saturated_nitrous(tank_state):
    """The physical point of the two-phase model: flashing chokes the orifice."""
    result = Injector(AREA, CD).mass_flow(tank_state, 2.6e6)
    assert result.hem_mass_flow_kg_s < result.spi_mass_flow_kg_s


def test_kappa_is_exactly_one_for_a_saturated_tank(tank_state):
    """With p1 = p_v the Dyer parameter collapses to 1 for any backpressure."""
    for p2 in (1.0e6, 2.6e6, 4.0e6):
        result = Injector(AREA, CD).mass_flow(tank_state, p2)
        assert result.non_equilibrium_parameter == pytest.approx(1.0, rel=1e-14)


def test_dyer_is_the_equal_blend_when_kappa_is_one(tank_state):
    result = Injector(AREA, CD).mass_flow(tank_state, 2.6e6)
    assert result.mass_flow_kg_s == pytest.approx(
        0.5 * (result.spi_mass_flow_kg_s + result.hem_mass_flow_kg_s), rel=1e-14
    )


def test_dyer_matches_the_printed_weighting_with_independent_terms(tank_state):
    p2 = 2.6e6
    result = Injector(AREA, CD).mass_flow(tank_state, p2)
    spi = reference_spi(CD * AREA, tank_state.liquid_density_kg_m3, tank_state.pressure_pa - p2)
    hem = reference_hem(
        CD * AREA, tank_state.liquid_entropy_j_kg_k, tank_state.liquid_enthalpy_j_kg, p2
    )
    kappa = reference_kappa(tank_state.pressure_pa, p2, tank_state.pressure_pa)
    expected = (kappa / (1.0 + kappa)) * spi + (1.0 / (1.0 + kappa)) * hem
    assert result.mass_flow_kg_s == pytest.approx(expected, rel=1e-12)


def test_dyer_lies_between_hem_and_spi(tank_state):
    result = Injector(AREA, CD).mass_flow(tank_state, 2.6e6)
    assert result.hem_mass_flow_kg_s < result.mass_flow_kg_s < result.spi_mass_flow_kg_s


def test_kappa_becomes_infinite_without_flashing_potential():
    """p2 >= p_v means no flashing, so the Dyer weighting must go to pure SPI."""
    assert Injector.non_equilibrium_parameter(5.0e6, 3.0e6, 3.0e6) == math.inf
    assert Injector.non_equilibrium_parameter(5.0e6, 3.5e6, 3.0e6) == math.inf


def test_kappa_is_zero_without_a_pressure_drop():
    assert Injector.non_equilibrium_parameter(3.0e6, 3.0e6, 5.0e6) == 0.0


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_kappa_rejects_non_finite_pressures(bad):
    with pytest.raises(ValueError):
        Injector.non_equilibrium_parameter(bad, 1.0e6, 5.0e6)


def test_spi_falls_monotonically_with_rising_backpressure(tank_state):
    flows = [
        Injector(AREA, CD, InjectorModel.SPI).mass_flow(tank_state, p2).mass_flow_kg_s
        for p2 in (1.0e6, 2.0e6, 3.0e6, 4.0e6)
    ]
    assert all(later < earlier for earlier, later in itertools.pairwise(flows))


def test_hem_has_a_maximum_with_respect_to_backpressure(tank_state):
    """Waxman Eq. (5): the HEM flow maximises at the critical backpressure.

    This is the well-known critical-flow behaviour, and it means the HEM branch
    is NOT monotone in backpressure -- over most of the operating range it
    *increases* as the chamber pressure rises.  Recorded here because it is
    physics, not a defect, and because the coupled solver's uniqueness argument
    has to account for it.
    """
    pressures = [0.9e5, 1.0e6, 2.0e6, 3.0e6, 3.6e6, 4.0e6, 4.5e6]
    flows = [
        Injector(AREA, CD, InjectorModel.HEM).mass_flow(tank_state, p2).mass_flow_kg_s
        for p2 in pressures
    ]
    peak = flows.index(max(flows))
    assert 0 < peak < len(flows) - 1, "the maximum must be interior"
    assert all(later > earlier for earlier, later in itertools.pairwise(flows[: peak + 1]))
    assert all(later < earlier for earlier, later in itertools.pairwise(flows[peak:]))


def test_dyer_inherits_a_shallow_maximum_from_the_hem_branch(tank_state):
    """The Dyer blend is also not strictly monotone, by roughly 3 % of flow."""
    injector = Injector(AREA, CD)
    low = injector.mass_flow(tank_state, 0.9e5).mass_flow_kg_s
    peak = injector.mass_flow(tank_state, 1.0e6).mass_flow_kg_s
    high = injector.mass_flow(tank_state, 4.0e6).mass_flow_kg_s
    assert peak > low
    assert peak > high
    assert (peak - low) / peak < 0.05


def test_reference_injector_flow_hardcoded_literals(tank_state):
    """3.5 mm^2 at C_d 0.66, tank 50.525 bar, chamber 26.11 bar."""
    result = Injector(AREA, CD).mass_flow(tank_state, 26.11e5)
    assert result.spi_mass_flow_kg_s == pytest.approx(0.14303, rel=1e-3)
    assert result.hem_mass_flow_kg_s == pytest.approx(0.06009, rel=1e-3)
    assert result.mass_flow_kg_s == pytest.approx(0.10156, rel=1e-3)


# --- 23, 24. input validation ---------------------------------------------


@pytest.mark.parametrize("bad_area", [0.0, -1.0e-6, math.nan, math.inf])
def test_invalid_area_is_rejected(bad_area):
    with pytest.raises(ValueError):
        Injector(bad_area, CD)


@pytest.mark.parametrize("bad_cd", [0.0, -0.5, 1.5, math.nan, math.inf])
def test_invalid_discharge_coefficient_is_rejected(bad_cd):
    with pytest.raises(ValueError):
        Injector(AREA, bad_cd)


def test_invalid_model_is_rejected():
    with pytest.raises(TypeError):
        Injector(AREA, CD, "DYER")


@pytest.mark.parametrize("bad_pressure", [-1.0, math.nan, math.inf])
def test_invalid_chamber_pressure_is_rejected(bad_pressure, tank_state):
    with pytest.raises(ValueError):
        Injector(AREA, CD).mass_flow(tank_state, bad_pressure)


def test_invalid_thermodynamic_state_is_rejected():
    with pytest.raises(TypeError):
        Injector(AREA, CD).mass_flow("not a tank state", 2.6e6)


@pytest.mark.parametrize("bad_density", [0.0, -100.0, math.nan])
def test_spi_rejects_invalid_density(bad_density):
    with pytest.raises(ValueError):
        Injector(AREA, CD).spi_mass_flow_kg_s(bad_density, 1.0e6)


def test_no_liquid_gives_a_dedicated_status():
    """A tank state with no liquid must not produce liquid flow."""
    state = REFERENCE_TANK.initial_state(280.0, 0.5)
    empty = type(state)(
        **{**state.__dict__, "liquid_mass_kg": 0.0, "phase": state.phase}
    )
    result = Injector(AREA, CD).mass_flow(empty, 1.0e6)
    assert result.status is InjectorFlowStatus.NO_FLOW_NO_LIQUID
    assert result.mass_flow_kg_s == 0.0


def test_injector_is_immutable():
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        Injector(AREA, CD).effective_area_m2 = 1.0


def test_colder_tank_gives_a_denser_liquid_and_lower_vapour_pressure():
    cold = saturated_state(273.15)
    warm = saturated_state(303.15)
    assert cold.liquid_density_kg_m3 > warm.liquid_density_kg_m3
    assert cold.pressure_pa < warm.pressure_pa
