"""Independent verification of the nested feed / variable-chemistry closure.

Independence strategy
---------------------
Every one of the six simultaneous relations the coupled state must satisfy is
re-formed here from the reported solution using the frozen Milestone 1-4
components directly, rather than being read back from the solver.  The closed-form
O/F inversion that defines the search window is checked against the forward
regression chain it inverts.  The frozen Milestone 4 coupling is re-run in the
same session to prove it is untouched.
"""

from __future__ import annotations

import math

import pytest

from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.feed_system import FeedStatus, solve_feed_coupling
from hybrid_rocket_motor.geometry import REPRESENTATIVE_GRAIN
from hybrid_rocket_motor.injector import Injector
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.regression import REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
)
from hybrid_rocket_motor.thermochemistry import default_table
from hybrid_rocket_motor.variable_feed_system import (
    VariableFeedStatus,
    oxidizer_flow_for_mixture_ratio,
    solve_variable_feed_coupling,
)

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
NOZZLE = REFERENCE_NOZZLE
INJECTOR = Injector(3.5e-6)
PORT_RADIUS = 0.020


@pytest.fixture(scope="module")
def table():
    return default_table()


@pytest.fixture(scope="module")
def tank_state():
    return REFERENCE_TANK.initial_state(REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION)


@pytest.fixture(scope="module")
def solution(table, tank_state):
    return solve_variable_feed_coupling(
        tank_state, INJECTOR, GRAIN, LAW, table, NOZZLE, PORT_RADIUS
    )


# -- the closed-form O/F inversion -------------------------------------------------


@pytest.mark.parametrize("of_ratio", [1.2, 1.75, 2.5, 3.3, 4.0])
@pytest.mark.parametrize("radius", [0.020, 0.031, 0.0445])
def test_oxidizer_flow_inversion_round_trips_through_the_regression_chain(of_ratio, radius):
    """Invert O/F, then run the forward chain and recover the same O/F."""
    oxidizer = oxidizer_flow_for_mixture_ratio(GRAIN, LAW, radius, of_ratio)
    flux = oxidizer / (math.pi * radius * radius)
    rate = LAW.regression_rate_si(flux)
    fuel = GRAIN.fuel_density_kg_m3 * 2.0 * math.pi * radius * GRAIN.length_m * rate
    assert oxidizer / fuel == pytest.approx(of_ratio, rel=1e-12)


def test_oxidizer_flow_inversion_is_monotone_in_mixture_ratio():
    previous = 0.0
    for of_ratio in (1.2, 2.0, 3.0, 4.0):
        flow = oxidizer_flow_for_mixture_ratio(GRAIN, LAW, PORT_RADIUS, of_ratio)
        assert flow > previous
        previous = flow


@pytest.mark.parametrize("of_ratio", [0.0, -1.0, math.nan, math.inf])
def test_invalid_mixture_ratio_is_refused_by_the_inversion(of_ratio):
    with pytest.raises(ValueError):
        oxidizer_flow_for_mixture_ratio(GRAIN, LAW, PORT_RADIUS, of_ratio)


# -- all six relations hold simultaneously -------------------------------------------


def test_reference_point_is_flowing(solution):
    assert solution.status is VariableFeedStatus.FLOWING
    assert solution.is_flowing
    assert solution.oxidizer_mass_flow_kg_s > 0.0


def test_relation_1_injector_flow(solution, tank_state):
    """The oxidizer flow must be the flow the injector passes at this pressure."""
    injector_flow = INJECTOR.mass_flow(tank_state, solution.chamber_pressure_pa).mass_flow_kg_s
    assert injector_flow == pytest.approx(solution.oxidizer_mass_flow_kg_s, abs=1e-12)
    assert abs(solution.feed_residual_kg_s) < 1e-12


def test_relation_2_and_3_regression_law_and_fuel_geometry(solution):
    """Re-derive G_ox, r_dot and m_dot_f from the frozen Milestone 1 chain."""
    area = math.pi * PORT_RADIUS * PORT_RADIUS
    flux = solution.oxidizer_mass_flow_kg_s / area
    assert flux == pytest.approx(solution.oxidizer_mass_flux_si, rel=1e-14)
    rate = LAW.regression_rate_si(flux)
    assert rate == pytest.approx(solution.regression_rate_m_s, rel=1e-14)
    fuel = GRAIN.fuel_density_kg_m3 * 2.0 * math.pi * PORT_RADIUS * GRAIN.length_m * rate
    assert fuel == pytest.approx(solution.fuel_mass_flow_kg_s, rel=1e-14)


def test_relation_4_mixture_ratio_definition(solution):
    assert solution.mixture_ratio == pytest.approx(
        solution.oxidizer_mass_flow_kg_s / solution.fuel_mass_flow_kg_s, rel=1e-15
    )


def test_relation_5_and_6_chamber_closure_with_tabulated_cstar(solution, table):
    """``p_c = m_dot c*(O/F, p_c)/A_t`` with c* read back from the table."""
    state = table.evaluate(solution.mixture_ratio, solution.chamber_pressure_pa)
    implied = solution.total_mass_flow_kg_s * state.c_star_m_s / NOZZLE.throat_area_m2
    assert implied == pytest.approx(solution.chamber_pressure_pa, rel=1e-11)
    assert abs(solution.chamber_residual_pa) < 1e-6


def test_reported_state_is_inside_the_table(solution, table):
    assert solution.thermochemistry.is_within_table
    of_low, of_high = table.of_bounds
    p_low, p_high = table.pressure_bounds_pa
    assert of_low <= solution.mixture_ratio <= of_high
    assert p_low <= solution.chamber_pressure_pa <= p_high


def test_pressure_margin_is_consistent_and_positive(solution, tank_state):
    assert solution.tank_pressure_pa == pytest.approx(tank_state.pressure_pa, rel=1e-15)
    assert solution.pressure_margin_pa == pytest.approx(
        solution.tank_pressure_pa - solution.chamber_pressure_pa, rel=1e-12
    )
    assert solution.pressure_margin_pa > 0.0


# -- comparison with the frozen Milestone 4 coupling ------------------------------------


def test_frozen_milestone_4_coupling_is_untouched(tank_state):
    """Milestone 4 must still converge to its own answer, unchanged."""
    frozen = solve_feed_coupling(
        tank_state,
        INJECTOR,
        GRAIN,
        LAW,
        ILLUSTRATIVE_N2O_HTPB_COMBUSTION,
        NOZZLE,
        PORT_RADIUS,
    )
    assert frozen.status is FeedStatus.FLOWING
    injector_flow = INJECTOR.mass_flow(tank_state, frozen.chamber_pressure_pa).mass_flow_kg_s
    assert injector_flow == pytest.approx(frozen.oxidizer_mass_flow_kg_s, abs=1e-12)


def test_variable_model_lowers_pressure_and_raises_flow(solution, tank_state):
    """The lower equilibrium c* drops p_c, which lets the injector pass more.

    This coupled response is the whole reason Milestone 5 is not simply Milestone 4
    with different constants: the chemistry and the feed move each other.
    """
    frozen = solve_feed_coupling(
        tank_state,
        INJECTOR,
        GRAIN,
        LAW,
        ILLUSTRATIVE_N2O_HTPB_COMBUSTION,
        NOZZLE,
        PORT_RADIUS,
    )
    assert solution.chamber_pressure_pa < frozen.chamber_pressure_pa
    assert solution.oxidizer_mass_flow_kg_s > frozen.oxidizer_mass_flow_kg_s
    assert solution.mixture_ratio > frozen.mixture_ratio


# -- behaviour across the grain and the operating envelope --------------------------------


@pytest.mark.parametrize("radius", [0.0201, 0.025, 0.030, 0.035, 0.040, 0.0448])
def test_solution_closes_across_the_whole_grain(table, tank_state, radius):
    solution = solve_variable_feed_coupling(tank_state, INJECTOR, GRAIN, LAW, table, NOZZLE, radius)
    assert solution.is_flowing
    assert abs(solution.feed_residual_kg_s) < 1e-11
    assert abs(solution.chamber_residual_pa) < 1e-6
    assert solution.thermochemistry.is_within_table


def test_mixture_ratio_falls_as_the_port_opens(table, tank_state):
    """The Milestone 2 result must survive the new chemistry: O/F falls with burn."""
    previous = math.inf
    for radius in (0.0201, 0.025, 0.030, 0.035, 0.040, 0.0448):
        solution = solve_variable_feed_coupling(
            tank_state, INJECTOR, GRAIN, LAW, table, NOZZLE, radius
        )
        assert solution.mixture_ratio < previous
        previous = solution.mixture_ratio


@pytest.mark.parametrize("area", [2.5e-6, 3.5e-6, 4.5e-6])
def test_solution_closes_across_the_injector_range(table, tank_state, area):
    solution = solve_variable_feed_coupling(
        tank_state, Injector(area), GRAIN, LAW, table, NOZZLE, PORT_RADIUS
    )
    assert solution.is_flowing
    assert abs(solution.feed_residual_kg_s) < 1e-11


@pytest.mark.parametrize("temperature", [283.15, 293.15, 303.15])
def test_solution_closes_across_the_tank_temperature_range(table, temperature):
    state = REFERENCE_TANK.initial_state(temperature, REFERENCE_TANK_FILL_FRACTION)
    solution = solve_variable_feed_coupling(state, INJECTOR, GRAIN, LAW, table, NOZZLE, PORT_RADIUS)
    assert solution.is_flowing
    assert abs(solution.feed_residual_kg_s) < 1e-11


# -- explicit failure instead of clamping ------------------------------------------------


def test_depleted_liquid_is_reported_not_solved(table, tank_state):
    """Once the liquid is gone the model has nothing to say, and says so."""
    tank = REFERENCE_TANK
    # 5 % of the mass at 10 % of the energy is a fully vaporised tank: the flash
    # returns quality 1 with no liquid left, which is exactly the terminal state
    # the blowdown integrator stops on.
    empty = tank.state_from_mass_and_energy(
        tank_state.total_mass_kg * 0.05,
        tank_state.total_internal_energy_j * 0.10,
        quality_overshoot=0.6,
    )
    assert not empty.has_liquid
    solution = solve_variable_feed_coupling(empty, INJECTOR, GRAIN, LAW, table, NOZZLE, PORT_RADIUS)
    assert solution.status is VariableFeedStatus.LIQUID_DEPLETED
    assert not solution.is_flowing
    assert solution.oxidizer_mass_flow_kg_s == 0.0
    assert solution.thermochemistry is None


def test_a_tiny_injector_falls_off_the_bottom_of_the_table(table, tank_state):
    """Too little flow cannot reach the tabulated pressure floor; refuse it."""
    solution = solve_variable_feed_coupling(
        tank_state, Injector(1.0e-8), GRAIN, LAW, table, NOZZLE, PORT_RADIUS
    )
    assert not solution.is_flowing
    assert solution.status in {
        VariableFeedStatus.OF_OUTSIDE_TABLE,
        VariableFeedStatus.PRESSURE_OUTSIDE_TABLE,
    }
    assert solution.thermochemistry is None


def test_a_failed_solve_reports_no_thermochemical_state(table, tank_state):
    solution = solve_variable_feed_coupling(
        tank_state, Injector(1.0e-8), GRAIN, LAW, table, NOZZLE, PORT_RADIUS
    )
    assert solution.thermochemistry is None
    assert solution.injector is None
    assert solution.feed_residual_kg_s == 0.0


@pytest.mark.parametrize("radius", [0.0, -0.01, math.nan, math.inf])
def test_invalid_port_radius_is_refused(table, tank_state, radius):
    with pytest.raises(ValueError):
        solve_variable_feed_coupling(tank_state, INJECTOR, GRAIN, LAW, table, NOZZLE, radius)


def test_a_burnt_through_grain_is_refused(table, tank_state):
    outer_radius = GRAIN.outer_diameter_m / 2.0
    with pytest.raises(ValueError):
        solve_variable_feed_coupling(
            tank_state, INJECTOR, GRAIN, LAW, table, NOZZLE, 1.05 * outer_radius
        )


def test_burnout_overshoot_allowance_is_numerical_only(table, tank_state):
    """The allowance lets the integrator bracket burnout; it must not change a result."""
    outer_radius = GRAIN.outer_diameter_m / 2.0
    inside = solve_variable_feed_coupling(
        tank_state, INJECTOR, GRAIN, LAW, table, NOZZLE, 0.999 * outer_radius
    )
    with_allowance = solve_variable_feed_coupling(
        tank_state,
        INJECTOR,
        GRAIN,
        LAW,
        table,
        NOZZLE,
        0.999 * outer_radius,
        burnout_overshoot_fraction=0.05,
    )
    assert with_allowance.chamber_pressure_pa == pytest.approx(
        inside.chamber_pressure_pa, rel=1e-15
    )
