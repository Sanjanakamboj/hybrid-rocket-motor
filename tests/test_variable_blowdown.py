"""Independent verification of the coupled variable-property blowdown.

Independence strategy
---------------------
Closure residuals are re-formed here from the reported histories -- fuel mass from
the swept port geometry, oxidizer mass from the tank inventory, impulse from a
trapezoidal integral of the reported thrust -- rather than read from the solver's
own accumulators.  Individual samples are re-solved from scratch with the static
coupled solver.  The frozen Milestone 4 blowdown is re-run in the same session on
identical conditions to prove it is untouched and to isolate the chemistry.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from hybrid_rocket_motor.blowdown import BlowdownTermination, simulate_blowdown
from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.geometry import REPRESENTATIVE_GRAIN
from hybrid_rocket_motor.injector import Injector
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.performance import STANDARD_GRAVITY_M_S2
from hybrid_rocket_motor.regression import REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
    NitrousTank,
)
from hybrid_rocket_motor.thermochemistry import default_table
from hybrid_rocket_motor.variable_blowdown import (
    VariableBlowdownTermination,
    simulate_variable_blowdown,
)
from hybrid_rocket_motor.variable_feed_system import (
    VariableFeedStatus,
    solve_variable_feed_coupling,
)
from hybrid_rocket_motor.variable_nozzle import solve_variable_performance_point

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
NOZZLE = REFERENCE_NOZZLE
INJECTOR = Injector(3.5e-6)
SEA_LEVEL_PA = 101325.0


@pytest.fixture(scope="module")
def table():
    return default_table()


@pytest.fixture(scope="module")
def initial_state():
    return REFERENCE_TANK.initial_state(REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION)


@pytest.fixture(scope="module")
def nominal(table, initial_state):
    return simulate_variable_blowdown(
        REFERENCE_TANK, initial_state, INJECTOR, GRAIN, LAW, table, NOZZLE, SEA_LEVEL_PA
    )


@pytest.fixture(scope="module")
def frozen_nominal(initial_state):
    return simulate_blowdown(
        REFERENCE_TANK,
        initial_state,
        INJECTOR,
        GRAIN,
        LAW,
        ILLUSTRATIVE_N2O_HTPB_COMBUSTION,
        NOZZLE,
        SEA_LEVEL_PA,
    )


# -- the run terminates for a stated physical reason ---------------------------------


def test_nominal_case_burns_the_grain_through(nominal):
    assert nominal.termination is VariableBlowdownTermination.GRAIN_BURNOUT
    assert nominal.burn_time_s > 0.0


def test_nominal_case_never_leaves_the_table(nominal):
    assert nominal.stayed_within_table
    assert np.all(nominal.of_margin_low[nominal.firing] > 0.0)
    assert np.all(nominal.of_margin_high[nominal.firing] > 0.0)
    assert np.all(nominal.pressure_margin_low_pa[nominal.firing] > 0.0)
    assert np.all(nominal.pressure_margin_high_pa[nominal.firing] > 0.0)


def test_port_radius_is_clamped_at_the_outer_radius_on_burnout(nominal):
    outer = GRAIN.outer_diameter_m / 2.0
    assert np.all(nominal.port_radius_m <= outer + 1e-15)
    assert nominal.port_radius_m[-1] == pytest.approx(outer, rel=1e-9)


def test_oxidizer_limited_case_terminates_on_liquid_depletion(table):
    small_tank = NitrousTank(2.5e-3)
    state = small_tank.initial_state(REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION)
    result = simulate_variable_blowdown(
        small_tank, state, INJECTOR, GRAIN, LAW, table, NOZZLE, SEA_LEVEL_PA
    )
    assert result.termination is VariableBlowdownTermination.LIQUID_DEPLETED
    assert result.stayed_within_table


# -- closures re-formed independently of the solver's accumulators ---------------------


def test_fuel_mass_matches_the_swept_port_geometry(nominal):
    swept = (
        GRAIN.fuel_density_kg_m3
        * GRAIN.length_m
        * math.pi
        * (nominal.port_radius_m[-1] ** 2 - nominal.port_radius_m[0] ** 2)
    )
    assert nominal.fuel_consumed_kg == pytest.approx(swept, rel=1e-6)


def test_oxidizer_mass_matches_the_tank_inventory(nominal):
    drawn = nominal.tank_total_mass_kg[0] - nominal.tank_total_mass_kg[-1]
    assert nominal.oxidizer_consumed_kg == pytest.approx(drawn, abs=1e-12)
    assert np.max(np.abs(nominal.oxidizer_mass_closure_residual_kg)) < 1e-12


def test_impulse_matches_a_trapezoidal_integral_of_thrust(nominal):
    trapezoid = float(np.trapezoid(nominal.thrust_n, nominal.time_s))
    assert nominal.total_impulse_n_s == pytest.approx(trapezoid, rel=1e-5)


def test_equivalent_isp_matches_its_definition(nominal):
    expected = nominal.total_impulse_n_s / (nominal.propellant_consumed_kg * STANDARD_GRAVITY_M_S2)
    assert nominal.equivalent_specific_impulse_s == pytest.approx(expected, rel=1e-12)


def test_chamber_closure_holds_at_every_reported_sample(nominal):
    assert np.max(np.abs(nominal.c_star_identity_residual_pa)) < 1e-6


def test_mixture_ratio_identity_holds_at_every_reported_sample(nominal):
    assert np.max(np.abs(nominal.mixture_ratio_identity_residual)) < 1e-12


def test_feed_residual_is_at_machine_precision_at_every_sample(nominal):
    assert np.max(np.abs(nominal.feed_residual_kg_s)) < 1e-11


def test_gas_constant_matches_molar_mass_at_every_sample(nominal):
    firing = nominal.firing
    expected = 8.314462618 / nominal.effective_molar_mass_kg_mol[firing]
    assert np.allclose(nominal.specific_gas_constant_j_kg_k[firing], expected, rtol=1e-13)


# -- individual samples re-solved from scratch -------------------------------------------


@pytest.mark.parametrize("index", [0, 80, 160, 200])
def test_sampled_states_re_solve_to_the_same_operating_point(nominal, table, index):
    """Re-run the static coupled solver at a reported state and compare."""
    tank_state = REFERENCE_TANK.state_from_mass_and_energy(
        float(nominal.tank_total_mass_kg[index]),
        float(nominal.tank_internal_energy_j[index]),
        quality_overshoot=0.6,
    )
    solution = solve_variable_feed_coupling(
        tank_state,
        INJECTOR,
        GRAIN,
        LAW,
        table,
        NOZZLE,
        float(nominal.port_radius_m[index]),
        burnout_overshoot_fraction=0.05,
    )
    assert solution.status is VariableFeedStatus.FLOWING
    assert solution.oxidizer_mass_flow_kg_s == pytest.approx(
        float(nominal.oxidizer_mass_flow_kg_s[index]), rel=1e-9
    )
    assert solution.chamber_pressure_pa == pytest.approx(
        float(nominal.chamber_pressure_pa[index]), rel=1e-9
    )
    assert solution.mixture_ratio == pytest.approx(float(nominal.mixture_ratio[index]), rel=1e-9)


@pytest.mark.parametrize("index", [0, 120, 200])
def test_reported_thrust_re_derives_from_the_reported_state(nominal, table, index):
    """The frozen Milestone 3 nozzle, driven by the reported chemistry, must agree."""
    state = table.evaluate(
        float(nominal.mixture_ratio[index]), float(nominal.chamber_pressure_pa[index])
    )
    point = solve_variable_performance_point(
        state,
        NOZZLE,
        float(nominal.oxidizer_mass_flow_kg_s[index]),
        float(nominal.fuel_mass_flow_kg_s[index]),
        SEA_LEVEL_PA,
    )
    assert point.thrust_n == pytest.approx(float(nominal.thrust_n[index]), rel=1e-9)
    assert point.exit_state.mach == pytest.approx(float(nominal.exit_mach[index]), rel=1e-9)


def test_thrust_is_the_sum_of_its_two_reported_terms(nominal):
    firing = nominal.firing
    total = nominal.momentum_thrust_n[firing] + nominal.pressure_thrust_n[firing]
    assert np.allclose(total, nominal.thrust_n[firing], rtol=1e-12)


# -- physical behaviour of the histories ---------------------------------------------------


def test_mixture_ratio_and_chamber_pressure_both_fall_over_the_burn(nominal):
    firing = nominal.firing
    assert nominal.mixture_ratio[firing][-1] < nominal.mixture_ratio[firing][0]
    assert nominal.chamber_pressure_pa[firing][-1] < nominal.chamber_pressure_pa[firing][0]


def test_characteristic_velocity_falls_with_the_mixture_ratio(nominal):
    """On the fuel-rich side of the table a falling O/F means a falling c*.

    This is the Milestone 5 mechanism that Milestone 4 could not represent, so it
    is asserted rather than merely described.
    """
    firing = nominal.firing
    of_values = nominal.mixture_ratio[firing]
    c_star = nominal.c_star_m_s[firing]
    assert of_values[-1] < of_values[0]
    assert c_star[-1] < c_star[0]
    assert np.corrcoef(of_values, c_star)[0, 1] > 0.99


def test_tank_empties_monotonically(nominal):
    assert np.all(np.diff(nominal.tank_total_mass_kg) <= 1e-12)
    assert np.all(np.diff(nominal.tank_temperature_k) <= 1e-9)
    assert np.all(np.diff(nominal.cumulative_impulse_n_s) >= -1e-9)


def test_port_opens_monotonically(nominal):
    assert np.all(np.diff(nominal.port_radius_m) >= -1e-15)


# -- comparison against the untouched Milestone 4 model --------------------------------------


def test_frozen_milestone_4_model_is_unchanged(frozen_nominal):
    assert frozen_nominal.termination is BlowdownTermination.GRAIN_BURNOUT
    assert frozen_nominal.total_impulse_n_s > 0.0
    assert frozen_nominal.combustion is ILLUSTRATIVE_N2O_HTPB_COMBUSTION


def test_both_models_consume_the_same_fuel_because_the_grain_sets_it(nominal, frozen_nominal):
    """Burnout is geometric, so the fuel mass cannot depend on the chemistry."""
    assert nominal.fuel_consumed_kg == pytest.approx(frozen_nominal.fuel_consumed_kg, rel=1e-5)


def test_variable_chemistry_amplifies_the_thrust_decay_without_reversing_it(
    nominal, frozen_nominal
):
    """The Milestone 5 headline result, asserted so a regression would fail loudly."""
    frozen_decay = frozen_nominal.thrust_n[-1] / frozen_nominal.thrust_n[0] - 1.0
    variable_decay = nominal.thrust_n[-1] / nominal.thrust_n[0] - 1.0
    assert frozen_decay < 0.0, "the Milestone 4 result was a falling thrust"
    assert variable_decay < 0.0, "the sign of that result must survive"
    assert variable_decay < frozen_decay, "the decay should be amplified, not damped"
    assert 1.2 < variable_decay / frozen_decay < 2.0


def test_variable_chemistry_lowers_the_whole_thrust_curve(nominal, frozen_nominal):
    assert nominal.total_impulse_n_s < frozen_nominal.total_impulse_n_s
    assert nominal.peak_thrust_n < frozen_nominal.peak_thrust_n
    relative = nominal.total_impulse_n_s / frozen_nominal.total_impulse_n_s - 1.0
    assert -0.15 < relative < -0.05


# -- efficiency and table variants -----------------------------------------------------------


def test_efficiency_moves_chamber_pressure_but_barely_moves_impulse(table, initial_state):
    """The pressure-fed injector self-compensates; the study explains why."""
    runs = {}
    for efficiency in (0.92, 1.00):
        runs[efficiency] = simulate_variable_blowdown(
            REFERENCE_TANK,
            initial_state,
            INJECTOR,
            GRAIN,
            LAW,
            table,
            NOZZLE,
            SEA_LEVEL_PA,
            c_star_efficiency=efficiency,
            n_report=51,
        )
    low, high = runs[0.92], runs[1.00]
    pressure_change = high.chamber_pressure_pa[0] / low.chamber_pressure_pa[0] - 1.0
    impulse_change = high.total_impulse_n_s / low.total_impulse_n_s - 1.0
    assert pressure_change > 0.05, "a higher c* must raise chamber pressure"
    assert abs(impulse_change) < 0.005, "yet total impulse is nearly unchanged"
    assert high.oxidizer_mass_flow_kg_s[0] < low.oxidizer_mass_flow_kg_s[0]


def test_a_coarser_table_changes_the_answer_only_slightly(table, initial_state, nominal):
    coarse = simulate_variable_blowdown(
        REFERENCE_TANK,
        initial_state,
        INJECTOR,
        GRAIN,
        LAW,
        table.subsample(of_stride=8, pressure_stride=4),
        NOZZLE,
        SEA_LEVEL_PA,
        n_report=51,
    )
    relative = coarse.total_impulse_n_s / nominal.total_impulse_n_s - 1.0
    assert abs(relative) < 1e-3


# -- argument validation -------------------------------------------------------------------------


@pytest.mark.parametrize("n_report", [0, 1, -5])
def test_invalid_report_count_is_refused(table, initial_state, n_report):
    with pytest.raises(ValueError):
        simulate_variable_blowdown(
            REFERENCE_TANK,
            initial_state,
            INJECTOR,
            GRAIN,
            LAW,
            table,
            NOZZLE,
            SEA_LEVEL_PA,
            n_report=n_report,
        )


@pytest.mark.parametrize("ambient", [-1.0, math.nan, math.inf])
def test_invalid_ambient_pressure_is_refused(table, initial_state, ambient):
    with pytest.raises(ValueError):
        simulate_variable_blowdown(
            REFERENCE_TANK, initial_state, INJECTOR, GRAIN, LAW, table, NOZZLE, ambient
        )


def test_result_carries_its_inputs_for_provenance(nominal, table):
    assert nominal.table is table
    assert nominal.geometry is GRAIN
    assert nominal.law is LAW
    assert nominal.nozzle is NOZZLE
    assert nominal.ambient_pressure_pa == SEA_LEVEL_PA
