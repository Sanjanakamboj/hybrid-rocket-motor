"""Independent verification of the coupled N2O blowdown model.

Independence strategy
---------------------
Mass, geometry, impulse and energy closures are re-formed here from the reported
histories rather than read back from the solver's own accumulators, and the
initial coupled state is re-solved from scratch.  The frozen Milestone 3
prescribed-flow reference is re-run in the same session to prove it is untouched.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from hybrid_rocket_motor.blowdown import (
    BlowdownTermination,
    simulate_blowdown,
)
from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.feed_system import FeedStatus, solve_feed_coupling
from hybrid_rocket_motor.geometry import REPRESENTATIVE_GRAIN
from hybrid_rocket_motor.injector import Injector
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.performance import (
    STANDARD_GRAVITY_M_S2,
    couple_transient_to_performance,
    solve_performance_point,
)
from hybrid_rocket_motor.regression import REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
)
from hybrid_rocket_motor.transient import ConstantOxidizerFlow, simulate_transient

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
COMBUSTION = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
NOZZLE = REFERENCE_NOZZLE
TANK = REFERENCE_TANK
AREA = 3.5e-6
SEA_LEVEL_PA = 101325.0


@pytest.fixture(scope="module")
def initial_state():
    return TANK.initial_state(REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION)


@pytest.fixture(scope="module")
def nominal(initial_state):
    return simulate_blowdown(
        TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=201,
    )


# --- 33. t = 0 matches an independently solved static state ---------------


def test_initial_state_matches_a_standalone_feed_solve(nominal, initial_state):
    standalone = solve_feed_coupling(
        initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, 0.020
    )
    assert nominal.time_s[0] == 0.0
    assert nominal.oxidizer_mass_flow_kg_s[0] == pytest.approx(
        standalone.oxidizer_mass_flow_kg_s, rel=1e-10
    )
    assert nominal.fuel_mass_flow_kg_s[0] == pytest.approx(
        standalone.fuel_mass_flow_kg_s, rel=1e-10
    )
    assert nominal.chamber_pressure_pa[0] == pytest.approx(
        standalone.chamber_pressure_pa, rel=1e-10
    )
    assert nominal.tank_pressure_pa[0] == pytest.approx(initial_state.pressure_pa, rel=1e-9)


def test_initial_thrust_matches_a_standalone_performance_solve(nominal):
    point = solve_performance_point(
        COMBUSTION,
        NOZZLE,
        float(nominal.oxidizer_mass_flow_kg_s[0]),
        float(nominal.fuel_mass_flow_kg_s[0]),
        SEA_LEVEL_PA,
    )
    assert nominal.thrust_n[0] == pytest.approx(point.thrust_n, rel=1e-12)
    assert nominal.exit_mach[0] == pytest.approx(point.exit_state.mach, rel=1e-12)


def test_reference_case_hardcoded_literals(nominal):
    """10 L tank, 80 % fill, 293.15 K, 3.5 mm^2 Dyer injector."""
    assert nominal.tank_pressure_pa[0] / 1e5 == pytest.approx(50.525, rel=1e-4)
    assert nominal.oxidizer_mass_flow_kg_s[0] == pytest.approx(0.10132, rel=1e-3)
    assert nominal.chamber_pressure_pa[0] / 1e5 == pytest.approx(26.396, rel=1e-3)
    assert nominal.thrust_n[0] == pytest.approx(313.73, rel=1e-3)
    assert nominal.termination is BlowdownTermination.GRAIN_BURNOUT
    assert nominal.burn_time_s == pytest.approx(42.917, rel=1e-3)


# --- 34, 35. monotonic behaviour ------------------------------------------


def test_port_radius_is_nondecreasing(nominal):
    assert np.all(np.diff(nominal.port_radius_m) >= 0.0)
    assert np.all(np.diff(nominal.port_diameter_m) >= 0.0)


def test_tank_mass_is_nonincreasing(nominal):
    assert np.all(np.diff(nominal.tank_total_mass_kg) <= 0.0)
    assert np.all(np.diff(nominal.tank_liquid_mass_kg) <= 0.0)


def test_tank_pressure_and_temperature_fall_monotonically(nominal):
    assert np.all(np.diff(nominal.tank_pressure_pa) < 0.0)
    assert np.all(np.diff(nominal.tank_temperature_k) < 0.0)


def test_cumulative_quantities_are_nondecreasing(nominal):
    for array in (
        nominal.cumulative_oxidizer_mass_kg,
        nominal.cumulative_fuel_mass_kg,
        nominal.cumulative_impulse_n_s,
    ):
        assert np.all(np.diff(array) >= 0.0)


def test_blowdown_cools_the_tank_and_lowers_the_feed_pressure(nominal):
    assert nominal.tank_temperature_k[-1] < nominal.tank_temperature_k[0]
    assert nominal.tank_pressure_pa[-1] < nominal.tank_pressure_pa[0]
    assert nominal.oxidizer_mass_flow_kg_s[-1] < nominal.oxidizer_mass_flow_kg_s[0]


# --- 36, 37, 38. closure checks -------------------------------------------


def test_oxidizer_mass_closure(nominal):
    """Oxidizer integrated at the injector must equal that drawn from the tank."""
    drawn = nominal.tank_total_mass_kg[0] - nominal.tank_total_mass_kg
    assert np.allclose(nominal.cumulative_oxidizer_mass_kg, drawn, rtol=0.0, atol=1e-9)
    assert np.max(np.abs(nominal.oxidizer_mass_closure_residual_kg)) < 1e-9


def test_fuel_mass_closure_against_the_port_geometry(nominal):
    """Integrated fuel flow vs rho_f L pi (r^2 - r_0^2), formed here."""
    geometric = 930.0 * 0.40 * math.pi * (
        nominal.port_radius_m**2 - nominal.port_radius_m[0] ** 2
    )
    assert np.allclose(nominal.cumulative_fuel_mass_kg, geometric, rtol=0.0, atol=1e-6)
    assert np.max(np.abs(nominal.fuel_mass_closure_residual_kg)) < 1e-6


def test_burnout_consumes_the_whole_grain(nominal):
    assert nominal.fuel_consumed_kg == pytest.approx(GRAIN.initial_fuel_mass_kg, rel=1e-5)
    assert nominal.port_diameter_m[-1] == pytest.approx(GRAIN.outer_diameter_m, rel=1e-9)


def test_impulse_closure(nominal):
    """Integrated impulse state vs an independent trapezoidal integral of thrust."""
    manual = float(
        np.sum(
            0.5
            * (nominal.thrust_n[1:] + nominal.thrust_n[:-1])
            * np.diff(nominal.time_s)
        )
    )
    assert nominal.total_impulse_n_s == pytest.approx(manual, rel=1e-5)
    assert abs(nominal.impulse_closure_residual_n_s) < 1e-2


def test_energy_closure_of_the_tank(nominal):
    """The adiabatic energy balance residual, normalised by the initial energy."""
    residual = np.max(np.abs(nominal.energy_closure_residual_j))
    assert residual / abs(nominal.tank_internal_energy_j[0]) < 1e-6


def test_equivalent_specific_impulse_identity(nominal):
    expected = nominal.total_impulse_n_s / (
        nominal.propellant_consumed_kg * STANDARD_GRAVITY_M_S2
    )
    assert nominal.equivalent_specific_impulse_s == pytest.approx(expected, rel=1e-15)
    assert 100.0 < nominal.equivalent_specific_impulse_s < 400.0


def test_feed_residual_is_negligible_throughout(nominal):
    assert np.max(np.abs(nominal.feed_residual_kg_s)) < 1e-9


# --- 39, 40, 41. finite, positive firing states ---------------------------


def test_all_firing_states_are_finite_and_positive(nominal):
    for array in (
        nominal.oxidizer_mass_flow_kg_s,
        nominal.fuel_mass_flow_kg_s,
        nominal.chamber_pressure_pa,
        nominal.thrust_n,
        nominal.tank_pressure_pa,
        nominal.tank_temperature_k,
        nominal.exit_velocity_m_s,
        nominal.specific_impulse_s,
    ):
        assert np.all(np.isfinite(array))
        assert np.all(array > 0.0)


def test_pressure_margin_stays_positive(nominal):
    assert np.all(nominal.pressure_margin_pa > 0.0)
    assert nominal.minimum_pressure_margin_pa > 0.0


def test_exit_mach_is_supersonic_and_constant(nominal):
    assert np.all(nominal.exit_mach > 1.0)
    assert np.max(nominal.exit_mach) - np.min(nominal.exit_mach) < 1e-9


def test_every_reported_sample_is_flowing(nominal):
    assert all(status is FeedStatus.FLOWING for status in nominal.feed_status)


def test_thrust_decomposition_is_exact(nominal):
    assert np.allclose(
        nominal.thrust_n,
        nominal.momentum_thrust_n + nominal.pressure_thrust_n,
        rtol=0.0,
        atol=1e-12,
    )


# --- 42, 43, 44. terminal events ------------------------------------------


def test_nominal_case_terminates_on_grain_burnout(nominal):
    assert nominal.termination is BlowdownTermination.GRAIN_BURNOUT
    assert nominal.port_radius_m[-1] == pytest.approx(GRAIN.outer_diameter_m / 2.0, rel=1e-9)
    assert nominal.tank_liquid_mass_kg[-1] > 0.0


def test_no_state_is_integrated_past_burnout(nominal):
    assert np.all(nominal.port_radius_m <= GRAIN.outer_diameter_m / 2.0 + 1e-12)


def test_a_small_tank_terminates_on_liquid_depletion(initial_state):
    """Too little oxidizer to burn the grain through: liquid runs out first."""
    from hybrid_rocket_motor.tank import NitrousTank

    small = NitrousTank(volume_m3=0.0025)
    state = small.initial_state(REFERENCE_TANK_TEMPERATURE_K, 0.80)
    result = simulate_blowdown(
        small, state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=101,
    )
    assert result.termination is BlowdownTermination.LIQUID_DEPLETED
    assert result.tank_liquid_mass_kg[-1] < 1e-3
    assert result.port_diameter_m[-1] < GRAIN.outer_diameter_m
    assert np.all(result.tank_liquid_mass_kg >= 0.0)


def test_a_short_horizon_completes_the_time_window(initial_state):
    result = simulate_blowdown(
        TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, t_end_s=5.0, n_report=51,
    )
    assert result.termination is BlowdownTermination.COMPLETED_TIME_WINDOW
    assert result.burn_time_s == pytest.approx(5.0, rel=1e-9)
    assert result.port_diameter_m[-1] < GRAIN.outer_diameter_m
    assert result.tank_liquid_mass_kg[-1] > 0.0


# --- 45. convergence -------------------------------------------------------


def test_results_converge_as_the_ode_tolerance_tightens(initial_state):
    results = []
    for rtol in (1e-6, 1e-8, 1e-10):
        result = simulate_blowdown(
            TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
            SEA_LEVEL_PA, n_report=51, rtol=rtol, atol=rtol * 1e-2,
        )
        results.append(
            (result.burn_time_s, result.tank_pressure_pa[-1],
             result.chamber_pressure_pa[-1], result.total_impulse_n_s)
        )
    tight = results[-1]
    for loose in results[:-1]:
        for a, b in zip(loose, tight, strict=True):
            assert a == pytest.approx(b, rel=1e-3)


def test_results_are_insensitive_to_the_reporting_grid(initial_state):
    coarse = simulate_blowdown(
        TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=26,
    )
    fine = simulate_blowdown(
        TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=401,
    )
    assert coarse.burn_time_s == pytest.approx(fine.burn_time_s, rel=1e-12)
    assert coarse.total_impulse_n_s == pytest.approx(fine.total_impulse_n_s, rel=1e-9)


def test_results_converge_as_the_feed_root_tolerance_tightens(initial_state):
    reference = simulate_blowdown(
        TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=51, feed_xtol=1e-14,
    )
    loose = simulate_blowdown(
        TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=51, feed_xtol=1e-9,
    )
    assert loose.burn_time_s == pytest.approx(reference.burn_time_s, rel=1e-6)
    assert loose.total_impulse_n_s == pytest.approx(reference.total_impulse_n_s, rel=1e-6)


# --- 46, 47. consistency and determinism ----------------------------------


def test_history_matches_standalone_solves_sample_by_sample(nominal):
    for index in (0, 60, 140):
        tank_state = TANK.state_from_mass_and_energy(
            float(nominal.tank_total_mass_kg[index]),
            float(nominal.tank_internal_energy_j[index]),
        )
        standalone = solve_feed_coupling(
            tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
            float(nominal.port_radius_m[index]),
        )
        assert nominal.oxidizer_mass_flow_kg_s[index] == pytest.approx(
            standalone.oxidizer_mass_flow_kg_s, rel=1e-9
        )
        assert nominal.chamber_pressure_pa[index] == pytest.approx(
            standalone.chamber_pressure_pa, rel=1e-9
        )


def test_all_history_arrays_share_one_length(nominal):
    length = nominal.time_s.size
    for array in (
        nominal.port_radius_m, nominal.tank_temperature_k, nominal.tank_pressure_pa,
        nominal.tank_liquid_mass_kg, nominal.oxidizer_mass_flow_kg_s,
        nominal.fuel_mass_flow_kg_s, nominal.chamber_pressure_pa, nominal.thrust_n,
        nominal.cumulative_impulse_n_s, nominal.specific_impulse_s,
    ):
        assert array.size == length
    assert len(nominal.feed_status) == length


def test_repeated_runs_are_bit_identical(initial_state):
    first = simulate_blowdown(
        TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=51,
    )
    second = simulate_blowdown(
        TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=51,
    )
    assert np.array_equal(first.thrust_n, second.thrust_n)
    assert np.array_equal(first.tank_pressure_pa, second.tank_pressure_pa)
    assert first.burn_time_s == second.burn_time_s


# --- 48. the Milestone 3 reference must be untouched ----------------------


def test_milestone_3_prescribed_flow_reference_is_unchanged():
    """The frozen comparison baseline must still produce its published numbers."""
    transient = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(0.100, 120.0), n_report=801
    )
    history = couple_transient_to_performance(transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA)
    assert transient.burnout_time_s == pytest.approx(41.740618, rel=1e-6)
    assert history.thrust_n[0] == pytest.approx(310.048, rel=1e-5)
    assert history.peak_thrust_n == pytest.approx(333.528, rel=1e-5)
    assert history.total_impulse_n_s == pytest.approx(13525.65, rel=1e-5)
    assert history.equivalent_specific_impulse_s == pytest.approx(227.103, rel=1e-5)


def test_coupling_reverses_the_milestone_3_thrust_trend(nominal):
    """The headline Milestone 4 result, pinned as a test.

    Milestone 3 held the oxidizer flow constant and predicted thrust RISING
    through the burn.  With the tank coupled, the falling feed pressure wins and
    thrust FALLS instead.  This is a qualitative sign reversal, not a small
    correction.
    """
    transient = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(0.100, 120.0), n_report=201
    )
    m3 = couple_transient_to_performance(transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA)
    assert m3.thrust_n[-1] > m3.thrust_n[0]
    assert nominal.thrust_n[-1] < nominal.thrust_n[0]


# --- input validation -----------------------------------------------------


@pytest.mark.parametrize("bad", [1, 0, -5])
def test_invalid_report_count_is_rejected(bad, initial_state):
    with pytest.raises(ValueError):
        simulate_blowdown(
            TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
            SEA_LEVEL_PA, n_report=bad,
        )


@pytest.mark.parametrize("bad", [-1.0, math.nan, math.inf])
def test_invalid_ambient_pressure_is_rejected(bad, initial_state):
    with pytest.raises(ValueError):
        simulate_blowdown(
            TANK, initial_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, bad,
        )


def test_result_is_immutable(nominal):
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        nominal.termination = BlowdownTermination.LIQUID_DEPLETED


def test_higher_tank_temperature_gives_more_initial_thrust():
    thrusts = []
    for temperature in (283.15, 293.15, 303.15):
        state = TANK.initial_state(temperature, 0.80)
        result = simulate_blowdown(
            TANK, state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE,
            SEA_LEVEL_PA, t_end_s=1.0, n_report=3,
        )
        thrusts.append(result.thrust_n[0])
    assert all(later > earlier for earlier, later in itertools.pairwise(thrusts))
