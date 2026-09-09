"""Independent verification of the coupled chamber + nozzle + thrust model.

Independence strategy
---------------------
* The reference thrust is rebuilt here from raw ``math`` primitives: the area-Mach
  relation is inverted with an independent bisection, the isentropic ratios are
  written out longhand, and thrust is assembled from the NASA Glenn equation.
* Two independent thrust routes are compared -- the exit-state form
  ``F = m_dot V_e + (p_e - p_a) A_e`` and the coefficient form
  ``F = C_F p_c A_t`` -- neither of which is used to produce the other.
* A third, exact structural identity is used: with fixed ``epsilon``, ``gamma``
  and ``T_c``, thrust is affine in the total mass flow.
* Transient chamber pressures are checked against the Milestone 2 *closed-form*
  radius solution, re-derived here, so the numerical integrator is bypassed.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from hybrid_rocket_motor.chamber import (
    ILLUSTRATIVE_N2O_HTPB_COMBUSTION,
    ChamberStatus,
    CombustionProperties,
)
from hybrid_rocket_motor.geometry import REPRESENTATIVE_GRAIN
from hybrid_rocket_motor.nozzle import (
    REFERENCE_NOZZLE,
    ExpansionRegime,
    NozzleGeometry,
)
from hybrid_rocket_motor.operating_point import OperatingPoint
from hybrid_rocket_motor.performance import (
    STANDARD_GRAVITY_M_S2,
    couple_transient_to_performance,
    solve_performance_point,
)
from hybrid_rocket_motor.regression import REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.transient import (
    ConstantOxidizerFlow,
    PiecewiseConstantOxidizerFlow,
    simulate_transient,
)

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
COMBUSTION = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
NOZZLE = REFERENCE_NOZZLE
SEA_LEVEL_PA = 101325.0

M_DOT_OX = 0.100
M_DOT_F = 0.039776644

# Milestone 1/2 constants, re-declared independently.
R0 = 0.020
R_OUTER = 0.045
RHO_F = 930.0
LENGTH = 0.40
A_SI = 1.0e-3 * 0.3977 * math.pow(10.0, -0.3667)
N_EXPONENT = 0.3667


# --- independent reference implementation ---------------------------------


def independent_area_ratio(mach: float, gamma: float) -> float:
    k = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (((gamma + 1.0) / 2.0) ** (-k)) / mach * (1.0 + mach * mach * (gamma - 1.0) / 2.0) ** k


def independent_exit_mach(eps: float, gamma: float) -> float:
    """Plain bisection on the supersonic branch -- no scipy, no production code."""
    low, high = 1.0, 2.0
    while independent_area_ratio(high, gamma) < eps:
        high *= 2.0
    for _ in range(400):
        mid = 0.5 * (low + high)
        if independent_area_ratio(mid, gamma) < eps:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def independent_performance(m_ox: float, m_f: float, p_a: float, nozzle=NOZZLE,
                            combustion=COMBUSTION) -> dict[str, float]:
    """Rebuild the whole chain longhand from the published relations."""
    gamma = combustion.gamma
    r_gas = combustion.gas_constant_j_kg_k
    a_t = math.pi * nozzle.throat_diameter_m**2 / 4.0
    a_e = math.pi * nozzle.exit_diameter_m**2 / 4.0
    eps = a_e / a_t
    m_total = m_ox + m_f
    p_c = m_total * combustion.c_star_m_s / a_t
    mach = independent_exit_mach(eps, gamma)
    p_e = p_c * math.pow(1.0 + (gamma - 1.0) * mach * mach / 2.0, -gamma / (gamma - 1.0))
    t_e = combustion.chamber_temperature_k / (1.0 + (gamma - 1.0) * mach * mach / 2.0)
    v_e = mach * math.sqrt(gamma * r_gas * t_e)
    momentum = m_total * v_e
    pressure = (p_e - p_a) * a_e
    thrust = momentum + pressure
    return {
        "m_total": m_total,
        "p_c": p_c,
        "mach": mach,
        "p_e": p_e,
        "t_e": t_e,
        "v_e": v_e,
        "momentum": momentum,
        "pressure": pressure,
        "thrust": thrust,
        "cf": thrust / (p_c * a_t),
        "isp": thrust / (m_total * STANDARD_GRAVITY_M_S2),
    }


def analytic_radius(elapsed_s, m_dot_ox: float, r_initial: float = R0):
    """Milestone 2 closed form, re-derived here."""
    k = A_SI * (m_dot_ox / math.pi) ** N_EXPONENT
    power = 2.0 * N_EXPONENT + 1.0
    return (r_initial**power + power * k * np.asarray(elapsed_s, dtype=float)) ** (1.0 / power)


def analytic_fuel_flow(radius, m_dot_ox: float):
    """m_dot_f = rho_f (2 pi r L) a (m_dot_ox/(pi r^2))^n."""
    radius = np.asarray(radius, dtype=float)
    flux = m_dot_ox / (math.pi * radius**2)
    r_dot = A_SI * flux**N_EXPONENT
    return RHO_F * 2.0 * math.pi * radius * LENGTH * r_dot


# ---------------------------------------------------------------------------
# 31. Milestone 1 instantaneous state couples correctly
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def reference_point():
    return solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)


def test_reference_point_matches_the_independent_reimplementation(reference_point):
    expected = independent_performance(M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    assert reference_point.total_mass_flow_kg_s == pytest.approx(expected["m_total"], rel=1e-15)
    assert reference_point.chamber_pressure_pa == pytest.approx(expected["p_c"], rel=1e-14)
    assert reference_point.exit_state.mach == pytest.approx(expected["mach"], rel=1e-11)
    assert reference_point.exit_state.pressure_pa == pytest.approx(expected["p_e"], rel=1e-10)
    assert reference_point.exit_state.temperature_k == pytest.approx(expected["t_e"], rel=1e-11)
    assert reference_point.exit_state.velocity_m_s == pytest.approx(expected["v_e"], rel=1e-11)
    assert reference_point.momentum_thrust_n == pytest.approx(expected["momentum"], rel=1e-11)
    assert reference_point.pressure_thrust_n == pytest.approx(expected["pressure"], rel=1e-9)
    assert reference_point.thrust_n == pytest.approx(expected["thrust"], rel=1e-11)
    assert reference_point.thrust_coefficient == pytest.approx(expected["cf"], rel=1e-11)
    assert reference_point.specific_impulse_s == pytest.approx(expected["isp"], rel=1e-11)


def test_reference_point_uses_the_milestone_1_fuel_flow():
    """The M3 reference point is the M1 operating point, unchanged."""
    m1 = OperatingPoint(GRAIN, M_DOT_OX)
    assert m1.fuel_mass_flow_kg_s(LAW) == pytest.approx(M_DOT_F, rel=1e-8)
    point = solve_performance_point(
        COMBUSTION, NOZZLE, M_DOT_OX, m1.fuel_mass_flow_kg_s(LAW), SEA_LEVEL_PA
    )
    assert point.chamber.mixture_ratio == pytest.approx(m1.mixture_ratio(LAW), rel=1e-12)


def test_reference_point_hardcoded_literals(reference_point):
    """p_c 2611424.85 Pa, M_e 2.619447, V_e 2190.506 m/s, F 310.048 N."""
    assert reference_point.chamber_pressure_pa == pytest.approx(2611424.85, rel=1e-8)
    assert reference_point.exit_state.mach == pytest.approx(2.6194467766658, rel=1e-11)
    assert reference_point.exit_state.pressure_pa == pytest.approx(113631.85, rel=1e-7)
    assert reference_point.exit_state.temperature_k == pytest.approx(1541.974, rel=1e-6)
    assert reference_point.exit_state.velocity_m_s == pytest.approx(2190.506, rel=1e-6)
    assert reference_point.thrust_n == pytest.approx(310.0479, rel=1e-6)
    assert reference_point.thrust_coefficient == pytest.approx(1.511685, rel=1e-6)
    assert reference_point.specific_impulse_s == pytest.approx(226.190, rel=1e-5)


# --- 18, 19, 20, 21, 22. thrust decomposition and identities --------------


def test_thrust_is_exactly_momentum_plus_pressure(reference_point):
    assert reference_point.thrust_n == (
        reference_point.momentum_thrust_n + reference_point.pressure_thrust_n
    )


def test_momentum_thrust_equals_mass_flow_times_exit_velocity(reference_point):
    assert reference_point.momentum_thrust_n == pytest.approx(
        reference_point.total_mass_flow_kg_s * reference_point.exit_state.velocity_m_s, rel=1e-15
    )


def test_pressure_thrust_equals_pressure_difference_times_exit_area(reference_point):
    assert reference_point.pressure_thrust_n == pytest.approx(
        (reference_point.exit_state.pressure_pa - SEA_LEVEL_PA) * NOZZLE.exit_area_m2, rel=1e-14
    )


def test_route_a_and_route_b_thrust_agree(reference_point):
    """Route A: exit-state form. Route B: F = C_F p_c A_t."""
    route_b = (
        reference_point.thrust_coefficient
        * reference_point.chamber_pressure_pa
        * NOZZLE.throat_area_m2
    )
    assert route_b == pytest.approx(reference_point.thrust_n, rel=1e-13)


def test_thrust_coefficient_identity(reference_point):
    assert reference_point.thrust_coefficient == pytest.approx(
        reference_point.thrust_n / (reference_point.chamber_pressure_pa * NOZZLE.throat_area_m2),
        rel=1e-15,
    )
    assert 1.0 < reference_point.thrust_coefficient < 2.0


def test_specific_impulse_identity(reference_point):
    assert reference_point.specific_impulse_s == pytest.approx(
        reference_point.thrust_n
        / (reference_point.total_mass_flow_kg_s * STANDARD_GRAVITY_M_S2),
        rel=1e-15,
    )
    assert reference_point.specific_impulse_s > 0.0


def test_thrust_is_affine_in_total_mass_flow():
    """F = m_dot [V_e + c* eps (p_e/p_c)] - p_a A_e exactly, for fixed properties.

    A structural identity that follows from p_c being proportional to m_dot while
    M_e, T_e and V_e stay fixed. Verified against three separate flow levels.
    """
    ref = solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    slope = (
        ref.exit_state.velocity_m_s
        + COMBUSTION.c_star_m_s * NOZZLE.expansion_ratio * ref.exit_state.pressure_ratio
    )
    intercept = -SEA_LEVEL_PA * NOZZLE.exit_area_m2
    for m_ox, m_f in ((0.05, 0.02), (0.10, 0.04), (0.18, 0.05)):
        point = solve_performance_point(COMBUSTION, NOZZLE, m_ox, m_f, SEA_LEVEL_PA)
        assert point.thrust_n == pytest.approx(slope * (m_ox + m_f) + intercept, rel=1e-12)


def test_exit_velocity_is_independent_of_mass_flow():
    low = solve_performance_point(COMBUSTION, NOZZLE, 0.05, 0.02, SEA_LEVEL_PA)
    high = solve_performance_point(COMBUSTION, NOZZLE, 0.20, 0.06, SEA_LEVEL_PA)
    assert low.exit_state.velocity_m_s == pytest.approx(high.exit_state.velocity_m_s, rel=1e-15)
    assert low.exit_state.temperature_k == pytest.approx(high.exit_state.temperature_k, rel=1e-15)
    assert low.exit_state.mach == pytest.approx(high.exit_state.mach, rel=1e-15)


# --- 26. ambient-pressure behaviour ---------------------------------------


def test_lower_ambient_pressure_increases_thrust_at_a_fixed_chamber_state():
    thrusts = [
        solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, p_a).thrust_n
        for p_a in (101325.0, 54000.0, 20000.0, 0.0)
    ]
    assert all(later > earlier for earlier, later in itertools.pairwise(thrusts))


def test_vacuum_thrust_exceeds_sea_level_thrust_by_pa_times_exit_area():
    sea = solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    vacuum = solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, 0.0)
    assert vacuum.thrust_n - sea.thrust_n == pytest.approx(
        SEA_LEVEL_PA * NOZZLE.exit_area_m2, rel=1e-12
    )


def test_momentum_thrust_is_unaffected_by_ambient_pressure():
    sea = solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    vacuum = solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, 0.0)
    assert sea.momentum_thrust_n == pytest.approx(vacuum.momentum_thrust_n, rel=1e-15)


def test_reference_nozzle_pressure_thrust_is_positive_at_sea_level(reference_point):
    """Slightly underexpanded at the reference condition, so the term helps."""
    assert reference_point.pressure_thrust_n > 0.0
    assert reference_point.expansion_regime is ExpansionRegime.UNDEREXPANDED


def test_a_large_expansion_ratio_gives_negative_pressure_thrust_at_sea_level():
    wide = NozzleGeometry.from_throat_and_expansion_ratio(0.010, 20.0)
    point = solve_performance_point(COMBUSTION, wide, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    assert point.pressure_thrust_n < 0.0
    assert point.expansion_regime is ExpansionRegime.OVEREXPANDED


def test_invalid_ambient_pressure_is_rejected():
    for bad in (-1.0, math.nan, math.inf):
        with pytest.raises(ValueError):
            solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, bad)


# --- idle handling ---------------------------------------------------------


def test_idle_point_produces_no_thrust_and_no_nozzle_state():
    idle = solve_performance_point(COMBUSTION, NOZZLE, 0.0, 0.0, SEA_LEVEL_PA)
    assert idle.chamber.status is ChamberStatus.IDLE_NO_FLOW
    assert not idle.is_firing
    assert idle.chamber_pressure_pa == 0.0
    assert idle.thrust_n == 0.0
    assert idle.momentum_thrust_n == 0.0
    assert idle.pressure_thrust_n == 0.0
    assert idle.exit_state is None
    assert idle.expansion_regime is None
    assert math.isnan(idle.thrust_coefficient)
    assert math.isnan(idle.specific_impulse_s)


def test_idle_thrust_is_zero_not_minus_ambient_times_exit_area():
    """A nozzle full of still ambient gas produces no thrust, not -p_a A_e."""
    idle = solve_performance_point(COMBUSTION, NOZZLE, 0.0, 0.0, SEA_LEVEL_PA)
    assert idle.thrust_n == 0.0
    assert idle.thrust_n != pytest.approx(-SEA_LEVEL_PA * NOZZLE.exit_area_m2)


# ---------------------------------------------------------------------------
# Transient coupling
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def nominal_history():
    transient = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M_DOT_OX, 120.0), n_report=401
    )
    return couple_transient_to_performance(transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA)


PIECEWISE_STEPS = ((10.0, 0.100), (6.0, 0.0), (10.0, 0.150), (14.0, 0.100))


@pytest.fixture(scope="module")
def piecewise_history():
    schedule = PiecewiseConstantOxidizerFlow.from_durations(PIECEWISE_STEPS)
    transient = simulate_transient(GRAIN, LAW, schedule, n_report=301)
    return couple_transient_to_performance(transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA)


# --- 32. transient t = 0 reproduces the static result ---------------------


def test_transient_start_reproduces_the_static_reference(nominal_history, reference_point):
    """The reference point uses the *published* (rounded) M1 fuel flow, 0.039776644,
    while the transient carries the exact value 0.03977664393870721.  Tolerances
    here are therefore set by that ~1e-9 difference in the input, not by the model."""
    assert nominal_history.time_s[0] == 0.0
    assert nominal_history.total_mass_flow_kg_s[0] == pytest.approx(
        reference_point.total_mass_flow_kg_s, rel=1e-9
    )
    assert nominal_history.chamber_pressure_pa[0] == pytest.approx(
        reference_point.chamber_pressure_pa, rel=1e-9
    )
    assert nominal_history.thrust_n[0] == pytest.approx(reference_point.thrust_n, rel=1e-9)
    assert nominal_history.thrust_coefficient[0] == pytest.approx(
        reference_point.thrust_coefficient, rel=1e-8
    )
    assert nominal_history.specific_impulse_s[0] == pytest.approx(
        reference_point.specific_impulse_s, rel=1e-9
    )


# --- 33. numerical transient vs analytic-radius chamber pressure ----------


def test_transient_chamber_pressure_matches_the_analytic_radius_route(nominal_history):
    """Bypass the ODE integrator: build p_c(t) from the M2 closed-form radius."""
    times = nominal_history.time_s
    radius = analytic_radius(times, M_DOT_OX)
    fuel = analytic_fuel_flow(radius, M_DOT_OX)
    expected = (M_DOT_OX + fuel) * COMBUSTION.c_star_m_s / NOZZLE.throat_area_m2
    assert np.allclose(nominal_history.chamber_pressure_pa, expected, rtol=1e-8)


def test_transient_thrust_matches_the_analytic_radius_route(nominal_history):
    times = nominal_history.time_s
    fuel = analytic_fuel_flow(analytic_radius(times, M_DOT_OX), M_DOT_OX)
    expected = np.array(
        [
            independent_performance(M_DOT_OX, float(f), SEA_LEVEL_PA)["thrust"]
            for f in fuel[:: max(1, len(fuel) // 25)]
        ]
    )
    actual = nominal_history.thrust_n[:: max(1, len(fuel) // 25)]
    assert np.allclose(actual, expected, rtol=1e-8)


def test_exact_c_star_identity_holds_across_the_whole_history(nominal_history):
    residual = (
        nominal_history.chamber_pressure_pa * NOZZLE.throat_area_m2
        - nominal_history.total_mass_flow_kg_s * COMBUSTION.c_star_m_s
    )
    assert np.max(np.abs(residual)) < 1.0e-9


# --- 34. constant-flow trends ---------------------------------------------


def test_constant_flow_chamber_pressure_rises_with_fuel_flow(nominal_history):
    """m_dot_f rises through the burn, so p_c must rise with it."""
    assert np.all(np.diff(nominal_history.fuel_mass_flow_kg_s) > 0.0)
    assert np.all(np.diff(nominal_history.total_mass_flow_kg_s) > 0.0)
    assert np.all(np.diff(nominal_history.chamber_pressure_pa) > 0.0)
    assert np.all(np.diff(nominal_history.thrust_n) > 0.0)


def test_exit_velocity_and_mach_are_constant_through_the_burn(nominal_history):
    for array in (nominal_history.exit_mach, nominal_history.exit_velocity_m_s,
                  nominal_history.exit_temperature_k):
        assert np.max(array) - np.min(array) < 1.0e-9


def test_exit_pressure_tracks_chamber_pressure_by_a_constant_ratio(nominal_history):
    ratio = nominal_history.exit_pressure_pa / nominal_history.chamber_pressure_pa
    assert np.max(ratio) - np.min(ratio) < 1.0e-12


def test_all_firing_states_are_physically_sensible(nominal_history):
    assert np.all(nominal_history.total_mass_flow_kg_s > 0.0)
    assert np.all(nominal_history.chamber_pressure_pa > 0.0)
    assert np.all(nominal_history.exit_mach > 1.0)
    assert np.all(nominal_history.exit_velocity_m_s > 0.0)
    assert np.all(nominal_history.momentum_thrust_n > 0.0)
    assert np.all(nominal_history.thrust_n > 0.0)
    assert np.all(nominal_history.specific_impulse_s > 0.0)


# --- 38. no NaN/inf in valid firing states --------------------------------


def test_no_nan_or_inf_in_firing_states(nominal_history):
    for array in (
        nominal_history.chamber_pressure_pa,
        nominal_history.exit_mach,
        nominal_history.exit_pressure_pa,
        nominal_history.exit_temperature_k,
        nominal_history.exit_velocity_m_s,
        nominal_history.momentum_thrust_n,
        nominal_history.pressure_thrust_n,
        nominal_history.thrust_n,
        nominal_history.thrust_coefficient,
        nominal_history.specific_impulse_s,
    ):
        assert np.all(np.isfinite(array))


# --- 35, 36. total impulse and equivalent specific impulse ----------------


def test_total_impulse_matches_an_independent_trapezoidal_integral(nominal_history):
    times, thrust = nominal_history.time_s, nominal_history.thrust_n
    manual = float(
        np.sum(0.5 * (thrust[1:] + thrust[:-1]) * (times[1:] - times[:-1]))
    )
    assert nominal_history.total_impulse_n_s == pytest.approx(manual, rel=1e-13)
    assert nominal_history.total_impulse_n_s > 0.0


def test_equivalent_specific_impulse_identity(nominal_history):
    expected = nominal_history.total_impulse_n_s / (
        nominal_history.propellant_mass_kg * STANDARD_GRAVITY_M_S2
    )
    assert nominal_history.equivalent_specific_impulse_s == pytest.approx(expected, rel=1e-15)
    assert 100.0 < nominal_history.equivalent_specific_impulse_s < 400.0


def test_propellant_mass_is_oxidizer_plus_fuel(nominal_history):
    transient = nominal_history.transient
    assert nominal_history.propellant_mass_kg == pytest.approx(
        transient.cumulative_oxidizer_mass_kg[-1] + transient.cumulative_fuel_mass_kg[-1],
        rel=1e-15,
    )
    # The whole grain is consumed, and oxidizer is flow times burn time.
    assert transient.cumulative_fuel_mass_kg[-1] == pytest.approx(1.899092759, rel=1e-7)
    assert transient.cumulative_oxidizer_mass_kg[-1] == pytest.approx(
        M_DOT_OX * nominal_history.time_s[-1], rel=1e-9
    )


def test_mean_thrust_is_total_impulse_over_active_burn_time(nominal_history):
    assert nominal_history.active_burn_time_s == pytest.approx(
        nominal_history.time_s[-1], rel=1e-12
    )
    assert nominal_history.mean_thrust_n == pytest.approx(
        nominal_history.total_impulse_n_s / nominal_history.active_burn_time_s, rel=1e-14
    )
    assert (
        nominal_history.thrust_n[0]
        < nominal_history.mean_thrust_n
        < nominal_history.peak_thrust_n
    )


def test_peak_thrust_occurs_at_burnout_for_a_constant_flow_burn(nominal_history):
    assert nominal_history.peak_thrust_n == pytest.approx(
        nominal_history.thrust_n[-1], rel=1e-15
    )


# --- 37. piecewise zero-flow segment --------------------------------------


def test_zero_flow_segment_produces_zero_pressure_and_zero_thrust(piecewise_history):
    coast = piecewise_history.oxidizer_mass_flow_kg_s == 0.0
    assert coast.sum() > 5
    assert np.all(piecewise_history.chamber_pressure_pa[coast] == 0.0)
    assert np.all(piecewise_history.thrust_n[coast] == 0.0)
    assert np.all(piecewise_history.momentum_thrust_n[coast] == 0.0)
    assert np.all(piecewise_history.pressure_thrust_n[coast] == 0.0)
    assert np.all(~piecewise_history.firing[coast])


def test_zero_flow_segment_reports_no_nozzle_state(piecewise_history):
    coast = piecewise_history.oxidizer_mass_flow_kg_s == 0.0
    for array in (
        piecewise_history.exit_mach,
        piecewise_history.exit_pressure_pa,
        piecewise_history.exit_temperature_k,
        piecewise_history.exit_velocity_m_s,
        piecewise_history.thrust_coefficient,
        piecewise_history.specific_impulse_s,
    ):
        assert np.all(np.isnan(array[coast]))
    for index in np.flatnonzero(coast):
        assert piecewise_history.expansion_regime[index] is None


def test_piecewise_firing_states_remain_valid(piecewise_history):
    firing = piecewise_history.firing
    assert np.all(piecewise_history.chamber_pressure_pa[firing] > 0.0)
    assert np.all(piecewise_history.thrust_n[firing] > 0.0)
    assert np.all(np.isfinite(piecewise_history.thrust_n))


def test_piecewise_active_burn_time_excludes_the_coast(piecewise_history):
    assert piecewise_history.active_burn_time_s == pytest.approx(34.0, abs=1e-9)
    assert piecewise_history.time_s[-1] == pytest.approx(40.0, rel=1e-12)


def test_piecewise_total_impulse_is_positive_and_mean_thrust_is_consistent(piecewise_history):
    assert piecewise_history.total_impulse_n_s > 0.0
    assert piecewise_history.mean_thrust_n == pytest.approx(
        piecewise_history.total_impulse_n_s / piecewise_history.active_burn_time_s, rel=1e-14
    )


def test_higher_prescribed_flow_raises_chamber_pressure_and_thrust(piecewise_history):
    """The 0.150 kg/s segment must sit above the 0.100 kg/s segments."""
    times = piecewise_history.time_s
    mid = (times > 17.0) & (times < 25.0)
    early = (times > 1.0) & (times < 9.0)
    assert np.min(piecewise_history.chamber_pressure_pa[mid]) > np.max(
        piecewise_history.chamber_pressure_pa[early]
    )
    assert np.min(piecewise_history.thrust_n[mid]) > np.max(piecewise_history.thrust_n[early])


# --- 39. scalar / history consistency -------------------------------------


def test_history_matches_point_solves_sample_by_sample(nominal_history):
    for index in (0, 37, 150, 400):
        point = solve_performance_point(
            COMBUSTION,
            NOZZLE,
            float(nominal_history.oxidizer_mass_flow_kg_s[index]),
            float(nominal_history.fuel_mass_flow_kg_s[index]),
            SEA_LEVEL_PA,
        )
        assert nominal_history.thrust_n[index] == pytest.approx(point.thrust_n, rel=1e-14)
        assert nominal_history.chamber_pressure_pa[index] == pytest.approx(
            point.chamber_pressure_pa, rel=1e-14
        )


def test_all_history_arrays_share_one_length(nominal_history):
    length = nominal_history.time_s.size
    for array in (
        nominal_history.oxidizer_mass_flow_kg_s,
        nominal_history.fuel_mass_flow_kg_s,
        nominal_history.total_mass_flow_kg_s,
        nominal_history.mixture_ratio,
        nominal_history.chamber_pressure_pa,
        nominal_history.exit_mach,
        nominal_history.exit_pressure_pa,
        nominal_history.exit_temperature_k,
        nominal_history.exit_velocity_m_s,
        nominal_history.momentum_thrust_n,
        nominal_history.pressure_thrust_n,
        nominal_history.thrust_n,
        nominal_history.thrust_coefficient,
        nominal_history.specific_impulse_s,
        nominal_history.firing,
    ):
        assert array.size == length
    assert len(nominal_history.expansion_regime) == length


def test_coupling_is_deterministic(nominal_history):
    repeat = couple_transient_to_performance(
        nominal_history.transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA
    )
    assert np.array_equal(repeat.thrust_n, nominal_history.thrust_n)
    assert repeat.total_impulse_n_s == nominal_history.total_impulse_n_s


def test_history_is_immutable(nominal_history):
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        nominal_history.ambient_pressure_pa = 0.0


def test_coupling_rejects_invalid_ambient_pressure(nominal_history):
    for bad in (-1.0, math.nan):
        with pytest.raises(ValueError):
            couple_transient_to_performance(
                nominal_history.transient, COMBUSTION, NOZZLE, bad
            )


# --- 40. convergence -------------------------------------------------------


def test_integrated_performance_converges_with_solver_tolerance():
    results = []
    for rtol in (1e-5, 1e-7, 1e-9):
        transient = simulate_transient(
            GRAIN,
            LAW,
            ConstantOxidizerFlow(M_DOT_OX, 120.0),
            n_report=801,
            rtol=rtol,
            atol=rtol * 1e-4,
        )
        history = couple_transient_to_performance(transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA)
        results.append(
            (
                transient.burnout_time_s,
                history.peak_thrust_n,
                history.total_impulse_n_s,
                history.equivalent_specific_impulse_s,
            )
        )
    tight = results[-1]
    for loose in results[:-1]:
        for a, b in zip(loose, tight, strict=True):
            assert a == pytest.approx(b, rel=1e-6)


def test_integrated_performance_converges_with_reporting_grid():
    coarse = couple_transient_to_performance(
        simulate_transient(GRAIN, LAW, ConstantOxidizerFlow(M_DOT_OX, 120.0), n_report=201),
        COMBUSTION, NOZZLE, SEA_LEVEL_PA,
    )
    fine = couple_transient_to_performance(
        simulate_transient(GRAIN, LAW, ConstantOxidizerFlow(M_DOT_OX, 120.0), n_report=3201),
        COMBUSTION, NOZZLE, SEA_LEVEL_PA,
    )
    assert coarse.total_impulse_n_s == pytest.approx(fine.total_impulse_n_s, rel=1e-6)
    assert coarse.peak_thrust_n == pytest.approx(fine.peak_thrust_n, rel=1e-9)
    assert coarse.equivalent_specific_impulse_s == pytest.approx(
        fine.equivalent_specific_impulse_s, rel=1e-6
    )


# --- combustion-property sensitivity structure ----------------------------


def test_chamber_pressure_scales_directly_with_prescribed_c_star():
    hotter = COMBUSTION.replace(c_star_m_s=1.10 * COMBUSTION.c_star_m_s)
    base = solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    scaled = solve_performance_point(hotter, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    assert scaled.chamber_pressure_pa == pytest.approx(
        1.10 * base.chamber_pressure_pa, rel=1e-13
    )


def test_exit_velocity_scales_as_sqrt_chamber_temperature():
    hotter = COMBUSTION.replace(chamber_temperature_k=4.0 * COMBUSTION.chamber_temperature_k)
    base = solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    scaled = solve_performance_point(hotter, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    assert scaled.exit_state.velocity_m_s == pytest.approx(
        2.0 * base.exit_state.velocity_m_s, rel=1e-13
    )


def test_gamma_changes_the_exit_mach_number():
    softer = CombustionProperties(
        c_star_m_s=COMBUSTION.c_star_m_s,
        gamma=1.15,
        chamber_temperature_k=COMBUSTION.chamber_temperature_k,
        molar_mass_kg_mol=COMBUSTION.molar_mass_kg_mol,
        label="gamma sensitivity",
        reference="ILLUSTRATIVE one-factor sensitivity case",
    )
    base = solve_performance_point(COMBUSTION, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    varied = solve_performance_point(softer, NOZZLE, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    assert varied.exit_state.mach != pytest.approx(base.exit_state.mach, rel=1e-6)
    assert varied.chamber_pressure_pa == pytest.approx(base.chamber_pressure_pa, rel=1e-15)


# --- permanent scope guard -------------------------------------------------

DEFERRED_PHYSICS_TOKENS = (
    "tank",
    "blowdown",
    "vapour",
    "vapor",
    "saturation",
    "injector",
    "feedline",
    "feed_line",
    "feed_system",
    "cea",
    "equilibrium",
    "ignition",
    "instability",
    "structural",
    "contour",
)


def test_package_exports_no_deferred_physics():
    """Milestone 3 may compute chamber pressure, nozzle flow and thrust.

    It may NOT compute anything from the deferred list: tank state, blowdown,
    vapour pressure, injector flow, feed-system losses, equilibrium chemistry,
    ignition transients, combustion instability, structural sizing or nozzle
    contours.  This test pins that boundary so a later milestone cannot cross it
    by accident.
    """
    import hybrid_rocket_motor as package

    exported = {name.lower() for name in dir(package)}
    leaked = sorted(
        token for token in DEFERRED_PHYSICS_TOKENS if any(token in name for name in exported)
    )
    assert leaked == [], f"deferred physics leaked into the public API: {leaked}"


def test_oxidizer_flow_is_never_derived_from_chamber_pressure():
    """The coupling must stay one-way: p_c is an output, never an input to flow.

    Solving at two very different throat areas (hence very different chamber
    pressures) must leave the prescribed flows completely untouched.
    """
    narrow = NozzleGeometry.from_throat_and_expansion_ratio(0.006, 4.0)
    wide = NozzleGeometry.from_throat_and_expansion_ratio(0.016, 4.0)
    a = solve_performance_point(COMBUSTION, narrow, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    b = solve_performance_point(COMBUSTION, wide, M_DOT_OX, M_DOT_F, SEA_LEVEL_PA)
    assert a.chamber_pressure_pa > 5.0 * b.chamber_pressure_pa
    assert a.chamber.oxidizer_mass_flow_kg_s == b.chamber.oxidizer_mass_flow_kg_s == M_DOT_OX
    assert a.chamber.fuel_mass_flow_kg_s == b.chamber.fuel_mass_flow_kg_s == M_DOT_F
    assert a.total_mass_flow_kg_s == b.total_mass_flow_kg_s
