"""Independent verification of the Milestone 2 transient port-regression model.

Independence strategy
---------------------
* The closed-form constant-flow solution is **re-derived here** from the source
  coefficient printed in the paper (0.3977 mm/s, g/(cm^2 s)), converted with its
  own arithmetic, and never read back from the production law object.
* Piecewise-flow trajectories are checked against the analytic solution *chained
  segment by segment*, which the production solver never does -- it integrates
  numerically.
* The ``n = 0.5`` special case is checked against the published exact solution of
  Karabeyoglu, Cantwell & Zilliac (2007), in which the squared port diameter is
  linear in time.
* Mass conservation is checked against geometry formed independently of the
  integrator's own accumulator states.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from hybrid_rocket_motor.geometry import REPRESENTATIVE_GRAIN, GrainGeometry
from hybrid_rocket_motor.operating_point import OperatingPoint
from hybrid_rocket_motor.regression import (
    REZAEI_2018_N2O_HTPB,
    MassFluxUnit,
    PowerLawRegressionLaw,
    Provenance,
    RegressionRateUnit,
)
from hybrid_rocket_motor.transient import (
    CallableOxidizerFlow,
    ConstantOxidizerFlow,
    FlowSegment,
    FluxRangeStatus,
    PiecewiseConstantOxidizerFlow,
    TerminationReason,
    analytic_constant_flow_burnout_time,
    analytic_constant_flow_radius,
    analytic_constant_flow_speed_coefficient,
    remaining_fuel_mass_kg,
    simulate_transient,
)

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB

R0 = GRAIN.initial_port_diameter_m / 2.0          # 0.020 m
R_OUTER = GRAIN.outer_diameter_m / 2.0            # 0.045 m
RHO_F = GRAIN.fuel_density_kg_m3                  # 930 kg/m^3
LENGTH = GRAIN.length_m                           # 0.40 m

M_DOT_REPRESENTATIVE = 0.100

# --- independent re-derivation of the correlation in SI --------------------
# Rezaei et al. (2018) Eq. (10): r_dot [mm/s] = 0.3977 (G_ox [g/(cm^2 s)])^0.3667
SOURCE_A_MM_S = 0.3977
N_EXPONENT = 0.3667
A_SI_REFERENCE = 1.0e-3 * SOURCE_A_MM_S * math.pow(10.0, -N_EXPONENT)


def reference_speed_coefficient(m_dot_ox: float, exponent: float = N_EXPONENT,
                                a_si: float = A_SI_REFERENCE) -> float:
    """K in dr/dt = K r^(-2n), from G_ox = m_dot_ox/(pi r^2)."""
    return a_si * (m_dot_ox / math.pi) ** exponent


def reference_radius(r_initial: float, m_dot_ox: float, elapsed_s, exponent: float = N_EXPONENT,
                     a_si: float = A_SI_REFERENCE):
    """r(t) = [r_0^(2n+1) + (2n+1) K t]^(1/(2n+1))."""
    k = reference_speed_coefficient(m_dot_ox, exponent, a_si)
    power = 2.0 * exponent + 1.0
    return (r_initial**power + power * k * np.asarray(elapsed_s, dtype=float)) ** (1.0 / power)


def reference_burnout_time(r_initial: float, r_outer: float, m_dot_ox: float,
                           exponent: float = N_EXPONENT, a_si: float = A_SI_REFERENCE) -> float:
    """t_burn = [r_outer^(2n+1) - r_0^(2n+1)] / [(2n+1) K]."""
    k = reference_speed_coefficient(m_dot_ox, exponent, a_si)
    power = 2.0 * exponent + 1.0
    return (r_outer**power - r_initial**power) / (power * k)


def reference_radius_after_schedule(r_initial: float, steps) -> float:
    """Chain the closed-form solution across a piecewise-constant schedule."""
    radius = r_initial
    for duration, flow in steps:
        radius = float(reference_radius(radius, flow, duration))
    return radius


# ---------------------------------------------------------------------------
# Shared solutions (module-scoped so the heavier integrations run once)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def constant_run():
    """Representative constant-flow burn, integrated through to burnout."""
    return simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 120.0), n_report=401
    )


@pytest.fixture(scope="module")
def short_run():
    """A short window that stops well before burnout."""
    return simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 10.0), n_report=201
    )


ZERO_FLOW_STEPS = ((8.0, 0.100), (5.0, 0.0), (7.0, 0.150))


@pytest.fixture(scope="module")
def piecewise_run():
    """Nominal flow, a zero-flow coast, then a higher flow."""
    schedule = PiecewiseConstantOxidizerFlow.from_durations(ZERO_FLOW_STEPS)
    return simulate_transient(GRAIN, LAW, schedule, n_report=301)


# ---------------------------------------------------------------------------
# 1. initial state reproduces Milestone 1 instantaneous bookkeeping
# ---------------------------------------------------------------------------


def test_initial_state_reproduces_milestone_1_bookkeeping(constant_run):
    point = OperatingPoint(GRAIN, M_DOT_REPRESENTATIVE, GRAIN.initial_port_diameter_m)
    assert constant_run.time_s[0] == 0.0
    assert constant_run.port_radius_m[0] == pytest.approx(R0, rel=1e-15)
    assert constant_run.port_area_m2[0] == pytest.approx(point.port_area_m2, rel=1e-14)
    assert constant_run.oxidizer_mass_flux_si[0] == pytest.approx(
        point.oxidizer_mass_flux_si, rel=1e-14
    )
    assert constant_run.regression_rate_m_s[0] == pytest.approx(
        point.regression_rate_si(LAW), rel=1e-14
    )
    assert constant_run.burning_area_m2[0] == pytest.approx(point.burning_area_m2, rel=1e-14)
    assert constant_run.fuel_mass_flow_kg_s[0] == pytest.approx(
        point.fuel_mass_flow_kg_s(LAW), rel=1e-14
    )
    assert constant_run.mixture_ratio[0] == pytest.approx(point.mixture_ratio(LAW), rel=1e-14)


def test_initial_state_matches_milestone_1_published_headline_values(constant_run):
    assert constant_run.oxidizer_mass_flux_si[0] == pytest.approx(79.57747154594766, rel=1e-12)
    assert constant_run.regression_rate_mm_s[0] == pytest.approx(0.8508937503295065, rel=1e-12)
    assert constant_run.fuel_mass_flow_kg_s[0] == pytest.approx(0.039776643938707, rel=1e-12)
    assert constant_run.mixture_ratio[0] == pytest.approx(2.5140381414302, rel=1e-12)


def test_initial_remaining_fuel_matches_loaded_grain(constant_run):
    assert constant_run.remaining_fuel_mass_kg[0] == pytest.approx(
        GRAIN.initial_fuel_mass_kg, rel=1e-14
    )
    assert constant_run.cumulative_fuel_mass_kg[0] == 0.0
    assert constant_run.cumulative_oxidizer_mass_kg[0] == 0.0


def test_transient_fuel_mass_helper_agrees_with_milestone_1_geometry():
    """The M2 closed-interval formula must agree with M1's on the open interior."""
    for diameter in (0.040, 0.050, 0.070, 0.089):
        assert remaining_fuel_mass_kg(GRAIN, diameter / 2.0) == pytest.approx(
            GRAIN.fuel_mass_kg(diameter), rel=1e-13
        )


def test_transient_fuel_mass_helper_is_zero_at_burnout():
    assert remaining_fuel_mass_kg(GRAIN, R_OUTER) == pytest.approx(0.0, abs=1e-18)


# ---------------------------------------------------------------------------
# 2, 3. dr/dt = r_dot and dD/dt = 2 r_dot
# ---------------------------------------------------------------------------


def test_reported_regression_rate_equals_the_law_applied_to_the_reported_state(constant_run):
    """r_dot(t) must equal K r(t)^(-2n), formed independently of the solver."""
    k = reference_speed_coefficient(M_DOT_REPRESENTATIVE)
    expected = k * constant_run.port_radius_m ** (-2.0 * N_EXPONENT)
    assert np.allclose(constant_run.regression_rate_m_s, expected, rtol=1e-12, atol=0.0)


def test_radius_derivative_equals_regression_rate(constant_run):
    """Central-difference dr/dt against the reported r_dot."""
    t = constant_run.time_s
    r = constant_run.port_radius_m
    derivative = (r[2:] - r[:-2]) / (t[2:] - t[:-2])
    assert np.allclose(derivative, constant_run.regression_rate_m_s[1:-1], rtol=2e-5, atol=0.0)


def test_diameter_derivative_equals_twice_regression_rate(constant_run):
    t = constant_run.time_s
    d = constant_run.port_diameter_m
    derivative = (d[2:] - d[:-2]) / (t[2:] - t[:-2])
    assert np.allclose(
        derivative, 2.0 * constant_run.regression_rate_m_s[1:-1], rtol=2e-5, atol=0.0
    )


def test_diameter_is_exactly_twice_the_radius(constant_run):
    assert np.array_equal(constant_run.port_diameter_m, 2.0 * constant_run.port_radius_m)


def test_web_thickness_is_the_remaining_radial_fuel(constant_run):
    assert np.allclose(
        constant_run.web_thickness_m, R_OUTER - constant_run.port_radius_m, rtol=0.0, atol=1e-18
    )


# ---------------------------------------------------------------------------
# 4, 5. closed-form solution and burnout time
# ---------------------------------------------------------------------------


def test_production_speed_coefficient_matches_independent_derivation():
    for m_dot in (0.045, 0.100, 0.150):
        assert analytic_constant_flow_speed_coefficient(LAW, m_dot) == pytest.approx(
            reference_speed_coefficient(m_dot), rel=1e-13
        )


def test_analytic_speed_coefficient_hardcoded_literal():
    """K = 1.709446806e-4 * (0.1/pi)^0.3667 = 4.828928...e-5 SI."""
    assert analytic_constant_flow_speed_coefficient(LAW, 0.100) == pytest.approx(
        4.828928e-5, rel=1e-6
    )


def test_analytic_radius_matches_independent_derivation():
    for m_dot in (0.050, 0.100, 0.150):
        for t in (0.0, 1.0, 7.5, 30.0):
            assert analytic_constant_flow_radius(LAW, R0, m_dot, t) == pytest.approx(
                float(reference_radius(R0, m_dot, t)), rel=1e-13
            )


def test_analytic_radius_satisfies_the_state_equation_by_construction():
    """Differentiating the closed form must return K r^(-2n)."""
    m_dot, k = 0.100, reference_speed_coefficient(0.100)
    step = 1.0e-5
    for t in (2.0, 10.0, 30.0):
        forward = analytic_constant_flow_radius(LAW, R0, m_dot, t + step)
        backward = analytic_constant_flow_radius(LAW, R0, m_dot, t - step)
        derivative = (forward - backward) / (2.0 * step)
        radius = analytic_constant_flow_radius(LAW, R0, m_dot, t)
        assert derivative == pytest.approx(k * radius ** (-2.0 * N_EXPONENT), rel=1e-8)


def test_analytic_radius_starts_at_the_initial_radius():
    assert analytic_constant_flow_radius(LAW, R0, 0.100, 0.0) == pytest.approx(R0, rel=1e-15)


def test_analytic_burnout_time_matches_independent_derivation():
    for m_dot in (0.050, 0.100, 0.150):
        assert analytic_constant_flow_burnout_time(LAW, R0, R_OUTER, m_dot) == pytest.approx(
            reference_burnout_time(R0, R_OUTER, m_dot), rel=1e-13
        )


def test_analytic_burnout_time_hardcoded_literal():
    """(0.045^1.7334 - 0.020^1.7334) / (1.7334 * 4.828928e-5) = 41.7406 s."""
    assert analytic_constant_flow_burnout_time(LAW, R0, R_OUTER, 0.100) == pytest.approx(
        41.7406, rel=1e-4
    )


def test_analytic_radius_at_burnout_time_is_the_outer_radius():
    t_burn = analytic_constant_flow_burnout_time(LAW, R0, R_OUTER, 0.100)
    assert analytic_constant_flow_radius(LAW, R0, 0.100, t_burn) == pytest.approx(
        R_OUTER, rel=1e-13
    )


def test_zero_flow_never_burns_out():
    assert analytic_constant_flow_burnout_time(LAW, R0, R_OUTER, 0.0) == math.inf
    assert analytic_constant_flow_radius(LAW, R0, 0.0, 1000.0) == pytest.approx(R0, rel=1e-15)


def test_half_exponent_case_reduces_to_the_published_exact_solution():
    """Karabeyoglu et al. (2007): for n = 0.5 the squared diameter is linear in t.

    Their exact solution for a flux exponent of 0.5 with no length dependence has
    the form D(t) = [D_i^2 + C t]^(1/2).  Our general closed form must collapse to
    exactly that, which is an independent structural check on the derivation.
    """
    half = PowerLawRegressionLaw(
        coefficient_source_units=1.0e-4,
        exponent=0.5,
        flux_unit=MassFluxUnit.KG_PER_M2_S,
        rate_unit=RegressionRateUnit.M_PER_S,
        label="structural check, n = 0.5",
        provenance=Provenance.ILLUSTRATIVE,
        reference="n = 0.5 reduction check against Karabeyoglu et al. (2007)",
    )
    times = np.array([0.0, 5.0, 12.0, 25.0, 60.0])
    diameters = 2.0 * analytic_constant_flow_radius(half, R0, 0.100, times)
    squared = diameters**2
    slopes = np.diff(squared) / np.diff(times)
    assert np.allclose(slopes, slopes[0], rtol=1e-12)
    k = 1.0e-4 * (0.100 / math.pi) ** 0.5
    assert slopes[0] == pytest.approx(8.0 * k, rel=1e-12)


# ---------------------------------------------------------------------------
# 6, 7. numerical solution against the closed form
# ---------------------------------------------------------------------------


def test_numerical_radius_matches_the_analytic_radius(constant_run):
    expected = reference_radius(R0, M_DOT_REPRESENTATIVE, constant_run.time_s)
    assert np.allclose(constant_run.port_radius_m, expected, rtol=1e-9, atol=0.0)


def test_numerical_burnout_time_matches_the_analytic_burnout_time(constant_run):
    expected = reference_burnout_time(R0, R_OUTER, M_DOT_REPRESENTATIVE)
    assert constant_run.burnout_time_s == pytest.approx(expected, rel=1e-9)


def test_numerical_solution_matches_analytic_for_several_flows():
    for m_dot in (0.050, 0.150):
        run = simulate_transient(GRAIN, LAW, ConstantOxidizerFlow(m_dot, 15.0), n_report=51)
        expected = reference_radius(R0, m_dot, run.time_s)
        assert np.allclose(run.port_radius_m, expected, rtol=1e-9, atol=0.0), m_dot


# ---------------------------------------------------------------------------
# 8-12. monotonicity, positivity, no overshoot
# ---------------------------------------------------------------------------


def test_port_radius_is_nondecreasing(constant_run, piecewise_run):
    for run in (constant_run, piecewise_run):
        assert np.all(np.diff(run.port_radius_m) >= 0.0)


def test_port_radius_strictly_increases_under_positive_flow(constant_run):
    assert np.all(np.diff(constant_run.port_radius_m) > 0.0)


def test_port_diameter_is_nondecreasing(constant_run, piecewise_run):
    for run in (constant_run, piecewise_run):
        assert np.all(np.diff(run.port_diameter_m) >= 0.0)


def test_port_area_is_nondecreasing(constant_run, piecewise_run):
    for run in (constant_run, piecewise_run):
        assert np.all(np.diff(run.port_area_m2) >= 0.0)


def test_remaining_fuel_mass_is_nonincreasing(constant_run, piecewise_run):
    for run in (constant_run, piecewise_run):
        assert np.all(np.diff(run.remaining_fuel_mass_kg) <= 0.0)


def test_remaining_fuel_mass_is_never_materially_negative(constant_run, piecewise_run):
    for run in (constant_run, piecewise_run):
        assert np.all(run.remaining_fuel_mass_kg >= -1.0e-12)


def test_port_never_overshoots_the_outer_radius(constant_run, piecewise_run):
    for run in (constant_run, piecewise_run):
        assert np.all(run.port_radius_m <= R_OUTER)
        assert np.all(run.port_diameter_m <= GRAIN.outer_diameter_m)


def test_web_thickness_is_nonnegative_and_shrinking(constant_run):
    assert np.all(constant_run.web_thickness_m >= -1.0e-15)
    assert np.all(np.diff(constant_run.web_thickness_m) <= 0.0)


# ---------------------------------------------------------------------------
# 13, 14. constant flow: G_ox and r_dot fall as the port grows
# ---------------------------------------------------------------------------


def test_constant_flow_flux_decreases_as_the_port_grows(constant_run):
    assert np.all(np.diff(constant_run.oxidizer_mass_flux_si) < 0.0)


def test_constant_flow_regression_rate_decreases_as_the_port_grows(constant_run):
    assert np.all(np.diff(constant_run.regression_rate_m_s) < 0.0)


def test_constant_flow_flux_equals_flow_over_area(constant_run):
    expected = M_DOT_REPRESENTATIVE / (math.pi * constant_run.port_radius_m**2)
    assert np.allclose(constant_run.oxidizer_mass_flux_si, expected, rtol=1e-14)


def test_burning_area_grows_while_regression_rate_falls(constant_run):
    assert np.all(np.diff(constant_run.burning_area_m2) > 0.0)
    assert constant_run.regression_rate_m_s[-1] < constant_run.regression_rate_m_s[0]


# ---------------------------------------------------------------------------
# 15, 16. fuel mass flow and O/F
# ---------------------------------------------------------------------------


def test_fuel_mass_flow_equation_holds_at_every_reported_time(constant_run):
    expected = RHO_F * 2.0 * math.pi * constant_run.port_radius_m * LENGTH * (
        constant_run.regression_rate_m_s
    )
    assert np.allclose(constant_run.fuel_mass_flow_kg_s, expected, rtol=1e-14)


def test_fuel_mass_flow_matches_milestone_1_at_several_sampled_times(constant_run):
    for index in (0, 50, 150, 300):
        diameter = float(constant_run.port_diameter_m[index])
        if diameter >= GRAIN.outer_diameter_m:
            continue
        point = OperatingPoint(GRAIN, M_DOT_REPRESENTATIVE, diameter)
        assert constant_run.fuel_mass_flow_kg_s[index] == pytest.approx(
            point.fuel_mass_flow_kg_s(LAW), rel=1e-12
        )


def test_mixture_ratio_equation_holds_where_defined(constant_run):
    expected = constant_run.oxidizer_mass_flow_kg_s / constant_run.fuel_mass_flow_kg_s
    assert np.allclose(constant_run.mixture_ratio, expected, rtol=1e-14)


def test_mixture_ratio_matches_milestone_1_at_several_sampled_times(constant_run):
    for index in (0, 50, 150, 300):
        diameter = float(constant_run.port_diameter_m[index])
        if diameter >= GRAIN.outer_diameter_m:
            continue
        point = OperatingPoint(GRAIN, M_DOT_REPRESENTATIVE, diameter)
        assert constant_run.mixture_ratio[index] == pytest.approx(
            point.mixture_ratio(LAW), rel=1e-12
        )


def test_all_positive_flow_fuel_flows_are_positive_and_finite(constant_run):
    assert np.all(constant_run.fuel_mass_flow_kg_s > 0.0)
    assert np.all(np.isfinite(constant_run.fuel_mass_flow_kg_s))


def test_mixture_ratio_is_positive_and_finite_under_positive_flow(constant_run):
    assert np.all(constant_run.mixture_ratio > 0.0)
    assert np.all(np.isfinite(constant_run.mixture_ratio))


def test_no_nan_or_inf_in_positive_flow_states(constant_run):
    for array in (
        constant_run.port_radius_m,
        constant_run.port_area_m2,
        constant_run.oxidizer_mass_flux_si,
        constant_run.regression_rate_m_s,
        constant_run.fuel_mass_flow_kg_s,
        constant_run.mixture_ratio,
        constant_run.remaining_fuel_mass_kg,
        constant_run.cumulative_fuel_mass_kg,
        constant_run.cumulative_oxidizer_mass_kg,
    ):
        assert np.all(np.isfinite(array))


# ---------------------------------------------------------------------------
# 17, 18, 19. mass conservation
# ---------------------------------------------------------------------------


def test_geometric_and_integrated_fuel_mass_agree(constant_run):
    """rho_f L pi (r^2 - r_0^2) formed here, not read from the solver."""
    geometric = RHO_F * LENGTH * math.pi * (constant_run.port_radius_m**2 - R0**2)
    assert np.allclose(constant_run.cumulative_fuel_mass_kg, geometric, rtol=0.0, atol=1e-9)
    assert constant_run.max_absolute_fuel_closure_residual_kg < 1.0e-9
    assert constant_run.max_relative_fuel_closure_residual < 1.0e-8


def test_fuel_consumed_plus_remaining_equals_the_loaded_grain(constant_run):
    total = constant_run.remaining_fuel_mass_kg + constant_run.cumulative_fuel_mass_kg
    assert np.allclose(total, GRAIN.initial_fuel_mass_kg, rtol=0.0, atol=1e-9)


def test_burnout_consumes_the_whole_grain(constant_run):
    assert constant_run.cumulative_fuel_mass_kg[-1] == pytest.approx(
        GRAIN.initial_fuel_mass_kg, rel=1e-8
    )
    assert constant_run.remaining_fuel_mass_kg[-1] == pytest.approx(0.0, abs=1e-9)


def test_piecewise_mass_closure(piecewise_run):
    geometric = RHO_F * LENGTH * math.pi * (piecewise_run.port_radius_m**2 - R0**2)
    assert np.allclose(piecewise_run.cumulative_fuel_mass_kg, geometric, rtol=0.0, atol=1e-10)


def test_constant_flow_oxidizer_mass_is_flow_times_time(constant_run):
    expected = M_DOT_REPRESENTATIVE * constant_run.time_s
    assert np.allclose(constant_run.cumulative_oxidizer_mass_kg, expected, rtol=0.0, atol=1e-11)


def test_piecewise_oxidizer_mass_matches_an_independent_segment_sum(piecewise_run):
    """Integrate the schedule by hand: sum of flow x duration up to each time."""
    boundaries = np.cumsum([duration for duration, _ in ZERO_FLOW_STEPS])
    starts = np.concatenate(([0.0], boundaries[:-1]))
    flows = np.array([flow for _, flow in ZERO_FLOW_STEPS])

    def expected_mass(t: float) -> float:
        total = 0.0
        for start, end, flow in zip(starts, boundaries, flows, strict=True):
            total += flow * max(0.0, min(t, end) - start)
        return total

    expected = np.array([expected_mass(float(t)) for t in piecewise_run.time_s])
    assert np.allclose(piecewise_run.cumulative_oxidizer_mass_kg, expected, rtol=0.0, atol=1e-11)


def test_piecewise_radius_matches_the_chained_analytic_solution(piecewise_run):
    """The solver integrates numerically; here the closed form is chained by hand."""
    final = reference_radius_after_schedule(R0, ZERO_FLOW_STEPS)
    assert piecewise_run.port_radius_m[-1] == pytest.approx(final, rel=1e-9)


# ---------------------------------------------------------------------------
# 20, 21, 22. zero-flow behaviour
# ---------------------------------------------------------------------------


def _zero_flow_mask(run) -> np.ndarray:
    return run.oxidizer_mass_flow_kg_s == 0.0


def test_zero_flow_interval_produces_zero_regression(piecewise_run):
    mask = _zero_flow_mask(piecewise_run)
    assert mask.sum() > 5
    assert np.all(piecewise_run.oxidizer_mass_flux_si[mask] == 0.0)
    assert np.all(piecewise_run.regression_rate_m_s[mask] == 0.0)
    assert np.all(piecewise_run.fuel_mass_flow_kg_s[mask] == 0.0)


def test_zero_flow_interval_leaves_the_radius_unchanged(piecewise_run):
    mask = _zero_flow_mask(piecewise_run)
    radii = piecewise_run.port_radius_m[mask]
    assert np.max(radii) - np.min(radii) < 1.0e-15


def test_zero_flow_interval_leaves_fuel_mass_unchanged(piecewise_run):
    mask = _zero_flow_mask(piecewise_run)
    consumed = piecewise_run.cumulative_fuel_mass_kg[mask]
    assert np.max(consumed) - np.min(consumed) < 1.0e-12


def test_mixture_ratio_is_nan_during_zero_flow(piecewise_run):
    mask = _zero_flow_mask(piecewise_run)
    assert np.all(np.isnan(piecewise_run.mixture_ratio[mask]))


def test_mixture_ratio_is_defined_outside_the_zero_flow_interval(piecewise_run):
    mask = ~_zero_flow_mask(piecewise_run)
    assert np.all(np.isfinite(piecewise_run.mixture_ratio[mask]))
    assert np.all(piecewise_run.mixture_ratio[mask] > 0.0)


def test_oxidizer_mass_does_not_accumulate_during_zero_flow(piecewise_run):
    mask = _zero_flow_mask(piecewise_run)
    used = piecewise_run.cumulative_oxidizer_mass_kg[mask]
    assert np.max(used) - np.min(used) < 1.0e-12


def test_a_fully_zero_flow_run_never_regresses():
    run = simulate_transient(GRAIN, LAW, ConstantOxidizerFlow(0.0, 50.0), n_report=21)
    assert np.all(run.port_radius_m == R0)
    assert np.all(run.regression_rate_m_s == 0.0)
    assert np.all(np.isnan(run.mixture_ratio))
    assert run.termination_reason is TerminationReason.COMPLETED_TIME_WINDOW
    assert run.cumulative_fuel_mass_kg[-1] == pytest.approx(0.0, abs=1e-15)
    assert run.cumulative_oxidizer_mass_kg[-1] == pytest.approx(0.0, abs=1e-15)


# ---------------------------------------------------------------------------
# 23, 24. termination events
# ---------------------------------------------------------------------------


def test_burnout_terminates_the_integration_at_the_outer_boundary(constant_run):
    assert constant_run.termination_reason is TerminationReason.FUEL_DEPLETED
    assert constant_run.burnout_time_s is not None
    assert constant_run.port_radius_m[-1] == pytest.approx(R_OUTER, rel=1e-11)
    assert constant_run.time_s[-1] == pytest.approx(constant_run.burnout_time_s, rel=1e-15)


def test_burnout_stops_before_the_requested_horizon(constant_run):
    assert constant_run.time_s[-1] < 120.0


def test_short_run_completes_the_window_without_burnout(short_run):
    assert short_run.termination_reason is TerminationReason.COMPLETED_TIME_WINDOW
    assert short_run.burnout_time_s is None
    assert short_run.time_s[-1] == pytest.approx(10.0, rel=1e-15)
    assert short_run.port_radius_m[-1] < R_OUTER
    assert short_run.remaining_fuel_mass_kg[-1] > 0.0


def test_short_run_agrees_with_the_analytic_solution(short_run):
    expected = reference_radius(R0, M_DOT_REPRESENTATIVE, short_run.time_s)
    assert np.allclose(short_run.port_radius_m, expected, rtol=1e-10, atol=0.0)


def test_t_end_override_truncates_the_window():
    run = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 100.0), t_end_s=5.0, n_report=11
    )
    assert run.time_s[-1] == pytest.approx(5.0, rel=1e-15)
    assert run.termination_reason is TerminationReason.COMPLETED_TIME_WINDOW


# ---------------------------------------------------------------------------
# 25, 26, 27. input validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad_flow", [-1e-9, -0.1, math.nan, math.inf, -math.inf])
def test_constant_flow_rejects_invalid_mass_flow(bad_flow):
    with pytest.raises(ValueError):
        ConstantOxidizerFlow(bad_flow, 10.0)


@pytest.mark.parametrize("bad_flow", [-0.05, math.nan, math.inf])
def test_flow_segment_rejects_invalid_mass_flow(bad_flow):
    with pytest.raises(ValueError):
        FlowSegment(0.0, 1.0, bad_flow)


@pytest.mark.parametrize(("start", "end"), [(0.0, 0.0), (2.0, 1.0), (-1.0, 5.0)])
def test_flow_segment_rejects_invalid_time_range(start, end):
    with pytest.raises(ValueError):
        FlowSegment(start, end, 0.1)


def test_constant_flow_rejects_nonpositive_window():
    with pytest.raises(ValueError):
        ConstantOxidizerFlow(0.1, 0.0)
    with pytest.raises(ValueError):
        ConstantOxidizerFlow(0.1, t_end_s=1.0, t_start_s=1.0)


def test_constant_flow_rejects_negative_time():
    with pytest.raises(ValueError):
        ConstantOxidizerFlow(0.1, t_end_s=5.0, t_start_s=-1.0)


def test_schedule_rejects_a_gap_between_segments():
    with pytest.raises(ValueError, match="contiguous"):
        PiecewiseConstantOxidizerFlow(
            (FlowSegment(0.0, 5.0, 0.1), FlowSegment(6.0, 10.0, 0.1))
        )


def test_schedule_rejects_overlapping_segments():
    with pytest.raises(ValueError, match="contiguous"):
        PiecewiseConstantOxidizerFlow(
            (FlowSegment(0.0, 5.0, 0.1), FlowSegment(4.0, 10.0, 0.1))
        )


def test_schedule_rejects_out_of_order_segments():
    with pytest.raises(ValueError, match="contiguous"):
        PiecewiseConstantOxidizerFlow(
            (FlowSegment(5.0, 10.0, 0.1), FlowSegment(0.0, 5.0, 0.1))
        )


def test_schedule_rejects_an_empty_sequence():
    with pytest.raises(ValueError):
        PiecewiseConstantOxidizerFlow(())
    with pytest.raises(ValueError):
        PiecewiseConstantOxidizerFlow.from_durations([])


@pytest.mark.parametrize("bad_duration", [0.0, -1.0, math.nan, math.inf])
def test_from_durations_rejects_invalid_durations(bad_duration):
    with pytest.raises(ValueError):
        PiecewiseConstantOxidizerFlow.from_durations([(bad_duration, 0.1)])


def test_schedule_rejects_queries_outside_its_window():
    schedule = PiecewiseConstantOxidizerFlow.from_durations([(5.0, 0.1)])
    with pytest.raises(ValueError):
        schedule.mass_flow_kg_s(6.0)


def test_callable_flow_rejects_a_negative_sample():
    with pytest.raises(ValueError):
        CallableOxidizerFlow(lambda t: 0.1 - 0.02 * t, 10.0)


def test_callable_flow_rejects_a_non_finite_sample():
    with pytest.raises(ValueError):
        CallableOxidizerFlow(lambda t: math.nan, 10.0)


def test_callable_flow_rejects_a_non_callable():
    with pytest.raises(TypeError):
        CallableOxidizerFlow(0.1, 10.0)


def test_solver_rejects_an_invalid_reporting_grid():
    with pytest.raises(ValueError):
        simulate_transient(GRAIN, LAW, ConstantOxidizerFlow(0.1, 5.0), n_report=1)


def test_solver_rejects_a_t_end_outside_the_flow_window():
    with pytest.raises(ValueError):
        simulate_transient(GRAIN, LAW, ConstantOxidizerFlow(0.1, 5.0), t_end_s=9.0)


@pytest.mark.parametrize("bad_diameter", [0.0, -0.01, 0.090, 0.2])
def test_solver_rejects_an_invalid_initial_port_diameter(bad_diameter):
    with pytest.raises(ValueError):
        simulate_transient(
            GRAIN, LAW, ConstantOxidizerFlow(0.1, 5.0), initial_port_diameter_m=bad_diameter
        )


def test_remaining_fuel_helper_rejects_a_radius_beyond_the_grain():
    with pytest.raises(ValueError):
        remaining_fuel_mass_kg(GRAIN, R_OUTER * 1.05)
    with pytest.raises(ValueError):
        remaining_fuel_mass_kg(GRAIN, 0.0)


# ---------------------------------------------------------------------------
# 28. scalar / history consistency
# ---------------------------------------------------------------------------


def test_analytic_radius_scalar_and_array_agree():
    times = np.array([0.0, 3.0, 11.0, 27.5])
    vectorised = analytic_constant_flow_radius(LAW, R0, 0.100, times)
    assert isinstance(vectorised, np.ndarray)
    for t, value in zip(times.tolist(), vectorised.tolist(), strict=True):
        assert value == pytest.approx(analytic_constant_flow_radius(LAW, R0, 0.100, t), rel=1e-15)


def test_analytic_radius_scalar_input_returns_a_float():
    assert isinstance(analytic_constant_flow_radius(LAW, R0, 0.100, 4.0), float)


def test_remaining_fuel_helper_scalar_and_array_agree():
    radii = np.array([0.020, 0.030, 0.044])
    vectorised = remaining_fuel_mass_kg(GRAIN, radii)
    assert isinstance(vectorised, np.ndarray)
    for radius, value in zip(radii.tolist(), vectorised.tolist(), strict=True):
        assert value == pytest.approx(remaining_fuel_mass_kg(GRAIN, radius), rel=1e-15)


def test_history_arrays_all_share_one_length(constant_run):
    length = constant_run.time_s.size
    for array in (
        constant_run.port_radius_m,
        constant_run.port_diameter_m,
        constant_run.port_area_m2,
        constant_run.web_thickness_m,
        constant_run.remaining_fuel_mass_kg,
        constant_run.oxidizer_mass_flow_kg_s,
        constant_run.oxidizer_mass_flux_si,
        constant_run.regression_rate_m_s,
        constant_run.regression_rate_mm_s,
        constant_run.burning_area_m2,
        constant_run.fuel_mass_flow_kg_s,
        constant_run.mixture_ratio,
        constant_run.cumulative_oxidizer_mass_kg,
        constant_run.cumulative_fuel_mass_kg,
    ):
        assert array.size == length
    assert len(constant_run.flux_status) == length


def test_reported_times_are_strictly_increasing(constant_run, piecewise_run):
    for run in (constant_run, piecewise_run):
        assert np.all(np.diff(run.time_s) > 0.0)


def test_regression_rate_in_mm_per_second_is_exactly_1000x(constant_run):
    assert np.allclose(
        constant_run.regression_rate_mm_s, constant_run.regression_rate_m_s * 1.0e3, rtol=1e-15
    )


# ---------------------------------------------------------------------------
# 29. source-flux-range classification
# ---------------------------------------------------------------------------


def test_flux_classification_matches_the_recorded_source_range(constant_run):
    low, high = LAW.valid_flux_range_si
    for flux, flow, status in zip(
        constant_run.oxidizer_mass_flux_si,
        constant_run.oxidizer_mass_flow_kg_s,
        constant_run.flux_status,
        strict=True,
    ):
        if flow == 0.0:
            expected = FluxRangeStatus.ZERO_OXIDIZER_FLOW
        elif flux < low:
            expected = FluxRangeStatus.BELOW_SOURCE_FLUX_RANGE
        elif flux > high:
            expected = FluxRangeStatus.ABOVE_SOURCE_FLUX_RANGE
        else:
            expected = FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE
        assert status is expected


def test_constant_flow_run_starts_in_range_and_ends_below_it(constant_run):
    assert constant_run.flux_status[0] is FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE
    assert constant_run.flux_status[-1] is FluxRangeStatus.BELOW_SOURCE_FLUX_RANGE


def test_within_range_time_fraction_matches_the_analytic_exit_time(constant_run):
    """The trajectory leaves the source range when G_ox falls to 35 kg/(m^2 s)."""
    low = LAW.valid_flux_range_si[0]
    radius_at_exit = math.sqrt(M_DOT_REPRESENTATIVE / (low * math.pi))
    power = 2.0 * N_EXPONENT + 1.0
    k = reference_speed_coefficient(M_DOT_REPRESENTATIVE)
    t_exit = (radius_at_exit**power - R0**power) / (power * k)
    t_burn = reference_burnout_time(R0, R_OUTER, M_DOT_REPRESENTATIVE)
    fraction = constant_run.time_fraction_with_status(FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE)
    assert fraction == pytest.approx(t_exit / t_burn, abs=0.01)
    assert 0.30 < fraction < 0.37


def test_zero_flow_is_classified_separately_from_below_range(piecewise_run):
    statuses = set(piecewise_run.flux_status)
    assert FluxRangeStatus.ZERO_OXIDIZER_FLOW in statuses
    mask = _zero_flow_mask(piecewise_run)
    for index in np.flatnonzero(mask):
        assert piecewise_run.flux_status[index] is FluxRangeStatus.ZERO_OXIDIZER_FLOW


def test_a_law_without_a_source_range_is_reported_as_such():
    unranged = PowerLawRegressionLaw(
        coefficient_source_units=LAW.coefficient_si,
        exponent=LAW.exponent,
        flux_unit=MassFluxUnit.KG_PER_M2_S,
        rate_unit=RegressionRateUnit.M_PER_S,
        label="ILLUSTRATIVE, no declared range",
        provenance=Provenance.ILLUSTRATIVE,
        reference="classification test - no validity range recorded",
    )
    run = simulate_transient(GRAIN, unranged, ConstantOxidizerFlow(0.1, 4.0), n_report=11)
    assert all(s is FluxRangeStatus.NO_SOURCE_RANGE_DECLARED for s in run.flux_status)


def test_time_fractions_over_all_statuses_sum_to_one(constant_run, piecewise_run):
    for run in (constant_run, piecewise_run):
        total = sum(run.time_fraction_with_status(status) for status in FluxRangeStatus)
        assert total == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# 30. convergence
# ---------------------------------------------------------------------------


def test_burnout_time_converges_towards_the_analytic_value():
    expected = reference_burnout_time(R0, R_OUTER, M_DOT_REPRESENTATIVE)
    errors = []
    for rtol in (1e-5, 1e-7, 1e-9):
        run = simulate_transient(
            GRAIN,
            LAW,
            ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 120.0),
            n_report=11,
            rtol=rtol,
            atol=rtol * 1e-4,
        )
        errors.append(abs(run.burnout_time_s - expected) / expected)
    assert errors[-1] <= errors[0]
    assert errors[-1] < 1.0e-9


def test_results_are_insensitive_to_the_reporting_grid():
    """The reporting grid must not influence the physics, only the sampling."""
    coarse = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 120.0), n_report=11
    )
    fine = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 120.0), n_report=2001
    )
    assert coarse.burnout_time_s == pytest.approx(fine.burnout_time_s, rel=1e-12)
    assert coarse.cumulative_fuel_mass_kg[-1] == pytest.approx(
        fine.cumulative_fuel_mass_kg[-1], rel=1e-10
    )


def test_max_step_restriction_does_not_change_the_answer():
    reference = reference_burnout_time(R0, R_OUTER, M_DOT_REPRESENTATIVE)
    run = simulate_transient(
        GRAIN,
        LAW,
        ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 120.0),
        n_report=51,
        max_step=0.25,
    )
    assert run.burnout_time_s == pytest.approx(reference, rel=1e-9)


def test_an_alternative_integrator_reaches_the_same_burnout_time():
    reference = reference_burnout_time(R0, R_OUTER, M_DOT_REPRESENTATIVE)
    run = simulate_transient(
        GRAIN,
        LAW,
        ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 120.0),
        n_report=51,
        method="RK45",
        rtol=1e-11,
        atol=1e-14,
    )
    assert run.burnout_time_s == pytest.approx(reference, rel=1e-8)


def test_mass_closure_residual_shrinks_with_tighter_tolerance():
    residuals = []
    for rtol in (1e-5, 1e-9):
        run = simulate_transient(
            GRAIN,
            LAW,
            ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 30.0),
            n_report=51,
            rtol=rtol,
            atol=rtol * 1e-4,
        )
        residuals.append(run.max_absolute_fuel_closure_residual_kg)
    assert residuals[-1] <= residuals[0]
    assert residuals[-1] < 1.0e-9


# ---------------------------------------------------------------------------
# Callable flow history and miscellaneous
# ---------------------------------------------------------------------------


def test_callable_flow_history_integrates_and_conserves_oxidizer_mass():
    """m_dot_ox(t) = 0.08 + 0.002 t, so the oxidizer mass is 0.08t + 0.001 t^2."""
    history = CallableOxidizerFlow(lambda t: 0.08 + 0.002 * t, 12.0)
    run = simulate_transient(GRAIN, LAW, history, n_report=61)
    expected = 0.08 * run.time_s + 0.001 * run.time_s**2
    assert np.allclose(run.cumulative_oxidizer_mass_kg, expected, rtol=0.0, atol=1e-10)
    assert np.all(np.diff(run.port_radius_m) > 0.0)


def test_callable_flow_matches_constant_flow_when_it_is_constant():
    history = CallableOxidizerFlow(lambda t: M_DOT_REPRESENTATIVE, 12.0)
    run = simulate_transient(GRAIN, LAW, history, n_report=41)
    expected = reference_radius(R0, M_DOT_REPRESENTATIVE, run.time_s)
    assert np.allclose(run.port_radius_m, expected, rtol=1e-9)


def test_schedule_reports_its_interior_breakpoints():
    schedule = PiecewiseConstantOxidizerFlow.from_durations(ZERO_FLOW_STEPS)
    assert schedule.interior_breakpoints_s == (8.0, 13.0)
    assert schedule.t_start_s == 0.0
    assert schedule.t_end_s == 20.0


def test_schedule_flow_lookup_uses_half_open_segments():
    schedule = PiecewiseConstantOxidizerFlow.from_durations(ZERO_FLOW_STEPS)
    assert schedule.mass_flow_kg_s(0.0) == 0.100
    assert schedule.mass_flow_kg_s(7.999) == 0.100
    assert schedule.mass_flow_kg_s(8.0) == 0.0
    assert schedule.mass_flow_kg_s(13.0) == 0.150
    assert schedule.mass_flow_kg_s(20.0) == 0.150


def test_breakpoints_are_sampled_on_both_sides(piecewise_run):
    """A discontinuous prescribed flow must not be drawn as a ramp."""
    for breakpoint_s in (8.0, 13.0):
        index = int(np.searchsorted(piecewise_run.time_s, breakpoint_s))
        assert piecewise_run.time_s[index] == breakpoint_s
        assert piecewise_run.time_s[index] - piecewise_run.time_s[index - 1] < 1.0e-12


def test_a_larger_initial_port_burns_out_sooner():
    wide = simulate_transient(
        GRAIN,
        LAW,
        ConstantOxidizerFlow(M_DOT_REPRESENTATIVE, 200.0),
        initial_port_diameter_m=0.060,
        n_report=21,
    )
    expected = reference_burnout_time(0.030, R_OUTER, M_DOT_REPRESENTATIVE)
    assert wide.burnout_time_s == pytest.approx(expected, rel=1e-9)
    assert wide.burnout_time_s < reference_burnout_time(R0, R_OUTER, M_DOT_REPRESENTATIVE)


def test_higher_flow_burns_out_sooner_than_lower_flow():
    times = []
    for m_dot in (0.050, 0.100, 0.150):
        run = simulate_transient(GRAIN, LAW, ConstantOxidizerFlow(m_dot, 200.0), n_report=11)
        assert run.termination_reason is TerminationReason.FUEL_DEPLETED
        times.append(run.burnout_time_s)
    assert all(later < earlier for earlier, later in itertools.pairwise(times))


def test_a_different_grain_is_handled_consistently():
    grain = GrainGeometry(
        length_m=0.25,
        initial_port_diameter_m=0.030,
        outer_diameter_m=0.080,
        fuel_density_kg_m3=930.0,
    )
    run = simulate_transient(grain, LAW, ConstantOxidizerFlow(0.08, 200.0), n_report=51)
    expected = reference_burnout_time(0.015, 0.040, 0.08)
    assert run.burnout_time_s == pytest.approx(expected, rel=1e-9)
    geometric = 930.0 * 0.25 * math.pi * (run.port_radius_m**2 - 0.015**2)
    assert np.allclose(run.cumulative_fuel_mass_kg, geometric, rtol=0.0, atol=1e-9)


def test_result_object_is_immutable(constant_run):
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        constant_run.burnout_time_s = 1.0
