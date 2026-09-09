"""Independent verification of the implicit variable-property chamber solve.

Independence strategy
---------------------
The converged chamber pressure is checked by re-forming the defining relation
``p_c = m_dot c*(O/F, p_c) / A_t`` from the reported state rather than by trusting
the solver's own residual, and by re-solving the same closure with a completely
different method (a fine scan plus a secant refinement) in one case.  The frozen
Milestone 3 chamber is re-run in the same session to prove it is untouched.
"""

from __future__ import annotations

import itertools
import math

import pytest

from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION, solve_chamber
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.thermochemistry import DEFAULT_C_STAR_EFFICIENCY, default_table
from hybrid_rocket_motor.variable_chamber import (
    VariableChamberStatus,
    solve_variable_chamber,
)

THROAT_AREA = REFERENCE_NOZZLE.throat_area_m2
M3_OXIDIZER_FLOW = 0.100
M3_FUEL_FLOW = 0.039776644


@pytest.fixture(scope="module")
def table():
    return default_table()


# -- the closure actually closes -------------------------------------------------


def test_reference_point_solves(table):
    solution = solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)
    assert solution.status is VariableChamberStatus.SOLVED
    assert solution.is_solved


def test_converged_pressure_satisfies_the_defining_relation(table):
    """Re-form ``p_c = m_dot c*/A_t`` from the reported state, independently."""
    solution = solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)
    state = table.evaluate(solution.mixture_ratio, solution.chamber_pressure_pa)
    implied = solution.total_mass_flow_kg_s * state.c_star_m_s / THROAT_AREA
    assert implied == pytest.approx(solution.chamber_pressure_pa, rel=1e-12)


def test_mixture_ratio_identity_holds_exactly(table):
    solution = solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)
    assert solution.mixture_ratio == pytest.approx(M3_OXIDIZER_FLOW / M3_FUEL_FLOW, rel=1e-15)
    assert solution.total_mass_flow_kg_s == pytest.approx(
        M3_OXIDIZER_FLOW + M3_FUEL_FLOW, rel=1e-15
    )


def test_reported_residuals_are_at_machine_precision(table):
    solution = solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)
    assert abs(solution.residual_pa) < 1e-6
    assert abs(solution.c_star_identity_residual) < 1e-9


def test_independent_root_finder_agrees(table):
    """Re-solve the same closure by scan-and-secant instead of Brent."""
    solution = solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)
    total = M3_OXIDIZER_FLOW + M3_FUEL_FLOW
    of_ratio = M3_OXIDIZER_FLOW / M3_FUEL_FLOW

    def residual(pressure: float) -> float:
        state = table.evaluate(of_ratio, pressure)
        return total * state.c_star_m_s / THROAT_AREA - pressure

    low, high = table.pressure_bounds_pa
    scan = [low + (high - low) * i / 400.0 for i in range(401)]
    bracket = next((a, b) for a, b in itertools.pairwise(scan) if residual(a) * residual(b) <= 0.0)
    a, b = bracket
    for _ in range(200):
        fa, fb = residual(a), residual(b)
        if fb == fa:
            break
        c = b - fb * (b - a) / (fb - fa)
        a, b = b, c
    assert b == pytest.approx(solution.chamber_pressure_pa, rel=1e-9)


# -- the variable model differs from the frozen one in the expected direction ------


def test_frozen_milestone_3_chamber_is_untouched():
    """The Milestone 3 solver must still give its own, unchanged answer."""
    frozen = solve_chamber(
        ILLUSTRATIVE_N2O_HTPB_COMBUSTION, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW
    )
    expected = (
        (M3_OXIDIZER_FLOW + M3_FUEL_FLOW)
        * ILLUSTRATIVE_N2O_HTPB_COMBUSTION.c_star_m_s
        / THROAT_AREA
    )
    assert frozen.chamber_pressure_pa == pytest.approx(expected, rel=1e-15)


def test_variable_pressure_is_lower_than_the_frozen_one_at_this_point(table):
    """The equilibrium c* here is below the illustrative constant, so p_c falls."""
    variable = solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)
    frozen = solve_chamber(
        ILLUSTRATIVE_N2O_HTPB_COMBUSTION, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW
    )
    assert variable.chamber_pressure_pa < frozen.chamber_pressure_pa
    relative = variable.chamber_pressure_pa / frozen.chamber_pressure_pa - 1.0
    assert -0.10 < relative < -0.05


def test_pressure_scales_with_efficiency(table):
    """With the table's weak pressure dependence, p_c is near-proportional to eta."""
    low = solve_variable_chamber(
        table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW, c_star_efficiency=0.80
    )
    high = solve_variable_chamber(
        table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW, c_star_efficiency=1.00
    )
    assert high.chamber_pressure_pa > low.chamber_pressure_pa
    assert high.chamber_pressure_pa / low.chamber_pressure_pa == pytest.approx(1.25, rel=2e-3)


def test_pressure_rises_with_total_mass_flow(table):
    previous = 0.0
    for scale in (0.9, 1.0, 1.1, 1.2):
        solution = solve_variable_chamber(
            table, THROAT_AREA, scale * M3_OXIDIZER_FLOW, scale * M3_FUEL_FLOW
        )
        assert solution.is_solved
        assert solution.chamber_pressure_pa > previous
        previous = solution.chamber_pressure_pa


def test_mixture_ratio_is_unchanged_by_scaling_both_flows(table):
    """Only the ratio enters the chemistry, so scaling both flows must not move it."""
    base = solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)
    scaled = solve_variable_chamber(table, THROAT_AREA, 1.3 * M3_OXIDIZER_FLOW, 1.3 * M3_FUEL_FLOW)
    assert scaled.mixture_ratio == pytest.approx(base.mixture_ratio, rel=1e-14)


# -- explicit failure instead of clamping -----------------------------------------


def test_no_flow_is_reported_as_idle_not_solved(table):
    solution = solve_variable_chamber(table, THROAT_AREA, 0.0, 0.0)
    assert solution.status is VariableChamberStatus.IDLE_NO_FLOW
    assert not solution.is_solved
    assert solution.chamber_pressure_pa == 0.0
    assert solution.thermochemistry is None


def test_mixture_ratio_below_the_table_is_refused_not_clipped(table):
    """A very fuel-rich point must fail loudly rather than be pulled onto the edge."""
    of_low, _ = table.of_bounds
    fuel = 1.0
    oxidizer = 0.5 * of_low * fuel
    solution = solve_variable_chamber(table, THROAT_AREA, oxidizer, fuel)
    assert solution.status is VariableChamberStatus.OF_OUTSIDE_TABLE
    assert not solution.is_solved
    assert solution.mixture_ratio == pytest.approx(oxidizer / fuel)
    assert solution.mixture_ratio < of_low, "the reported O/F must not be clipped into the table"


def test_mixture_ratio_above_the_table_is_refused(table):
    _, of_high = table.of_bounds
    fuel = 0.001
    oxidizer = 2.0 * of_high * fuel
    solution = solve_variable_chamber(table, THROAT_AREA, oxidizer, fuel)
    assert solution.status is VariableChamberStatus.OF_OUTSIDE_TABLE
    assert solution.mixture_ratio > of_high


def test_pressure_above_the_table_is_refused(table):
    """A huge mass flow through this throat would need a pressure the table lacks."""
    solution = solve_variable_chamber(table, THROAT_AREA, 10.0, 4.0)
    assert solution.status is VariableChamberStatus.PRESSURE_OUTSIDE_TABLE
    assert not solution.is_solved


def test_pressure_below_the_table_is_refused(table):
    """A tiny mass flow cannot reach the tabulated pressure floor."""
    solution = solve_variable_chamber(table, THROAT_AREA, 0.0010, 0.00045)
    assert solution.status is VariableChamberStatus.PRESSURE_OUTSIDE_TABLE
    assert not solution.is_solved


def test_a_failed_solve_never_reports_a_thermochemical_state(table):
    for oxidizer, fuel in ((0.0, 0.0), (0.5, 1.0), (10.0, 4.0)):
        solution = solve_variable_chamber(table, THROAT_AREA, oxidizer, fuel)
        assert not solution.is_solved
        assert solution.thermochemistry is None


# -- argument validation ------------------------------------------------------------


@pytest.mark.parametrize("throat_area", [0.0, -1.0, math.nan, math.inf])
def test_invalid_throat_area_is_refused(table, throat_area):
    with pytest.raises(ValueError):
        solve_variable_chamber(table, throat_area, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)


@pytest.mark.parametrize("flow", [-1.0, math.nan, math.inf])
def test_invalid_mass_flows_are_refused(table, flow):
    with pytest.raises(ValueError):
        solve_variable_chamber(table, THROAT_AREA, flow, M3_FUEL_FLOW)
    with pytest.raises(ValueError):
        solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, flow)


def test_default_efficiency_is_the_shared_one(table):
    explicit = solve_variable_chamber(
        table,
        THROAT_AREA,
        M3_OXIDIZER_FLOW,
        M3_FUEL_FLOW,
        c_star_efficiency=DEFAULT_C_STAR_EFFICIENCY,
    )
    implicit = solve_variable_chamber(table, THROAT_AREA, M3_OXIDIZER_FLOW, M3_FUEL_FLOW)
    assert explicit.chamber_pressure_pa == pytest.approx(implicit.chamber_pressure_pa, rel=1e-15)
