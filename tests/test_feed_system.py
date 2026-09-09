"""Independent verification of the coupled injector/chamber closure.

Independence strategy
---------------------
The converged oxidizer flow is re-solved here by a plain bisection written in the
test file, and the fixed point is verified directly by pushing the solved flow
back through the frozen Milestone 1-3 chain by hand.
"""

from __future__ import annotations

import itertools
import math

import pytest

from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION, solve_chamber
from hybrid_rocket_motor.feed_system import FeedStatus, solve_feed_coupling
from hybrid_rocket_motor.geometry import REPRESENTATIVE_GRAIN
from hybrid_rocket_motor.injector import Injector, InjectorModel
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.regression import REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
)

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
COMBUSTION = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
NOZZLE = REFERENCE_NOZZLE
AREA = 3.5e-6
R0 = 0.020


@pytest.fixture(scope="module")
def tank_state():
    return REFERENCE_TANK.initial_state(
        REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION
    )


def hand_chamber_pressure(oxidizer_flow: float, radius: float) -> float:
    """The frozen M1 -> M2 -> M3 chain, rebuilt longhand."""
    flux = oxidizer_flow / (math.pi * radius * radius)
    rate = LAW.regression_rate_si(flux)
    fuel = 930.0 * 2.0 * math.pi * radius * 0.40 * rate
    return (oxidizer_flow + fuel) * COMBUSTION.c_star_m_s / NOZZLE.throat_area_m2


def hand_bisection(tank, injector, radius, iterations: int = 200) -> float:
    """Plain bisection on the same residual -- no scipy, no production solver."""
    low = 0.0
    high = injector.spi_mass_flow_kg_s(tank.liquid_density_kg_m3, tank.pressure_pa)

    def residual(flow: float) -> float:
        p_c = max(hand_chamber_pressure(flow, radius), 87837.4)
        return injector.mass_flow(tank, p_c).mass_flow_kg_s - flow

    for _ in range(iterations):
        mid = 0.5 * (low + high)
        if residual(mid) > 0.0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


# --- 25, 26. the solved flow is the fixed point ---------------------------


def test_solved_flow_is_the_injector_flow_at_the_solved_chamber_pressure(tank_state):
    solution = solve_feed_coupling(
        tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0
    )
    reproduced = Injector(AREA).mass_flow(tank_state, solution.chamber_pressure_pa)
    assert reproduced.mass_flow_kg_s == pytest.approx(
        solution.oxidizer_mass_flow_kg_s, rel=1e-9
    )


def test_solved_chamber_pressure_matches_the_hand_built_chain(tank_state):
    solution = solve_feed_coupling(
        tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0
    )
    assert solution.chamber_pressure_pa == pytest.approx(
        hand_chamber_pressure(solution.oxidizer_mass_flow_kg_s, R0), rel=1e-12
    )
    assert solution.chamber_pressure_pa == pytest.approx(
        solve_chamber(
            COMBUSTION,
            NOZZLE.throat_area_m2,
            solution.oxidizer_mass_flow_kg_s,
            solution.fuel_mass_flow_kg_s,
        ).chamber_pressure_pa,
        rel=1e-14,
    )


def test_residual_is_negligible(tank_state):
    for area in (2.0e-6, 3.5e-6, 6.0e-6):
        solution = solve_feed_coupling(
            tank_state, Injector(area), GRAIN, LAW, COMBUSTION, NOZZLE, R0
        )
        assert abs(solution.residual_kg_s) < 1.0e-10


def test_production_solve_matches_an_independent_bisection(tank_state):
    for area, model in ((3.5e-6, InjectorModel.DYER_NHNE), (3.5e-6, InjectorModel.SPI)):
        injector = Injector(area, model=model)
        solution = solve_feed_coupling(
            tank_state, injector, GRAIN, LAW, COMBUSTION, NOZZLE, R0
        )
        assert solution.oxidizer_mass_flow_kg_s == pytest.approx(
            hand_bisection(tank_state, injector, R0), rel=1e-9
        )


def test_reference_coupled_state_hardcoded_literals(tank_state):
    """3.5 mm^2 Dyer injector at r_p = 20 mm -> 0.1013 kg/s and 26.40 bar."""
    solution = solve_feed_coupling(
        tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0
    )
    assert solution.oxidizer_mass_flow_kg_s == pytest.approx(0.101316, rel=1e-4)
    assert solution.chamber_pressure_pa / 1e5 == pytest.approx(26.396, rel=1e-4)
    assert solution.mixture_ratio == pytest.approx(2.535, rel=1e-3)


# --- the uniqueness argument, verified numerically ------------------------


def test_the_coupled_residual_is_strictly_decreasing(tank_state):
    """The injector characteristic is not monotone, but the residual is.

    Checks the loop-gain argument in ``feed_system`` directly: even across the
    shallow Dyer maximum the residual falls monotonically, so the bracketed root
    is unique.
    """
    injector = Injector(AREA)
    upper = injector.spi_mass_flow_kg_s(tank_state.liquid_density_kg_m3, tank_state.pressure_pa)
    flows = [upper * i / 60.0 for i in range(1, 61)]
    residuals = [
        injector.mass_flow(
            tank_state, max(hand_chamber_pressure(flow, R0), 87837.4)
        ).mass_flow_kg_s
        - flow
        for flow in flows
    ]
    assert all(later < earlier for earlier, later in itertools.pairwise(residuals))


# --- 27, 28. pressure margin ---------------------------------------------


def test_tank_pressure_exceeds_chamber_pressure_while_flowing(tank_state):
    for area in (1.0e-6, 3.5e-6, 8.0e-6):
        solution = solve_feed_coupling(
            tank_state, Injector(area), GRAIN, LAW, COMBUSTION, NOZZLE, R0
        )
        assert solution.is_flowing
        assert solution.tank_pressure_pa > solution.chamber_pressure_pa
        assert solution.pressure_margin_pa > 0.0
        assert solution.oxidizer_mass_flow_kg_s > 0.0


def test_no_liquid_gives_no_flow_and_a_dedicated_status(tank_state):
    empty = type(tank_state)(**{**tank_state.__dict__, "liquid_mass_kg": 0.0})
    solution = solve_feed_coupling(
        empty, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0
    )
    assert solution.status is FeedStatus.LIQUID_DEPLETED
    assert solution.oxidizer_mass_flow_kg_s == 0.0
    assert solution.fuel_mass_flow_kg_s == 0.0
    assert solution.chamber_pressure_pa == 0.0
    assert math.isnan(solution.mixture_ratio)
    assert solution.injector is None
    assert not solution.is_flowing


# --- 31, 32. parametric direction ----------------------------------------


def test_larger_effective_area_raises_the_coupled_flow(tank_state):
    flows = [
        solve_feed_coupling(
            tank_state, Injector(area), GRAIN, LAW, COMBUSTION, NOZZLE, R0
        ).oxidizer_mass_flow_kg_s
        for area in (2.0e-6, 3.0e-6, 4.0e-6, 5.0e-6)
    ]
    assert all(later > earlier for earlier, later in itertools.pairwise(flows))


def test_larger_effective_area_raises_the_chamber_pressure(tank_state):
    pressures = [
        solve_feed_coupling(
            tank_state, Injector(area), GRAIN, LAW, COMBUSTION, NOZZLE, R0
        ).chamber_pressure_pa
        for area in (2.0e-6, 3.0e-6, 4.0e-6, 5.0e-6)
    ]
    assert all(later > earlier for earlier, later in itertools.pairwise(pressures))


def test_higher_tank_pressure_raises_the_coupled_flow():
    flows = []
    for temperature in (273.15, 283.15, 293.15, 303.15):
        state = REFERENCE_TANK.initial_state(temperature, 0.8)
        flows.append(
            solve_feed_coupling(
                state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0
            ).oxidizer_mass_flow_kg_s
        )
    assert all(later > earlier for earlier, later in itertools.pairwise(flows))


def test_a_wider_port_raises_the_chamber_pressure_and_backs_off_the_injector(tank_state):
    """The Milestone 2 exponent effect, carried through the coupling.

    A wider port lowers ``G_ox`` and hence ``r_dot``, but the burning area grows
    faster (the net exponent on the radius is ``1 - 2n = +0.267``), so the fuel
    flow and the total flow both rise and the chamber stiffens.  The injector
    then sees a smaller pressure drop and delivers slightly *less* oxidizer --
    a feedback that has no counterpart in the prescribed-flow Milestone 3 model.
    """
    solutions = [
        solve_feed_coupling(
            tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, radius
        )
        for radius in (0.020, 0.030, 0.040)
    ]
    pressures = [s.chamber_pressure_pa for s in solutions]
    fuel = [s.fuel_mass_flow_kg_s for s in solutions]
    oxidizer = [s.oxidizer_mass_flow_kg_s for s in solutions]
    assert all(later > earlier for earlier, later in itertools.pairwise(pressures))
    assert all(later > earlier for earlier, later in itertools.pairwise(fuel))
    assert all(later < earlier for earlier, later in itertools.pairwise(oxidizer))


def test_the_spi_model_predicts_more_flow_than_the_dyer_blend(tank_state):
    spi = solve_feed_coupling(
        tank_state, Injector(AREA, model=InjectorModel.SPI), GRAIN, LAW, COMBUSTION, NOZZLE, R0
    )
    dyer = solve_feed_coupling(
        tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0
    )
    assert spi.oxidizer_mass_flow_kg_s > dyer.oxidizer_mass_flow_kg_s
    assert spi.chamber_pressure_pa > dyer.chamber_pressure_pa


# --- 29, 30. determinism and tolerance -----------------------------------


def test_the_solve_is_deterministic(tank_state):
    first = solve_feed_coupling(tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0)
    second = solve_feed_coupling(tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0)
    assert first.oxidizer_mass_flow_kg_s == second.oxidizer_mass_flow_kg_s
    assert first.chamber_pressure_pa == second.chamber_pressure_pa


def test_the_solve_converges_as_the_root_tolerance_tightens(tank_state):
    reference = solve_feed_coupling(
        tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0, xtol=1e-14
    ).oxidizer_mass_flow_kg_s
    for xtol in (1e-6, 1e-9, 1e-12):
        solution = solve_feed_coupling(
            tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0, xtol=xtol
        )
        assert solution.oxidizer_mass_flow_kg_s == pytest.approx(reference, abs=max(xtol, 1e-12))


# --- input validation -----------------------------------------------------


@pytest.mark.parametrize("bad_radius", [0.0, -0.01, math.nan, math.inf])
def test_invalid_port_radius_is_rejected(bad_radius, tank_state):
    with pytest.raises(ValueError):
        solve_feed_coupling(
            tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, bad_radius
        )


def test_a_burnt_through_port_is_rejected_by_default(tank_state):
    with pytest.raises(ValueError, match="burnt through"):
        solve_feed_coupling(
            tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, 0.050
        )


def test_the_burnout_overshoot_allowance_is_bounded(tank_state):
    """A small allowance is accepted for event bracketing; a large one is not."""
    solution = solve_feed_coupling(
        tank_state,
        Injector(AREA),
        GRAIN,
        LAW,
        COMBUSTION,
        NOZZLE,
        0.0455,
        burnout_overshoot_fraction=0.05,
    )
    assert solution.is_flowing
    with pytest.raises(ValueError, match="burnt through"):
        solve_feed_coupling(
            tank_state,
            Injector(AREA),
            GRAIN,
            LAW,
            COMBUSTION,
            NOZZLE,
            0.060,
            burnout_overshoot_fraction=0.05,
        )


@pytest.mark.parametrize("bad", [-0.1, math.nan, math.inf])
def test_invalid_overshoot_fraction_is_rejected(bad, tank_state):
    with pytest.raises(ValueError):
        solve_feed_coupling(
            tank_state,
            Injector(AREA),
            GRAIN,
            LAW,
            COMBUSTION,
            NOZZLE,
            R0,
            burnout_overshoot_fraction=bad,
        )


def test_solution_is_immutable(tank_state):
    solution = solve_feed_coupling(tank_state, Injector(AREA), GRAIN, LAW, COMBUSTION, NOZZLE, R0)
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        solution.oxidizer_mass_flow_kg_s = 1.0
