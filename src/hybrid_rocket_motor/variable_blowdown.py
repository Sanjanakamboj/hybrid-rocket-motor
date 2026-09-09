"""Coupled blowdown with O/F-dependent thermochemistry (Milestone 5 top level).

This is the Milestone 4 model with one assumption removed.  Milestone 4 carried a
single prescribed ``(c*, gamma, T_c, M)`` set through the whole burn; here those
four numbers are interpolated from a frozen NASA CEA table at the *instantaneous*
mixture ratio and chamber pressure, and fed back into the same coupled loop.

Integrated state -- identical to Milestone 4
--------------------------------------------
::

    y = [ r_p , m_tank , U_tank , m_ox_cumulative , m_f_cumulative , impulse ]

    d r_p / dt    = r_dot
    d m_tank/dt   = -m_dot_ox
    d U_tank/dt   = -m_dot_ox h_l(T_tank)
    d m_ox_cum/dt = m_dot_ox
    d m_f_cum/dt  = m_dot_f
    d impulse/dt  = F

Nothing about the tank, the injector, the regression law or the nozzle relations
changes.  What changes is that each right-hand-side evaluation now solves a
*nested* closure (:mod:`hybrid_rocket_motor.variable_feed_system`) instead of a
single feed root, because ``c*`` depends on the operating point it helps set.

Termination
-----------
The four Milestone 4 terminal events are kept unchanged -- grain burnout, liquid
depletion, tank/chamber pressure equalisation, tank temperature leaving the
checked property band.  Milestone 5 adds four more, one per edge of the
thermochemical table:

* mixture ratio below / above the tabulated O/F span,
* chamber pressure below / above the tabulated pressure span.

The table is **never extrapolated**.  If the coupled state walks off the table
the run stops and reports which edge it left by, rather than quietly continuing
on an invented property.

How the table events are bracketed
----------------------------------
Each event returns the signed distance to its table edge while the coupled solve
succeeds, so it decreases smoothly to zero as the state approaches that edge.
Past the edge the solve has no valid answer at all and the event returns a fixed
negative value.  The jump sits exactly at the edge -- the last point where the
solve succeeds is the last point where the margin is non-negative -- so the root
located by the integrator is the crossing itself, not an artefact of the step
size.  These events are guards: the reference cases in the Milestone 5 study do
not approach any table edge, and a run that terminates this way says so in its
:class:`VariableBlowdownTermination`.

Units
-----
SI throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

import numpy as np
from scipy.integrate import solve_ivp

from .blowdown import BURNOUT_OVERSHOOT_FRACTION, DEPLETION_QUALITY_OVERSHOOT
from .geometry import GrainGeometry
from .injector import Injector
from .nitrous_properties import MIN_VERIFIED_TEMPERATURE_K
from .nozzle import NozzleGeometry
from .performance import STANDARD_GRAVITY_M_S2
from .regression import PowerLawRegressionLaw
from .tank import NitrousTank, TankState
from .thermochemistry import DEFAULT_C_STAR_EFFICIENCY, ThermochemistryTable
from .variable_feed_system import (
    VariableFeedSolution,
    VariableFeedStatus,
    solve_variable_feed_coupling,
)
from .variable_nozzle import solve_variable_performance_point

__all__ = [
    "VariableBlowdownResult",
    "VariableBlowdownTermination",
    "simulate_variable_blowdown",
]


class VariableBlowdownTermination(Enum):
    """Why the coupled variable-property integration stopped."""

    GRAIN_BURNOUT = "GRAIN_BURNOUT"
    LIQUID_DEPLETED = "LIQUID_DEPLETED"
    PRESSURE_EQUALIZED = "PRESSURE_EQUALIZED"
    PROPERTY_LIMIT = "PROPERTY_LIMIT"
    OF_BELOW_TABLE = "OF_BELOW_TABLE"
    OF_ABOVE_TABLE = "OF_ABOVE_TABLE"
    CHAMBER_PRESSURE_BELOW_TABLE = "CHAMBER_PRESSURE_BELOW_TABLE"
    CHAMBER_PRESSURE_ABOVE_TABLE = "CHAMBER_PRESSURE_ABOVE_TABLE"
    COMPLETED_TIME_WINDOW = "COMPLETED_TIME_WINDOW"


@dataclass(frozen=True)
class VariableBlowdownResult:
    """Coupled time histories from a variable-property blowdown simulation."""

    time_s: np.ndarray
    port_radius_m: np.ndarray
    port_diameter_m: np.ndarray
    tank_temperature_k: np.ndarray
    tank_pressure_pa: np.ndarray
    tank_total_mass_kg: np.ndarray
    tank_liquid_mass_kg: np.ndarray
    tank_vapour_mass_kg: np.ndarray
    tank_vapour_quality: np.ndarray
    tank_internal_energy_j: np.ndarray
    oxidizer_mass_flow_kg_s: np.ndarray
    fuel_mass_flow_kg_s: np.ndarray
    total_mass_flow_kg_s: np.ndarray
    mixture_ratio: np.ndarray
    regression_rate_m_s: np.ndarray
    oxidizer_mass_flux_si: np.ndarray
    chamber_pressure_pa: np.ndarray
    pressure_margin_pa: np.ndarray
    feed_residual_kg_s: np.ndarray
    chamber_residual_pa: np.ndarray
    # -- thermochemistry, the Milestone 5 additions ------------------------
    ideal_c_star_m_s: np.ndarray
    c_star_m_s: np.ndarray
    gamma: np.ndarray
    chamber_temperature_k: np.ndarray
    effective_molar_mass_kg_mol: np.ndarray
    specific_gas_constant_j_kg_k: np.ndarray
    of_margin_low: np.ndarray
    of_margin_high: np.ndarray
    pressure_margin_low_pa: np.ndarray
    pressure_margin_high_pa: np.ndarray
    # -- nozzle and thrust --------------------------------------------------
    exit_mach: np.ndarray
    exit_pressure_pa: np.ndarray
    exit_velocity_m_s: np.ndarray
    momentum_thrust_n: np.ndarray
    pressure_thrust_n: np.ndarray
    thrust_n: np.ndarray
    thrust_coefficient: np.ndarray
    specific_impulse_s: np.ndarray
    cumulative_oxidizer_mass_kg: np.ndarray
    cumulative_fuel_mass_kg: np.ndarray
    cumulative_impulse_n_s: np.ndarray
    feed_status: tuple[VariableFeedStatus, ...]
    termination: VariableBlowdownTermination
    tank: NitrousTank
    injector: Injector
    geometry: GrainGeometry
    law: PowerLawRegressionLaw
    table: ThermochemistryTable
    nozzle: NozzleGeometry
    ambient_pressure_pa: float
    c_star_efficiency: float
    initial_tank_state: TankState

    # -- summary ------------------------------------------------------------

    @property
    def burn_time_s(self) -> float:
        return float(self.time_s[-1])

    @property
    def total_impulse_n_s(self) -> float:
        return float(self.cumulative_impulse_n_s[-1])

    @property
    def oxidizer_consumed_kg(self) -> float:
        return float(self.cumulative_oxidizer_mass_kg[-1])

    @property
    def fuel_consumed_kg(self) -> float:
        return float(self.cumulative_fuel_mass_kg[-1])

    @property
    def propellant_consumed_kg(self) -> float:
        return self.oxidizer_consumed_kg + self.fuel_consumed_kg

    @property
    def equivalent_specific_impulse_s(self) -> float:
        mass = self.propellant_consumed_kg
        if mass <= 0.0:
            return math.nan
        return self.total_impulse_n_s / (mass * STANDARD_GRAVITY_M_S2)

    @property
    def peak_thrust_n(self) -> float:
        return float(np.max(self.thrust_n))

    @property
    def mean_thrust_n(self) -> float:
        return self.total_impulse_n_s / self.burn_time_s if self.burn_time_s > 0.0 else 0.0

    @property
    def firing(self) -> np.ndarray:
        """Boolean mask of samples at which propellant is actually flowing."""
        return np.array(
            [status is VariableFeedStatus.FLOWING for status in self.feed_status],
            dtype=bool,
        )

    @property
    def stayed_within_table(self) -> bool:
        """Whether every firing sample sat strictly inside the tabulated rectangle."""
        mask = self.firing
        margins = (
            self.of_margin_low[mask],
            self.of_margin_high[mask],
            self.pressure_margin_low_pa[mask],
            self.pressure_margin_high_pa[mask],
        )
        return bool(all(np.all(m >= 0.0) for m in margins))

    # -- independent closure checks ----------------------------------------

    @property
    def oxidizer_mass_closure_residual_kg(self) -> np.ndarray:
        """Oxidizer drawn from the tank minus oxidizer integrated at the injector."""
        drawn = self.tank_total_mass_kg[0] - self.tank_total_mass_kg
        return self.cumulative_oxidizer_mass_kg - drawn

    @property
    def fuel_mass_closure_residual_kg(self) -> np.ndarray:
        """Integrated fuel flow minus the fuel implied by the port geometry [kg]."""
        geometric = (
            self.geometry.fuel_density_kg_m3
            * self.geometry.length_m
            * math.pi
            * (self.port_radius_m**2 - self.port_radius_m[0] ** 2)
        )
        return self.cumulative_fuel_mass_kg - geometric

    @property
    def c_star_identity_residual_pa(self) -> np.ndarray:
        """``m_dot c* / A_t - p_c`` at every firing sample [Pa].

        Zero by construction at a converged state, so this measures how well the
        *reported* histories satisfy the chamber closure rather than re-deriving
        it.  Non-firing samples are reported as zero.
        """
        mask = self.firing
        residual = np.zeros_like(self.chamber_pressure_pa)
        residual[mask] = (
            self.total_mass_flow_kg_s[mask] * self.c_star_m_s[mask] / self.nozzle.throat_area_m2
            - self.chamber_pressure_pa[mask]
        )
        return residual

    @property
    def mixture_ratio_identity_residual(self) -> np.ndarray:
        """``O/F - m_dot_ox / m_dot_f`` at every firing sample [-]."""
        mask = self.firing
        residual = np.zeros_like(self.mixture_ratio)
        residual[mask] = (
            self.mixture_ratio[mask]
            - self.oxidizer_mass_flow_kg_s[mask] / self.fuel_mass_flow_kg_s[mask]
        )
        return residual

    @property
    def relative_impulse_closure_residual(self) -> float:
        """Integrated impulse state against a trapezoidal integral of thrust [-]."""
        total = self.total_impulse_n_s
        if total <= 0.0:
            return 0.0
        trapezoid = float(np.trapezoid(self.thrust_n, self.time_s))
        return abs(total - trapezoid) / total


def simulate_variable_blowdown(
    tank: NitrousTank,
    initial_tank_state: TankState,
    injector: Injector,
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    table: ThermochemistryTable,
    nozzle: NozzleGeometry,
    ambient_pressure_pa: float,
    *,
    c_star_efficiency: float = DEFAULT_C_STAR_EFFICIENCY,
    initial_port_diameter_m: float | None = None,
    t_end_s: float = 200.0,
    n_report: int = 201,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
    method: str = "LSODA",
    feed_xtol: float = 1.0e-12,
    pressure_margin_floor_pa: float = 1.0e4,
) -> VariableBlowdownResult:
    """Integrate the coupled tank/injector/motor system with tabulated chemistry.

    The signature mirrors :func:`hybrid_rocket_motor.blowdown.simulate_blowdown`
    with the prescribed ``combustion`` argument replaced by ``table`` plus the
    single ``c_star_efficiency`` that turns tabulated ideal ``c*`` into delivered
    ``c*``.  Every numerical control keeps its Milestone 4 default so the two
    models can be compared without a tolerance difference confounding the result.
    """
    if n_report < 2:
        raise ValueError(f"n_report must be at least 2, got {n_report!r}")
    p_a = float(ambient_pressure_pa)
    if not math.isfinite(p_a) or p_a < 0.0:
        raise ValueError(f"ambient_pressure_pa must be finite and non-negative, got {p_a!r}")

    diameter0 = (
        geometry.initial_port_diameter_m
        if initial_port_diameter_m is None
        else float(initial_port_diameter_m)
    )
    geometry.port_area_m2(diameter0)
    r0 = diameter0 / 2.0
    r_outer = geometry.outer_diameter_m / 2.0

    of_low, of_high = table.of_bounds
    pressure_low, pressure_high = table.pressure_bounds_pa

    cache: dict[tuple[float, float, float], tuple[TankState, VariableFeedSolution]] = {}

    def evaluate(
        radius: float, mass: float, energy: float
    ) -> tuple[TankState, VariableFeedSolution]:
        key = (radius, mass, energy)
        hit = cache.get(key)
        if hit is not None:
            return hit
        tank_state = tank.state_from_mass_and_energy(
            mass, energy, quality_overshoot=DEPLETION_QUALITY_OVERSHOOT
        )
        feed = solve_variable_feed_coupling(
            tank_state,
            injector,
            geometry,
            law,
            table,
            nozzle,
            radius,
            c_star_efficiency=c_star_efficiency,
            xtol=feed_xtol,
            burnout_overshoot_fraction=BURNOUT_OVERSHOOT_FRACTION,
        )
        if len(cache) > 4096:
            cache.clear()
        cache[key] = (tank_state, feed)
        return tank_state, feed

    def rhs(t: float, y: np.ndarray) -> list[float]:
        radius, mass, energy = float(y[0]), float(y[1]), float(y[2])
        tank_state, feed = evaluate(radius, mass, energy)
        m_ox = feed.oxidizer_mass_flow_kg_s
        thrust = 0.0
        if feed.is_flowing:
            thrust = solve_variable_performance_point(
                feed.thermochemistry, nozzle, m_ox, feed.fuel_mass_flow_kg_s, p_a
            ).thrust_n
        return [
            feed.regression_rate_m_s,
            -m_ox,
            -m_ox * tank_state.liquid_enthalpy_j_kg,
            m_ox,
            feed.fuel_mass_flow_kg_s,
            thrust,
        ]

    # -- Milestone 4 terminal events, unchanged ------------------------------

    def event_burnout(t: float, y: np.ndarray) -> float:
        return float(y[0]) - r_outer

    def event_liquid(t: float, y: np.ndarray) -> float:
        tank_state, _ = evaluate(float(y[0]), float(y[1]), float(y[2]))
        return tank_state.liquid_fraction_remaining

    def event_margin(t: float, y: np.ndarray) -> float:
        _, feed = evaluate(float(y[0]), float(y[1]), float(y[2]))
        return feed.pressure_margin_pa - pressure_margin_floor_pa

    def event_property(t: float, y: np.ndarray) -> float:
        tank_state, _ = evaluate(float(y[0]), float(y[1]), float(y[2]))
        return tank_state.temperature_k - MIN_VERIFIED_TEMPERATURE_K

    # -- Milestone 5 table-validity events -----------------------------------

    def table_margin(feed: VariableFeedSolution, which: str) -> float:
        """Signed distance to one table edge; see the module docstring."""
        if feed.is_flowing:
            state = feed.thermochemistry
            if which == "of_low":
                return state.of_ratio - of_low
            if which == "of_high":
                return of_high - state.of_ratio
            if which == "p_low":
                return state.chamber_pressure_pa - pressure_low
            return pressure_high - state.chamber_pressure_pa
        scale = of_low if which.startswith("of") else pressure_low
        failed_on_of = feed.status is VariableFeedStatus.OF_OUTSIDE_TABLE
        failed_on_pressure = feed.status is VariableFeedStatus.PRESSURE_OUTSIDE_TABLE
        # Only the edge that actually failed is driven negative; the others stay
        # positive so a liquid-depletion or equalisation stop is not mislabelled.
        if which == "of_low" and failed_on_of:
            return -scale
        if which == "p_low" and failed_on_pressure:
            return -scale
        return scale

    def make_table_event(which: str):
        def event(t: float, y: np.ndarray) -> float:
            _, feed = evaluate(float(y[0]), float(y[1]), float(y[2]))
            return table_margin(feed, which)

        return event

    event_of_low = make_table_event("of_low")
    event_of_high = make_table_event("of_high")
    event_pressure_low = make_table_event("p_low")
    event_pressure_high = make_table_event("p_high")

    events = [
        event_burnout,
        event_liquid,
        event_margin,
        event_property,
        event_of_low,
        event_of_high,
        event_pressure_low,
        event_pressure_high,
    ]
    directions = (1, -1, -1, -1, -1, -1, -1, -1)
    reasons = (
        VariableBlowdownTermination.GRAIN_BURNOUT,
        VariableBlowdownTermination.LIQUID_DEPLETED,
        VariableBlowdownTermination.PRESSURE_EQUALIZED,
        VariableBlowdownTermination.PROPERTY_LIMIT,
        VariableBlowdownTermination.OF_BELOW_TABLE,
        VariableBlowdownTermination.OF_ABOVE_TABLE,
        VariableBlowdownTermination.CHAMBER_PRESSURE_BELOW_TABLE,
        VariableBlowdownTermination.CHAMBER_PRESSURE_ABOVE_TABLE,
    )
    for event, direction in zip(events, directions, strict=True):
        event.terminal = True
        event.direction = direction

    y0 = np.array(
        [
            r0,
            initial_tank_state.total_mass_kg,
            initial_tank_state.total_internal_energy_j,
            0.0,
            0.0,
            0.0,
        ]
    )

    solution = solve_ivp(
        rhs,
        (0.0, float(t_end_s)),
        y0,
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=True,
        events=events,
    )
    if not solution.success:
        raise RuntimeError(f"variable blowdown integration failed: {solution.message}")

    termination = VariableBlowdownTermination.COMPLETED_TIME_WINDOW
    t_final = float(solution.t[-1])
    # Earliest event wins, not the first in list order: the table guards and the
    # physical stops are genuinely independent and either may come first.
    best: float | None = None
    for reason, times in zip(reasons, solution.t_events, strict=True):
        if times.size > 0 and (best is None or float(times[0]) < best):
            best = float(times[0])
            termination = reason
    if best is not None:
        t_final = best

    times = np.linspace(0.0, t_final, n_report)
    rows = [solution.sol(float(t)) for t in times]

    tank_states: list[TankState] = []
    feeds: list[VariableFeedSolution] = []
    for row in rows:
        tank_state, feed = evaluate(float(row[0]), float(row[1]), float(row[2]))
        tank_states.append(tank_state)
        feeds.append(feed)

    def thrust_fields(feed: VariableFeedSolution):
        if not feed.is_flowing:
            return (math.nan, math.nan, math.nan, 0.0, 0.0, 0.0, math.nan, math.nan)
        point = solve_variable_performance_point(
            feed.thermochemistry,
            nozzle,
            feed.oxidizer_mass_flow_kg_s,
            feed.fuel_mass_flow_kg_s,
            p_a,
        )
        exit_state = point.exit_state
        return (
            exit_state.mach,
            exit_state.pressure_pa,
            exit_state.velocity_m_s,
            point.momentum_thrust_n,
            point.pressure_thrust_n,
            point.thrust_n,
            point.thrust_coefficient,
            point.specific_impulse_s,
        )

    thrust_rows = np.array([thrust_fields(f) for f in feeds], dtype=float)

    def column(getter, source) -> np.ndarray:
        return np.array([getter(item) for item in source], dtype=float)

    def chemistry(getter) -> np.ndarray:
        return np.array(
            [getter(f.thermochemistry) if f.is_flowing else math.nan for f in feeds],
            dtype=float,
        )

    radii = np.array([row[0] for row in rows], dtype=float)
    if termination is VariableBlowdownTermination.GRAIN_BURNOUT:
        radii = np.minimum(radii, r_outer)

    return VariableBlowdownResult(
        time_s=times,
        port_radius_m=radii,
        port_diameter_m=2.0 * radii,
        tank_temperature_k=column(lambda s: s.temperature_k, tank_states),
        tank_pressure_pa=column(lambda s: s.pressure_pa, tank_states),
        tank_total_mass_kg=column(lambda s: s.total_mass_kg, tank_states),
        tank_liquid_mass_kg=column(lambda s: s.liquid_mass_kg, tank_states),
        tank_vapour_mass_kg=column(lambda s: s.vapour_mass_kg, tank_states),
        tank_vapour_quality=column(lambda s: s.vapour_quality, tank_states),
        tank_internal_energy_j=np.array([row[2] for row in rows], dtype=float),
        oxidizer_mass_flow_kg_s=column(lambda f: f.oxidizer_mass_flow_kg_s, feeds),
        fuel_mass_flow_kg_s=column(lambda f: f.fuel_mass_flow_kg_s, feeds),
        total_mass_flow_kg_s=column(lambda f: f.total_mass_flow_kg_s, feeds),
        mixture_ratio=column(lambda f: f.mixture_ratio, feeds),
        regression_rate_m_s=column(lambda f: f.regression_rate_m_s, feeds),
        oxidizer_mass_flux_si=column(lambda f: f.oxidizer_mass_flux_si, feeds),
        chamber_pressure_pa=column(lambda f: f.chamber_pressure_pa, feeds),
        pressure_margin_pa=column(lambda f: f.pressure_margin_pa, feeds),
        feed_residual_kg_s=column(lambda f: f.feed_residual_kg_s, feeds),
        chamber_residual_pa=column(lambda f: f.chamber_residual_pa, feeds),
        ideal_c_star_m_s=chemistry(lambda s: s.ideal_c_star_m_s),
        c_star_m_s=chemistry(lambda s: s.c_star_m_s),
        gamma=chemistry(lambda s: s.gamma),
        chamber_temperature_k=chemistry(lambda s: s.chamber_temperature_k),
        effective_molar_mass_kg_mol=chemistry(lambda s: s.effective_molar_mass_kg_mol),
        specific_gas_constant_j_kg_k=chemistry(lambda s: s.specific_gas_constant_j_kg_k),
        of_margin_low=column(lambda f: table_margin(f, "of_low"), feeds),
        of_margin_high=column(lambda f: table_margin(f, "of_high"), feeds),
        pressure_margin_low_pa=column(lambda f: table_margin(f, "p_low"), feeds),
        pressure_margin_high_pa=column(lambda f: table_margin(f, "p_high"), feeds),
        exit_mach=thrust_rows[:, 0],
        exit_pressure_pa=thrust_rows[:, 1],
        exit_velocity_m_s=thrust_rows[:, 2],
        momentum_thrust_n=thrust_rows[:, 3],
        pressure_thrust_n=thrust_rows[:, 4],
        thrust_n=thrust_rows[:, 5],
        thrust_coefficient=thrust_rows[:, 6],
        specific_impulse_s=thrust_rows[:, 7],
        cumulative_oxidizer_mass_kg=np.array([row[3] for row in rows], dtype=float),
        cumulative_fuel_mass_kg=np.array([row[4] for row in rows], dtype=float),
        cumulative_impulse_n_s=np.array([row[5] for row in rows], dtype=float),
        feed_status=tuple(f.status for f in feeds),
        termination=termination,
        tank=tank,
        injector=injector,
        geometry=geometry,
        law=law,
        table=table,
        nozzle=nozzle,
        ambient_pressure_pa=p_a,
        c_star_efficiency=float(c_star_efficiency),
        initial_tank_state=initial_tank_state,
    )
