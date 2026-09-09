"""Coupled N2O blowdown: tank + injector + regression + chamber + nozzle.

This is the Milestone 4 top-level model.  It removes the prescribed oxidizer mass
flow of Milestones 1-3 and replaces it with a tank/injector closure, so the
oxidizer flow, chamber pressure and thrust are all outputs of one coupled system.

Integrated state
----------------
::

    y = [ r_p , m_tank , U_tank , m_ox_cumulative , m_f_cumulative , impulse ]

with, at every evaluation:

1. the tank state recovered from ``(m_tank, U_tank)`` by the saturated
   equilibrium flash (:mod:`hybrid_rocket_motor.tank`),
2. the coupled oxidizer flow solved from the injector/chamber/regression loop
   (:mod:`hybrid_rocket_motor.feed_system`),
3. the *frozen* Milestone 1 regression law giving ``r_dot`` and ``m_dot_f``,
4. the *frozen* Milestone 3 chamber and nozzle model giving ``p_c`` and thrust.

The derivatives are then::

    d r_p / dt   = r_dot
    d m_tank/dt  = -m_dot_ox
    d U_tank/dt  = -m_dot_ox h_l(T_tank)
    d m_ox_cum/dt = m_dot_ox
    d m_f_cum/dt  = m_dot_f
    d impulse/dt  = F

Nothing in Milestones 1-3 is modified: this module composes them.

Termination
-----------
Integration stops on the first physically meaningful event -- grain burnout,
liquid depletion, tank/chamber pressure equalisation, or the tank temperature
leaving the checked property band -- rather than integrating blindly into a state
where the model does not hold.  The vapour-only discharge tail after liquid
depletion is **not modelled**; the run terminates there and says so.

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

from .chamber import CombustionProperties
from .feed_system import FeedSolution, FeedStatus, solve_feed_coupling
from .geometry import GrainGeometry
from .injector import Injector
from .nitrous_properties import MIN_VERIFIED_TEMPERATURE_K
from .nozzle import NozzleGeometry
from .performance import STANDARD_GRAVITY_M_S2, solve_performance_point
from .regression import PowerLawRegressionLaw
from .tank import NitrousTank, TankState

__all__ = [
    "BURNOUT_OVERSHOOT_FRACTION",
    "DEPLETION_QUALITY_OVERSHOOT",
    "BlowdownResult",
    "BlowdownTermination",
    "simulate_blowdown",
]


#: How far past the outer grain radius the integrator is permitted to step while
#: bracketing the terminal burnout event.  Purely numerical: the event root is
#: located inside this margin and every reported radius is clamped back to the
#: outer radius, so no overshoot appears in any result.
BURNOUT_OVERSHOOT_FRACTION = 0.05

#: How far past liquid depletion the tank flash is allowed to be evaluated while
#: the terminal depletion event is bracketed, expressed as vapour quality beyond
#: 1.  Purely numerical, exactly like :data:`BURNOUT_OVERSHOOT_FRACTION`: the
#: reported phase masses stay clamped at zero and no sample past the located
#: event is ever reported.
DEPLETION_QUALITY_OVERSHOOT = 0.6


class BlowdownTermination(Enum):
    """Why the coupled integration stopped."""

    GRAIN_BURNOUT = "GRAIN_BURNOUT"
    LIQUID_DEPLETED = "LIQUID_DEPLETED"
    PRESSURE_EQUALIZED = "PRESSURE_EQUALIZED"
    PROPERTY_LIMIT = "PROPERTY_LIMIT"
    COMPLETED_TIME_WINDOW = "COMPLETED_TIME_WINDOW"


@dataclass(frozen=True)
class BlowdownResult:
    """Coupled time histories from a blowdown simulation."""

    time_s: np.ndarray
    port_radius_m: np.ndarray
    port_diameter_m: np.ndarray
    tank_temperature_k: np.ndarray
    tank_pressure_pa: np.ndarray
    tank_total_mass_kg: np.ndarray
    tank_liquid_mass_kg: np.ndarray
    tank_vapour_mass_kg: np.ndarray
    tank_liquid_volume_m3: np.ndarray
    tank_vapour_volume_m3: np.ndarray
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
    spi_mass_flow_kg_s: np.ndarray
    hem_mass_flow_kg_s: np.ndarray
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
    feed_status: tuple[FeedStatus, ...]
    termination: BlowdownTermination
    tank: NitrousTank
    injector: Injector
    geometry: GrainGeometry
    law: PowerLawRegressionLaw
    combustion: CombustionProperties
    nozzle: NozzleGeometry
    ambient_pressure_pa: float
    initial_tank_state: TankState

    # -- integrated performance -------------------------------------------

    @property
    def burn_time_s(self) -> float:
        return float(self.time_s[-1])

    @property
    def total_impulse_n_s(self) -> float:
        """Total impulse from the integrated impulse state [N s]."""
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
        """``I_total / (m_propellant g_0)`` [s]."""
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
    def minimum_pressure_margin_pa(self) -> float:
        return float(np.min(self.pressure_margin_pa))

    # -- independent closure checks ---------------------------------------

    @property
    def firing(self) -> np.ndarray:
        """Boolean mask of samples at which propellant is actually flowing."""
        return np.array(
            [status is FeedStatus.FLOWING for status in self.feed_status], dtype=bool
        )

    @property
    def oxidizer_mass_closure_residual_kg(self) -> np.ndarray:
        """Oxidizer drawn from the tank minus oxidizer integrated at the injector.

        Zero in exact arithmetic (the tank mass and the cumulative-oxidizer state
        are driven by the same derivative with opposite signs), so this is a pure
        measure of integration error.
        """
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
    def impulse_closure_residual_n_s(self) -> float:
        """Integrated impulse state minus a trapezoidal integral of thrust [N s].

        A coarse reconstruction check.  When a run ends on liquid depletion the
        thrust steps discontinuously to zero at the terminal sample, and a
        trapezoid cannot represent that step; the residual is then dominated by
        the discontinuity rather than by integration error.  The impulse state
        carried by the integrator is the accurate quantity.  Use
        :attr:`relative_impulse_closure_residual` to judge the size.
        """
        return self.total_impulse_n_s - float(np.trapezoid(self.thrust_n, self.time_s))

    @property
    def relative_impulse_closure_residual(self) -> float:
        """:attr:`impulse_closure_residual_n_s` normalised by the total impulse [-]."""
        total = self.total_impulse_n_s
        return abs(self.impulse_closure_residual_n_s) / total if total > 0.0 else 0.0

    @property
    def energy_closure_residual_j(self) -> np.ndarray:
        """Tank internal energy minus its integrated balance [J].

        ``U(t) = U_0 - integral m_dot_ox h_l dt``.  Recomputed here by trapezoidal
        integration of the reported ``m_dot_ox h_l`` product, so it is an
        independent check on the energy state carried by the integrator.

        Evaluated over the **firing** samples only.  A run terminated by liquid
        depletion ends with one non-firing sample at which the outflow drops
        discontinuously to zero; a trapezoid across that step misrepresents the
        integral by orders of magnitude more than the true integration error, so
        including it would measure the discontinuity rather than the accuracy.
        """
        mask = self.firing
        times = self.time_s[mask]
        energies = self.tank_internal_energy_j[mask]
        outflow_enthalpy_rate = self.oxidizer_mass_flow_kg_s[mask] * np.array(
            [
                self.tank.state_from_mass_and_energy(
                    m, u, quality_overshoot=DEPLETION_QUALITY_OVERSHOOT
                ).liquid_enthalpy_j_kg
                for m, u in zip(
                    self.tank_total_mass_kg[mask], energies, strict=True
                )
            ]
        )
        integrated = np.concatenate(
            ([0.0], np.cumsum(
                0.5
                * (outflow_enthalpy_rate[1:] + outflow_enthalpy_rate[:-1])
                * np.diff(times)
            ))
        )
        return energies - (energies[0] - integrated)

    @property
    def relative_energy_closure_residual(self) -> float:
        """Worst tank energy residual normalised by the initial internal energy [-]."""
        return float(np.max(np.abs(self.energy_closure_residual_j))) / abs(
            self.tank_internal_energy_j[0]
        )


def simulate_blowdown(
    tank: NitrousTank,
    initial_tank_state: TankState,
    injector: Injector,
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    combustion: CombustionProperties,
    nozzle: NozzleGeometry,
    ambient_pressure_pa: float,
    *,
    initial_port_diameter_m: float | None = None,
    t_end_s: float = 200.0,
    n_report: int = 201,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
    method: str = "LSODA",
    feed_xtol: float = 1.0e-12,
    pressure_margin_floor_pa: float = 1.0e4,
) -> BlowdownResult:
    """Integrate the coupled tank/injector/motor system.

    Parameters
    ----------
    pressure_margin_floor_pa:
        Terminal event threshold on ``p_tank - p_c``.  The coupled solve loses
        meaning as the margin vanishes, so the run stops while the solution is
        still well posed rather than chasing the singular point.
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

    cache: dict[tuple[float, float, float], tuple[TankState, FeedSolution]] = {}

    def evaluate(radius: float, mass: float, energy: float) -> tuple[TankState, FeedSolution]:
        key = (radius, mass, energy)
        hit = cache.get(key)
        if hit is not None:
            return hit
        tank_state = tank.state_from_mass_and_energy(
            mass, energy, quality_overshoot=DEPLETION_QUALITY_OVERSHOOT
        )
        feed = solve_feed_coupling(
            tank_state,
            injector,
            geometry,
            law,
            combustion,
            nozzle,
            radius,
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
            thrust = solve_performance_point(
                combustion, nozzle, m_ox, feed.fuel_mass_flow_kg_s, p_a
            ).thrust_n
        return [
            feed.regression_rate_m_s,
            -m_ox,
            -m_ox * tank_state.liquid_enthalpy_j_kg,
            m_ox,
            feed.fuel_mass_flow_kg_s,
            thrust,
        ]

    def event_burnout(t: float, y: np.ndarray) -> float:
        return float(y[0]) - r_outer

    def event_liquid(t: float, y: np.ndarray) -> float:
        # The *unclamped* liquid fraction is used so the root can be bracketed:
        # the reported liquid mass is clamped at zero and would go flat instead
        # of changing sign.
        tank_state, _ = evaluate(float(y[0]), float(y[1]), float(y[2]))
        return tank_state.liquid_fraction_remaining

    def event_margin(t: float, y: np.ndarray) -> float:
        _, feed = evaluate(float(y[0]), float(y[1]), float(y[2]))
        return feed.pressure_margin_pa - pressure_margin_floor_pa

    def event_property(t: float, y: np.ndarray) -> float:
        tank_state, _ = evaluate(float(y[0]), float(y[1]), float(y[2]))
        return tank_state.temperature_k - MIN_VERIFIED_TEMPERATURE_K

    for event, direction in (
        (event_burnout, 1),
        (event_liquid, -1),
        (event_margin, -1),
        (event_property, -1),
    ):
        event.terminal = True
        event.direction = direction

    events = [event_burnout, event_liquid, event_margin, event_property]
    reasons = (
        BlowdownTermination.GRAIN_BURNOUT,
        BlowdownTermination.LIQUID_DEPLETED,
        BlowdownTermination.PRESSURE_EQUALIZED,
        BlowdownTermination.PROPERTY_LIMIT,
    )

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
        raise RuntimeError(f"blowdown integration failed: {solution.message}")

    termination = BlowdownTermination.COMPLETED_TIME_WINDOW
    t_final = float(solution.t[-1])
    for reason, times in zip(reasons, solution.t_events, strict=True):
        if times.size > 0:
            termination = reason
            t_final = float(times[0])
            break

    times = np.linspace(0.0, t_final, n_report)
    rows = [solution.sol(float(t)) for t in times]

    tank_states: list[TankState] = []
    feeds: list[FeedSolution] = []
    for row in rows:
        tank_state, feed = evaluate(float(row[0]), float(row[1]), float(row[2]))
        tank_states.append(tank_state)
        feeds.append(feed)

    def thrust_fields(feed: FeedSolution):
        if not feed.is_flowing:
            return (math.nan, math.nan, math.nan, 0.0, 0.0, 0.0, math.nan, math.nan)
        point = solve_performance_point(
            combustion, nozzle, feed.oxidizer_mass_flow_kg_s, feed.fuel_mass_flow_kg_s, p_a
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

    radii = np.array([row[0] for row in rows], dtype=float)
    if termination is BlowdownTermination.GRAIN_BURNOUT:
        radii = np.minimum(radii, r_outer)

    return BlowdownResult(
        time_s=times,
        port_radius_m=radii,
        port_diameter_m=2.0 * radii,
        tank_temperature_k=column(lambda s: s.temperature_k, tank_states),
        tank_pressure_pa=column(lambda s: s.pressure_pa, tank_states),
        tank_total_mass_kg=column(lambda s: s.total_mass_kg, tank_states),
        tank_liquid_mass_kg=column(lambda s: s.liquid_mass_kg, tank_states),
        tank_vapour_mass_kg=column(lambda s: s.vapour_mass_kg, tank_states),
        tank_liquid_volume_m3=column(lambda s: s.liquid_volume_m3, tank_states),
        tank_vapour_volume_m3=column(lambda s: s.vapour_volume_m3, tank_states),
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
        feed_residual_kg_s=column(lambda f: f.residual_kg_s, feeds),
        spi_mass_flow_kg_s=column(
            lambda f: f.injector.spi_mass_flow_kg_s if f.injector else math.nan, feeds
        ),
        hem_mass_flow_kg_s=column(
            lambda f: f.injector.hem_mass_flow_kg_s if f.injector else math.nan, feeds
        ),
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
        combustion=combustion,
        nozzle=nozzle,
        ambient_pressure_pa=p_a,
        initial_tank_state=initial_tank_state,
    )
