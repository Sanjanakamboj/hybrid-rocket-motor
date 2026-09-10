"""Coupled feed / chamber closure with O/F-dependent combustion properties.

Milestone 4 solved one scalar root: the oxidizer flow that makes the injector and
the chamber agree, with a constant ``c*``.  Milestone 5 adds a second implicit
relation, because ``c*`` now depends on the operating point:

    m_dot_ox = injector(p_tank, p_c)              feed relation
    m_dot_f  = rho_f (2 pi r_p L) a G_ox^n        frozen Milestone 1/2 regression
    O/F      = m_dot_ox / m_dot_f
    p_c      = m_dot_total c*(O/F, p_c) / A_t     implicit chamber closure

Architecture
------------
A **nested** solve, both levels bracketed:

* **outer** -- Brent on ``residual(m_dot_ox) = injector_flow(p_c) - m_dot_ox``;
* **inner** -- for each trial oxidizer flow, the fuel flow and hence O/F follow
  directly, and :func:`~hybrid_rocket_motor.variable_chamber.solve_variable_chamber`
  solves the implicit pressure.

No fixed-point iteration is used at either level.

The outer bracket
-----------------
Milestone 4 could bracket on ``[0, m_dot_SPI(p_2 = 0)]``.  That no longer works,
because the chemistry table is valid only on a finite ``(O/F, p_c)`` rectangle,
and both coordinates run out of the table near the ends of that Milestone 4
bracket.  The outer bracket is therefore intersected with the oxidizer-flow
window on which the *inner* chamber solve is guaranteed to have a root inside
the table.  Two constraints define that window.

**Mixture ratio.**  As ``m_dot_ox -> 0`` the mixture ratio tends to zero as well,
because the regression law gives ``m_dot_f ~ m_dot_ox^n`` and hence
``O/F ~ m_dot_ox^(1-n)``.  The table's O/F bounds invert in closed form:

    m_dot_ox(O/F) = [ (O/F) rho_f (2 pi r_p L) a_SI / (pi r_p^2)^n ]^(1/(1-n))

**Chamber pressure.**  Writing the inner residual as ``F(p_c) - p_c`` with
``F(p_c) = m_dot_total c*(O/F, p_c) / A_t``, a root exists inside the tabulated
pressure span exactly when ``F(p_min) >= p_min`` and ``F(p_max) <= p_max``.  Both
tests evaluate the table *on* its boundary, so they need no extrapolation, and
both left-hand sides increase with oxidizer flow (the total mass flow grows
roughly linearly while ``c*`` varies by only a few per cent).  Each therefore
contributes one more bracketed root that trims an end off the window.

If no sign change survives inside the trimmed window the solver reports an
explicit failure status.  The mixture ratio is never clipped into the table, and
the chamber pressure is never clamped to the table edge.

Units
-----
SI throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from scipy.optimize import brentq

from .geometry import GrainGeometry
from .injector import Injector, InjectorResult
from .nitrous_properties import triple_point_pressure_pa
from .nozzle import NozzleGeometry
from .regression import PowerLawRegressionLaw
from .tank import TankState
from .thermochemistry import (
    DEFAULT_C_STAR_EFFICIENCY,
    ThermochemicalState,
    ThermochemistryTable,
)
from .variable_chamber import VariableChamberStatus, solve_variable_chamber

__all__ = [
    "VariableFeedSolution",
    "VariableFeedStatus",
    "oxidizer_flow_for_mixture_ratio",
    "solve_variable_feed_coupling",
]


class VariableFeedStatus(Enum):
    """Outcome of the coupled variable-property feed solve."""

    FLOWING = "FLOWING"
    PRESSURE_EQUALIZED = "PRESSURE_EQUALIZED"
    LIQUID_DEPLETED = "LIQUID_DEPLETED"
    NO_ROOT = "NO_ROOT"
    OF_OUTSIDE_TABLE = "OF_OUTSIDE_TABLE"
    PRESSURE_OUTSIDE_TABLE = "PRESSURE_OUTSIDE_TABLE"


@dataclass(frozen=True)
class VariableFeedSolution:
    """The converged coupled operating point at one instant."""

    oxidizer_mass_flow_kg_s: float
    fuel_mass_flow_kg_s: float
    total_mass_flow_kg_s: float
    mixture_ratio: float
    chamber_pressure_pa: float
    tank_pressure_pa: float
    pressure_margin_pa: float
    regression_rate_m_s: float
    oxidizer_mass_flux_si: float
    thermochemistry: ThermochemicalState | None
    injector: InjectorResult | None
    feed_residual_kg_s: float
    chamber_residual_pa: float
    status: VariableFeedStatus

    @property
    def is_flowing(self) -> bool:
        return self.status is VariableFeedStatus.FLOWING


def _fuel_flow_and_rate(
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    port_radius_m: float,
    oxidizer_mass_flow_kg_s: float,
) -> tuple[float, float, float]:
    """Frozen Milestone 1/2 chain: ``m_dot_ox, r_p -> G_ox, r_dot, m_dot_f``."""
    area = math.pi * port_radius_m * port_radius_m
    flux = oxidizer_mass_flow_kg_s / area
    rate = law.regression_rate_si(flux)
    fuel = geometry.fuel_density_kg_m3 * 2.0 * math.pi * port_radius_m * geometry.length_m * rate
    return flux, rate, fuel


def oxidizer_flow_for_mixture_ratio(
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    port_radius_m: float,
    mixture_ratio: float,
) -> float:
    """Invert ``O/F(m_dot_ox)`` at fixed port radius [kg/s].

    Substituting the regression law into ``O/F = m_dot_ox / m_dot_f`` gives
    ``O/F ~ m_dot_ox^(1-n)``, which inverts exactly:

        m_dot_ox = [ (O/F) rho_f (2 pi r L) a_SI / (pi r^2)^n ]^(1/(1-n))
    """
    of_ratio = float(mixture_ratio)
    if not math.isfinite(of_ratio) or of_ratio <= 0.0:
        raise ValueError(f"mixture_ratio must be finite and positive, got {mixture_ratio!r}")
    radius = float(port_radius_m)
    exponent = law.exponent
    if exponent >= 1.0:
        raise ValueError(f"the inversion requires a regression exponent below 1, got {exponent!r}")
    burning_area = 2.0 * math.pi * radius * geometry.length_m
    port_area = math.pi * radius * radius
    coefficient = (
        of_ratio
        * geometry.fuel_density_kg_m3
        * burning_area
        * law.coefficient_si
        * port_area ** (-exponent)
    )
    return coefficient ** (1.0 / (1.0 - exponent))


def solve_variable_feed_coupling(
    tank_state: TankState,
    injector: Injector,
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    table: ThermochemistryTable,
    nozzle: NozzleGeometry,
    port_radius_m: float,
    *,
    c_star_efficiency: float = DEFAULT_C_STAR_EFFICIENCY,
    xtol: float = 1.0e-12,
    chamber_xtol: float = 1.0e-9,
    burnout_overshoot_fraction: float = 0.0,
) -> VariableFeedSolution:
    """Solve the injector / regression / variable-chemistry chamber loop."""
    radius = float(port_radius_m)
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError(f"port_radius_m must be finite and positive, got {port_radius_m!r}")
    overshoot = float(burnout_overshoot_fraction)
    if not math.isfinite(overshoot) or overshoot < 0.0:
        raise ValueError(
            f"burnout_overshoot_fraction must be finite and non-negative, got {overshoot!r}"
        )
    if radius > geometry.outer_diameter_m / 2.0 * (1.0 + max(overshoot, 1.0e-12)):
        raise ValueError(
            f"port_radius_m {radius!r} exceeds the outer grain radius; the grain is burnt through"
        )

    throat_area = nozzle.throat_area_m2
    floor_pressure = triple_point_pressure_pa()

    def not_flowing(status: VariableFeedStatus) -> VariableFeedSolution:
        return VariableFeedSolution(
            oxidizer_mass_flow_kg_s=0.0,
            fuel_mass_flow_kg_s=0.0,
            total_mass_flow_kg_s=0.0,
            mixture_ratio=math.nan,
            chamber_pressure_pa=0.0,
            tank_pressure_pa=tank_state.pressure_pa,
            pressure_margin_pa=tank_state.pressure_pa,
            regression_rate_m_s=0.0,
            oxidizer_mass_flux_si=0.0,
            thermochemistry=None,
            injector=None,
            feed_residual_kg_s=0.0,
            chamber_residual_pa=0.0,
            status=status,
        )

    if not tank_state.has_liquid:
        return not_flowing(VariableFeedStatus.LIQUID_DEPLETED)

    free_flow = injector.spi_mass_flow_kg_s(tank_state.liquid_density_kg_m3, tank_state.pressure_pa)
    if free_flow <= 0.0:
        return not_flowing(VariableFeedStatus.PRESSURE_EQUALIZED)

    # --- window on which the mixture ratio stays inside the table -------------
    of_low, of_high = table.of_bounds
    flow_low = oxidizer_flow_for_mixture_ratio(geometry, law, radius, of_low)
    flow_high = min(oxidizer_flow_for_mixture_ratio(geometry, law, radius, of_high), free_flow)
    if not (flow_low < flow_high):
        return not_flowing(VariableFeedStatus.OF_OUTSIDE_TABLE)

    def chamber_pressure_at_edge(oxidizer_flow: float, edge_pressure_pa: float) -> float:
        """``F(p_edge) - p_edge`` for the inner residual, evaluated on the table edge."""
        _, _, fuel = _fuel_flow_and_rate(geometry, law, radius, oxidizer_flow)
        state = table.evaluate(
            (oxidizer_flow / fuel) if fuel > 0.0 else math.inf,
            edge_pressure_pa,
            c_star_efficiency=c_star_efficiency,
        )
        return (oxidizer_flow + fuel) * state.c_star_m_s / throat_area - edge_pressure_pa

    # --- trim the window to where the inner chamber root is inside the table ---
    # Which constraint sets each edge is remembered, so that a feed root landing
    # outside the window is reported with the reason that actually bound it. The
    # mixture ratio and the chamber pressure can each bind independently, and a
    # root can sit comfortably inside the O/F span while still demanding a chamber
    # pressure the table does not cover.
    low_edge_reason = VariableFeedStatus.OF_OUTSIDE_TABLE
    # The upper edge starts as whichever of the two candidates min() actually
    # picked. If the injector's free-discharge flow is the binding one then a
    # root above the window is not a table limitation at all, and NO_ROOT is the
    # honest answer.
    high_edge_reason = (
        VariableFeedStatus.OF_OUTSIDE_TABLE if flow_high < free_flow else VariableFeedStatus.NO_ROOT
    )
    pressure_low, pressure_high = table.pressure_bounds_pa
    for edge_pressure, want_non_negative in (
        (pressure_low, True),
        (pressure_high, False),
    ):
        margin_low = chamber_pressure_at_edge(flow_low, edge_pressure)
        margin_high = chamber_pressure_at_edge(flow_high, edge_pressure)
        if want_non_negative:
            # Need F(p_min) >= p_min; the margin increases with flow, so trim below.
            if margin_low >= 0.0:
                continue
            if margin_high < 0.0:
                return not_flowing(VariableFeedStatus.PRESSURE_OUTSIDE_TABLE)
            flow_low = float(
                brentq(
                    chamber_pressure_at_edge,
                    flow_low,
                    flow_high,
                    args=(edge_pressure,),
                    xtol=xtol,
                    maxiter=200,
                )
            )
            low_edge_reason = VariableFeedStatus.PRESSURE_OUTSIDE_TABLE
        else:
            # Need F(p_max) <= p_max; trim above.
            if margin_high <= 0.0:
                continue
            if margin_low > 0.0:
                return not_flowing(VariableFeedStatus.PRESSURE_OUTSIDE_TABLE)
            flow_high = float(
                brentq(
                    chamber_pressure_at_edge,
                    flow_low,
                    flow_high,
                    args=(edge_pressure,),
                    xtol=xtol,
                    maxiter=200,
                )
            )
            high_edge_reason = VariableFeedStatus.PRESSURE_OUTSIDE_TABLE
        if not (flow_low < flow_high):
            return not_flowing(VariableFeedStatus.PRESSURE_OUTSIDE_TABLE)

    # The Brent endpoints land within ``xtol`` of the boundary and may fall a
    # rounding step outside it, so nudge them inwards before the outer solve.
    span = flow_high - flow_low
    flow_low += 1.0e-9 * span
    flow_high -= 1.0e-9 * span

    def chamber_at(oxidizer_flow: float):
        _, _, fuel = _fuel_flow_and_rate(geometry, law, radius, oxidizer_flow)
        return solve_variable_chamber(
            table,
            throat_area,
            oxidizer_flow,
            fuel,
            c_star_efficiency=c_star_efficiency,
            xtol=chamber_xtol,
        )

    invalid: list[VariableChamberStatus] = []

    def residual(oxidizer_flow: float) -> float:
        chamber = chamber_at(oxidizer_flow)
        if not chamber.is_solved:
            invalid.append(chamber.status)
            raise RuntimeError(
                "the inner chamber solve failed inside the pre-screened valid window "
                f"at m_dot_ox = {oxidizer_flow!r} kg/s with status {chamber.status.value}"
            )
        pressure = max(chamber.chamber_pressure_pa, floor_pressure)
        return injector.mass_flow(tank_state, pressure).mass_flow_kg_s - oxidizer_flow

    residual_low = residual(flow_low)
    residual_high = residual(flow_high)

    if residual_low <= 0.0:
        # The feed root sits at a lower oxidizer flow than the window allows, so
        # the operating point would need chemistry the table does not carry. The
        # reason reported is whichever constraint set that edge.
        return not_flowing(low_edge_reason)
    if residual_high > 0.0:
        # The root sits above the window: the operating point would need chemistry
        # the table does not carry, reported with the constraint that bound it.
        return not_flowing(high_edge_reason)

    oxidizer_flow = float(brentq(residual, flow_low, flow_high, xtol=xtol, maxiter=200))

    chamber = chamber_at(oxidizer_flow)
    if not chamber.is_solved:
        return not_flowing(
            VariableFeedStatus.PRESSURE_OUTSIDE_TABLE
            if chamber.status is VariableChamberStatus.PRESSURE_OUTSIDE_TABLE
            else VariableFeedStatus.OF_OUTSIDE_TABLE
        )

    flux, rate, fuel = _fuel_flow_and_rate(geometry, law, radius, oxidizer_flow)
    pressure = chamber.chamber_pressure_pa
    injector_result = injector.mass_flow(tank_state, max(pressure, floor_pressure))

    return VariableFeedSolution(
        oxidizer_mass_flow_kg_s=oxidizer_flow,
        fuel_mass_flow_kg_s=fuel,
        total_mass_flow_kg_s=chamber.total_mass_flow_kg_s,
        mixture_ratio=chamber.mixture_ratio,
        chamber_pressure_pa=pressure,
        tank_pressure_pa=tank_state.pressure_pa,
        pressure_margin_pa=tank_state.pressure_pa - pressure,
        regression_rate_m_s=rate,
        oxidizer_mass_flux_si=flux,
        thermochemistry=chamber.thermochemistry,
        injector=injector_result,
        feed_residual_kg_s=injector_result.mass_flow_kg_s - oxidizer_flow,
        chamber_residual_pa=chamber.residual_pa,
        status=VariableFeedStatus.FLOWING,
    )
