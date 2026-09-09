"""Coupled injector / chamber closure.

The problem
-----------
Milestones 1-3 took the oxidizer mass flow as a *prescribed input*.  Milestone 4
removes that assumption, and in doing so creates an implicit loop:

* the injector flow depends on the chamber pressure it discharges into,
* the chamber pressure depends on the total propellant flow,
* the fuel flow depends on the oxidizer flow through the regression law.

Concretely, with the frozen Milestone 1-3 physics::

    G_ox   = m_dot_ox / A_port(r_p)
    r_dot  = a G_ox^n                                  (Milestone 1, unchanged)
    m_dot_f = rho_f (2 pi r_p L) r_dot                 (Milestone 2, unchanged)
    p_c    = (m_dot_ox + m_dot_f) c* / A_t             (Milestone 3, unchanged)
    m_dot_ox = injector(tank_state, p_c)               (Milestone 4, new)

Solution method
---------------
This is a scalar fixed point, solved as a bracketed root of

    residual(m_dot_ox) = injector_flow(tank, p_c(m_dot_ox, r_p)) - m_dot_ox

``p_c`` increases monotonically with ``m_dot_ox`` (both the oxidizer and the fuel
term increase).  The injector characteristic is **not** monotone in ``p_c``: the
homogeneous-equilibrium branch has a maximum with respect to backpressure -- the
critical-flow behaviour of Waxman et al. Eq. (5) -- and the Dyer blend inherits a
shallow maximum from it (about 3 % of the flow, near 10 bar for the reference
tank).  The uniqueness argument is therefore made on the loop gain rather than on
injector monotonicity:

    d(residual)/d(m_dot_ox) = (d m_dot_inj / d p_c)(d p_c / d m_dot_ox) - 1

``d p_c / d m_dot_ox`` is about ``c*/A_t`` ~ 1.9e7 Pa per kg/s, while
``d m_dot_inj / d p_c`` reaches at most ~ +4e-9 kg/s per Pa on the rising side of
the Dyer maximum.  Their product is ~ 0.07, so the derivative stays close to
``-1`` and the residual is strictly decreasing in practice even where the
injector characteristic is not.  ``test_feed_system.py`` verifies this
numerically rather than assuming it.

The bracket is available in closed form and needs no searching:

* at ``m_dot_ox = 0`` the chamber pressure is zero, so the injector delivers its
  free-discharge flow and ``residual > 0``;
* at ``m_dot_ox = m_dot_SPI(p_2 = 0)`` -- the largest flow any of the models can
  produce, since HEM and the Dyer blend are both bounded above by SPI at zero
  backpressure -- the chamber pressure is positive, so the injector delivers less
  than that and ``residual < 0``.

Brent's method is then guaranteed to converge.  A fixed-point iteration
``m_dot <- injector(p_c(m_dot))`` is deliberately *not* used: it has no
convergence guarantee here and can oscillate when the loop gain approaches unity.

Units
-----
SI throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from scipy.optimize import brentq

from .chamber import CombustionProperties, solve_chamber
from .geometry import GrainGeometry
from .injector import Injector, InjectorResult
from .nitrous_properties import triple_point_pressure_pa
from .nozzle import NozzleGeometry
from .regression import PowerLawRegressionLaw
from .tank import TankState

__all__ = [
    "CHAMBER_PRESSURE_FLOOR_PA",
    "FeedSolution",
    "FeedStatus",
    "solve_feed_coupling",
]

#: The injector's downstream pressure is floored at the N2O triple-point pressure
#: (87 837 Pa) purely so that the equation-of-state call at the hypothetical
#: ``m_dot_ox = 0`` bracket end is well posed.  The homogeneous equilibrium model
#: expands the fluid isentropically to the downstream pressure, and below the
#: triple point that isentrope leaves the liquid-vapour region entirely, so the
#: model has no meaning there.  Chamber pressures of interest here are 20-40 bar,
#: between one and two orders of magnitude above this floor, so it regularises
#: the bracket end and cannot influence the root [Pa].
CHAMBER_PRESSURE_FLOOR_PA = triple_point_pressure_pa()


class FeedStatus(Enum):
    """Outcome of the coupled feed/chamber solve."""

    FLOWING = "FLOWING"
    PRESSURE_EQUALIZED = "PRESSURE_EQUALIZED"
    LIQUID_DEPLETED = "LIQUID_DEPLETED"
    NO_ROOT = "NO_ROOT"
    PROPERTY_LIMIT = "PROPERTY_LIMIT"


@dataclass(frozen=True)
class FeedSolution:
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
    injector: InjectorResult | None
    residual_kg_s: float
    status: FeedStatus

    @property
    def is_flowing(self) -> bool:
        return self.status is FeedStatus.FLOWING


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


def solve_feed_coupling(
    tank_state: TankState,
    injector: Injector,
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    combustion: CombustionProperties,
    nozzle: NozzleGeometry,
    port_radius_m: float,
    *,
    xtol: float = 1.0e-12,
    rtol: float = 8.881784197001252e-16,
    burnout_overshoot_fraction: float = 0.0,
) -> FeedSolution:
    """Solve the injector/chamber/regression loop for the oxidizer mass flow.

    Returns a non-flowing solution -- with every flow exactly zero -- when the
    tank holds no liquid or when the tank and chamber pressures have equalised.
    No impossible coupled state is ever clipped into a valid-looking one.

    Parameters
    ----------
    burnout_overshoot_fraction:
        Fraction of the outer radius by which ``port_radius_m`` may exceed the
        grain's outer radius before being rejected.  Zero for ordinary callers,
        which is the physical constraint.  A time integrator must be allowed to
        evaluate *slightly* past the burnout radius in order to bracket the
        terminal burnout event, and the closure remains well posed there because
        it depends only on the port and burning areas -- both plain functions of
        the radius -- and never on the remaining fuel mass.  Radii reported by
        :func:`~hybrid_rocket_motor.blowdown.simulate_blowdown` are clamped back
        to the outer radius, so no overshoot ever reaches a result.
    """
    radius = float(port_radius_m)
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError(f"port_radius_m must be finite and positive, got {port_radius_m!r}")
    overshoot = float(burnout_overshoot_fraction)
    if not math.isfinite(overshoot) or overshoot < 0.0:
        raise ValueError(
            f"burnout_overshoot_fraction must be finite and non-negative, got {overshoot!r}"
        )
    radius_limit = geometry.outer_diameter_m / 2.0 * (1.0 + max(overshoot, 1.0e-12))
    if radius > radius_limit:
        raise ValueError(
            f"port_radius_m {radius!r} exceeds the outer grain radius; the grain is burnt through"
        )

    throat_area = nozzle.throat_area_m2

    def chamber_pressure(oxidizer_flow: float) -> float:
        _, _, fuel = _fuel_flow_and_rate(geometry, law, radius, oxidizer_flow)
        return solve_chamber(combustion, throat_area, oxidizer_flow, fuel).chamber_pressure_pa

    def residual(oxidizer_flow: float) -> float:
        p_c = max(chamber_pressure(oxidizer_flow), CHAMBER_PRESSURE_FLOOR_PA)
        return injector.mass_flow(tank_state, p_c).mass_flow_kg_s - oxidizer_flow

    def not_flowing(status: FeedStatus) -> FeedSolution:
        return FeedSolution(
            oxidizer_mass_flow_kg_s=0.0,
            fuel_mass_flow_kg_s=0.0,
            total_mass_flow_kg_s=0.0,
            mixture_ratio=math.nan,
            chamber_pressure_pa=0.0,
            tank_pressure_pa=tank_state.pressure_pa,
            pressure_margin_pa=tank_state.pressure_pa,
            regression_rate_m_s=0.0,
            oxidizer_mass_flux_si=0.0,
            injector=None,
            residual_kg_s=0.0,
            status=status,
        )

    if not tank_state.has_liquid:
        return not_flowing(FeedStatus.LIQUID_DEPLETED)

    # Upper bracket: the SPI free-discharge flow bounds every model from above.
    upper = injector.spi_mass_flow_kg_s(tank_state.liquid_density_kg_m3, tank_state.pressure_pa)
    if upper <= 0.0:
        return not_flowing(FeedStatus.PRESSURE_EQUALIZED)

    residual_low = residual(0.0)
    if residual_low <= 0.0:
        # The injector cannot deliver anything even against a vacuum chamber.
        return not_flowing(FeedStatus.PRESSURE_EQUALIZED)

    residual_high = residual(upper)
    if residual_high > 0.0:
        return not_flowing(FeedStatus.NO_ROOT)

    oxidizer_flow = float(brentq(residual, 0.0, upper, xtol=xtol, rtol=rtol, maxiter=200))

    flux, rate, fuel = _fuel_flow_and_rate(geometry, law, radius, oxidizer_flow)
    chamber = solve_chamber(combustion, throat_area, oxidizer_flow, fuel)
    p_c = chamber.chamber_pressure_pa
    injector_result = injector.mass_flow(tank_state, max(p_c, CHAMBER_PRESSURE_FLOOR_PA))

    return FeedSolution(
        oxidizer_mass_flow_kg_s=oxidizer_flow,
        fuel_mass_flow_kg_s=fuel,
        total_mass_flow_kg_s=chamber.total_mass_flow_kg_s,
        mixture_ratio=chamber.mixture_ratio,
        chamber_pressure_pa=p_c,
        tank_pressure_pa=tank_state.pressure_pa,
        pressure_margin_pa=tank_state.pressure_pa - p_c,
        regression_rate_m_s=rate,
        oxidizer_mass_flux_si=flux,
        injector=injector_result,
        residual_kg_s=injector_result.mass_flow_kg_s - oxidizer_flow,
        status=FeedStatus.FLOWING,
    )
