"""Instantaneous operating-point bookkeeping for a single-port hybrid grain.

Scope (Milestone 1)
-------------------
The oxidizer mass flow is a **prescribed input only**.  Nothing in this module
derives it from tank pressure, injector characteristics, feed-system losses or
chamber pressure, and nothing here computes chamber pressure, c*, nozzle flow or
thrust.  The mixture ratio produced below is pure bookkeeping: it is *not* fed
into any combustion or performance calculation at this milestone.

Units
-----
SI throughout: kg/s, m, m^2, kg/(m^2 s), m/s, kg/m^3.  ``O/F`` is dimensionless.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .geometry import GrainGeometry
from .regression import PowerLawRegressionLaw

__all__ = [
    "OperatingPoint",
    "fuel_mass_flow",
    "mixture_ratio",
    "oxidizer_mass_flux",
]


def oxidizer_mass_flux(oxidizer_mass_flow_kg_s: float, port_area_m2: float) -> float:
    """Oxidizer mass flux ``G_ox = m_dot_ox / A_port`` [kg/(m^2 s)].

    The port area is the *oxidizer* flow area; by the usual convention the fuel
    added along the port is excluded, so this is an oxidizer flux and not a total
    mass flux.
    """
    m_dot = float(oxidizer_mass_flow_kg_s)
    area = float(port_area_m2)
    if not math.isfinite(m_dot):
        raise ValueError(f"oxidizer_mass_flow_kg_s must be finite, got {m_dot!r}")
    if m_dot < 0.0:
        raise ValueError(f"oxidizer_mass_flow_kg_s must be non-negative, got {m_dot!r}")
    if not math.isfinite(area) or area <= 0.0:
        raise ValueError(f"port_area_m2 must be finite and positive, got {area!r}")
    return m_dot / area


def fuel_mass_flow(
    fuel_density_kg_m3: float, burning_area_m2: float, regression_rate_m_s: float
) -> float:
    """Fuel mass flow ``m_dot_f = rho_f * A_burn * r_dot`` [kg/s].

    This is the instantaneous rate at which solid fuel is converted to gas by a
    surface receding at ``r_dot`` over an area ``A_burn``.
    """
    rho = float(fuel_density_kg_m3)
    area = float(burning_area_m2)
    rate = float(regression_rate_m_s)
    if not math.isfinite(rho) or rho <= 0.0:
        raise ValueError(f"fuel_density_kg_m3 must be finite and positive, got {rho!r}")
    if not math.isfinite(area) or area <= 0.0:
        raise ValueError(f"burning_area_m2 must be finite and positive, got {area!r}")
    if not math.isfinite(rate):
        raise ValueError(f"regression_rate_m_s must be finite, got {rate!r}")
    if rate < 0.0:
        raise ValueError(f"regression_rate_m_s must be non-negative, got {rate!r}")
    return rho * area * rate


def mixture_ratio(oxidizer_mass_flow_kg_s: float, fuel_mass_flow_kg_s: float) -> float:
    """Instantaneous mixture ratio ``O/F = m_dot_ox / m_dot_f`` [-].

    Undefined when there is no fuel flow, so a zero (or non-positive) fuel flow
    raises rather than returning ``inf``/``nan``.
    """
    m_ox = float(oxidizer_mass_flow_kg_s)
    m_f = float(fuel_mass_flow_kg_s)
    if not math.isfinite(m_ox) or m_ox < 0.0:
        raise ValueError(f"oxidizer_mass_flow_kg_s must be finite and non-negative, got {m_ox!r}")
    if not math.isfinite(m_f):
        raise ValueError(f"fuel_mass_flow_kg_s must be finite, got {m_f!r}")
    if m_f <= 0.0:
        raise ValueError(
            "mixture ratio is undefined for non-positive fuel mass flow "
            f"(got {m_f!r}); no fuel is being produced"
        )
    return m_ox / m_f


@dataclass(frozen=True)
class OperatingPoint:
    """One instantaneous, prescribed-oxidizer-flow operating point.

    Parameters
    ----------
    geometry:
        The :class:`~hybrid_rocket_motor.geometry.GrainGeometry` being operated.
    oxidizer_mass_flow_kg_s:
        Prescribed ``m_dot_ox`` [kg/s].  An input, never a derived quantity at
        this milestone.
    port_diameter_m:
        Current port diameter ``D_p`` [m].  Defaults to the grain's initial port
        diameter.

    Notes
    -----
    This object is a frozen snapshot.  It carries no time and performs no
    integration; port evolution over a burn is deferred to Milestone 2.
    """

    geometry: GrainGeometry
    oxidizer_mass_flow_kg_s: float
    port_diameter_m: float | None = None

    def __post_init__(self) -> None:
        if self.port_diameter_m is None:
            object.__setattr__(self, "port_diameter_m", self.geometry.initial_port_diameter_m)
        # Delegates range/positivity checks to the geometry.
        self.geometry.port_area_m2(self.port_diameter_m)
        m_dot = float(self.oxidizer_mass_flow_kg_s)
        if not math.isfinite(m_dot):
            raise ValueError(f"oxidizer_mass_flow_kg_s must be finite, got {m_dot!r}")
        if m_dot < 0.0:
            raise ValueError(f"oxidizer_mass_flow_kg_s must be non-negative, got {m_dot!r}")

    # -- geometry at this operating point ----------------------------------

    @property
    def port_area_m2(self) -> float:
        """Port flow area ``A_port`` [m^2]."""
        return self.geometry.port_area_m2(self.port_diameter_m)

    @property
    def port_perimeter_m(self) -> float:
        """Burning perimeter ``P_port`` [m]."""
        return self.geometry.port_perimeter_m(self.port_diameter_m)

    @property
    def burning_area_m2(self) -> float:
        """Burning surface area ``A_burn`` [m^2]."""
        return self.geometry.burning_area_m2(self.port_diameter_m)

    @property
    def remaining_fuel_mass_kg(self) -> float:
        """Solid fuel still present at this port diameter [kg]."""
        return self.geometry.fuel_mass_kg(self.port_diameter_m)

    # -- flow bookkeeping ---------------------------------------------------

    @property
    def oxidizer_mass_flux_si(self) -> float:
        """Oxidizer mass flux ``G_ox`` [kg/(m^2 s)]."""
        return oxidizer_mass_flux(self.oxidizer_mass_flow_kg_s, self.port_area_m2)

    def regression_rate_si(self, law: PowerLawRegressionLaw) -> float:
        """Regression rate [m/s] from a power-law correlation at this flux."""
        return law.regression_rate_si(self.oxidizer_mass_flux_si)

    def fuel_mass_flow_kg_s(self, law: PowerLawRegressionLaw) -> float:
        """Fuel mass flow ``m_dot_f`` [kg/s] implied by ``law`` at this point."""
        return fuel_mass_flow(
            self.geometry.fuel_density_kg_m3,
            self.burning_area_m2,
            self.regression_rate_si(law),
        )

    def mixture_ratio(self, law: PowerLawRegressionLaw) -> float:
        """Instantaneous ``O/F`` [-] implied by ``law`` at this point.

        Bookkeeping only.  Milestone 1 does not use ``O/F`` for flame
        temperature, c*, chamber pressure, nozzle flow or thrust.
        """
        return mixture_ratio(self.oxidizer_mass_flow_kg_s, self.fuel_mass_flow_kg_s(law))
