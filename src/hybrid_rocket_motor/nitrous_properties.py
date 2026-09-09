"""Saturated nitrous oxide thermophysical properties.

Property backend
----------------
Properties come from **CoolProp**, whose nitrous oxide equation of state is the
short fundamental Helmholtz equation of state of

    E. W. Lemmon and R. Span, "Short Fundamental Equations of State for 20
    Industrial Fluids", J. Chem. Eng. Data 51(3) (2006) 785-850.
    (CoolProp reports this as BibTeX key ``Lemmon-JCED-2006``.)

This is the same equation of state the NIST Chemistry WebBook cites for nitrous
oxide, so the WebBook tabulation in ``tests/test_nitrous_properties.py`` is a
check on *this module's units and call structure*, not an independent check of
the equation of state itself.  That distinction is stated deliberately: agreeing
with NIST here does not validate the EOS, it validates the plumbing.

No property correlation is fitted, invented or transcribed anywhere in this
project.  Option A of the Milestone 4 property strategy (hand-transcribed ESDU
91022 correlations) was rejected because ESDU 91022 is paywalled and its
coefficients could not be verified from a primary source.

Scope
-----
Only *saturation-line* properties and the isentropic downstream states needed by
the homogeneous-equilibrium injector model are exposed.  There is no
supercritical model, no metastable/superheated-liquid model, no transport
property, and no solid phase.

Units
-----
SI throughout: K, Pa, kg/m^3, J/kg, J/(kg K).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

import CoolProp
import CoolProp.CoolProp as CP

__all__ = [
    "FLUID",
    "MAX_VERIFIED_TEMPERATURE_K",
    "MIN_VERIFIED_TEMPERATURE_K",
    "PropertyRangeStatus",
    "SaturationState",
    "critical_pressure_pa",
    "critical_temperature_k",
    "isentropic_downstream_enthalpy_j_kg",
    "property_range_status",
    "saturated_liquid_density_kg_m3",
    "saturated_liquid_enthalpy_j_kg",
    "saturated_liquid_entropy_j_kg_k",
    "saturated_liquid_internal_energy_j_kg",
    "saturated_state",
    "saturated_vapour_density_kg_m3",
    "saturated_vapour_enthalpy_j_kg",
    "saturated_vapour_internal_energy_j_kg",
    "saturation_pressure_pa",
    "triple_point_pressure_pa",
    "triple_point_temperature_k",
]

#: CoolProp fluid name.
FLUID = "NitrousOxide"

# CoolProp's low-level ``AbstractState`` interface is used instead of repeated
# ``PropsSI`` string calls.  One ``update`` fixes the thermodynamic state and every
# property is then a cheap getter, which matters because the tank temperature
# solve and the coupled feed solve each evaluate saturation properties tens of
# times per right-hand-side evaluation.  The equation of state is identical --
# this is purely a call-overhead optimisation, and the test suite checks the two
# interfaces agree.
_SATURATED = CoolProp.AbstractState("HEOS", FLUID)
_ISENTROPIC = CoolProp.AbstractState("HEOS", FLUID)


def _update_saturated(temperature_k: float, quality: float) -> CoolProp.AbstractState:
    """Fix the shared state object on the saturation line and return it."""
    _SATURATED.update(CoolProp.QT_INPUTS, quality, temperature_k)
    return _SATURATED

#: Lower bound of the temperature band over which this module is checked against
#: independently tabulated NIST Chemistry WebBook data [K].
MIN_VERIFIED_TEMPERATURE_K = 200.0

#: Upper bound of the checked band [K].  Deliberately below the 309.52 K critical
#: temperature: saturated liquid and vapour properties converge there and the
#: two-phase tank model degenerates well before the critical point is reached.
#: Set at 305 K so that a warm-day tank start (303.15 K) sits inside the band that
#: is actually checked against independently tabulated data.
MAX_VERIFIED_TEMPERATURE_K = 305.0


class PropertyRangeStatus(Enum):
    """Where a temperature sits relative to this module's checked band."""

    VERIFIED = "VERIFIED"
    BELOW_VERIFIED_RANGE = "BELOW_VERIFIED_RANGE"
    ABOVE_VERIFIED_RANGE = "ABOVE_VERIFIED_RANGE"


def critical_temperature_k() -> float:
    """Critical temperature of N2O [K] (309.52 K from the Lemmon-Span EOS)."""
    return float(CP.PropsSI("Tcrit", FLUID))


def critical_pressure_pa() -> float:
    """Critical pressure of N2O [Pa]."""
    return float(CP.PropsSI("pcrit", FLUID))


def triple_point_temperature_k() -> float:
    """Triple-point temperature of N2O [K]."""
    return float(CP.PropsSI("Ttriple", FLUID))


def triple_point_pressure_pa() -> float:
    """Triple-point pressure of N2O [Pa] (87 837 Pa from the Lemmon-Span EOS).

    The lower bound of the fluid domain: an isentropic expansion taken below this
    pressure leaves the liquid-vapour region entirely, so the homogeneous
    equilibrium injector model is not defined there.
    """
    return float(CP.PropsSI("ptriple", FLUID))


def _check_temperature(temperature_k: float) -> float:
    """Validate a saturation temperature and return it as a float.

    Rejects non-finite input, anything at or below the triple point, and anything
    at or above the critical temperature -- a saturated two-phase state does not
    exist outside that interval, and letting saturation equations run past the
    critical point silently is exactly the failure mode this guard exists to
    prevent.
    """
    t = float(temperature_k)
    if not math.isfinite(t):
        raise ValueError(f"temperature_k must be finite, got {temperature_k!r}")
    t_triple = triple_point_temperature_k()
    t_crit = critical_temperature_k()
    if t <= t_triple:
        raise ValueError(
            f"temperature_k must be above the N2O triple point {t_triple:.3f} K, got {t!r}"
        )
    if t >= t_crit:
        raise ValueError(
            f"temperature_k must be below the N2O critical point {t_crit:.4f} K "
            f"(no saturated two-phase state exists there), got {t!r}"
        )
    return t


def property_range_status(temperature_k: float) -> PropertyRangeStatus:
    """Classify a temperature against the checked band.

    Evaluation is *not* blocked outside the band -- the underlying equation of
    state remains valid up to the critical point -- but callers are given an
    explicit flag so that extrapolation beyond what this project has checked is
    never silent.
    """
    t = _check_temperature(temperature_k)
    if t < MIN_VERIFIED_TEMPERATURE_K:
        return PropertyRangeStatus.BELOW_VERIFIED_RANGE
    if t > MAX_VERIFIED_TEMPERATURE_K:
        return PropertyRangeStatus.ABOVE_VERIFIED_RANGE
    return PropertyRangeStatus.VERIFIED


def saturation_pressure_pa(temperature_k: float) -> float:
    """Saturation (vapour) pressure ``p_sat(T)`` [Pa]."""
    return float(_update_saturated(_check_temperature(temperature_k), 0.0).p())


def saturated_liquid_density_kg_m3(temperature_k: float) -> float:
    """Saturated liquid density ``rho_l(T)`` [kg/m^3]."""
    return float(_update_saturated(_check_temperature(temperature_k), 0.0).rhomass())


def saturated_vapour_density_kg_m3(temperature_k: float) -> float:
    """Saturated vapour density ``rho_v(T)`` [kg/m^3]."""
    return float(_update_saturated(_check_temperature(temperature_k), 1.0).rhomass())


def saturated_liquid_enthalpy_j_kg(temperature_k: float) -> float:
    """Saturated liquid specific enthalpy ``h_l(T)`` [J/kg]."""
    return float(_update_saturated(_check_temperature(temperature_k), 0.0).hmass())


def saturated_vapour_enthalpy_j_kg(temperature_k: float) -> float:
    """Saturated vapour specific enthalpy ``h_v(T)`` [J/kg]."""
    return float(_update_saturated(_check_temperature(temperature_k), 1.0).hmass())


def saturated_liquid_internal_energy_j_kg(temperature_k: float) -> float:
    """Saturated liquid specific internal energy ``u_l(T)`` [J/kg]."""
    return float(_update_saturated(_check_temperature(temperature_k), 0.0).umass())


def saturated_vapour_internal_energy_j_kg(temperature_k: float) -> float:
    """Saturated vapour specific internal energy ``u_v(T)`` [J/kg]."""
    return float(_update_saturated(_check_temperature(temperature_k), 1.0).umass())


def saturated_liquid_entropy_j_kg_k(temperature_k: float) -> float:
    """Saturated liquid specific entropy ``s_l(T)`` [J/(kg K)]."""
    return float(_update_saturated(_check_temperature(temperature_k), 0.0).smass())


def isentropic_downstream_enthalpy_j_kg(
    upstream_entropy_j_kg_k: float, downstream_pressure_pa: float
) -> tuple[float, float]:
    """Enthalpy and density after an isentropic expansion to a given pressure.

    Returns ``(h2, rho2)`` for the state at ``downstream_pressure_pa`` having the
    same specific entropy as upstream.  This is the state the homogeneous
    equilibrium injector model needs (Waxman et al., Eqns. 3-4); for nitrous
    oxide expanded below its saturation pressure the result is a two-phase
    mixture, which the equation of state handles directly.
    """
    s1 = float(upstream_entropy_j_kg_k)
    p2 = float(downstream_pressure_pa)
    if not math.isfinite(s1):
        raise ValueError(f"upstream_entropy_j_kg_k must be finite, got {upstream_entropy_j_kg_k!r}")
    if not math.isfinite(p2) or p2 <= 0.0:
        raise ValueError(
            f"downstream_pressure_pa must be finite and positive, got {downstream_pressure_pa!r}"
        )
    _ISENTROPIC.update(CoolProp.PSmass_INPUTS, p2, s1)
    return float(_ISENTROPIC.hmass()), float(_ISENTROPIC.rhomass())


@dataclass(frozen=True)
class SaturationState:
    """Complete saturated-line state at one temperature."""

    temperature_k: float
    pressure_pa: float
    liquid_density_kg_m3: float
    vapour_density_kg_m3: float
    liquid_enthalpy_j_kg: float
    vapour_enthalpy_j_kg: float
    liquid_internal_energy_j_kg: float
    vapour_internal_energy_j_kg: float
    liquid_entropy_j_kg_k: float
    range_status: PropertyRangeStatus

    @property
    def latent_heat_j_kg(self) -> float:
        """Latent heat of vaporisation ``h_v - h_l`` [J/kg]."""
        return self.vapour_enthalpy_j_kg - self.liquid_enthalpy_j_kg

    @property
    def is_within_verified_range(self) -> bool:
        return self.range_status is PropertyRangeStatus.VERIFIED


def saturated_state(temperature_k: float) -> SaturationState:
    """Assemble the full saturated state at a temperature (one guard, two updates)."""
    t = _check_temperature(temperature_k)
    liquid = _update_saturated(t, 0.0)
    pressure = float(liquid.p())
    liquid_density = float(liquid.rhomass())
    liquid_enthalpy = float(liquid.hmass())
    liquid_energy = float(liquid.umass())
    liquid_entropy = float(liquid.smass())
    vapour = _update_saturated(t, 1.0)
    return SaturationState(
        temperature_k=t,
        pressure_pa=pressure,
        liquid_density_kg_m3=liquid_density,
        vapour_density_kg_m3=float(vapour.rhomass()),
        liquid_enthalpy_j_kg=liquid_enthalpy,
        vapour_enthalpy_j_kg=float(vapour.hmass()),
        liquid_internal_energy_j_kg=liquid_energy,
        vapour_internal_energy_j_kg=float(vapour.umass()),
        liquid_entropy_j_kg_k=liquid_entropy,
        range_status=property_range_status(t),
    )
