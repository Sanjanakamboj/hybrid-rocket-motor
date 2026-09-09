"""Rigid, adiabatic, saturated-equilibrium nitrous oxide tank.

Model
-----
The tank is a rigid control volume holding a *saturated* two-phase nitrous oxide
mixture that is assumed spatially uniform and in thermal equilibrium.  This is
the standard "equilibrium model" baseline for self-pressurising propellant tanks.
Liquid is withdrawn from the bottom; the remaining propellant cools because the
latent heat of the vapour it generates is drawn from the fluid itself, and the
tank pressure falls with it.

Governing equations
-------------------
Mass and energy conservation for a rigid, adiabatic control volume with a single
liquid outflow (the textbook open-system first law with ``Q = 0``, ``W = 0``):

    dm/dt = -m_dot_out
    dU/dt = -m_dot_out h_l(T)

closed by the volume constraint and the saturation relations:

    V     = m_l / rho_l(T) + m_v / rho_v(T)
    m     = m_l + m_v
    U     = m_l u_l(T) + m_v u_v(T)
    p     = p_sat(T)

Given ``(m, U)`` and the fixed volume ``V``, the temperature follows from a
single scalar root solve (see :meth:`NitrousTank.state_from_mass_and_energy`).
Nothing here prescribes a pressure-versus-time curve: the pressure history is an
*output* of the energy balance.

Assumptions, all of them consequential
--------------------------------------
* **Saturated equilibrium at all times** -- no thermal stratification, no
  metastable superheat, no bubble-nucleation delay.
* **Adiabatic walls** -- no heat transfer with the tank structure or ambient, and
  no thermal mass of the tank itself.
* **Spatially uniform** -- a 0-D lumped model with no ullage dynamics.
* **Saturated liquid leaves the tank** -- pure liquid outflow at ``h_l(T)``.

Known limitations
-----------------
Zimmerman (Stanford PhD, 2015) shows experimentally that self-pressurising tanks
exhibit an initial transient with rapid pressure fluctuations and bubble
nucleation that equilibrium models do not reproduce, followed by a
quasi-steady regime.  The equilibrium model is also reported to over-predict tank
pressure relative to experiment, and is outperformed by the Zilliac-Karabeyoglu
and Casalino-Pastrone models.  It is used here because it is the simplest
formulation that is a genuine thermodynamic closure rather than a fitted decay
curve, and because it needs saturation properties only.

Units
-----
SI throughout: m^3, kg, J, K, Pa.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from scipy.optimize import brentq

from .nitrous_properties import (
    PropertyRangeStatus,
    _update_saturated,
    critical_temperature_k,
    property_range_status,
    saturated_state,
    triple_point_temperature_k,
)
from .nitrous_properties import _check_temperature as _check_saturation_temperature

__all__ = [
    "REFERENCE_TANK",
    "SINGLE_PHASE_QUALITY_TOLERANCE",
    "NitrousTank",
    "TankPhase",
    "TankState",
]


#: How far the vapour quality may stray outside ``[0, 1]`` before the contents
#: are declared single-phase and the saturated model refused.  A small margin is
#: required so that a time integrator can bracket the terminal liquid-depletion
#: event, exactly as the burnout event needs a little room past the outer radius.
#: Reported phase masses are still clamped to be non-negative.
SINGLE_PHASE_QUALITY_TOLERANCE = 0.02


class TankPhase(Enum):
    """Phase condition of the tank contents."""

    TWO_PHASE = "TWO_PHASE"
    VAPOUR_ONLY = "VAPOUR_ONLY"
    DEPLETED = "DEPLETED"


def _check_positive(name: str, value: float) -> float:
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if v <= 0.0:
        raise ValueError(f"{name} must be strictly positive, got {v!r}")
    return v


@dataclass(frozen=True)
class TankState:
    """Instantaneous tank state."""

    temperature_k: float
    pressure_pa: float
    total_mass_kg: float
    liquid_mass_kg: float
    vapour_mass_kg: float
    liquid_volume_m3: float
    vapour_volume_m3: float
    vapour_quality: float
    liquid_fill_fraction: float
    total_internal_energy_j: float
    liquid_enthalpy_j_kg: float
    liquid_density_kg_m3: float
    liquid_entropy_j_kg_k: float
    phase: TankPhase
    range_status: PropertyRangeStatus
    unclamped_vapour_quality: float

    @property
    def liquid_fraction_remaining(self) -> float:
        """``1 - x`` using the *unclamped* quality [-].

        Crosses zero smoothly at liquid depletion and continues slightly negative
        beyond it, which is what a terminal-event root finder needs.  Use
        :attr:`liquid_mass_kg` for any physical quantity; that one is clamped.
        """
        return 1.0 - self.unclamped_vapour_quality

    @property
    def has_liquid(self) -> bool:
        return self.phase is TankPhase.TWO_PHASE and self.liquid_mass_kg > 0.0

    @property
    def is_within_verified_property_range(self) -> bool:
        return self.range_status is PropertyRangeStatus.VERIFIED


@dataclass(frozen=True)
class NitrousTank:
    """A rigid tank of fixed internal volume holding saturated nitrous oxide."""

    volume_m3: float

    def __post_init__(self) -> None:
        _check_positive("volume_m3", self.volume_m3)

    # -- construction ------------------------------------------------------

    def initial_state(self, temperature_k: float, liquid_fill_fraction: float) -> TankState:
        """Build the initial saturated state from temperature and liquid fill.

        ``liquid_fill_fraction`` is the fraction of the tank *volume* occupied by
        saturated liquid, which is how a tank is filled in practice.  It must lie
        strictly inside ``(0, 1)``: a tank with no ullage cannot self-pressurise
        and a tank with no liquid is outside this model's intended regime.
        """
        fill = float(liquid_fill_fraction)
        if not math.isfinite(fill):
            raise ValueError(f"liquid_fill_fraction must be finite, got {liquid_fill_fraction!r}")
        if not (0.0 < fill < 1.0):
            raise ValueError(
                f"liquid_fill_fraction must lie strictly in (0, 1), got {fill!r}"
            )
        sat = saturated_state(temperature_k)
        liquid_volume = fill * self.volume_m3
        vapour_volume = self.volume_m3 - liquid_volume
        liquid_mass = liquid_volume * sat.liquid_density_kg_m3
        vapour_mass = vapour_volume * sat.vapour_density_kg_m3
        total_mass = liquid_mass + vapour_mass
        total_energy = (
            liquid_mass * sat.liquid_internal_energy_j_kg
            + vapour_mass * sat.vapour_internal_energy_j_kg
        )
        return TankState(
            temperature_k=sat.temperature_k,
            pressure_pa=sat.pressure_pa,
            total_mass_kg=total_mass,
            liquid_mass_kg=liquid_mass,
            vapour_mass_kg=vapour_mass,
            liquid_volume_m3=liquid_volume,
            vapour_volume_m3=vapour_volume,
            vapour_quality=vapour_mass / total_mass,
            liquid_fill_fraction=fill,
            total_internal_energy_j=total_energy,
            liquid_enthalpy_j_kg=sat.liquid_enthalpy_j_kg,
            liquid_density_kg_m3=sat.liquid_density_kg_m3,
            liquid_entropy_j_kg_k=sat.liquid_entropy_j_kg_k,
            phase=TankPhase.TWO_PHASE,
            range_status=sat.range_status,
            unclamped_vapour_quality=vapour_mass / total_mass,
        )

    # -- the equilibrium flash --------------------------------------------

    @staticmethod
    def _flash_properties(temperature_k: float) -> tuple[float, float, float, float]:
        """``(v_l, v_v, u_l, u_v)`` at a saturation temperature.

        A lean four-call helper used inside the temperature root solve; the full
        :func:`~hybrid_rocket_motor.nitrous_properties.saturated_state` assembles
        eight properties and is needlessly expensive to call per iteration.
        """
        t = _check_saturation_temperature(temperature_k)
        liquid = _update_saturated(t, 0.0)
        v_l = 1.0 / float(liquid.rhomass())
        u_l = float(liquid.umass())
        vapour = _update_saturated(t, 1.0)
        return v_l, 1.0 / float(vapour.rhomass()), u_l, float(vapour.umass())

    def _quality_at(self, temperature_k: float, specific_volume_m3_kg: float) -> float:
        """Vapour quality of a saturated mixture of given specific volume.

        ``x = (v - v_l) / (v_v - v_l)``.  Since ``v_l`` rises and ``v_v`` falls
        with temperature, ``x`` increases monotonically with temperature at fixed
        ``v`` -- which is what makes the temperature solve well behaved.
        """
        sat = saturated_state(temperature_k)
        v_l = 1.0 / sat.liquid_density_kg_m3
        v_v = 1.0 / sat.vapour_density_kg_m3
        return (specific_volume_m3_kg - v_l) / (v_v - v_l)

    def _specific_internal_energy_at(
        self, temperature_k: float, specific_volume_m3_kg: float
    ) -> float:
        v_l, v_v, u_l, u_v = self._flash_properties(temperature_k)
        quality = (specific_volume_m3_kg - v_l) / (v_v - v_l)
        return u_l + quality * (u_v - u_l)

    def state_from_mass_and_energy(
        self,
        total_mass_kg: float,
        total_internal_energy_j: float,
        *,
        xtol: float = 1.0e-8,
        quality_overshoot: float = SINGLE_PHASE_QUALITY_TOLERANCE,
    ) -> TankState:
        """Solve the saturated equilibrium state from conserved ``(m, U)``.

        The temperature is found by a bracketed Brent solve on

            residual(T) = u_mix(T, v) - U/m ,   v = V/m

        which is monotonically increasing in ``T`` (see :meth:`_quality_at`), so
        the root is unique.  A bracketed solve is used rather than a fixed-point
        or Newton iteration because it cannot wander out of the two-phase dome.

        Raises if no two-phase root exists in the bracket -- that means the state
        is single-phase (liquid-full or vapour-only) and the caller must handle it
        rather than have a saturated equation applied where it does not hold.

        Parameters
        ----------
        quality_overshoot:
            How far the vapour quality may stray outside ``[0, 1]`` before the
            state is refused.  The default is the physical tolerance.  A time
            integrator passes a larger value so that it can step past liquid
            depletion far enough to bracket the terminal event; the reported
            phase masses remain clamped to be non-negative, and
            :func:`~hybrid_rocket_motor.blowdown.simulate_blowdown` never reports
            a sample beyond the event it locates.
        """
        mass = _check_positive("total_mass_kg", total_mass_kg)
        energy = float(total_internal_energy_j)
        if not math.isfinite(energy):
            raise ValueError(
                f"total_internal_energy_j must be finite, got {total_internal_energy_j!r}"
            )

        specific_volume = self.volume_m3 / mass
        specific_energy = energy / mass

        low = triple_point_temperature_k() + 1.0e-3
        high = critical_temperature_k() - 1.0e-3

        def residual(temperature: float) -> float:
            return self._specific_internal_energy_at(temperature, specific_volume) - specific_energy

        residual_low = residual(low)
        residual_high = residual(high)
        if residual_low > 0.0 or residual_high < 0.0:
            raise ValueError(
                "no saturated two-phase equilibrium state exists for "
                f"m = {mass!r} kg, U = {energy!r} J in V = {self.volume_m3!r} m^3 "
                "(the contents are single-phase); the caller must handle this case"
            )

        temperature = float(brentq(residual, low, high, xtol=xtol, maxiter=200))
        sat = saturated_state(temperature)
        unclamped_quality = self._quality_at(temperature, specific_volume)
        overshoot = float(quality_overshoot)
        if not math.isfinite(overshoot) or overshoot < 0.0:
            raise ValueError(
                f"quality_overshoot must be finite and non-negative, got {quality_overshoot!r}"
            )
        if unclamped_quality > 1.0 + overshoot:
            raise ValueError(
                f"the tank contents are single-phase vapour (quality {unclamped_quality:.4f}); "
                "the saturated two-phase model does not apply and no vapour-only "
                "discharge model exists in this project"
            )
        if unclamped_quality < -overshoot:
            raise ValueError(
                f"the tank contents are compressed liquid (quality {unclamped_quality:.4f}); "
                "the saturated two-phase model does not apply"
            )
        quality = min(max(unclamped_quality, 0.0), 1.0)

        vapour_mass = quality * mass
        liquid_mass = mass - vapour_mass
        liquid_volume = liquid_mass / sat.liquid_density_kg_m3
        vapour_volume = self.volume_m3 - liquid_volume

        phase = TankPhase.VAPOUR_ONLY if liquid_mass <= 0.0 else TankPhase.TWO_PHASE

        return TankState(
            temperature_k=temperature,
            pressure_pa=sat.pressure_pa,
            total_mass_kg=mass,
            liquid_mass_kg=liquid_mass,
            vapour_mass_kg=vapour_mass,
            liquid_volume_m3=liquid_volume,
            vapour_volume_m3=vapour_volume,
            vapour_quality=quality,
            liquid_fill_fraction=liquid_volume / self.volume_m3,
            total_internal_energy_j=energy,
            liquid_enthalpy_j_kg=sat.liquid_enthalpy_j_kg,
            liquid_density_kg_m3=sat.liquid_density_kg_m3,
            liquid_entropy_j_kg_k=sat.liquid_entropy_j_kg_k,
            phase=phase,
            range_status=property_range_status(temperature),
            unclamped_vapour_quality=unclamped_quality,
        )

    def liquid_mass_at(self, total_mass_kg: float, total_internal_energy_j: float) -> float:
        """Liquid mass for a conserved ``(m, U)`` pair [kg], for event detection."""
        return self.state_from_mass_and_energy(
            total_mass_kg, total_internal_energy_j
        ).liquid_mass_kg


#: Reference *illustrative* tank for the Milestone 4 study.  Round generic
#: numbers: a 10 litre rigid volume, 80 % liquid fill by volume, at 293.15 K.
#: Not a real, certified or commercially available tank, and not sized for any
#: pressure, material or safety requirement -- this project performs no
#: structural analysis whatsoever.
REFERENCE_TANK = NitrousTank(volume_m3=0.010)
REFERENCE_TANK_TEMPERATURE_K = 293.15
REFERENCE_TANK_FILL_FRACTION = 0.80
