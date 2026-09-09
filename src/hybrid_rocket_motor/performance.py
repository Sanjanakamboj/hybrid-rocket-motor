"""Coupled chamber + nozzle performance and thrust prediction.

Scope (Milestone 3)
-------------------
Couples the Milestone 1/2 propellant flows to the quasi-steady chamber closure
(:mod:`hybrid_rocket_motor.chamber`) and the 1-D isentropic nozzle
(:mod:`hybrid_rocket_motor.nozzle`) to produce thrust.

The oxidizer mass flow is still **prescribed**: chamber pressure does not feed
back into it.  Combustion properties are **prescribed constants** -- ``c*``,
``gamma`` and ``T_c`` do not vary with O/F, and Milestone 3 contains no
equilibrium chemistry.  Consequently:

    The baseline Milestone 3 transient changes thrust because the total
    propellant flow evolves; ``c*``, ``gamma`` and ``T_c`` remain prescribed
    constants.

A structural consequence worth knowing
--------------------------------------
With ``epsilon`` and ``gamma`` fixed, the exit Mach number and both isentropic
ratios are constant in time.  Because ``T_c`` is also prescribed constant, the
exit temperature and therefore the **exit velocity are constant throughout a
burn**.  Since ``p_c`` is proportional to the total mass flow, thrust is exactly
*affine* in the total mass flow:

    F = m_dot_total [ V_e + c* epsilon (p_e/p_c) ] - p_a A_e

This identity is used as an independent check in the test suite.

Units
-----
SI throughout: s, kg/s, Pa, K, m/s, N, N·s.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .chamber import ChamberState, ChamberStatus, CombustionProperties, solve_chamber
from .nozzle import ExpansionRegime, NozzleExitState, NozzleGeometry, solve_nozzle_exit
from .transient import TransientResult

__all__ = [
    "STANDARD_GRAVITY_M_S2",
    "PerformancePoint",
    "ThrustHistory",
    "couple_transient_to_performance",
    "solve_performance_point",
]

#: Standard gravity, the defined constant used in the definition of specific
#: impulse in seconds (``I_sp = F / (m_dot g_0)``).
STANDARD_GRAVITY_M_S2 = 9.80665


@dataclass(frozen=True)
class PerformancePoint:
    """Coupled chamber + nozzle + thrust solution at one instant."""

    chamber: ChamberState
    exit_state: NozzleExitState | None
    nozzle: NozzleGeometry
    combustion: CombustionProperties
    ambient_pressure_pa: float
    momentum_thrust_n: float
    pressure_thrust_n: float
    thrust_n: float
    thrust_coefficient: float
    specific_impulse_s: float

    @property
    def is_firing(self) -> bool:
        return self.chamber.is_firing

    @property
    def chamber_pressure_pa(self) -> float:
        return self.chamber.chamber_pressure_pa

    @property
    def total_mass_flow_kg_s(self) -> float:
        return self.chamber.total_mass_flow_kg_s

    @property
    def expansion_regime(self) -> ExpansionRegime | None:
        return None if self.exit_state is None else self.exit_state.expansion_regime


def solve_performance_point(
    combustion: CombustionProperties,
    nozzle: NozzleGeometry,
    oxidizer_mass_flow_kg_s: float,
    fuel_mass_flow_kg_s: float,
    ambient_pressure_pa: float,
) -> PerformancePoint:
    """Solve chamber, nozzle exit and thrust for one instantaneous flow pair.

    Thrust follows the NASA Glenn rocket thrust equation

        F = m_dot V_e + (p_e - p_a) A_e

    with ``C_F = F / (p_c A_t)`` and ``I_sp = F / (m_dot g_0)``.

    Idle convention
    ---------------
    With no propellant flow the chamber is not firing.  Thrust is then exactly
    zero -- **not** ``-p_a A_e``, which would be the (meaningless) result of
    applying the pressure term to a nozzle that is simply full of ambient gas.
    No exit state is reported at all, so no nozzle numbers can be mistaken for a
    firing condition, and ``C_F`` and ``I_sp`` are ``NaN`` because both are
    ``0/0``.
    """
    p_a = float(ambient_pressure_pa)
    if not math.isfinite(p_a) or p_a < 0.0:
        raise ValueError(f"ambient_pressure_pa must be finite and non-negative, got {p_a!r}")

    chamber = solve_chamber(
        combustion, nozzle.throat_area_m2, oxidizer_mass_flow_kg_s, fuel_mass_flow_kg_s
    )

    if chamber.status is ChamberStatus.IDLE_NO_FLOW:
        return PerformancePoint(
            chamber=chamber,
            exit_state=None,
            nozzle=nozzle,
            combustion=combustion,
            ambient_pressure_pa=p_a,
            momentum_thrust_n=0.0,
            pressure_thrust_n=0.0,
            thrust_n=0.0,
            thrust_coefficient=math.nan,
            specific_impulse_s=math.nan,
        )

    exit_state = solve_nozzle_exit(
        nozzle,
        chamber.chamber_pressure_pa,
        combustion.chamber_temperature_k,
        combustion.gamma,
        combustion.gas_constant_j_kg_k,
        p_a,
    )

    momentum = chamber.total_mass_flow_kg_s * exit_state.velocity_m_s
    pressure = (exit_state.pressure_pa - p_a) * nozzle.exit_area_m2
    thrust = momentum + pressure

    return PerformancePoint(
        chamber=chamber,
        exit_state=exit_state,
        nozzle=nozzle,
        combustion=combustion,
        ambient_pressure_pa=p_a,
        momentum_thrust_n=momentum,
        pressure_thrust_n=pressure,
        thrust_n=thrust,
        thrust_coefficient=thrust / (chamber.chamber_pressure_pa * nozzle.throat_area_m2),
        specific_impulse_s=thrust / (chamber.total_mass_flow_kg_s * STANDARD_GRAVITY_M_S2),
    )


@dataclass(frozen=True)
class ThrustHistory:
    """Time histories from coupling a Milestone 2 transient to the M3 model.

    Exit-state arrays are ``NaN`` at idle samples so that no nozzle number can be
    read as if the motor were firing there; thrust is exactly ``0`` at those
    samples, and ``C_F`` / ``I_sp`` are ``NaN`` (both are ``0/0``).
    """

    time_s: np.ndarray
    oxidizer_mass_flow_kg_s: np.ndarray
    fuel_mass_flow_kg_s: np.ndarray
    total_mass_flow_kg_s: np.ndarray
    mixture_ratio: np.ndarray
    chamber_pressure_pa: np.ndarray
    exit_mach: np.ndarray
    exit_pressure_pa: np.ndarray
    exit_temperature_k: np.ndarray
    exit_velocity_m_s: np.ndarray
    momentum_thrust_n: np.ndarray
    pressure_thrust_n: np.ndarray
    thrust_n: np.ndarray
    thrust_coefficient: np.ndarray
    specific_impulse_s: np.ndarray
    firing: np.ndarray
    expansion_regime: tuple[ExpansionRegime | None, ...]
    ambient_pressure_pa: float
    nozzle: NozzleGeometry
    combustion: CombustionProperties
    transient: TransientResult

    # -- integrated performance -------------------------------------------

    @property
    def total_impulse_n_s(self) -> float:
        """``I_total = integral F dt`` over the simulated window [N s]."""
        return float(np.trapezoid(self.thrust_n, self.time_s))

    @property
    def active_burn_time_s(self) -> float:
        """Elapsed time with propellant actually flowing [s]."""
        return float(np.trapezoid(self.firing.astype(float), self.time_s))

    @property
    def mean_thrust_n(self) -> float:
        """Total impulse divided by the active burn time [N]."""
        active = self.active_burn_time_s
        return self.total_impulse_n_s / active if active > 0.0 else 0.0

    @property
    def peak_thrust_n(self) -> float:
        return float(np.max(self.thrust_n))

    @property
    def propellant_mass_kg(self) -> float:
        """Total oxidizer plus fuel consumed over the window [kg]."""
        return float(
            self.transient.cumulative_oxidizer_mass_kg[-1]
            + self.transient.cumulative_fuel_mass_kg[-1]
        )

    @property
    def equivalent_specific_impulse_s(self) -> float:
        """``I_sp,eq = I_total / (m_propellant g_0)`` [s].

        Preferred over averaging the instantaneous ``I_sp``: this has a clean
        integral definition and is well defined even when some samples are idle.
        """
        mass = self.propellant_mass_kg
        if mass <= 0.0:
            return math.nan
        return self.total_impulse_n_s / (mass * STANDARD_GRAVITY_M_S2)


def couple_transient_to_performance(
    transient: TransientResult,
    combustion: CombustionProperties,
    nozzle: NozzleGeometry,
    ambient_pressure_pa: float,
) -> ThrustHistory:
    """Apply the quasi-steady chamber/nozzle model to every transient sample.

    The coupling is one-way and instantaneous: each Milestone 2 state supplies
    ``m_dot_ox`` and ``m_dot_f``, and the chamber/nozzle model is solved at that
    state.  Nothing flows back into the regression or flow model.
    """
    p_a = float(ambient_pressure_pa)
    if not math.isfinite(p_a) or p_a < 0.0:
        raise ValueError(f"ambient_pressure_pa must be finite and non-negative, got {p_a!r}")

    points = [
        solve_performance_point(combustion, nozzle, float(m_ox), float(m_f), p_a)
        for m_ox, m_f in zip(
            transient.oxidizer_mass_flow_kg_s, transient.fuel_mass_flow_kg_s, strict=True
        )
    ]

    def gather(getter, default: float = math.nan) -> np.ndarray:
        return np.array(
            [default if p.exit_state is None else getter(p.exit_state) for p in points],
            dtype=float,
        )

    return ThrustHistory(
        time_s=transient.time_s,
        oxidizer_mass_flow_kg_s=transient.oxidizer_mass_flow_kg_s,
        fuel_mass_flow_kg_s=transient.fuel_mass_flow_kg_s,
        total_mass_flow_kg_s=np.array(
            [p.chamber.total_mass_flow_kg_s for p in points], dtype=float
        ),
        mixture_ratio=np.array([p.chamber.mixture_ratio for p in points], dtype=float),
        chamber_pressure_pa=np.array([p.chamber_pressure_pa for p in points], dtype=float),
        exit_mach=gather(lambda e: e.mach),
        exit_pressure_pa=gather(lambda e: e.pressure_pa),
        exit_temperature_k=gather(lambda e: e.temperature_k),
        exit_velocity_m_s=gather(lambda e: e.velocity_m_s),
        momentum_thrust_n=np.array([p.momentum_thrust_n for p in points], dtype=float),
        pressure_thrust_n=np.array([p.pressure_thrust_n for p in points], dtype=float),
        thrust_n=np.array([p.thrust_n for p in points], dtype=float),
        thrust_coefficient=np.array([p.thrust_coefficient for p in points], dtype=float),
        specific_impulse_s=np.array([p.specific_impulse_s for p in points], dtype=float),
        firing=np.array([p.is_firing for p in points], dtype=bool),
        expansion_regime=tuple(p.expansion_regime for p in points),
        ambient_pressure_pa=p_a,
        nozzle=nozzle,
        combustion=combustion,
        transient=transient,
    )
