"""Quasi-steady chamber closure with O/F- and pressure-dependent properties.

Milestones 3 and 4 used a single prescribed ``c*``, which made the chamber
closure explicit::

    p_c = (m_dot_ox + m_dot_f) c* / A_t

Milestone 5 makes ``c*`` a function of the operating point, so the same balance
becomes **implicit** in the chamber pressure::

    p_c = m_dot_total c*_delivered(O/F, p_c) / A_t

The mixture ratio is fixed once the two flows are known, so only the pressure is
implicit here; the wider feed coupling is handled in
:mod:`hybrid_rocket_motor.variable_feed_system`.

Solution method
---------------
A bracketed Brent solve on

    residual(p_c) = m_dot_total c*_delivered(O/F, p_c) / A_t - p_c

over the table's own pressure range.  The derivative is
``(m_dot/A_t) dc*/dp_c - 1``; for this propellant ``c*`` varies by about 0.03 %
across the whole 10-45 bar span, so the first term is ~1e-3 and the residual is
strictly decreasing.  The root is therefore unique inside the table.

Nothing is clamped: if the implied pressure lies outside the tabulated range, or
the mixture ratio lies outside it, the solver reports an explicit failure status
rather than extrapolating the chemistry.

Units
-----
SI throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from scipy.optimize import brentq

from .thermochemistry import (
    DEFAULT_C_STAR_EFFICIENCY,
    TableStatus,
    ThermochemicalState,
    ThermochemistryTable,
)

__all__ = [
    "VariableChamberSolution",
    "VariableChamberStatus",
    "solve_variable_chamber",
]


class VariableChamberStatus(Enum):
    """Outcome of the implicit chamber solve."""

    SOLVED = "SOLVED"
    IDLE_NO_FLOW = "IDLE_NO_FLOW"
    OF_OUTSIDE_TABLE = "OF_OUTSIDE_TABLE"
    PRESSURE_OUTSIDE_TABLE = "PRESSURE_OUTSIDE_TABLE"
    NO_ROOT = "NO_ROOT"


@dataclass(frozen=True)
class VariableChamberSolution:
    """A converged variable-property chamber state."""

    oxidizer_mass_flow_kg_s: float
    fuel_mass_flow_kg_s: float
    total_mass_flow_kg_s: float
    mixture_ratio: float
    chamber_pressure_pa: float
    throat_area_m2: float
    thermochemistry: ThermochemicalState | None
    residual_pa: float
    status: VariableChamberStatus

    @property
    def is_solved(self) -> bool:
        return self.status is VariableChamberStatus.SOLVED

    @property
    def c_star_m_s(self) -> float:
        """Delivered characteristic velocity [m/s]."""
        return math.nan if self.thermochemistry is None else self.thermochemistry.c_star_m_s

    @property
    def c_star_identity_residual(self) -> float:
        """``p_c A_t - m_dot_total c*_delivered`` [N], zero in exact arithmetic."""
        if self.thermochemistry is None:
            return 0.0
        return (
            self.chamber_pressure_pa * self.throat_area_m2
            - self.total_mass_flow_kg_s * self.thermochemistry.c_star_m_s
        )


def _check_non_negative(name: str, value: float) -> float:
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if v < 0.0:
        raise ValueError(f"{name} must be non-negative, got {v!r}")
    return v


def solve_variable_chamber(
    table: ThermochemistryTable,
    throat_area_m2: float,
    oxidizer_mass_flow_kg_s: float,
    fuel_mass_flow_kg_s: float,
    *,
    c_star_efficiency: float = DEFAULT_C_STAR_EFFICIENCY,
    xtol: float = 1.0e-9,
) -> VariableChamberSolution:
    """Solve the implicit chamber pressure for a given pair of propellant flows.

    With no propellant flow at all the chamber is reported idle at zero pressure,
    matching the Milestone 3 convention; the mixture ratio is ``NaN`` whenever
    there is no fuel flow to form a ratio with.
    """
    area = float(throat_area_m2)
    if not math.isfinite(area) or area <= 0.0:
        raise ValueError(f"throat_area_m2 must be finite and positive, got {throat_area_m2!r}")
    m_ox = _check_non_negative("oxidizer_mass_flow_kg_s", oxidizer_mass_flow_kg_s)
    m_f = _check_non_negative("fuel_mass_flow_kg_s", fuel_mass_flow_kg_s)

    total = m_ox + m_f
    of_ratio = m_ox / m_f if m_f > 0.0 else math.nan

    def failed(status: VariableChamberStatus, pressure: float = 0.0) -> VariableChamberSolution:
        return VariableChamberSolution(
            oxidizer_mass_flow_kg_s=m_ox,
            fuel_mass_flow_kg_s=m_f,
            total_mass_flow_kg_s=total,
            mixture_ratio=of_ratio,
            chamber_pressure_pa=pressure,
            throat_area_m2=area,
            thermochemistry=None,
            residual_pa=math.nan,
            status=status,
        )

    if total == 0.0:
        return failed(VariableChamberStatus.IDLE_NO_FLOW)
    if not math.isfinite(of_ratio):
        return failed(VariableChamberStatus.OF_OUTSIDE_TABLE)

    of_low, of_high = table.of_bounds
    if not (of_low <= of_ratio <= of_high):
        return failed(VariableChamberStatus.OF_OUTSIDE_TABLE)

    pressure_low, pressure_high = table.pressure_bounds_pa

    def residual(pressure: float) -> float:
        state = table.evaluate(of_ratio, pressure, c_star_efficiency=c_star_efficiency)
        return total * state.c_star_m_s / area - pressure

    residual_low = residual(pressure_low)
    residual_high = residual(pressure_high)
    if residual_low < 0.0:
        # The implied pressure sits below the tabulated range.
        return failed(VariableChamberStatus.PRESSURE_OUTSIDE_TABLE, pressure_low)
    if residual_high > 0.0:
        return failed(VariableChamberStatus.PRESSURE_OUTSIDE_TABLE, pressure_high)

    pressure = float(brentq(residual, pressure_low, pressure_high, xtol=xtol, maxiter=200))
    state = table.evaluate(of_ratio, pressure, c_star_efficiency=c_star_efficiency)
    if state.status is not TableStatus.WITHIN_TABLE:  # pragma: no cover - guarded above
        return failed(VariableChamberStatus.PRESSURE_OUTSIDE_TABLE, pressure)

    return VariableChamberSolution(
        oxidizer_mass_flow_kg_s=m_ox,
        fuel_mass_flow_kg_s=m_f,
        total_mass_flow_kg_s=total,
        mixture_ratio=of_ratio,
        chamber_pressure_pa=pressure,
        throat_area_m2=area,
        thermochemistry=state,
        residual_pa=total * state.c_star_m_s / area - pressure,
        status=VariableChamberStatus.SOLVED,
    )
