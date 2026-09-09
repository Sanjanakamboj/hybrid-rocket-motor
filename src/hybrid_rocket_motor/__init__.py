"""Reduced-order model of a generic N2O/HTPB hybrid rocket motor.

This package is an educational, reduced-order engineering study of a *generic*
nitrous-oxide / HTPB hybrid rocket motor.  It is not a motor design, a
fabrication guide, a test procedure, or a flight-hardware model.

Milestone 1 covers grain geometry, prescribed-oxidizer operating-point
bookkeeping, the empirical regression law ``r_dot = a * G_ox ** n``, and the
resulting fuel-flow and O/F bookkeeping.

Milestone 2 adds transient port evolution: integration of ``dr_p/dt = r_dot``
under a *prescribed* oxidizer mass-flow history, with fuel-depletion event
handling and time histories of regression rate, fuel flow and O/F.

Both milestones deliberately do **not** model N2O tank thermodynamics, injector
flow, chamber pressure, combustion efficiency, equilibrium chemistry, c*,
nozzle flow, thrust, structures, thermal response, ignition or fabrication.
"""

from __future__ import annotations

from .geometry import REPRESENTATIVE_GRAIN, GrainGeometry
from .operating_point import (
    OperatingPoint,
    fuel_mass_flow,
    mixture_ratio,
    oxidizer_mass_flux,
)
from .regression import (
    ILLUSTRATIVE_ANCHOR_FLUX_SI,
    ILLUSTRATIVE_SENSITIVITY_LAWS,
    REZAEI_2018_N2O_HTPB,
    MassFluxUnit,
    PowerLawRegressionLaw,
    Provenance,
    RegressionRateUnit,
    illustrative_exponent_variant,
)
from .transient import (
    CallableOxidizerFlow,
    ConstantOxidizerFlow,
    FlowSegment,
    FluxRangeStatus,
    OxidizerFlowHistory,
    PiecewiseConstantOxidizerFlow,
    TerminationReason,
    TransientResult,
    analytic_constant_flow_burnout_time,
    analytic_constant_flow_radius,
    analytic_constant_flow_speed_coefficient,
    remaining_fuel_mass_kg,
    simulate_transient,
)

__version__ = "0.2.0"

__all__ = [
    "ILLUSTRATIVE_ANCHOR_FLUX_SI",
    "ILLUSTRATIVE_SENSITIVITY_LAWS",
    "REPRESENTATIVE_GRAIN",
    "REZAEI_2018_N2O_HTPB",
    "CallableOxidizerFlow",
    "ConstantOxidizerFlow",
    "FlowSegment",
    "FluxRangeStatus",
    "GrainGeometry",
    "MassFluxUnit",
    "OperatingPoint",
    "OxidizerFlowHistory",
    "PiecewiseConstantOxidizerFlow",
    "PowerLawRegressionLaw",
    "Provenance",
    "RegressionRateUnit",
    "TerminationReason",
    "TransientResult",
    "__version__",
    "analytic_constant_flow_burnout_time",
    "analytic_constant_flow_radius",
    "analytic_constant_flow_speed_coefficient",
    "fuel_mass_flow",
    "illustrative_exponent_variant",
    "mixture_ratio",
    "oxidizer_mass_flux",
    "remaining_fuel_mass_kg",
    "simulate_transient",
]
