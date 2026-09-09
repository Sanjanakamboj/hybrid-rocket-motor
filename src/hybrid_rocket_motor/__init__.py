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

Milestone 3 adds a quasi-steady chamber mass balance with a *prescribed*
characteristic velocity and a 1-D isentropic converging-diverging nozzle,
producing chamber pressure, exit state, thrust, thrust coefficient and specific
impulse -- still under a prescribed oxidizer mass flow.

Milestone 4 removes the prescribed oxidizer flow.  A saturated-equilibrium N2O
tank, a sourced injector model (SPI / HEM / Dyer-NHNE) and the Milestone 1-3
motor are solved as one coupled system, so oxidizer flow, chamber pressure and
thrust all fall out of a self-pressurising blowdown.

All four milestones deliberately do **not** model equilibrium combustion
chemistry, combustion instability, ignition or chamber-filling transients,
nozzle contours, shocks or flow separation, feed-line pressure losses,
structural design or thermal design.
"""

from __future__ import annotations

from .blowdown import (
    BURNOUT_OVERSHOOT_FRACTION,
    BlowdownResult,
    BlowdownTermination,
    simulate_blowdown,
)
from .chamber import (
    ILLUSTRATIVE_N2O_HTPB_COMBUSTION,
    UNIVERSAL_GAS_CONSTANT_J_MOL_K,
    ChamberState,
    ChamberStatus,
    CombustionProperties,
    ideal_characteristic_velocity,
    solve_chamber,
)
from .feed_system import (
    CHAMBER_PRESSURE_FLOOR_PA,
    FeedSolution,
    FeedStatus,
    solve_feed_coupling,
)
from .geometry import REPRESENTATIVE_GRAIN, GrainGeometry
from .injector import (
    DYER_REFERENCE_DISCHARGE_COEFFICIENT,
    Injector,
    InjectorFlowStatus,
    InjectorModel,
    InjectorResult,
)
from .nitrous_properties import (
    MAX_VERIFIED_TEMPERATURE_K,
    MIN_VERIFIED_TEMPERATURE_K,
    PropertyRangeStatus,
    SaturationState,
    saturated_state,
    saturation_pressure_pa,
)
from .nozzle import (
    NEAR_IDEAL_PRESSURE_TOLERANCE,
    REFERENCE_NOZZLE,
    SUMMERFIELD_SEPARATION_RATIO,
    ExpansionRegime,
    NozzleExitState,
    NozzleGeometry,
    area_mach_ratio,
    classify_expansion_regime,
    isentropic_pressure_ratio,
    isentropic_temperature_ratio,
    solve_nozzle_exit,
    solve_supersonic_exit_mach,
    speed_of_sound,
)
from .operating_point import (
    OperatingPoint,
    fuel_mass_flow,
    mixture_ratio,
    oxidizer_mass_flux,
)
from .performance import (
    STANDARD_GRAVITY_M_S2,
    PerformancePoint,
    ThrustHistory,
    couple_transient_to_performance,
    solve_performance_point,
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
from .tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
    NitrousTank,
    TankPhase,
    TankState,
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

__version__ = "0.4.0"

__all__ = [
    "BURNOUT_OVERSHOOT_FRACTION",
    "CHAMBER_PRESSURE_FLOOR_PA",
    "DYER_REFERENCE_DISCHARGE_COEFFICIENT",
    "ILLUSTRATIVE_ANCHOR_FLUX_SI",
    "ILLUSTRATIVE_N2O_HTPB_COMBUSTION",
    "ILLUSTRATIVE_SENSITIVITY_LAWS",
    "MAX_VERIFIED_TEMPERATURE_K",
    "MIN_VERIFIED_TEMPERATURE_K",
    "NEAR_IDEAL_PRESSURE_TOLERANCE",
    "REFERENCE_NOZZLE",
    "REFERENCE_TANK",
    "REFERENCE_TANK_FILL_FRACTION",
    "REFERENCE_TANK_TEMPERATURE_K",
    "REPRESENTATIVE_GRAIN",
    "REZAEI_2018_N2O_HTPB",
    "STANDARD_GRAVITY_M_S2",
    "SUMMERFIELD_SEPARATION_RATIO",
    "UNIVERSAL_GAS_CONSTANT_J_MOL_K",
    "BlowdownResult",
    "BlowdownTermination",
    "CallableOxidizerFlow",
    "ChamberState",
    "ChamberStatus",
    "CombustionProperties",
    "ConstantOxidizerFlow",
    "ExpansionRegime",
    "FeedSolution",
    "FeedStatus",
    "FlowSegment",
    "FluxRangeStatus",
    "GrainGeometry",
    "Injector",
    "InjectorFlowStatus",
    "InjectorModel",
    "InjectorResult",
    "MassFluxUnit",
    "NitrousTank",
    "NozzleExitState",
    "NozzleGeometry",
    "OperatingPoint",
    "OxidizerFlowHistory",
    "PerformancePoint",
    "PiecewiseConstantOxidizerFlow",
    "PowerLawRegressionLaw",
    "PropertyRangeStatus",
    "Provenance",
    "RegressionRateUnit",
    "SaturationState",
    "TankPhase",
    "TankState",
    "TerminationReason",
    "ThrustHistory",
    "TransientResult",
    "__version__",
    "analytic_constant_flow_burnout_time",
    "analytic_constant_flow_radius",
    "analytic_constant_flow_speed_coefficient",
    "area_mach_ratio",
    "classify_expansion_regime",
    "couple_transient_to_performance",
    "fuel_mass_flow",
    "ideal_characteristic_velocity",
    "illustrative_exponent_variant",
    "isentropic_pressure_ratio",
    "isentropic_temperature_ratio",
    "mixture_ratio",
    "oxidizer_mass_flux",
    "remaining_fuel_mass_kg",
    "saturated_state",
    "saturation_pressure_pa",
    "simulate_blowdown",
    "simulate_transient",
    "solve_chamber",
    "solve_feed_coupling",
    "solve_nozzle_exit",
    "solve_performance_point",
    "solve_supersonic_exit_mach",
    "speed_of_sound",
]
