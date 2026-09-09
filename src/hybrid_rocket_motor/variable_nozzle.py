"""Adapter that drives the frozen Milestone 3 nozzle with tabulated properties.

Milestone 3 expanded the flow through a fixed-area nozzle using a *prescribed*
:class:`~hybrid_rocket_motor.chamber.CombustionProperties` -- one ``gamma``, one
``T_c``, one molar mass for the whole burn.  Milestone 5 changes nothing about
that expansion: it changes where the numbers come from.  The isentropic relations
are re-used unmodified, and this module is only the translation layer that turns
an interpolated :class:`~hybrid_rocket_motor.thermochemistry.ThermochemicalState`
into the object the frozen code already accepts.

Why the nozzle cares about ``gamma`` at all
-------------------------------------------
The exit Mach number at a fixed area ratio is the root of the area--Mach relation

    eps = (1/M) [ (2/(g+1)) (1 + (g-1)/2 M^2) ]^((g+1)/(2(g-1)))

which contains ``gamma`` and nothing else.  So a mixture ratio that shifts
``gamma`` shifts ``M_e``, and with it the exit pressure ratio, the exit velocity
and the pressure-thrust term -- *even at constant chamber pressure and constant
area ratio*.  This is a genuinely separate channel from the ``c*`` effect on
chamber pressure, and the Milestone 5 study reports the two separately rather
than only their sum.

Efficiency bookkeeping
----------------------
The table stores **ideal** ``c*``; the single ``eta_c*`` of Milestone 3 is
applied once, in :mod:`hybrid_rocket_motor.thermochemistry`, and the delivered
value is carried through here.  Nothing in this module reapplies it, so the loss
is not double-counted.  There is still **no nozzle efficiency**: the expansion
remains ideal, one-dimensional and shock-free, exactly as in Milestone 3.

Units
-----
SI throughout.
"""

from __future__ import annotations

from .chamber import CombustionProperties
from .nozzle import NozzleGeometry
from .performance import PerformancePoint, solve_performance_point
from .regression import Provenance
from .thermochemistry import ThermochemicalState

__all__ = [
    "combustion_properties_from_state",
    "solve_variable_performance_point",
]


def combustion_properties_from_state(
    state: ThermochemicalState,
    *,
    label: str | None = None,
    reference: str | None = None,
) -> CombustionProperties:
    """Package one interpolated state as frozen-Milestone-3 combustion properties.

    The delivered ``c*`` is carried across unchanged, so feeding the result back
    through :func:`~hybrid_rocket_motor.chamber.solve_chamber` at the same total
    mass flow reproduces the chamber pressure the variable-property solver
    converged on.  That identity is asserted in the tests.

    The provenance is marked :attr:`~hybrid_rocket_motor.regression.Provenance.SOURCED`
    rather than ``ILLUSTRATIVE``: unlike the Milestone 3 constants, these numbers
    come from a NASA CEA equilibrium solution.  Sourced still does not mean
    validated against a firing.
    """
    described = label or f"CEA equilibrium at O/F = {state.of_ratio:.4f}"
    return CombustionProperties(
        c_star_m_s=state.c_star_m_s,
        gamma=state.gamma,
        chamber_temperature_k=state.chamber_temperature_k,
        molar_mass_kg_mol=state.effective_molar_mass_kg_mol,
        label=described,
        reference=(
            reference
            or "NASA CEA (Gordon & McBride, NASA RP-1311) via the frozen project table; "
            "see data/n2o_htpb_equilibrium_metadata.json"
        ),
        provenance=Provenance.SOURCED,
    )


def solve_variable_performance_point(
    state: ThermochemicalState,
    nozzle: NozzleGeometry,
    oxidizer_mass_flow_kg_s: float,
    fuel_mass_flow_kg_s: float,
    ambient_pressure_pa: float,
) -> PerformancePoint:
    """Frozen Milestone 3 thrust solve, driven by one interpolated state."""
    return solve_performance_point(
        combustion_properties_from_state(state),
        nozzle,
        oxidizer_mass_flow_kg_s,
        fuel_mass_flow_kg_s,
        ambient_pressure_pa,
    )
