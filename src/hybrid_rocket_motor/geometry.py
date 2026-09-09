"""Geometry of a generic single-port cylindrical hybrid fuel grain.

Scope
-----
This module describes a *generic, illustrative* cylindrical fuel grain with one
central circular port.  It is a reduced-order engineering abstraction used for a
regression-rate study; it is **not** a manufacturing drawing, a case/liner design,
or a structural model.  No wall thickness, liner, insulation, bulkhead, injector
or nozzle geometry is represented.

Units
-----
All quantities are SI: metres, square metres, cubic metres, kilograms,
kilograms per cubic metre.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "REPRESENTATIVE_GRAIN",
    "GrainGeometry",
]


@dataclass(frozen=True)
class GrainGeometry:
    """Immutable single-port cylindrical hybrid fuel grain.

    Parameters
    ----------
    length_m:
        Axial length ``L`` of the fuel grain [m].
    initial_port_diameter_m:
        Initial bore diameter ``D_p,0`` of the single circular port [m].
    outer_diameter_m:
        Outer diameter ``D_o`` of the fuel grain [m].
    fuel_density_kg_m3:
        Solid fuel density ``rho_f`` [kg/m^3].

    Notes
    -----
    The port is assumed to remain circular, concentric and axially uniform.  The
    geometry object itself is time-invariant: the *current* port diameter is
    passed explicitly to each method, so a single geometry can be queried at any
    hypothetical port size without mutation.
    """

    length_m: float
    initial_port_diameter_m: float
    outer_diameter_m: float
    fuel_density_kg_m3: float

    def __post_init__(self) -> None:
        for name in (
            "length_m",
            "initial_port_diameter_m",
            "outer_diameter_m",
            "fuel_density_kg_m3",
        ):
            value = getattr(self, name)
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite, got {value!r}")
            if value <= 0.0:
                raise ValueError(f"{name} must be strictly positive, got {value!r}")
        if self.initial_port_diameter_m >= self.outer_diameter_m:
            raise ValueError(
                "initial_port_diameter_m must be strictly less than outer_diameter_m "
                f"(got {self.initial_port_diameter_m} >= {self.outer_diameter_m})"
            )

    # -- validation helper -------------------------------------------------

    def _check_port_diameter(self, port_diameter_m: float) -> float:
        """Validate a candidate port diameter and return it as a float."""
        d = float(port_diameter_m)
        if not math.isfinite(d):
            raise ValueError(f"port_diameter_m must be finite, got {port_diameter_m!r}")
        if d <= 0.0:
            raise ValueError(f"port_diameter_m must be strictly positive, got {d!r}")
        if d >= self.outer_diameter_m:
            raise ValueError(
                "port_diameter_m must be strictly less than outer_diameter_m "
                f"(got {d} >= {self.outer_diameter_m}); the grain would be burnt through"
            )
        return d

    # -- derived geometry --------------------------------------------------

    @property
    def initial_web_thickness_m(self) -> float:
        """Initial radial fuel web ``(D_o - D_p,0) / 2`` [m]."""
        return 0.5 * (self.outer_diameter_m - self.initial_port_diameter_m)

    def port_area_m2(self, port_diameter_m: float) -> float:
        """Cross-sectional flow area of the circular port, ``A_port = pi D_p^2 / 4`` [m^2]."""
        d = self._check_port_diameter(port_diameter_m)
        return math.pi * d * d / 4.0

    def port_perimeter_m(self, port_diameter_m: float) -> float:
        """Wetted perimeter of the circular port, ``P_port = pi D_p`` [m]."""
        d = self._check_port_diameter(port_diameter_m)
        return math.pi * d

    def burning_area_m2(self, port_diameter_m: float) -> float:
        """Instantaneous burning surface area, ``A_burn = P_port * L`` [m^2].

        Only the cylindrical port wall is treated as burning surface; grain end
        faces are assumed inhibited.
        """
        return self.port_perimeter_m(port_diameter_m) * self.length_m

    def fuel_cross_section_area_m2(self, port_diameter_m: float) -> float:
        """Remaining solid-fuel annular cross-section, ``pi (D_o^2 - D_p^2) / 4`` [m^2]."""
        d = self._check_port_diameter(port_diameter_m)
        return math.pi * (self.outer_diameter_m**2 - d * d) / 4.0

    def fuel_volume_m3(self, port_diameter_m: float) -> float:
        """Remaining solid-fuel volume [m^3]."""
        return self.fuel_cross_section_area_m2(port_diameter_m) * self.length_m

    def fuel_mass_kg(self, port_diameter_m: float) -> float:
        """Remaining solid-fuel mass, ``rho_f * V_fuel`` [kg]."""
        return self.fuel_density_kg_m3 * self.fuel_volume_m3(port_diameter_m)

    # -- convenience -------------------------------------------------------

    @property
    def initial_port_area_m2(self) -> float:
        """Port area at the initial port diameter [m^2]."""
        return self.port_area_m2(self.initial_port_diameter_m)

    @property
    def initial_burning_area_m2(self) -> float:
        """Burning area at the initial port diameter [m^2]."""
        return self.burning_area_m2(self.initial_port_diameter_m)

    @property
    def initial_fuel_mass_kg(self) -> float:
        """Loaded fuel mass at the initial port diameter [kg]."""
        return self.fuel_mass_kg(self.initial_port_diameter_m)


#: Representative *illustrative* grain used throughout the Milestone 1 study.
#: These dimensions are chosen to be a plausible small laboratory-scale hybrid.
#: They are NOT taken from, nor intended to represent, any real motor.
REPRESENTATIVE_GRAIN = GrainGeometry(
    length_m=0.40,
    initial_port_diameter_m=0.040,
    outer_diameter_m=0.090,
    fuel_density_kg_m3=930.0,
)
