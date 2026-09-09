"""One-dimensional steady isentropic converging-diverging nozzle.

Scope (Milestone 3)
-------------------
A *thermodynamic performance* model only.  There is no wall contour, no
material, no cooling, no manufacturing information, and no structural or thermal
sizing.  The flow model is 1-D, steady, adiabatic and isentropic with a
calorically perfect gas; shocks, boundary layers, flow separation, divergence
losses and two-phase effects are **not** modelled.

Equations and their sources
---------------------------
All relations below are the standard compressible-flow results as published by
NASA Glenn Research Center (see ``DESIGN.md`` §24 for the full audit):

* Area-Mach relation (``astar.html``)::

      A/A* = (1/M) [ (2/(g+1)) (1 + (g-1) M^2 / 2) ]^((g+1)/(2(g-1)))

* Isentropic static-to-total pressure ratio (``isentrop.html``, Eq. 6)::

      p/p_t = [1 + (g-1) M^2 / 2]^(-g/(g-1))

* Isentropic static-to-total temperature ratio (``isentrop.html``, Eq. 7)::

      T/T_t = [1 + (g-1) M^2 / 2]^(-1)

* Speed of sound and velocity (``mflchk.html``)::

      V = M a = M sqrt(g R T)

* Rocket thrust (``rockth.html``)::

      F = m_dot V_e + (p_e - p_a) A_e

Units
-----
SI throughout: m, m^2, Pa, K, kg/s, m/s, J/(kg K).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

import numpy as np
from scipy.optimize import brentq

__all__ = [
    "NEAR_IDEAL_PRESSURE_TOLERANCE",
    "REFERENCE_NOZZLE",
    "SUMMERFIELD_SEPARATION_RATIO",
    "ExpansionRegime",
    "NozzleExitState",
    "NozzleGeometry",
    "area_mach_ratio",
    "classify_expansion_regime",
    "isentropic_pressure_ratio",
    "isentropic_temperature_ratio",
    "solve_nozzle_exit",
    "solve_supersonic_exit_mach",
    "speed_of_sound",
]

#: ``|p_e/p_a - 1|`` below this counts as "near ideal" expansion.  A reporting
#: convention only; nothing in the physics depends on it.
NEAR_IDEAL_PRESSURE_TOLERANCE = 0.05

#: Summerfield's rule-of-thumb separation criterion: an overexpanded nozzle is
#: expected to separate when the exit pressure falls below roughly 0.4 of
#: ambient.  Used **only** to flag that a reported attached-flow result is
#: unlikely to be physical; this project does not model separation.
SUMMERFIELD_SEPARATION_RATIO = 0.4


class ExpansionRegime(Enum):
    """Qualitative nozzle expansion regime at the reported ambient pressure."""

    UNDEREXPANDED = "UNDEREXPANDED"
    NEAR_IDEAL_EXPANSION = "NEAR_IDEAL_EXPANSION"
    OVEREXPANDED = "OVEREXPANDED"
    VACUUM = "VACUUM"


def _check_gamma(gamma: float) -> float:
    g = float(gamma)
    if not math.isfinite(g):
        raise ValueError(f"gamma must be finite, got {gamma!r}")
    if g <= 1.0:
        raise ValueError(f"gamma must be strictly greater than 1, got {g!r}")
    return g


def _check_positive(name: str, value: float) -> float:
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if v <= 0.0:
        raise ValueError(f"{name} must be strictly positive, got {v!r}")
    return v


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NozzleGeometry:
    """Axisymmetric converging-diverging nozzle, described by two diameters.

    Only throat and exit areas matter to a 1-D model; no contour, length or
    divergence half-angle is represented, and none is needed.
    """

    throat_diameter_m: float
    exit_diameter_m: float

    def __post_init__(self) -> None:
        d_t = _check_positive("throat_diameter_m", self.throat_diameter_m)
        d_e = _check_positive("exit_diameter_m", self.exit_diameter_m)
        if d_e < d_t:
            raise ValueError(
                "exit_diameter_m must be at least the throat diameter "
                f"(got {d_e} < {d_t}); the diverging section cannot contract"
            )

    @classmethod
    def from_areas(cls, throat_area_m2: float, exit_area_m2: float) -> NozzleGeometry:
        """Build from areas rather than diameters."""
        a_t = _check_positive("throat_area_m2", throat_area_m2)
        a_e = _check_positive("exit_area_m2", exit_area_m2)
        if a_e < a_t:
            raise ValueError(f"exit_area_m2 must be at least throat_area_m2, got {a_e} < {a_t}")
        return cls(
            throat_diameter_m=math.sqrt(4.0 * a_t / math.pi),
            exit_diameter_m=math.sqrt(4.0 * a_e / math.pi),
        )

    @classmethod
    def from_throat_and_expansion_ratio(
        cls, throat_diameter_m: float, expansion_ratio: float
    ) -> NozzleGeometry:
        """Build from a throat diameter and an area expansion ratio ``A_e/A_t``."""
        d_t = _check_positive("throat_diameter_m", throat_diameter_m)
        eps = float(expansion_ratio)
        if not math.isfinite(eps):
            raise ValueError(f"expansion_ratio must be finite, got {expansion_ratio!r}")
        if eps < 1.0:
            raise ValueError(f"expansion_ratio must be at least 1, got {eps!r}")
        return cls(throat_diameter_m=d_t, exit_diameter_m=d_t * math.sqrt(eps))

    @property
    def throat_area_m2(self) -> float:
        """``A_t = pi D_t^2 / 4`` [m^2]."""
        return math.pi * self.throat_diameter_m**2 / 4.0

    @property
    def exit_area_m2(self) -> float:
        """``A_e = pi D_e^2 / 4`` [m^2]."""
        return math.pi * self.exit_diameter_m**2 / 4.0

    @property
    def expansion_ratio(self) -> float:
        """``epsilon = A_e/A_t`` [-]."""
        return self.exit_area_m2 / self.throat_area_m2


#: Reference *illustrative* nozzle used throughout Milestone 3.  A round 10 mm
#: throat with a round area expansion ratio of 4.  Chosen so that the chamber
#: pressure implied by the Milestone 1 mass flow lands inside the range measured
#: by the same experimental study the regression law comes from -- see
#: ``DESIGN.md`` §26 for the documented selection.  It is **not** optimised and
#: **not** taken from any real motor.
REFERENCE_NOZZLE = NozzleGeometry.from_throat_and_expansion_ratio(0.010, 4.0)


# ---------------------------------------------------------------------------
# Compressible-flow relations
# ---------------------------------------------------------------------------


def area_mach_ratio(mach, gamma: float):
    """Isentropic area ratio ``A/A*`` at a given Mach number [-].

    ``A/A* = (1/M) [ (2/(g+1)) (1 + (g-1) M^2/2) ]^((g+1)/(2(g-1)))``

    Equals exactly 1 at ``M = 1`` for any ``gamma``.  Accepts a scalar or an
    array-like and mirrors the input kind, like the rest of this package.
    """
    g = _check_gamma(gamma)
    m = np.asarray(mach, dtype=float)
    if np.any(~np.isfinite(m)):
        raise ValueError("mach must be finite")
    if np.any(m <= 0.0):
        raise ValueError("mach must be strictly positive")
    exponent = (g + 1.0) / (2.0 * (g - 1.0))
    ratio = (1.0 / m) * ((2.0 / (g + 1.0)) * (1.0 + (g - 1.0) * m * m / 2.0)) ** exponent
    if np.isscalar(mach) or np.ndim(mach) == 0:
        return float(ratio)
    return ratio


def isentropic_pressure_ratio(mach, gamma: float):
    """Static-to-total pressure ratio ``p/p_t = [1 + (g-1)M^2/2]^(-g/(g-1))`` [-]."""
    g = _check_gamma(gamma)
    m = np.asarray(mach, dtype=float)
    if np.any(~np.isfinite(m)) or np.any(m < 0.0):
        raise ValueError("mach must be finite and non-negative")
    ratio = (1.0 + (g - 1.0) * m * m / 2.0) ** (-g / (g - 1.0))
    if np.isscalar(mach) or np.ndim(mach) == 0:
        return float(ratio)
    return ratio


def isentropic_temperature_ratio(mach, gamma: float):
    """Static-to-total temperature ratio ``T/T_t = [1 + (g-1)M^2/2]^(-1)`` [-]."""
    g = _check_gamma(gamma)
    m = np.asarray(mach, dtype=float)
    if np.any(~np.isfinite(m)) or np.any(m < 0.0):
        raise ValueError("mach must be finite and non-negative")
    ratio = 1.0 / (1.0 + (g - 1.0) * m * m / 2.0)
    if np.isscalar(mach) or np.ndim(mach) == 0:
        return float(ratio)
    return ratio


def speed_of_sound(gamma: float, gas_constant_j_kg_k: float, temperature_k: float) -> float:
    """``a = sqrt(g R T)`` [m/s] for a calorically perfect gas."""
    g = _check_gamma(gamma)
    r = _check_positive("gas_constant_j_kg_k", gas_constant_j_kg_k)
    t = _check_positive("temperature_k", temperature_k)
    return math.sqrt(g * r * t)


def solve_supersonic_exit_mach(
    expansion_ratio: float,
    gamma: float,
    *,
    xtol: float = 1.0e-14,
    rtol: float = 8.881784197001252e-16,
    max_mach: float = 1.0e4,
) -> float:
    """Solve the **supersonic** branch of the area-Mach relation for ``M_e``.

    The area-Mach relation has two roots for any ``A/A* > 1`` -- one subsonic and
    one supersonic.  A converging-diverging nozzle running full and choked sits
    on the supersonic branch, which is monotonically increasing in ``M`` for
    ``M > 1``.  The root is therefore bracketed on ``[1, M_hi]`` and found with
    Brent's method, which is guaranteed to converge on a sign-changing bracket --
    unlike an unbracketed Newton iteration, which can wander onto the subsonic
    branch or diverge for large expansion ratios.

    ``epsilon = 1`` returns exactly ``M_e = 1`` (the throat limit).
    """
    g = _check_gamma(gamma)
    eps = float(expansion_ratio)
    if not math.isfinite(eps):
        raise ValueError(f"expansion_ratio must be finite, got {expansion_ratio!r}")
    if eps < 1.0:
        raise ValueError(f"expansion_ratio must be at least 1, got {eps!r}")
    if eps == 1.0:
        return 1.0

    def residual(mach: float) -> float:
        return area_mach_ratio(mach, g) - eps

    upper = 2.0
    while residual(upper) < 0.0:
        upper *= 2.0
        if upper > max_mach:
            raise ValueError(
                f"expansion_ratio {eps!r} is too large to bracket below Mach {max_mach}"
            )
    return float(brentq(residual, 1.0, upper, xtol=xtol, rtol=rtol, maxiter=200))


def classify_expansion_regime(
    exit_pressure_pa: float,
    ambient_pressure_pa: float,
    *,
    tolerance: float = NEAR_IDEAL_PRESSURE_TOLERANCE,
) -> ExpansionRegime:
    """Classify the nozzle expansion regime from the exit/ambient pressure pair.

    ``p_a = 0`` (a vacuum) is reported as :attr:`ExpansionRegime.VACUUM` rather
    than "infinitely underexpanded", because the relative comparison used for the
    other regimes is undefined there.
    """
    p_e = float(exit_pressure_pa)
    p_a = float(ambient_pressure_pa)
    if not math.isfinite(p_e) or not math.isfinite(p_a):
        raise ValueError("exit and ambient pressures must be finite")
    if p_a < 0.0:
        raise ValueError(f"ambient_pressure_pa must be non-negative, got {p_a!r}")
    if p_a == 0.0:
        return ExpansionRegime.VACUUM
    if abs(p_e / p_a - 1.0) <= tolerance:
        return ExpansionRegime.NEAR_IDEAL_EXPANSION
    return ExpansionRegime.UNDEREXPANDED if p_e > p_a else ExpansionRegime.OVEREXPANDED


# ---------------------------------------------------------------------------
# Exit state
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NozzleExitState:
    """Solved 1-D isentropic exit state of a choked converging-diverging nozzle."""

    mach: float
    pressure_pa: float
    temperature_k: float
    speed_of_sound_m_s: float
    velocity_m_s: float
    pressure_ratio: float
    temperature_ratio: float
    expansion_regime: ExpansionRegime
    ambient_pressure_pa: float
    area_mach_residual: float
    separation_expected: bool

    @property
    def is_separation_expected(self) -> bool:
        """True if Summerfield's rule of thumb suggests the flow would separate.

        This model assumes fully attached flow regardless; the flag exists so a
        physically doubtful result is reported as doubtful rather than presented
        as a validated nozzle-flow prediction.
        """
        return self.separation_expected


def solve_nozzle_exit(
    geometry: NozzleGeometry,
    chamber_pressure_pa: float,
    chamber_temperature_k: float,
    gamma: float,
    gas_constant_j_kg_k: float,
    ambient_pressure_pa: float,
    *,
    near_ideal_tolerance: float = NEAR_IDEAL_PRESSURE_TOLERANCE,
    xtol: float = 1.0e-14,
) -> NozzleExitState:
    """Solve the supersonic exit state of a choked, isentropic 1-D nozzle.

    Chamber conditions are treated as stagnation (total) conditions, which is the
    standard assumption for a low-Mach combustion chamber.
    """
    g = _check_gamma(gamma)
    p_c = _check_positive("chamber_pressure_pa", chamber_pressure_pa)
    t_c = _check_positive("chamber_temperature_k", chamber_temperature_k)
    r_gas = _check_positive("gas_constant_j_kg_k", gas_constant_j_kg_k)
    p_a = float(ambient_pressure_pa)
    if not math.isfinite(p_a) or p_a < 0.0:
        raise ValueError(f"ambient_pressure_pa must be finite and non-negative, got {p_a!r}")

    eps = geometry.expansion_ratio
    mach = solve_supersonic_exit_mach(eps, g, xtol=xtol)
    p_ratio = isentropic_pressure_ratio(mach, g)
    t_ratio = isentropic_temperature_ratio(mach, g)
    p_e = p_c * p_ratio
    t_e = t_c * t_ratio
    a_e = speed_of_sound(g, r_gas, t_e)
    v_e = mach * a_e

    regime = classify_expansion_regime(p_e, p_a, tolerance=near_ideal_tolerance)
    separation = bool(p_a > 0.0 and p_e < SUMMERFIELD_SEPARATION_RATIO * p_a)

    return NozzleExitState(
        mach=mach,
        pressure_pa=p_e,
        temperature_k=t_e,
        speed_of_sound_m_s=a_e,
        velocity_m_s=v_e,
        pressure_ratio=p_ratio,
        temperature_ratio=t_ratio,
        expansion_regime=regime,
        ambient_pressure_pa=p_a,
        area_mach_residual=area_mach_ratio(mach, g) - eps,
        separation_expected=separation,
    )
