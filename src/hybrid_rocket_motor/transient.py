"""Transient port regression under a prescribed oxidizer mass-flow history.

Scope (Milestone 2)
-------------------
This module integrates the growth of a single circular port in time.  The
oxidizer mass flow remains a **prescribed input**: nothing here derives it from
N2O tank state, vapour pressure, two-phase blowdown, feed-line pressure drop,
injector characteristics or chamber pressure.  Nothing here computes chamber
pressure, combustion equilibrium, c*, nozzle flow, thrust or specific impulse.

State equation
--------------
The regression rate ``r_dot`` is the *recession velocity of the fuel surface*
(Marquardt & Majdalani 2020, Eq. 1: ``Q_w = rho_f r_dot h_v``, so ``rho_f r_dot``
is the fuel mass flux leaving the surface).  For a circular port the surface
normal is radial, so the port radius grows at exactly the regression rate:

    dr_p/dt = r_dot              equivalently   dD_p/dt = 2 r_dot

Karabeyoglu, Cantwell & Zilliac (JPP 23(4), 2007) state this directly --
"Based on the definition of regression rate, ``2 r_dot = dD/dt``" -- immediately
before their port-diameter equation (their Eq. 51), derived for the same case of
a single circular port and under an explicit quasi-steady assumption.

Augmented state
---------------
Cumulative masses are carried as extra ODE states so that they inherit the
integrator's error control rather than being post-processed with a cruder rule:

    dr_p/dt      = r_dot(G_ox(t, r_p))
    dm_ox,cum/dt = m_dot_ox(t)
    dm_f,cum/dt  = rho_f * A_burn(r_p) * r_dot

Because ``m_dot_f = rho_f (2 pi r_p L) dr_p/dt = rho_f L pi d(r_p^2)/dt``, the
integrated fuel mass must equal the geometric fuel loss
``rho_f L pi (r_p^2 - r_p,0^2)`` *exactly* in exact arithmetic.  Any difference is
purely numerical, which makes the mass-closure residual a sharp, independent
check on the integration -- the same global mass balance Karabeyoglu et al.
recommend for checking numerical solutions (their Eq. 63).

Units
-----
SI throughout: s, m, m^2, kg, kg/s, kg/(m^2 s), m/s.  Regression rate is
additionally exposed in mm/s for readability only.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from enum import Enum

import numpy as np
from scipy.integrate import solve_ivp

from .geometry import GrainGeometry
from .regression import PowerLawRegressionLaw

__all__ = [
    "CallableOxidizerFlow",
    "ConstantOxidizerFlow",
    "FlowSegment",
    "FluxRangeStatus",
    "OxidizerFlowHistory",
    "PiecewiseConstantOxidizerFlow",
    "TerminationReason",
    "TransientResult",
    "analytic_constant_flow_burnout_time",
    "analytic_constant_flow_radius",
    "analytic_constant_flow_speed_coefficient",
    "remaining_fuel_mass_kg",
    "simulate_transient",
]


# ---------------------------------------------------------------------------
# Enumerations
# ---------------------------------------------------------------------------


class TerminationReason(Enum):
    """Why the integration stopped."""

    COMPLETED_TIME_WINDOW = "COMPLETED_TIME_WINDOW"
    FUEL_DEPLETED = "FUEL_DEPLETED"


class FluxRangeStatus(Enum):
    """Where an evaluated oxidizer mass flux sits relative to the source data.

    ``ZERO_OXIDIZER_FLOW`` is kept separate from ``BELOW_SOURCE_FLUX_RANGE``
    deliberately: at ``m_dot_ox = 0`` the model is not extrapolating an empirical
    correlation at all, it is returning the exact ``r_dot = 0`` limit.  Lumping
    the two together would overstate how much of a burn is extrapolated.
    """

    WITHIN_SOURCE_FLUX_RANGE = "WITHIN_SOURCE_FLUX_RANGE"
    BELOW_SOURCE_FLUX_RANGE = "BELOW_SOURCE_FLUX_RANGE"
    ABOVE_SOURCE_FLUX_RANGE = "ABOVE_SOURCE_FLUX_RANGE"
    ZERO_OXIDIZER_FLOW = "ZERO_OXIDIZER_FLOW"
    NO_SOURCE_RANGE_DECLARED = "NO_SOURCE_RANGE_DECLARED"


# ---------------------------------------------------------------------------
# Prescribed oxidizer-flow histories
# ---------------------------------------------------------------------------


def _validate_time(name: str, value: float) -> float:
    t = float(value)
    if not math.isfinite(t):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if t < 0.0:
        raise ValueError(f"{name} must be non-negative, got {t!r}")
    return t


def _validate_flow(name: str, value: float) -> float:
    m = float(value)
    if not math.isfinite(m):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if m < 0.0:
        raise ValueError(f"{name} must be non-negative, got {m!r}")
    return m


class OxidizerFlowHistory:
    """Base class for a prescribed ``m_dot_ox(t)``.

    Subclasses expose the flow at a time, the window over which the history is
    defined, and the interior breakpoints at which the flow is discontinuous.
    The solver restarts integration at every breakpoint, so discontinuous
    schedules are integrated exactly rather than smeared by an adaptive step.
    """

    t_start_s: float
    t_end_s: float

    def mass_flow_kg_s(self, time_s: float) -> float:  # pragma: no cover - interface
        raise NotImplementedError

    @property
    def interior_breakpoints_s(self) -> tuple[float, ...]:
        """Interior times at which the prescribed flow is discontinuous."""
        return ()

    def flow_on_segment(self, seg_start_s: float, seg_end_s: float) -> Callable[[float], float]:
        """Return the flow law valid *strictly inside* one integration segment.

        The solver integrates segment by segment between breakpoints, and an
        integrator legitimately evaluates the right-hand side at a segment's
        right endpoint.  A plain half-open lookup would there return the *next*
        segment's flow and contaminate the final step, so the solver always asks
        for a segment-restricted accessor instead.  Segments are cut at the
        history's own breakpoints, so the flow is continuous inside one and
        clamping the query time is exact.
        """

        def flow(time_s: float) -> float:
            clamped = min(max(float(time_s), seg_start_s), seg_end_s)
            return self.mass_flow_kg_s(clamped)

        return flow

    def _check_window(self) -> None:
        if self.t_end_s <= self.t_start_s:
            raise ValueError(
                f"t_end_s must be strictly greater than t_start_s, got "
                f"{self.t_start_s!r} and {self.t_end_s!r}"
            )


@dataclass(frozen=True)
class ConstantOxidizerFlow(OxidizerFlowHistory):
    """A constant prescribed oxidizer mass flow over ``[t_start_s, t_end_s]``."""

    mass_flow_kg_s_value: float
    t_end_s: float
    t_start_s: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "mass_flow_kg_s_value", _validate_flow("mass_flow_kg_s_value",
                                                         self.mass_flow_kg_s_value)
        )
        object.__setattr__(self, "t_start_s", _validate_time("t_start_s", self.t_start_s))
        object.__setattr__(self, "t_end_s", _validate_time("t_end_s", self.t_end_s))
        self._check_window()

    def mass_flow_kg_s(self, time_s: float) -> float:
        return self.mass_flow_kg_s_value


@dataclass(frozen=True)
class FlowSegment:
    """One interval of a piecewise-constant prescribed flow schedule."""

    t_start_s: float
    t_end_s: float
    mass_flow_kg_s: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "t_start_s", _validate_time("t_start_s", self.t_start_s))
        object.__setattr__(self, "t_end_s", _validate_time("t_end_s", self.t_end_s))
        object.__setattr__(
            self, "mass_flow_kg_s", _validate_flow("mass_flow_kg_s", self.mass_flow_kg_s)
        )
        if self.t_end_s <= self.t_start_s:
            raise ValueError(
                f"segment must have positive duration, got [{self.t_start_s}, {self.t_end_s}]"
            )


@dataclass(frozen=True)
class PiecewiseConstantOxidizerFlow(OxidizerFlowHistory):
    """A contiguous, strictly increasing schedule of constant-flow segments.

    The schedule must be gapless and non-overlapping: segment ``i+1`` must start
    exactly where segment ``i`` ends.  Each segment is half-open ``[start, end)``
    so the flow is single-valued everywhere except at the final time, which takes
    the last segment's value.
    """

    segments: tuple[FlowSegment, ...]
    t_start_s: float = field(init=False)
    t_end_s: float = field(init=False)

    def __post_init__(self) -> None:
        segments = tuple(self.segments)
        if not segments:
            raise ValueError("schedule must contain at least one segment")
        for earlier, later in itertools.pairwise(segments):
            if later.t_start_s != earlier.t_end_s:
                raise ValueError(
                    "schedule segments must be contiguous and strictly increasing: "
                    f"segment ending at {earlier.t_end_s} is followed by one starting at "
                    f"{later.t_start_s}"
                )
        object.__setattr__(self, "segments", segments)
        object.__setattr__(self, "t_start_s", segments[0].t_start_s)
        object.__setattr__(self, "t_end_s", segments[-1].t_end_s)
        self._check_window()

    @classmethod
    def from_durations(
        cls,
        steps: Sequence[tuple[float, float]],
        t_start_s: float = 0.0,
    ) -> PiecewiseConstantOxidizerFlow:
        """Build from ``(duration_s, mass_flow_kg_s)`` pairs."""
        if not steps:
            raise ValueError("at least one (duration, mass flow) step is required")
        cursor = _validate_time("t_start_s", t_start_s)
        segments: list[FlowSegment] = []
        for duration, flow in steps:
            span = float(duration)
            if not math.isfinite(span) or span <= 0.0:
                raise ValueError(f"each duration must be finite and positive, got {duration!r}")
            segments.append(FlowSegment(cursor, cursor + span, flow))
            cursor += span
        return cls(tuple(segments))

    def mass_flow_kg_s(self, time_s: float) -> float:
        t = float(time_s)
        for segment in self.segments:
            if segment.t_start_s <= t < segment.t_end_s:
                return segment.mass_flow_kg_s
        if t == self.t_end_s:
            return self.segments[-1].mass_flow_kg_s
        raise ValueError(
            f"time {t!r} is outside the schedule window "
            f"[{self.t_start_s}, {self.t_end_s}]"
        )

    @property
    def interior_breakpoints_s(self) -> tuple[float, ...]:
        return tuple(segment.t_start_s for segment in self.segments[1:])

    def flow_on_segment(self, seg_start_s: float, seg_end_s: float) -> Callable[[float], float]:
        """Bind the single constant that applies across this segment.

        The midpoint unambiguously identifies the owning segment regardless of
        the half-open convention at either end, so the returned law is exactly
        constant over the whole closed segment -- including its endpoints.
        """
        value = self.mass_flow_kg_s(0.5 * (seg_start_s + seg_end_s))

        def flow(time_s: float) -> float:
            return value

        return flow


@dataclass(frozen=True)
class CallableOxidizerFlow(OxidizerFlowHistory):
    """A user-supplied ``m_dot_ox(t)`` callable.

    The callable is validated by sampling ``n_validation_samples`` points across
    the window at construction; every sample must be finite and non-negative.
    Sampling cannot prove the function is well behaved everywhere, so the
    integrator also re-validates on every right-hand-side evaluation.

    The callable is integrated as a **single smooth segment**.  A history with
    genuine discontinuities should use :class:`PiecewiseConstantOxidizerFlow`
    instead, so that the solver restarts at each jump.
    """

    function: Callable[[float], float]
    t_end_s: float
    t_start_s: float = 0.0
    n_validation_samples: int = 64

    def __post_init__(self) -> None:
        if not callable(self.function):
            raise TypeError("function must be callable")
        object.__setattr__(self, "t_start_s", _validate_time("t_start_s", self.t_start_s))
        object.__setattr__(self, "t_end_s", _validate_time("t_end_s", self.t_end_s))
        self._check_window()
        samples = int(self.n_validation_samples)
        if samples < 2:
            raise ValueError(f"n_validation_samples must be at least 2, got {samples!r}")
        for t in np.linspace(self.t_start_s, self.t_end_s, samples):
            _validate_flow(f"function({t})", self.function(float(t)))

    def mass_flow_kg_s(self, time_s: float) -> float:
        return _validate_flow(f"function({time_s})", self.function(float(time_s)))


# ---------------------------------------------------------------------------
# Geometry helper valid on the closed interval up to burnout
# ---------------------------------------------------------------------------


def remaining_fuel_mass_kg(geometry: GrainGeometry, port_radius_m: float | np.ndarray):
    """Remaining solid fuel, ``rho_f L pi (r_outer^2 - r_p^2)`` [kg].

    Milestone 1's :meth:`GrainGeometry.fuel_mass_kg` deliberately rejects
    ``D_p >= D_o`` because a fully burnt-through grain is a degenerate state for
    instantaneous bookkeeping.  A transient integration must be able to *reach*
    that boundary exactly, so this function evaluates the same geometry on the
    closed interval ``r_p <= r_outer``, returning exactly zero at burnout.  It is
    written independently of the M1 routine, which also makes it a useful
    cross-check (see ``tests/test_transient.py``).
    """
    r = np.asarray(port_radius_m, dtype=float)
    r_outer = geometry.outer_diameter_m / 2.0
    if np.any(~np.isfinite(r)):
        raise ValueError("port_radius_m must be finite")
    if np.any(r <= 0.0):
        raise ValueError("port_radius_m must be strictly positive")
    if np.any(r > r_outer * (1.0 + 1.0e-12)):
        raise ValueError(
            f"port_radius_m must not exceed the outer radius {r_outer!r}; "
            "the grain would be burnt through"
        )
    mass = geometry.fuel_density_kg_m3 * geometry.length_m * math.pi * (r_outer**2 - r**2)
    if np.isscalar(port_radius_m) or np.ndim(port_radius_m) == 0:
        return float(mass)
    return mass


# ---------------------------------------------------------------------------
# Closed-form constant-flow solution
# ---------------------------------------------------------------------------


def analytic_constant_flow_speed_coefficient(
    law: PowerLawRegressionLaw, oxidizer_mass_flow_kg_s: float
) -> float:
    """``K = a_SI (m_dot_ox / pi)^n`` for ``dr/dt = K r^(-2n)`` [SI].

    Substituting ``G_ox = m_dot_ox / (pi r^2)`` into ``r_dot = a_SI G_ox^n``:

        dr/dt = a_SI (m_dot_ox / (pi r^2))^n = a_SI (m_dot_ox/pi)^n r^(-2n)
    """
    m_dot = _validate_flow("oxidizer_mass_flow_kg_s", oxidizer_mass_flow_kg_s)
    return law.coefficient_si * (m_dot / math.pi) ** law.exponent


def analytic_constant_flow_radius(
    law: PowerLawRegressionLaw,
    initial_port_radius_m: float,
    oxidizer_mass_flow_kg_s: float,
    elapsed_time_s,
):
    """Closed-form port radius under constant prescribed oxidizer flow [m].

    Separating ``dr/dt = K r^(-2n)`` gives ``r^(2n) dr = K dt``, hence

        r(t)^(2n+1) = r_0^(2n+1) + (2n+1) K t
        r(t)        = [ r_0^(2n+1) + (2n+1) K t ]^(1/(2n+1))

    For ``n = 0.5`` this reduces to ``r^2 = r_0^2 + 2Kt``, i.e. ``D^2`` linear in
    ``t`` -- the same form as the exact solution published by Karabeyoglu,
    Cantwell & Zilliac (2007) for a flux exponent of 0.5 with no length
    dependence.  This function does **not** know about burnout; it is the
    unbounded solution of the state equation.
    """
    r0 = float(initial_port_radius_m)
    if not math.isfinite(r0) or r0 <= 0.0:
        raise ValueError(f"initial_port_radius_m must be finite and positive, got {r0!r}")
    t = np.asarray(elapsed_time_s, dtype=float)
    if np.any(~np.isfinite(t)):
        raise ValueError("elapsed_time_s must be finite")
    if np.any(t < 0.0):
        raise ValueError("elapsed_time_s must be non-negative")
    k = analytic_constant_flow_speed_coefficient(law, oxidizer_mass_flow_kg_s)
    exponent = 2.0 * law.exponent + 1.0
    radius = (r0**exponent + exponent * k * t) ** (1.0 / exponent)
    if np.isscalar(elapsed_time_s) or np.ndim(elapsed_time_s) == 0:
        return float(radius)
    return radius


def analytic_constant_flow_burnout_time(
    law: PowerLawRegressionLaw,
    initial_port_radius_m: float,
    outer_radius_m: float,
    oxidizer_mass_flow_kg_s: float,
) -> float:
    """Closed-form burnout time under constant prescribed oxidizer flow [s].

    Inverting the radius solution at ``r = r_outer``:

        t_burn = [ r_outer^(2n+1) - r_0^(2n+1) ] / [ (2n+1) K ]

    Returns ``inf`` when the prescribed flow is zero, because the port then never
    regresses and burnout is never reached.
    """
    r0 = float(initial_port_radius_m)
    r_outer = float(outer_radius_m)
    if not math.isfinite(r0) or r0 <= 0.0:
        raise ValueError(f"initial_port_radius_m must be finite and positive, got {r0!r}")
    if not math.isfinite(r_outer) or r_outer <= r0:
        raise ValueError(
            f"outer_radius_m must be finite and strictly greater than the initial radius, "
            f"got {r_outer!r} with r_0 = {r0!r}"
        )
    k = analytic_constant_flow_speed_coefficient(law, oxidizer_mass_flow_kg_s)
    if k == 0.0:
        return math.inf
    exponent = 2.0 * law.exponent + 1.0
    return (r_outer**exponent - r0**exponent) / (exponent * k)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TransientResult:
    """Time histories from a prescribed-flow transient integration.

    All arrays share the reporting grid ``time_s``.  At an interior schedule
    breakpoint the grid carries two nearly coincident samples so that a
    discontinuous prescribed flow is represented by its left and right limits
    rather than by a spurious ramp.

    ``mixture_ratio`` is ``NaN`` wherever the fuel mass flow is zero: the mixture
    ratio is genuinely undefined there and is not fabricated.
    """

    time_s: np.ndarray
    port_radius_m: np.ndarray
    port_diameter_m: np.ndarray
    port_area_m2: np.ndarray
    web_thickness_m: np.ndarray
    remaining_fuel_mass_kg: np.ndarray
    oxidizer_mass_flow_kg_s: np.ndarray
    oxidizer_mass_flux_si: np.ndarray
    regression_rate_m_s: np.ndarray
    regression_rate_mm_s: np.ndarray
    burning_area_m2: np.ndarray
    fuel_mass_flow_kg_s: np.ndarray
    mixture_ratio: np.ndarray
    cumulative_oxidizer_mass_kg: np.ndarray
    cumulative_fuel_mass_kg: np.ndarray
    flux_status: tuple[FluxRangeStatus, ...]
    termination_reason: TerminationReason
    burnout_time_s: float | None
    geometry: GrainGeometry
    law: PowerLawRegressionLaw
    initial_port_radius_m: float

    # -- mass closure ------------------------------------------------------

    @property
    def geometric_fuel_consumed_kg(self) -> np.ndarray:
        """Fuel consumed inferred from the port geometry alone [kg]."""
        return (
            self.geometry.fuel_density_kg_m3
            * self.geometry.length_m
            * math.pi
            * (self.port_radius_m**2 - self.initial_port_radius_m**2)
        )

    @property
    def fuel_mass_closure_residual_kg(self) -> np.ndarray:
        """Integrated fuel mass minus geometric fuel loss [kg].

        Zero in exact arithmetic, so this is a pure measure of integration error.
        """
        return self.cumulative_fuel_mass_kg - self.geometric_fuel_consumed_kg

    @property
    def max_absolute_fuel_closure_residual_kg(self) -> float:
        return float(np.max(np.abs(self.fuel_mass_closure_residual_kg)))

    @property
    def max_relative_fuel_closure_residual(self) -> float:
        """Worst residual relative to the fuel consumed at that time [-]."""
        consumed = self.geometric_fuel_consumed_kg
        mask = consumed > 0.0
        if not np.any(mask):
            return 0.0
        return float(
            np.max(np.abs(self.fuel_mass_closure_residual_kg[mask] / consumed[mask]))
        )

    # -- validity bookkeeping ---------------------------------------------

    def time_fraction_with_status(self, status: FluxRangeStatus) -> float:
        """Fraction of the simulated *duration* spent in a flux-range status.

        Measured by trapezoidal weighting of the reporting grid, so it reflects
        elapsed time rather than the number of samples.
        """
        total = self.time_s[-1] - self.time_s[0]
        if total <= 0.0:
            return 0.0
        flags = np.array([s is status for s in self.flux_status], dtype=float)
        return float(np.trapezoid(flags, self.time_s) / total)

    @property
    def final_index(self) -> int:
        return int(self.time_s.size - 1)


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------


def _classify_flux(law: PowerLawRegressionLaw, flux_si: float, flow_kg_s: float) -> FluxRangeStatus:
    if flow_kg_s == 0.0:
        return FluxRangeStatus.ZERO_OXIDIZER_FLOW
    if law.valid_flux_range_si is None:
        return FluxRangeStatus.NO_SOURCE_RANGE_DECLARED
    low, high = law.valid_flux_range_si
    if flux_si < low:
        return FluxRangeStatus.BELOW_SOURCE_FLUX_RANGE
    if flux_si > high:
        return FluxRangeStatus.ABOVE_SOURCE_FLUX_RANGE
    return FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE


def simulate_transient(
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    flow_history: OxidizerFlowHistory,
    *,
    initial_port_diameter_m: float | None = None,
    t_end_s: float | None = None,
    n_report: int = 401,
    rtol: float = 1.0e-10,
    atol: float = 1.0e-14,
    max_step: float = np.inf,
    method: str = "DOP853",
) -> TransientResult:
    """Integrate port growth under a prescribed oxidizer mass-flow history.

    The integration is restarted at every interior breakpoint of the flow
    history, so a piecewise-constant schedule is integrated exactly segment by
    segment rather than having a step discontinuity smeared across it.  A
    terminal event stops the integration the moment the port reaches the outer
    grain radius, so the solution never enters negative fuel thickness.

    Parameters
    ----------
    geometry, law:
        The Milestone 1 grain and regression correlation, used unchanged.
    flow_history:
        Prescribed ``m_dot_ox(t)``.
    initial_port_diameter_m:
        Defaults to the grain's initial port diameter.
    t_end_s:
        Optional earlier end time; must lie inside the flow history's window.
    n_report:
        Number of evenly spaced reporting samples.  Segment breakpoints and the
        termination time are always added on top of these.
    rtol, atol, max_step, method:
        Passed through to :func:`scipy.integrate.solve_ivp`.
    """
    if n_report < 2:
        raise ValueError(f"n_report must be at least 2, got {n_report!r}")

    diameter0 = (
        geometry.initial_port_diameter_m
        if initial_port_diameter_m is None
        else float(initial_port_diameter_m)
    )
    # Delegate positivity / inside-the-grain validation to the M1 geometry.
    geometry.port_area_m2(diameter0)
    r0 = diameter0 / 2.0
    r_outer = geometry.outer_diameter_m / 2.0

    t_start = flow_history.t_start_s
    t_final_requested = flow_history.t_end_s if t_end_s is None else _validate_time(
        "t_end_s", t_end_s
    )
    if not (t_start < t_final_requested <= flow_history.t_end_s):
        raise ValueError(
            f"t_end_s must lie in ({t_start}, {flow_history.t_end_s}], got {t_final_requested!r}"
        )

    rho_f = geometry.fuel_density_kg_m3
    length = geometry.length_m

    def make_rhs(segment_flow: Callable[[float], float]):
        def rhs(t: float, y: np.ndarray) -> list[float]:
            radius = float(y[0])
            flow = segment_flow(t)
            area = math.pi * radius * radius
            flux = flow / area
            rate = law.regression_rate_si(flux)
            return [rate, flow, rho_f * 2.0 * math.pi * radius * length * rate]

        return rhs

    def burnout_event(t: float, y: np.ndarray) -> float:
        return float(y[0]) - r_outer

    burnout_event.terminal = True
    burnout_event.direction = 1

    # Integration boundaries: window ends plus interior flow discontinuities.
    boundaries = [t_start, *flow_history.interior_breakpoints_s, t_final_requested]
    boundaries = sorted({b for b in boundaries if t_start <= b <= t_final_requested})

    state = np.array([r0, 0.0, 0.0], dtype=float)
    solutions: list[tuple[float, float, object]] = []
    termination = TerminationReason.COMPLETED_TIME_WINDOW
    burnout_time: float | None = None
    t_final = t_final_requested

    for seg_start, seg_end in itertools.pairwise(boundaries):
        solution = solve_ivp(
            make_rhs(flow_history.flow_on_segment(seg_start, seg_end)),
            (seg_start, seg_end),
            state,
            method=method,
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=True,
            events=burnout_event,
        )
        if not solution.success:
            raise RuntimeError(
                f"integration failed on [{seg_start}, {seg_end}]: {solution.message}"
            )
        solutions.append((seg_start, float(solution.t[-1]), solution.sol))
        state = solution.y[:, -1].copy()
        if solution.t_events[0].size > 0:
            termination = TerminationReason.FUEL_DEPLETED
            burnout_time = float(solution.t_events[0][0])
            t_final = burnout_time
            break
        t_final = float(solution.t[-1])

    # -- reporting grid ----------------------------------------------------
    grid = set(np.linspace(t_start, t_final, n_report).tolist())
    grid.add(t_start)
    grid.add(t_final)
    for breakpoint_s in flow_history.interior_breakpoints_s:
        if t_start < breakpoint_s < t_final:
            # Two nearly coincident samples capture the left and right limits of
            # the discontinuous prescribed flow.
            grid.add(float(np.nextafter(breakpoint_s, -np.inf)))
            grid.add(float(breakpoint_s))
    times = np.array(sorted(grid), dtype=float)

    def evaluate(t: float) -> np.ndarray:
        for seg_start, seg_end, dense in solutions:
            if seg_start <= t <= seg_end:
                return np.asarray(dense(t), dtype=float)
        return np.asarray(solutions[-1][2](t), dtype=float)

    raw = np.array([evaluate(float(t)) for t in times])
    radius = raw[:, 0]
    # The event root is found to solver precision; clip only the last sample's
    # round-off so that r <= r_outer holds exactly at burnout.
    if termination is TerminationReason.FUEL_DEPLETED:
        radius = np.minimum(radius, r_outer)
    cumulative_ox = raw[:, 1]
    cumulative_fuel = raw[:, 2]

    flows = np.array([flow_history.mass_flow_kg_s(float(t)) for t in times])
    areas = math.pi * radius**2
    fluxes = flows / areas
    rates = law.regression_rate_si(fluxes)
    burning_areas = 2.0 * math.pi * radius * length
    fuel_flows = rho_f * burning_areas * rates

    with np.errstate(divide="ignore", invalid="ignore"):
        mixture = np.where(fuel_flows > 0.0, flows / fuel_flows, np.nan)

    statuses = tuple(
        _classify_flux(law, float(flux), float(flow))
        for flux, flow in zip(fluxes, flows, strict=True)
    )

    return TransientResult(
        time_s=times,
        port_radius_m=radius,
        port_diameter_m=2.0 * radius,
        port_area_m2=areas,
        web_thickness_m=r_outer - radius,
        remaining_fuel_mass_kg=remaining_fuel_mass_kg(geometry, radius),
        oxidizer_mass_flow_kg_s=flows,
        oxidizer_mass_flux_si=fluxes,
        regression_rate_m_s=rates,
        regression_rate_mm_s=rates * 1.0e3,
        burning_area_m2=burning_areas,
        fuel_mass_flow_kg_s=fuel_flows,
        mixture_ratio=mixture,
        cumulative_oxidizer_mass_kg=cumulative_ox,
        cumulative_fuel_mass_kg=cumulative_fuel,
        flux_status=statuses,
        termination_reason=termination,
        burnout_time_s=burnout_time,
        geometry=geometry,
        law=law,
        initial_port_radius_m=r0,
    )
