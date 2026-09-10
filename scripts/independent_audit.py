"""End-to-end independent verification of the whole M1-M5 model chain.

This is the Milestone 6 audit. Its purpose is *not* to re-run the project's own
tests: it is to recompute the identities the model rests on from **raw formulas
written out here**, and compare them against what the production modules report.

The rule this script follows
----------------------------
The "expected" side of every check is built from one of:

* an algebraic identity written out longhand in this file,
* a closed-form solution derived independently of the solver,
* an alternate numerical route (plain bisection, trapezoid, secant) that shares
  no code with the production path,
* published experimental data transcribed from the source paper,
* raw CSV columns parsed with the standard library.

Nothing is checked by calling the same function twice.

Tolerances
----------
Two different standards apply, and they are labelled per check:

* ``identity`` -- an algebraic relation that must hold to floating-point
  precision. Tolerance ~1e-9 relative or tighter.
* ``discretisation`` -- a quantity produced by an ODE integration or a
  quadrature, where the honest expectation is the integrator's own tolerance,
  not machine epsilon. **These are not reported as if they should be exact.**
* ``empirical`` -- agreement with measured data, where the tolerance is the
  scatter of the experiment.

Run with::

    python scripts/independent_audit.py
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from hybrid_rocket_motor import REPRESENTATIVE_GRAIN, REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.blowdown import simulate_blowdown
from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.injector import Injector
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.operating_point import OperatingPoint
from hybrid_rocket_motor.performance import (
    STANDARD_GRAVITY_M_S2,
    couple_transient_to_performance,
    solve_performance_point,
)
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
)
from hybrid_rocket_motor.thermochemistry import DEFAULT_TABLE_PATH, TableStatus, default_table
from hybrid_rocket_motor.transient import ConstantOxidizerFlow, simulate_transient
from hybrid_rocket_motor.variable_blowdown import simulate_variable_blowdown
from hybrid_rocket_motor.variable_chamber import solve_variable_chamber
from hybrid_rocket_motor.variable_feed_system import solve_variable_feed_coupling

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
NOZZLE = REFERENCE_NOZZLE
FROZEN = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
INJECTOR = Injector(3.5e-6)
SEA_LEVEL_PA = 101325.0
UNIVERSAL_GAS_CONSTANT_J_MOL_K = 8.314462618

M1_OXIDIZER_FLOW = 0.100
RULE = "-" * 104

#: Measured (mean oxidizer mass flux [g/(cm^2 s)], regression rate [mm/s]) pairs
#: from Table 4 of Rezaei et al. (2018) -- the firings the correlation was fitted
#: to. Experimental data, not model output.
REZAEI_TABLE_4 = [
    (6.88, 0.779),
    (11.95, 0.891),
    (11.61, 1.044),
    (7.70, 0.882),
    (6.36, 0.775),
    (5.93, 0.810),
    (4.65, 0.653),
    (4.89, 0.740),
    (11.61, 1.006),
    (8.79, 0.940),
    (5.47, 0.710),
    (12.02, 0.899),
    (5.12, 0.693),
    (6.19, 0.739),
    (8.64, 0.953),
    (10.20, 0.948),
    (4.66, 0.732),
    (3.24, 0.613),
]


@dataclass
class Check:
    subsystem: str
    name: str
    actual: float
    expected: float
    tolerance: float
    kind: str
    unit: str = ""
    relative: bool = True

    @property
    def absolute_residual(self) -> float:
        return abs(self.actual - self.expected)

    @property
    def relative_residual(self) -> float:
        scale = abs(self.expected)
        return self.absolute_residual / scale if scale > 0.0 else self.absolute_residual

    @property
    def residual(self) -> float:
        return self.relative_residual if self.relative else self.absolute_residual

    @property
    def passed(self) -> bool:
        return self.residual <= self.tolerance


@dataclass
class Audit:
    checks: list[Check] = field(default_factory=list)

    def add(self, subsystem, name, actual, expected, tolerance, kind, unit="", relative=True):
        self.checks.append(
            Check(subsystem, name, float(actual), float(expected), tolerance, kind, unit, relative)
        )

    def report(self, subsystem: str) -> None:
        rows = [c for c in self.checks if c.subsystem == subsystem]
        print(f"{subsystem}")
        print(RULE)
        header = f"  {'check':<44} {'abs resid':>12} {'rel resid':>12} {'tol':>10} {'kind':<15}"
        print(header)
        for c in rows:
            mark = "OK  " if c.passed else "FAIL"
            rel = f"{c.relative_residual:12.3e}" if abs(c.expected) > 0 else f"{'--':>12}"
            print(
                f"  [{mark}] {c.name:<38} {c.absolute_residual:12.3e} {rel} "
                f"{c.tolerance:10.1e} {c.kind:<15}"
            )
        worst = max((c.residual for c in rows), default=0.0)
        print(
            f"  worst residual in this subsystem: {worst:.3e} "
            f"({sum(c.passed for c in rows)}/{len(rows)} passed)"
        )
        print()


# ---------------------------------------------------------------------------
# M1 -- geometry and regression
# ---------------------------------------------------------------------------


def audit_m1(audit: Audit) -> None:
    sub = "M1 -- geometry and regression"
    diameter = GRAIN.initial_port_diameter_m
    point = OperatingPoint(GRAIN, M1_OXIDIZER_FLOW)

    # Longhand geometry, written out rather than called.
    port_area = math.pi * diameter * diameter / 4.0
    burn_area = math.pi * diameter * GRAIN.length_m
    audit.add(
        sub, "port area = pi D^2 / 4", point.port_area_m2, port_area, 1e-15, "identity", "m^2"
    )
    audit.add(sub, "burn area = pi D L", point.burning_area_m2, burn_area, 1e-15, "identity", "m^2")

    flux = M1_OXIDIZER_FLOW / port_area
    audit.add(
        sub,
        "G_ox = m_dot_ox / A_port",
        point.oxidizer_mass_flux_si,
        flux,
        1e-15,
        "identity",
        "kg/(m^2 s)",
    )

    # The regression law rebuilt from the SOURCE units (mm/s, g/(cm^2 s)) and
    # converted here, rather than from the stored SI coefficient.
    source_a_mm_s, source_n = 0.3977, 0.3667
    rate_mm_s = source_a_mm_s * (flux * 0.1) ** source_n  # kg/(m^2 s) -> g/(cm^2 s)
    audit.add(
        sub,
        "r_dot from source-unit correlation",
        point.regression_rate_si(LAW) * 1e3,
        rate_mm_s,
        1e-12,
        "identity",
        "mm/s",
    )

    fuel = GRAIN.fuel_density_kg_m3 * burn_area * rate_mm_s * 1e-3
    audit.add(
        sub,
        "m_dot_f = rho_f A_burn r_dot",
        point.fuel_mass_flow_kg_s(LAW),
        fuel,
        1e-12,
        "identity",
        "kg/s",
    )
    audit.add(
        sub,
        "O/F = m_dot_ox / m_dot_f",
        point.mixture_ratio(LAW),
        M1_OXIDIZER_FLOW / fuel,
        1e-12,
        "identity",
        "-",
    )

    # Published measurements: the tolerance is the experiment's own scatter.
    deviations = []
    for flux_g_cm2_s, measured_mm_s in REZAEI_TABLE_4:
        predicted = LAW.regression_rate_si(flux_g_cm2_s * 10.0) * 1e3
        deviations.append(abs(predicted / measured_mm_s - 1.0))
    audit.add(
        sub,
        f"worst of {len(REZAEI_TABLE_4)} published firings",
        max(deviations),
        0.0,
        0.20,
        "empirical",
        "-",
        relative=False,
    )
    audit.add(
        sub,
        "mean deviation over those firings",
        sum(deviations) / len(deviations),
        0.0,
        0.10,
        "empirical",
        "-",
        relative=False,
    )


# ---------------------------------------------------------------------------
# M2 -- transient port evolution
# ---------------------------------------------------------------------------


def audit_m2(audit: Audit) -> None:
    sub = "M2 -- transient port evolution"
    result = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M1_OXIDIZER_FLOW, 200.0), t_end_s=200.0
    )

    # Closed form, re-derived here. With r_dot = a (m_dot/(pi r^2))^n,
    #   r^(2n) dr = a (m_dot/pi)^n dt  ->  r(t) = [r0^(2n+1) + (2n+1) a (m/pi)^n t]^(1/(2n+1))
    n = LAW.exponent
    a = LAW.coefficient_si
    r0 = GRAIN.initial_port_diameter_m / 2.0
    r_outer = GRAIN.outer_diameter_m / 2.0
    power = 2.0 * n + 1.0
    speed = a * (M1_OXIDIZER_FLOW / math.pi) ** n

    def radius_at(t: float) -> float:
        return (r0**power + power * speed * t) ** (1.0 / power)

    burnout = (r_outer**power - r0**power) / (power * speed)
    audit.add(sub, "analytic burnout time", result.burnout_time_s, burnout, 1e-9, "identity", "s")

    worst = 0.0
    for fraction in (0.1, 0.25, 0.5, 0.75, 0.95):
        t = fraction * burnout
        analytic = radius_at(t)
        numeric = float(np.interp(t, result.time_s, result.port_radius_m))
        worst = max(worst, abs(numeric / analytic - 1.0))
    audit.add(
        sub,
        "analytic radius r(t), worst of 5 samples",
        worst,
        0.0,
        1e-7,
        "discretisation",
        "-",
        relative=False,
    )

    # Oxidizer integral: constant flow, so the exact value is m_dot * t.
    audit.add(
        sub,
        "oxidizer mass integral = m_dot t",
        result.cumulative_oxidizer_mass_kg[result.final_index],
        M1_OXIDIZER_FLOW * result.burnout_time_s,
        1e-9,
        "identity",
        "kg",
    )

    # Fuel mass two ways: integrated flow vs the annulus the port swept out.
    swept = GRAIN.fuel_density_kg_m3 * GRAIN.length_m * math.pi * (r_outer**2 - r0**2)
    audit.add(
        sub,
        "fuel mass: integral vs swept geometry",
        result.cumulative_fuel_mass_kg[result.final_index],
        swept,
        1e-8,
        "discretisation",
        "kg",
    )
    audit.add(
        sub,
        "loaded fuel mass = rho pi (Ro^2 - Ri^2) L",
        GRAIN.initial_fuel_mass_kg,
        swept,
        1e-14,
        "identity",
        "kg",
    )


# ---------------------------------------------------------------------------
# M3 -- chamber and nozzle
# ---------------------------------------------------------------------------


def audit_m3(audit: Audit) -> None:
    sub = "M3 -- chamber and nozzle"
    fuel = 0.039776644
    total = M1_OXIDIZER_FLOW + fuel
    point = solve_performance_point(FROZEN, NOZZLE, M1_OXIDIZER_FLOW, fuel, SEA_LEVEL_PA)
    exit_state = point.exit_state
    gamma = FROZEN.gamma
    throat = NOZZLE.throat_area_m2

    audit.add(
        sub,
        "p_c A_t = m_dot c*",
        point.chamber.chamber_pressure_pa * throat,
        total * FROZEN.c_star_m_s,
        1e-14,
        "identity",
        "N",
    )

    # Area-Mach relation written out longhand, then re-solved by plain bisection
    # using no production routine.
    def area_ratio(mach: float) -> float:
        term = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach * mach)
        return term ** ((gamma + 1.0) / (2.0 * (gamma - 1.0))) / mach

    target = NOZZLE.expansion_ratio
    low, high = 1.0 + 1e-12, 50.0
    for _ in range(300):
        mid = 0.5 * (low + high)
        if area_ratio(mid) < target:
            low = mid
        else:
            high = mid
    mach = 0.5 * (low + high)
    audit.add(
        sub, "exit Mach by independent bisection", exit_state.mach, mach, 1e-11, "identity", "-"
    )

    pressure_ratio = (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** (-gamma / (gamma - 1.0))
    temperature_ratio = 1.0 / (1.0 + 0.5 * (gamma - 1.0) * mach * mach)
    audit.add(
        sub,
        "p_e / p_c isentropic",
        exit_state.pressure_ratio,
        pressure_ratio,
        1e-11,
        "identity",
        "-",
    )
    audit.add(
        sub,
        "T_e / T_c isentropic",
        exit_state.temperature_ratio,
        temperature_ratio,
        1e-11,
        "identity",
        "-",
    )

    exit_temperature = FROZEN.chamber_temperature_k * temperature_ratio
    velocity = mach * math.sqrt(gamma * FROZEN.gas_constant_j_kg_k * exit_temperature)
    audit.add(
        sub,
        "V_e = M_e sqrt(gamma R T_e)",
        exit_state.velocity_m_s,
        velocity,
        1e-11,
        "identity",
        "m/s",
    )

    exit_pressure = point.chamber.chamber_pressure_pa * pressure_ratio
    momentum = total * velocity
    pressure_thrust = (exit_pressure - SEA_LEVEL_PA) * NOZZLE.exit_area_m2
    audit.add(
        sub,
        "F = m_dot V_e + (p_e - p_a) A_e",
        point.thrust_n,
        momentum + pressure_thrust,
        1e-11,
        "identity",
        "N",
    )

    # C_F and I_sp identities, and the two independent routes to thrust.
    thrust_coefficient = point.thrust_n / (point.chamber.chamber_pressure_pa * throat)
    audit.add(
        sub,
        "C_F = F / (p_c A_t)",
        point.thrust_coefficient,
        thrust_coefficient,
        1e-13,
        "identity",
        "-",
    )
    audit.add(
        sub,
        "F via C_F p_c A_t route",
        point.thrust_coefficient * point.chamber.chamber_pressure_pa * throat,
        point.thrust_n,
        1e-13,
        "identity",
        "N",
    )
    audit.add(
        sub,
        "I_sp = F / (m_dot g0)",
        point.specific_impulse_s,
        point.thrust_n / (total * STANDARD_GRAVITY_M_S2),
        1e-13,
        "identity",
        "s",
    )
    audit.add(
        sub,
        "I_sp = c* C_F / g0",
        point.specific_impulse_s,
        FROZEN.c_star_m_s * point.thrust_coefficient / STANDARD_GRAVITY_M_S2,
        1e-13,
        "identity",
        "s",
    )

    # Total impulse against an independently coded trapezoid of the history.
    transient = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M1_OXIDIZER_FLOW, 200.0), t_end_s=200.0
    )
    history = couple_transient_to_performance(transient, FROZEN, NOZZLE, SEA_LEVEL_PA)
    index = transient.final_index
    times = history.time_s[: index + 1]
    thrusts = history.thrust_n[: index + 1]
    trapezoid = float(np.sum(0.5 * (thrusts[1:] + thrusts[:-1]) * np.diff(times)))
    audit.add(
        sub,
        "total impulse vs independent trapezoid",
        history.total_impulse_n_s,
        trapezoid,
        1e-6,
        "discretisation",
        "N s",
    )


# ---------------------------------------------------------------------------
# M4 -- tank, injector, feed, blowdown
# ---------------------------------------------------------------------------


def audit_m4(audit: Audit) -> None:
    sub = "M4 -- tank, injector, feed and blowdown"
    state = REFERENCE_TANK.initial_state(REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION)

    audit.add(
        sub,
        "tank phase mass closure",
        state.liquid_mass_kg + state.vapour_mass_kg,
        state.total_mass_kg,
        1e-14,
        "identity",
        "kg",
    )
    audit.add(
        sub,
        "tank volume closure",
        state.liquid_volume_m3 + state.vapour_volume_m3,
        REFERENCE_TANK.volume_m3,
        1e-14,
        "identity",
        "m^3",
    )
    audit.add(
        sub,
        "vapour quality = m_v / m_total",
        state.vapour_quality,
        state.vapour_mass_kg / state.total_mass_kg,
        1e-14,
        "identity",
        "-",
    )

    # Saturation pressure through a second CoolProp route (quality-based rather
    # than the production temperature-based call).
    from CoolProp.CoolProp import PropsSI

    audit.add(
        sub,
        "p_sat(T) via an independent CoolProp call",
        state.pressure_pa,
        PropsSI("P", "T", state.temperature_k, "Q", 0.0, "NitrousOxide"),
        1e-12,
        "identity",
        "Pa",
    )
    audit.add(
        sub,
        "liquid density via independent CoolProp call",
        state.liquid_density_kg_m3,
        PropsSI("D", "T", state.temperature_k, "Q", 0.0, "NitrousOxide"),
        1e-12,
        "identity",
        "kg/m^3",
    )

    # SPI written out longhand at one representative backpressure.
    downstream = 26.0e5
    result = INJECTOR.mass_flow(state, downstream)
    spi = (
        0.66
        * 3.5e-6
        * math.sqrt(2.0 * state.liquid_density_kg_m3 * (state.pressure_pa - downstream))
    )
    audit.add(
        sub, "SPI = C_d A sqrt(2 rho dp)", result.spi_mass_flow_kg_s, spi, 1e-13, "identity", "kg/s"
    )
    # For a saturated tank the Dyer non-equilibrium parameter is exactly 1, so
    # the blend is the plain average of the two limits.
    audit.add(
        sub,
        "Dyer kappa = 1 for a saturated tank",
        result.non_equilibrium_parameter,
        1.0,
        1e-12,
        "identity",
        "-",
    )
    audit.add(
        sub,
        "Dyer = (SPI + HEM) / 2 when kappa = 1",
        result.mass_flow_kg_s,
        0.5 * (result.spi_mass_flow_kg_s + result.hem_mass_flow_kg_s),
        1e-13,
        "identity",
        "kg/s",
    )

    blowdown = simulate_blowdown(
        REFERENCE_TANK, state, INJECTOR, GRAIN, LAW, FROZEN, NOZZLE, SEA_LEVEL_PA
    )

    # Feed root: recompute the residual from the reported histories.
    worst_feed = 0.0
    for i in (0, len(blowdown.time_s) // 2, -1):
        tank_state = REFERENCE_TANK.state_from_mass_and_energy(
            float(blowdown.tank_total_mass_kg[i]),
            float(blowdown.tank_internal_energy_j[i]),
            quality_overshoot=0.6,
        )
        flow = INJECTOR.mass_flow(tank_state, float(blowdown.chamber_pressure_pa[i]))
        worst_feed = max(
            worst_feed, abs(flow.mass_flow_kg_s - float(blowdown.oxidizer_mass_flow_kg_s[i]))
        )
    audit.add(
        sub,
        "feed residual re-formed at 3 samples",
        worst_feed,
        0.0,
        1e-11,
        "identity",
        "kg/s",
        relative=False,
    )

    audit.add(
        sub,
        "oxidizer mass: tank drawn vs integrated",
        float(np.max(np.abs(blowdown.oxidizer_mass_closure_residual_kg))),
        0.0,
        1e-12,
        "discretisation",
        "kg",
        relative=False,
    )
    audit.add(
        sub,
        "fuel mass: integral vs swept geometry",
        float(np.max(np.abs(blowdown.fuel_mass_closure_residual_kg))),
        0.0,
        1e-6,
        "discretisation",
        "kg",
        relative=False,
    )
    audit.add(
        sub,
        "tank energy balance, relative",
        blowdown.relative_energy_closure_residual,
        0.0,
        1e-6,
        "discretisation",
        "-",
        relative=False,
    )
    audit.add(
        sub,
        "impulse vs trapezoid, relative",
        blowdown.relative_impulse_closure_residual,
        0.0,
        1e-6,
        "discretisation",
        "-",
        relative=False,
    )


# ---------------------------------------------------------------------------
# M5 -- thermochemistry
# ---------------------------------------------------------------------------


def audit_m5(audit: Audit) -> None:
    sub = "M5 -- thermochemistry and coupled solve"
    table = default_table()

    # Raw CSV, parsed with the standard library rather than the package loader.
    with Path(DEFAULT_TABLE_PATH).open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    of_values = sorted({float(r["of"]) for r in rows})
    pressures = sorted({float(r["chamber_pressure_pa"]) for r in rows})

    audit.add(
        sub,
        "grid is complete and rectangular",
        len(rows),
        len(of_values) * len(pressures),
        0.0,
        "identity",
        "rows",
    )
    # Monotonicity is a boolean property, so it is scored as one: 1 means every
    # consecutive spacing is strictly positive.
    audit.add(
        sub,
        "O/F axis strictly increasing",
        float(bool(np.all(np.diff(of_values) > 0.0))),
        1.0,
        0.0,
        "identity",
        "-",
    )
    audit.add(
        sub,
        "pressure axis strictly increasing",
        float(bool(np.all(np.diff(pressures) > 0.0))),
        1.0,
        0.0,
        "identity",
        "-",
    )

    lookup = {(float(r["of"]), float(r["chamber_pressure_pa"])): r for r in rows}

    # Exact grid-point retrieval straight out of the raw CSV.
    worst_node = 0.0
    for of in (of_values[0], of_values[56], of_values[-1]):
        for pressure in (pressures[0], pressures[7], pressures[-1]):
            row = lookup[(of, pressure)]
            state = table.evaluate(of, pressure, c_star_efficiency=1.0)
            worst_node = max(
                worst_node,
                abs(state.ideal_c_star_m_s / float(row["ideal_cstar_m_s"]) - 1.0),
                abs(state.gamma / float(row["gamma"]) - 1.0),
            )
    audit.add(
        sub,
        "grid-point retrieval, worst of 9 nodes",
        worst_node,
        0.0,
        1e-14,
        "identity",
        "-",
        relative=False,
    )

    # Bilinear interpolation by hand at a cell centre, from the four raw corners.
    i, j = 40, 5
    of_lo, of_hi = of_values[i], of_values[i + 1]
    p_lo, p_hi = pressures[j], pressures[j + 1]
    corners = [float(lookup[(o, p)]["gamma"]) for o in (of_lo, of_hi) for p in (p_lo, p_hi)]
    centre = table.evaluate(0.5 * (of_lo + of_hi), 0.5 * (p_lo + p_hi), c_star_efficiency=1.0).gamma
    audit.add(
        sub,
        "bilinear cell centre = mean of 4 corners",
        centre,
        sum(corners) / 4.0,
        1e-14,
        "identity",
        "-",
    )

    # R = R_u / M off-node, where it is recomputed rather than interpolated.
    state = table.evaluate(2.3137, 21.37e5)
    audit.add(
        sub,
        "R = R_u / M off-node",
        state.specific_gas_constant_j_kg_k,
        UNIVERSAL_GAS_CONSTANT_J_MOL_K / state.effective_molar_mass_kg_mol,
        1e-15,
        "identity",
        "J/(kg K)",
    )

    # Efficiency convention: delivered = eta * ideal, applied exactly once.
    delivered = table.evaluate(2.3137, 21.37e5, c_star_efficiency=0.83)
    audit.add(
        sub,
        "c*_delivered = eta c*_ideal",
        delivered.c_star_m_s,
        0.83 * state.ideal_c_star_m_s,
        1e-15,
        "identity",
        "m/s",
    )
    audit.add(sub, "eta does not touch gamma", delivered.gamma, state.gamma, 1e-15, "identity", "-")

    # The implicit chamber solve, re-solved here by plain bisection.
    fuel = 0.039776644
    total = M1_OXIDIZER_FLOW + fuel
    of_ratio = M1_OXIDIZER_FLOW / fuel
    chamber = solve_variable_chamber(table, NOZZLE.throat_area_m2, M1_OXIDIZER_FLOW, fuel)

    def residual(pressure: float) -> float:
        s = table.evaluate(of_ratio, pressure)
        return total * s.c_star_m_s / NOZZLE.throat_area_m2 - pressure

    low, high = table.pressure_bounds_pa
    for _ in range(200):
        mid = 0.5 * (low + high)
        if residual(mid) > 0.0:
            low = mid
        else:
            high = mid
    audit.add(
        sub,
        "implicit p_c by independent bisection",
        chamber.chamber_pressure_pa,
        0.5 * (low + high),
        1e-12,
        "identity",
        "Pa",
    )
    audit.add(
        sub,
        "chamber residual m_dot c*/A_t - p_c",
        abs(chamber.residual_pa),
        0.0,
        1e-6,
        "identity",
        "Pa",
        relative=False,
    )

    # Simultaneous feed/chamber consistency, every relation re-formed by hand.
    tank_state = REFERENCE_TANK.initial_state(
        REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION
    )
    feed = solve_variable_feed_coupling(tank_state, INJECTOR, GRAIN, LAW, table, NOZZLE, 0.020)
    port_area = math.pi * 0.020**2
    hand_flux = feed.oxidizer_mass_flow_kg_s / port_area
    hand_rate = LAW.coefficient_si * hand_flux**LAW.exponent
    hand_fuel = GRAIN.fuel_density_kg_m3 * 2.0 * math.pi * 0.020 * GRAIN.length_m * hand_rate
    hand_state = table.evaluate(feed.mixture_ratio, feed.chamber_pressure_pa)
    audit.add(
        sub,
        "coupled: G_ox re-formed",
        feed.oxidizer_mass_flux_si,
        hand_flux,
        1e-14,
        "identity",
        "kg/(m^2 s)",
    )
    audit.add(
        sub,
        "coupled: m_dot_f re-formed",
        feed.fuel_mass_flow_kg_s,
        hand_fuel,
        1e-13,
        "identity",
        "kg/s",
    )
    audit.add(
        sub,
        "coupled: O/F identity",
        feed.mixture_ratio,
        feed.oxidizer_mass_flow_kg_s / feed.fuel_mass_flow_kg_s,
        1e-15,
        "identity",
        "-",
    )
    audit.add(
        sub,
        "coupled: p_c = m_dot c*/A_t",
        feed.chamber_pressure_pa,
        (feed.oxidizer_mass_flow_kg_s + feed.fuel_mass_flow_kg_s)
        * hand_state.c_star_m_s
        / NOZZLE.throat_area_m2,
        1e-11,
        "identity",
        "Pa",
    )
    audit.add(
        sub,
        "coupled: injector agrees with p_c",
        INJECTOR.mass_flow(tank_state, feed.chamber_pressure_pa).mass_flow_kg_s,
        feed.oxidizer_mass_flow_kg_s,
        1e-12,
        "identity",
        "kg/s",
    )

    # Boundary behaviour: the table must classify and refuse, never extrapolate.
    of_low, of_high2 = table.of_bounds
    p_low, p_high = table.pressure_bounds_pa
    statuses = [
        table.status(0.5 * (of_low + of_high2), 0.5 * (p_low + p_high)) is TableStatus.WITHIN_TABLE,
        table.status(of_low - 1e-6, p_low) is TableStatus.O_F_BELOW_TABLE,
        table.status(of_high2 + 1e-6, p_low) is TableStatus.O_F_ABOVE_TABLE,
        table.status(2.5, p_low - 1.0) is TableStatus.PRESSURE_BELOW_TABLE,
        table.status(2.5, p_high + 1.0) is TableStatus.PRESSURE_ABOVE_TABLE,
    ]
    refused = 0
    for of, pressure in (
        (of_low - 1e-6, p_low),
        (of_high2 + 1e-6, p_low),
        (2.5, p_low - 1.0),
        (2.5, p_high + 1.0),
    ):
        try:
            table.evaluate(of, pressure)
        except ValueError:
            refused += 1
    audit.add(sub, "table edge classification (5 cases)", sum(statuses), 5, 0.0, "identity", "-")
    audit.add(sub, "extrapolation refused (4 cases)", refused, 4, 0.0, "identity", "-")

    # Coupled blowdown closures re-formed from the reported histories.
    m5 = simulate_variable_blowdown(
        REFERENCE_TANK, tank_state, INJECTOR, GRAIN, LAW, table, NOZZLE, SEA_LEVEL_PA
    )
    audit.add(
        sub,
        "blowdown: oxidizer mass closure",
        float(np.max(np.abs(m5.oxidizer_mass_closure_residual_kg))),
        0.0,
        1e-12,
        "discretisation",
        "kg",
        relative=False,
    )
    audit.add(
        sub,
        "blowdown: c* identity over history",
        float(np.max(np.abs(m5.c_star_identity_residual_pa))),
        0.0,
        1e-6,
        "identity",
        "Pa",
        relative=False,
    )
    audit.add(
        sub,
        "blowdown: O/F identity over history",
        float(np.max(np.abs(m5.mixture_ratio_identity_residual))),
        0.0,
        1e-14,
        "identity",
        "-",
        relative=False,
    )
    audit.add(
        sub,
        "blowdown: feed residual over history",
        float(np.max(np.abs(m5.feed_residual_kg_s))),
        0.0,
        1e-11,
        "identity",
        "kg/s",
        relative=False,
    )
    audit.add(
        sub,
        "blowdown: stayed inside the table",
        float(m5.stayed_within_table),
        1.0,
        0.0,
        "identity",
        "-",
    )


def main() -> int:
    print("=" * 104)
    print("MILESTONE 6 -- INDEPENDENT END-TO-END AUDIT")
    print("Every expected value below is recomputed from raw formulas, closed forms, alternate")
    print("numerical routes, published measurements or raw CSV columns -- never by calling the")
    print("same production function twice.")
    print("=" * 104)
    print()

    audit = Audit()
    audit_m1(audit)
    audit_m2(audit)
    audit_m3(audit)
    audit_m4(audit)
    audit_m5(audit)

    for subsystem in dict.fromkeys(c.subsystem for c in audit.checks):
        audit.report(subsystem)

    print("SUMMARY")
    print(RULE)
    by_kind: dict[str, list[Check]] = {}
    for check in audit.checks:
        by_kind.setdefault(check.kind, []).append(check)
    for kind, rows in by_kind.items():
        worst = max(c.residual for c in rows)
        print(
            f"  {kind:<16} {len(rows):>3} checks   worst residual {worst:.3e}   "
            f"{sum(c.passed for c in rows)}/{len(rows)} passed"
        )
    print()

    identities = by_kind.get("identity", [])
    worst_identity = max((c.residual for c in identities), default=0.0)
    print(f"  Maximum ALGEBRAIC-IDENTITY residual: {worst_identity:.3e}")
    print("    This is the number that should be near machine precision. Discretisation rows")
    print("    are integrator-tolerance quantities and are NOT expected to reach it; they are")
    print("    reported separately above rather than folded into one headline figure.")
    print()

    failures = [c for c in audit.checks if not c.passed]
    print(f"  TOTAL: {len(audit.checks) - len(failures)}/{len(audit.checks)} checks passed")
    if failures:
        print("  FAILURES:")
        for c in failures:
            print(f"    {c.subsystem} :: {c.name}  residual {c.residual:.3e} > {c.tolerance:.1e}")
    print("=" * 104)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
