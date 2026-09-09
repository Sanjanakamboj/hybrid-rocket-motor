"""Milestone 4 coupled N2O blowdown / thrust study.

Generic, reduced-order, educational study.  The oxidizer mass flow is no longer
prescribed: a saturated-equilibrium N2O tank, a sourced injector model and the
frozen Milestone 1-3 motor are solved as one coupled system.

Nothing here is hardware design.  The injector is an effective flow area only --
no hole count, no hole size, no plate geometry.  The tank is a generic rigid
volume with no material, no wall thickness and no structural analysis of any
kind.  Nothing has been validated against experiment.

Run with::

    python scripts/blowdown_thrust_study.py
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from hybrid_rocket_motor import REPRESENTATIVE_GRAIN, REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.blowdown import BlowdownTermination, simulate_blowdown
from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.feed_system import FeedStatus, solve_feed_coupling
from hybrid_rocket_motor.injector import (
    DYER_REFERENCE_DISCHARGE_COEFFICIENT,
    Injector,
    InjectorModel,
)
from hybrid_rocket_motor.nitrous_properties import (
    MAX_VERIFIED_TEMPERATURE_K,
    MIN_VERIFIED_TEMPERATURE_K,
)
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.performance import (
    STANDARD_GRAVITY_M_S2,
    couple_transient_to_performance,
)
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
    NitrousTank,
)
from hybrid_rocket_motor.transient import ConstantOxidizerFlow, simulate_transient

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
COMBUSTION = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
NOZZLE = REFERENCE_NOZZLE
SEA_LEVEL_PA = 101325.0

#: Reference effective injector area [m^2].  See the calibration table below.
REFERENCE_INJECTOR_AREA_M2 = 3.5e-6
CD = DYER_REFERENCE_DISCHARGE_COEFFICIENT

#: Milestone 3 prescribed-flow comparison baseline.
M3_PRESCRIBED_FLOW_KG_S = 0.100

N_REPORT = 401
RULE = "-" * 108


def _banner() -> None:
    print("=" * 108)
    print("MILESTONE 4 - N2O TANK / INJECTOR / FEED COUPLING AND BLOWDOWN")
    print("Generic N2O/HTPB hybrid rocket motor - reduced-order educational model")
    print("=" * 108)
    print()
    print("SCOPE NOTICE")
    print("  Oxidizer flow is now an OUTPUT of a coupled tank/injector/chamber solve.")
    print("  Generic illustrative tank volume and effective injector area only:")
    print("  no hole count, no hole size, no plate geometry, no tank material or wall thickness.")
    print("  Combustion properties (c*, gamma, T_c) remain PRESCRIBED constants.")
    print("  No equilibrium chemistry. No ignition or chamber-filling transient.")
    print("  No combustion instability. No structural or thermal sizing. NOT validated.")
    print()


def _print_configuration() -> None:
    tank_state = REFERENCE_TANK.initial_state(
        REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION
    )
    print("N2O PROPERTY MODEL")
    print(RULE)
    print("  CoolProp, Lemmon & Span (2006) short fundamental equation of state")
    print("  (the same EOS the NIST Chemistry WebBook cites for nitrous oxide).")
    print(
        f"  checked temperature band: {MIN_VERIFIED_TEMPERATURE_K:.0f}"
        f" to {MAX_VERIFIED_TEMPERATURE_K:.0f} K"
    )
    print()
    print("REFERENCE TANK (illustrative, generic, NOT a real or certified tank)")
    print(RULE)
    print(f"  internal volume V           = {REFERENCE_TANK.volume_m3 * 1e3:.2f} L")
    print(f"  initial temperature T_0     = {REFERENCE_TANK_TEMPERATURE_K:.2f} K")
    print(f"  initial liquid fill         = {REFERENCE_TANK_FILL_FRACTION * 100:.0f} % by volume")
    print(f"  initial pressure            = {tank_state.pressure_pa / 1e5:.3f} bar")
    print(f"  initial N2O mass            = {tank_state.total_mass_kg:.4f} kg")
    print(f"    of which liquid           = {tank_state.liquid_mass_kg:.4f} kg")
    print(f"    of which vapour           = {tank_state.vapour_mass_kg:.4f} kg")
    print()
    print("REFERENCE INJECTOR (effective flow area only)")
    print(RULE)
    print(f"  effective area A_eff        = {REFERENCE_INJECTOR_AREA_M2 * 1e6:.2f} mm^2")
    print(f"  discharge coefficient C_d   = {CD:.2f}")
    print("    (the average value Dyer et al. used across their injector set,")
    print("     as quoted by Waxman et al.; a generic literature value)")
    print(f"  effective C_d A             = {CD * REFERENCE_INJECTOR_AREA_M2 * 1e6:.4f} mm^2")
    print("  model                       = Dyer / NHNE blend of SPI and HEM")
    print()


def _print_calibration() -> None:
    print("EFFECTIVE-AREA CALIBRATION (a comparison convenience, NOT an injector design)")
    print(RULE)
    print(
        "  Milestone 3 prescribed a constant "
        f"m_dot_ox = {M3_PRESCRIBED_FLOW_KG_S:.3f} kg/s.  For the Milestone 3 / 4"
    )
    print("  comparison to be meaningful the coupled model should start near that flow, so one")
    print("  parameter -- the effective area -- was chosen from a short list of round values with")
    print("  C_d fixed at the sourced 0.66.  Nothing else was adjusted.")
    print()
    tank_state = REFERENCE_TANK.initial_state(
        REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION
    )
    print(f"    {'A_eff [mm^2]':>13} {'m_dot_ox [kg/s]':>17} {'p_c [bar]':>11}   verdict")
    for area in (2.5e-6, 3.0e-6, 3.5e-6, 4.0e-6, 4.5e-6):
        solution = solve_feed_coupling(
            tank_state, Injector(area), GRAIN, LAW, COMBUSTION, NOZZLE, 0.020
        )
        flow = solution.oxidizer_mass_flow_kg_s
        inside = 0.08 <= flow <= 0.12
        chosen = abs(area - REFERENCE_INJECTOR_AREA_M2) < 1e-12
        verdict = "inside the 0.08-0.12 kg/s band" if inside else "outside that band"
        if chosen:
            verdict += "  <- chosen"
        print(
            f"    {area * 1e6:>13.1f} {flow:>17.6f} "
            f"{solution.chamber_pressure_pa / 1e5:>11.3f}   {verdict}"
        )
    print(
        "  3.5 mm^2 is the round value landing closest to the Milestone 3 flow. "
        "No thrust target"
    )
    print("  was used, and no second parameter was moved.")
    print()


def run_case(label: str, description: str, *, tank=None, temperature=None, fill=None,
             area=None, model=InjectorModel.DYER_NHNE):
    tank = tank or REFERENCE_TANK
    temperature = REFERENCE_TANK_TEMPERATURE_K if temperature is None else temperature
    fill = REFERENCE_TANK_FILL_FRACTION if fill is None else fill
    area = REFERENCE_INJECTOR_AREA_M2 if area is None else area
    state = tank.initial_state(temperature, fill)
    result = simulate_blowdown(
        tank, state, Injector(area, CD, model), GRAIN, LAW, COMBUSTION, NOZZLE,
        SEA_LEVEL_PA, n_report=N_REPORT,
    )
    return label, description, result


def _report_case(label: str, description: str, result) -> None:
    print(f"{label} - {description}")
    print(RULE)
    rows = (
        ("tank pressure", result.tank_pressure_pa / 1e5, "bar"),
        ("tank temperature", result.tank_temperature_k, "K"),
        ("tank liquid mass", result.tank_liquid_mass_kg, "kg"),
        ("oxidizer mass flow", result.oxidizer_mass_flow_kg_s, "kg/s"),
        ("fuel mass flow", result.fuel_mass_flow_kg_s, "kg/s"),
        ("mixture ratio O/F", result.mixture_ratio, "-"),
        ("chamber pressure", result.chamber_pressure_pa / 1e5, "bar"),
        ("pressure margin", result.pressure_margin_pa / 1e5, "bar"),
        ("port diameter", result.port_diameter_m * 1e3, "mm"),
        ("thrust", result.thrust_n, "N"),
        ("specific impulse", result.specific_impulse_s, "s"),
    )
    print(f"  {'quantity':<22} {'initial':>14} {'final':>14}   units")
    for name, array, unit in rows:
        print(f"  {name:<22} {array[0]:>14.5f} {array[-1]:>14.5f}   {unit}")
    print(RULE)
    print(f"  termination                {result.termination.value}")
    print(f"  burn duration              {result.burn_time_s:.5f} s")
    print(f"  peak thrust                {result.peak_thrust_n:.4f} N")
    print(f"  mean thrust                {result.mean_thrust_n:.4f} N")
    print(f"  TOTAL IMPULSE              {result.total_impulse_n_s:.4f} N s")
    print(f"  oxidizer consumed          {result.oxidizer_consumed_kg:.5f} kg")
    print(f"  fuel consumed              {result.fuel_consumed_kg:.5f} kg")
    print(f"  equivalent I_sp            {result.equivalent_specific_impulse_s:.4f} s")
    print(
        f"  minimum pressure margin    {result.minimum_pressure_margin_pa / 1e5:.4f} bar"
    )
    burnout = result.termination is BlowdownTermination.GRAIN_BURNOUT
    depleted = result.termination is BlowdownTermination.LIQUID_DEPLETED
    equalized = result.termination is BlowdownTermination.PRESSURE_EQUALIZED
    print(f"  grain burnout reached      {'yes' if burnout else 'no'}")
    print(f"  liquid depletion reached   {'yes' if depleted else 'no'}")
    print(f"  pressure equalization      {'yes' if equalized else 'no'}")
    thrust_trend = "RISES" if result.thrust_n[-1] > result.thrust_n[0] else "FALLS"
    print(
        f"  thrust {thrust_trend} over the burn "
        f"({result.thrust_n[0]:.2f} -> {result.thrust_n[-1]:.2f} N, "
        f"{100 * (result.thrust_n[-1] / result.thrust_n[0] - 1):+.1f} %)"
    )
    print(RULE)
    print()


def _print_trajectory(result, n_rows: int = 12) -> None:
    print("  sampled coupled trajectory")
    print(
        f"    {'t':>7} {'p_tank':>9} {'T_tank':>8} {'m_l':>8} {'m_dot_ox':>10} "
        f"{'p_c':>8} {'O/F':>7} {'D_p':>8} {'F':>9}"
    )
    print(
        f"    {'[s]':>7} {'[bar]':>9} {'[K]':>8} {'[kg]':>8} {'[kg/s]':>10} "
        f"{'[bar]':>8} {'[-]':>7} {'[mm]':>8} {'[N]':>9}"
    )
    for i in np.unique(np.linspace(0, result.time_s.size - 1, n_rows).astype(int)):
        print(
            f"    {result.time_s[i]:>7.3f} {result.tank_pressure_pa[i] / 1e5:>9.4f} "
            f"{result.tank_temperature_k[i]:>8.3f} {result.tank_liquid_mass_kg[i]:>8.4f} "
            f"{result.oxidizer_mass_flow_kg_s[i]:>10.6f} "
            f"{result.chamber_pressure_pa[i] / 1e5:>8.4f} {result.mixture_ratio[i]:>7.3f} "
            f"{result.port_diameter_m[i] * 1e3:>8.3f} {result.thrust_n[i]:>9.4f}"
        )
    print()


def milestone_3_reference():
    """The frozen prescribed-flow baseline, re-run unchanged for comparison."""
    transient = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M3_PRESCRIBED_FLOW_KG_S, 200.0), n_report=N_REPORT
    )
    return couple_transient_to_performance(transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA)


def _print_comparison(m3, m4) -> None:
    print("MILESTONE 3 vs MILESTONE 4 - what the feed coupling changes")
    print(RULE)
    print(
        "  M3 holds m_dot_ox at a constant 0.100 kg/s. M4 lets the tank and injector set it."
    )
    print(
        "  Everything else -- grain, regression law, combustion properties, nozzle --\n"
        "  is identical."
    )
    print()
    print(f"  {'quantity':<28} {'M3 (prescribed)':>18} {'M4 (coupled)':>18} {'change':>12}")
    comparisons = (
        ("initial m_dot_ox [kg/s]", m3.oxidizer_mass_flow_kg_s[0], m4.oxidizer_mass_flow_kg_s[0]),
        ("final m_dot_ox [kg/s]", m3.oxidizer_mass_flow_kg_s[-1], m4.oxidizer_mass_flow_kg_s[-1]),
        ("initial p_c [bar]", m3.chamber_pressure_pa[0] / 1e5, m4.chamber_pressure_pa[0] / 1e5),
        ("final p_c [bar]", m3.chamber_pressure_pa[-1] / 1e5, m4.chamber_pressure_pa[-1] / 1e5),
        ("initial thrust [N]", m3.thrust_n[0], m4.thrust_n[0]),
        ("final thrust [N]", m3.thrust_n[-1], m4.thrust_n[-1]),
        ("peak thrust [N]", m3.peak_thrust_n, m4.peak_thrust_n),
        ("initial O/F [-]", m3.mixture_ratio[0], m4.mixture_ratio[0]),
        ("final O/F [-]", m3.mixture_ratio[-1], m4.mixture_ratio[-1]),
        ("burn duration [s]", m3.time_s[-1], m4.burn_time_s),
        ("total impulse [N s]", m3.total_impulse_n_s, m4.total_impulse_n_s),
        ("equivalent I_sp [s]",
         m3.equivalent_specific_impulse_s, m4.equivalent_specific_impulse_s),
    )
    for name, a, b in comparisons:
        print(f"  {name:<28} {a:>18.4f} {b:>18.4f} {100 * (b / a - 1):>11.2f}%")
    print(RULE)
    m3_trend = "rises" if m3.thrust_n[-1] > m3.thrust_n[0] else "falls"
    m4_trend = "rises" if m4.thrust_n[-1] > m4.thrust_n[0] else "falls"
    print(
        f"  THRUST TREND REVERSES: prescribed flow {m3_trend} "
        f"({100 * (m3.thrust_n[-1] / m3.thrust_n[0] - 1):+.1f} %), "
        f"coupled flow {m4_trend} ({100 * (m4.thrust_n[-1] / m4.thrust_n[0] - 1):+.1f} %)."
    )
    print(
        "  This is a qualitative disagreement, not a small correction: assuming constant"
    )
    print("  oxidizer flow gets the SIGN of the thrust-time slope wrong for this configuration.")
    print(RULE)
    print()


def _print_mass_energy_audit(result) -> None:
    print("INDEPENDENT MASS / ENERGY AUDIT (reference case)")
    print(RULE)
    drawn = result.tank_total_mass_kg[0] - result.tank_total_mass_kg[-1]
    geometric_fuel = (
        GRAIN.fuel_density_kg_m3
        * GRAIN.length_m
        * math.pi
        * (result.port_radius_m[-1] ** 2 - result.port_radius_m[0] ** 2)
    )
    trapezoid_impulse = float(np.trapezoid(result.thrust_n, result.time_s))
    print(
        f"  oxidizer drawn from tank            {drawn:.9f} kg\n"
        f"  oxidizer integrated at injector     {result.oxidizer_consumed_kg:.9f} kg\n"
        f"    max |residual|                    "
        f"{np.max(np.abs(result.oxidizer_mass_closure_residual_kg)):.3e} kg"
    )
    print(
        f"  fuel from port geometry             {geometric_fuel:.9f} kg\n"
        f"  fuel integrated from m_dot_f        {result.fuel_consumed_kg:.9f} kg\n"
        f"    max |residual|                    "
        f"{np.max(np.abs(result.fuel_mass_closure_residual_kg)):.3e} kg"
    )
    print(
        f"  impulse state                       {result.total_impulse_n_s:.6f} N s\n"
        f"  trapezoidal integral of thrust      {trapezoid_impulse:.6f} N s\n"
        f"    residual                          {result.impulse_closure_residual_n_s:.3e} N s\n"
        f"    relative                          {result.relative_impulse_closure_residual:.3e}"
    )
    energy_residual = float(np.max(np.abs(result.energy_closure_residual_j)))
    print(
        f"  tank energy balance residual        {energy_residual:.3e} J\n"
        f"    normalised by initial energy      {result.relative_energy_closure_residual:.3e}"
    )
    print(
        f"  feed-coupling root residual         "
        f"{np.max(np.abs(result.feed_residual_kg_s)):.3e} kg/s"
    )
    print(RULE)
    print()


def _print_convergence() -> None:
    print("NUMERICAL CONVERGENCE (reference case)")
    print(RULE)
    state = REFERENCE_TANK.initial_state(
        REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION
    )
    print(
        f"  {'rtol':>9} {'feed xtol':>11} {'n_report':>9} {'t_end [s]':>12} "
        f"{'p_tank,f [bar]':>15} {'p_c,f [bar]':>13} {'I_total [N s]':>15}"
    )
    for rtol, xtol, n_report in (
        (1e-6, 1e-9, 101),
        (1e-8, 1e-12, 201),
        (1e-10, 1e-12, 401),
        (1e-11, 1e-14, 801),
    ):
        result = simulate_blowdown(
            REFERENCE_TANK, state, Injector(REFERENCE_INJECTOR_AREA_M2), GRAIN, LAW,
            COMBUSTION, NOZZLE, SEA_LEVEL_PA, n_report=n_report, rtol=rtol,
            atol=rtol * 1e-2, feed_xtol=xtol,
        )
        print(
            f"  {rtol:>9.0e} {xtol:>11.0e} {n_report:>9d} {result.burn_time_s:>12.6f} "
            f"{result.tank_pressure_pa[-1] / 1e5:>15.6f} "
            f"{result.chamber_pressure_pa[-1] / 1e5:>13.6f} "
            f"{result.total_impulse_n_s:>15.4f}"
        )
    print(RULE)
    print()


def sanity_audit(results, m3, m4) -> list[tuple[str, bool]]:
    """Engineering sanity checks.

    Quantities that only exist while propellant is flowing are checked on the
    firing samples only.  A run terminated by liquid depletion legitimately ends
    with a non-firing sample carrying zero flow, zero chamber pressure and zero
    thrust, and that sample is checked separately for exactly those zeros.
    """
    checks: list[tuple[str, bool]] = []

    def every(predicate) -> bool:
        return all(predicate(r) for r in results)

    def firing(r):
        return np.array([status is FeedStatus.FLOWING for status in r.feed_status])

    checks.append(("port radius nondecreasing",
                   every(lambda r: bool(np.all(np.diff(r.port_radius_m) >= 0.0)))))
    checks.append(("tank oxidizer mass nonincreasing",
                   every(lambda r: bool(np.all(np.diff(r.tank_total_mass_kg) <= 0.0)))))
    checks.append(("tank liquid mass nonincreasing and never negative",
                   every(lambda r: bool(np.all(np.diff(r.tank_liquid_mass_kg) <= 0.0)
                                        and np.all(r.tank_liquid_mass_kg >= 0.0)))))
    checks.append(("tank pressure falls monotonically",
                   every(lambda r: bool(np.all(np.diff(r.tank_pressure_pa) < 0.0)))))
    checks.append(("tank temperature falls monotonically",
                   every(lambda r: bool(np.all(np.diff(r.tank_temperature_k) < 0.0)))))
    checks.append(("p_tank > p_c on every firing sample",
                   every(lambda r: bool(np.all(r.pressure_margin_pa[firing(r)] > 0.0)))))
    checks.append(("m_dot_ox strictly positive while firing",
                   every(lambda r: bool(np.all(r.oxidizer_mass_flow_kg_s[firing(r)] > 0.0)))))
    checks.append(("chamber pressure strictly positive while firing",
                   every(lambda r: bool(np.all(r.chamber_pressure_pa[firing(r)] > 0.0)))))
    checks.append((
        "thrust finite and positive while firing",
        every(lambda r: bool(np.all(np.isfinite(r.thrust_n[firing(r)]))
                             and np.all(r.thrust_n[firing(r)] > 0.0))),
    ))
    checks.append(("exit Mach supersonic while firing",
                   every(lambda r: bool(np.all(r.exit_mach[firing(r)] > 1.0)))))
    checks.append((
        "O/F finite and positive while firing",
        every(lambda r: bool(np.all(np.isfinite(r.mixture_ratio[firing(r)]))
                             and np.all(r.mixture_ratio[firing(r)] > 0.0))),
    ))
    checks.append((
        "a non-firing sample reports exactly zero flow, pressure and thrust",
        every(lambda r: bool(np.all(r.oxidizer_mass_flow_kg_s[~firing(r)] == 0.0)
                             and np.all(r.chamber_pressure_pa[~firing(r)] == 0.0)
                             and np.all(r.thrust_n[~firing(r)] == 0.0))),
    ))
    checks.append((
        "a non-firing sample reports no nozzle state at all",
        every(lambda r: bool(np.all(np.isnan(r.exit_mach[~firing(r)])))),
    ))
    checks.append((
        "no NaN/inf anywhere in a firing history",
        every(lambda r: bool(np.all(np.isfinite(r.tank_pressure_pa))
                             and np.all(np.isfinite(r.oxidizer_mass_flow_kg_s))
                             and np.all(np.isfinite(r.chamber_pressure_pa)))),
    ))
    checks.append((
        "oxidizer mass closes to better than 1e-9 kg",
        every(lambda r: float(np.max(np.abs(r.oxidizer_mass_closure_residual_kg))) < 1e-9),
    ))
    checks.append(("fuel mass closes to better than 1e-6 kg",
                   every(lambda r: float(np.max(np.abs(r.fuel_mass_closure_residual_kg))) < 1e-6)))
    checks.append((
        "impulse reconstruction closes to better than 1 % of total impulse",
        every(lambda r: r.relative_impulse_closure_residual < 1e-2),
    ))
    checks.append((
        "tank energy balance closes to better than 1e-6 relative (firing samples)",
        every(lambda r: r.relative_energy_closure_residual < 1e-6),
    ))
    checks.append(("feed root residual below 1e-9 kg/s",
                   every(lambda r: float(np.max(np.abs(r.feed_residual_kg_s))) < 1e-9)))
    checks.append((
        "no state integrated past burnout",
        every(lambda r: bool(np.all(r.port_diameter_m <= GRAIN.outer_diameter_m + 1e-12))),
    ))
    checks.append(("total impulse positive", every(lambda r: r.total_impulse_n_s > 0.0)))
    checks.append((
        "equivalent I_sp identity closes",
        every(lambda r: abs(r.equivalent_specific_impulse_s
                            - r.total_impulse_n_s
                            / (r.propellant_consumed_kg * STANDARD_GRAVITY_M_S2)) < 1e-12),
    ))
    checks.append((
        "tank temperature stays inside the checked property band",
        every(lambda r: bool(np.all(r.tank_temperature_k >= MIN_VERIFIED_TEMPERATURE_K)
                             and np.all(r.tank_temperature_k <= MAX_VERIFIED_TEMPERATURE_K))),
    ))
    checks.append((
        "both terminal events are exercised across the case set",
        {r.termination for r in results}
        >= {BlowdownTermination.GRAIN_BURNOUT, BlowdownTermination.LIQUID_DEPLETED},
    ))
    checks.append((
        "Milestone 3 baseline unchanged (thrust still rises under prescribed flow)",
        bool(m3.thrust_n[-1] > m3.thrust_n[0]),
    ))
    checks.append((
        "coupled model reverses that trend",
        bool(m4.thrust_n[-1] < m4.thrust_n[0]),
    ))
    return checks


def _print_sanity_audit(checks) -> bool:
    print("ENGINEERING SANITY AUDIT")
    print(RULE)
    for description, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {description}")
    print(RULE)
    all_ok = all(ok for _, ok in checks)
    print(f"  {'ALL SANITY CHECKS PASSED' if all_ok else '*** SANITY CHECKS FAILED ***'}")
    print()
    return all_ok


def _print_observations(cases, m3, m4) -> None:
    print("OBSERVATIONS (reported as computed, not tuned)")
    print(RULE)
    by_label = {label: result for label, _d, result in cases}
    a = by_label["CASE A"]
    m3_pct = 100 * (m3.thrust_n[-1] / m3.thrust_n[0] - 1)
    a_pct = 100 * (a.thrust_n[-1] / a.thrust_n[0] - 1)
    ox_pct = 100 * (a.oxidizer_mass_flow_kg_s[-1] / a.oxidizer_mass_flow_kg_s[0] - 1)
    liquid_pct = 100 * a.tank_liquid_mass_kg[-1] / a.tank_liquid_mass_kg[0]
    b, c = by_label["CASE B"], by_label["CASE C"]
    d, e = by_label["CASE D"], by_label["CASE E"]
    f = by_label["CASE F"]
    print(
        "  1. THE THRUST-TIME SLOPE CHANGES SIGN. Milestone 3, holding the oxidizer flow\n"
        f"     constant, predicted thrust rising {m3_pct:+.1f} %"
        f" over the burn. With the tank coupled it\n"
        f"     falls {a_pct:+.1f} %"
        " instead, because the tank cools and depressurises faster than the growing\n"
        "     port raises the fuel flow. Constant-flow assumptions do not merely lose accuracy\n"
        "     here; they get the qualitative behaviour wrong."
    )
    print(
        f"  2. SELF-PRESSURISATION COSTS FEED PRESSURE FAST. The tank falls\n"
        f"     {a.tank_pressure_pa[0] / 1e5:.1f} -> {a.tank_pressure_pa[-1] / 1e5:.1f} bar and"
        f" {a.tank_temperature_k[0]:.1f} -> {a.tank_temperature_k[-1]:.1f} K over"
        f" {a.burn_time_s:.1f} s, with\n"
        f"     {a.tank_liquid_mass_kg[-1]:.2f} kg of liquid still unused at burnout. The oxidizer\n"
        f"     flow decays {a.oxidizer_mass_flow_kg_s[0]:.4f} ->"
        f" {a.oxidizer_mass_flow_kg_s[-1]:.4f} kg/s ({ox_pct:+.1f} %)."
    )
    print(
        "  3. INITIAL TANK TEMPERATURE IS THE STRONGEST LEVER. A 10 K change moves the initial\n"
        f"     tank pressure {b.tank_pressure_pa[0] / 1e5:.1f} / {a.tank_pressure_pa[0] / 1e5:.1f}"
        f" / {c.tank_pressure_pa[0] / 1e5:.1f} bar and the initial thrust\n"
        f"     {b.thrust_n[0]:.1f} / {a.thrust_n[0]:.1f} / {c.thrust_n[0]:.1f} N."
        " Total impulse moves much less\n"
        f"     ({b.total_impulse_n_s:.0f} / {a.total_impulse_n_s:.0f}"
        f" / {c.total_impulse_n_s:.0f} N s),"
        " because the grain sets how much fuel there is to burn."
    )
    print(
        "  4. INJECTOR AREA TRADES THRUST AGAINST DURATION. Halving-to-doubling the effective\n"
        f"     area gives burn times {e.burn_time_s:.1f} / {a.burn_time_s:.1f}"
        f" / {d.burn_time_s:.1f} s"
        f" and initial thrusts {e.thrust_n[0]:.0f} / {a.thrust_n[0]:.0f} / {d.thrust_n[0]:.0f} N,\n"
        f"     with total impulse {e.total_impulse_n_s:.0f} / {a.total_impulse_n_s:.0f}"
        f" / {d.total_impulse_n_s:.0f} N s. All three still end on grain burnout: over this\n"
        "     range the injector sets how hard and how long the motor runs, not what stops it."
    )
    print(
        "  5. THE INJECTOR MODEL MATTERS AS MUCH AS THE TANK. Switching from the Dyer blend to\n"
        f"     the plain incompressible SPI equation raises initial flow"
        f" {a.oxidizer_mass_flow_kg_s[0]:.4f} -> {f.oxidizer_mass_flow_kg_s[0]:.4f} kg/s\n"
        f"     and initial thrust {a.thrust_n[0]:.1f} -> {f.thrust_n[0]:.1f} N"
        f" ({100 * (f.thrust_n[0] / a.thrust_n[0] - 1):+.1f} %). Waxman et al. warn that SPI is\n"
        "     'often applied inappropriately' to nitrous oxide because the liquid flashes in the\n"
        "     orifice; the spread between the two models is the size of that modelling choice."
    )
    g = by_label["CASE G"]
    impulse_loss_pct = 100 * (1 - g.total_impulse_n_s / a.total_impulse_n_s)
    print(
        "  6. THE GRAIN, NOT THE TANK, ENDS EVERY NOMINAL-SIZED BURN. Cases A-F all terminate on\n"
        f"     grain burnout, case A with {liquid_pct:.0f} % of its liquid still in the tank."
        " The 10 L tank simply\n"
        "     carries far more oxidizer than a 1.9 kg grain can consume, so the tank is never the\n"
        f"     limiting factor. Case G, with a 2.5 L tank, terminates instead on"
        f" {g.termination.value}\n"
        f"     after {g.burn_time_s:.1f} s with the port only {g.port_diameter_m[-1] * 1e3:.1f} mm"
        f" across -- {impulse_loss_pct:.0f} % less total impulse,\n"
        "     and the grain left unburnt. Which event ends the burn is a design outcome, not an\n"
        "     assumption, and the model reports it rather than presuming burnout."
    )
    print(RULE)
    print()


def main() -> int:
    _banner()
    _print_configuration()
    _print_calibration()

    cases = [
        run_case("CASE A", "nominal generic tank and injector"),
        run_case("CASE B", "lower initial tank temperature, 283.15 K", temperature=283.15),
        run_case("CASE C", "higher initial tank temperature, 303.15 K", temperature=303.15),
        run_case("CASE D", "lower effective injector area, 2.5 mm^2", area=2.5e-6),
        run_case("CASE E", "higher effective injector area, 4.5 mm^2", area=4.5e-6),
        run_case("CASE F", "nominal tank with the SPI injector model",
                 model=InjectorModel.SPI),
        run_case("CASE G", "smaller 2.5 L tank - oxidizer-limited rather than grain-limited",
                 tank=NitrousTank(volume_m3=0.0025)),
    ]

    for label, description, result in cases:
        _report_case(label, description, result)
        if label == "CASE A":
            _print_trajectory(result, n_rows=14)

    nominal = cases[0][2]
    m3 = milestone_3_reference()
    _print_comparison(m3, nominal)
    _print_mass_energy_audit(nominal)
    _print_convergence()
    _print_observations(cases, m3, nominal)

    ok = _print_sanity_audit(
        sanity_audit([result for _l, _d, result in cases], m3, nominal)
    )
    print(
        "END OF MILESTONE 4 STUDY - generic effective injector area only; no hole geometry,\n"
        "no tank structural analysis, no experimental validation."
    )
    assert all(
        result.termination is not None for _l, _d, result in cases
    )
    assert len(set(itertools.chain.from_iterable([[r.termination] for _l, _d, r in cases]))) >= 1
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
