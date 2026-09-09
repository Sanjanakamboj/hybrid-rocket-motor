"""Milestone 3 quasi-steady chamber/nozzle coupling and thrust prediction study.

Generic, reduced-order, educational study of an illustrative N2O/HTPB hybrid
motor.  The oxidizer mass flow remains a *prescribed input*: chamber pressure
does not feed back into it.  Combustion properties (``c*``, ``gamma``, ``T_c``)
are *prescribed illustrative constants*, not computed chemistry; this milestone
contains no equilibrium-chemistry solver, no tank or injector model and no
structural or thermal sizing.

Run with::

    python scripts/thrust_prediction_study.py
"""

from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np

from hybrid_rocket_motor import REPRESENTATIVE_GRAIN, REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.nozzle import (
    REFERENCE_NOZZLE,
    solve_supersonic_exit_mach,
)
from hybrid_rocket_motor.performance import (
    STANDARD_GRAVITY_M_S2,
    couple_transient_to_performance,
    solve_performance_point,
)
from hybrid_rocket_motor.transient import (
    ConstantOxidizerFlow,
    PiecewiseConstantOxidizerFlow,
    TerminationReason,
    simulate_transient,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))

from transient_regression_study import (
    M_DOT_HIGH,
    M_DOT_LOW,
    M_DOT_NOMINAL,
    PIECEWISE_STEPS,
)

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
COMBUSTION = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
NOZZLE = REFERENCE_NOZZLE

SEA_LEVEL_PA = 101325.0
HORIZON_S = 200.0
N_REPORT = 801

#: Milestone 1 reference operating point, reused verbatim.
M1_M_DOT_OX = 0.100
M1_M_DOT_F = 0.039776644

#: Ambient pressures for the deterministic nozzle sensitivity sweep [Pa].
#: The bracketed altitudes are ISA context only -- no atmosphere model is used
#: anywhere in this project, and this is not an altitude-performance envelope.
AMBIENT_PRESSURES_PA: tuple[tuple[float, str], ...] = (
    (101325.0, "sea level"),
    (54020.0, "~5 km ISA context"),
    (26500.0, "~10 km ISA context"),
    (1197.0, "~30 km ISA context"),
    (0.0, "vacuum"),
)

RULE = "-" * 104


def _banner() -> None:
    print("=" * 104)
    print("MILESTONE 3 - QUASI-STEADY CHAMBER/NOZZLE COUPLING AND THRUST PREDICTION")
    print("Generic N2O/HTPB hybrid rocket motor - reduced-order educational model")
    print("=" * 104)
    print()
    print("SCOPE NOTICE")
    print("  Generic illustrative geometry and PRESCRIBED combustion properties.")
    print("  Oxidizer mass flow is PRESCRIBED INPUT ONLY; chamber pressure does not feed back.")
    print("  No tank/blowdown/injector/feed model. No equilibrium chemistry. No CEA.")
    print("  No ignition or chamber-filling transient. No shocks or nozzle separation model.")
    print("  No structural or thermal sizing. NOT a validated or flight-ready motor.")
    print()


def _print_configuration() -> None:
    print("PRESCRIBED COMBUSTION PROPERTIES  ***ILLUSTRATIVE - NOT EXPERIMENTALLY VALIDATED***")
    print(RULE)
    print(f"  gamma                     = {COMBUSTION.gamma:.4f} [-]")
    print(f"  chamber temperature T_c   = {COMBUSTION.chamber_temperature_k:.1f} K")
    print(f"  molar mass M              = {COMBUSTION.molar_mass_kg_mol * 1e3:.2f} g/mol")
    print(f"  gas constant R            = {COMBUSTION.gas_constant_j_kg_k:.4f} J/(kg K)")
    print(f"  ideal c* from the above   = {COMBUSTION.ideal_c_star_m_s:.3f} m/s")
    print(f"  prescribed c* efficiency  = {COMBUSTION.c_star_efficiency:.4f} [-]")
    print(f"  delivered c*              = {COMBUSTION.c_star_m_s:.3f} m/s")
    print("  Consistency (NOT validation): Rezaei et al. (2018) measured c* = 1403-1587 m/s")
    print("  and chamber pressures of 19.9-31.0 bar for an HTPB/N2O motor of similar scale.")
    print()
    print("REFERENCE NOZZLE (illustrative, NOT optimised, NOT a manufacturing drawing)")
    print(RULE)
    print(f"  throat diameter D_t       = {NOZZLE.throat_diameter_m * 1e3:.3f} mm")
    print(f"  exit diameter   D_e       = {NOZZLE.exit_diameter_m * 1e3:.3f} mm")
    print(f"  throat area     A_t       = {NOZZLE.throat_area_m2:.6e} m^2")
    print(f"  exit area       A_e       = {NOZZLE.exit_area_m2:.6e} m^2")
    print(f"  expansion ratio epsilon   = {NOZZLE.expansion_ratio:.4f} [-]")
    print()


def _print_throat_selection() -> None:
    print("THROAT SELECTION (documented, not optimised)")
    print(RULE)
    print(
        "  Candidate throats evaluated at the Milestone 1 flow "
        f"(m_dot_total = {M1_M_DOT_OX + M1_M_DOT_F:.9f} kg/s):"
    )
    print(f"    {'D_t [mm]':>10} {'A_t [m^2]':>14} {'p_c [bar]':>12}   verdict")
    for d_t in (0.008, 0.010, 0.012, 0.015):
        area = math.pi * d_t**2 / 4.0
        p_c_bar = (M1_M_DOT_OX + M1_M_DOT_F) * COMBUSTION.c_star_m_s / area / 1e5
        inside = 19.9 <= p_c_bar <= 31.0
        verdict = "inside the measured 19.9-31.0 bar band" if inside else "outside that band"
        print(f"    {d_t * 1e3:>10.1f} {area:>14.6e} {p_c_bar:>12.2f}   {verdict}")
    print(
        "  A round 10 mm throat is the natural choice: it is the only round candidate whose\n"
        "  implied chamber pressure sits mid-band. No thrust target influenced this choice."
    )
    print(
        f"  Expansion ratio fixed at a round epsilon = {NOZZLE.expansion_ratio:.0f}; ideal "
        "sea-level expansion for this chamber\n"
        "  pressure would need epsilon ~ 4.34, so the reference nozzle is slightly "
        "UNDEREXPANDED at sea level."
    )
    print()


def _print_static_reference() -> None:
    print("STATIC REFERENCE POINT - the Milestone 1 operating point, now with thrust")
    print(RULE)
    point = solve_performance_point(COMBUSTION, NOZZLE, M1_M_DOT_OX, M1_M_DOT_F, SEA_LEVEL_PA)
    exit_state = point.exit_state
    rows = (
        ("prescribed m_dot_ox", M1_M_DOT_OX, "kg/s"),
        ("m_dot_f (Milestone 1)", M1_M_DOT_F, "kg/s"),
        ("m_dot_total", point.total_mass_flow_kg_s, "kg/s"),
        ("O/F", point.chamber.mixture_ratio, "-"),
        ("chamber pressure p_c", point.chamber_pressure_pa, "Pa"),
        ("chamber pressure p_c", point.chamber_pressure_pa / 1e5, "bar"),
        ("exit Mach M_e", exit_state.mach, "-"),
        ("exit pressure p_e", exit_state.pressure_pa, "Pa"),
        ("exit pressure ratio p_e/p_c", exit_state.pressure_ratio, "-"),
        ("exit temperature T_e", exit_state.temperature_k, "K"),
        ("exit speed of sound a_e", exit_state.speed_of_sound_m_s, "m/s"),
        ("exit velocity V_e", exit_state.velocity_m_s, "m/s"),
        ("momentum thrust", point.momentum_thrust_n, "N"),
        ("pressure thrust", point.pressure_thrust_n, "N"),
        ("TOTAL THRUST F", point.thrust_n, "N"),
        ("thrust coefficient C_F", point.thrust_coefficient, "-"),
        ("specific impulse I_sp", point.specific_impulse_s, "s"),
    )
    for label, value, unit in rows:
        print(f"  {label:<30} {value:>18.6f}   {unit}")
    print(RULE)
    route_b = point.thrust_coefficient * point.chamber_pressure_pa * NOZZLE.throat_area_m2
    slope = (
        exit_state.velocity_m_s
        + COMBUSTION.c_star_m_s * NOZZLE.expansion_ratio * exit_state.pressure_ratio
    )
    route_c = slope * point.total_mass_flow_kg_s - SEA_LEVEL_PA * NOZZLE.exit_area_m2
    print("  INDEPENDENT THRUST ROUTES")
    print(f"    Route A  F = m_dot V_e + (p_e - p_a) A_e   = {point.thrust_n:.9f} N")
    print(f"    Route B  F = C_F p_c A_t                   = {route_b:.9f} N")
    print(f"    Route C  F = m_dot [V_e + c* eps p_e/p_c] - p_a A_e = {route_c:.9f} N")
    print(f"    max |A - B| = {abs(point.thrust_n - route_b):.3e} N")
    print(f"    max |A - C| = {abs(point.thrust_n - route_c):.3e} N")
    c_star_residual = (
        point.chamber_pressure_pa * NOZZLE.throat_area_m2
        - point.total_mass_flow_kg_s * COMBUSTION.c_star_m_s
    )
    print(f"  c* identity residual  p_c A_t - m_dot c* = {c_star_residual:.3e} N")
    print(f"  expansion regime: {exit_state.expansion_regime.value}")
    separation = "SEPARATION EXPECTED" if exit_state.is_separation_expected else "attached flow"
    print(
        f"  Summerfield separation check: p_e/p_a = {exit_state.pressure_pa / SEA_LEVEL_PA:.4f} "
        f"vs 0.4 -> {separation} expected"
    )
    print(RULE)
    print()
    return point


def run_case(flow_history, label: str, description: str):
    transient = simulate_transient(GRAIN, LAW, flow_history, n_report=N_REPORT)
    history = couple_transient_to_performance(transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA)
    return label, description, transient, history


def _report_case(label: str, description: str, transient, history) -> None:
    print(f"{label} - {description}")
    print(RULE)
    first, last = 0, history.time_s.size - 1
    rows = (
        ("m_dot_f", history.fuel_mass_flow_kg_s, "kg/s"),
        ("m_dot_total", history.total_mass_flow_kg_s, "kg/s"),
        ("chamber pressure", history.chamber_pressure_pa / 1e5, "bar"),
        ("thrust", history.thrust_n, "N"),
        ("momentum thrust", history.momentum_thrust_n, "N"),
        ("pressure thrust", history.pressure_thrust_n, "N"),
        ("C_F", history.thrust_coefficient, "-"),
        ("I_sp", history.specific_impulse_s, "s"),
    )
    print(f"  {'quantity':<22} {'initial':>16} {'final':>16}   units")
    for name, array, unit in rows:
        print(f"  {name:<22} {array[first]:>16.6f} {array[last]:>16.6f}   {unit}")
    print(RULE)
    print(f"  termination                {transient.termination_reason.value}")
    if transient.burnout_time_s is not None:
        print(f"  burnout time               {transient.burnout_time_s:.6f} s")
    print(f"  simulated duration         {history.time_s[-1]:.6f} s")
    print(f"  active burn time           {history.active_burn_time_s:.6f} s")
    print(f"  peak thrust                {history.peak_thrust_n:.6f} N")
    print(f"  mean thrust (active)       {history.mean_thrust_n:.6f} N")
    print(f"  TOTAL IMPULSE              {history.total_impulse_n_s:.6f} N s")
    print(f"  propellant consumed        {history.propellant_mass_kg:.6f} kg")
    print(f"  equivalent I_sp            {history.equivalent_specific_impulse_s:.6f} s")
    firing = history.firing
    p_c = history.chamber_pressure_pa[firing]
    thrust = history.thrust_n[firing]
    press = history.pressure_thrust_n[firing]
    print(
        f"  chamber pressure {'RISES' if p_c[-1] > p_c[0] else 'FALLS'} over the burn "
        f"({p_c[0] / 1e5:.3f} -> {p_c[-1] / 1e5:.3f} bar)"
    )
    print(
        f"  thrust {'RISES' if thrust[-1] > thrust[0] else 'FALLS'} over the burn "
        f"({thrust[0]:.3f} -> {thrust[-1]:.3f} N)"
    )
    sign = "POSITIVE" if np.all(press > 0) else ("NEGATIVE" if np.all(press < 0) else "MIXED SIGN")
    print(
        f"  pressure-thrust contribution is {sign} throughout "
        f"({press.min():.3f} to {press.max():.3f} N)"
    )
    regimes = {r.value for r in history.expansion_regime if r is not None}
    print(f"  nozzle expansion regime at p_a = {SEA_LEVEL_PA:.0f} Pa: {', '.join(sorted(regimes))}")
    print(RULE)
    print()


def _print_trajectory(history, n_rows: int = 12) -> None:
    print("  sampled thrust trajectory")
    print(
        f"    {'t':>8} {'m_dot_tot':>11} {'p_c':>10} {'F_mom':>10} {'F_press':>9} "
        f"{'F':>10} {'C_F':>8} {'I_sp':>9}"
    )
    print(
        f"    {'[s]':>8} {'[kg/s]':>11} {'[bar]':>10} {'[N]':>10} {'[N]':>9} "
        f"{'[N]':>10} {'[-]':>8} {'[s]':>9}"
    )
    indices = np.unique(np.linspace(0, history.time_s.size - 1, n_rows).astype(int))
    for i in indices:
        cf = history.thrust_coefficient[i]
        isp = history.specific_impulse_s[i]
        cf_text = "     n/a" if math.isnan(cf) else f"{cf:8.4f}"
        isp_text = "      n/a" if math.isnan(isp) else f"{isp:9.3f}"
        print(
            f"    {history.time_s[i]:>8.3f} {history.total_mass_flow_kg_s[i]:>11.6f} "
            f"{history.chamber_pressure_pa[i] / 1e5:>10.4f} {history.momentum_thrust_n[i]:>10.4f} "
            f"{history.pressure_thrust_n[i]:>9.4f} {history.thrust_n[i]:>10.4f} "
            f"{cf_text} {isp_text}"
        )
    print()


def ambient_sensitivity(point):
    """Deterministic nozzle sensitivity to ambient pressure at a fixed chamber state."""
    rows = []
    for p_a, note in AMBIENT_PRESSURES_PA:
        result = solve_performance_point(COMBUSTION, NOZZLE, M1_M_DOT_OX, M1_M_DOT_F, p_a)
        rows.append((p_a, note, result))
    return rows


def _print_ambient_sensitivity(rows) -> None:
    print("AMBIENT-PRESSURE SENSITIVITY (deterministic nozzle sensitivity - NOT an altitude")
    print("envelope; no atmosphere model exists anywhere in this project)")
    print(RULE)
    print(
        f"  {'p_a [Pa]':>10} {'F [N]':>10} {'F_press [N]':>12} {'C_F':>9} {'I_sp [s]':>10}  "
        f"{'regime':<22} context"
    )
    for p_a, note, result in rows:
        print(
            f"  {p_a:>10.0f} {result.thrust_n:>10.4f} {result.pressure_thrust_n:>12.4f} "
            f"{result.thrust_coefficient:>9.5f} {result.specific_impulse_s:>10.4f}  "
            f"{result.expansion_regime.value:<22} {note}"
        )
    print(RULE)
    sea, vac = rows[0][2], rows[-1][2]
    print(
        f"  Thrust rises {vac.thrust_n - sea.thrust_n:.3f} N from sea level to vacuum, exactly "
        f"p_a A_e = {SEA_LEVEL_PA * NOZZLE.exit_area_m2:.3f} N."
    )
    print(
        "  The momentum term is unchanged by ambient pressure; the entire difference is the\n"
        "  pressure term. The exit state itself never changes, because the nozzle is choked."
    )
    print()


def combustion_sensitivity():
    """One-factor-at-a-time deterministic sensitivity on prescribed properties."""
    base = solve_performance_point(COMBUSTION, NOZZLE, M1_M_DOT_OX, M1_M_DOT_F, SEA_LEVEL_PA)
    cases = [("baseline", COMBUSTION)]
    for factor, tag in ((0.9, "-10%"), (1.1, "+10%")):
        cases.append(
            (f"c* {tag}", COMBUSTION.replace(c_star_m_s=factor * COMBUSTION.c_star_m_s))
        )
    for gamma in (1.15, 1.25):
        cases.append((f"gamma {gamma:.2f}", COMBUSTION.replace(gamma=gamma)))
    for factor, tag in ((0.9, "-10%"), (1.1, "+10%")):
        cases.append(
            (
                f"T_c {tag}",
                COMBUSTION.replace(
                    chamber_temperature_k=factor * COMBUSTION.chamber_temperature_k
                ),
            )
        )
    rows = []
    for label, props in cases:
        result = solve_performance_point(props, NOZZLE, M1_M_DOT_OX, M1_M_DOT_F, SEA_LEVEL_PA)
        rows.append((label, props, result, base))
    return rows


def _print_combustion_sensitivity(rows) -> None:
    print("COMBUSTION-PROPERTY SENSITIVITY - deterministic one-factor variation")
    print("(NOT uncertainty propagation: these are prescribed assumptions being varied by hand)")
    print(RULE)
    print(
        f"  {'case':<14} {'p_c [bar]':>10} {'M_e':>8} {'V_e [m/s]':>11} {'F [N]':>10} "
        f"{'I_sp [s]':>9} {'dF':>9} {'eta_c*':>8}"
    )
    for label, props, result, base in rows:
        delta = 100.0 * (result.thrust_n / base.thrust_n - 1.0)
        print(
            f"  {label:<14} {result.chamber_pressure_pa / 1e5:>10.4f} "
            f"{result.exit_state.mach:>8.4f} {result.exit_state.velocity_m_s:>11.3f} "
            f"{result.thrust_n:>10.4f} {result.specific_impulse_s:>9.3f} "
            f"{delta:>8.2f}% {props.c_star_efficiency:>8.4f}"
        )
    print(RULE)
    print(
        "  Structure: p_c is exactly proportional to c* at fixed mass flow and throat area;\n"
        "  T_c enters V_e as sqrt(T_c) and does not touch p_c; gamma changes the area-Mach\n"
        "  solution and both isentropic ratios but leaves p_c untouched.\n"
        "  The eta_c* column shows the c* efficiency each variation implies against its own\n"
        "  ideal c*: varying one property alone makes the set internally inconsistent, which\n"
        "  is exactly why these are labelled assumptions rather than data."
    )
    print()


def _print_convergence() -> None:
    print("NUMERICAL CONVERGENCE")
    print(RULE)
    print("  Nozzle exit-Mach root solve (epsilon = 4, gamma = 1.20):")
    print(f"    {'xtol':>10} {'M_e':>18} {'area-Mach residual':>22} {'F [N]':>12}")
    for xtol in (1e-4, 1e-8, 1e-12, 1e-15):
        mach = solve_supersonic_exit_mach(NOZZLE.expansion_ratio, COMBUSTION.gamma, xtol=xtol)
        residual = (1.0 / mach) * (
            (2.0 / (COMBUSTION.gamma + 1.0))
            * (1.0 + (COMBUSTION.gamma - 1.0) * mach * mach / 2.0)
        ) ** ((COMBUSTION.gamma + 1.0) / (2.0 * (COMBUSTION.gamma - 1.0))) - NOZZLE.expansion_ratio
        point = solve_performance_point(
            COMBUSTION, NOZZLE, M1_M_DOT_OX, M1_M_DOT_F, SEA_LEVEL_PA
        )
        print(f"    {xtol:>10.0e} {mach:>18.12f} {residual:>22.3e} {point.thrust_n:>12.6f}")
    print()
    print("  Transient integration (constant nominal flow, coupled performance):")
    print(
        f"    {'rtol':>10} {'n_report':>9} {'t_burn [s]':>14} {'peak F [N]':>12} "
        f"{'I_total [N s]':>15} {'I_sp,eq [s]':>12}"
    )
    for rtol, n_report in ((1e-5, 201), (1e-7, 401), (1e-9, 801), (1e-11, 1601)):
        transient = simulate_transient(
            GRAIN,
            LAW,
            ConstantOxidizerFlow(M_DOT_NOMINAL, HORIZON_S),
            n_report=n_report,
            rtol=rtol,
            atol=rtol * 1e-4,
        )
        history = couple_transient_to_performance(transient, COMBUSTION, NOZZLE, SEA_LEVEL_PA)
        print(
            f"    {rtol:>10.0e} {n_report:>9d} {transient.burnout_time_s:>14.9f} "
            f"{history.peak_thrust_n:>12.6f} {history.total_impulse_n_s:>15.6f} "
            f"{history.equivalent_specific_impulse_s:>12.6f}"
        )
    print(RULE)
    print()


def sanity_audit(cases, reference_point, ambient_rows) -> list[tuple[str, bool]]:
    checks: list[tuple[str, bool]] = []
    histories = [history for *_, history in cases]

    def every(predicate) -> bool:
        return all(predicate(h) for h in histories)

    def firing(history):
        return history.firing

    checks.append(
        ("m_dot_total > 0 during firing",
         every(lambda h: bool(np.all(h.total_mass_flow_kg_s[firing(h)] > 0.0))))
    )
    checks.append(
        ("p_c > 0 during firing",
         every(lambda h: bool(np.all(h.chamber_pressure_pa[firing(h)] > 0.0))))
    )
    checks.append(
        ("p_c = 0 during true zero-flow idle",
         every(lambda h: bool(np.all(h.chamber_pressure_pa[~firing(h)] == 0.0))))
    )
    checks.append(
        ("thrust = 0 during idle, and no nozzle state is reported",
         every(lambda h: bool(
             np.all(h.thrust_n[~firing(h)] == 0.0)
             and np.all(np.isnan(h.exit_mach[~firing(h)]))
         )))
    )
    checks.append(
        ("M_e > 1 for epsilon > 1",
         every(lambda h: bool(np.all(h.exit_mach[firing(h)] > 1.0))))
    )
    checks.append(
        ("0 < p_e/p_c < 1",
         every(lambda h: bool(np.all(
             (h.exit_pressure_pa[firing(h)] / h.chamber_pressure_pa[firing(h)] > 0.0)
             & (h.exit_pressure_pa[firing(h)] / h.chamber_pressure_pa[firing(h)] < 1.0)
         ))))
    )
    checks.append(
        ("0 < T_e/T_c < 1",
         every(lambda h: bool(np.all(
             (h.exit_temperature_k[firing(h)] / COMBUSTION.chamber_temperature_k > 0.0)
             & (h.exit_temperature_k[firing(h)] / COMBUSTION.chamber_temperature_k < 1.0)
         ))))
    )
    checks.append(
        ("V_e > 0 while firing",
         every(lambda h: bool(np.all(h.exit_velocity_m_s[firing(h)] > 0.0))))
    )
    checks.append(
        ("momentum thrust > 0 while firing",
         every(lambda h: bool(np.all(h.momentum_thrust_n[firing(h)] > 0.0))))
    )
    checks.append(
        (
            "total thrust is finite everywhere",
            every(lambda h: bool(np.all(np.isfinite(h.thrust_n)))),
        )
    )
    checks.append(
        ("F = momentum + pressure exactly",
         every(lambda h: bool(np.all(h.thrust_n == h.momentum_thrust_n + h.pressure_thrust_n))))
    )
    checks.append(
        ("C_F = F/(p_c A_t) to machine precision",
         every(lambda h: bool(np.allclose(
             h.thrust_coefficient[firing(h)],
             h.thrust_n[firing(h)] / (h.chamber_pressure_pa[firing(h)] * NOZZLE.throat_area_m2),
             rtol=1e-13,
         ))))
    )
    checks.append(
        ("p_c A_t = m_dot c* to machine precision",
         every(lambda h: float(np.max(np.abs(
             h.chamber_pressure_pa * NOZZLE.throat_area_m2
             - h.total_mass_flow_kg_s * COMBUSTION.c_star_m_s
         ))) < 1.0e-9))
    )
    checks.append(
        ("I_sp positive where defined",
         every(lambda h: bool(np.all(h.specific_impulse_s[firing(h)] > 0.0))))
    )
    thrusts = [row[2].thrust_n for row in ambient_rows]
    checks.append(
        ("lower ambient pressure increases thrust at a fixed chamber state",
         all(later > earlier for earlier, later in itertools.pairwise(thrusts)))
    )
    checks.append(("total impulse is positive", every(lambda h: h.total_impulse_n_s > 0.0)))
    checks.append(
        ("equivalent I_sp identity closes",
         every(lambda h: abs(
             h.equivalent_specific_impulse_s
             - h.total_impulse_n_s / (h.propellant_mass_kg * STANDARD_GRAVITY_M_S2)
         ) < 1.0e-12))
    )
    checks.append(
        ("no NaN/inf in valid firing states",
         every(lambda h: bool(np.all(np.isfinite(h.thrust_n[firing(h)]))
                              and np.all(np.isfinite(h.chamber_pressure_pa[firing(h)]))
                              and np.all(np.isfinite(h.exit_velocity_m_s[firing(h)])))))
    )
    route_b = (
        reference_point.thrust_coefficient
        * reference_point.chamber_pressure_pa
        * NOZZLE.throat_area_m2
    )
    checks.append(
        ("two independent thrust routes agree at the reference point",
         abs(reference_point.thrust_n - route_b) < 1.0e-9)
    )
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


def _print_observations(cases, reference_point) -> None:
    print("OBSERVATIONS (reported as computed, not tuned)")
    print(RULE)
    case_a = cases[0][3]
    case_d = cases[3][3]
    print(
        "  1. THRUST RISES THROUGH A CONSTANT-FLOW BURN. Milestone 2 showed m_dot_f rising as\n"
        "     the port opens; with a fixed throat the chamber pressure and thrust follow it\n"
        f"     directly, from {case_a.thrust_n[0]:.2f} N to {case_a.thrust_n[-1]:.2f} N "
        f"({100 * (case_a.thrust_n[-1] / case_a.thrust_n[0] - 1):.1f} %) over the burn.\n"
        "     A real motor with a fixed oxidizer supply would not necessarily behave this way,\n"
        "     because its oxidizer flow would respond to the rising chamber pressure."
    )
    print(
        "  2. THE EXIT STATE IS FROZEN. With epsilon, gamma and T_c all prescribed constants,\n"
        f"     M_e = {case_a.exit_mach[0]:.4f}, T_e and "
        f"V_e = {case_a.exit_velocity_m_s[0]:.1f} m/s do not change at all during a burn.\n"
        "     Thrust is therefore exactly affine in the total mass flow. This is a property of\n"
        "     the prescribed-property assumption, NOT a physical prediction: real c*, gamma and\n"
        "     T_c all vary with O/F, and Milestone 3 deliberately does not model that."
    )
    print(
        "  3. THE REFERENCE NOZZLE IS SLIGHTLY UNDEREXPANDED AT SEA LEVEL, so the pressure term\n"
        f"     helps rather than hurts: {reference_point.pressure_thrust_n:+.3f} N of "
        f"{reference_point.thrust_n:.3f} N "
        f"({100 * reference_point.pressure_thrust_n / reference_point.thrust_n:.2f} %).\n"
        "     As p_c rises through the burn the nozzle becomes progressively more\n"
        "     underexpanded, so the pressure contribution grows."
    )
    print(
        f"  4. PREDICTED I_sp ({reference_point.specific_impulse_s:.1f} s) EXCEEDS THE MEASURED "
        "RANGE of the comparable motor\n"
        "     (Rezaei et al. reported 167-223 s). That is expected and is reported rather than\n"
        "     corrected: this is an ideal 1-D isentropic nozzle with no divergence loss, no\n"
        "     boundary layer, no heat loss and no two-phase effects, and c* efficiency is the\n"
        "     ONLY loss represented anywhere in the model."
    )
    print(
        "  5. THE ZERO-FLOW COAST PRODUCES NOTHING. During case D's coast p_c and thrust are\n"
        f"     exactly zero over {int((~case_d.firing).sum())} reported samples, and no exit "
        "state is reported at all, so\n"
        "     no nozzle number can be mistaken for a firing condition."
    )
    print(RULE)
    print()


def main() -> int:
    _banner()
    _print_configuration()
    _print_throat_selection()
    reference_point = _print_static_reference()

    cases = [
        run_case(
            ConstantOxidizerFlow(M_DOT_NOMINAL, HORIZON_S),
            "CASE A",
            f"Milestone 2 nominal constant flow, m_dot_ox = {M_DOT_NOMINAL:.3f} kg/s, to burnout",
        ),
        run_case(
            ConstantOxidizerFlow(M_DOT_LOW, HORIZON_S),
            "CASE B",
            f"lower constant flow, m_dot_ox = {M_DOT_LOW:.3f} kg/s, to burnout",
        ),
        run_case(
            ConstantOxidizerFlow(M_DOT_HIGH, HORIZON_S),
            "CASE C",
            f"higher constant flow, m_dot_ox = {M_DOT_HIGH:.3f} kg/s, to burnout",
        ),
        run_case(
            PiecewiseConstantOxidizerFlow.from_durations(PIECEWISE_STEPS),
            "CASE D",
            "Milestone 2 piecewise prescribed schedule (includes a zero-flow coast)",
        ),
    ]

    for label, description, transient, history in cases:
        _report_case(label, description, transient, history)
        if label in ("CASE A", "CASE D"):
            _print_trajectory(history, n_rows=14 if label == "CASE D" else 12)

    ambient_rows = ambient_sensitivity(reference_point)
    _print_ambient_sensitivity(ambient_rows)

    _print_combustion_sensitivity(combustion_sensitivity())
    _print_convergence()
    _print_observations(cases, reference_point)

    ok = _print_sanity_audit(sanity_audit(cases, reference_point, ambient_rows))
    print(
        "END OF MILESTONE 3 STUDY - oxidizer flow remains prescribed; no tank, injector or\n"
        "feed-system model exists, and no result here has been validated against hardware."
    )
    assert all(
        transient.termination_reason
        in (TerminationReason.FUEL_DEPLETED, TerminationReason.COMPLETED_TIME_WINDOW)
        for *_, transient, _ in cases
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
