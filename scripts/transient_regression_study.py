"""Milestone 2 transient regression study under prescribed oxidizer-flow histories.

Generic, reduced-order, educational study of an illustrative N2O/HTPB hybrid
grain.  The oxidizer mass flow is a *prescribed input* throughout: nothing here
is derived from tank state, vapour pressure, blowdown, feed lines, injectors or
chamber pressure, and no chamber pressure, c*, nozzle flow or thrust is
computed.

Run with::

    python scripts/transient_regression_study.py
"""

from __future__ import annotations

import math

import numpy as np

from hybrid_rocket_motor import (
    REPRESENTATIVE_GRAIN,
    REZAEI_2018_N2O_HTPB,
)
from hybrid_rocket_motor.transient import (
    ConstantOxidizerFlow,
    FluxRangeStatus,
    PiecewiseConstantOxidizerFlow,
    TerminationReason,
    analytic_constant_flow_burnout_time,
    analytic_constant_flow_radius,
    simulate_transient,
)

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
R0 = GRAIN.initial_port_diameter_m / 2.0
R_OUTER = GRAIN.outer_diameter_m / 2.0

#: Representative prescribed oxidizer mass flow [kg/s] (same as Milestone 1).
M_DOT_NOMINAL = 0.100
M_DOT_LOW = 0.050
M_DOT_HIGH = 0.150

#: Illustrative piecewise schedule: nominal, a zero-flow coast that exercises the
#: undefined-O/F convention, a higher flow, then back to nominal.  It is not a
#: valve, controller or feed-system design, and it is not tuned to shape O/F.
PIECEWISE_STEPS: tuple[tuple[float, float], ...] = (
    (10.0, M_DOT_NOMINAL),
    (6.0, 0.0),
    (10.0, M_DOT_HIGH),
    (14.0, M_DOT_NOMINAL),
)

RULE = "-" * 100
HORIZON_S = 200.0


def _banner() -> None:
    print("=" * 100)
    print("MILESTONE 2 - TRANSIENT PORT REGRESSION UNDER PRESCRIBED OXIDIZER-FLOW HISTORIES")
    print("Generic N2O/HTPB hybrid rocket motor - reduced-order educational model")
    print("=" * 100)
    print()
    print("SCOPE NOTICE")
    print("  Generic illustrative geometry - not a real motor, not a manufacturing drawing.")
    print("  Oxidizer mass flow is PRESCRIBED INPUT ONLY.")
    print("  This milestone does NOT model chamber pressure, c*, nozzle flow or thrust.")
    print()
    print("STATE EQUATION")
    print("  dr_p/dt = r_dot          (equivalently  dD_p/dt = 2 r_dot)")
    print("  G_ox(t) = m_dot_ox(t) / (pi r_p(t)^2)")
    print(f"  r_dot   = {LAW.si_equation}")
    lo, hi = LAW.valid_flux_range_si
    print(f"  source data range: G_ox = {lo:.1f} to {hi:.1f} kg/(m^2 s)")
    print()
    print("GRAIN (illustrative)")
    print(
        f"  L = {GRAIN.length_m:.3f} m   D_p,0 = {GRAIN.initial_port_diameter_m:.3f} m   "
        f"D_o = {GRAIN.outer_diameter_m:.3f} m   rho_f = {GRAIN.fuel_density_kg_m3:.0f} kg/m^3"
    )
    print(
        f"  initial web = {GRAIN.initial_web_thickness_m * 1e3:.1f} mm   "
        f"loaded fuel mass = {GRAIN.initial_fuel_mass_kg:.6f} kg"
    )
    print()


def _status_fractions(result) -> dict[FluxRangeStatus, float]:
    return {status: result.time_fraction_with_status(status) for status in FluxRangeStatus}


def _report_case(name: str, description: str, result, *, analytic_burnout: float | None) -> None:
    print(f"{name} - {description}")
    print(RULE)
    first, last = 0, result.final_index
    print(
        f"  {'quantity':<34} {'initial':>16} {'final':>16}   units"
    )
    rows = (
        ("port diameter D_p", result.port_diameter_m * 1e3, "mm"),
        ("remaining web", result.web_thickness_m * 1e3, "mm"),
        ("oxidizer mass flow m_dot_ox", result.oxidizer_mass_flow_kg_s, "kg/s"),
        ("oxidizer mass flux G_ox", result.oxidizer_mass_flux_si, "kg/(m^2 s)"),
        ("regression rate r_dot", result.regression_rate_mm_s, "mm/s"),
        ("burning area A_burn", result.burning_area_m2, "m^2"),
        ("fuel mass flow m_dot_f", result.fuel_mass_flow_kg_s, "kg/s"),
        ("mixture ratio O/F", result.mixture_ratio, "-"),
        ("remaining fuel mass", result.remaining_fuel_mass_kg, "kg"),
    )
    for label, array, unit in rows:
        print(f"  {label:<34} {array[first]:>16.6f} {array[last]:>16.6f}   {unit}")
    print(RULE)
    print(f"  simulated duration                 {result.time_s[-1]:.4f} s")
    print(f"  termination                        {result.termination_reason.value}")
    if result.burnout_time_s is not None:
        print(f"  burnout time (numerical)           {result.burnout_time_s:.6f} s")
    if analytic_burnout is not None and math.isfinite(analytic_burnout):
        print(f"  burnout time (closed form)         {analytic_burnout:.6f} s")
        if result.burnout_time_s is not None:
            error = abs(result.burnout_time_s - analytic_burnout) / analytic_burnout
            print(f"  relative difference                {error:.3e}")
    print(f"  oxidizer consumed                  {result.cumulative_oxidizer_mass_kg[-1]:.6f} kg")
    print(f"  fuel consumed (integrated)         {result.cumulative_fuel_mass_kg[-1]:.6f} kg")
    print(
        f"  fuel consumed (from geometry)      "
        f"{result.geometric_fuel_consumed_kg[-1]:.6f} kg"
    )
    print(
        f"  mass-closure residual              "
        f"{result.max_absolute_fuel_closure_residual_kg:.3e} kg (abs), "
        f"{result.max_relative_fuel_closure_residual:.3e} (rel)"
    )
    fractions = _status_fractions(result)
    print("  time spent by source-flux status:")
    for status, fraction in fractions.items():
        if fraction > 0.0:
            print(f"      {status.value:<28} {100.0 * fraction:6.2f} %")
    within = fractions[FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE]
    if result.termination_reason is TerminationReason.FUEL_DEPLETED and within < 0.5:
        print(
            f"  *** Burnout prediction relies materially on EXTRAPOLATION: only "
            f"{100.0 * within:.1f} % of the burn"
        )
        print("      lies inside the correlation's measured flux range. ***")
    print(RULE)
    print()


def _print_trajectory(result, n_rows: int = 12) -> None:
    print("  sampled trajectory")
    print(
        f"    {'t':>8} {'D_p':>9} {'m_dot_ox':>10} {'G_ox':>11} {'r_dot':>9} "
        f"{'m_dot_f':>10} {'O/F':>8} {'fuel left':>10}  status"
    )
    print(
        f"    {'[s]':>8} {'[mm]':>9} {'[kg/s]':>10} {'[kg/m^2 s]':>11} {'[mm/s]':>9} "
        f"{'[kg/s]':>10} {'[-]':>8} {'[kg]':>10}"
    )
    indices = np.unique(np.linspace(0, result.final_index, n_rows).astype(int))
    for index in indices:
        of_value = result.mixture_ratio[index]
        of_text = "undefined" if math.isnan(of_value) else f"{of_value:8.3f}"
        print(
            f"    {result.time_s[index]:>8.3f} {result.port_diameter_m[index] * 1e3:>9.4f} "
            f"{result.oxidizer_mass_flow_kg_s[index]:>10.4f} "
            f"{result.oxidizer_mass_flux_si[index]:>11.3f} "
            f"{result.regression_rate_mm_s[index]:>9.5f} "
            f"{result.fuel_mass_flow_kg_s[index]:>10.6f} {of_text:>8} "
            f"{result.remaining_fuel_mass_kg[index]:>10.6f}  "
            f"{result.flux_status[index].value}"
        )
    print()


def run_constant_case(m_dot: float):
    result = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(m_dot, HORIZON_S), n_report=801
    )
    analytic = analytic_constant_flow_burnout_time(LAW, R0, R_OUTER, m_dot)
    return result, analytic


def run_piecewise_case():
    schedule = PiecewiseConstantOxidizerFlow.from_durations(PIECEWISE_STEPS)
    return simulate_transient(GRAIN, LAW, schedule, n_report=801)


def _convergence_study() -> None:
    print("NUMERICAL CONVERGENCE - representative constant-flow case")
    print(
        "  Reference is the closed-form solution, not another numerical run: "
        "t_burn = [r_o^(2n+1) - r_0^(2n+1)] / [(2n+1) K]."
    )
    print(RULE)
    analytic_t = analytic_constant_flow_burnout_time(LAW, R0, R_OUTER, M_DOT_NOMINAL)
    analytic_r10 = analytic_constant_flow_radius(LAW, R0, M_DOT_NOMINAL, 10.0)
    print(f"  closed-form burnout time        = {analytic_t:.9f} s")
    print(f"  closed-form radius at t = 10 s  = {analytic_r10:.12f} m")
    print(RULE)
    print(
        f"  {'rtol':>9} {'method':>8} {'t_burn [s]':>16} {'rel. err':>11} "
        f"{'r(10 s) rel. err':>17} {'fuel [kg]':>11} {'closure [kg]':>13}"
    )
    for rtol, method in ((1e-5, "DOP853"), (1e-7, "DOP853"), (1e-9, "DOP853"),
                         (1e-11, "DOP853"), (1e-11, "RK45")):
        run = simulate_transient(
            GRAIN,
            LAW,
            ConstantOxidizerFlow(M_DOT_NOMINAL, HORIZON_S),
            n_report=201,
            rtol=rtol,
            atol=rtol * 1.0e-4,
            method=method,
        )
        t_error = abs(run.burnout_time_s - analytic_t) / analytic_t
        short = simulate_transient(
            GRAIN,
            LAW,
            ConstantOxidizerFlow(M_DOT_NOMINAL, 10.0),
            n_report=11,
            rtol=rtol,
            atol=rtol * 1.0e-4,
            method=method,
        )
        r_error = abs(short.port_radius_m[-1] - analytic_r10) / analytic_r10
        print(
            f"  {rtol:>9.0e} {method:>8} {run.burnout_time_s:>16.9f} {t_error:>11.2e} "
            f"{r_error:>17.2e} {run.cumulative_fuel_mass_kg[-1]:>11.6f} "
            f"{run.max_absolute_fuel_closure_residual_kg:>13.3e}"
        )
    print(RULE)
    print(
        "  Reporting-grid independence: burnout time and consumed fuel are unchanged when the\n"
        "  reporting grid is refined, because the grid only samples the dense solution."
    )
    coarse = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M_DOT_NOMINAL, HORIZON_S), n_report=11
    )
    fine = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M_DOT_NOMINAL, HORIZON_S), n_report=4001
    )
    print(
        f"    n_report =   11 -> t_burn = {coarse.burnout_time_s:.9f} s, "
        f"fuel = {coarse.cumulative_fuel_mass_kg[-1]:.9f} kg"
    )
    print(
        f"    n_report = 4001 -> t_burn = {fine.burnout_time_s:.9f} s, "
        f"fuel = {fine.cumulative_fuel_mass_kg[-1]:.9f} kg"
    )
    print(RULE)
    print()


def sanity_audit(results) -> list[tuple[str, bool]]:
    """Explicit engineering sanity checks across every simulated case."""
    checks: list[tuple[str, bool]] = []

    def every(predicate) -> bool:
        return all(predicate(result) for result in results)

    checks.append(("port radius is nondecreasing",
                   every(lambda r: bool(np.all(np.diff(r.port_radius_m) >= 0.0)))))
    checks.append(("port diameter is nondecreasing",
                   every(lambda r: bool(np.all(np.diff(r.port_diameter_m) >= 0.0)))))
    checks.append(("remaining fuel mass is nonincreasing",
                   every(lambda r: bool(np.all(np.diff(r.remaining_fuel_mass_kg) <= 0.0)))))
    checks.append(("remaining fuel mass is never materially negative",
                   every(lambda r: bool(np.all(r.remaining_fuel_mass_kg >= -1.0e-12)))))
    checks.append(("port never exceeds the outer grain radius",
                   every(lambda r: bool(np.all(r.port_radius_m <= R_OUTER + 1.0e-15)))))

    constant_results = results[:3]
    checks.append((
        "constant flow: G_ox decreases as the port grows",
        all(bool(np.all(np.diff(r.oxidizer_mass_flux_si) < 0.0)) for r in constant_results),
    ))
    checks.append((
        "constant flow: r_dot decreases as the port grows",
        all(bool(np.all(np.diff(r.regression_rate_m_s) < 0.0)) for r in constant_results),
    ))

    def positive_flow_fuel(result) -> bool:
        mask = result.oxidizer_mass_flow_kg_s > 0.0
        return bool(np.all(result.fuel_mass_flow_kg_s[mask] > 0.0))

    checks.append(("m_dot_f > 0 wherever oxidizer flow is positive", every(positive_flow_fuel)))

    def of_defined(result) -> bool:
        mask = result.fuel_mass_flow_kg_s > 0.0
        values = result.mixture_ratio[mask]
        return bool(np.all(np.isfinite(values)) and np.all(values > 0.0))

    checks.append(("O/F is positive and finite where defined", every(of_defined)))

    def zero_flow_frozen(result) -> bool:
        mask = result.oxidizer_mass_flow_kg_s == 0.0
        if not np.any(mask):
            return True
        radii = result.port_radius_m[mask]
        return bool(
            np.all(result.regression_rate_m_s[mask] == 0.0)
            and np.all(np.isnan(result.mixture_ratio[mask]))
            and float(np.max(radii) - np.min(radii)) <= 1.0e-15
        )

    checks.append(("zero-flow intervals cause no regression and leave O/F undefined",
                   every(zero_flow_frozen)))

    def burnout_at_boundary(result) -> bool:
        if result.termination_reason is not TerminationReason.FUEL_DEPLETED:
            return True
        return (
            abs(result.port_radius_m[-1] - R_OUTER) < 1.0e-9
            and abs(result.remaining_fuel_mass_kg[-1]) < 1.0e-9
        )

    checks.append(("burnout occurs at the geometric outer boundary", every(burnout_at_boundary)))

    analytic_ok = True
    for result, m_dot in zip(constant_results, (M_DOT_NOMINAL, M_DOT_LOW, M_DOT_HIGH), strict=True):
        expected = analytic_constant_flow_burnout_time(LAW, R0, R_OUTER, m_dot)
        analytic_ok = analytic_ok and abs(result.burnout_time_s - expected) / expected < 1.0e-8
    checks.append(("numerical burnout agrees with the closed-form solution", analytic_ok))

    checks.append((
        "fuel-mass conservation closes to better than 1e-8 kg",
        every(lambda r: r.max_absolute_fuel_closure_residual_kg < 1.0e-8),
    ))

    ox_ok = all(
        abs(r.cumulative_oxidizer_mass_kg[-1] - m_dot * r.time_s[-1]) < 1.0e-9
        for r, m_dot in zip(constant_results, (M_DOT_NOMINAL, M_DOT_LOW, M_DOT_HIGH), strict=True)
    )
    checks.append(("constant-flow oxidizer-mass conservation closes", ox_ok))

    def no_bad_values(result) -> bool:
        mask = result.oxidizer_mass_flow_kg_s > 0.0
        arrays = (
            result.port_radius_m,
            result.oxidizer_mass_flux_si,
            result.regression_rate_m_s,
            result.fuel_mass_flow_kg_s,
            result.mixture_ratio,
            result.cumulative_fuel_mass_kg,
            result.cumulative_oxidizer_mass_kg,
        )
        return all(bool(np.all(np.isfinite(array[mask]))) for array in arrays)

    checks.append(("no NaN/inf in positive-flow states", every(no_bad_values)))

    def classification_honest(result) -> bool:
        low, high = LAW.valid_flux_range_si
        for flux, flow, status in zip(
            result.oxidizer_mass_flux_si,
            result.oxidizer_mass_flow_kg_s,
            result.flux_status,
            strict=True,
        ):
            if flow == 0.0:
                expected = FluxRangeStatus.ZERO_OXIDIZER_FLOW
            elif flux < low:
                expected = FluxRangeStatus.BELOW_SOURCE_FLUX_RANGE
            elif flux > high:
                expected = FluxRangeStatus.ABOVE_SOURCE_FLUX_RANGE
            else:
                expected = FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE
            if status is not expected:
                return False
        return True

    checks.append(("validity/extrapolation classification is honest", every(classification_honest)))
    return checks


def _print_sanity_audit(results) -> bool:
    print("ENGINEERING SANITY AUDIT")
    print(RULE)
    checks = sanity_audit(results)
    for description, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {description}")
    print(RULE)
    all_ok = all(ok for _, ok in checks)
    print(f"  {'ALL SANITY CHECKS PASSED' if all_ok else '*** SANITY CHECKS FAILED ***'}")
    print()
    return all_ok


def _print_observations(case_a, case_b, case_c, case_d) -> None:
    print("OBSERVATIONS (reported as computed, not tuned)")
    print(RULE)
    within_a = case_a.time_fraction_with_status(FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE)
    within_b = case_b.time_fraction_with_status(FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE)
    within_c = case_c.time_fraction_with_status(FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE)
    print(
        "  1. EXTRAPOLATION DOMINATES THE LONG BURNS. At constant oxidizer flow the port area\n"
        "     grows, so G_ox falls monotonically and eventually leaves the correlation's measured\n"
        f"     range. Only {100 * within_a:.1f} % (case A), {100 * within_b:.1f} % (case B) and "
        f"{100 * within_c:.1f} % (case C) of each burn lies inside\n"
        "     35-120 kg/(m^2 s). The predicted burnout times therefore rest mostly on\n"
        "     extrapolation below the data. This is reported, not clamped or discarded."
    )
    print(
        "  2. LOWER FLOW EXTRAPOLATES SOONER, NOT LATER. Case B starts at "
        f"G_ox = {case_b.oxidizer_mass_flux_si[0]:.1f} kg/(m^2 s),\n"
        "     barely inside the band, and leaves it almost immediately; the higher-flow case C\n"
        "     stays inside longest. Lower flow is not the 'safer' choice for model validity here."
    )
    ratio_a = case_a.fuel_mass_flow_kg_s[-1] / case_a.fuel_mass_flow_kg_s[0]
    print(
        "  3. FUEL FLOW RISES WHILE REGRESSION RATE FALLS. Over case A the regression rate drops\n"
        f"     from {case_a.regression_rate_mm_s[0]:.4f} to "
        f"{case_a.regression_rate_mm_s[-1]:.4f} mm/s, yet m_dot_f *rises* by a factor "
        f"{ratio_a:.3f},\n"
        "     because the burning area grows as r_p while r_dot falls only as r_p^(-2n) with\n"
        "     n = 0.3667. The net exponent on r_p is 1 - 2n = +0.267. O/F therefore drifts\n"
        f"     downward from {case_a.mixture_ratio[0]:.3f} to {case_a.mixture_ratio[-1]:.3f} "
        "rather than upward."
    )
    print(
        "  4. THE O/F EXCURSION IS MILD. Across the whole case-A burn O/F moves by only about\n"
        f"     {100 * (1 - case_a.mixture_ratio[-1] / case_a.mixture_ratio[0]):.1f} %. That is a "
        "direct consequence of the low sourced flux exponent; a\n"
        "     steeper exponent would give a much larger drift. Milestone 1 documents that\n"
        "     n = 0.3667 sits below the 0.5-0.8 band reported for most hybrid systems."
    )
    coast_mask = case_d.oxidizer_mass_flow_kg_s == 0.0
    print(
        "  5. THE ZERO-FLOW COAST IS EXACTLY FROZEN. During case D's coast the port radius is\n"
        f"     unchanged to machine precision over {int(coast_mask.sum())} reported samples, "
        "r_dot and m_dot_f are exactly\n"
        "     zero, and O/F is reported as undefined (NaN) rather than fabricated."
    )
    worst_closure = max(
        r.max_absolute_fuel_closure_residual_kg for r in (case_a, case_b, case_c, case_d)
    )
    print(
        "  6. MASS CLOSURE IS TIGHT. The integrated fuel mass and the fuel mass inferred from the\n"
        f"     port geometry agree to {worst_closure:.2e} kg worst case across all four cases, "
        "on a\n"
        f"     {GRAIN.initial_fuel_mass_kg:.3f} kg grain."
    )
    print(RULE)
    print()


def main() -> int:
    _banner()

    case_a, analytic_a = run_constant_case(M_DOT_NOMINAL)
    _report_case(
        "CASE A",
        f"constant representative flow, m_dot_ox = {M_DOT_NOMINAL:.3f} kg/s, run to burnout",
        case_a,
        analytic_burnout=analytic_a,
    )
    _print_trajectory(case_a)

    case_b, analytic_b = run_constant_case(M_DOT_LOW)
    _report_case(
        "CASE B",
        f"lower constant flow, m_dot_ox = {M_DOT_LOW:.3f} kg/s, run to burnout",
        case_b,
        analytic_burnout=analytic_b,
    )

    case_c, analytic_c = run_constant_case(M_DOT_HIGH)
    _report_case(
        "CASE C",
        f"higher constant flow, m_dot_ox = {M_DOT_HIGH:.3f} kg/s, run to burnout",
        case_c,
        analytic_burnout=analytic_c,
    )

    case_d = run_piecewise_case()
    schedule_text = ", ".join(
        f"{duration:g} s at {flow:.3f} kg/s" for duration, flow in PIECEWISE_STEPS
    )
    _report_case(
        "CASE D",
        f"piecewise prescribed schedule ({schedule_text})",
        case_d,
        analytic_burnout=None,
    )
    _print_trajectory(case_d, n_rows=16)

    _convergence_study()
    _print_observations(case_a, case_b, case_c, case_d)

    ok = _print_sanity_audit([case_a, case_b, case_c, case_d])
    print(
        "END OF MILESTONE 2 STUDY - no chamber pressure, nozzle or thrust result is produced."
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
