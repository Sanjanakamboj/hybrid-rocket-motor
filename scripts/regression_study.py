"""Milestone 1 prescribed-oxidizer-flow regression study.

Generic, reduced-order, educational study of an illustrative N2O/HTPB hybrid
grain.  Oxidizer mass flow is a *prescribed input*: nothing here is derived from
tank thermodynamics, injector flow or chamber pressure, and no chamber pressure,
c*, nozzle flow or thrust is computed at this milestone.

Run with::

    python scripts/regression_study.py
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from hybrid_rocket_motor import (
    ILLUSTRATIVE_ANCHOR_FLUX_SI,
    ILLUSTRATIVE_SENSITIVITY_LAWS,
    REPRESENTATIVE_GRAIN,
    REZAEI_2018_N2O_HTPB,
    OperatingPoint,
)

# --- study configuration (shared with the figure script) -------------------

#: Prescribed oxidizer mass flows [kg/s].  Chosen so that the resulting fluxes
#: span the range over which the sourced correlation reports data; this is a
#: study range, not a claimed hardware operating envelope.
OXIDIZER_MASS_FLOWS_KG_S: tuple[float, ...] = tuple(
    round(float(v), 6) for v in np.linspace(0.045, 0.150, 8)
)

#: Representative prescribed oxidizer mass flow [kg/s].
REPRESENTATIVE_M_DOT_OX_KG_S = 0.100

#: Hypothetical instantaneous port diameters for the port-size sensitivity [m].
#: These are *instantaneous* what-if sizes, not a burn history: no time
#: integration of port growth happens at Milestone 1.
PORT_DIAMETERS_M: tuple[float, ...] = (0.040, 0.050, 0.060, 0.070)

SOURCED_LAW = REZAEI_2018_N2O_HTPB
GRAIN = REPRESENTATIVE_GRAIN

RULE = "-" * 96


def _banner() -> None:
    print("=" * 96)
    print("MILESTONE 1 - PRESCRIBED-OXIDIZER-FLOW REGRESSION STUDY")
    print("Generic N2O/HTPB hybrid rocket motor - reduced-order educational model")
    print("=" * 96)
    print()
    print("SCOPE NOTICE")
    print("  Generic illustrative geometry - not a real motor, not a manufacturing drawing.")
    print("  Oxidizer mass flow is PRESCRIBED INPUT ONLY.")
    print("  This milestone does NOT model chamber pressure, c*, nozzle flow or thrust.")
    print()


def _print_configuration() -> None:
    print("GRAIN GEOMETRY (illustrative)")
    print(f"  grain length            L        = {GRAIN.length_m:.4f} m")
    print(f"  initial port diameter   D_p,0    = {GRAIN.initial_port_diameter_m:.4f} m")
    print(f"  outer grain diameter    D_o      = {GRAIN.outer_diameter_m:.4f} m")
    print(f"  fuel density            rho_f    = {GRAIN.fuel_density_kg_m3:.1f} kg/m^3")
    print(f"  initial web thickness            = {GRAIN.initial_web_thickness_m:.4f} m")
    print(f"  initial port area       A_port,0 = {GRAIN.initial_port_area_m2:.6e} m^2")
    print(f"  initial burning area    A_burn,0 = {GRAIN.initial_burning_area_m2:.6e} m^2")
    print(f"  loaded fuel mass        m_f,0    = {GRAIN.initial_fuel_mass_kg:.6f} kg")
    print()
    print("REGRESSION LAW (SOURCED)")
    print(f"  as published : {SOURCED_LAW.source_equation}")
    print(f"  converted    : {SOURCED_LAW.si_equation}")
    lo, hi = SOURCED_LAW.valid_flux_range_si
    print(f"  source data range : G_ox = {lo:.1f} to {hi:.1f} kg/(m^2 s)")
    print(f"  reference    : {SOURCED_LAW.reference}")
    print()


def _flag(law, flux_si: float) -> str:
    return "in-range" if law.is_within_validated_flux_range(flux_si) else "EXTRAPOLATED"


def study_oxidizer_flow_sweep() -> list[dict[str, float | str]]:
    """Sweep prescribed oxidizer mass flow at the initial port diameter."""
    rows: list[dict[str, float | str]] = []
    for m_dot_ox in OXIDIZER_MASS_FLOWS_KG_S:
        point = OperatingPoint(GRAIN, m_dot_ox, GRAIN.initial_port_diameter_m)
        flux = point.oxidizer_mass_flux_si
        rows.append(
            {
                "m_dot_ox_kg_s": m_dot_ox,
                "G_ox_si": flux,
                "G_ox_g_cm2_s": flux / 10.0,
                "r_dot_mm_s": point.regression_rate_si(SOURCED_LAW) * 1.0e3,
                "m_dot_f_kg_s": point.fuel_mass_flow_kg_s(SOURCED_LAW),
                "of_ratio": point.mixture_ratio(SOURCED_LAW),
                "range_flag": _flag(SOURCED_LAW, flux),
            }
        )
    return rows


def study_port_diameter_sensitivity(
    m_dot_ox: float = REPRESENTATIVE_M_DOT_OX_KG_S,
) -> list[dict[str, float | str]]:
    """Instantaneous port-size sensitivity at a single fixed oxidizer mass flow."""
    rows: list[dict[str, float | str]] = []
    for diameter in PORT_DIAMETERS_M:
        point = OperatingPoint(GRAIN, m_dot_ox, diameter)
        flux = point.oxidizer_mass_flux_si
        rows.append(
            {
                "d_port_m": diameter,
                "A_port_m2": point.port_area_m2,
                "G_ox_si": flux,
                "r_dot_mm_s": point.regression_rate_si(SOURCED_LAW) * 1.0e3,
                "A_burn_m2": point.burning_area_m2,
                "m_dot_f_kg_s": point.fuel_mass_flow_kg_s(SOURCED_LAW),
                "of_ratio": point.mixture_ratio(SOURCED_LAW),
                "fuel_left_kg": point.remaining_fuel_mass_kg,
                "range_flag": _flag(SOURCED_LAW, flux),
            }
        )
    return rows


def _print_flow_sweep(rows: list[dict[str, float | str]]) -> None:
    print("TABLE 1 - PRESCRIBED OXIDIZER MASS FLOW SWEEP (at initial port D_p = "
          f"{GRAIN.initial_port_diameter_m * 1e3:.0f} mm)")
    print(RULE)
    print(
        f"{'m_dot_ox':>10} {'G_ox':>12} {'G_ox':>12} {'r_dot':>10} "
        f"{'m_dot_f':>11} {'O/F':>8}  {'source range':<13}"
    )
    print(
        f"{'[kg/s]':>10} {'[kg/m^2 s]':>12} {'[g/cm^2 s]':>12} {'[mm/s]':>10} "
        f"{'[kg/s]':>11} {'[-]':>8}  {'':<13}"
    )
    print(RULE)
    for r in rows:
        print(
            f"{r['m_dot_ox_kg_s']:>10.4f} {r['G_ox_si']:>12.3f} {r['G_ox_g_cm2_s']:>12.4f} "
            f"{r['r_dot_mm_s']:>10.4f} {r['m_dot_f_kg_s']:>11.6f} {r['of_ratio']:>8.3f}  "
            f"{r['range_flag']:<13}"
        )
    print(RULE)
    print()


def _print_port_sensitivity(rows: list[dict[str, float | str]], m_dot_ox: float) -> None:
    print(
        "TABLE 2 - INSTANTANEOUS PORT-DIAMETER SENSITIVITY at fixed "
        f"m_dot_ox = {m_dot_ox:.4f} kg/s"
    )
    print("         (what-if port sizes at one instant; NOT a burn-time integration)")
    print(RULE)
    print(
        f"{'D_p':>7} {'A_port':>12} {'G_ox':>11} {'r_dot':>10} {'A_burn':>10} "
        f"{'m_dot_f':>11} {'O/F':>8}  {'source range':<13}"
    )
    print(
        f"{'[mm]':>7} {'[m^2]':>12} {'[kg/m^2 s]':>11} {'[mm/s]':>10} {'[m^2]':>10} "
        f"{'[kg/s]':>11} {'[-]':>8}  {'':<13}"
    )
    print(RULE)
    for r in rows:
        print(
            f"{r['d_port_m'] * 1e3:>7.1f} {r['A_port_m2']:>12.6e} {r['G_ox_si']:>11.3f} "
            f"{r['r_dot_mm_s']:>10.4f} {r['A_burn_m2']:>10.6f} {r['m_dot_f_kg_s']:>11.6f} "
            f"{r['of_ratio']:>8.3f}  {r['range_flag']:<13}"
        )
    print(RULE)
    first, last = rows[0], rows[-1]
    print(
        f"  Coupling check: D_p {first['d_port_m'] * 1e3:.0f} -> {last['d_port_m'] * 1e3:.0f} mm "
        f"raises A_port {first['A_port_m2']:.3e} -> {last['A_port_m2']:.3e} m^2, "
        f"lowers G_ox {first['G_ox_si']:.1f} -> {last['G_ox_si']:.1f} kg/(m^2 s), "
        f"lowers r_dot {first['r_dot_mm_s']:.4f} -> {last['r_dot_mm_s']:.4f} mm/s."
    )
    print(
        "  Note: m_dot_f does NOT fall in step with r_dot, because the burning area grows "
        "with D_p at the same time."
    )
    print()


def _print_exponent_sensitivity(m_dot_ox: float) -> None:
    point = OperatingPoint(GRAIN, m_dot_ox, GRAIN.initial_port_diameter_m)
    flux = point.oxidizer_mass_flux_si
    print(
        "TABLE 3 - FLUX-EXPONENT SENSITIVITY  "
        "***ILLUSTRATIVE COEFFICIENTS - NOT MEASURED DATA***"
    )
    print(
        f"         at m_dot_ox = {m_dot_ox:.4f} kg/s, "
        f"D_p = {GRAIN.initial_port_diameter_m * 1e3:.0f} mm, "
        f"G_ox = {flux:.3f} kg/(m^2 s)"
    )
    print(
        "         Illustrative laws share the sourced law's value at the disclosed anchor "
        f"G_ox = {ILLUSTRATIVE_ANCHOR_FLUX_SI:.1f} kg/(m^2 s);"
    )
    print("         only the exponent is varied.  They are NOT fits to N2O/HTPB measurements.")
    print(RULE)
    print(
        f"{'law':<34} {'n':>7} {'a_SI':>13} {'r_dot':>10} {'m_dot_f':>11} {'O/F':>8}  "
        f"{'status':<12}"
    )
    print(f"{'':<34} {'[-]':>7} {'[SI]':>13} {'[mm/s]':>10} {'[kg/s]':>11} {'[-]':>8}  {'':<12}")
    print(RULE)
    for law in (SOURCED_LAW, *ILLUSTRATIVE_SENSITIVITY_LAWS):
        label = law.label.replace("$_2$", "2").replace("$n=", "n=").replace("$", "")
        status = "ILLUSTRATIVE" if law.is_illustrative else "sourced"
        print(
            f"{label:<34} {law.exponent:>7.4f} {law.coefficient_si:>13.6e} "
            f"{point.regression_rate_si(law) * 1e3:>10.4f} "
            f"{point.fuel_mass_flow_kg_s(law):>11.6f} "
            f"{point.mixture_ratio(law):>8.3f}  {status:<12}"
        )
    print(RULE)
    print()


def sanity_audit() -> list[tuple[str, bool]]:
    """Explicit engineering sanity checks over the whole study domain."""
    checks: list[tuple[str, bool]] = []

    flow_rows = study_oxidizer_flow_sweep()
    port_rows = study_port_diameter_sensitivity()
    all_rows = flow_rows + port_rows

    checks.append(("A_port > 0 everywhere", all(r["A_port_m2"] > 0 for r in port_rows)))
    checks.append(
        ("remaining fuel mass > 0 everywhere", all(r["fuel_left_kg"] > 0 for r in port_rows))
    )
    fluxes = [float(r["G_ox_si"]) for r in all_rows]
    checks.append(
        (
            "G_ox finite and in a reasonable small-hybrid band (1-500 kg/m^2 s)",
            all(math.isfinite(g) and 1.0 < g < 500.0 for g in fluxes),
        )
    )
    checks.append(
        ("r_dot > 0 for every positive G_ox", all(float(r["r_dot_mm_s"]) > 0 for r in all_rows))
    )
    rates = [float(r["r_dot_mm_s"]) for r in flow_rows]
    checks.append(
        (
            "r_dot increases monotonically with G_ox",
            all(b > a for a, b in itertools.pairwise(rates)),
        )
    )
    port_fluxes = [float(r["G_ox_si"]) for r in port_rows]
    checks.append(
        (
            "larger port at fixed m_dot_ox lowers G_ox",
            all(b < a for a, b in itertools.pairwise(port_fluxes)),
        )
    )
    port_rates = [float(r["r_dot_mm_s"]) for r in port_rows]
    checks.append(
        (
            "larger port at fixed m_dot_ox lowers r_dot",
            all(b < a for a, b in itertools.pairwise(port_rates)),
        )
    )
    checks.append(
        (
            "m_dot_f > 0 wherever regression is positive",
            all(r["m_dot_f_kg_s"] > 0 for r in all_rows),
        )
    )
    ofs = [float(r["of_ratio"]) for r in all_rows]
    checks.append(
        ("O/F finite and positive where defined", all(math.isfinite(o) and o > 0 for o in ofs))
    )
    numeric = [v for r in all_rows for v in r.values() if isinstance(v, float)]
    checks.append(
        ("no NaN/inf in the normal study domain", all(math.isfinite(v) for v in numeric))
    )
    return checks


def _print_sanity_audit() -> bool:
    print("ENGINEERING SANITY AUDIT")
    print(RULE)
    checks = sanity_audit()
    for description, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {description}")
    print(RULE)
    all_ok = all(ok for _, ok in checks)
    print(f"  {'ALL SANITY CHECKS PASSED' if all_ok else '*** SANITY CHECKS FAILED ***'}")
    print()
    return all_ok


def _print_observations(flow_rows, port_rows) -> None:
    rep = OperatingPoint(GRAIN, REPRESENTATIVE_M_DOT_OX_KG_S, GRAIN.initial_port_diameter_m)
    print("OBSERVATIONS (reported as computed, not tuned)")
    print(RULE)
    print(
        f"  Representative point: m_dot_ox = {REPRESENTATIVE_M_DOT_OX_KG_S:.3f} kg/s, "
        f"G_ox = {rep.oxidizer_mass_flux_si:.2f} kg/(m^2 s) "
        f"({rep.oxidizer_mass_flux_si / 10:.3f} g/(cm^2 s)),"
    )
    print(
        f"    r_dot = {rep.regression_rate_si(SOURCED_LAW) * 1e3:.4f} mm/s, "
        f"m_dot_f = {rep.fuel_mass_flow_kg_s(SOURCED_LAW):.6f} kg/s, "
        f"O/F = {rep.mixture_ratio(SOURCED_LAW):.3f}."
    )
    print(
        "  1. r_dot of order 0.7-1.0 mm/s over the sweep is consistent with the ~1 mm/s scale "
        "reported for HTPB/N2O."
    )
    print(
        "  2. The sourced exponent n = 0.3667 is BELOW the 0.5-0.8 band reported for most hybrid "
        "systems, so r_dot"
    )
    print(
        "     responds only weakly to flux here: a 3.3x rise in m_dot_ox lifts r_dot by just "
        f"{float(flow_rows[-1]['r_dot_mm_s']) / float(flow_rows[0]['r_dot_mm_s']):.2f}x."
    )
    print(
        "     Table 3 shows what a steeper exponent would imply. This discrepancy is reported, "
        "not tuned away."
    )
    print(
        f"  3. O/F over the sweep runs {min(float(r['of_ratio']) for r in flow_rows):.2f} to "
        f"{max(float(r['of_ratio']) for r in flow_rows):.2f}. This is below the 2.9-4.6 "
        "measured by the source,"
    )
    print(
        "     as expected: this illustrative grain is 400 mm long against the source's 250 mm, so "
        "it exposes ~1.6x"
    )
    print("     the burning area at a comparable port size and therefore makes more fuel.")
    print(
        f"  4. The largest sensitivity port ({PORT_DIAMETERS_M[-1] * 1e3:.0f} mm) drops G_ox to "
        f"{float(port_rows[-1]['G_ox_si']):.1f} kg/(m^2 s), below the source's"
    )
    print("     35 kg/(m^2 s) lower data bound; that row is flagged EXTRAPOLATED above.")
    print(RULE)
    print()


def main() -> int:
    _banner()
    _print_configuration()

    flow_rows = study_oxidizer_flow_sweep()
    _print_flow_sweep(flow_rows)

    port_rows = study_port_diameter_sensitivity()
    _print_port_sensitivity(port_rows, REPRESENTATIVE_M_DOT_OX_KG_S)

    _print_exponent_sensitivity(REPRESENTATIVE_M_DOT_OX_KG_S)
    _print_observations(flow_rows, port_rows)

    ok = _print_sanity_audit()
    print("END OF MILESTONE 1 STUDY - no chamber pressure, nozzle or thrust result is produced.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
