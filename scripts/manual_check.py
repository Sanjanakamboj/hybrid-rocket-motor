"""Hand-arithmetic cross-check of the Milestone 1 model chain.

Every quantity is recomputed here from raw ``math`` primitives, written out
longhand and in a different algebraic arrangement from the production code, and
then compared against the package result.  Nothing in this file imports the
production formulas; it only imports the geometry/law *constants* and the
objects whose output is being audited.

Run with::

    python scripts/manual_check.py
"""

from __future__ import annotations

import math

from hybrid_rocket_motor import (
    REPRESENTATIVE_GRAIN,
    REZAEI_2018_N2O_HTPB,
    OperatingPoint,
)

TOLERANCE = 1.0e-12
RULE = "-" * 86

# Prescribed inputs for the audited point.
M_DOT_OX_KG_S = 0.100
D_PORT_M = REPRESENTATIVE_GRAIN.initial_port_diameter_m
L_M = REPRESENTATIVE_GRAIN.length_m
D_OUTER_M = REPRESENTATIVE_GRAIN.outer_diameter_m
RHO_F = REPRESENTATIVE_GRAIN.fuel_density_kg_m3

# Coefficients exactly as printed in the source publication.
A_SOURCE_MM_S = 0.3977
N_SOURCE = 0.3667


def _compare(name: str, hand: float, package: float, unit: str) -> bool:
    rel = abs(hand - package) / abs(hand) if hand != 0.0 else abs(package)
    ok = rel <= TOLERANCE
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<34} hand = {hand:>16.9e} {unit:<12}")
    print(f"         {'':<34} pkg  = {package:>16.9e} {unit:<12} (rel. diff {rel:.2e})")
    return ok


def main() -> int:
    print("=" * 86)
    print("MILESTONE 1 MANUAL CHECK - longhand arithmetic vs package output")
    print("Generic illustrative N2O/HTPB grain. Prescribed oxidizer flow, Milestone 1 chain.")
    print("=" * 86)
    print()
    print("PRESCRIBED INPUTS")
    print(f"  m_dot_ox = {M_DOT_OX_KG_S} kg/s   D_p = {D_PORT_M} m   L = {L_M} m")
    print(f"  D_o      = {D_OUTER_M} m         rho_f = {RHO_F} kg/m^3")
    print()

    point = OperatingPoint(REPRESENTATIVE_GRAIN, M_DOT_OX_KG_S, D_PORT_M)
    results: list[bool] = []

    print("STEP 1  port area:  A_port = pi r^2,  r = D_p/2 = 0.020 m")
    radius = D_PORT_M / 2.0
    a_port = math.pi * radius * radius
    results.append(_compare("A_port", a_port, point.port_area_m2, "m^2"))
    print()

    print("STEP 2  burning perimeter:  P_port = 2 pi r")
    perimeter = 2.0 * math.pi * radius
    results.append(_compare("P_port", perimeter, point.port_perimeter_m, "m"))
    print()

    print("STEP 3  burning area:  A_burn = P_port * L")
    a_burn = perimeter * L_M
    results.append(_compare("A_burn", a_burn, point.burning_area_m2, "m^2"))
    print()

    print("STEP 4  remaining fuel mass:  m_f = rho_f * (pi/4)(D_o^2 - D_p^2) * L")
    fuel_mass = RHO_F * (math.pi / 4.0) * (D_OUTER_M**2 - D_PORT_M**2) * L_M
    results.append(_compare("m_f remaining", fuel_mass, point.remaining_fuel_mass_kg, "kg"))
    print()

    print("STEP 5  oxidizer mass flux:  G_ox = m_dot_ox / A_port = 4 m_dot_ox / (pi D_p^2)")
    flux = 4.0 * M_DOT_OX_KG_S / (math.pi * D_PORT_M**2)
    results.append(_compare("G_ox", flux, point.oxidizer_mass_flux_si, "kg/(m^2 s)"))
    print(f"         in source units: {flux / 10.0:.6f} g/(cm^2 s)")
    print()

    print("STEP 6  regression rate, evaluated in the SOURCE's units then converted:")
    print(f"         r_dot [mm/s] = {A_SOURCE_MM_S} * (G_ox [g/(cm^2 s)])^{N_SOURCE}")
    flux_g_cm2_s = flux / 10.0
    rate_mm_s = A_SOURCE_MM_S * flux_g_cm2_s**N_SOURCE
    rate_m_s = rate_mm_s / 1000.0
    print(f"         G_ox = {flux_g_cm2_s:.9f} g/(cm^2 s)  ->  r_dot = {rate_mm_s:.9f} mm/s")
    results.append(
        _compare("r_dot", rate_m_s, point.regression_rate_si(REZAEI_2018_N2O_HTPB), "m/s")
    )
    print()

    print("STEP 7  fuel mass flow:  m_dot_f = rho_f * A_burn * r_dot")
    m_dot_f = RHO_F * a_burn * rate_m_s
    results.append(
        _compare("m_dot_f", m_dot_f, point.fuel_mass_flow_kg_s(REZAEI_2018_N2O_HTPB), "kg/s")
    )
    print()

    print("STEP 8  mixture ratio:  O/F = m_dot_ox / m_dot_f")
    of_ratio = M_DOT_OX_KG_S / m_dot_f
    results.append(_compare("O/F", of_ratio, point.mixture_ratio(REZAEI_2018_N2O_HTPB), "-"))
    print()

    print("STEP 9  unit-conversion audit: collapsing the source coefficient into SI once")
    print("         a_SI = (1e-3 m/s per mm/s) * a_src * (10 kg/m^2 s per g/cm^2 s)^(-n)")
    a_si = 1.0e-3 * A_SOURCE_MM_S * math.pow(10.0, -N_SOURCE)
    results.append(_compare("a_SI", a_si, REZAEI_2018_N2O_HTPB.coefficient_si, "SI"))
    print()

    print("STEP 10 port-growth coupling at fixed m_dot_ox (sign check, not a burn history)")
    print(RULE)
    print(f"  {'D_p [mm]':>9} {'A_port [m^2]':>14} {'G_ox [kg/m^2 s]':>17} {'r_dot [mm/s]':>14}")
    print(RULE)
    previous_flux = math.inf
    previous_rate = math.inf
    monotonic = True
    for diameter in (0.040, 0.050, 0.060, 0.070):
        probe = OperatingPoint(REPRESENTATIVE_GRAIN, M_DOT_OX_KG_S, diameter)
        g = probe.oxidizer_mass_flux_si
        r = probe.regression_rate_si(REZAEI_2018_N2O_HTPB) * 1e3
        print(f"  {diameter * 1e3:>9.1f} {probe.port_area_m2:>14.6e} {g:>17.4f} {r:>14.4f}")
        monotonic = monotonic and g < previous_flux and r < previous_rate
        previous_flux, previous_rate = g, r
    print(RULE)
    print(f"  [{'OK  ' if monotonic else 'FAIL'}] larger port -> lower G_ox -> lower r_dot")
    results.append(monotonic)
    print()

    print("STEP 11 scope guard: confirm no DEFERRED physics is present in the model")
    # Scope-boundary note (Milestones 3, 4 and 5).
    # This guard originally forbade chamber-pressure, nozzle and thrust symbols
    # (correct through Milestone 2), then tank/injector/blowdown symbols (correct
    # through Milestone 3), then CEA and equilibrium-chemistry symbols (correct
    # through Milestone 4). Each milestone was chartered to add exactly what the
    # guard then forbade, so it is retargeted each time at the physics that is
    # STILL deferred. A documented scope change, not a correction: no Milestone 1
    # physics, coefficient or reported value is affected, and every numeric check
    # above is unchanged. See DESIGN.md.
    #
    # Equilibrium chemistry is now in scope. Equilibrium gives instantaneous,
    # complete reaction; it does not give reaction RATES, so finite-rate chemistry
    # and kinetics join the forbidden set rather than being treated as covered.
    forbidden = (
        "finite_rate",
        "kinetics",
        "instability",
        "ignition",
        "contour",
        "structural",
        "stress",
        "thermal_sizing",
        "flight_dynamics",
        "feedline",
        "feed_line",
        "heat_transfer",
    )
    import hybrid_rocket_motor as package_module

    exported = {name.lower() for name in dir(package_module)}
    leaked = sorted(name for name in forbidden if any(name in e for e in exported))
    clean = not leaked
    verdict = "OK  " if clean else "FAIL"
    print(
        f"  [{verdict}] no finite-rate / ignition / structural / thermal "
        "symbol exported"
    )
    if leaked:
        print(f"         leaked: {leaked}")
    results.append(clean)
    print()

    print("=" * 86)
    passed = all(results)
    print(f"MANUAL CHECK: {sum(results)}/{len(results)} checks passed"
          f"{'' if passed else '  *** FAILURES PRESENT ***'}")
    print("=" * 86)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
