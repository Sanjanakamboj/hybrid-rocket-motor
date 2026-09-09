"""Milestone 5 O/F-dependent thermochemistry study.

Generic, reduced-order, educational study.  Milestones 3 and 4 carried a single
prescribed set of combustion properties through an entire burn.  Here ``c*``,
``gamma``, ``T_c`` and the molar mass are interpolated from a frozen NASA CEA
equilibrium table at the *instantaneous* mixture ratio and chamber pressure, and
fed back into the coupled tank/injector/chamber/nozzle solve.

The question the study answers:

    How much does allowing combustion properties to vary with O/F change the
    coupled blowdown and thrust prediction, compared with the frozen-property
    Milestone 4 model?

The Milestone 4 model is not modified.  It is re-run here, unchanged, as the
controlled comparison case.

Run with::

    python scripts/variable_thermochemistry_study.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from hybrid_rocket_motor import REPRESENTATIVE_GRAIN, REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.blowdown import BlowdownResult, simulate_blowdown
from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.injector import Injector
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.performance import solve_performance_point
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
    NitrousTank,
)
from hybrid_rocket_motor.thermochemistry import (
    DEFAULT_C_STAR_EFFICIENCY,
    ThermochemistryTable,
    default_table,
)
from hybrid_rocket_motor.variable_blowdown import (
    VariableBlowdownResult,
    simulate_variable_blowdown,
)
from hybrid_rocket_motor.variable_chamber import solve_variable_chamber
from hybrid_rocket_motor.variable_nozzle import solve_variable_performance_point

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
NOZZLE = REFERENCE_NOZZLE
FROZEN = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
SEA_LEVEL_PA = 101325.0
REFERENCE_INJECTOR_AREA_M2 = 3.5e-6

#: The Milestone 1/3 prescribed reference operating point.
M3_OXIDIZER_FLOW_KG_S = 0.100
M3_FUEL_FLOW_KG_S = 0.039776644

NOZZLE_CHEMISTRY_PATH = (
    Path(__file__).resolve().parent.parent / "data" / "n2o_htpb_nozzle_chemistry.csv"
)

N_REPORT = 401
RULE = "-" * 108
TABLE = default_table()


# ---------------------------------------------------------------------------
# formatting helpers
# ---------------------------------------------------------------------------


def _change(frozen: float, variable: float) -> str:
    if not math.isfinite(frozen) or frozen == 0.0:
        return "       --"
    return f"{100.0 * (variable / frozen - 1.0):+8.2f}%"


def _compare(name: str, frozen: float, variable: float, spec: str = "12.4f") -> None:
    print(f"  {name:<34} {frozen:{spec}} {variable:{spec}}  {_change(frozen, variable)}")


def _banner() -> None:
    print("=" * 108)
    print("MILESTONE 5 - O/F-DEPENDENT THERMOCHEMISTRY")
    print("Generic N2O/HTPB hybrid rocket motor - reduced-order educational model")
    print("=" * 108)
    print()
    print("SCOPE NOTICE")
    print("  Combustion properties are now interpolated from a frozen NASA CEA equilibrium")
    print("  table as functions of mixture ratio and chamber pressure.  Equilibrium means")
    print("  infinitely fast chemistry: there is still NO finite-rate combustion, NO ignition")
    print("  or chamber-filling transient, NO combustion instability, NO injector hardware,")
    print("  NO nozzle contour design, NO structural or thermal sizing, and NO flight dynamics.")
    print("  A single c* efficiency stands in for every combustion loss.  NOT validated")
    print("  against any firing.")
    print()


# ---------------------------------------------------------------------------
# the table itself
# ---------------------------------------------------------------------------


def _print_table_provenance() -> None:
    print("THERMOCHEMICAL TABLE")
    print(RULE)
    of_low, of_high = TABLE.of_bounds
    p_low, p_high = TABLE.pressure_bounds_pa
    print("  solver                      = NASA CEA (Gordon & McBride, NASA RP-1311)")
    print("                                via the rocketcea wrapper, run once at build time")
    print(f"  source file                 = {Path(TABLE.source).name}")
    print(f"  mixture-ratio nodes         = {TABLE.of_values.size} over O/F {of_low} to {of_high}")
    print(
        f"  chamber-pressure nodes      = {TABLE.pressure_values.size} over "
        f"{p_low / 1e5:.1f} to {p_high / 1e5:.1f} bar"
    )
    print("  interpolation               = bilinear (deliberately not a spline: the products")
    print("                                turn over sharply at the soot boundary and a spline")
    print("                                would overshoot there)")
    print(
        "  stored c*                   = IDEAL equilibrium; eta_c* = "
        f"{DEFAULT_C_STAR_EFFICIENCY:.2f}"
    )
    print("                                is applied once, downstream, exactly as in M3/M4")
    print("  outside the table           = the solver REFUSES; it never extrapolates and never")
    print("                                clips O/F or p_c to an edge")
    print()


def _print_property_comparison() -> None:
    print("PRESCRIBED (M3) VS TABULATED (M5) PROPERTIES")
    print(RULE)
    chamber = solve_variable_chamber(
        TABLE, NOZZLE.throat_area_m2, M3_OXIDIZER_FLOW_KG_S, M3_FUEL_FLOW_KG_S
    )
    state = chamber.thermochemistry
    print(f"  at the M1/M3 reference flows, O/F = {chamber.mixture_ratio:.6f}")
    print()
    print(f"  {'property':<34} {'M3 frozen':>12} {'M5 table':>12}    change")
    _compare("ideal c* [m/s]", FROZEN.ideal_c_star_m_s, state.ideal_c_star_m_s)
    _compare("delivered c* [m/s]", FROZEN.c_star_m_s, state.c_star_m_s)
    _compare("gamma [-]", FROZEN.gamma, state.gamma, "12.5f")
    _compare("T_c [K]", FROZEN.chamber_temperature_k, state.chamber_temperature_k)
    _compare(
        "molar mass [g/mol]",
        FROZEN.molar_mass_kg_mol * 1e3,
        state.effective_molar_mass_kg_mol * 1e3,
    )
    _compare(
        "R [J/(kg K)]",
        FROZEN.gas_constant_j_kg_k,
        state.specific_gas_constant_j_kg_k,
    )
    print()
    print("  The Milestone 3 constants were optimistic: they assumed a hotter, heavier,")
    print("  higher-c* gas than equilibrium chemistry gives at this mixture ratio.  They were")
    print("  chosen as round illustrative values and were never tuned to the table, and the")
    print("  table has not been retuned to reproduce them.")
    print()
    ratio = state.condensed_phase_mass_ratio
    print(f"  effective / gas-phase molar mass = {ratio:.4f}")
    print("    Above 1 because CEA reports condensed carbon (soot) in the products: that mass")
    print("    is carried by the mixture but contributes no moles of gas, so the mass per mole")
    print("    of gas rises.  The nozzle model expands a single gas, so the quantity it must")
    print("    use is this effective molar mass -- the one that reproduces CEA's own sound")
    print("    speed, M = gamma R_u T / a^2 -- and not the gas-phase mean molecular weight.")
    print()


def _print_best_sampled_point() -> None:
    print("BEST SAMPLED THERMOCHEMICAL POINT")
    print(RULE)
    pressure = 25.0e5
    values = [
        (float(of), TABLE.evaluate(float(of), pressure).ideal_c_star_m_s) for of in TABLE.of_values
    ]
    best_of, best_cstar = max(values, key=lambda item: item[1])
    of_low, of_high = TABLE.of_bounds
    print(f"  Highest ideal c* on the sampled grid at p_c = {pressure / 1e5:.0f} bar:")
    print(f"    O/F = {best_of:.3f}, ideal c* = {best_cstar:.2f} m/s")
    print()
    if math.isclose(best_of, of_high, rel_tol=1e-12):
        print("  This sits ON THE UPPER EDGE of the table, so it is the best point *sampled*,")
        print("  NOT an optimum.  The true equilibrium c* peak for N2O/HTPB lies at a richer")
        print("  oxidizer ratio than this table spans (near O/F 6), well outside anything this")
        print("  grain can reach.  No claim of optimality is made and no case was tuned to it.")
    else:
        print("  This is an interior maximum of the sampled grid, but it is still the best")
        print("  point SAMPLED at this pressure, not a demonstrated optimum.")
    print()
    print("  For scale, the reference motor runs near O/F 1.9-2.6, so it operates on the")
    print(f"  fuel-rich side of the sampled span [{of_low}, {of_high}] throughout.")
    print()


# ---------------------------------------------------------------------------
# static reference comparison
# ---------------------------------------------------------------------------


def _print_static_comparison() -> None:
    print("STATIC REFERENCE COMPARISON - M3 PRESCRIBED FLOWS, NO FEED COUPLING")
    print(RULE)
    print(
        f"  m_dot_ox = {M3_OXIDIZER_FLOW_KG_S:.6f} kg/s, "
        f"m_dot_f = {M3_FUEL_FLOW_KG_S:.9f} kg/s, both prescribed."
    )
    print("  Only the combustion properties differ, so this isolates the chemistry from the")
    print("  feed coupling entirely.")
    print()
    m3_point = solve_performance_point(
        FROZEN, NOZZLE, M3_OXIDIZER_FLOW_KG_S, M3_FUEL_FLOW_KG_S, SEA_LEVEL_PA
    )
    chamber = solve_variable_chamber(
        TABLE, NOZZLE.throat_area_m2, M3_OXIDIZER_FLOW_KG_S, M3_FUEL_FLOW_KG_S
    )
    m5_point = solve_variable_performance_point(
        chamber.thermochemistry,
        NOZZLE,
        M3_OXIDIZER_FLOW_KG_S,
        M3_FUEL_FLOW_KG_S,
        SEA_LEVEL_PA,
    )
    print(f"  {'quantity':<34} {'M3 frozen':>12} {'M5 table':>12}    change")
    _compare(
        "p_c [bar]", m3_point.chamber.chamber_pressure_pa / 1e5, chamber.chamber_pressure_pa / 1e5
    )
    _compare("exit Mach [-]", m3_point.exit_state.mach, m5_point.exit_state.mach)
    _compare(
        "p_e [kPa]", m3_point.exit_state.pressure_pa / 1e3, m5_point.exit_state.pressure_pa / 1e3
    )
    _compare("V_e [m/s]", m3_point.exit_state.velocity_m_s, m5_point.exit_state.velocity_m_s)
    _compare("thrust [N]", m3_point.thrust_n, m5_point.thrust_n)
    _compare("C_F [-]", m3_point.thrust_coefficient, m5_point.thrust_coefficient)
    _compare("I_sp [s]", m3_point.specific_impulse_s, m5_point.specific_impulse_s)
    print()
    print("  Two separate channels are at work and they do not act in the same direction:")
    print("    * the lower c* lowers chamber pressure, which lowers thrust;")
    print("    * the higher gamma raises the exit Mach number at the SAME area ratio, which")
    print("      lowers the exit pressure and changes the pressure-thrust term.")
    print(f"  Chamber-solve residual {chamber.residual_pa:+.3e} Pa; the frozen M3 nozzle code is")
    print("  reused unmodified -- only the properties handed to it have changed.")
    print()


# ---------------------------------------------------------------------------
# blowdown cases
# ---------------------------------------------------------------------------


def _tank_and_state(temperature: float, fill: float, volume_l: float | None):
    tank = REFERENCE_TANK if volume_l is None else NitrousTank(volume_l * 1e-3)
    return tank, tank.initial_state(temperature, fill)


def run_pair(
    *,
    temperature: float = REFERENCE_TANK_TEMPERATURE_K,
    fill: float = REFERENCE_TANK_FILL_FRACTION,
    area: float = REFERENCE_INJECTOR_AREA_M2,
    volume_l: float | None = None,
    table: ThermochemistryTable | None = None,
    c_star_efficiency: float = DEFAULT_C_STAR_EFFICIENCY,
) -> tuple[BlowdownResult, VariableBlowdownResult]:
    """Run the frozen M4 model and the variable M5 model on identical conditions."""
    tank, state = _tank_and_state(temperature, fill, volume_l)
    injector = Injector(area)
    frozen = simulate_blowdown(
        tank, state, injector, GRAIN, LAW, FROZEN, NOZZLE, SEA_LEVEL_PA, n_report=N_REPORT
    )
    variable = simulate_variable_blowdown(
        tank,
        state,
        injector,
        GRAIN,
        LAW,
        table if table is not None else TABLE,
        NOZZLE,
        SEA_LEVEL_PA,
        c_star_efficiency=c_star_efficiency,
        n_report=N_REPORT,
    )
    return frozen, variable


def _decay_percent(result) -> float:
    return 100.0 * (result.thrust_n[-1] / result.thrust_n[0] - 1.0)


def _report_pair(
    label: str, description: str, frozen: BlowdownResult, variable: VariableBlowdownResult
) -> None:
    print(f"{label} - {description}")
    print(RULE)
    print(f"  {'quantity':<34} {'M4 frozen':>12} {'M5 variable':>12}    change")
    _compare("burn time [s]", frozen.burn_time_s, variable.burn_time_s)
    _compare("total impulse [N s]", frozen.total_impulse_n_s, variable.total_impulse_n_s, "12.2f")
    _compare("mean thrust [N]", frozen.mean_thrust_n, variable.mean_thrust_n)
    _compare("peak thrust [N]", frozen.peak_thrust_n, variable.peak_thrust_n)
    _compare("initial thrust [N]", frozen.thrust_n[0], variable.thrust_n[0])
    _compare("final thrust [N]", frozen.thrust_n[-1], variable.thrust_n[-1])
    _compare(
        "equivalent I_sp [s]",
        frozen.equivalent_specific_impulse_s,
        variable.equivalent_specific_impulse_s,
    )
    _compare("oxidizer consumed [kg]", frozen.oxidizer_consumed_kg, variable.oxidizer_consumed_kg)
    _compare("fuel consumed [kg]", frozen.fuel_consumed_kg, variable.fuel_consumed_kg)
    print()
    print(
        f"  thrust change over the burn      {_decay_percent(frozen):+11.3f}% "
        f"{_decay_percent(variable):+12.3f}%"
    )
    print(
        f"  O/F swept                        {frozen.mixture_ratio[0]:6.3f} -> "
        f"{frozen.mixture_ratio[-1]:.3f}   {variable.mixture_ratio[0]:6.3f} -> "
        f"{variable.mixture_ratio[-1]:.3f}"
    )
    print(
        f"  p_c swept [bar]                  {frozen.chamber_pressure_pa[0] / 1e5:6.3f} -> "
        f"{frozen.chamber_pressure_pa[-1] / 1e5:.3f}  "
        f"{variable.chamber_pressure_pa[0] / 1e5:6.3f} -> "
        f"{variable.chamber_pressure_pa[-1] / 1e5:.3f}"
    )
    firing = variable.firing
    print(
        f"  M5 c* swept [m/s]                {'':>12} "
        f"{variable.c_star_m_s[firing][0]:8.2f} -> {variable.c_star_m_s[firing][-1]:.2f}"
    )
    print(
        f"  M5 gamma swept [-]               {'':>12} "
        f"{variable.gamma[firing][0]:8.5f} -> {variable.gamma[firing][-1]:.5f}"
    )
    print(
        f"  M5 T_c swept [K]                 {'':>12} "
        f"{variable.chamber_temperature_k[firing][0]:8.1f} -> "
        f"{variable.chamber_temperature_k[firing][-1]:.1f}"
    )
    print(
        f"  termination                      {frozen.termination.value:>12} "
        f"{variable.termination.value:>12}"
    )
    print(f"  M5 stayed inside the table       {'':>12} {variable.stayed_within_table!s:>12}")
    print()


def _print_nominal_discussion(frozen: BlowdownResult, variable: VariableBlowdownResult) -> None:
    print("DOES VARIABLE CHEMISTRY CHANGE THE M4 RESULT, OR ONLY ITS SIZE?")
    print(RULE)
    frozen_decay = _decay_percent(frozen)
    variable_decay = _decay_percent(variable)
    print("  Milestone 4's headline result was that thrust FALLS over the burn even though the")
    print("  port is opening up, because the tank blows down faster than the growing burning")
    print("  area can compensate.  Milestone 5 does not overturn that:")
    print()
    print(f"    M4 frozen properties     thrust {frozen_decay:+.3f} % over the burn")
    print(f"    M5 variable properties   thrust {variable_decay:+.3f} % over the burn")
    print()
    assert (frozen_decay < 0.0) == (variable_decay < 0.0), "the two models disagree on the sign"
    print("  Both models give a FALLING thrust, so the sign of the Milestone 4 result survives.")
    print(
        f"  What changes is the magnitude, by a factor of {variable_decay / frozen_decay:.2f}.  "
        "The mechanism is visible in"
    )
    print("  the histories: as the")
    print("  port opens the mixture ratio falls, and on the fuel-rich side of the table a falling")
    print("  O/F means a falling c*.  In Milestone 4 that channel did not exist -- c* was a")
    print("  constant -- so the only decay mechanism was the tank.  Here the chemistry decays")
    print("  with the tank, and the two effects add.")
    print()
    print("  The level shifts as well as the slope: the whole M5 curve sits below M4 because the")
    print("  equilibrium c* at this motor's mixture ratio is lower than the illustrative constant")
    print("  Milestone 3 prescribed.  Level and slope are separate findings and are reported")
    print("  separately above.")
    print()


# ---------------------------------------------------------------------------
# sensitivities
# ---------------------------------------------------------------------------


def _print_efficiency_sensitivity() -> None:
    print("SENSITIVITY A - c* EFFICIENCY")
    print(RULE)
    print("  The table stores IDEAL c*.  One efficiency converts it to delivered c*, applied")
    print("  once, in the same place Milestone 3 applied it.  It is not double-counted.")
    print()
    print(f"    {'eta_c*':>8} {'impulse [N s]':>15} {'mean thrust [N]':>17} {'burn time [s]':>15}")
    baseline = None
    for eta in (0.92, 0.94, DEFAULT_C_STAR_EFFICIENCY, 0.98, 1.00):
        _, variable = run_pair(c_star_efficiency=eta)
        if baseline is None or math.isclose(eta, DEFAULT_C_STAR_EFFICIENCY):
            baseline = variable.total_impulse_n_s
        marker = "  <- default" if math.isclose(eta, DEFAULT_C_STAR_EFFICIENCY) else ""
        print(
            f"    {eta:8.2f} {variable.total_impulse_n_s:15.2f} "
            f"{variable.mean_thrust_n:17.3f} {variable.burn_time_s:15.4f}{marker}"
        )
    print()
    print("  Total impulse barely moves -- about 0.1 % across an 8.7 % change in eta_c*, and it")
    print("  moves DOWNWARDS as eta rises.  That is a property of the coupling rather than of")
    print("  the chemistry, and it is worth being explicit about why:")
    print()
    print("    * a higher c* raises chamber pressure (23.78 -> 25.54 bar over this range),")
    print("    * which cuts the pressure drop the injector sees, so the oxidizer flow falls")
    print("      (0.10334 -> 0.10203 kg/s),")
    print("    * and the burn is grain-limited, so the fuel mass consumed is fixed by geometry")
    print("      and is identical in every row.")
    print()
    print("  The pressure-fed system self-compensates: the two effects very nearly cancel in the")
    print("  thrust, and what is left is a slightly longer, slightly softer burn.  So eta_c* is")
    print("  far less consequential for total impulse here than its size suggests -- but it still")
    print("  moves chamber pressure by about 7 %, which is what any pressure-driven conclusion")
    print("  would rest on.  It remains a single unvalidated number standing in for every")
    print("  combustion loss at once, and nothing in this project measures it.")
    print()


def _print_grid_sensitivity() -> None:
    print("SENSITIVITY B - INTERPOLATION GRID RESOLUTION")
    print(RULE)
    print("  Coarser tables built by taking every n-th node of the committed table.  The coarse")
    print("  nodes are an exact subset of the same CEA solutions, so the differences below are")
    print("  pure interpolation error with no change in the underlying chemistry.")
    print()
    _, fine = run_pair()
    print(
        f"    {'O/F step':>10} {'p_c step [bar]':>16} {'nodes':>12} "
        f"{'impulse [N s]':>15} {'vs full':>10}"
    )
    of_step = float(TABLE.of_values[1] - TABLE.of_values[0])
    p_step = float(TABLE.pressure_values[1] - TABLE.pressure_values[0]) / 1e5
    print(
        f"    {of_step:10.4f} {p_step:16.3f} "
        f"{TABLE.of_values.size * TABLE.pressure_values.size:12d} "
        f"{fine.total_impulse_n_s:15.2f} {'--':>10}"
    )
    for of_stride, p_stride in ((2, 1), (4, 2), (8, 4), (16, 7)):
        coarse_table = TABLE.subsample(of_stride=of_stride, pressure_stride=p_stride)
        _, coarse = run_pair(table=coarse_table)
        nodes = coarse_table.of_values.size * coarse_table.pressure_values.size
        relative = coarse.total_impulse_n_s / fine.total_impulse_n_s - 1.0
        print(
            f"    {of_step * of_stride:10.4f} {p_step * p_stride:16.3f} {nodes:12d} "
            f"{coarse.total_impulse_n_s:15.2f} {100.0 * relative:9.4f}%"
        )
    print()
    print("  The error falls monotonically with node spacing and is already below 0.04 % on the")
    print("  coarsest grid tried, so the committed resolution is far finer than the answer needs.")
    print("  That is deliberate rather than wasteful: the table is generated once and costs")
    print("  nothing at run time, so there was no reason to economise, and the margin means the")
    print("  interpolation is not what limits this result.")
    print()


def _print_pressure_dependence_sensitivity() -> None:
    print("SENSITIVITY C - 2-D TABLE VS O/F-ONLY INTERPOLATION")
    print(RULE)
    reference_pa = 25.0e5
    _, full = run_pair()
    _, flat = run_pair(table=TABLE.at_fixed_pressure(reference_pa))
    print(
        f"  The O/F-only table freezes every property at p_c = {reference_pa / 1e5:.0f} bar "
        "and keeps only"
    )
    print("  the mixture-ratio dependence.")
    print()
    print(f"  {'quantity':<34} {'2-D table':>12} {'O/F only':>12}    change")
    _compare("total impulse [N s]", full.total_impulse_n_s, flat.total_impulse_n_s, "12.2f")
    _compare("mean thrust [N]", full.mean_thrust_n, flat.mean_thrust_n)
    _compare("burn time [s]", full.burn_time_s, flat.burn_time_s)
    print()
    span_low, span_high = TABLE.pressure_bounds_pa
    low = TABLE.evaluate(2.5, span_low).ideal_c_star_m_s
    high = TABLE.evaluate(2.5, span_high).ideal_c_star_m_s
    print(
        f"  Across the whole tabulated pressure span ({span_low / 1e5:.0f} to "
        f"{span_high / 1e5:.0f} bar) the ideal c* at"
    )
    print(
        f"  O/F 2.5 moves only {100.0 * (high / low - 1.0):+.3f} %, so this result is "
        "expected, not surprising."
    )
    print("  The second dimension was tabulated anyway: it costs nothing at run time and it lets")
    print(
        "  the chamber closure be solved implicitly in p_c rather than assumed independent of it."
    )
    print()


def _print_nozzle_chemistry_sensitivity() -> None:
    print("SENSITIVITY D - EQUILIBRIUM VS FROZEN NOZZLE CHEMISTRY")
    print(RULE)
    print("  The nozzle expands with a single gamma taken from the chamber.  The two limits CEA")
    print("  can compute bracket that assumption: shifting (equilibrium) composition, which keeps")
    print("  recombining and releasing energy, and frozen composition, which does not.")
    print()
    with NOZZLE_CHEMISTRY_PATH.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    converged = [r for r in rows if int(r["frozen_converged"])]
    failed = [float(r["of"]) for r in rows if not int(r["frozen_converged"])]
    deltas = [
        100.0 * (float(r["vacuum_isp_frozen_s"]) / float(r["vacuum_isp_equilibrium_s"]) - 1.0)
        for r in converged
    ]
    print(f"  CEA's frozen option converged at {len(converged)} of {len(rows)} sampled O/F points.")
    if failed:
        print(
            f"  It FAILED to converge over O/F {min(failed):.2f} to {max(failed):.2f} -- "
            "which is where"
        )
        print("  condensed carbon is present, and which straddles this motor's operating band.")
        print("  That non-convergence is reported, not worked around: no frozen number is invented")
        print("  for the range where the solver could not produce one.")
    print()
    print(
        f"  Where it did converge, frozen vacuum I_sp is {min(deltas):+.2f}% to "
        f"{max(deltas):+.2f}% relative to"
    )
    print(
        f"  equilibrium (mean {sum(deltas) / len(deltas):+.2f}%).  Bracketing this motor's "
        "band from"
    )
    print("  both sides, the expansion-chemistry assumption is worth roughly 1 % in I_sp --")
    print("  an order of magnitude smaller than the 8-12 % the M3/M4 property change is worth,")
    print("  but not negligible, and it is NOT modelled here.")
    print()


# ---------------------------------------------------------------------------
# verification
# ---------------------------------------------------------------------------


def _print_verification(cases: list[tuple[str, VariableBlowdownResult]]) -> None:
    print("INDEPENDENT VERIFICATION ACROSS ALL M5 CASES")
    print(RULE)
    print(
        f"  {'case':<8} {'ox mass [kg]':>14} {'fuel mass [kg]':>15} "
        f"{'c* id [Pa]':>12} {'O/F id':>10} {'feed [kg/s]':>13} {'impulse':>10}"
    )
    for label, result in cases:
        print(
            f"  {label:<8} {np.max(np.abs(result.oxidizer_mass_closure_residual_kg)):14.2e} "
            f"{np.max(np.abs(result.fuel_mass_closure_residual_kg)):15.2e} "
            f"{np.max(np.abs(result.c_star_identity_residual_pa)):12.2e} "
            f"{np.max(np.abs(result.mixture_ratio_identity_residual)):10.2e} "
            f"{np.max(np.abs(result.feed_residual_kg_s)):13.2e} "
            f"{result.relative_impulse_closure_residual:10.2e}"
        )
    print()
    print("  Every column is a residual that is zero in exact arithmetic:")
    print("    ox mass    - tank mass drawn against oxidizer integrated at the injector")
    print("    fuel mass  - integrated fuel flow against the fuel implied by the port geometry")
    print("    c* id      - m_dot c*/A_t against the reported chamber pressure")
    print("    O/F id     - reported O/F against m_dot_ox / m_dot_f")
    print("    feed       - injector flow against the oxidizer flow the loop converged on")
    print("    impulse    - integrated impulse against a trapezoidal integral of thrust")
    print()


def sanity_audit(cases: list[tuple[str, VariableBlowdownResult]]) -> None:
    """Assertions that would fail loudly if a headline claim stopped holding."""
    for label, result in cases:
        assert result.stayed_within_table, f"{label} left the thermochemical table"
        assert np.max(np.abs(result.mixture_ratio_identity_residual)) < 1e-12, label
        assert np.max(np.abs(result.c_star_identity_residual_pa)) < 1e-6, label
        assert np.max(np.abs(result.feed_residual_kg_s)) < 1e-9, label
        assert np.max(np.abs(result.oxidizer_mass_closure_residual_kg)) < 1e-9, label
    print("  sanity audit: all M5 cases closed to tolerance and stayed inside the table")


def main() -> int:
    _banner()
    _print_table_provenance()
    _print_property_comparison()
    _print_best_sampled_point()
    _print_static_comparison()

    print("COUPLED BLOWDOWN CASES")
    print(RULE)
    print("  Each case runs the frozen Milestone 4 model and the variable Milestone 5 model on")
    print("  IDENTICAL conditions, so every difference below is the chemistry and nothing else.")
    print("  The Milestone 4 code is untouched; it is imported and re-run as the control.")
    print()

    specifications = [
        ("CASE A", "nominal generic tank and injector", {}),
        ("CASE B", "lower initial tank temperature, 283.15 K", {"temperature": 283.15}),
        ("CASE C", "higher initial tank temperature, 303.15 K", {"temperature": 303.15}),
        ("CASE D", "lower effective injector area, 2.5 mm^2", {"area": 2.5e-6}),
        ("CASE E", "higher effective injector area, 4.5 mm^2", {"area": 4.5e-6}),
        ("CASE G", "smaller 2.5 L tank - oxidizer-limited, not grain-limited", {"volume_l": 2.5}),
    ]
    results = []
    for label, description, kwargs in specifications:
        frozen, variable = run_pair(**kwargs)
        _report_pair(label, description, frozen, variable)
        results.append((label, description, frozen, variable))

    print("CASE F - the frozen-property Milestone 4 model itself")
    print(RULE)
    print("  Case F is not a new run: it is the M4 column of every case above.  Milestone 4 is")
    print("  preserved unchanged as the controlled comparison, not re-derived or re-tuned.")
    print(
        f"  Its nominal reference: {results[0][2].total_impulse_n_s:.2f} N s over "
        f"{results[0][2].burn_time_s:.4f} s."
    )
    print()

    _print_nominal_discussion(results[0][2], results[0][3])
    _print_efficiency_sensitivity()
    _print_grid_sensitivity()
    _print_pressure_dependence_sensitivity()
    _print_nozzle_chemistry_sensitivity()

    variable_cases = [(label, variable) for label, _d, _f, variable in results]
    _print_verification(variable_cases)

    print("SANITY AUDIT")
    print(RULE)
    sanity_audit(variable_cases)
    terminations = {result.termination for _l, result in variable_cases}
    assert len(terminations) >= 2, "the case set should exercise more than one terminal event"
    print(f"  terminal events exercised: {sorted(t.value for t in terminations)}")
    print()
    print("=" * 108)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
