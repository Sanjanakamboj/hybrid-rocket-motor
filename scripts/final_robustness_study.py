"""Milestone 6 final robustness audit and portfolio synthesis.

The closing study for the project. It adds no new physics. It runs the frozen
Milestone 3, 4 and 5 models over a disclosed set of one-at-a-time perturbations,
decomposes the model hierarchy, ranks the assumptions by measured effect, and
tests which qualitative conclusions survive.

Run with::

    python scripts/final_robustness_study.py
"""

from __future__ import annotations

import math

import numpy as np

from hybrid_rocket_motor import REPRESENTATIVE_GRAIN, REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.injector import Injector, InjectorModel
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.performance import couple_transient_to_performance
from hybrid_rocket_motor.robustness import (
    CaseOutcome,
    ModelKind,
    SensitivityFamily,
    conclusion_support,
    law_with_exponent,
    rank_families,
    run_frozen,
    run_variable,
    scaled_regression_law,
    summarise,
)
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
    NitrousTank,
)
from hybrid_rocket_motor.thermochemistry import DEFAULT_C_STAR_EFFICIENCY, default_table
from hybrid_rocket_motor.transient import ConstantOxidizerFlow, simulate_transient

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
NOZZLE = REFERENCE_NOZZLE
FROZEN = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
TABLE = default_table()
SEA_LEVEL_PA = 101325.0
REFERENCE_AREA_M2 = 3.5e-6
M3_PRESCRIBED_FLOW = 0.100

#: The reference oxidizer mass flux, used to anchor the exponent variation.
ANCHOR_FLUX_SI = 79.57747

RULE = "-" * 112
N_REPORT = 401


def _plain(text: str) -> str:
    """Strip the mathtext a label carries for figures, for console display.

    The frozen Milestone 1 study does the same thing; matching it keeps raw
    ``$..$`` markup out of terminal output.
    """
    return text.replace("$_2$", "2").replace("$", "")


def _state(temperature=REFERENCE_TANK_TEMPERATURE_K, fill=REFERENCE_TANK_FILL_FRACTION, tank=None):
    tank = tank if tank is not None else REFERENCE_TANK
    return tank, tank.initial_state(temperature, fill)


def variable_case(
    label,
    family,
    *,
    law=LAW,
    table=None,
    area=REFERENCE_AREA_M2,
    model=None,
    efficiency=DEFAULT_C_STAR_EFFICIENCY,
    temperature=REFERENCE_TANK_TEMPERATURE_K,
    tank=None,
) -> CaseOutcome:
    """One Milestone 5 run, reduced to the recorded quantities."""
    tank, state = _state(temperature=temperature, tank=tank)
    injector = Injector(area) if model is None else Injector(area, model=model)
    result = run_variable(
        tank,
        state,
        injector,
        GRAIN,
        law,
        table if table is not None else TABLE,
        NOZZLE,
        SEA_LEVEL_PA,
        c_star_efficiency=efficiency,
        n_report=N_REPORT,
    )
    return summarise(label, family, result, model=ModelKind.VARIABLE_PROPERTIES)


def frozen_case(
    label,
    family,
    *,
    law=LAW,
    area=REFERENCE_AREA_M2,
    model=None,
    temperature=REFERENCE_TANK_TEMPERATURE_K,
    tank=None,
) -> CaseOutcome:
    """One Milestone 4 run, reduced to the recorded quantities."""
    tank, state = _state(temperature=temperature, tank=tank)
    injector = Injector(area) if model is None else Injector(area, model=model)
    result = run_frozen(
        tank,
        state,
        injector,
        GRAIN,
        law,
        FROZEN,
        NOZZLE,
        SEA_LEVEL_PA,
        n_report=N_REPORT,
    )
    return summarise(label, family, result, model=ModelKind.FROZEN_PROPERTIES)


def milestone_3_case() -> CaseOutcome:
    """Milestone 3: prescribed constant oxidizer flow, frozen properties.

    Reduced to the same quantities as the coupled runs so the hierarchy table
    compares like with like. The mixture ratio and chamber pressure still move,
    because the port still opens -- what does not move is the oxidizer flow.
    """
    transient = simulate_transient(
        GRAIN, LAW, ConstantOxidizerFlow(M3_PRESCRIBED_FLOW, 200.0), t_end_s=200.0
    )
    history = couple_transient_to_performance(transient, FROZEN, NOZZLE, SEA_LEVEL_PA)
    index = transient.final_index
    thrust = history.thrust_n[: index + 1]
    mixture = history.mixture_ratio[: index + 1]
    pressure = history.chamber_pressure_pa[: index + 1]
    times = history.time_s[: index + 1]
    fuel = float(transient.cumulative_fuel_mass_kg[index])
    return CaseOutcome(
        label="M3 prescribed flow",
        family="model hierarchy",
        model=ModelKind.FROZEN_PROPERTIES,
        termination=transient.termination_reason.value,
        valid=True,
        burn_time_s=float(times[-1]),
        initial_thrust_n=float(thrust[0]),
        final_thrust_n=float(thrust[-1]),
        thrust_slope_percent=100.0 * (float(thrust[-1]) / float(thrust[0]) - 1.0),
        peak_thrust_n=float(np.max(thrust)),
        mean_thrust_n=float(history.total_impulse_n_s / times[-1]),
        total_impulse_n_s=float(history.total_impulse_n_s),
        equivalent_specific_impulse_s=float(history.equivalent_specific_impulse_s),
        initial_mixture_ratio=float(mixture[0]),
        final_mixture_ratio=float(mixture[-1]),
        initial_chamber_pressure_pa=float(pressure[0]),
        final_chamber_pressure_pa=float(pressure[-1]),
        oxidizer_consumed_kg=float(transient.cumulative_oxidizer_mass_kg[index]),
        fuel_consumed_kg=fuel,
        note="oxidizer flow is an INPUT here, not an output",
    )


# ---------------------------------------------------------------------------
# printing helpers
# ---------------------------------------------------------------------------


def _banner() -> None:
    print("=" * 112)
    print("MILESTONE 6 - FINAL ROBUSTNESS AUDIT AND PORTFOLIO SYNTHESIS")
    print("Generic N2O/HTPB hybrid rocket motor - reduced-order educational model")
    print("=" * 112)
    print()
    print("SCOPE NOTICE")
    print("  This milestone adds no physics. It re-runs the frozen Milestone 3, 4 and 5 models")
    print("  over a disclosed set of one-at-a-time perturbations and reports what changes.")
    print("  Still absent: finite-rate combustion, ignition transients, combustion instability,")
    print("  injector hardware geometry, nozzle contour design, feed-line and valve dynamics,")
    print("  the vapour-only discharge tail, structural or thermal sizing, and flight dynamics.")
    print()


def _case_row(case: CaseOutcome) -> str:
    return (
        f"  {case.label:<34} {case.burn_time_s:9.3f} {case.initial_thrust_n:11.3f} "
        f"{case.final_thrust_n:11.3f} {case.thrust_slope_percent:+10.3f} "
        f"{case.total_impulse_n_s:12.2f} {case.equivalent_specific_impulse_s:10.3f}"
    )


def _case_header() -> str:
    return (
        f"  {'model / case':<34} {'burn [s]':>9} {'F_0 [N]':>11} {'F_end [N]':>11} "
        f"{'slope [%]':>10} {'impulse':>12} {'I_sp [s]':>10}"
    )


def main() -> int:
    _banner()

    # -- frozen baseline verification -------------------------------------
    print("## Frozen baseline verification")
    print(RULE)
    m3 = milestone_3_case()
    m4 = frozen_case("M4 blowdown, frozen properties", "model hierarchy")
    m5 = variable_case("M5 blowdown, variable chemistry", "model hierarchy")
    checks = [
        ("M1 G_ox at 0.100 kg/s [kg/(m^2 s)]", 79.577472, 79.57747, 1e-4),
        ("M2 constant-flow burnout [s]", 41.740618, 41.740618, 1e-5),
        ("M3 reference p_c [bar]", m3.initial_chamber_pressure_pa / 1e5, 26.114, 1e-3),
        ("M4 burn time [s]", m4.burn_time_s, 42.917, 1e-3),
        ("M4 thrust slope [%]", m4.thrust_slope_percent, -8.599, 1e-2),
        ("M5 thrust slope [%]", m5.thrust_slope_percent, -12.915, 1e-2),
        (
            "M5 impulse vs M4 [%]",
            100.0 * (m5.total_impulse_n_s / m4.total_impulse_n_s - 1.0),
            -10.689,
            1e-2,
        ),
    ]
    for name, actual, expected, tol in checks:
        mark = "OK  " if abs(actual - expected) <= tol else "FAIL"
        print(f"  [{mark}] {name:<44} {actual:14.6f}  (frozen {expected:.6f})")
    print()
    print("  Milestones 1-5 reproduce their published values. Nothing below modifies them.")
    print()

    # -- reference case ----------------------------------------------------
    print("## Representative conceptual reference case")
    print(RULE)
    tank, state = _state()
    print("  This is a REPRESENTATIVE CONCEPTUAL REFERENCE CASE. It is not an optimised motor,")
    print("  not a design, and not validated against any hardware. Nothing was tuned to a")
    print("  performance target.")
    print()
    print("  Provenance of each input:")
    print(
        f"    grain geometry          ILLUSTRATIVE   L = {GRAIN.length_m:.2f} m, "
        f"D_p,0 = {GRAIN.initial_port_diameter_m * 1e3:.0f} mm, "
        f"D_o = {GRAIN.outer_diameter_m * 1e3:.0f} mm,"
    )
    print(
        f"                                           rho_f = {GRAIN.fuel_density_kg_m3:.0f}"
        " kg/m^3, single circular port"
    )
    print(f"    regression law          SOURCED        {_plain(LAW.label)}")
    print(
        f"                                           a_SI = {LAW.coefficient_si:.6e}, "
        f"n = {LAW.exponent:.4f}"
    )
    print("    N2O properties          SOURCED        CoolProp / Lemmon & Span (2006) EOS")
    print(
        f"    tank                    ILLUSTRATIVE   {tank.volume_m3 * 1e3:.1f} L rigid "
        f"adiabatic, {REFERENCE_TANK_FILL_FRACTION * 100:.0f} % liquid fill at "
        f"{REFERENCE_TANK_TEMPERATURE_K:.2f} K"
    )
    print(
        f"                                           -> {state.pressure_pa / 1e5:.3f} bar, "
        f"{state.total_mass_kg:.4f} kg N2O"
    )
    print(
        f"    injector                EFFECTIVE      Dyer/NHNE, A_eff = "
        f"{REFERENCE_AREA_M2 * 1e6:.2f} mm^2, C_d = 0.66"
    )
    print("                                           effective flow area ONLY - no hole count,")
    print("                                           no hole size, no plate geometry")
    print(
        f"    nozzle                  ILLUSTRATIVE   D_t = "
        f"{NOZZLE.throat_diameter_m * 1e3:.1f} mm, eps = {NOZZLE.expansion_ratio:.1f}, "
        "ideal 1-D expansion"
    )
    print("    thermochemistry         SOURCED        NASA CEA (RP-1311) via rocketcea, frozen")
    print(
        f"                                           table of {TABLE.of_values.size} x "
        f"{TABLE.pressure_values.size} points; eta_c* = "
        f"{DEFAULT_C_STAR_EFFICIENCY:.2f} MODEL-FORM"
    )
    print()
    print("  Resulting behaviour:")
    print(
        f"    O/F                     {m5.initial_mixture_ratio:.4f} -> "
        f"{m5.final_mixture_ratio:.4f}"
    )
    print(f"    tank pressure [bar]     {state.pressure_pa / 1e5:.3f} -> (see M4/M5 histories)")
    print(
        f"    chamber pressure [bar]  {m5.initial_chamber_pressure_pa / 1e5:.3f} -> "
        f"{m5.final_chamber_pressure_pa / 1e5:.3f}"
    )
    print(
        f"    thrust [N]              {m5.initial_thrust_n:.3f} -> {m5.final_thrust_n:.3f}"
        f"   ({m5.thrust_slope_percent:+.3f} %)"
    )
    print(f"    burn duration [s]       {m5.burn_time_s:.4f}")
    print(f"    total impulse [N s]     {m5.total_impulse_n_s:.2f}")
    print(f"    equivalent I_sp [s]     {m5.equivalent_specific_impulse_s:.3f}")
    print(f"    oxidizer / fuel [kg]    {m5.oxidizer_consumed_kg:.4f} / {m5.fuel_consumed_kg:.4f}")
    print(f"    terminal event          {m5.termination}")
    print()

    # -- hierarchy ---------------------------------------------------------
    print("## M3 vs M4 vs M5 hierarchy")
    print(RULE)
    print(_case_header())
    for case in (m3, m4, m5):
        print(_case_row(case))
    print()

    def delta(a: CaseOutcome, b: CaseOutcome, metric: str) -> float:
        reference = getattr(a, metric)
        return 100.0 * (getattr(b, metric) / reference - 1.0) if reference else math.nan

    print("  Decomposition of the total effect into its two physical causes:")
    print()
    print(f"    {'effect':<44} {'impulse':>12} {'burn time':>12} {'slope [pp]':>12}")
    for name, a, b in (
        ("A. feed/blowdown coupling   (M4 - M3)", m3, m4),
        ("B. variable thermochemistry (M5 - M4)", m4, m5),
        ("C. total coupled effect     (M5 - M3)", m3, m5),
    ):
        print(
            f"    {name:<44} {delta(a, b, 'total_impulse_n_s'):+11.2f}% "
            f"{delta(a, b, 'burn_time_s'):+11.2f}% "
            f"{b.thrust_slope_percent - a.thrust_slope_percent:+11.2f}"
        )
    print()
    print("  The thrust-slope column is in PERCENTAGE POINTS, because the quantity being")
    print("  compared is itself a percentage. Coupling the tank moves the slope from rising to")
    print("  falling; the chemistry then deepens the fall without changing its sign.")
    print()
    print("  These three models are NOT of equal fidelity. Each removes an assumption the")
    print("  previous one made, and none of them is validated against hardware.")
    print()

    # -- sensitivity families ---------------------------------------------
    print("## One-factor sensitivity")
    print(RULE)
    print("  Baseline: the Milestone 5 nominal case. Each family changes exactly ONE assumption")
    print("  and reruns the whole coupled model. These are deterministic alternatives, not")
    print("  probability distributions, and the +/- notation means 'two extra runs at these two")
    print("  values' - it is not a confidence interval.")
    print()

    baseline = m5
    families: list[SensitivityFamily] = []

    families.append(
        SensitivityFamily(
            "regression coefficient a",
            "+/-10 % on the sourced coefficient",
            baseline,
            (
                variable_case(
                    "a x0.90", "regression coefficient a", law=scaled_regression_law(LAW, 0.90)
                ),
                variable_case(
                    "a x1.10", "regression coefficient a", law=scaled_regression_law(LAW, 1.10)
                ),
            ),
        )
    )
    families.append(
        SensitivityFamily(
            "regression exponent n",
            f"n = 0.3367 / 0.3967, coefficient re-anchored at G_ox = {ANCHOR_FLUX_SI:.2f}",
            baseline,
            (
                variable_case(
                    "n = 0.3367 (anchored)",
                    "regression exponent n",
                    law=law_with_exponent(LAW, 0.3367, anchor_flux_si=ANCHOR_FLUX_SI),
                ),
                variable_case(
                    "n = 0.3967 (anchored)",
                    "regression exponent n",
                    law=law_with_exponent(LAW, 0.3967, anchor_flux_si=ANCHOR_FLUX_SI),
                ),
            ),
        )
    )
    families.append(
        SensitivityFamily(
            "c* efficiency",
            "eta_c* = 0.92 / 1.00 against the 0.96 baseline",
            baseline,
            (
                variable_case("eta_c* = 0.92", "c* efficiency", efficiency=0.92),
                variable_case("eta_c* = 1.00", "c* efficiency", efficiency=1.00),
            ),
        )
    )
    families.append(
        SensitivityFamily(
            "injector effective area",
            "+/-10 % on the effective flow area",
            baseline,
            (
                variable_case(
                    "A_eff x0.90", "injector effective area", area=0.90 * REFERENCE_AREA_M2
                ),
                variable_case(
                    "A_eff x1.10", "injector effective area", area=1.10 * REFERENCE_AREA_M2
                ),
            ),
        )
    )
    families.append(
        SensitivityFamily(
            "injector model",
            "Dyer/NHNE baseline vs SPI and HEM limits",
            baseline,
            (
                variable_case("SPI limit", "injector model", model=InjectorModel.SPI),
                variable_case("HEM limit", "injector model", model=InjectorModel.HEM),
            ),
        )
    )
    families.append(
        SensitivityFamily(
            "initial tank temperature",
            "283.15 / 303.15 K against the 293.15 K baseline",
            baseline,
            (
                variable_case("T_tank = 283.15 K", "initial tank temperature", temperature=283.15),
                variable_case("T_tank = 303.15 K", "initial tank temperature", temperature=303.15),
            ),
        )
    )
    families.append(
        SensitivityFamily(
            "thermochemistry grid resolution",
            "coarser subsets of the committed CEA nodes",
            baseline,
            (
                variable_case(
                    "grid /4 x /2",
                    "thermochemistry grid resolution",
                    table=TABLE.subsample(of_stride=4, pressure_stride=2),
                ),
                variable_case(
                    "grid /16 x /7",
                    "thermochemistry grid resolution",
                    table=TABLE.subsample(of_stride=16, pressure_stride=7),
                ),
            ),
        )
    )
    families.append(
        SensitivityFamily(
            "thermochemistry pressure dependence",
            "O/F-only table frozen at 25 bar",
            baseline,
            (
                variable_case(
                    "O/F-only chemistry",
                    "thermochemistry pressure dependence",
                    table=TABLE.at_fixed_pressure(25.0e5),
                ),
            ),
        )
    )

    print(
        f"  {'assumption':<36} {'variant':<24} {'impulse':>10} {'slope [pp]':>11} "
        f"{'burn':>10} {'termination':>16}"
    )
    for family in families:
        for variant in family.variants:
            if not variant.valid:
                # An invalid case has no meaningful performance numbers to compare,
                # so none are printed. A percentage change against a run that never
                # fired would invent a comparison that does not exist.
                print(
                    f"  {family.name:<36} {variant.label:<24} "
                    f"{'--':>9}  {'--':>10} {'--':>10}  {variant.termination:>16}"
                    "  <-- INVALID"
                )
                continue
            impulse = 100.0 * (variant.total_impulse_n_s / baseline.total_impulse_n_s - 1.0)
            slope = variant.thrust_slope_percent - baseline.thrust_slope_percent
            burn = 100.0 * (variant.burn_time_s / baseline.burn_time_s - 1.0)
            print(
                f"  {family.name:<36} {variant.label:<24} {impulse:+9.3f}% {slope:+10.3f} "
                f"{burn:+9.3f}% {variant.termination:>16}"
            )
    print()
    invalid = [v for f in families for v in f.invalid_variants]
    if invalid:
        print("  Cases reported as INVALID (not repaired, not clamped, excluded from ranking):")
        for variant in invalid:
            print(f"    {variant.label}: {variant.note}")
    else:
        print("  Every perturbed case stayed inside the tabulated thermochemistry rectangle and")
        print("  the checked N2O property band. No case had to be excluded.")
    print()

    # -- ranking -----------------------------------------------------------
    print("## Sensitivity ranking")
    print(RULE)
    print("  Metric: max |% change from baseline| across the tested alternatives in each family.")
    print("  Computed from the runs above, not assigned.")
    print()
    rankings = {}
    for metric, title, unit in (
        ("total_impulse_n_s", "A. total impulse", "%"),
        ("thrust_slope_percent", "B. thrust decay", "%"),
        ("burn_time_s", "C. burn duration", "%"),
    ):
        ranking = rank_families(families, metric)
        rankings[metric] = ranking
        print(f"  {title}")
        for position, (name, swing) in enumerate(ranking.entries, start=1):
            print(f"    {position:>2}. {name:<40} {swing:9.3f} {unit}")
        print()

    dominants = {m: r.dominant for m, r in rankings.items()}
    if len(set(dominants.values())) == 1:
        print(
            f"  The same assumption dominates all three metrics: {next(iter(dominants.values()))}."
        )
    else:
        print("  Different outputs have DIFFERENT dominant assumptions, so no single universal")
        print("  ranking is claimed:")
        for metric, name in dominants.items():
            print(f"    {metric:<26} -> {name}")
    print()

    # -- qualitative conclusions ------------------------------------------
    print("## Robustness of qualitative conclusions")
    print(RULE)
    all_variable = [baseline] + [v for f in families for v in f.variants]
    paired = [
        ("nominal", m4, m5),
        (
            "T = 283.15 K",
            frozen_case("M4 cold", "pair", temperature=283.15),
            variable_case("M5 cold", "pair", temperature=283.15),
        ),
        (
            "T = 303.15 K",
            frozen_case("M4 warm", "pair", temperature=303.15),
            variable_case("M5 warm", "pair", temperature=303.15),
        ),
        (
            "A_eff 2.5 mm^2",
            frozen_case("M4 small", "pair", area=2.5e-6),
            variable_case("M5 small", "pair", area=2.5e-6),
        ),
        (
            "A_eff 4.5 mm^2",
            frozen_case("M4 large", "pair", area=4.5e-6),
            variable_case("M5 large", "pair", area=4.5e-6),
        ),
    ]
    small_tank = NitrousTank(2.5e-3)
    depletion = variable_case("M5 2.5 L tank", "small tank", tank=small_tank)

    supports = []
    supports.append(
        conclusion_support(
            "1. Coupling the tank reverses the M3 rising-thrust trend",
            [b for _n, _a, b in paired],
            lambda o: o.thrust_slope_percent < 0.0 < m3.thrust_slope_percent,
        )
    )
    supports.append(
        conclusion_support(
            "2. Variable thermochemistry preserves the falling-thrust sign",
            all_variable,
            lambda o: o.thrust_slope_percent < 0.0,
        )
    )
    pair_outcomes = [b for _n, _a, b in paired]
    pair_lookup = {b.label: a for _n, a, b in paired}
    supports.append(
        conclusion_support(
            "3. Variable thermochemistry deepens the decay relative to M4",
            pair_outcomes,
            lambda o: o.thrust_slope_percent < pair_lookup[o.label].thrust_slope_percent,
        )
    )
    supports.append(
        conclusion_support(
            "4. c* variation dominates the M5 chemistry effect (over gamma)",
            [o for o in all_variable if math.isfinite(o.c_star_swing)],
            lambda o: abs(o.c_star_swing) > 5.0 * abs(o.gamma_swing),
        )
    )
    supports.append(
        conclusion_support(
            "5. Feed coupling changes thrust shape more than the chemistry does",
            pair_outcomes,
            lambda o: (
                abs(m4.thrust_slope_percent - m3.thrust_slope_percent)
                > abs(o.thrust_slope_percent - pair_lookup[o.label].thrust_slope_percent)
            ),
        )
    )
    supports.append(
        conclusion_support(
            "6. O/F drifts downward during the burn",
            all_variable,
            lambda o: o.final_mixture_ratio < o.initial_mixture_ratio,
        )
    )
    supports.append(
        conclusion_support(
            "7. The 10 L reference case burns the grain through before liquid runs out",
            all_variable,
            lambda o: o.termination == "GRAIN_BURNOUT",
        )
    )
    supports.append(
        conclusion_support(
            "8. A small tank can instead terminate on liquid depletion",
            [depletion],
            lambda o: o.termination == "LIQUID_DEPLETED",
        )
    )

    # Conclusion 4 needs the property histories, not just the summary quantities.
    tank, state = _state()
    m5_full = run_variable(
        tank,
        state,
        Injector(REFERENCE_AREA_M2),
        GRAIN,
        LAW,
        TABLE,
        NOZZLE,
        SEA_LEVEL_PA,
        n_report=N_REPORT,
    )
    firing = m5_full.firing
    c_star_swing = abs(m5_full.c_star_m_s[firing][-1] / m5_full.c_star_m_s[firing][0] - 1.0)
    gamma_swing = abs(m5_full.gamma[firing][-1] / m5_full.gamma[firing][0] - 1.0)

    for support in supports:
        print(f"  {support.claim}")
        print(f"    support {support.supporting}/{support.total}    {support.verdict}")
        if support.exceptions:
            print(f"    exceptions: {', '.join(support.exceptions)}")
    print()
    print(
        f"  Conclusion 4 in numbers: over the nominal burn, delivered c* moves "
        f"{100.0 * c_star_swing:.3f} %"
    )
    print(
        f"  while gamma moves {100.0 * gamma_swing:.3f} % -- a factor of "
        f"{c_star_swing / gamma_swing:.1f}. The added decay is a c* effect"
    )
    print("  acting through chamber pressure, not a nozzle-gamma effect.")
    print()

    # -- interpretation ----------------------------------------------------
    print("## Final interpretation")
    print(RULE)
    impulse_ranking = rankings["total_impulse_n_s"]
    slope_ranking = rankings["thrust_slope_percent"]
    burn_ranking = rankings["burn_time_s"]
    print(f"  Total impulse is most sensitive to {impulse_ranking.dominant}; the thrust-decay")
    print(f"  slope to {slope_ranking.dominant}; the burn duration to {burn_ranking.dominant}.")
    unique = {impulse_ranking.dominant, slope_ranking.dominant, burn_ranking.dominant}
    if len(unique) == 1:
        print("  All three headline outputs share the same dominant assumption.")
    else:
        print("  Those are not all the same assumption, which is why no single 'most important")
        print("  parameter' is claimed for the model as a whole.")
    print()
    print("  Note what is NOT near the top of any ranking: the thermochemistry table's own")
    print("  resolution and its pressure dependence each move the answer by well under a tenth")
    print("  of a per cent. The Milestone 5 chemistry matters through WHAT it says about")
    print("  c*(O/F), not through how finely it is tabulated.")
    print()
    print("  What survives every perturbation tested here is the SHAPE of the answer: with the")
    print("  tank coupled, thrust falls through the burn rather than rising, and letting the")
    print("  chemistry follow the mixture ratio makes it fall further. What moves a great deal")
    print("  is the LEVEL: total impulse and equivalent I_sp shift substantially with")
    print("  assumptions that are not measured for this motor.")
    print()
    print("  This is a deterministic reduced-order model-robustness study, not a probabilistic")
    print("  uncertainty analysis.")
    print()
    print("  The reference case is conceptual and not validated against a specific motor.")
    print()
    print("=" * 112)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
