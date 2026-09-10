"""Deterministic generation of the Milestone 6 final portfolio figures.

Produces five PNGs in ``figures/``:

* ``final_model_hierarchy.png``          - M3 / M4 / M5 flow, pressure and thrust
* ``final_thermochemistry_evolution.png`` - which property actually drives the M5 decay
* ``final_sensitivity_ranking.png``      - deterministic one-factor tornado
* ``final_reference_case.png``           - the representative conceptual reference case
* ``final_coupling_decomposition.png``   - feed effect vs chemistry effect, signed

Regeneration is byte-identical for a given matplotlib/FreeType build. PNG glyph
rasterisation depends on the FreeType version, so hashes are not expected to match
across environments with different builds; the content is unchanged. Verified with
matplotlib 3.10.9 / FreeType 2.6.1.

Run with::

    python scripts/generate_final_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from hybrid_rocket_motor import REPRESENTATIVE_GRAIN, REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.blowdown import simulate_blowdown
from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.injector import Injector, InjectorModel
from hybrid_rocket_motor.nozzle import REFERENCE_NOZZLE
from hybrid_rocket_motor.performance import couple_transient_to_performance
from hybrid_rocket_motor.robustness import (
    ModelKind,
    SensitivityFamily,
    law_with_exponent,
    rank_families,
    scaled_regression_law,
    summarise,
)
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
)
from hybrid_rocket_motor.thermochemistry import default_table
from hybrid_rocket_motor.transient import ConstantOxidizerFlow, simulate_transient
from hybrid_rocket_motor.variable_blowdown import simulate_variable_blowdown

sys.path.insert(0, str(Path(__file__).resolve().parent))

from final_robustness_study import ANCHOR_FLUX_SI, variable_case

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
NOZZLE = REFERENCE_NOZZLE
FROZEN = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
TABLE = default_table()
SEA_LEVEL_PA = 101325.0
AREA_M2 = 3.5e-6
M3_FLOW = 0.100
N_REPORT = 401

M3_COLOUR = "#7570b3"
M4_COLOUR = "#d95f02"
M5_COLOUR = "#0b5394"
ACCENT = "#a11b1b"
SECONDARY = "#1b9e77"
BAND = "#bcd7ee"

DISCLAIMER = (
    "Generic coupled N₂O/HTPB reduced-order model — NOT experimentally validated. Conceptual "
    "reference case, not an optimised motor.\n"
    "Injector represented as an effective flow area only. No ignition transient, no finite-rate "
    "chemistry, no structural or thermal design."
)


def _configure_matplotlib() -> None:
    plt.rcdefaults()
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 150,
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linewidth": 0.6,
            "legend.fontsize": 8.5,
            "legend.framealpha": 0.95,
            "lines.linewidth": 1.8,
            "mathtext.fontset": "dejavusans",
            "svg.hashsalt": "hybrid-rocket-motor-m6",
            "path.simplify": False,
            "axes.axisbelow": True,
        }
    )


def _readable(termination: str) -> str:
    """Turn an enum name into prose, so no raw code-style label reaches a figure."""
    return termination.replace("_", " ").lower()


def _save(fig, name: str) -> Path:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / name
    fig.savefig(path, format="png", metadata={"Software": None})
    plt.close(fig)
    return path


def _load_models():
    """The three model levels on the same nominal configuration."""
    tank, state = (
        REFERENCE_TANK,
        REFERENCE_TANK.initial_state(REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION),
    )
    injector = Injector(AREA_M2)
    transient = simulate_transient(GRAIN, LAW, ConstantOxidizerFlow(M3_FLOW, 200.0), t_end_s=200.0)
    m3 = couple_transient_to_performance(transient, FROZEN, NOZZLE, SEA_LEVEL_PA)
    m3_index = transient.final_index
    m4 = simulate_blowdown(
        tank, state, injector, GRAIN, LAW, FROZEN, NOZZLE, SEA_LEVEL_PA, n_report=N_REPORT
    )
    m5 = simulate_variable_blowdown(
        tank, state, injector, GRAIN, LAW, TABLE, NOZZLE, SEA_LEVEL_PA, n_report=N_REPORT
    )
    return transient, m3, m3_index, m4, m5, state


# ---------------------------------------------------------------------------
# Figure 1 -- model hierarchy
# ---------------------------------------------------------------------------


def figure_model_hierarchy(transient, m3, m3_index, m4, m5) -> Path:
    fig, axes = plt.subplots(3, 1, figsize=(9.2, 10.2), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.965, top=0.925, bottom=0.165, hspace=0.16)
    fig.suptitle(
        "Final Figure 1 — model hierarchy: what each removed assumption changes", fontsize=12.5
    )
    ax_flow, ax_pressure, ax_thrust = axes

    m3_time = m3.time_s[: m3_index + 1]
    series = (
        ("M3  prescribed flow, frozen properties", M3_COLOUR, ":", m3_time),
        ("M4  tank/feed coupled, frozen properties", M4_COLOUR, "--", m4.time_s),
        ("M5  tank/feed coupled, CEA chemistry", M5_COLOUR, "-", m5.time_s),
    )

    flows = (
        transient.oxidizer_mass_flow_kg_s[: m3_index + 1],
        m4.oxidizer_mass_flow_kg_s,
        m5.oxidizer_mass_flow_kg_s,
    )
    pressures = (
        m3.chamber_pressure_pa[: m3_index + 1] / 1e5,
        m4.chamber_pressure_pa / 1e5,
        m5.chamber_pressure_pa / 1e5,
    )
    thrusts = (m3.thrust_n[: m3_index + 1], m4.thrust_n, m5.thrust_n)

    for axis, values, label in (
        (ax_flow, flows, r"oxidizer flow  $\dot{m}_{ox}$   [kg/s]"),
        (ax_pressure, pressures, r"chamber pressure  $p_c$   [bar]"),
        (ax_thrust, thrusts, "thrust  $F$   [N]"),
    ):
        for (name, colour, style, time), value in zip(series, values, strict=True):
            axis.plot(time, value, color=colour, linestyle=style, zorder=4, label=name)
            axis.plot([time[-1]], [value[-1]], marker="o", color=colour, markersize=5, zorder=5)
        axis.set_ylabel(label)

    ax_flow.set_ylim(0.0, 0.135)
    ax_flow.legend(loc="lower left")
    ax_flow.annotate(
        "M3 holds this constant BY ASSUMPTION;\nin M4/M5 it is an output of the tank",
        xy=(0.985, 0.94),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=8.4,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": "0.7"},
    )
    ax_thrust.set_xlabel("time  $t$   [s]")
    ax_thrust.set_xlim(0.0, max(m3_time[-1], m4.time_s[-1], m5.time_s[-1]) * 1.02)
    ax_thrust.set_ylim(0.0, 380.0)
    ax_thrust.annotate(
        f"M3 {100 * (thrusts[0][-1] / thrusts[0][0] - 1):+.1f} %   rising\n"
        f"M4 {100 * (thrusts[1][-1] / thrusts[1][0] - 1):+.1f} %   falling\n"
        f"M5 {100 * (thrusts[2][-1] / thrusts[2][0] - 1):+.1f} %   falling further",
        xy=(0.985, 0.30),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=8.6,
        bbox={"boxstyle": "round,pad=0.42", "facecolor": "white", "edgecolor": "0.7"},
    )
    fig.text(
        0.105,
        0.115,
        "Filled markers are the terminal event of each run. The three models are NOT of equal "
        "fidelity: each removes an assumption\n"
        "the previous one made, and none is validated against hardware. Coupling the tank "
        "reverses the sign of the thrust-time slope;\n"
        "letting the chemistry follow the mixture ratio then deepens the fall without changing "
        "its sign.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.105, 0.040, DISCLAIMER, fontsize=7.5, color="0.42", va="top")
    return _save(fig, "final_model_hierarchy.png")


# ---------------------------------------------------------------------------
# Figure 2 -- thermochemistry evolution
# ---------------------------------------------------------------------------


def figure_thermochemistry_evolution(m5) -> Path:
    firing = m5.firing
    time = m5.time_s[firing]

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.8), sharex=True)
    fig.subplots_adjust(left=0.085, right=0.975, top=0.915, bottom=0.245, hspace=0.20, wspace=0.25)
    fig.suptitle(
        "Final Figure 2 — which property drives the added Milestone 5 thrust decay",
        fontsize=12.5,
    )
    ax_of, ax_cstar, ax_tc, ax_gamma = axes[0][0], axes[0][1], axes[1][0], axes[1][1]

    ax_of.plot(time, m5.mixture_ratio[firing], color=M5_COLOUR, zorder=4)
    ax_of.set_ylabel(r"mixture ratio  $O/F$   [-]")

    ax_cstar.plot(
        time,
        m5.ideal_c_star_m_s[firing],
        color=SECONDARY,
        linestyle="--",
        zorder=4,
        label=r"ideal  $c^*$",
    )
    ax_cstar.plot(
        time, m5.c_star_m_s[firing], color=M5_COLOUR, zorder=5, label=r"delivered  $\eta_{c^*} c^*$"
    )
    ax_cstar.set_ylabel(r"characteristic velocity  $c^*$   [m/s]")
    ax_cstar.legend(loc="upper right")

    ax_tc.plot(time, m5.chamber_temperature_k[firing], color=M4_COLOUR, zorder=4)
    ax_tc.set_ylabel(r"chamber temperature  $T_c$   [K]")
    ax_tc.set_xlabel("time  $t$   [s]")

    ax_gamma.plot(time, m5.gamma[firing], color=ACCENT, zorder=4)
    ax_gamma.set_ylabel(r"ratio of specific heats  $\gamma$   [-]")
    ax_gamma.set_xlabel("time  $t$   [s]")

    for axis in (ax_of, ax_cstar, ax_tc, ax_gamma):
        axis.set_xlim(0.0, float(time[-1]))

    def swing(values):
        return 100.0 * (values[-1] / values[0] - 1.0)

    ax_gamma.annotate(
        f"{swing(m5.gamma[firing]):+.2f} % over the burn",
        xy=(0.5, 0.12),
        xycoords="axes fraction",
        ha="center",
        fontsize=8.6,
        bbox={"boxstyle": "round,pad=0.36", "facecolor": "white", "edgecolor": "0.7"},
    )
    ax_cstar.annotate(
        f"{swing(m5.c_star_m_s[firing]):+.2f} % over the burn",
        xy=(0.04, 0.10),
        xycoords="axes fraction",
        ha="left",
        fontsize=8.6,
        bbox={"boxstyle": "round,pad=0.36", "facecolor": "white", "edgecolor": "0.7"},
    )

    ratio = abs(swing(m5.c_star_m_s[firing]) / swing(m5.gamma[firing]))
    fig.text(
        0.085,
        0.185,
        "As the port opens the mixture ratio falls, and on the fuel-rich side of the CEA table "
        "a falling $O/F$ means a falling $c^*$\n"
        f"and a falling $T_c$. Over this burn $c^*$ moves {abs(swing(m5.c_star_m_s[firing])):.2f} %"
        f" while $\\gamma$ moves only {abs(swing(m5.gamma[firing])):.2f} % — a factor of "
        f"{ratio:.0f}.\n"
        "The added decay is therefore a $c^*$ effect acting through chamber pressure, not a "
        "nozzle-$\\gamma$ effect. The ideal and\n"
        r"delivered curves differ only by the constant $\eta_{c^*}$, applied once.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.085, 0.048, DISCLAIMER, fontsize=7.5, color="0.42", va="top")
    return _save(fig, "final_thermochemistry_evolution.png")


# ---------------------------------------------------------------------------
# Figure 3 -- sensitivity ranking
# ---------------------------------------------------------------------------


def _build_families(baseline):
    return [
        SensitivityFamily(
            "regression coefficient $a$",
            "",
            baseline,
            (
                variable_case("a x0.90", "", law=scaled_regression_law(LAW, 0.90)),
                variable_case("a x1.10", "", law=scaled_regression_law(LAW, 1.10)),
            ),
        ),
        SensitivityFamily(
            "regression exponent $n$",
            "",
            baseline,
            (
                variable_case(
                    "n 0.3367",
                    "",
                    law=law_with_exponent(LAW, 0.3367, anchor_flux_si=ANCHOR_FLUX_SI),
                ),
                variable_case(
                    "n 0.3967",
                    "",
                    law=law_with_exponent(LAW, 0.3967, anchor_flux_si=ANCHOR_FLUX_SI),
                ),
            ),
        ),
        SensitivityFamily(
            r"$c^*$ efficiency  $\eta_{c^*}$",
            "",
            baseline,
            (
                variable_case("eta 0.92", "", efficiency=0.92),
                variable_case("eta 1.00", "", efficiency=1.00),
            ),
        ),
        SensitivityFamily(
            "injector effective area  $A_{eff}$",
            "",
            baseline,
            (
                variable_case("A x0.90", "", area=0.90 * AREA_M2),
                variable_case("A x1.10", "", area=1.10 * AREA_M2),
            ),
        ),
        SensitivityFamily(
            "injector model (SPI limit)",
            "",
            baseline,
            (variable_case("SPI", "", model=InjectorModel.SPI),),
        ),
        SensitivityFamily(
            "initial tank temperature  $T_{tank}$",
            "",
            baseline,
            (
                variable_case("283.15 K", "", temperature=283.15),
                variable_case("303.15 K", "", temperature=303.15),
            ),
        ),
        SensitivityFamily(
            "chemistry grid resolution",
            "",
            baseline,
            (
                variable_case("grid /4", "", table=TABLE.subsample(of_stride=4, pressure_stride=2)),
                variable_case(
                    "grid /16", "", table=TABLE.subsample(of_stride=16, pressure_stride=7)
                ),
            ),
        ),
        SensitivityFamily(
            r"chemistry $p_c$ dependence",
            "",
            baseline,
            (variable_case("O/F only", "", table=TABLE.at_fixed_pressure(25.0e5)),),
        ),
    ]


def figure_sensitivity_ranking(families) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 6.4))
    fig.subplots_adjust(left=0.215, right=0.975, top=0.865, bottom=0.275, wspace=0.52)
    fig.suptitle(
        "Final Figure 3 — deterministic one-factor sensitivity — NOT confidence intervals",
        fontsize=12.5,
    )

    for axis, metric, title, unit in (
        (axes[0], "total_impulse_n_s", "effect on total impulse", "% change from baseline"),
        (axes[1], "thrust_slope_percent", "effect on thrust decay", "% change in the decay figure"),
    ):
        ranking = rank_families(families, metric)
        names = [name for name, _s in ranking.entries]
        lookup = {f.name: f for f in families}
        lows = [lookup[n].signed_extremes(metric)[0] for n in names]
        highs = [lookup[n].signed_extremes(metric)[1] for n in names]
        positions = np.arange(len(names))[::-1]

        for position, low, high in zip(positions, lows, highs, strict=True):
            left, right = min(low, 0.0), max(high, 0.0)
            axis.barh(
                position,
                right - left,
                left=left,
                height=0.62,
                color=BAND,
                edgecolor=M5_COLOUR,
                linewidth=1.1,
                zorder=3,
            )
        axis.axvline(0.0, color="0.25", linewidth=1.2, zorder=4)
        axis.set_yticks(positions, names, fontsize=9)
        axis.set_xlabel(unit)
        axis.set_title(title, fontsize=10.5)
        span = max(max(abs(v) for v in lows), max(abs(v) for v in highs)) * 1.48
        axis.set_xlim(-span, span)
        for position, low, high in zip(positions, lows, highs, strict=True):
            widest = high if abs(high) >= abs(low) else low
            offset = 6 if widest >= 0 else -6
            axis.annotate(
                f"{widest:+.2f}",
                xy=(widest, position),
                fontsize=8.0,
                xytext=(offset, 0),
                textcoords="offset points",
                ha="left" if widest >= 0 else "right",
                va="center",
            )

    fig.text(
        0.045,
        0.195,
        "Each bar spans the range produced by the tested alternatives for ONE assumption, with "
        "everything else held at the reference case.\n"
        "These are deterministic re-runs at stated values — a bar is a range of tried options, "
        "not a probability distribution, not a standard\n"
        "deviation and not a confidence interval. No uncertainty distribution is assumed for any "
        "input anywhere in this project.\n"
        "Bars are ordered by the largest absolute change in each panel, computed from the runs "
        "themselves.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.045, 0.045, DISCLAIMER, fontsize=7.5, color="0.42", va="top")
    return _save(fig, "final_sensitivity_ranking.png")


# ---------------------------------------------------------------------------
# Figure 4 -- reference case
# ---------------------------------------------------------------------------


def figure_reference_case(m5, state) -> Path:
    firing = m5.firing
    time = m5.time_s

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2))
    fig.subplots_adjust(left=0.072, right=0.935, top=0.885, bottom=0.315, hspace=0.34, wspace=0.30)
    fig.suptitle("Final Figure 4 — representative conceptual reference case", fontsize=13, y=0.965)
    fig.text(
        0.5,
        0.925,
        "Generic coupled N₂O/HTPB reduced-order model — not experimentally validated",
        ha="center",
        fontsize=9.6,
        color="0.3",
    )

    ax_p, ax_flow, ax_thrust, ax_mass = axes[0][0], axes[0][1], axes[1][0], axes[1][1]

    ax_p.plot(
        time,
        m5.tank_pressure_pa / 1e5,
        color=M5_COLOUR,
        zorder=4,
        label=r"tank  $p_{tank} = p_{sat}(T)$",
    )
    ax_p.plot(
        time,
        m5.chamber_pressure_pa / 1e5,
        color=M4_COLOUR,
        linestyle="--",
        zorder=4,
        label=r"chamber  $p_c$",
    )
    ax_p.fill_between(
        time,
        m5.chamber_pressure_pa / 1e5,
        m5.tank_pressure_pa / 1e5,
        color=BAND,
        alpha=0.5,
        zorder=1,
        label="injector pressure margin",
    )
    ax_p.set_ylabel("pressure   [bar]")
    ax_p.set_xlabel("time  $t$   [s]")
    ax_p.set_ylim(0.0, 56.0)
    ax_p.legend(loc="lower left")

    ax_flow.plot(
        time, m5.oxidizer_mass_flow_kg_s, color=M5_COLOUR, zorder=4, label=r"$\dot{m}_{ox}$"
    )
    ax_flow.plot(time, m5.fuel_mass_flow_kg_s, color=SECONDARY, zorder=4, label=r"$\dot{m}_f$")
    ax_flow.set_ylabel("mass flow   [kg/s]")
    ax_flow.set_xlabel("time  $t$   [s]")
    ax_flow.set_ylim(0.0, 0.125)
    ax_flow.legend(loc="center left")
    ax_of = ax_flow.twinx()
    ax_of.plot(
        time[firing], m5.mixture_ratio[firing], color=ACCENT, linestyle=":", zorder=4, label="$O/F$"
    )
    ax_of.set_ylabel(r"mixture ratio  $O/F$   [-]", color=ACCENT)
    ax_of.tick_params(axis="y", colors=ACCENT)
    ax_of.grid(False)
    ax_of.legend(loc="upper right")

    ax_thrust.plot(time, m5.thrust_n, color=M5_COLOUR, zorder=4)
    ax_thrust.set_ylabel("thrust  $F$   [N]")
    ax_thrust.set_xlabel("time  $t$   [s]")
    ax_thrust.set_ylim(0.0, 330.0)

    ax_mass.plot(
        time, m5.cumulative_impulse_n_s, color=M5_COLOUR, zorder=4, label="cumulative impulse"
    )
    ax_mass.set_ylabel("cumulative impulse   [N s]")
    ax_mass.set_xlabel("time  $t$   [s]")
    ax_mass.legend(loc="upper left")
    ax_consumed = ax_mass.twinx()
    ax_consumed.plot(
        time,
        m5.cumulative_oxidizer_mass_kg,
        color=M4_COLOUR,
        linestyle="--",
        zorder=4,
        label="oxidizer consumed",
    )
    ax_consumed.plot(
        time,
        m5.cumulative_fuel_mass_kg,
        color=SECONDARY,
        linestyle="-.",
        zorder=4,
        label="fuel consumed",
    )
    ax_consumed.set_ylabel("mass consumed   [kg]")
    ax_consumed.grid(False)
    ax_consumed.legend(loc="lower right")

    summary = (
        f"initial tank pressure      {m5.tank_pressure_pa[0] / 1e5:8.3f} bar\n"
        f"initial chamber pressure   {m5.chamber_pressure_pa[0] / 1e5:8.3f} bar\n"
        f"thrust                     {m5.thrust_n[0]:8.1f} → {m5.thrust_n[firing][-1]:.1f} N\n"
        f"burn duration              {m5.burn_time_s:8.3f} s\n"
        f"total impulse              {m5.total_impulse_n_s:8.0f} N s\n"
        f"equivalent $I_{{sp}}$              "
        f"{m5.equivalent_specific_impulse_s:8.1f} s\n"
        f"terminal event             {_readable(m5.termination.value)}"
    )
    fig.text(
        0.078,
        0.205,
        summary,
        fontsize=8.6,
        va="top",
        family="DejaVu Sans",
        bbox={"boxstyle": "round,pad=0.6", "facecolor": "#f4f7fa", "edgecolor": "0.75"},
    )
    fig.text(
        0.40,
        0.205,
        "Illustrative geometry: 0.40 m grain, 40 → 90 mm port, single circular port.\n"
        "Sourced physics: Rezaei et al. (2018) regression law; Lemmon & Span (2006) N₂O EOS;\n"
        "NASA CEA (RP-1311) equilibrium chemistry from a committed table.\n"
        "Effective/calibrated: injector effective flow area 3.50 mm² at $C_d$ = 0.66 — no hole\n"
        "count, no hole size, no plate geometry.\n"
        r"Model-form assumption: a single $\eta_{c^*}$ = 0.96 stands in for every combustion loss."
        "\n"
        "Nothing here was tuned to a performance target, and nothing is validated against a\n"
        "specific motor.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.072, 0.040, DISCLAIMER, fontsize=7.5, color="0.42", va="top")
    return _save(fig, "final_reference_case.png")


# ---------------------------------------------------------------------------
# Figure 5 -- coupling decomposition
# ---------------------------------------------------------------------------


def figure_coupling_decomposition(m3_case, m4_case, m5_case) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 5.4))
    fig.subplots_adjust(left=0.068, right=0.982, top=0.855, bottom=0.315, wspace=0.33)
    fig.suptitle(
        "Final Figure 5 — decomposing the total effect into feed coupling and chemistry",
        fontsize=12.5,
    )

    labels = ["feed coupling\nM3 → M4", "thermochemistry\nM4 → M5", "net\nM3 → M5"]
    colours = [M4_COLOUR, M5_COLOUR, "0.35"]

    def percent(a, b, metric):
        return 100.0 * (getattr(b, metric) / getattr(a, metric) - 1.0)

    panels = (
        (
            axes[0],
            [
                percent(m3_case, m4_case, "total_impulse_n_s"),
                percent(m4_case, m5_case, "total_impulse_n_s"),
                percent(m3_case, m5_case, "total_impulse_n_s"),
            ],
            "total impulse",
            "change   [%]",
            "{:+.2f}",
        ),
        (
            axes[1],
            [
                m4_case.thrust_slope_percent - m3_case.thrust_slope_percent,
                m5_case.thrust_slope_percent - m4_case.thrust_slope_percent,
                m5_case.thrust_slope_percent - m3_case.thrust_slope_percent,
            ],
            "thrust-time slope",
            "change   [percentage points]",
            "{:+.2f}",
        ),
        (
            axes[2],
            [
                percent(m3_case, m4_case, "burn_time_s"),
                percent(m4_case, m5_case, "burn_time_s"),
                percent(m3_case, m5_case, "burn_time_s"),
            ],
            "burn duration",
            "change   [%]",
            "{:+.2f}",
        ),
    )
    for axis, values, title, ylabel, fmt in panels:
        axis.bar(range(3), values, color=colours, width=0.6, zorder=3)
        axis.axhline(0.0, color="0.25", linewidth=1.2, zorder=4)
        axis.set_xticks(range(3), labels, fontsize=8.8)
        axis.set_ylabel(ylabel)
        axis.set_title(title, fontsize=10.5)
        span = max(abs(v) for v in values) * 1.35
        axis.set_ylim(-span, span)
        for position, value in enumerate(values):
            axis.annotate(
                fmt.format(value),
                xy=(position, value),
                ha="center",
                va="bottom" if value >= 0 else "top",
                xytext=(0, 5 if value >= 0 else -5),
                textcoords="offset points",
                fontsize=8.6,
            )

    fig.text(
        0.068,
        0.235,
        "Signed bars: each shows what one modelling step does to the reference case. Coupling "
        "the tank (M3 → M4) is what flips the thrust-time\n"
        "slope from rising to falling — a 16 percentage-point swing — while barely moving total "
        "impulse. Letting the chemistry follow the mixture\n"
        "ratio (M4 → M5) does the opposite: it is the larger effect on total impulse and the "
        "smaller one on slope. The two steps are therefore\n"
        "not interchangeable, and neither alone accounts for the net change. Slope is quoted in "
        "percentage points because it is itself a percentage.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.068, 0.052, DISCLAIMER, fontsize=7.5, color="0.42", va="top")
    return _save(fig, "final_coupling_decomposition.png")


def main() -> int:
    _configure_matplotlib()
    print("Generating Milestone 6 final figures...")
    transient, m3, m3_index, m4, m5, state = _load_models()

    from final_robustness_study import milestone_3_case

    m3_case = milestone_3_case()
    m4_case = summarise("M4", "hierarchy", m4, model=ModelKind.FROZEN_PROPERTIES)
    m5_case = summarise("M5", "hierarchy", m5, model=ModelKind.VARIABLE_PROPERTIES)

    paths = [
        figure_model_hierarchy(transient, m3, m3_index, m4, m5),
        figure_thermochemistry_evolution(m5),
    ]
    families = _build_families(m5_case)
    paths.append(figure_sensitivity_ranking(families))
    paths.append(figure_reference_case(m5, state))
    paths.append(figure_coupling_decomposition(m3_case, m4_case, m5_case))

    for path in paths:
        print(f"  wrote {path.relative_to(FIGURES_DIR.parent)} ({path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
