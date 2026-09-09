"""Deterministic generation of the Milestone 3 figures.

Produces six PNGs in ``figures/``:

* ``fig_m3_a_chamber_pressure.png``   - chamber pressure history
* ``fig_m3_b_thrust.png``             - thrust history (the flagship M3 figure)
* ``fig_m3_c_nozzle_state.png``       - nozzle pressure ratios and expansion regimes
* ``fig_m3_d_ambient_sensitivity.png``- thrust and C_F versus ambient pressure
* ``fig_m3_e_property_sensitivity.png`` - prescribed-property one-factor sensitivity
* ``fig_m3_f_piecewise_thrust.png``   - thrust response to the piecewise schedule

Regeneration is byte-identical for a given matplotlib/FreeType build: the Agg
backend is forced, every rendering rcParam is pinned, no timestamp metadata is
written, and every input is a fixed constant.  PNG glyph rasterisation depends on
the FreeType version, so hashes are not expected to match across environments
with different builds; the figure content is unchanged.  Verified with
matplotlib 3.10.9 / FreeType 2.6.1.

Run with::

    python scripts/generate_m3_figures.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION
from hybrid_rocket_motor.nozzle import (
    REFERENCE_NOZZLE,
    SUMMERFIELD_SEPARATION_RATIO,
    NozzleGeometry,
    area_mach_ratio,
    isentropic_pressure_ratio,
    solve_supersonic_exit_mach,
)
from hybrid_rocket_motor.performance import solve_performance_point

sys.path.insert(0, str(Path(__file__).resolve().parent))

from thrust_prediction_study import (
    AMBIENT_PRESSURES_PA,
    M1_M_DOT_F,
    M1_M_DOT_OX,
    PIECEWISE_STEPS,
    SEA_LEVEL_PA,
    combustion_sensitivity,
    run_case,
)
from transient_regression_study import M_DOT_NOMINAL

from hybrid_rocket_motor.transient import (  # isort: skip
    ConstantOxidizerFlow,
    PiecewiseConstantOxidizerFlow,
)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

COMBUSTION = ILLUSTRATIVE_N2O_HTPB_COMBUSTION
NOZZLE = REFERENCE_NOZZLE

PRIMARY = "#0b5394"
SECONDARY = "#1b9e77"
TERTIARY = "#d95f02"
NEUTRAL = "#7570b3"
ACCENT = "#a11b1b"
BAND_COLOUR = "#bcd7ee"

DISCLAIMER = (
    "Generic reduced-order educational model — NOT experimentally validated, not an optimised "
    "nozzle, not a flight-ready motor.\n"
    "Oxidizer flow is prescribed and does not respond to chamber pressure; "
    r"$c^*$, $\gamma$ and $T_c$ are prescribed constants."
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
            "svg.hashsalt": "hybrid-rocket-motor-m3",
            "path.simplify": False,
            "axes.axisbelow": True,
        }
    )


def _save(fig, name: str) -> Path:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / name
    fig.savefig(path, format="png", metadata={"Software": None})
    plt.close(fig)
    return path


def _load():
    _, _, transient_a, history_a = run_case(
        ConstantOxidizerFlow(M_DOT_NOMINAL, 200.0), "A", "nominal"
    )
    _, _, transient_d, history_d = run_case(
        PiecewiseConstantOxidizerFlow.from_durations(PIECEWISE_STEPS), "D", "piecewise"
    )
    return transient_a, history_a, transient_d, history_d


# ---------------------------------------------------------------------------
# Figure A - chamber pressure history
# ---------------------------------------------------------------------------


def figure_a_chamber_pressure(history, transient) -> Path:
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    fig.subplots_adjust(left=0.10, right=0.975, top=0.870, bottom=0.290)

    ax.plot(
        history.time_s,
        history.chamber_pressure_pa / 1e5,
        color=PRIMARY,
        zorder=4,
        label=r"$p_c = (\dot{m}_{ox} + \dot{m}_f)\,c^*/A_t$",
    )
    ax.axhspan(
        19.9, 31.0, color=BAND_COLOUR, alpha=0.40, zorder=0,
        label=(
            "19.9–31.0 bar measured by Rezaei et al. (2018)\n"
            "(context for plausibility, not validation)"
        ),
    )
    burnout = transient.burnout_time_s
    ax.plot(
        burnout,
        history.chamber_pressure_pa[-1] / 1e5,
        "o",
        color=ACCENT,
        markersize=8,
        markeredgecolor="white",
        markeredgewidth=1.0,
        zorder=6,
        label=f"burnout, {burnout:.2f} s",
    )
    ax.axvline(burnout, color=ACCENT, linestyle="--", linewidth=1.2, zorder=2)
    ax.annotate(
        f"{history.chamber_pressure_pa[0] / 1e5:.3f} bar",
        xy=(0.0, history.chamber_pressure_pa[0] / 1e5),
        xytext=(2.0, 25.2),
        fontsize=9,
        color=PRIMARY,
        arrowprops={"arrowstyle": "-", "color": PRIMARY, "linewidth": 0.9},
    )
    ax.annotate(
        f"{history.chamber_pressure_pa[-1] / 1e5:.3f} bar",
        xy=(burnout, history.chamber_pressure_pa[-1] / 1e5),
        xytext=(burnout - 11.0, 29.3),
        fontsize=9,
        color=PRIMARY,
        arrowprops={"arrowstyle": "-", "color": PRIMARY, "linewidth": 0.9},
    )

    ax.set_xlabel(r"time  $t$   [s]")
    ax.set_ylabel(r"chamber pressure  $p_c$   [bar]")
    ax.set_xlim(0.0, 46.0)
    ax.set_ylim(18.5, 32.5)
    ax.legend(loc="lower left")
    ax.set_title(
        r"Figure M3-A — Chamber pressure, constant prescribed $\dot{m}_{ox} = "
        rf"{M_DOT_NOMINAL:.3f}$ kg/s",
        pad=10,
    )
    fig.text(
        0.5,
        0.150,
        r"Quasi-steady closure with a PRESCRIBED $c^* = "
        rf"{COMBUSTION.c_star_m_s:.1f}$ m/s and a fixed throat "
        rf"$D_t = {NOZZLE.throat_diameter_m * 1e3:.0f}$ mm."
        "\n"
        r"$p_c$ rises only because $\dot{m}_f$ rises as the port opens; "
        "there is no chamber-filling or ignition transient.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.045, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m3_a_chamber_pressure.png")


# ---------------------------------------------------------------------------
# Figure B - thrust history (flagship)
# ---------------------------------------------------------------------------


def figure_b_thrust(history, transient) -> Path:
    """Two panels: the momentum term dominates, so the small pressure term needs
    its own axis to be legible at all."""
    fig, axes = plt.subplots(
        2, 1, figsize=(9.4, 7.2), sharex=True, gridspec_kw={"height_ratios": [2.4, 1.0]}
    )
    fig.subplots_adjust(left=0.100, right=0.975, top=0.885, bottom=0.215, hspace=0.13)
    ax, ax_press = axes
    burnout = transient.burnout_time_s

    ax.plot(history.time_s, history.thrust_n, color=PRIMARY, linewidth=2.6, zorder=5,
            label=r"total thrust  $F = \dot{m}\,V_e + (p_e - p_a)A_e$")
    ax.plot(history.time_s, history.momentum_thrust_n, color=SECONDARY, linestyle="--",
            zorder=4, label=r"momentum thrust  $\dot{m}\,V_e$")
    ax.plot(burnout, history.thrust_n[-1], "o", color=ACCENT, markersize=8,
            markeredgecolor="white", markeredgewidth=1.0, zorder=6,
            label=f"burnout, {burnout:.2f} s")
    for axis in axes:
        axis.axvline(burnout, color=ACCENT, linestyle="--", linewidth=1.2, zorder=2)

    ax.annotate(
        f"{history.thrust_n[0]:.1f} N",
        xy=(0.0, history.thrust_n[0]),
        xytext=(2.6, 302.5),
        fontsize=9, color=PRIMARY,
        arrowprops={"arrowstyle": "-", "color": PRIMARY, "linewidth": 0.9},
    )
    ax.annotate(
        f"peak {history.peak_thrust_n:.1f} N",
        xy=(burnout, history.thrust_n[-1]),
        xytext=(29.0, 338.0),
        fontsize=9, color=PRIMARY,
        arrowprops={"arrowstyle": "-", "color": PRIMARY, "linewidth": 0.9},
    )
    ax.text(
        0.885,
        0.045,
        f"total impulse  {history.total_impulse_n_s:.0f} N\u00b7s\n"
        f"mean thrust  {history.mean_thrust_n:.1f} N\n"
        f"equivalent $I_{{sp}}$  {history.equivalent_specific_impulse_s:.1f} s",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.42", "facecolor": "white", "edgecolor": PRIMARY,
              "alpha": 0.95, "linewidth": 0.9},
    )
    ax.set_ylabel(r"thrust   [N]")
    ax.set_ylim(297.0, 343.0)
    ax.legend(loc="upper left")
    ax.set_title(
        r"Figure M3-B — Predicted thrust, constant prescribed $\dot{m}_{ox} = "
        rf"{M_DOT_NOMINAL:.3f}$ kg/s at $p_a = 101325$ Pa",
        pad=10,
    )

    ax_press.plot(history.time_s, history.pressure_thrust_n, color=TERTIARY, linewidth=2.2,
                  zorder=4, label=r"pressure thrust  $(p_e - p_a)A_e$")
    ax_press.set_xlabel(r"time  $t$   [s]")
    ax_press.set_ylabel(r"[N]")
    ax_press.set_xlim(0.0, 46.0)
    ax_press.set_ylim(3.2, 7.0)
    ax_press.legend(loc="upper left")

    fig.text(
        0.5,
        0.128,
        "Thrust rises through the burn only because the total propellant flow rises: with "
        r"$\varepsilon$, $\gamma$ and $T_c$ prescribed constant,"
        "\n"
        r"$M_e$, $T_e$ and $V_e$ are frozen, so $F$ is exactly affine in $\dot{m}_{total}$. "
        "This is a modelling consequence, not a validated prediction.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.038, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m3_b_thrust.png")


# ---------------------------------------------------------------------------
# Figure C - nozzle state
# ---------------------------------------------------------------------------


def figure_c_nozzle_state(history) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.6))
    fig.subplots_adjust(left=0.075, right=0.975, top=0.835, bottom=0.320, wspace=0.26)
    ax_time, ax_eps = axes

    # -- left: pressure ratios through the burn ---------------------------
    p_e_over_p_c = history.exit_pressure_pa / history.chamber_pressure_pa
    p_a_over_p_c = SEA_LEVEL_PA / history.chamber_pressure_pa
    ax_time.plot(history.time_s, p_e_over_p_c, color=PRIMARY, zorder=4,
                 label=r"$p_e/p_c$  (fixed by $\varepsilon$ and $\gamma$)")
    ax_time.plot(history.time_s, p_a_over_p_c, color=TERTIARY, linestyle="--", zorder=4,
                 label=r"$p_a/p_c$  (falls as $p_c$ rises)")
    ax_time.fill_between(history.time_s, p_a_over_p_c, p_e_over_p_c,
                         where=p_e_over_p_c > p_a_over_p_c, color=BAND_COLOUR, alpha=0.55,
                         zorder=1, label="underexpanded margin")
    ax_time.set_xlabel(r"time  $t$   [s]")
    ax_time.set_ylabel("pressure ratio   [–]")
    ax_time.set_xlim(0.0, 42.0)
    ax_time.set_ylim(0.030, 0.046)
    ax_time.legend(loc="lower left")
    ax_time.set_title("Through the burn", pad=8)

    # -- right: static sweep over expansion ratio -------------------------
    ratios = np.linspace(1.02, 30.0, 400)
    machs = np.array([solve_supersonic_exit_mach(e, COMBUSTION.gamma) for e in ratios])
    p_e = history.chamber_pressure_pa[0] * isentropic_pressure_ratio(machs, COMBUSTION.gamma)
    cf = np.array([
        solve_performance_point(
            COMBUSTION,
            NozzleGeometry.from_throat_and_expansion_ratio(NOZZLE.throat_diameter_m, e),
            M1_M_DOT_OX,
            M1_M_DOT_F,
            SEA_LEVEL_PA,
        ).thrust_coefficient
        for e in ratios
    ])
    ax_eps.plot(ratios, cf, color=PRIMARY, zorder=5, label=r"$C_F$ at the reference chamber state")
    # Ideal expansion solved exactly by inverting the isentropic pressure
    # relation for p_e = p_a, then evaluating the area-Mach relation there.
    # (An np.interp over p_e - p_a would be wrong: that array is DECREASING.)
    p_c_ref = history.chamber_pressure_pa[0]
    gamma = COMBUSTION.gamma
    mach_ideal = math.sqrt(
        ((p_c_ref / SEA_LEVEL_PA) ** ((gamma - 1.0) / gamma) - 1.0) * 2.0 / (gamma - 1.0)
    )
    ideal = area_mach_ratio(mach_ideal, gamma)
    ax_eps.axvline(NOZZLE.expansion_ratio, color=ACCENT, linewidth=1.4, zorder=4,
                   label=rf"reference $\varepsilon = {NOZZLE.expansion_ratio:.0f}$")
    ax_eps.axvline(ideal, color=SECONDARY, linestyle="--", linewidth=1.2, zorder=3,
                   label=rf"ideal expansion, $\varepsilon \approx {ideal:.2f}$")
    mask_sep = p_e < SUMMERFIELD_SEPARATION_RATIO * SEA_LEVEL_PA
    if np.any(mask_sep):
        ax_eps.axvspan(float(ratios[mask_sep][0]), 30.0, color=ACCENT, alpha=0.10, zorder=0,
                       label="Summerfield separation expected\n(NOT modelled here)")
    ax_eps.set_xlabel(r"expansion ratio  $\varepsilon = A_e/A_t$   [–]")
    ax_eps.set_ylabel(r"thrust coefficient  $C_F$   [–]")
    ax_eps.set_xlim(1.0, 30.0)
    ax_eps.set_ylim(1.20, 1.75)
    ax_eps.legend(loc="upper right")
    ax_eps.set_title(r"Static sweep at $p_a = 101325$ Pa", pad=8)

    fig.suptitle(
        "Figure M3-C — Nozzle state: pressure ratios and expansion regime", fontsize=12.5, y=0.955
    )
    fig.text(
        0.5,
        0.170,
        r"The reference nozzle is slightly UNDEREXPANDED at sea level ($p_e > p_a$), so the "
        "pressure term adds thrust. As $p_c$ rises\n"
        "through the burn the nozzle becomes progressively more underexpanded. "
        "The shaded band on the right marks where a\n"
        "1-D attached-flow result would be physically doubtful; "
        "shocks and separation are NOT modelled anywhere in this project.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.050, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m3_c_nozzle_state.png")


# ---------------------------------------------------------------------------
# Figure D - ambient-pressure sensitivity
# ---------------------------------------------------------------------------


def figure_d_ambient_sensitivity() -> Path:
    pressures = np.linspace(0.0, 101325.0, 400)
    points = [
        solve_performance_point(COMBUSTION, NOZZLE, M1_M_DOT_OX, M1_M_DOT_F, float(p))
        for p in pressures
    ]
    thrust = np.array([p.thrust_n for p in points])
    cf = np.array([p.thrust_coefficient for p in points])

    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    fig.subplots_adjust(left=0.095, right=0.885, top=0.870, bottom=0.300)

    ax.plot(pressures / 1e3, thrust, color=PRIMARY, zorder=4, label=r"thrust $F$")
    ax_cf = ax.twinx()
    ax_cf.plot(pressures / 1e3, cf, color=TERTIARY, linestyle="--", zorder=4,
               label=r"thrust coefficient $C_F$")
    ax_cf.set_ylabel(r"thrust coefficient  $C_F$   [–]", color=TERTIARY)
    ax_cf.tick_params(axis="y", labelcolor=TERTIARY)
    ax_cf.set_ylim(1.48, 1.70)
    ax_cf.grid(False)

    # Several sampled points sit within a fraction of a newton of each other, so
    # they are listed in one block rather than annotated individually.
    lines = []
    for p_a, note in AMBIENT_PRESSURES_PA:
        point = solve_performance_point(COMBUSTION, NOZZLE, M1_M_DOT_OX, M1_M_DOT_F, p_a)
        ax.plot(p_a / 1e3, point.thrust_n, "o", color=PRIMARY, markersize=6,
                markeredgecolor="white", markeredgewidth=0.8, zorder=6)
        lines.append(f"{p_a / 1e3:>7.2f} kPa   {point.thrust_n:>7.2f} N   ({note})")
    ax.text(
        0.025,
        0.045,
        "sampled points\n" + "\n".join(lines),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        family="DejaVu Sans Mono",
        color=PRIMARY,
        bbox={"boxstyle": "round,pad=0.42", "facecolor": "white", "edgecolor": PRIMARY,
              "alpha": 0.95, "linewidth": 0.9},
    )

    ax.set_xlabel(r"ambient pressure  $p_a$   [kPa]")
    ax.set_ylabel(r"thrust  $F$   [N]", color=PRIMARY)
    ax.tick_params(axis="y", labelcolor=PRIMARY)
    ax.set_xlim(-3.0, 112.0)
    ax.set_ylim(303.0, 350.0)
    handles = ax.get_lines()[:1] + ax_cf.get_lines()[:1]
    ax.legend(handles, [h.get_label() for h in handles], loc="upper right")
    ax.set_title(
        "Figure M3-D — Deterministic nozzle sensitivity to ambient pressure "
        "at the fixed reference chamber state",
        pad=10,
    )
    fig.text(
        0.5,
        0.160,
        r"The chamber and exit states never change: the nozzle is choked, so only the "
        r"$(p_e - p_a)A_e$ term moves."
        "\n"
        r"Sea level to vacuum gains exactly $p_a A_e = "
        rf"{SEA_LEVEL_PA * NOZZLE.exit_area_m2:.2f}$ N."
        "\n"
        "This is a pressure sweep, NOT an altitude-performance envelope — "
        "no atmosphere model exists anywhere here.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.048, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m3_d_ambient_sensitivity.png")


# ---------------------------------------------------------------------------
# Figure E - prescribed-property sensitivity
# ---------------------------------------------------------------------------


def figure_e_property_sensitivity() -> Path:
    """Thrust against chamber pressure.

    Plotting the change in specific impulse alongside the change in thrust would
    be redundant: at fixed mass flow ``I_sp = F/(m_dot g_0)``, so the two
    percentage changes are identical by construction.  The chamber pressure is
    plotted instead, because it responds very differently -- a 10 % change in
    ``c*`` moves ``p_c`` by exactly 10 % but thrust by only about 1 %.
    """
    rows = combustion_sensitivity()
    base = rows[0][2]
    labels, thrust_delta, pressure_delta = [], [], []
    for label, _props, result, _b in rows[1:]:
        labels.append(label)
        thrust_delta.append(100.0 * (result.thrust_n / base.thrust_n - 1.0))
        pressure_delta.append(
            100.0 * (result.chamber_pressure_pa / base.chamber_pressure_pa - 1.0)
        )

    order = np.argsort(np.abs(thrust_delta))
    labels = [labels[i] for i in order]
    thrust_delta = [thrust_delta[i] for i in order]
    pressure_delta = [pressure_delta[i] for i in order]

    y = np.arange(len(labels), dtype=float)
    height = 0.38

    fig, ax = plt.subplots(figsize=(9.6, 6.2))
    fig.subplots_adjust(left=0.155, right=0.975, top=0.880, bottom=0.300)

    ax.barh(y + height / 2.0, thrust_delta, height=height, color=PRIMARY, zorder=3,
            label=r"change in thrust $F$  (and, identically, in $I_{sp}$)")
    ax.barh(y - height / 2.0, pressure_delta, height=height, color=TERTIARY, zorder=3,
            label=r"change in chamber pressure $p_c$")
    ax.axvline(0.0, color="black", linewidth=1.0, zorder=4)

    for value, position in zip(thrust_delta, y + height / 2.0, strict=True):
        ax.annotate(f"{value:+.2f}%", xy=(value, position),
                    xytext=(4 if value >= 0 else -4, 0), textcoords="offset points",
                    va="center", ha="left" if value >= 0 else "right", fontsize=8.5,
                    color=PRIMARY)
    for value, position in zip(pressure_delta, y - height / 2.0, strict=True):
        ax.annotate(f"{value:+.1f}%", xy=(value, position),
                    xytext=(4 if value >= 0 else -4, 0), textcoords="offset points",
                    va="center", ha="left" if value >= 0 else "right", fontsize=8.5,
                    color=TERTIARY)

    pretty = {
        "c* -10%": r"$c^*$  $-10\%$",
        "c* +10%": r"$c^*$  $+10\%$",
        "gamma 1.15": r"$\gamma = 1.15$",
        "gamma 1.25": r"$\gamma = 1.25$",
        "T_c -10%": r"$T_c$  $-10\%$",
        "T_c +10%": r"$T_c$  $+10\%$",
    }
    ax.set_yticks(y)
    ax.set_yticklabels([pretty[label] for label in labels])
    ax.set_xlabel("change from the baseline prescribed properties   [%]")
    ax.set_xlim(-14.5, 14.5)
    ax.set_ylim(-0.68, len(labels) - 0.32)
    ax.legend(loc="lower right")
    ax.set_title(
        "Figure M3-E — One-factor sensitivity to the PRESCRIBED combustion properties",
        pad=10,
    )
    fig.text(
        0.5,
        0.165,
        r"Deterministic one-at-a-time variation, NOT uncertainty propagation. $c^*$ moves "
        r"$p_c$ proportionally yet barely moves $F$,"
        "\n"
        r"because $V_e$ is unchanged and only the small pressure term scales; $T_c$ moves "
        r"$V_e$ as $\sqrt{T_c}$ and dominates thrust;"
        "\n"
        r"$\gamma$ shifts the area–Mach solution but leaves $p_c$ untouched. "
        "Mass flow and nozzle geometry are held fixed throughout.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.048, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m3_e_property_sensitivity.png")


# ---------------------------------------------------------------------------
# Figure F - piecewise thrust response
# ---------------------------------------------------------------------------


def figure_f_piecewise(history) -> Path:
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.4), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.975, top=0.900, bottom=0.190, hspace=0.16)
    ax_flow, ax_pc, ax_f = axes

    coast = ~history.firing
    coast_start = float(history.time_s[coast][0])
    coast_end = float(history.time_s[coast][-1])
    for axis in axes:
        axis.axvspan(coast_start, coast_end, color=ACCENT, alpha=0.10, zorder=1,
                     label="zero-flow coast")

    ax_flow.plot(history.time_s, history.oxidizer_mass_flow_kg_s, color=PRIMARY, zorder=4,
                 label=r"prescribed $\dot{m}_{ox}(t)$")
    ax_flow.plot(history.time_s, history.total_mass_flow_kg_s, color=SECONDARY, linestyle="--",
                 zorder=4, label=r"total $\dot{m}_{ox} + \dot{m}_f$")
    ax_flow.set_ylabel(r"$\dot{m}$   [kg/s]")
    ax_flow.set_ylim(-0.018, 0.285)
    ax_flow.legend(loc="upper left", ncols=3)
    ax_flow.set_title(
        "Figure M3-F — Thrust response to the piecewise prescribed oxidizer schedule", pad=10
    )

    ax_pc.plot(history.time_s, history.chamber_pressure_pa / 1e5, color=TERTIARY, zorder=4,
               label=r"chamber pressure $p_c$")
    ax_pc.set_ylabel(r"$p_c$   [bar]")
    ax_pc.set_ylim(-2.5, 46.0)
    ax_pc.legend(loc="upper left", ncols=2)

    ax_f.plot(history.time_s, history.thrust_n, color=PRIMARY, linewidth=2.4, zorder=5,
              label=r"total thrust $F$")
    ax_f.plot(history.time_s, history.momentum_thrust_n, color=SECONDARY, linestyle="--",
              zorder=4, label=r"momentum $\dot{m}\,V_e$")
    ax_f.plot(history.time_s, history.pressure_thrust_n, color=NEUTRAL, linestyle=":",
              linewidth=2.0, zorder=4, label=r"pressure $(p_e - p_a)A_e$")
    ax_f.set_xlabel(r"time  $t$   [s]")
    ax_f.set_ylabel(r"thrust   [N]")
    ax_f.set_xlim(0.0, 40.0)
    ax_f.set_ylim(-30.0, 570.0)
    ax_f.legend(loc="upper left", ncols=4)

    fig.text(
        0.5,
        0.108,
        r"During the coast $\dot{m}_{total} = 0$, so $p_c$ and $F$ are exactly zero and no "
        "nozzle exit state is reported at all.\n"
        "The steps are instantaneous because the model has no chamber-filling, ignition or "
        "tail-off dynamics —\n"
        "a real motor would not respond this way.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.030, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m3_f_piecewise_thrust.png")


def main() -> int:
    _configure_matplotlib()
    print("Generating Milestone 3 figures (deterministic)...")
    transient_a, history_a, _transient_d, history_d = _load()
    outputs = [
        figure_a_chamber_pressure(history_a, transient_a),
        figure_b_thrust(history_a, transient_a),
        figure_c_nozzle_state(history_a),
        figure_d_ambient_sensitivity(),
        figure_e_property_sensitivity(),
        figure_f_piecewise(history_d),
    ]
    for path in outputs:
        print(f"  wrote {path.relative_to(FIGURES_DIR.parent)}  ({path.stat().st_size} bytes)")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
