"""Deterministic generation of the Milestone 2 figures.

Produces four PNGs in ``figures/``:

* ``fig_m2_a_port_growth.png``      - port diameter and remaining web versus time
* ``fig_m2_b_flux_and_rate.png``    - oxidizer mass flux and regression rate histories
* ``fig_m2_c_fuel_flow_and_of.png`` - fuel mass flow and mixture ratio histories
* ``fig_m2_d_mass_closure.png``     - fuel mass bookkeeping and conservation residual
* ``fig_m2_e_piecewise.png``        - response to the piecewise prescribed schedule

Regeneration is byte-identical for a given matplotlib/FreeType build: the Agg
backend is forced, every rendering rcParam is pinned, no timestamp metadata is
written, and every input is a fixed constant.  PNG glyph rasterisation depends on
the FreeType version, so hashes are not expected to match across environments
with different builds; the figure content is unchanged.  Verified with
matplotlib 3.10.9 / FreeType 2.6.1.

Run with::

    python scripts/generate_m2_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from hybrid_rocket_motor import REPRESENTATIVE_GRAIN, REZAEI_2018_N2O_HTPB
from hybrid_rocket_motor.transient import FluxRangeStatus

sys.path.insert(0, str(Path(__file__).resolve().parent))

from transient_regression_study import (
    M_DOT_HIGH,
    M_DOT_LOW,
    M_DOT_NOMINAL,
    PIECEWISE_STEPS,
    run_constant_case,
    run_piecewise_case,
)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB
R_OUTER = GRAIN.outer_diameter_m / 2.0

CASE_COLOURS = ("#0b5394", "#1b9e77", "#d95f02")
BAND_COLOUR = "#bcd7ee"
ACCENT = "#a11b1b"
NEUTRAL = "#7570b3"

DISCLAIMER_MODEL = (
    "Generic reduced-order educational model — no validation against any real motor.\n"
    "Milestone 2: prescribed oxidizer flow only; no chamber-pressure, nozzle or thrust "
    "result is computed."
)


def _configure_matplotlib() -> None:
    """Pin every rcParam that could make output non-deterministic or unreadable."""
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
            "svg.hashsalt": "hybrid-rocket-motor-m2",
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


def _in_range_mask(result) -> np.ndarray:
    return np.array(
        [s is FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE for s in result.flux_status], dtype=bool
    )


def _masked(values: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """Copy with the discarded samples set to NaN so the line breaks cleanly."""
    out = np.array(values, dtype=float)
    out[~keep] = np.nan
    return out


def _grown(mask: np.ndarray) -> np.ndarray:
    """Extend a mask by one sample on each side.

    Used for the out-of-range strand so that the solid and dotted parts of a
    curve meet instead of leaving a visible gap at the crossing.
    """
    return (
        mask
        | np.concatenate(([False], mask[:-1]))
        | np.concatenate((mask[1:], [False]))
    )


def _load_cases():
    case_a, _ = run_constant_case(M_DOT_NOMINAL)
    case_b, _ = run_constant_case(M_DOT_LOW)
    case_c, _ = run_constant_case(M_DOT_HIGH)
    case_d = run_piecewise_case()
    return case_a, case_b, case_c, case_d


# ---------------------------------------------------------------------------
# Figure A - port growth and remaining web
# ---------------------------------------------------------------------------


def figure_a_port_growth(cases) -> Path:
    case_a, case_b, case_c, _ = cases
    runs = ((case_a, M_DOT_NOMINAL), (case_b, M_DOT_LOW), (case_c, M_DOT_HIGH))
    # Distinct vertical offsets: the burnout markers are close together in time,
    # so a single offset would overlap the labels.
    label_offsets = (-17, -17, -33)

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.4), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.975, top=0.885, bottom=0.185, hspace=0.14)
    ax_d, ax_web = axes

    ax_d.axhline(
        GRAIN.outer_diameter_m * 1e3,
        color=ACCENT,
        linestyle="--",
        linewidth=1.4,
        zorder=2,
        label=rf"outer grain diameter $D_o = {GRAIN.outer_diameter_m * 1e3:.0f}$ mm (burnout)",
    )
    for (run, m_dot), colour, dy in zip(runs, CASE_COLOURS, label_offsets, strict=True):
        ax_d.plot(
            run.time_s,
            run.port_diameter_m * 1e3,
            color=colour,
            zorder=4,
            label=rf"$\dot{{m}}_{{ox}} = {m_dot:.3f}$ kg/s",
        )
        ax_d.plot(
            run.burnout_time_s,
            run.port_diameter_m[-1] * 1e3,
            "o",
            color=colour,
            markersize=7,
            markeredgecolor="white",
            markeredgewidth=1.0,
            zorder=6,
        )
        ax_d.annotate(
            rf"burnout {run.burnout_time_s:.2f} s",
            xy=(run.burnout_time_s, run.port_diameter_m[-1] * 1e3),
            xytext=(-6, dy),
            textcoords="offset points",
            ha="right",
            fontsize=8.5,
            color=colour,
        )
    ax_d.set_ylabel(r"port diameter  $D_p$   [mm]")
    ax_d.set_ylim(36.0, 97.0)
    ax_d.legend(loc="lower right")
    ax_d.set_title(
        "Figure M2-A — Port growth under constant prescribed oxidizer flow", pad=10
    )

    for (run, m_dot), colour in zip(runs, CASE_COLOURS, strict=True):
        ax_web.plot(
            run.time_s,
            run.web_thickness_m * 1e3,
            color=colour,
            zorder=4,
            label=rf"$\dot{{m}}_{{ox}} = {m_dot:.3f}$ kg/s",
        )
    ax_web.axhline(0.0, color=ACCENT, linestyle="--", linewidth=1.4, zorder=2,
                   label="fuel exhausted (zero web)")
    ax_web.set_xlabel(r"time  $t$   [s]")
    ax_web.set_ylabel(r"remaining radial web   [mm]")
    ax_web.set_xlim(0.0, 56.0)
    ax_web.set_ylim(-1.6, 27.0)
    ax_web.legend(loc="upper right")

    fig.text(
        0.5,
        0.105,
        r"$\mathrm{d}r_p/\mathrm{d}t = \dot{r}$, so $\mathrm{d}D_p/\mathrm{d}t = 2\dot{r}$. "
        "Integration terminates exactly at the outer grain boundary;\n"
        "the solution is never continued into negative fuel thickness.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.030, DISCLAIMER_MODEL, ha="center", va="center", fontsize=8.5, color="0.25")
    return _save(fig, "fig_m2_a_port_growth.png")


# ---------------------------------------------------------------------------
# Figure B - flux and regression-rate history
# ---------------------------------------------------------------------------


def figure_b_flux_and_rate(cases) -> Path:
    case_a, case_b, case_c, _ = cases
    runs = ((case_a, M_DOT_NOMINAL), (case_b, M_DOT_LOW), (case_c, M_DOT_HIGH))
    low, high = LAW.valid_flux_range_si

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.6), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.975, top=0.885, bottom=0.205, hspace=0.14)
    ax_flux, ax_rate = axes

    ax_flux.axhspan(
        low,
        high,
        color=BAND_COLOUR,
        alpha=0.45,
        zorder=0,
        label=rf"source data range, {low:.0f}–{high:.0f} kg/(m$^2\,$s)",
    )
    ax_flux.axhline(low, color=ACCENT, linestyle="--", linewidth=1.1, zorder=2)

    for (run, m_dot), colour in zip(runs, CASE_COLOURS, strict=True):
        inside = _in_range_mask(run)
        ax_flux.plot(
            run.time_s,
            _masked(run.oxidizer_mass_flux_si, inside),
            color=colour,
            zorder=4,
            label=rf"$\dot{{m}}_{{ox}} = {m_dot:.3f}$ kg/s (within source range)",
        )
        ax_flux.plot(
            run.time_s,
            _masked(run.oxidizer_mass_flux_si, _grown(~inside)),
            color=colour,
            linestyle=":",
            linewidth=1.5,
            zorder=3,
        )
        ax_rate.plot(run.time_s, _masked(run.regression_rate_mm_s, inside), color=colour, zorder=4)
        ax_rate.plot(
            run.time_s,
            _masked(run.regression_rate_mm_s, _grown(~inside)),
            color=colour,
            linestyle=":",
            linewidth=1.5,
            zorder=3,
        )

    ax_flux.plot([], [], color="0.35", linestyle=":", linewidth=1.5,
                 label="extrapolated below the source range")
    ax_flux.set_ylabel(r"oxidizer mass flux  $G_{ox}$   [kg/(m$^2\,$s)]")
    ax_flux.set_ylim(0.0, 132.0)
    ax_flux.legend(loc="upper right")
    ax_flux.set_title(
        "Figure M2-B — Oxidizer mass flux and regression rate histories", pad=10
    )

    ax_rate.set_xlabel(r"time  $t$   [s]")
    ax_rate.set_ylabel(r"regression rate  $\dot{r}$   [mm/s]")
    ax_rate.set_xlim(0.0, 56.0)
    ax_rate.set_ylim(0.30, 1.02)
    ax_rate.plot([], [], color="0.35", linewidth=1.8, label="within source range")
    ax_rate.plot([], [], color="0.35", linestyle=":", linewidth=1.5,
                 label="extrapolated (dotted)")
    ax_rate.legend(loc="upper right")

    fractions = [
        100.0 * run.time_fraction_with_status(FluxRangeStatus.WITHIN_SOURCE_FLUX_RANGE)
        for run, _ in runs
    ]
    fig.text(
        0.5,
        0.118,
        "Dotted segments lie OUTSIDE the flux range over which the correlation was measured and "
        "are NOT validated.\n"
        f"Fraction of each burn inside the source range: Case A {fractions[0]:.1f} %, "
        f"Case B {fractions[1]:.1f} %, Case C {fractions[2]:.1f} %.",
        ha="center",
        va="center",
        fontsize=9.5,
        color=ACCENT,
        fontweight="bold",
    )
    fig.text(0.5, 0.030, DISCLAIMER_MODEL, ha="center", va="center", fontsize=8.5, color="0.25")
    return _save(fig, "fig_m2_b_flux_and_rate.png")


# ---------------------------------------------------------------------------
# Figure C - fuel mass flow and mixture ratio
# ---------------------------------------------------------------------------


def figure_c_fuel_flow_and_of(cases) -> Path:
    case_a, case_b, case_c, case_d = cases
    runs = (
        (case_a, rf"Case A: $\dot{{m}}_{{ox}} = {M_DOT_NOMINAL:.3f}$ kg/s", CASE_COLOURS[0], "-"),
        (case_b, rf"Case B: $\dot{{m}}_{{ox}} = {M_DOT_LOW:.3f}$ kg/s", CASE_COLOURS[1], "-"),
        (case_c, rf"Case C: $\dot{{m}}_{{ox}} = {M_DOT_HIGH:.3f}$ kg/s", CASE_COLOURS[2], "-"),
        (case_d, "Case D: piecewise schedule", NEUTRAL, "--"),
    )

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.6), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.975, top=0.885, bottom=0.195, hspace=0.14)
    ax_fuel, ax_of = axes

    for run, label, colour, style in runs:
        ax_fuel.plot(
            run.time_s, run.fuel_mass_flow_kg_s, color=colour, linestyle=style, zorder=4,
            label=label
        )
        ax_of.plot(run.time_s, run.mixture_ratio, color=colour, linestyle=style, zorder=4)

    # Mark the interval where O/F is genuinely undefined; nothing is drawn there.
    coast = case_d.oxidizer_mass_flow_kg_s == 0.0
    coast_start = float(case_d.time_s[coast][0])
    coast_end = float(case_d.time_s[coast][-1])
    for axis in (ax_fuel, ax_of):
        axis.axvspan(coast_start, coast_end, color=ACCENT, alpha=0.10, zorder=1)
    ax_of.annotate(
        "Case D zero-flow coast:\n"
        r"$\dot{m}_f = 0$, so O/F is undefined"
        "\nand nothing is plotted here",
        xy=(0.5 * (coast_start + coast_end), 1.92),
        xytext=(18.0, 1.56),
        fontsize=8.5,
        color=ACCENT,
        ha="left",
        arrowprops={"arrowstyle": "-", "color": ACCENT, "linewidth": 0.9},
        bbox={"boxstyle": "round,pad=0.34", "facecolor": "white", "edgecolor": ACCENT,
              "alpha": 0.93, "linewidth": 0.8},
    )

    ax_fuel.set_ylabel(r"fuel mass flow  $\dot{m}_f$   [kg/s]")
    ax_fuel.set_ylim(0.0, 0.062)
    ax_fuel.legend(loc="lower right", ncols=2)
    ax_fuel.set_title("Figure M2-C — Fuel mass flow and mixture-ratio histories", pad=10)

    ax_of.set_xlabel(r"time  $t$   [s]")
    ax_of.set_ylabel(r"mixture ratio  O/F   [–]")
    ax_of.set_xlim(0.0, 56.0)
    ax_of.set_ylim(1.20, 3.40)

    fig.text(
        0.5,
        0.112,
        r"At constant $\dot{m}_{ox}$ the burning area grows as $r_p$ while $\dot{r}$ falls only as "
        r"$r_p^{-2n}$ with $n = 0.3667$,"
        "\n"
        r"so $\dot{m}_f$ rises even as $\dot{r}$ falls and O/F drifts downward. "
        "Bookkeeping only — O/F is not used to compute performance.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.030, DISCLAIMER_MODEL, ha="center", va="center", fontsize=8.5, color="0.25")
    return _save(fig, "fig_m2_c_fuel_flow_and_of.png")


# ---------------------------------------------------------------------------
# Figure D - mass bookkeeping and closure
# ---------------------------------------------------------------------------


def figure_d_mass_closure(cases) -> Path:
    case_a = cases[0]

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.4), sharex=True,
                             gridspec_kw={"height_ratios": [2.0, 1.0]})
    fig.subplots_adjust(left=0.125, right=0.975, top=0.885, bottom=0.185, hspace=0.16)
    ax_mass, ax_res = axes

    ax_mass.plot(
        case_a.time_s,
        case_a.remaining_fuel_mass_kg,
        color=CASE_COLOURS[0],
        zorder=4,
        label=r"remaining fuel mass, $\rho_f L \pi (r_o^2 - r_p^2)$",
    )
    ax_mass.plot(
        case_a.time_s,
        case_a.geometric_fuel_consumed_kg,
        color=CASE_COLOURS[2],
        zorder=4,
        label=r"fuel consumed from geometry, $\rho_f L \pi (r_p^2 - r_{p,0}^2)$",
    )
    ax_mass.plot(
        case_a.time_s,
        case_a.cumulative_fuel_mass_kg,
        color="black",
        linestyle="--",
        linewidth=1.4,
        zorder=5,
        label=r"fuel consumed from integration, $\int_0^t \dot{m}_f\,\mathrm{d}t$",
    )
    ax_mass.axhline(
        GRAIN.initial_fuel_mass_kg,
        color=NEUTRAL,
        linestyle=":",
        linewidth=1.3,
        zorder=2,
        label=rf"loaded fuel mass, {GRAIN.initial_fuel_mass_kg:.4f} kg",
    )
    ax_mass.set_ylabel(r"fuel mass   [kg]")
    ax_mass.set_ylim(-0.08, 2.82)
    ax_mass.legend(loc="upper right")
    ax_mass.set_title(
        rf"Figure M2-D — Fuel mass bookkeeping and conservation, Case A "
        rf"($\dot{{m}}_{{ox}} = {M_DOT_NOMINAL:.3f}$ kg/s)",
        pad=10,
    )

    residual = np.abs(case_a.fuel_mass_closure_residual_kg)
    floor = 1.0e-16
    ax_res.semilogy(
        case_a.time_s,
        np.maximum(residual, floor),
        color=ACCENT,
        zorder=4,
        label=r"$\left|\int_0^t \dot{m}_f\,\mathrm{d}t - \rho_f L \pi (r_p^2 - r_{p,0}^2)\right|$",
    )
    ax_res.set_xlabel(r"time  $t$   [s]")
    ax_res.set_ylabel("closure residual   [kg]")
    ax_res.set_xlim(0.0, 42.5)
    ax_res.set_ylim(1.0e-16, 1.0e-8)
    ax_res.legend(loc="lower right")
    ax_res.annotate(
        f"worst residual {case_a.max_absolute_fuel_closure_residual_kg:.2e} kg\n"
        f"({case_a.max_relative_fuel_closure_residual:.1e} relative)",
        xy=(0.03, 0.80),
        xycoords="axes fraction",
        fontsize=8.5,
        color=ACCENT,
        va="top",
    )

    fig.text(
        0.5,
        0.104,
        r"$\dot{m}_f = \rho_f (2\pi r_p L)\,\mathrm{d}r_p/\mathrm{d}t "
        r"= \rho_f L \pi\,\mathrm{d}(r_p^2)/\mathrm{d}t$, so the two curves coincide exactly in "
        "exact arithmetic.\n"
        "The residual is therefore a pure measure of integration error, not a modelling "
        "approximation.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.030, DISCLAIMER_MODEL, ha="center", va="center", fontsize=8.5, color="0.25")
    return _save(fig, "fig_m2_d_mass_closure.png")


# ---------------------------------------------------------------------------
# Figure E - piecewise schedule response
# ---------------------------------------------------------------------------


def figure_e_piecewise(cases) -> Path:
    case_d = cases[3]

    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.6), sharex=True)
    fig.subplots_adjust(left=0.100, right=0.900, top=0.905, bottom=0.215, hspace=0.16)
    ax_flow, ax_port, ax_of = axes

    coast = case_d.oxidizer_mass_flow_kg_s == 0.0
    coast_start = float(case_d.time_s[coast][0])
    coast_end = float(case_d.time_s[coast][-1])
    for axis in axes:
        axis.axvspan(coast_start, coast_end, color=ACCENT, alpha=0.10, zorder=1,
                     label="zero-flow coast")

    ax_flow.plot(case_d.time_s, case_d.oxidizer_mass_flow_kg_s, color=CASE_COLOURS[0], zorder=4,
                 label=r"prescribed $\dot{m}_{ox}(t)$")
    ax_flow.set_ylabel(r"$\dot{m}_{ox}$   [kg/s]")
    ax_flow.set_ylim(-0.012, 0.176)
    handles, labels = ax_flow.get_legend_handles_labels()
    ax_flow.legend(handles[::-1], labels[::-1], loc="lower right", ncols=2)
    ax_flow.set_title(
        "Figure M2-E — Response to a piecewise prescribed oxidizer-flow schedule", pad=10
    )

    ax_port.plot(case_d.time_s, case_d.port_diameter_m * 1e3, color=CASE_COLOURS[2], zorder=4,
                 label=r"port diameter $D_p$")
    ax_port.axhline(GRAIN.outer_diameter_m * 1e3, color=ACCENT, linestyle="--", linewidth=1.3,
                    zorder=2, label=rf"$D_o = {GRAIN.outer_diameter_m * 1e3:.0f}$ mm")
    ax_port.set_ylabel(r"$D_p$   [mm]")
    ax_port.set_ylim(36.0, 95.0)
    ax_port.legend(loc="lower right", ncols=2)

    ax_rate = ax_of
    ax_rate.plot(case_d.time_s, case_d.regression_rate_mm_s, color=CASE_COLOURS[1], zorder=4,
                 label=r"$\dot{r}$   [mm/s]")
    ax_rate.set_ylabel(r"$\dot{r}$   [mm/s]", color=CASE_COLOURS[1])
    ax_rate.tick_params(axis="y", labelcolor=CASE_COLOURS[1])
    ax_rate.set_ylim(-0.06, 1.02)
    ax_rate.set_xlabel(r"time  $t$   [s]")
    ax_rate.set_xlim(0.0, 40.0)

    ax_twin = ax_rate.twinx()
    ax_twin.plot(case_d.time_s, case_d.mixture_ratio, color=NEUTRAL, linestyle="--", zorder=4,
                 label="O/F   [–]")
    ax_twin.set_ylabel(r"mixture ratio  O/F   [–]", color=NEUTRAL)
    ax_twin.tick_params(axis="y", labelcolor=NEUTRAL)
    ax_twin.set_ylim(1.9, 3.25)
    ax_twin.grid(False)

    lines = ax_rate.get_lines()[:1] + ax_twin.get_lines()[:1]
    ax_rate.legend(lines, [line.get_label() for line in lines], loc="lower right", ncols=2)

    schedule_text = ",  ".join(
        f"{duration:g} s at {flow:.3f} kg/s" for duration, flow in PIECEWISE_STEPS
    )
    fig.text(
        0.5,
        0.108,
        f"Illustrative schedule: {schedule_text}.\n"
        "Not a valve, controller or feed-system design, and not tuned to shape O/F. During the "
        r"coast $\dot{r} = 0$ exactly, the port is frozen,"
        "\nand O/F is undefined so no value is drawn.",
        ha="center",
        va="center",
        fontsize=9.0,
    )
    fig.text(0.5, 0.028, DISCLAIMER_MODEL, ha="center", va="center", fontsize=8.5, color="0.25")
    return _save(fig, "fig_m2_e_piecewise.png")


def main() -> int:
    _configure_matplotlib()
    print("Generating Milestone 2 figures (deterministic)...")
    cases = _load_cases()
    builders = (
        figure_a_port_growth,
        figure_b_flux_and_rate,
        figure_c_fuel_flow_and_of,
        figure_d_mass_closure,
        figure_e_piecewise,
    )
    for builder in builders:
        path = builder(cases)
        print(f"  wrote {path.relative_to(FIGURES_DIR.parent)}  ({path.stat().st_size} bytes)")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
