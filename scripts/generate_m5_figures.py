"""Deterministic generation of the Milestone 5 figures.

Produces five PNGs in ``figures/``:

* ``fig_m5_a_thermochemistry_table.png`` - the CEA table itself, with the M3 constants
* ``fig_m5_b_thrust_comparison.png``     - the flagship: coupled thrust, M4 vs M5
* ``fig_m5_c_property_histories.png``    - how the chemistry moves during a burn
* ``fig_m5_d_case_comparison.png``       - the case set and its closure residuals
* ``fig_m5_e_sensitivities.png``         - efficiency, grid resolution, frozen vs equilibrium

Regeneration is byte-identical for a given matplotlib/FreeType build.  PNG glyph
rasterisation depends on the FreeType version, so hashes are not expected to
match across environments with different builds; the content is unchanged.
Verified with matplotlib 3.10.9 / FreeType 2.6.1.

Run with::

    python scripts/generate_m5_figures.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION

sys.path.insert(0, str(Path(__file__).resolve().parent))

from variable_thermochemistry_study import (
    NOZZLE_CHEMISTRY_PATH,
    TABLE,
    run_pair,
)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

PRIMARY = "#0b5394"
SECONDARY = "#1b9e77"
TERTIARY = "#d95f02"
NEUTRAL = "#7570b3"
ACCENT = "#a11b1b"
BAND = "#bcd7ee"

#: Reference pressure for the table slices [Pa].  Mid-envelope for this motor.
SLICE_PRESSURE_PA = 25.0e5

DISCLAIMER = (
    "Generic reduced-order educational model — NOT experimentally validated. Equilibrium "
    "chemistry only: no finite-rate combustion,\n"
    "no ignition transient, no combustion instability. A single "
    r"$\eta_{c^*}$ stands in for every combustion loss."
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
            "svg.hashsalt": "hybrid-rocket-motor-m5",
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


def _table_slice(field: str) -> np.ndarray:
    return np.array(
        [
            getattr(TABLE.evaluate(float(of), SLICE_PRESSURE_PA, c_star_efficiency=1.0), field)
            for of in TABLE.of_values
        ]
    )


# ---------------------------------------------------------------------------
# Figure A - the table
# ---------------------------------------------------------------------------


def figure_a_table(operating_band: tuple[float, float]) -> Path:
    of_values = TABLE.of_values
    ideal_c_star = _table_slice("ideal_c_star_m_s")
    gamma = _table_slice("gamma")
    temperature = _table_slice("chamber_temperature_k")
    soot = _table_slice("condensed_phase_mass_ratio")
    # The boundary is where equilibrium stops producing condensed carbon.  It is
    # located from the largest step in gamma rather than from a tolerance on the
    # soot ratio: the step always exists, and its position is exactly what the
    # bilinear interpolation has to cope with.
    step = int(np.argmax(np.abs(np.diff(gamma))))
    boundary = float(0.5 * (of_values[step] + of_values[step + 1]))

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.6))
    fig.subplots_adjust(left=0.085, right=0.975, top=0.915, bottom=0.255, hspace=0.30, wspace=0.24)
    fig.suptitle(
        "Figure M5-A — NASA CEA equilibrium properties vs mixture ratio "
        f"(N₂O/HTPB, $p_c$ = {SLICE_PRESSURE_PA / 1e5:.0f} bar)",
        fontsize=12,
    )

    panels = (
        (
            axes[0][0],
            ideal_c_star,
            r"ideal $c^*$   [m/s]",
            PRIMARY,
            ILLUSTRATIVE_N2O_HTPB_COMBUSTION.ideal_c_star_m_s,
        ),
        (axes[0][1], gamma, r"$\gamma$   [-]", SECONDARY, ILLUSTRATIVE_N2O_HTPB_COMBUSTION.gamma),
        (
            axes[1][0],
            temperature,
            r"$T_c$   [K]",
            TERTIARY,
            ILLUSTRATIVE_N2O_HTPB_COMBUSTION.chamber_temperature_k,
        ),
        (axes[1][1], soot, r"$M_{eff} / MW_{gas}$   [-]", NEUTRAL, None),
    )
    for axis, values, label, colour, prescribed in panels:
        axis.axvspan(
            *operating_band, color=BAND, alpha=0.55, zorder=1, label="operating band, case A"
        )
        axis.axvline(
            boundary,
            color=ACCENT,
            linestyle=":",
            linewidth=1.4,
            zorder=3,
            label=f"condensed-carbon boundary, O/F ≈ {boundary:.2f}",
        )
        axis.plot(of_values, values, color=colour, zorder=4, label="CEA equilibrium")
        if prescribed is not None:
            axis.axhline(
                prescribed,
                color="0.35",
                linestyle="--",
                linewidth=1.3,
                zorder=3,
                label="M3 prescribed constant",
            )
        axis.set_xlabel("mixture ratio  O/F   [-]")
        axis.set_ylabel(label)
        axis.set_xlim(float(of_values[0]), float(of_values[-1]))

    axes[1][1].set_ylim(0.98, 1.48)
    axes[0][0].legend(loc="lower right")

    fig.text(
        0.085,
        0.165,
        "Where equilibrium stops producing condensed carbon, γ jumps by about 4.5 % in a single "
        "grid step while $c^*$ and $T_c$ merely kink.\n"
        "Interpolation is bilinear, not spline, precisely because a spline would overshoot at "
        "that step. This motor never reaches the boundary: it runs fuel-rich of it throughout,\n"
        "so the best point SAMPLED here is the upper table edge and is not an optimum.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.085, 0.052, DISCLAIMER, fontsize=7.6, color="0.4", va="top")
    return _save(fig, "fig_m5_a_thermochemistry_table.png")


# ---------------------------------------------------------------------------
# Figure B - the flagship thrust comparison
# ---------------------------------------------------------------------------


def figure_b_thrust(frozen, variable) -> Path:
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.4), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.965, top=0.910, bottom=0.265, hspace=0.15)
    fig.suptitle(
        "Figure M5-B — coupled thrust: prescribed properties (M4) vs O/F-dependent chemistry (M5)",
        fontsize=12,
    )
    ax_f, ax_i = axes

    ax_f.plot(
        frozen.time_s,
        frozen.thrust_n,
        color="0.35",
        linestyle="--",
        zorder=4,
        label=r"M4, frozen $c^*$, $\gamma$, $T_c$",
    )
    ax_f.plot(
        variable.time_s,
        variable.thrust_n,
        color=PRIMARY,
        zorder=5,
        label="M5, CEA table at the instantaneous O/F",
    )
    ax_f.fill_between(
        variable.time_s,
        variable.thrust_n,
        np.interp(variable.time_s, frozen.time_s, frozen.thrust_n),
        color=BAND,
        alpha=0.55,
        zorder=1,
        label="difference",
    )
    ax_f.set_ylabel("thrust  $F$   [N]")
    ax_f.set_ylim(0.0, 360.0)
    ax_f.legend(loc="lower left")

    frozen_decay = 100.0 * (frozen.thrust_n[-1] / frozen.thrust_n[0] - 1.0)
    variable_decay = 100.0 * (variable.thrust_n[-1] / variable.thrust_n[0] - 1.0)
    ax_f.annotate(
        f"M4 decay {frozen_decay:+.1f} %\nM5 decay {variable_decay:+.1f} %\n"
        f"same sign, {variable_decay / frozen_decay:.2f}× the magnitude",
        xy=(0.985, 0.62),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=8.6,
        bbox={"boxstyle": "round,pad=0.42", "facecolor": "white", "edgecolor": "0.7"},
    )

    ax_i.plot(
        frozen.time_s,
        frozen.cumulative_impulse_n_s,
        color="0.35",
        linestyle="--",
        zorder=4,
        label=f"M4, total {frozen.total_impulse_n_s:.0f} N s",
    )
    ax_i.plot(
        variable.time_s,
        variable.cumulative_impulse_n_s,
        color=SECONDARY,
        zorder=5,
        label=f"M5, total {variable.total_impulse_n_s:.0f} N s",
    )
    ax_i.set_xlabel("time  $t$   [s]")
    ax_i.set_ylabel("cumulative impulse   [N s]")
    ax_i.set_xlim(0.0, max(frozen.time_s[-1], variable.time_s[-1]))
    ax_i.legend(loc="upper left")

    fig.text(
        0.105,
        0.175,
        "Both runs use identical tank, injector, grain and nozzle, and both terminate on grain "
        "burnout. The only difference is where the\n"
        "combustion properties come from. The curve moves DOWN (the equilibrium $c^*$ here is "
        "below the illustrative constant) and it also\n"
        "steepens (a falling O/F now drags $c^*$ down with it). Level and slope are separate "
        "effects.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.105, 0.055, DISCLAIMER, fontsize=7.6, color="0.4", va="top")
    return _save(fig, "fig_m5_b_thrust_comparison.png")


# ---------------------------------------------------------------------------
# Figure C - property histories
# ---------------------------------------------------------------------------


def figure_c_histories(frozen, variable) -> Path:
    firing = variable.firing
    time = variable.time_s[firing]

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.6), sharex=True)
    fig.subplots_adjust(left=0.085, right=0.975, top=0.915, bottom=0.265, hspace=0.20, wspace=0.24)
    fig.suptitle("Figure M5-C — what the chemistry does during a burn (case A)", fontsize=12)

    ax_of, ax_c, ax_g, ax_t = axes[0][0], axes[0][1], axes[1][0], axes[1][1]

    ax_of.plot(
        frozen.time_s, frozen.mixture_ratio, color="0.35", linestyle="--", zorder=4, label="M4"
    )
    ax_of.plot(time, variable.mixture_ratio[firing], color=PRIMARY, zorder=5, label="M5")
    ax_of.set_ylabel("mixture ratio  O/F   [-]")
    ax_of.legend(loc="upper right")

    ax_c.axhline(
        ILLUSTRATIVE_N2O_HTPB_COMBUSTION.c_star_m_s,
        color="0.35",
        linestyle="--",
        zorder=4,
        label="M4, constant by assumption",
    )
    ax_c.plot(time, variable.c_star_m_s[firing], color=SECONDARY, zorder=5, label="M5, follows O/F")
    ax_c.set_ylabel(r"delivered $c^*$   [m/s]")
    ax_c.legend(loc="center right")

    ax_g.axhline(
        ILLUSTRATIVE_N2O_HTPB_COMBUSTION.gamma, color="0.35", linestyle="--", zorder=4, label="M4"
    )
    ax_g.plot(time, variable.gamma[firing], color=TERTIARY, zorder=5, label="M5")
    ax_g.set_ylabel(r"$\gamma$   [-]")
    ax_g.set_xlabel("time  $t$   [s]")
    ax_g.legend(loc="center right")

    ax_t.axhline(
        ILLUSTRATIVE_N2O_HTPB_COMBUSTION.chamber_temperature_k,
        color="0.35",
        linestyle="--",
        zorder=4,
        label="M4",
    )
    ax_t.plot(time, variable.chamber_temperature_k[firing], color=NEUTRAL, zorder=5, label="M5")
    ax_t.set_ylabel(r"$T_c$   [K]")
    ax_t.set_xlabel("time  $t$   [s]")
    ax_t.legend(loc="upper right")

    for axis in (ax_of, ax_c, ax_g, ax_t):
        axis.set_xlim(0.0, float(time[-1]))

    fig.text(
        0.085,
        0.175,
        "The dashed lines are Milestone 4: flat by assumption, not by physics. As the port opens "
        "the mixture ratio falls, and on the\n"
        "fuel-rich side of the table a falling O/F means a falling $c^*$ and a falling $T_c$. "
        r"$\gamma$ barely moves over this range, so the extra thrust"
        "\n"
        "decay is almost entirely a $c^*$ effect acting through chamber pressure.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.085, 0.055, DISCLAIMER, fontsize=7.6, color="0.4", va="top")
    return _save(fig, "fig_m5_c_property_histories.png")


# ---------------------------------------------------------------------------
# Figure D - case comparison and closure
# ---------------------------------------------------------------------------


def figure_d_cases(cases) -> Path:
    labels = [label for label, _f, _v in cases]
    frozen_impulse = np.array([f.total_impulse_n_s for _l, f, _v in cases])
    variable_impulse = np.array([v.total_impulse_n_s for _l, _f, v in cases])
    positions = np.arange(len(labels))

    fig, axes = plt.subplots(1, 3, figsize=(12.4, 5.3))
    fig.subplots_adjust(left=0.065, right=0.985, top=0.885, bottom=0.320, wspace=0.30)
    fig.suptitle(
        "Figure M5-D — the case set: how much the chemistry is worth, and how well it closes",
        fontsize=12,
    )
    ax_bars, ax_change, ax_residual = axes

    ax_bars.bar(
        positions - 0.19, frozen_impulse, width=0.36, color="0.55", label="M4 frozen", zorder=3
    )
    ax_bars.bar(
        positions + 0.19, variable_impulse, width=0.36, color=PRIMARY, label="M5 variable", zorder=3
    )
    ax_bars.set_xticks(positions, labels)
    ax_bars.set_ylabel("total impulse   [N s]")
    ax_bars.set_ylim(0.0, 1.22 * float(frozen_impulse.max()))
    ax_bars.set_title("total impulse", fontsize=10)
    ax_bars.legend(loc="upper right")

    change = 100.0 * (variable_impulse / frozen_impulse - 1.0)
    ax_change.bar(positions, change, width=0.55, color=TERTIARY, zorder=3)
    ax_change.axhline(0.0, color="0.3", linewidth=1.0, zorder=4)
    ax_change.set_xticks(positions, labels)
    ax_change.set_ylabel("change from M4   [%]")
    ax_change.set_title("every case loses impulse", fontsize=10)
    for position, value in zip(positions, change, strict=True):
        ax_change.annotate(
            f"{value:.1f}",
            xy=(position, value),
            ha="center",
            va="top",
            xytext=(0, -4),
            textcoords="offset points",
            fontsize=8.0,
        )

    residuals = {
        r"feed  [kg/s]": [np.max(np.abs(v.feed_residual_kg_s)) for _l, _f, v in cases],
        r"chamber  [Pa]": [np.max(np.abs(v.c_star_identity_residual_pa)) for _l, _f, v in cases],
        r"ox mass  [kg]": [
            np.max(np.abs(v.oxidizer_mass_closure_residual_kg)) for _l, _f, v in cases
        ],
    }
    markers = ("o", "s", "^")
    colours = (PRIMARY, SECONDARY, ACCENT)
    for (name, values), marker, colour in zip(residuals.items(), markers, colours, strict=True):
        ax_residual.semilogy(
            positions,
            np.maximum(values, 1e-18),
            marker=marker,
            color=colour,
            linestyle="none",
            markersize=7,
            label=name,
            zorder=4,
        )
    ax_residual.set_xticks(positions, labels)
    ax_residual.set_ylabel("worst residual over the run")
    ax_residual.set_ylim(1e-16, 1e-4)
    ax_residual.set_title("closure residuals (log scale)", fontsize=10)
    ax_residual.legend(loc="upper right")

    fig.text(
        0.065,
        0.215,
        "Cases: A nominal · B tank at 283.15 K · C tank at 303.15 K · D injector 2.5 mm² · "
        "E injector 4.5 mm² · G 2.5 L tank (oxidizer-limited).\n"
        "Every case is run twice on identical conditions, so the bars differ only by the "
        "combustion properties. The residuals are re-formed from\n"
        "the reported histories, not read back from the solver, and every case stayed inside "
        "the tabulated O/F and pressure rectangle.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.065, 0.070, DISCLAIMER, fontsize=7.6, color="0.4", va="top")
    return _save(fig, "fig_m5_d_case_comparison.png")


# ---------------------------------------------------------------------------
# Figure E - sensitivities
# ---------------------------------------------------------------------------


def figure_e_sensitivities(efficiency_runs, grid_runs, reference) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 5.3))
    fig.subplots_adjust(left=0.065, right=0.985, top=0.885, bottom=0.345, wspace=0.32)
    fig.suptitle("Figure M5-E — what the remaining assumptions are worth", fontsize=12)
    ax_eta, ax_grid, ax_frozen = axes

    etas = np.array([eta for eta, _r in efficiency_runs])
    impulses = np.array([r.total_impulse_n_s for _e, r in efficiency_runs])
    pressures = np.array([r.chamber_pressure_pa[0] / 1e5 for _e, r in efficiency_runs])
    ax_eta.plot(
        etas,
        100.0 * (impulses / impulses[0] - 1.0),
        marker="o",
        color=PRIMARY,
        zorder=4,
        label="total impulse",
    )
    ax_eta.plot(
        etas,
        100.0 * (pressures / pressures[0] - 1.0),
        marker="s",
        color=TERTIARY,
        zorder=4,
        label=r"initial $p_c$",
    )
    ax_eta.set_xlabel(r"$\eta_{c^*}$   [-]")
    ax_eta.set_ylabel(r"change from $\eta_{c^*}$ = 0.92   [%]")
    ax_eta.set_title(r"$c^*$ efficiency", fontsize=10)
    ax_eta.legend(loc="upper left")

    nodes = np.array([n for n, _r in grid_runs])
    grid_impulse = np.array([r.total_impulse_n_s for _n, r in grid_runs])
    ax_grid.loglog(
        nodes,
        np.abs(grid_impulse / reference.total_impulse_n_s - 1.0),
        marker="o",
        color=SECONDARY,
        zorder=4,
    )
    ax_grid.set_xlabel("table nodes   [-]")
    ax_grid.set_ylabel("|relative impulse error|   [-]")
    ax_grid.set_title("interpolation grid resolution", fontsize=10)

    with NOZZLE_CHEMISTRY_PATH.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    of_values = np.array([float(r["of"]) for r in rows])
    equilibrium = np.array([float(r["vacuum_isp_equilibrium_s"]) for r in rows])
    converged = np.array([bool(int(r["frozen_converged"])) for r in rows])
    frozen_isp = np.array([float(r["vacuum_isp_frozen_s"]) for r in rows])
    delta = np.full(of_values.shape, np.nan)
    delta[converged] = 100.0 * (frozen_isp[converged] / equilibrium[converged] - 1.0)
    gap = of_values[~converged]
    ax_frozen.axvspan(
        gap.min(), gap.max(), color="0.88", zorder=1, label="CEA frozen did not converge"
    )
    ax_frozen.plot(
        of_values, delta, marker="o", color=NEUTRAL, zorder=4, label="frozen vs equilibrium"
    )
    ax_frozen.axhline(0.0, color="0.3", linewidth=1.0, zorder=3)
    ax_frozen.set_xlabel("mixture ratio  O/F   [-]")
    ax_frozen.set_ylabel(r"vacuum $I_{sp}$ change   [%]")
    ax_frozen.set_title("expansion chemistry", fontsize=10)
    ax_frozen.legend(loc="lower right")

    fig.text(
        0.065,
        0.240,
        r"Left: raising $\eta_{c^*}$ lifts chamber pressure but barely moves total impulse — a "
        "higher $c^*$ raises $p_c$, which throttles the\n"
        "pressure-fed injector, and the burn is grain-limited so the fuel mass is fixed. "
        "Middle: the committed table is far finer than needed.\n"
        "Right: the single-γ expansion sits between CEA's two limits; the shaded band is where "
        "CEA could not run frozen at all, because of\n"
        "condensed carbon — and that band straddles this motor's operating range.",
        fontsize=8.2,
        color="0.25",
        va="top",
    )
    fig.text(0.065, 0.070, DISCLAIMER, fontsize=7.6, color="0.4", va="top")
    return _save(fig, "fig_m5_e_sensitivities.png")


def main() -> int:
    _configure_matplotlib()
    print("Generating Milestone 5 figures...")

    frozen_a, variable_a = run_pair()
    band = (float(np.nanmin(variable_a.mixture_ratio)), float(np.nanmax(variable_a.mixture_ratio)))

    paths = [
        figure_a_table(band),
        figure_b_thrust(frozen_a, variable_a),
        figure_c_histories(frozen_a, variable_a),
    ]

    specifications = [
        ("A", {}),
        ("B", {"temperature": 283.15}),
        ("C", {"temperature": 303.15}),
        ("D", {"area": 2.5e-6}),
        ("E", {"area": 4.5e-6}),
        ("G", {"volume_l": 2.5}),
    ]
    cases = []
    for label, kwargs in specifications:
        if label == "A":
            cases.append((label, frozen_a, variable_a))
            continue
        frozen, variable = run_pair(**kwargs)
        cases.append((label, frozen, variable))
    paths.append(figure_d_cases(cases))

    efficiency_runs = [(0.92, None), (0.96, None), (1.00, None)]
    efficiency_runs = [
        (eta, variable_a if eta == 0.96 else run_pair(c_star_efficiency=eta)[1])
        for eta, _ in efficiency_runs
    ]
    grid_runs = []
    for of_stride, pressure_stride in ((2, 1), (4, 2), (8, 4), (16, 7)):
        coarse = TABLE.subsample(of_stride=of_stride, pressure_stride=pressure_stride)
        _, result = run_pair(table=coarse)
        grid_runs.append((coarse.of_values.size * coarse.pressure_values.size, result))
    paths.append(figure_e_sensitivities(efficiency_runs, grid_runs, variable_a))

    for path in paths:
        print(f"  wrote {path.relative_to(FIGURES_DIR.parent)} ({path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
