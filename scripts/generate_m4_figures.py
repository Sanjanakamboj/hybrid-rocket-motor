"""Deterministic generation of the Milestone 4 figures.

Produces five PNGs in ``figures/``:

* ``fig_m4_a_tank_blowdown.png``   - tank pressure, temperature and liquid mass
* ``fig_m4_b_feed_coupling.png``   - oxidizer flow and chamber pressure, M3 vs M4
* ``fig_m4_c_thrust_comparison.png`` - the flagship: thrust, M3 vs M4
* ``fig_m4_d_injector_models.png`` - SPI / HEM / Dyer characteristics
* ``fig_m4_e_case_comparison.png`` - parametric case set and closure residuals

Regeneration is byte-identical for a given matplotlib/FreeType build.  PNG glyph
rasterisation depends on the FreeType version, so hashes are not expected to
match across environments with different builds; the content is unchanged.
Verified with matplotlib 3.10.9 / FreeType 2.6.1.

Run with::

    python scripts/generate_m4_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from hybrid_rocket_motor.injector import Injector, InjectorModel
from hybrid_rocket_motor.nitrous_properties import triple_point_pressure_pa
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))

from blowdown_thrust_study import (
    CD,
    REFERENCE_INJECTOR_AREA_M2,
    milestone_3_reference,
    run_case,
)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

PRIMARY = "#0b5394"
SECONDARY = "#1b9e77"
TERTIARY = "#d95f02"
NEUTRAL = "#7570b3"
ACCENT = "#a11b1b"
BAND = "#bcd7ee"

DISCLAIMER = (
    "Generic reduced-order educational model — NOT experimentally validated. Effective injector "
    "area only: no hole geometry.\n"
    "Generic tank volume with no material, wall thickness or structural analysis. "
    r"$c^*$, $\gamma$ and $T_c$ remain prescribed constants."
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
            "svg.hashsalt": "hybrid-rocket-motor-m4",
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
    cases = {
        "A": run_case("CASE A", "nominal")[2],
        "B": run_case("CASE B", "cold", temperature=283.15)[2],
        "C": run_case("CASE C", "warm", temperature=303.15)[2],
        "D": run_case("CASE D", "small area", area=2.5e-6)[2],
        "E": run_case("CASE E", "large area", area=4.5e-6)[2],
        "F": run_case("CASE F", "SPI", model=InjectorModel.SPI)[2],
    }
    return cases, milestone_3_reference()


# ---------------------------------------------------------------------------
# Figure A - tank blowdown
# ---------------------------------------------------------------------------


def figure_a_tank(case) -> Path:
    liquid_left_pct = 100 * case.tank_liquid_mass_kg[-1] / case.tank_liquid_mass_kg[0]
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.4), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.90, top=0.905, bottom=0.180, hspace=0.15)
    ax_p, ax_t, ax_m = axes

    ax_p.plot(case.time_s, case.tank_pressure_pa / 1e5, color=PRIMARY, zorder=4,
              label=r"tank pressure $p_{tank} = p_{sat}(T_{tank})$")
    ax_p.plot(case.time_s, case.chamber_pressure_pa / 1e5, color=TERTIARY, linestyle="--",
              zorder=4, label=r"chamber pressure $p_c$")
    ax_p.fill_between(case.time_s, case.chamber_pressure_pa / 1e5,
                      case.tank_pressure_pa / 1e5, color=BAND, alpha=0.5, zorder=1,
                      label="injector pressure margin")
    ax_p.set_ylabel(r"pressure   [bar]")
    ax_p.set_ylim(0.0, 56.0)
    ax_p.legend(loc="lower left", ncols=3)
    ax_p.set_title(
        "Figure M4-A — Self-pressurising N$_2$O blowdown, reference case", pad=10
    )

    ax_t.plot(case.time_s, case.tank_temperature_k, color=SECONDARY, zorder=4,
              label=r"tank temperature $T_{tank}$")
    ax_t.set_ylabel(r"$T_{tank}$   [K]")
    ax_t.set_ylim(278.0, 295.0)
    ax_t.legend(loc="lower left")

    ax_m.plot(case.time_s, case.tank_liquid_mass_kg, color=PRIMARY, zorder=4,
              label="liquid N$_2$O in tank")
    ax_m.plot(case.time_s, case.tank_vapour_mass_kg, color=NEUTRAL, linestyle="--",
              zorder=4, label="vapour N$_2$O in tank")
    ax_m.plot(case.time_s, case.cumulative_oxidizer_mass_kg, color=TERTIARY,
              linestyle=":", linewidth=2.0, zorder=4, label="oxidizer delivered")
    ax_m.set_xlabel(r"time  $t$   [s]")
    ax_m.set_ylabel(r"mass   [kg]")
    ax_m.set_xlim(0.0, 45.0)
    ax_m.set_ylim(-0.3, 7.0)
    ax_m.legend(loc="upper right", ncols=3)

    fig.text(
        0.5,
        0.095,
        "The tank cools because the latent heat of the vapour it generates comes from the "
        "propellant itself,\n"
        r"so $p_{tank} = p_{sat}(T_{tank})$ falls throughout. "
        f"{liquid_left_pct:.0f} % of the liquid is "
        "still in the tank when the grain burns through.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.028, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m4_a_tank_blowdown.png")


# ---------------------------------------------------------------------------
# Figure B - feed coupling versus prescribed flow
# ---------------------------------------------------------------------------


def figure_b_feed(case, m3) -> Path:
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.4), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.975, top=0.885, bottom=0.200, hspace=0.14)
    ax_flow, ax_pc = axes

    ax_flow.plot(m3.time_s, m3.oxidizer_mass_flow_kg_s, color=NEUTRAL, linestyle="--",
                 zorder=3, label=r"M3: prescribed $\dot{m}_{ox}$ = 0.100 kg/s")
    ax_flow.plot(case.time_s, case.oxidizer_mass_flow_kg_s, color=PRIMARY, linewidth=2.4,
                 zorder=5, label=r"M4: coupled $\dot{m}_{ox}$ from the tank and injector")
    ax_flow.plot(case.time_s, case.fuel_mass_flow_kg_s, color=SECONDARY, linestyle=":",
                 linewidth=2.0, zorder=4, label=r"M4: $\dot{m}_f$ from the regression law")
    ax_flow.set_ylabel(r"mass flow   [kg/s]")
    ax_flow.set_ylim(0.030, 0.115)
    ax_flow.legend(loc="center left")
    ax_flow.set_title(
        "Figure M4-B — What the feed coupling changes: flow and chamber pressure", pad=10
    )

    ax_pc.plot(m3.time_s, m3.chamber_pressure_pa / 1e5, color=NEUTRAL, linestyle="--",
               zorder=3, label="M3: prescribed flow")
    ax_pc.plot(case.time_s, case.chamber_pressure_pa / 1e5, color=PRIMARY, linewidth=2.4,
               zorder=5, label="M4: coupled")
    ax_pc.set_xlabel(r"time  $t$   [s]")
    ax_pc.set_ylabel(r"chamber pressure  $p_c$   [bar]")
    ax_pc.set_xlim(0.0, 45.0)
    ax_pc.set_ylim(23.5, 28.5)
    ax_pc.legend(loc="lower left")

    fig.text(
        0.5,
        0.108,
        r"With the flow prescribed, $\dot{m}_f$ rises as the port opens and $p_c$ rises with it. "
        "With the tank coupled, the\n"
        r"falling feed pressure drags $\dot{m}_{ox}$ down "
        f"{100 * (case.oxidizer_mass_flow_kg_s[-1] / case.oxidizer_mass_flow_kg_s[0] - 1):+.0f} %"
        r" and $p_c$ falls instead of rising.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.030, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m4_b_feed_coupling.png")


# ---------------------------------------------------------------------------
# Figure C - thrust comparison (flagship)
# ---------------------------------------------------------------------------


def figure_c_thrust(case, m3) -> Path:
    fig, ax = plt.subplots(figsize=(9.4, 6.4))
    fig.subplots_adjust(left=0.095, right=0.975, top=0.875, bottom=0.230)

    ax.plot(m3.time_s, m3.thrust_n, color=NEUTRAL, linestyle="--", linewidth=2.2, zorder=4,
            label="Milestone 3: prescribed constant oxidizer flow")
    ax.plot(case.time_s, case.thrust_n, color=PRIMARY, linewidth=2.8, zorder=5,
            label="Milestone 4: coupled N$_2$O tank, injector and motor")

    for history, colour, marker_time in (
        (m3, NEUTRAL, m3.time_s[-1]),
        (case, PRIMARY, case.time_s[-1]),
    ):
        thrust = history.thrust_n[-1]
        ax.plot(marker_time, thrust, "o", color=colour, markersize=8,
                markeredgecolor="white", markeredgewidth=1.0, zorder=6)

    ax.annotate(
        f"M3 ends at {m3.thrust_n[-1]:.0f} N\n(rises +7.6 %)",
        xy=(m3.time_s[-1], m3.thrust_n[-1]),
        xytext=(28.0, 348.0),
        fontsize=9, color=NEUTRAL,
        arrowprops={"arrowstyle": "-", "color": NEUTRAL, "linewidth": 0.9},
    )
    ax.annotate(
        f"M4 ends at {case.thrust_n[-1]:.0f} N\n(falls −8.6 %)",
        xy=(case.time_s[-1], case.thrust_n[-1]),
        xytext=(28.0, 262.0),
        fontsize=9, color=PRIMARY,
        arrowprops={"arrowstyle": "-", "color": PRIMARY, "linewidth": 0.9},
    )
    ax.text(
        0.025,
        0.05,
        f"total impulse   M3 {m3.total_impulse_n_s:.0f} N·s"
        f"   vs   M4 {case.total_impulse_n_s:.0f} N·s"
        f"   ({100 * (case.total_impulse_n_s / m3.total_impulse_n_s - 1):+.1f} %)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.42", "facecolor": "white", "edgecolor": PRIMARY,
              "alpha": 0.95, "linewidth": 0.9},
    )

    ax.set_xlabel(r"time  $t$   [s]")
    ax.set_ylabel(r"thrust  $F$   [N]")
    ax.set_xlim(0.0, 46.0)
    ax.set_ylim(255.0, 365.0)
    ax.legend(loc="upper left")
    ax.set_title(
        "Figure M4-C — Predicted thrust: prescribed flow versus coupled N$_2$O blowdown",
        pad=10,
    )
    fig.text(
        0.5,
        0.125,
        "THE SIGN OF THE THRUST-TIME SLOPE REVERSES. Assuming a constant oxidizer flow does not "
        "merely lose accuracy\nhere — it predicts a rising thrust curve where the coupled model "
        "gives a falling one. Neither curve is validated.",
        ha="center",
        va="center",
        fontsize=9.5,
        color=ACCENT,
        fontweight="bold",
    )
    fig.text(0.5, 0.033, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m4_c_thrust_comparison.png")


# ---------------------------------------------------------------------------
# Figure D - injector characteristics
# ---------------------------------------------------------------------------


def figure_d_injector(case_a, case_f) -> Path:
    state = REFERENCE_TANK.initial_state(
        REFERENCE_TANK_TEMPERATURE_K, REFERENCE_TANK_FILL_FRACTION
    )
    injector = Injector(REFERENCE_INJECTOR_AREA_M2, CD)
    pressures = np.linspace(triple_point_pressure_pa(), 0.995 * state.pressure_pa, 220)
    spi, hem, dyer = [], [], []
    for p2 in pressures:
        result = injector.mass_flow(state, float(p2))
        spi.append(result.spi_mass_flow_kg_s)
        hem.append(result.hem_mass_flow_kg_s)
        dyer.append(result.mass_flow_kg_s)

    fig, ax = plt.subplots(figsize=(9.4, 6.4))
    fig.subplots_adjust(left=0.095, right=0.975, top=0.875, bottom=0.245)

    ax.plot(pressures / 1e5, spi, color=TERTIARY, linestyle="--", zorder=4,
            label=r"SPI  $\dot{m} = C_d A \sqrt{2\rho\,\Delta p}$")
    ax.plot(pressures / 1e5, hem, color=SECONDARY, linestyle=":", linewidth=2.2, zorder=4,
            label=r"HEM  $\dot{m} = C_d A \rho_2 \sqrt{2(h_1 - h_2)}$,  $s_2 = s_1$")
    ax.plot(pressures / 1e5, dyer, color=PRIMARY, linewidth=2.6, zorder=5,
            label=r"Dyer / NHNE  $\frac{\kappa}{1+\kappa}\dot{m}_{SPI} + "
                  r"\frac{1}{1+\kappa}\dot{m}_{HEM}$")

    ax.plot(case_a.chamber_pressure_pa[0] / 1e5, case_a.oxidizer_mass_flow_kg_s[0], "o",
            color=PRIMARY, markersize=9, markeredgecolor="white", markeredgewidth=1.0,
            zorder=7, label="coupled operating point, Dyer")
    ax.plot(case_f.chamber_pressure_pa[0] / 1e5, case_f.oxidizer_mass_flow_kg_s[0], "s",
            color=TERTIARY, markersize=8, markeredgecolor="white", markeredgewidth=1.0,
            zorder=7, label="coupled operating point, SPI")

    peak = int(np.argmax(hem))
    ax.annotate(
        "HEM critical-flow maximum\n(Waxman et al. Eq. 5)",
        xy=(pressures[peak] / 1e5, hem[peak]),
        xytext=(23.0, 0.016),
        fontsize=8.5, color=SECONDARY,
        arrowprops={"arrowstyle": "-", "color": SECONDARY, "linewidth": 0.9},
    )

    ax.set_xlabel(r"downstream (chamber) pressure  $p_2$   [bar]")
    ax.set_ylabel(r"oxidizer mass flow  $\dot{m}_{ox}$   [kg/s]")
    ax.set_xlim(0.0, 52.0)
    ax.set_ylim(0.0, 0.225)
    ax.legend(loc="upper right")
    ax.set_title(
        r"Figure M4-D — Injector characteristics at the initial tank state "
        rf"($p_{{tank}}$ = {state.pressure_pa / 1e5:.1f} bar)",
        pad=10,
    )
    fig.text(
        0.5,
        0.135,
        "Saturated nitrous flashes inside the orifice, so the incompressible SPI equation "
        r"over-predicts flow. For a saturated"
        "\ntank the Dyer parameter is exactly "
        r"$\kappa = 1$, making the blend an equal average of the two limits. "
        "The coupled operating\npoint is where each characteristic meets the chamber's own "
        r"$p_c(\dot{m}_{ox})$ line — the SPI choice shifts it by about 24 %.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.033, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m4_d_injector_models.png")


# ---------------------------------------------------------------------------
# Figure E - case comparison and closure
# ---------------------------------------------------------------------------


def figure_e_cases(cases) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.8))
    fig.subplots_adjust(left=0.075, right=0.975, top=0.845, bottom=0.290, wspace=0.24)
    ax_thrust, ax_closure = axes

    styles = (
        ("A", "nominal, 293 K, 3.5 mm$^2$", PRIMARY, "-"),
        ("B", "cold tank, 283 K", SECONDARY, "--"),
        ("C", "warm tank, 303 K", TERTIARY, "--"),
        ("D", "small area, 2.5 mm$^2$", NEUTRAL, ":"),
        ("E", "large area, 4.5 mm$^2$", ACCENT, ":"),
    )
    for key, label, colour, style in styles:
        case = cases[key]
        ax_thrust.plot(case.time_s, case.thrust_n, color=colour, linestyle=style, zorder=4,
                       label=f"Case {key}: {label}")
    ax_thrust.set_xlabel(r"time  $t$   [s]")
    ax_thrust.set_ylabel(r"thrust  $F$   [N]")
    ax_thrust.set_xlim(0.0, 50.0)
    ax_thrust.set_ylim(215.0, 432.0)
    ax_thrust.legend(loc="upper right")
    ax_thrust.set_title("Parametric case set", pad=8)

    labels, oxidizer, fuel, energy = [], [], [], []
    for key in ("A", "B", "C", "D", "E", "F"):
        case = cases[key]
        labels.append(key)
        oxidizer.append(max(float(np.max(np.abs(case.oxidizer_mass_closure_residual_kg))), 1e-18))
        fuel.append(max(float(np.max(np.abs(case.fuel_mass_closure_residual_kg))), 1e-18))
        energy.append(max(case.relative_energy_closure_residual, 1e-18))

    positions = np.arange(len(labels), dtype=float)
    width = 0.27
    ax_closure.bar(positions - width, oxidizer, width, color=PRIMARY, zorder=3,
                   label="oxidizer mass residual [kg]")
    ax_closure.bar(positions, fuel, width, color=SECONDARY, zorder=3,
                   label="fuel mass residual [kg]")
    ax_closure.bar(positions + width, energy, width, color=TERTIARY, zorder=3,
                   label="tank energy residual [relative]")
    ax_closure.set_yscale("log")
    ax_closure.set_xticks(positions)
    ax_closure.set_xticklabels([f"Case {label}" for label in labels])
    ax_closure.set_ylabel("worst closure residual")
    ax_closure.set_ylim(1e-17, 1e-4)
    ax_closure.legend(loc="upper left")
    ax_closure.set_title("Independent conservation closure", pad=8)

    fig.suptitle(
        "Figure M4-E — Parametric cases and independent mass / energy closure",
        fontsize=12.5, y=0.955,
    )
    fig.text(
        0.5,
        0.145,
        "Left: initial tank temperature sets how hard the motor starts; effective injector area "
        "trades thrust against duration.\n"
        "Right: oxidizer drawn from the tank against oxidizer integrated at the injector; fuel "
        "from port geometry against integrated\n"
        r"$\dot{m}_f$; and the adiabatic tank energy balance. All are reconstructed independently "
        "of the solver's own accumulators.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.033, DISCLAIMER, ha="center", va="center", fontsize=8, color="0.25")
    return _save(fig, "fig_m4_e_case_comparison.png")


def main() -> int:
    _configure_matplotlib()
    print("Generating Milestone 4 figures (deterministic)...")
    cases, m3 = _load()
    outputs = [
        figure_a_tank(cases["A"]),
        figure_b_feed(cases["A"], m3),
        figure_c_thrust(cases["A"], m3),
        figure_d_injector(cases["A"], cases["F"]),
        figure_e_cases(cases),
    ]
    for path in outputs:
        print(f"  wrote {path.relative_to(FIGURES_DIR.parent)}  ({path.stat().st_size} bytes)")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
