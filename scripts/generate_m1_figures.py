"""Deterministic generation of the Milestone 1 figures.

Produces three PNGs in ``figures/``:

* ``fig_a_grain_geometry.png``      - generic grain schematic
* ``fig_b_regression_vs_flux.png``  - regression rate versus oxidizer mass flux
* ``fig_c_port_sensitivity.png``    - port-size sensitivity at fixed oxidizer flow

Regeneration is byte-identical *for a given matplotlib/FreeType build*: the Agg
backend is forced, all rcParams that affect rendering are pinned, no timestamp
metadata is written, and every input is a fixed constant.  PNG glyph
rasterisation does depend on the FreeType version, so hashes are not expected to
match across environments with different matplotlib or FreeType builds; the
figure content is unchanged.  Verified with matplotlib 3.10.9 / FreeType 2.6.1.

Run with::

    python scripts/generate_m1_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle

from hybrid_rocket_motor import (
    ILLUSTRATIVE_ANCHOR_FLUX_SI,
    ILLUSTRATIVE_SENSITIVITY_LAWS,
    REPRESENTATIVE_GRAIN,
    REZAEI_2018_N2O_HTPB,
    OperatingPoint,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))

from regression_study import (
    OXIDIZER_MASS_FLOWS_KG_S,
    PORT_DIAMETERS_M,
    REPRESENTATIVE_M_DOT_OX_KG_S,
)

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

GRAIN = REPRESENTATIVE_GRAIN
LAW = REZAEI_2018_N2O_HTPB

# Colour-blind-safe, print-legible palette.
FUEL_COLOUR = "#c8b18f"
PORT_COLOUR = "#f7f7f7"
SOURCED_COLOUR = "#0b5394"
ILLUSTRATIVE_COLOURS = ("#d95f02", "#7570b3", "#1b9e77")
BAND_COLOUR = "#bcd7ee"
ACCENT = "#a11b1b"

DISCLAIMER_GEOMETRY = "Generic illustrative geometry — not a manufacturing drawing."
DISCLAIMER_MODEL = (
    "Generic reduced-order educational model — no validation against any real motor.\n"
    "Milestone 1: no chamber-pressure, nozzle or thrust result is computed."
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
            "svg.hashsalt": "hybrid-rocket-motor-m1",
            "path.simplify": False,
            "axes.axisbelow": True,
        }
    )


def _save(fig, name: str) -> Path:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / name
    # metadata={"Software": None} suppresses the only non-content PNG tag
    # matplotlib writes, keeping regeneration byte-identical.
    fig.savefig(path, format="png", metadata={"Software": None})
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Figure A - grain geometry schematic
# ---------------------------------------------------------------------------


def figure_a_geometry() -> Path:
    length_mm = GRAIN.length_m * 1e3
    d_outer_mm = GRAIN.outer_diameter_m * 1e3
    d_port_mm = GRAIN.initial_port_diameter_m * 1e3
    web_mm = GRAIN.initial_web_thickness_m * 1e3
    r_outer, r_port = d_outer_mm / 2.0, d_port_mm / 2.0

    fig = plt.figure(figsize=(10.6, 5.1))
    gs = fig.add_gridspec(
        1, 2, width_ratios=[2.5, 1.0], left=0.045, right=0.985, top=0.795, bottom=0.205, wspace=0.16
    )
    ax = fig.add_subplot(gs[0, 0])
    ax_end = fig.add_subplot(gs[0, 1])

    # -- longitudinal section ------------------------------------------------
    # Vertical scale is exaggerated (aspect left as "auto") so that a 400 mm x
    # 90 mm grain does not collapse into an unreadable sliver.  Every dimension
    # is annotated numerically, so no length is read off the drawing.
    for y_bottom in (r_port, -r_outer):
        ax.add_patch(
            Rectangle(
                (0.0, y_bottom),
                length_mm,
                r_outer - r_port,
                facecolor=FUEL_COLOUR,
                edgecolor="black",
                linewidth=1.3,
                hatch="///",
                zorder=2,
            )
        )
    ax.add_patch(
        Rectangle(
            (0.0, -r_port), length_mm, d_port_mm, facecolor=PORT_COLOUR, edgecolor="none", zorder=1
        )
    )

    for sign in (1, -1):
        ax.plot(
            [0.0, length_mm],
            [sign * r_port, sign * r_port],
            color=ACCENT,
            linewidth=2.8,
            zorder=4,
            solid_capstyle="butt",
        )
        for x in np.linspace(0.10 * length_mm, 0.90 * length_mm, 5):
            ax.annotate(
                "",
                xy=(x, sign * (r_port + 9.0)),
                xytext=(x, sign * r_port),
                arrowprops={"arrowstyle": "-|>", "color": ACCENT, "linewidth": 1.2},
                zorder=5,
            )
    ax.text(
        0.5 * length_mm,
        r_outer + 5.0,
        r"port wall regresses outward at $\dot{r}$",
        ha="center",
        va="bottom",
        color=ACCENT,
        fontsize=9,
    )

    ax.annotate(
        "",
        xy=(0.72 * length_mm, 0.0),
        xytext=(0.05 * length_mm, 0.0),
        arrowprops={"arrowstyle": "-|>", "color": SOURCED_COLOUR, "linewidth": 2.2},
        zorder=6,
    )
    ax.text(
        0.385 * length_mm,
        4.5,
        r"prescribed $\dot{m}_{ox}$ through the port",
        ha="center",
        va="bottom",
        color=SOURCED_COLOUR,
        fontsize=9.5,
    )

    # Grain length dimension, below the section.
    y_dim = -r_outer - 17.0
    ax.annotate(
        "",
        xy=(length_mm, y_dim),
        xytext=(0.0, y_dim),
        arrowprops={"arrowstyle": "<|-|>", "color": "black", "linewidth": 1.1},
    )
    ax.text(
        length_mm / 2.0,
        y_dim - 4.0,
        rf"grain length  $L = {length_mm:.0f}$ mm",
        ha="center",
        va="top",
    )
    for x in (0.0, length_mm):
        ax.plot([x, x], [-r_outer, y_dim - 2.0], color="black", linewidth=0.7, linestyle=":")

    # Outer diameter, dimensioned clear of the right-hand end face.
    x_do = length_mm + 34.0
    ax.annotate(
        "",
        xy=(x_do, r_outer),
        xytext=(x_do, -r_outer),
        arrowprops={"arrowstyle": "<|-|>", "color": "black", "linewidth": 1.1},
    )
    ax.text(x_do + 9.0, 0.0, rf"$D_o = {d_outer_mm:.0f}$ mm", ha="left", va="center", rotation=90)
    for y in (r_outer, -r_outer):
        ax.plot([length_mm, x_do + 3.0], [y, y], color="black", linewidth=0.7, linestyle=":")

    # Port diameter, dimensioned on the opposite (left) side so the two
    # dimensions cannot collide with each other or with the grain body.
    x_dp = -20.0
    ax.annotate(
        "",
        xy=(x_dp, r_port),
        xytext=(x_dp, -r_port),
        arrowprops={"arrowstyle": "<|-|>", "color": ACCENT, "linewidth": 1.1},
    )
    ax.text(
        x_dp - 9.0,
        0.0,
        rf"$D_{{p,0}} = {d_port_mm:.0f}$ mm",
        ha="right",
        va="center",
        rotation=90,
        color=ACCENT,
    )
    for y in (r_port, -r_port):
        ax.plot([x_dp - 3.0, 0.0], [y, y], color=ACCENT, linewidth=0.7, linestyle=":")

    ax.set_xlim(-78.0, length_mm + 84.0)
    ax.set_ylim(-r_outer - 40.0, r_outer + 34.0)
    ax.set_title("Longitudinal section  (vertical scale exaggerated)", pad=9)
    ax.axis("off")

    ax.plot([], [], color=ACCENT, linewidth=2.8, label="burning surface (regressing port wall)")
    ax.add_patch(
        Rectangle((0, 0), 0, 0, facecolor=FUEL_COLOUR, edgecolor="black", hatch="///",
                  label="solid fuel")
    )
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.035), ncols=2, frameon=False)

    # -- end view ------------------------------------------------------------
    ax_end.add_patch(
        Circle((0.0, 0.0), r_outer, facecolor=FUEL_COLOUR, edgecolor="black", linewidth=1.3,
               hatch="///", zorder=2)
    )
    ax_end.add_patch(
        Circle((0.0, 0.0), r_port, facecolor=PORT_COLOUR, edgecolor=ACCENT, linewidth=2.8, zorder=3)
    )
    ax_end.annotate(
        "",
        xy=(r_port, 0.0),
        xytext=(-r_port, 0.0),
        arrowprops={"arrowstyle": "<|-|>", "color": ACCENT, "linewidth": 1.1},
        zorder=4,
    )
    ax_end.text(
        0.0,
        r_outer + 5.0,
        rf"single circular port" "\n" rf"$D_{{p,0}} = {d_port_mm:.0f}$ mm",
        ha="center",
        va="bottom",
        color=ACCENT,
        fontsize=9.5,
    )

    y_do = -r_outer - 13.0
    ax_end.annotate(
        "",
        xy=(r_outer, y_do),
        xytext=(-r_outer, y_do),
        arrowprops={"arrowstyle": "<|-|>", "color": "black", "linewidth": 1.1},
    )
    ax_end.text(
        0.0, y_do - 4.0, rf"$D_o = {d_outer_mm:.0f}$ mm", ha="center", va="top", fontsize=9.5
    )
    for x in (-r_outer, r_outer):
        ax_end.plot([x, x], [0.0, y_do - 2.0], color="black", linewidth=0.7, linestyle=":")

    ax_end.set_xlim(-r_outer - 16.0, r_outer + 16.0)
    ax_end.set_ylim(-r_outer - 40.0, r_outer + 34.0)
    ax_end.set_aspect("equal")
    ax_end.set_title("End view", pad=9)
    ax_end.axis("off")

    # -- titles and annotation ----------------------------------------------
    fig.suptitle(
        "Figure A — Generic single-port cylindrical hybrid fuel grain (HTPB-like)",
        fontsize=13,
        y=0.955,
    )
    fig.text(
        0.5,
        0.885,
        rf"$\rho_f = {GRAIN.fuel_density_kg_m3:.0f}$ kg/m$^3$     "
        rf"initial web $=(D_o-D_{{p,0}})/2 = {web_mm:.0f}$ mm     "
        rf"$A_{{port}} = \pi D_p^2/4$     $A_{{burn}} = \pi D_p L$",
        ha="center",
        va="center",
        fontsize=10,
    )
    fig.text(
        0.5,
        0.105,
        DISCLAIMER_GEOMETRY,
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="bold",
        color=ACCENT,
    )
    fig.text(
        0.5,
        0.042,
        "Schematic only: no wall thickness, liner, insulation, injector, igniter or nozzle "
        "is represented or sized.",
        ha="center",
        va="center",
        fontsize=8.5,
        color="0.25",
    )

    return _save(fig, "fig_a_grain_geometry.png")


# ---------------------------------------------------------------------------
# Figure B - regression rate versus oxidizer mass flux
# ---------------------------------------------------------------------------


def figure_b_regression() -> Path:
    flux_lo, flux_hi = 20.0, 140.0
    flux = np.linspace(flux_lo, flux_hi, 601)

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    fig.subplots_adjust(left=0.085, right=0.975, top=0.845, bottom=0.215)

    band_lo, band_hi = LAW.valid_flux_range_si
    ax.axvspan(
        band_lo,
        band_hi,
        color=BAND_COLOUR,
        alpha=0.55,
        zorder=0,
        label=rf"source data range, {band_lo:.0f}–{band_hi:.0f} kg/(m$^2\,$s)",
    )

    ax.plot(
        flux,
        LAW.regression_rate_mm_s(flux),
        color=SOURCED_COLOUR,
        linewidth=2.6,
        zorder=4,
        label=(
            "SOURCED: Rezaei et al. (2018) N$_2$O/HTPB\n"
            r"$\dot{r}\,[\mathrm{mm/s}] = 0.3977\,(G_{ox}\,[\mathrm{g/(cm^2 s)}])^{0.3667}$"
        ),
    )

    for law, colour in zip(ILLUSTRATIVE_SENSITIVITY_LAWS, ILLUSTRATIVE_COLOURS, strict=True):
        ax.plot(
            flux,
            law.regression_rate_mm_s(flux),
            color=colour,
            linewidth=1.6,
            linestyle="--",
            zorder=3,
            label=rf"ILLUSTRATIVE ONLY: $n = {law.exponent:g}$",
        )

    # Study sample points on the sourced law.
    study_flux = np.array(
        [
            OperatingPoint(GRAIN, m, GRAIN.initial_port_diameter_m).oxidizer_mass_flux_si
            for m in OXIDIZER_MASS_FLOWS_KG_S
        ]
    )
    ax.plot(
        study_flux,
        LAW.regression_rate_mm_s(study_flux),
        "o",
        color=SOURCED_COLOUR,
        markersize=5.5,
        markeredgecolor="white",
        markeredgewidth=0.8,
        zorder=6,
        label=r"study points, $\dot{m}_{ox} = 0.045$–$0.150$ kg/s at $D_p = 40$ mm",
    )

    rep = OperatingPoint(GRAIN, REPRESENTATIVE_M_DOT_OX_KG_S, GRAIN.initial_port_diameter_m)
    rep_flux = rep.oxidizer_mass_flux_si
    rep_rate = rep.regression_rate_si(LAW) * 1e3
    ax.plot(rep_flux, rep_rate, "*", color=ACCENT, markersize=16, zorder=7,
            label="representative point")
    ax.annotate(
        rf"$\dot{{m}}_{{ox}} = {REPRESENTATIVE_M_DOT_OX_KG_S:.3f}$ kg/s"
        "\n"
        rf"$G_{{ox}} = {rep_flux:.1f}$ kg/(m$^2\,$s)"
        "\n"
        rf"$\dot{{r}} = {rep_rate:.3f}$ mm/s",
        xy=(rep_flux, rep_rate),
        xytext=(rep_flux + 18.0, rep_rate - 0.27),
        fontsize=9,
        color=ACCENT,
        arrowprops={"arrowstyle": "-", "color": ACCENT, "linewidth": 0.9},
        bbox={"boxstyle": "round,pad=0.34", "facecolor": "white", "edgecolor": ACCENT,
              "alpha": 0.92, "linewidth": 0.8},
    )

    ax.axvline(ILLUSTRATIVE_ANCHOR_FLUX_SI, color="0.45", linewidth=1.0, linestyle=":", zorder=2)
    ax.text(
        ILLUSTRATIVE_ANCHOR_FLUX_SI - 1.5,
        0.335,
        rf"illustrative anchor, $G_{{ox}} = {ILLUSTRATIVE_ANCHOR_FLUX_SI:.0f}$",
        rotation=90,
        ha="right",
        va="bottom",
        fontsize=8,
        color="0.35",
    )

    ax.set_xlabel(r"oxidizer mass flux  $G_{ox} = \dot{m}_{ox}/A_{port}$   [kg/(m$^2\,$s)]")
    ax.set_ylabel(r"fuel regression rate  $\dot{r}$   [mm/s]")
    ax.set_title(
        "Figure B — Regression rate versus oxidizer mass flux,  "
        r"$\dot{r} = a\,G_{ox}^{\,n}$",
        pad=10,
    )
    ax.set_xlim(flux_lo, flux_hi)
    ax.set_ylim(0.30, 1.62)
    ax.legend(loc="upper left", ncols=1)

    fig.text(
        0.5,
        0.105,
        "Dashed curves use ILLUSTRATIVE coefficients: the exponent is varied over the 0.5–0.8 band "
        "reported for hybrids and the\ncoefficient is fixed only by pivoting about the marked "
        "anchor flux. They are NOT fits to N$_2$O/HTPB measurements.",
        ha="center",
        va="center",
        fontsize=9,
        color=ACCENT,
        fontweight="bold",
    )
    fig.text(0.5, 0.028, DISCLAIMER_MODEL, ha="center", va="center", fontsize=8.5, color="0.25")

    return _save(fig, "fig_b_regression_vs_flux.png")


# ---------------------------------------------------------------------------
# Figure C - port-size sensitivity at fixed oxidizer mass flow
# ---------------------------------------------------------------------------


def figure_c_port_sensitivity() -> Path:
    m_dot_ox = REPRESENTATIVE_M_DOT_OX_KG_S
    d_lo, d_hi = 0.035, 0.075
    diameters = np.linspace(d_lo, d_hi, 401)
    areas = np.pi * diameters**2 / 4.0
    fluxes = m_dot_ox / areas
    rates_mm_s = LAW.regression_rate_mm_s(fluxes)

    marker_fluxes = np.array(
        [OperatingPoint(GRAIN, m_dot_ox, d).oxidizer_mass_flux_si for d in PORT_DIAMETERS_M]
    )
    marker_rates = np.array(
        [OperatingPoint(GRAIN, m_dot_ox, d).regression_rate_si(LAW) * 1e3 for d in PORT_DIAMETERS_M]
    )
    marker_d_mm = np.array(PORT_DIAMETERS_M) * 1e3

    fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.6), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.975, top=0.865, bottom=0.185, hspace=0.16)
    ax_flux, ax_rate = axes

    band_lo, band_hi = LAW.valid_flux_range_si

    # -- upper panel: flux ---------------------------------------------------
    ax_flux.axhspan(band_lo, band_hi, color=BAND_COLOUR, alpha=0.32, zorder=0,
                    label="source data range")
    ax_flux.axhline(
        band_lo,
        color=ACCENT,
        linestyle="--",
        linewidth=1.1,
        zorder=2,
        label=rf"source data lower bound, {band_lo:.0f} kg/(m$^2\,$s)",
    )
    ax_flux.plot(diameters * 1e3, fluxes, color=SOURCED_COLOUR, zorder=3,
                 label=r"$G_{ox} = 4\dot{m}_{ox}/(\pi D_p^2)$")
    ax_flux.plot(marker_d_mm, marker_fluxes, "o", color=SOURCED_COLOUR, markersize=6,
                 markeredgecolor="white", markeredgewidth=0.8, zorder=5,
                 label="study port diameters")
    for d_mm, value in zip(marker_d_mm, marker_fluxes, strict=True):
        ax_flux.annotate(
            f"{value:.1f}",
            xy=(d_mm, value),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=8.5,
            color=SOURCED_COLOUR,
        )
    ax_flux.set_ylabel(r"$G_{ox}$   [kg/(m$^2\,$s)]")
    ax_flux.set_ylim(15.0, 110.0)
    ax_flux.legend(loc="upper right")
    ax_flux.set_title(
        rf"Figure C — Instantaneous port-size sensitivity at fixed $\dot{{m}}_{{ox}} = "
        rf"{m_dot_ox:.3f}$ kg/s",
        pad=10,
    )

    # -- lower panel: regression rate ---------------------------------------
    in_range = (fluxes >= band_lo) & (fluxes <= band_hi)
    ax_rate.plot(diameters * 1e3, rates_mm_s, color=SOURCED_COLOUR, zorder=3,
                 label=r"$\dot{r} = a\,G_{ox}^{\,n}$  (sourced coefficients)")
    ax_rate.fill_between(
        diameters * 1e3,
        0.0,
        rates_mm_s,
        where=~in_range,
        color=ACCENT,
        alpha=0.10,
        zorder=1,
        label="extrapolated beyond source data range",
    )
    ax_rate.plot(marker_d_mm, marker_rates, "o", color=SOURCED_COLOUR, markersize=6,
                 markeredgecolor="white", markeredgewidth=0.8, zorder=5,
                 label="study port diameters")
    for d_mm, value in zip(marker_d_mm, marker_rates, strict=True):
        ax_rate.annotate(
            f"{value:.3f}",
            xy=(d_mm, value),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=8.5,
            color=SOURCED_COLOUR,
        )
    ax_rate.set_xlabel(r"port diameter  $D_p$   [mm]")
    ax_rate.set_ylabel(r"$\dot{r}$   [mm/s]")
    ax_rate.set_xlim(d_lo * 1e3, d_hi * 1e3)
    ax_rate.set_ylim(0.44, 1.00)
    ax_rate.legend(loc="upper right")

    fig.text(
        0.5,
        0.098,
        r"Larger port $\rightarrow$ larger $A_{port}$ $\rightarrow$ lower $G_{ox}$ "
        r"$\rightarrow$ lower $\dot{r}$, at fixed prescribed $\dot{m}_{ox}$."
        "\nInstantaneous what-if port sizes only — this is not a burn-time integration; "
        "port evolution is deferred to Milestone 2.",
        ha="center",
        va="center",
        fontsize=9.5,
    )
    fig.text(0.5, 0.028, DISCLAIMER_MODEL, ha="center", va="center", fontsize=8.5, color="0.25")

    return _save(fig, "fig_c_port_sensitivity.png")


def main() -> int:
    _configure_matplotlib()
    print("Generating Milestone 1 figures (deterministic)...")
    for builder in (figure_a_geometry, figure_b_regression, figure_c_port_sensitivity):
        path = builder()
        print(f"  wrote {path.relative_to(FIGURES_DIR.parent)}  ({path.stat().st_size} bytes)")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
