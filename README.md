# Hybrid Rocket Motor — N₂O/HTPB Regression Rate Study

A rigorous, reduced-order engineering model of a **generic** nitrous-oxide /
HTPB hybrid rocket motor.

> **Scope.** This is an educational, reduced-order study. It is **not** a motor
> design, a fabrication guide, a test procedure, or a flight-hardware model. The
> geometry is illustrative and is not taken from any real motor. Nothing here has
> been validated against hardware.

---

## Engineering question

> For a generic single-port cylindrical HTPB grain fed a **prescribed** nitrous-oxide
> mass flow, what oxidizer mass flux does the port see, what fuel regression rate
> does an experimentally sourced N₂O/HTPB power law predict, and what fuel mass
> flow and mixture ratio follow from it?

## Milestone 1 scope

Milestone 1 builds the foundation only:

* validated grain geometry,
* prescribed-oxidizer operating-point bookkeeping,
* the empirical regression law `ṙ = a · G_ox ⁿ` with rigorous unit handling,
* fuel mass-flow and O/F bookkeeping,
* a prescribed-flow regression study with port-size sensitivity.

**Milestone 1 does not predict thrust.** It does not model N₂O tank
thermodynamics, injector flow, chamber pressure, combustion efficiency,
equilibrium chemistry, `c*`, nozzle flow, thrust, structures, thermal response,
ignition or fabrication. See [DESIGN.md](DESIGN.md) §10.

## Model chain

```
prescribed ṁ_ox  ─┐
                  ├─►  G_ox = ṁ_ox / A_port          [kg/(m²·s)]
port diameter D_p ┘           A_port = π D_p² / 4

G_ox  ─────────────►  ṙ = a · G_ox ⁿ                 [m/s]     (sourced correlation)

ṙ, D_p, L, ρ_f ────►  ṁ_f = ρ_f · A_burn · ṙ         [kg/s]    A_burn = π D_p L

ṁ_ox, ṁ_f ─────────►  O/F = ṁ_ox / ṁ_f               [–]       (bookkeeping only)
```

The chain stops there. O/F is **not** carried into flame temperature, `c*`,
chamber pressure or nozzle performance at this milestone.

## Representative geometry (illustrative)

| Parameter | Symbol | Value |
| --- | --- | --- |
| Grain length | `L` | 0.40 m |
| Initial port diameter | `D_p,0` | 0.040 m |
| Outer grain diameter | `D_o` | 0.090 m |
| Ports | — | 1, circular, central |
| Fuel density | `ρ_f` | 930 kg/m³ (HTPB-like nominal) |
| Initial web | — | 0.025 m |
| Initial port area | `A_port,0` | 1.2566 × 10⁻³ m² |
| Initial burning area | `A_burn,0` | 5.0265 × 10⁻² m² |
| Loaded fuel mass | — | 1.8991 kg |

*Generic illustrative geometry — not a manufacturing drawing.*

## Regression law

```
ṙ = a · G_ox ⁿ
```

### Coefficient status: **SOURCED** (with a clearly labelled illustrative supplement)

The primary coefficients are taken **verbatim** from a peer-reviewed experimental
study of an HTPB/N₂O motor. No fitting or calibration was performed by this
project.

| | |
| --- | --- |
| Correlation as published | `ṙ = 0.3977 · G_o^0.3667` |
| Source units | `ṙ` in **mm/s**, `G_o` in **g/(cm²·s)** |
| Converted, SI | `ṙ [m/s] = 1.709446806 × 10⁻⁴ · (G_ox [kg/(m²·s)])^0.3667` |
| Propellant | **N₂O / HTPB** — validated for this specific combination |
| Source data range | `G_ox` = 3.5–12 g/(cm²·s) = **35–120 kg/(m²·s)** |
| Reference | H. Rezaei, M. R. Soltani, A. R. Mohammadi, *Experimental study of fuel regression rate in an HTPB/N₂O hybrid rocket motor*, **Scientia Iranica B** 25(1) (2018) 253–265, Eq. (10). DOI [10.24200/sci.2017.4317](https://doi.org/10.24200/sci.2017.4317) |

**Two discrepancies are reported rather than tuned away** (full detail in
[DESIGN.md](DESIGN.md) §4.6):

1. The sourced exponent `n = 0.3667` sits **below** the 0.5–0.8 band reported for
   most hybrid systems, so `ṙ` responds only weakly to flux in this model. Three
   **ILLUSTRATIVE** laws with `n = 0.5, 0.62, 0.8` are provided to show what a
   steeper exponent would imply; their coefficients come from pivoting about a
   single disclosed anchor flux and are **not** fits to N₂O/HTPB data.
2. The correlation was fitted on a 250 mm grain, and the same paper reports a
   strong grain-length dependence. For this 400 mm grain that implies predictions
   may be ~15 % low. The correction is **not** applied — 400 mm is itself outside
   the range the length fit was derived over.

Source audit, unit derivation and the full limitations list: [DESIGN.md](DESIGN.md) §4.

## Representative results

At the representative operating point — prescribed `ṁ_ox = 0.100 kg/s`,
`D_p = 40 mm`:

| Quantity | Value |
| --- | --- |
| `A_port` | 1.2566 × 10⁻³ m² |
| `G_ox` | **79.58 kg/(m²·s)** (7.958 g/(cm²·s)) — inside the source data range |
| `ṙ` | **0.8509 mm/s** |
| `A_burn` | 5.0265 × 10⁻² m² |
| `ṁ_f` | **0.039777 kg/s** |
| `O/F` | **2.514** |

Prescribed oxidizer mass-flow sweep at `D_p = 40 mm`:

| `ṁ_ox` [kg/s] | `G_ox` [kg/(m²·s)] | `ṙ` [mm/s] | `ṁ_f` [kg/s] | `O/F` |
| ---: | ---: | ---: | ---: | ---: |
| 0.045 | 35.81 | 0.6349 | 0.02968 | 1.516 |
| 0.060 | 47.75 | 0.7055 | 0.03298 | 1.819 |
| 0.075 | 59.68 | 0.7657 | 0.03579 | 2.095 |
| 0.090 | 71.62 | 0.8186 | 0.03827 | 2.352 |
| 0.105 | 83.56 | 0.8663 | 0.04050 | 2.593 |
| 0.120 | 95.49 | 0.9097 | 0.04253 | 2.822 |
| 0.135 | 107.43 | 0.9499 | 0.04440 | 3.040 |
| 0.150 | 119.37 | 0.9873 | 0.04615 | 3.250 |

All rows are inside the source's reported flux range.

Port-size sensitivity at fixed `ṁ_ox = 0.100 kg/s` — the expected coupling
**larger port → larger `A_port` → lower `G_ox` → lower `ṙ`**:

| `D_p` [mm] | `A_port` [m²] | `G_ox` [kg/(m²·s)] | `ṙ` [mm/s] | Source range |
| ---: | ---: | ---: | ---: | --- |
| 40 | 1.2566 × 10⁻³ | 79.58 | 0.8509 | in-range |
| 50 | 1.9635 × 10⁻³ | 50.93 | 0.7224 | in-range |
| 60 | 2.8274 × 10⁻³ | 35.37 | 0.6320 | in-range |
| 70 | 3.8485 × 10⁻³ | 25.98 | 0.5645 | **EXTRAPOLATED** |

Note that `ṁ_f` does **not** fall in step with `ṙ`: the burning area grows with
`D_p` at the same time, and with `n = 0.3667` the net effect is a slight *rise*
in fuel flow. This is explained in [DESIGN.md](DESIGN.md) §5.

Regression rates of order 0.6–1.0 mm/s are consistent with the ~1 mm/s scale
reported for HTPB/N₂O in the literature.

## Installation

```bash
pip install -e ".[dev]"
```

Requires Python ≥ 3.11. The `dev` extra installs `pytest` and `ruff`, which is
everything the verification workflow below needs.

## Running the study

```bash
python scripts/regression_study.py
```

Prints the configuration, the oxidizer-flow sweep, the port-size sensitivity, the
illustrative exponent sensitivity, the observations, and an engineering sanity
audit. Exits non-zero if any sanity check fails.

Longhand cross-check of the whole model chain against hand arithmetic:

```bash
python scripts/manual_check.py
```

Regenerate the figures (deterministic — byte-identical on repeated runs within a
given environment; see note below):

```bash
python scripts/generate_m1_figures.py
```

All figure inputs are fixed constants, the Agg backend is forced, every rendering
rcParam is pinned and no timestamp metadata is written, so repeated runs produce
byte-identical PNGs. Hashes are **not** expected to match across environments with
different matplotlib or FreeType builds, because PNG glyph rasterisation depends on
the FreeType version; the figure content is unchanged. Verified with
matplotlib 3.10.9 / FreeType 2.6.1.

## Tests

```bash
pytest -W error -q
ruff check .
```

167 tests. Expected values are written independently of the production code —
independently arranged algebra, hand-computed literals, an independent unit
conversion route, and all 18 published `(G_ox, ṙ)` measurements from the source
paper's own data table.

## Figures

| Figure | File |
| --- | --- |
| **A** — Generic grain geometry schematic | [`figures/fig_a_grain_geometry.png`](figures/fig_a_grain_geometry.png) |
| **B** — Regression rate vs oxidizer mass flux | [`figures/fig_b_regression_vs_flux.png`](figures/fig_b_regression_vs_flux.png) |
| **C** — Port-size sensitivity at fixed `ṁ_ox` | [`figures/fig_c_port_sensitivity.png`](figures/fig_c_port_sensitivity.png) |

![Generic hybrid grain geometry](figures/fig_a_grain_geometry.png)

![Regression rate versus oxidizer mass flux](figures/fig_b_regression_vs_flux.png)

![Port-size sensitivity](figures/fig_c_port_sensitivity.png)

## Repository layout

```
src/hybrid_rocket_motor/
    geometry.py          grain geometry, areas, volumes, masses, validation
    operating_point.py   prescribed-ṁ_ox snapshot: G_ox, ṁ_f, O/F bookkeeping
    regression.py        ṙ = a G_ox ⁿ with explicit source-unit handling
tests/                   independent verification (167 tests)
scripts/
    manual_check.py      longhand arithmetic cross-check + scope guard
    regression_study.py  the Milestone 1 study and sanity audit
    generate_m1_figures.py   deterministic figure generation
figures/                 three PNG figures
DESIGN.md                conventions, source audit, verification strategy
```

## Limitations

* Coefficients come from a **single** experimental campaign on a **single** motor,
  with the two discrepancies described above.
* The study grain (400 mm) is longer than the grain the correlation was fitted on
  (250 mm).
* Extrapolation outside 35–120 kg/(m²·s) is **flagged but not prevented**; the
  70 mm port case genuinely sits outside the source data.
* Regression is treated as spatially uniform, diffusion-limited and
  pressure-independent — a lumped space–time averaged approximation.
* No uncertainty quantification; point values only.
* Quasi-static snapshots only — no time integration, no transient, no burnout.
* **Nothing here has been validated against hardware.**

## Next milestone

**Milestone 2 — transient port-radius evolution under prescribed oxidizer-flow
histories:** integrate `dD_p/dt = 2 ṙ` for a prescribed `ṁ_ox(t)`, handle fuel
depletion / burnout when `D_p → D_o`, and produce time histories of `ṙ`, `ṁ_f`
and `O/F` — still without chamber-pressure or nozzle coupling.

## License

MIT
