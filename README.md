# Hybrid Rocket Motor — Reduced-Order N₂O/HTPB Performance Study

A coupled, reduced-order engineering model of a **generic** nitrous-oxide / HTPB
hybrid rocket motor: fuel regression, transient port evolution, chamber and
nozzle, self-pressurising tank blowdown, and O/F-dependent equilibrium
thermochemistry — built in layers, each removing one assumption made by the last,
and finished with a deterministic robustness audit of what survives.

> **Scope.** This is an educational, reduced-order study of a **generic conceptual
> motor**. It is **not experimentally validated**, **not hardware design guidance**,
> and not a fabrication or test procedure. The injector is represented as an
> **effective flow area only** — no hole count, hole size or plate geometry. The
> thermochemistry comes from a **frozen generated table**, so the **runtime does
> not require CEA**. There is **no ignition, structural or thermal design** and
> **no flight prediction**. Nothing here has been checked against a firing.

---

## Engineering question

> For a generic single-port HTPB grain fed by a self-pressurising nitrous-oxide
> tank, what thrust history does a coupled reduced-order model predict — and
> which of its conclusions are robust to the modelling assumptions rather than
> artefacts of them?

## Main result

**Two modelling assumptions each change the answer in a different way, and one of
them changes its sign.**

| Model | Thrust | Slope | Total impulse |
| --- | --- | --- | --- |
| Prescribed oxidizer flow | 310.0 → 333.5 N | **+7.57 %** (rising) | 13 525.7 N·s |
| Coupled tank blowdown | 313.7 → 286.8 N | **−8.60 %** (falling) | 13 018.6 N·s |
| Coupled + O/F-dependent chemistry | 290.8 → 253.2 N | **−12.91 %** (falling further) | 11 627.0 N·s |

* **Coupling the tank reverses the trend.** Holding oxidizer flow constant
  predicts thrust *rising* through the burn. With the tank coupled it *falls* —
  the tank blows down faster than the growing burning area can compensate. Total
  impulse barely moves (−3.75 %) while the shape of the curve inverts: a model can
  get the integral roughly right and the trend backwards.
* **Variable chemistry deepens the fall without changing its sign,** to 1.50× the
  magnitude, and drops the whole curve (−10.69 % impulse). As the port opens the
  mixture ratio falls, and on the fuel-rich side of the CEA table a falling `O/F`
  drags `c*` down with it. `c*` moves **21×** more than `γ` over the burn, so this
  is a `c*` effect acting through chamber pressure, not a nozzle-`γ` effect.
* **What is robust is the shape, not the level.** Across every perturbation tested,
  the falling-thrust sign never flips. But total impulse moves ±13.6 % on the
  injector-model choice alone. The sign is a conclusion; the number is not.

Full numbers in [RESULTS.md](RESULTS.md); the verification behind them in
[VERIFICATION.md](VERIFICATION.md).

## Representative conceptual reference case

Deliberately **not** called an optimised motor. Nothing was tuned to a performance
target. Each input is labelled by what kind of thing it is:

| Input | Kind | Value |
| --- | --- | --- |
| Grain | **illustrative** | 0.40 m long, 40 → 90 mm single circular port, `ρ_f` = 930 kg/m³ |
| Regression law | **sourced** | Rezaei et al. (2018): `a_SI` = 1.709447 × 10⁻⁴, `n` = 0.3667 |
| N₂O properties | **sourced** | CoolProp / Lemmon & Span (2006) equation of state |
| Tank | **illustrative** | 10 L rigid adiabatic, 80 % liquid fill at 293.15 K → 50.525 bar, 6.5968 kg |
| Injector | **effective parameter** | Dyer/NHNE, `A_eff` = 3.50 mm², `C_d` = 0.66 |
| Nozzle | **illustrative** | `D_t` = 10 mm, ε = 4.0, ideal 1-D expansion |
| Thermochemistry | **sourced** | NASA CEA (RP-1311), frozen 113 × 15 table |
| `η_c*` | **model-form assumption** | 0.96, standing in for every combustion loss |

| Output | Value |
| --- | --- |
| `O/F` | 2.557 → 1.856 |
| Tank pressure | 50.53 → 38.09 bar |
| Chamber pressure | 24.67 → 21.85 bar |
| Thrust | 290.8 → 253.2 N (−12.91 %) |
| Burn duration | 42.52 s |
| Total impulse | 11 627 N·s |
| Equivalent `I_sp` | 199.5 s |
| Terminal event | grain burnout |

---

## Model architecture

### Regression

The empirical power law, taken verbatim from a peer-reviewed N₂O/HTPB campaign
and converted into SI with the unit chain made explicit:

```
ṙ = a · G_ox ⁿ          a = 0.3977 mm/s, n = 0.3667   [source units: g/(cm²·s)]
G_ox = ṁ_ox / A_port
ṁ_f  = ρ_f · A_burn · ṙ
O/F  = ṁ_ox / ṁ_f
```

The coefficient is **sourced**, not fitted here. The implementation reproduces all
18 published `(G_ox, ṙ)` firings from the source paper's own data table to 10.9 %
worst-case and 5.0 % mean — that is agreement with the *scatter of that
experiment*, not a validated predictive accuracy for any other motor.

The correlation's validated flux band is tracked, and operation outside it is
flagged rather than silently extrapolated.

### Port evolution

The grain state is one scalar — the port radius — driven by the same law:

```
dr_p/dt = ṙ(G_ox(r_p, t))
```

integrated with a terminal burnout event at the outer radius. A closed-form
solution exists for constant oxidizer flow and is used as an independent check:

```
r(t) = [ r₀^(2n+1) + (2n+1) a (ṁ/π)ⁿ t ]^(1/(2n+1))
```

### Chamber and nozzle

A quasi-steady chamber mass balance closed by the definition of `c*`, and a 1-D
isentropic converging–diverging nozzle:

```
p_c = (ṁ_ox + ṁ_f) c* / A_t
ε   = (1/M_e)[(2/(γ+1))(1 + (γ−1)/2 · M_e²)]^((γ+1)/(2(γ−1)))
F   = ṁ V_e + (p_e − p_a) A_e
```

There is **no nozzle efficiency**: the expansion is ideal, shock-free and
one-dimensional, with a single `γ` from chamber to exit. Flow separation is
flagged by the Summerfield criterion, not modelled.

### N₂O tank and feed

A rigid, adiabatic, spatially uniform **saturated-equilibrium** tank whose
pressure is `p_sat(T)` and whose temperature falls as liquid is withdrawn:

```
dm/dt = −ṁ_ox
dU/dt = −ṁ_ox h_l(T)
```

closed by the volume constraint. Nothing prescribes a pressure-versus-time curve —
the pressure history is an output.

Three sourced injector models, all quoted from the same Stanford/NASA Ames paper:

```
SPI   ṁ = C_d A √(2 ρ Δp)                          their Eq. (2)
HEM   ṁ = C_d A ρ₂ √(2 (h₁ − h₂))                  their Eq. (3)
Dyer  ṁ = κ/(1+κ)·ṁ_SPI + 1/(1+κ)·ṁ_HEM            their Eq. (9)
```

For a saturated tank the non-equilibrium parameter `κ` is exactly 1, so the Dyer
blend is an equal average of the two limits. The authors warn that SPI "is often
applied inappropriately" to nitrous; the sensitivity study quantifies what that
choice is worth, and it is the single largest lever in the model.

The oxidizer flow is then the root of

```
residual(ṁ_ox) = injector_flow(tank, p_c(ṁ_ox, r_p)) − ṁ_ox
```

solved by Brent's method on a closed-form bracket.

### Thermochemistry

`c*`, `γ`, `T_c` and molar mass are interpolated from a **frozen NASA CEA
equilibrium table** (RP-1311, via the `rocketcea` wrapper) at the instantaneous
mixture ratio and chamber pressure — 1695 points over `O/F` 1.20–4.00 and
10–45 bar, committed to the repository.

Because `c*` now depends on the pressure it helps set, the chamber relation is
implicit and is solved as a bracketed root:

```
p_c = ṁ c*(O/F, p_c) / A_t
```

nested inside the outer feed root. No fixed-point iteration appears anywhere.

Three details carry most of the weight:

* **Bilinear interpolation, deliberately not a spline.** At the soot boundary
  (`O/F` ≈ 3.01) `γ` steps ~4.5 % between adjacent nodes; a spline would overshoot
  there, inventing properties the solver never produced.
* **`M, (1/n)`, not `MW`.** Where equilibrium produces condensed carbon, CEA prints
  a gas-phase molecular weight *and* a total-mass-per-mole-of-gas. A single-gas
  nozzle needs the latter — the one reproducing CEA's own sound speed,
  `M = γ R_u T / a²`. At the reference point that is 21.057 vs 19.543 g/mol,
  7.75 % apart.
* **The model refuses to extrapolate.** Outside the tabulated rectangle the solver
  fails with a named status and the integration stops. `O/F` is never clipped into
  the table and `p_c` is never clamped to an edge.

`η_c*` is applied **once**, downstream, to the ideal tabulated `c*` — it is not
double-counted.

---

## Model hierarchy

The three levels differ only in which assumption they make. They are **not of
equal fidelity**, and none is validated against hardware.

### M3 — prescribed flow

Oxidizer flow is an **input**; combustion properties are prescribed constants.
Thrust *rises* +7.57 % through the burn as the port opens.

### M4 — coupled blowdown

Oxidizer flow becomes an **output** of the tank/injector/chamber loop. Properties
are still frozen. Thrust *falls* −8.60 %. Burn ends on grain burnout at 42.917 s
with 28 % of the liquid unused; a 2.5 L variant instead runs out of liquid at
16.035 s.

### M5 — variable thermochemistry

Properties follow the operating point. Thrust falls −12.91 % — the same sign,
1.50× the magnitude — and the whole curve drops.

**Decomposition of the total effect:**

| Effect | Impulse | Burn time | Slope [pp] |
| --- | --- | --- | --- |
| **A.** feed/blowdown coupling (M4 − M3) | −3.75 % | +2.82 % | **−16.17** |
| **B.** variable thermochemistry (M5 − M4) | **−10.69 %** | −0.92 % | −4.32 |
| **C.** total (M5 − M3) | −14.04 % | +1.87 % | −20.49 |

The two steps are not interchangeable: feed coupling flips the slope while barely
moving impulse, and the chemistry does the opposite.

---

## Sensitivity / robustness

**Deterministic, one factor at a time.** Each family changes exactly one
assumption to stated alternative values and reruns the whole coupled model. No
probability distribution is assumed for any input, and no sampling is performed —
this is **not** an uncertainty quantification.

| Rank | Total impulse | Thrust decay | Burn duration |
| --- | --- | --- | --- |
| 1 | injector model **13.59 %** | injector model **52.77 %** | regression coefficient `a` **11.40 %** |
| 2 | regression coefficient `a` 9.03 % | regression coefficient `a` 18.31 % | injector model 6.85 % |
| 3 | initial tank temperature 7.35 % | initial tank temperature 18.19 % | initial tank temperature 4.32 % |
| 4 | injector effective area 5.38 % | injector effective area 14.07 % | regression exponent `n` 3.27 % |
| 5 | regression exponent `n` 2.61 % | regression exponent `n` 9.31 % | injector effective area 3.12 % |
| 6 | `c*` efficiency 0.07 % | `c*` efficiency 1.86 % | `c*` efficiency 0.30 % |
| 7 | chemistry grid resolution 0.04 % | chemistry `p_c` dependence 0.06 % | chemistry grid resolution 0.004 % |
| 8 | chemistry `p_c` dependence 0.003 % | chemistry grid resolution 0.04 % | chemistry `p_c` dependence 0.000 % |

The dominant assumption is **not the same for every output**, so no single
universal ranking is claimed.

Three findings worth stating plainly:

* **The injector model outranks every propellant property.** Which reduced-order
  injector correlation you pick matters more here than the chemistry does.
* **`η_c*` barely moves total impulse (0.07 %)** despite being the largest
  unvalidated number in the model — the pressure-fed system self-compensates,
  because a higher `c*` raises `p_c`, which throttles the injector, while the
  grain-limited burn fixes the fuel mass. It still moves chamber pressure ~7 %.
* **The table's own resolution is irrelevant** at the committed fidelity (under
  0.04 % everywhere). The chemistry matters through *what* it says about `c*(O/F)`,
  not how finely it is tabulated.

**Qualitative conclusions, with support counts** — "robust within tested
deterministic family" means every case in *this* set agreed. It is not a
probability and says nothing about untested regions.

| Conclusion | Support |
| --- | --- |
| Coupling the tank reverses the rising-thrust trend | 5/5 |
| Variable thermochemistry preserves the falling-thrust sign | 15/15 |
| Variable thermochemistry deepens the decay | 5/5 |
| `c*` variation dominates the chemistry effect over `γ` | 15/15 |
| Feed coupling changes thrust shape more than chemistry | 5/5 |
| `O/F` drifts downward during the burn | 15/15 |
| The 10 L case burns through before liquid runs out | 15/15 |
| A small tank can instead deplete liquid first | 1/1 |

One perturbation was **invalid and is reported, not repaired**: the bare HEM limit
gives too little flow for the chamber to reach the table's pressure floor, so no
firing state is established. It is excluded from the rankings and was not clamped
back into validity.

---

## Verification

**59 independent checks, all passing.** The expected side of each is built from
longhand algebra, closed-form solutions, alternate numerical routes (plain
bisection, secant, trapezoid), published measurements, or raw CSV columns parsed
with the standard library — never by calling the same function twice.

| Kind | Checks | Worst residual |
| --- | --- | --- |
| Algebraic identity | 49 | 9.313 × 10⁻¹⁰ Pa (≈ 4 × 10⁻¹⁶ relative) |
| Discretisation | 8 | 2.911 × 10⁻⁷ relative |
| Empirical (vs published firings) | 2 | 1.085 × 10⁻¹ |

Discretisation residuals are integrator-tolerance quantities and are reported
separately rather than folded into one headline figure.

**857 tests** (167 M1, 119 M2, 215 M3, 211 M4, 142 M5, 3 M6 regression), zero
warnings under `-W error`, `ruff check` clean. All 29 figures and both committed
CEA data files regenerate byte-identically.

Details, including the one defect found and fixed during the final audit, are in
[VERIFICATION.md](VERIFICATION.md).

---

## Featured figures

| Figure | File |
| --- | --- |
| **Model hierarchy** — flow, pressure and thrust for M3 / M4 / M5 | [`figures/final_model_hierarchy.png`](figures/final_model_hierarchy.png) |
| **Thermochemistry evolution** — which property drives the added decay | [`figures/final_thermochemistry_evolution.png`](figures/final_thermochemistry_evolution.png) |
| **Sensitivity ranking** — deterministic one-factor tornado | [`figures/final_sensitivity_ranking.png`](figures/final_sensitivity_ranking.png) |
| **Reference case** — the full conceptual reference case | [`figures/final_reference_case.png`](figures/final_reference_case.png) |
| **Coupling decomposition** — feed effect vs chemistry effect | [`figures/final_coupling_decomposition.png`](figures/final_coupling_decomposition.png) |

Twenty-four further figures cover each milestone in detail — grain geometry,
regression sweeps, port growth, mass closure, chamber and nozzle states, ambient
and property sensitivities, tank blowdown, injector characteristics, the CEA table
itself and its interpolation convergence. See [`figures/`](figures/).

---

## Reproducing the analysis

```bash
pip install -e ".[dev]"
```

Python ≥ 3.11. Runtime dependencies are `numpy`, `scipy`, `CoolProp` and
`matplotlib`. **The runtime does not require CEA** — the thermochemistry table is
committed. Only *regenerating* it needs `pip install -e ".[cea]"`, which pulls in
`rocketcea` and therefore a Fortran compiler.

Verification gate:

```bash
pytest -W error -q && ruff check .
```

The two audit entry points:

```bash
python scripts/independent_audit.py
```

```bash
python scripts/final_robustness_study.py
```

Per-milestone studies:

```bash
python scripts/regression_study.py
```

```bash
python scripts/transient_regression_study.py
```

```bash
python scripts/thrust_prediction_study.py
```

```bash
python scripts/blowdown_thrust_study.py
```

```bash
python scripts/variable_thermochemistry_study.py
```

Longhand arithmetic cross-check of the Milestone 1 chain:

```bash
python scripts/manual_check.py
```

Regenerate every figure (deterministic within one environment):

```bash
python scripts/generate_final_figures.py
```

Confirm the committed CEA table is what the generator produces (needs the `cea`
extra):

```bash
python scripts/generate_thermochemistry_table.py --check
```

All figure inputs are fixed constants, the Agg backend is forced, every rendering
rcParam is pinned and no timestamp metadata is written, so repeated runs produce
byte-identical PNGs. Hashes are **not** expected to match across environments with
different matplotlib or FreeType builds, because PNG glyph rasterisation depends
on the FreeType version; the figure content is unchanged. Verified with
matplotlib 3.10.9 / FreeType 2.6.1.

---

## Repository structure

```
src/hybrid_rocket_motor/
    geometry.py            grain geometry, areas, volumes, masses, validation
    operating_point.py     prescribed-ṁ_ox snapshot: G_ox, ṁ_f, O/F bookkeeping
    regression.py          ṙ = a G_ox ⁿ with explicit source-unit handling
    transient.py           flow histories, ODE solver, burnout event, closed form
    nozzle.py              geometry, area–Mach solve, isentropic exit state
    chamber.py             prescribed combustion properties, quasi-steady p_c
    performance.py         chamber + nozzle coupling, thrust, total impulse
    nitrous_properties.py  saturated N₂O properties (CoolProp / Lemmon–Span)
    tank.py                rigid adiabatic saturated-equilibrium N₂O tank
    injector.py            SPI / HEM / Dyer-NHNE injector flow models
    feed_system.py         coupled injector–chamber–regression root solve
    blowdown.py            coupled time integration with terminal events
    thermochemistry.py     frozen CEA table, bilinear interpolation, refusal
    variable_chamber.py    implicit p_c = ṁ c*(O/F, p_c) / A_t closure
    variable_feed_system.py  nested feed / chemistry root solve
    variable_nozzle.py     adapter driving the frozen nozzle from the table
    variable_blowdown.py   coupled integration with table-validity events
    robustness.py          deterministic one-factor sensitivity analysis
data/
    n2o_htpb_equilibrium.csv            1695 CEA points, committed, no timestamp
    n2o_htpb_equilibrium_metadata.json  solver and version provenance
    n2o_htpb_nozzle_chemistry.csv       equilibrium vs frozen expansion
tests/                     independent verification (857 tests)
scripts/
    manual_check.py                    longhand arithmetic cross-check
    regression_study.py                Milestone 1 study and sanity audit
    transient_regression_study.py      Milestone 2 study and convergence
    thrust_prediction_study.py         Milestone 3 study and sensitivity
    blowdown_thrust_study.py           Milestone 4 study and comparison
    variable_thermochemistry_study.py  Milestone 5 study and sensitivity
    independent_audit.py               end-to-end independent verification
    final_robustness_study.py          final robustness audit and synthesis
    generate_thermochemistry_table.py  one-time CEA table generation
    generate_m{1..5}_figures.py        deterministic per-milestone figures
    generate_final_figures.py          deterministic final figures
figures/                   twenty-nine PNG figures
RESULTS.md                 quantitative findings
VERIFICATION.md            independent verification and residuals
DESIGN.md                  derivations, source audit, verification strategy
```

---

## Limitations

* **Generic conceptual geometry.** The grain, tank and nozzle are illustrative and
  taken from no real motor. Nothing was optimised.
* **Coefficients from a single campaign.** The regression law comes from one
  experimental study on one motor; its scatter is ~11 % worst-case.
* **Effective injector area only** — no hole count, hole size or plate geometry,
  and `C_d` = 0.66 is a literature average for *other* injectors.
* **Reduced-order injector correlations**, which the sensitivity study shows to be
  the single largest lever on the answer.
* **Saturated-equilibrium tank model**, with **no vapour-only discharge tail**
  after liquid depletion, **no feed-line dynamics** and **no valve dynamics**.
* **Equilibrium chemistry from a frozen tabulated CEA run.** Equilibrium assumes
  instantaneous, complete reaction: there is **no finite-rate chemistry** and no
  mixing model. Outside the tabulated rectangle the model refuses to run.
* **One `η_c*` for every combustion loss**, unmeasured.
* **Ideal 1-D nozzle** — no losses, no contour design, no shock or separation
  modelling, single `γ` through the expansion.
* **No ignition or chamber-filling transient. No combustion instability.**
* **No structural or thermal design. No hardware validation. No flight
  prediction.**

Being more coupled is not being more validated. Each layer replaced one assumption
with several new ones, and none of them has been checked against a firing.

---

## Development history

The model was built in six milestones, each frozen before the next began and each
verified against the previous one's published results:

| Milestone | What it removed | Headline |
| --- | --- | --- |
| **1** | — | Grain geometry, sourced regression law, prescribed-flow O/F bookkeeping |
| **2** | Static operating point | Transient port evolution, burnout events, mass closure |
| **3** | No thrust prediction | Quasi-steady chamber, 1-D nozzle, thrust — flow still prescribed |
| **4** | Prescribed oxidizer flow | Coupled N₂O tank/injector blowdown; **the thrust-time slope reverses sign** |
| **5** | Prescribed combustion properties | O/F-dependent CEA thermochemistry; **decay deepens 1.50×** |
| **6** | — | Independent audit, deterministic robustness study, portfolio synthesis |

Development stops at Milestone 6. Derivations, source audits and the
verification strategy for every milestone are in [DESIGN.md](DESIGN.md).

## License

MIT
