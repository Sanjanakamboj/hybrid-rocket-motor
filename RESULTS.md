# Results

Quantitative summary of the coupled N₂O/HTPB reduced-order study. Every number
here is reproduced by `scripts/final_robustness_study.py` and
`scripts/independent_audit.py`; the verification detail lives in
[VERIFICATION.md](VERIFICATION.md).

---

## Engineering objective

> For a generic single-port HTPB grain fed by a self-pressurising nitrous-oxide
> tank, what thrust history does a coupled reduced-order model predict, which of
> its conclusions survive reasonable changes to the modelling assumptions, and
> which are artefacts of the assumptions themselves?

The project builds that model in layers, each removing one assumption made by the
previous layer, and finishes by measuring how much each remaining assumption is
worth.

## Reference case

A **representative conceptual reference case** — not an optimised motor, not a
design, and not validated against any hardware. Nothing was tuned to a
performance target.

| Input | Provenance | Value |
| --- | --- | --- |
| Grain | **illustrative** | 0.40 m long, 40 → 90 mm single circular port, ρ_f = 930 kg/m³ |
| Regression law | **sourced** | Rezaei et al. (2018), `a_SI` = 1.709447e-4, `n` = 0.3667 |
| N₂O properties | **sourced** | CoolProp / Lemmon & Span (2006) equation of state |
| Tank | **illustrative** | 10 L rigid adiabatic, 80 % liquid fill at 293.15 K → 50.525 bar, 6.5968 kg |
| Injector | **effective parameter** | Dyer/NHNE, `A_eff` = 3.50 mm², `C_d` = 0.66 — no hole count, no hole size, no plate geometry |
| Nozzle | **illustrative** | `D_t` = 10 mm, ε = 4.0, ideal 1-D expansion |
| Thermochemistry | **sourced** | NASA CEA (RP-1311), frozen 113 × 15 table |
| `η_c*` | **model-form assumption** | 0.96, standing in for every combustion loss |

Resulting behaviour:

| Quantity | Value |
| --- | --- |
| `O/F` | 2.5569 → 1.8560 |
| Tank pressure | 50.525 → 38.09 bar |
| Chamber pressure | 24.671 → 21.846 bar |
| Thrust | 290.76 → 253.21 N (**−12.91 %**) |
| Burn duration | 42.522 s |
| Total impulse | 11 627 N·s |
| Equivalent `I_sp` | 199.5 s |
| Oxidizer / fuel consumed | 4.045 / 1.899 kg |
| Terminal event | grain burnout |

## Regression and port evolution

At the reference oxidizer flow of 0.100 kg/s through the initial 40 mm port:

| Quantity | Value |
| --- | --- |
| `G_ox` | 79.577 kg/(m²·s) |
| `ṙ` | 0.8509 mm/s |
| `ṁ_f` | 0.039777 kg/s |
| `O/F` | 2.514 |

Under a constant 0.100 kg/s the port opens 40 → 90 mm and burns through at
**41.741 s**, matching a closed-form solution re-derived independently to
1.6 × 10⁻¹⁰ s. The correlation reproduces all 18 published firings of its source
paper to 10.9 % worst-case and 5.0 % mean — agreement with the *scatter of that
experiment*, not a validated predictive accuracy.

## Chamber and nozzle reference

Prescribed flows, frozen properties:

| Quantity | Value |
| --- | --- |
| `p_c` | 26.114 bar |
| Exit Mach | 2.619 |
| Thrust | 310.048 N |
| `C_F` | 1.512 |
| `I_sp` | 226.190 s |
| Total impulse over the burn | 13 525.65 N·s |
| Equivalent `I_sp` | 227.103 s |

## Blowdown coupling result

Making the oxidizer flow an output of the tank rather than an input **reverses the
sign of the thrust-time slope**:

| | M3 prescribed flow | M4 coupled blowdown |
| --- | --- | --- |
| Thrust | 310.05 → 333.53 N | 313.73 → 286.76 N |
| Slope | **+7.57 %** (rising) | **−8.60 %** (falling) |
| Burn duration | 41.741 s | 42.917 s |
| Total impulse | 13 525.65 N·s | 13 018.56 N·s |

The tank blows down faster than the growing burning area can compensate. Note
that integrated quantities barely move (−3.75 % impulse) while the *shape* of the
curve inverts — a model can get the total roughly right and the trend backwards.

A 2.5 L tank instead runs out of liquid at 16.035 s with the grain unburnt, so
which terminal event occurs is a configuration property, not a fixed feature.

## Variable thermochemistry result

Letting `c*`, `γ`, `T_c` and molar mass follow the instantaneous `O/F` and `p_c`
does **not** overturn the falling-thrust result — it deepens it:

| | M4 frozen | M5 variable | Change |
| --- | --- | --- | --- |
| Thrust slope | −8.60 % | **−12.91 %** | 1.50× the magnitude |
| Total impulse | 13 018.56 | 11 627.03 N·s | **−10.69 %** |
| Equivalent `I_sp` | 225.73 s | 199.45 s | −11.64 % |
| Burn duration | 42.917 s | 42.522 s | −0.92 % |
| Fuel consumed | 1.8991 kg | 1.8991 kg | 0.00 % |

Fuel consumed is identical because burnout is geometric. Oxidizer consumed *rises*
1.59 %, because the lower chamber pressure lets the injector pass more.

**The mechanism is new.** In M4 the only decay channel was the tank. Here the
falling mixture ratio drags `c*` down with it, and the two add:

| Property | Over the burn |
| --- | --- |
| Delivered `c*` | 1356.2 → 1279.8 m/s (**−5.63 %**) |
| `γ` | 1.23163 → 1.23490 (+0.27 %) |
| `T_c` | 2186.0 → 1925.8 K |

`c*` moves **21×** more than `γ`, so the added decay is a `c*` effect acting
through chamber pressure, not a nozzle-`γ` effect.

The Milestone 3 prescribed constants turned out **optimistic**: at the reference
point the CEA table gives ideal `c*` −7.84 %, `T_c` −16.45 %, `γ` +2.67 %,
`R` +4.48 %.

## Model hierarchy: M3 → M4 → M5

| Model | Burn [s] | `F₀` [N] | `F_end` [N] | Slope | Impulse [N·s] | `I_sp` [s] |
| --- | --- | --- | --- | --- | --- | --- |
| M3 prescribed flow | 41.741 | 310.05 | 333.53 | **+7.57 %** | 13 525.65 | 227.10 |
| M4 coupled, frozen properties | 42.917 | 313.73 | 286.76 | **−8.60 %** | 13 018.56 | 225.73 |
| M5 coupled, CEA chemistry | 42.522 | 290.76 | 253.21 | **−12.91 %** | 11 627.03 | 199.45 |

Decomposed into the two physical causes:

| Effect | Impulse | Burn time | Slope [pp] |
| --- | --- | --- | --- |
| **A.** feed/blowdown coupling (M4 − M3) | −3.75 % | +2.82 % | **−16.17** |
| **B.** variable thermochemistry (M5 − M4) | **−10.69 %** | −0.92 % | −4.32 |
| **C.** total coupled effect (M5 − M3) | −14.04 % | +1.87 % | −20.49 |

The two steps are not interchangeable and neither alone accounts for the net
change. Feed coupling is what flips the slope while barely moving impulse; the
chemistry does the opposite. Slope is quoted in percentage points because it is
itself a percentage.

These three models are **not of equal fidelity**. Each removes an assumption the
previous made, and none is validated against hardware.

## Sensitivity ranking

Deterministic, one factor at a time, ranked by max |% change from baseline| across
the tested alternatives. **Not** a probabilistic analysis.

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

**The dominant assumption is not the same for every output** — burn duration is
led by the regression coefficient while impulse and slope are led by the injector
model — so no single universal ranking is claimed.

Three results worth stating plainly:

* **The injector model is the largest single lever**, ahead of every propellant
  property. Which reduced-order injector correlation you pick matters more to this
  prediction than the chemistry does.
* **`η_c*` barely moves total impulse (0.07 %)** despite being the largest
  unvalidated number in the model. The pressure-fed system self-compensates: a
  higher `c*` raises `p_c`, which throttles the injector, while the grain-limited
  burn fixes the fuel mass. It still moves chamber pressure ~7 %.
* **The thermochemistry table's own resolution is irrelevant** at the committed
  fidelity — under 0.04 % on every metric. The M5 chemistry matters through *what*
  it says about `c*(O/F)`, not through how finely it is tabulated.

One perturbation was **invalid and is reported, not repaired**: the bare HEM
injector limit gives an oxidizer flow too low for the chamber to reach the table's
10 bar floor, so no firing state is established. It is excluded from the rankings
and was not clamped back into validity.

## Robust conclusions

Each tested across the deterministic family. "Robust within tested deterministic
family" means every case in *this* set agreed — it is not a probability and says
nothing about untested regions.

| Conclusion | Support |
| --- | --- |
| Coupling the tank reverses the M3 rising-thrust trend | 5/5 |
| Variable thermochemistry preserves the falling-thrust sign | 15/15 |
| Variable thermochemistry deepens the decay relative to M4 | 5/5 |
| `c*` variation dominates the M5 chemistry effect over `γ` | 15/15 |
| Feed coupling changes thrust shape more than the chemistry does | 5/5 |
| `O/F` drifts downward during the burn | 15/15 |
| The 10 L reference case burns through before liquid runs out | 15/15 |
| A small tank can instead terminate on liquid depletion | 1/1 |

No exceptions were found in any of the eight.

## Model-sensitive conclusions

What is **not** robust is the *level* of every performance number:

* **Total impulse** moves ±13.6 % across the injector-model choice alone, and
  ±9 % across a ±10 % change in the regression coefficient. Any absolute impulse
  or `I_sp` figure from this model should be read as order-of-magnitude.
* **Thrust decay magnitude** ranges from −2.5 % to −17.8 % across the tested set,
  even though its *sign* never changes. The sign is a robust conclusion; the
  number is not.
* **Which terminal event occurs** depends on tank sizing — 10 L burns the grain
  through, 2.5 L runs out of liquid first.
* **Absolute `I_sp`** carries the full weight of the unmeasured `η_c*` plus an
  ideal, loss-free nozzle. It is not a performance prediction for any real motor.

## Verification

59 independent checks, all passing. The expected side of each is built from
longhand algebra, closed forms, alternate numerical routes, published measurements
or raw CSV columns — never by calling the same function twice.

| Kind | Checks | Worst residual |
| --- | --- | --- |
| Algebraic identity | 49 | 9.313e-10 Pa (≈ 4e-16 relative) |
| Discretisation | 8 | 2.911e-07 relative |
| Empirical (vs published firings) | 2 | 1.085e-01 |

Discretisation residuals are integrator-tolerance quantities and are reported
separately rather than folded into one headline figure.

Test suite: **857 passing**, zero warnings under `-W error`, `ruff check` clean.
All 29 figures regenerate byte-identically. Both committed CEA data files
regenerate byte-identically.

One genuine defect was found and fixed this milestone: the Milestone 5 coupled
solver reported the wrong *reason* when a perturbed case fell off the
thermochemistry table. Diagnostic only — no published number changed, verified by
byte-identical study output. See [VERIFICATION.md](VERIFICATION.md).

## Limitations

The model is a coupled reduced-order study, not a design tool.

* **Generic conceptual geometry.** The grain, tank and nozzle are illustrative
  and taken from no real motor.
* **Effective injector area only** — no hole count, hole size or plate geometry,
  and `C_d` = 0.66 is a literature average for *other* injectors.
* **Reduced-order injector correlations** (SPI / HEM / Dyer-NHNE), which the
  sensitivity study shows to be the single largest lever on the answer.
* **Equilibrium tank model**, with no vapour-only discharge tail after liquid
  depletion, no feed-line dynamics and no valve dynamics.
* **Equilibrium chemistry, from a frozen table.** Equilibrium assumes reaction is
  instantaneous and complete; there is no finite-rate chemistry and no mixing
  model. Outside the tabulated rectangle the model refuses to run rather than
  extrapolating.
* **One `η_c*` for every combustion loss**, unmeasured.
* **Ideal 1-D nozzle** with no losses, no contour, no shocks or separation
  modelling, and a single chamber `γ` through the expansion.
* **No ignition or chamber-filling transient**, no combustion instability.
* **No structural or thermal design**, no hardware validation, no flight
  prediction.

Being more coupled is not being more validated. Each milestone replaced one
assumption with several new ones, and none of them has been checked against a
firing.
