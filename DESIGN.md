# DESIGN — Milestones 1, 2 and 3

* **Part I (§1–§11)** — Milestone 1: instantaneous prescribed-flow regression bookkeeping.
* **Part II (§12–§23)** — Milestone 2: transient prescribed-flow port evolution.
* **Part III (§24–§33)** — Milestone 3: quasi-steady chamber/nozzle performance and thrust.

Parts I and II model no chamber pressure, nozzle flow or thrust. Part III adds
all three, but the oxidizer flow remains **prescribed** in every part: chamber
pressure never feeds back into it, and no feed-system closure exists anywhere.

---

# Part I — Milestone 1

Reduced-order model of a **generic** N₂O / HTPB hybrid rocket motor.

> **Scope statement.** This repository is an educational, reduced-order engineering
> study. It is **not** a motor design, a fabrication guide, a test procedure, or a
> flight-hardware model. The geometry is illustrative and is not taken from any
> real motor. Nothing here should be used to build, load, or fire hardware.

---

## 1. Conventions

### 1.1 Units

**SI internally, everywhere, without exception.**

| Quantity | Symbol | Internal unit |
| --- | --- | --- |
| Length, diameter | `L`, `D_p`, `D_o` | m |
| Area | `A_port`, `A_burn` | m² |
| Volume | `V_f` | m³ |
| Mass | `m_f` | kg |
| Density | `ρ_f` | kg/m³ |
| Oxidizer mass flow | `ṁ_ox` | kg/s |
| Fuel mass flow | `ṁ_f` | kg/s |
| Oxidizer mass flux | `G_ox` | kg/(m²·s) |
| Regression rate | `ṙ` | m/s |
| Mixture ratio | `O/F` | – (dimensionless) |
| Time | `t` | s (Milestone 2) |

Non-SI units appear in exactly two places, both deliberate and both explicit:

1. **Correlation input.** A `PowerLawRegressionLaw` stores its coefficient in the
   units of its source publication and declares those units through the
   `MassFluxUnit` / `RegressionRateUnit` enumerations. Conversion to SI happens
   **once**, at construction.
2. **Human-readable output.** Regression rates are printed in mm/s and fluxes are
   additionally printed in g/(cm²·s) so that results can be compared against the
   literature by eye. These are presentation-layer conversions only.

`mm/s` and `m/s` are never implicitly mixed; neither are `g/(cm²·s)` and
`kg/(m²·s)`. Every public evaluation function is SI-in, SI-out.

### 1.2 Naming

Function and attribute names carry their unit as a suffix (`_m`, `_m2`, `_kg_s`,
`_si`, `_mm_s`). Coefficients that are **not** sourced from the literature are
named and labelled `ILLUSTRATIVE` in code, docstrings, printed output and figures.

### 1.3 Immutability

`GrainGeometry`, `PowerLawRegressionLaw` and `OperatingPoint` are frozen
dataclasses. Geometry carries no time state: the *current* port diameter is
passed to each query, so a single geometry can be evaluated at any hypothetical
port size without mutation. This keeps Milestone 1 explicitly quasi-static.

---

## 2. Geometry equations

A single, central, circular port through a cylindrical grain. The port is assumed
to stay circular, concentric and axially uniform; grain end faces are assumed
inhibited, so only the port wall burns.

```
A_port(D_p)   = π D_p² / 4                       [m²]   port flow area
P_port(D_p)   = π D_p                            [m]    burning perimeter
A_burn(D_p)   = P_port · L = π D_p L             [m²]   burning surface area
A_fuel(D_p)   = π (D_o² − D_p²) / 4              [m²]   remaining fuel cross-section
V_fuel(D_p)   = A_fuel · L                       [m³]   remaining fuel volume
m_fuel(D_p)   = ρ_f · V_fuel                     [kg]   remaining fuel mass
web₀          = (D_o − D_p,0) / 2                [m]    initial radial web
```

### 2.1 Rejected inputs

`GrainGeometry` raises `ValueError` for:

* any non-finite or non-positive `L`, `D_p,0`, `D_o` or `ρ_f`;
* `D_p,0 ≥ D_o` (there would be no fuel);
* any queried port diameter that is non-finite, non-positive, or `≥ D_o`
  (the grain would be burnt through).

### 2.2 Representative geometry (illustrative)

| Parameter | Value | Note |
| --- | --- | --- |
| Grain length `L` | 0.40 m | illustrative |
| Initial port diameter `D_p,0` | 0.040 m | illustrative |
| Outer grain diameter `D_o` | 0.090 m | illustrative |
| Number of ports | 1 | circular, central |
| Fuel density `ρ_f` | 930 kg/m³ | representative HTPB-like nominal value; **not** traced to a specific measurement |
| Initial web | 0.025 m | derived |
| `A_port,0` | 1.256637 × 10⁻³ m² | derived |
| `A_burn,0` | 5.026548 × 10⁻² m² | derived |
| Loaded fuel mass | 1.8991 kg | derived |

These are chosen to be a plausible small laboratory-scale hybrid, sized so that
the study fluxes land inside the range over which the sourced correlation
reports data (§4.4). They are **not** copied from any real motor.

---

## 3. Oxidizer mass-flux definition

```
G_ox = ṁ_ox / A_port          [kg/(m²·s)]
```

By the usual hybrid convention this is an **oxidizer** flux: the port area is the
oxidizer flow area, and fuel added along the port is *not* included. The
alternative "total mass flux" `G_tot = (ṁ_ox + ṁ_f)/A_port` is **not** used at
this milestone, because the sourced correlation (§4) is itself fitted against
oxidizer flux.

For a circular port this is equivalently `G_ox = 4 ṁ_ox / (π D_p²)`, which is the
independent form used by the test suite and by `scripts/manual_check.py`.

**`ṁ_ox` is a prescribed input at Milestone 1.** It is never derived from tank
pressure, N₂O thermodynamics, injector characteristics, feed-system losses or
chamber pressure.

---

## 4. Regression source audit

### 4.1 What was consulted

| # | Source | What it was used for | Access |
| --- | --- | --- | --- |
| S1 | H. Rezaei, M. R. Soltani, A. R. Mohammadi, *Experimental study of fuel regression rate in an HTPB/N₂O hybrid rocket motor*, **Scientia Iranica B** 25(1) (2018) 253–265, DOI [10.24200/sci.2017.4317](https://doi.org/10.24200/sci.2017.4317) | The N₂O/HTPB coefficient pair actually used | Full text read |
| S2 | T. A. Marquardt, J. Majdalani, *A Primer on Classical Regression Rate Modeling in Hybrid Rockets*, **AIAA 2020-3758**, AIAA Propulsion and Energy 2020 Forum, DOI 10.2514/6.2020-3758 | Confirmation of the power-law structure and of what the theory does and does not predict | Full text read |
| S3 | M. A. Karabeyoglu, B. J. Cantwell, G. Zilliac, *Development of Scalable Space–Time Averaged Regression Rate Expressions for Hybrid Rockets*, **Journal of Propulsion and Power** 23(4) (2007) | Reported exponent range; scalability limitations of dimensional power-law fits | Full text read |

### 4.2 What the literature directly supports

* **The power-law structure is standard and well founded.** S2 reconstructs
  Marxman's diffusion-limited theory and arrives at the local law
  `ṙ = a₁ G^0.8 x^(−0.2)` (S2 Eq. 21), then states that the *space–time averaged*
  equivalent `ṙ = a₂ Gⁿ` (S2 Eq. 22) "is also often used by researchers reporting
  experimental data measurements". S3 opens from the same local law
  `ṙ = a Gⁿ x^m` (S3 Eq. 1). S1 states the same reduced form as its Eq. (9).
  **This project uses the space–time averaged form `ṙ = a G_ox ⁿ`.**
* **`G_ox = ṁ_ox / A_port`.** All three sources work in terms of the oxidizer port
  mass flux; S3 defines `A_port` as the local port area in its nomenclature, and
  S1 tabulates `G_ox` alongside `ṁ_ox` and port diameter for each firing.
* **Expected dependence on oxidizer mass flux.** S3: the classical mass-flux and
  length exponents are 0.8 and −0.2 respectively, derived for a fully turbulent
  boundary layer on a flat plate; in reality "the flux exponent is reported to be
  in the range of **0.5–0.8** for most hybrid systems". S2 concurs: for classical,
  non-metallised hybrids the regression rate theoretically depends on `G^0.8` and
  **not** on chamber pressure.
* **N₂O/HTPB-specific coefficients exist and are published.** S1 fitted 18 firings
  of an HTPB grain with N₂O and reports (its Eq. 10):

  ```
  ṙ = 0.3977 · G_o^0.3667
  ```

  S1 states the units for its correlations explicitly: `ṙ`, `G_o` and `L` are in
  **mm/s**, **g/(cm²·s)** and **mm** respectively.

### 4.3 Coefficient decision

**The primary correlation is SOURCED, not invented and not calibrated by this
project.** It is S1 Eq. (10), used exactly as published:

| Property | Value |
| --- | --- |
| Coefficient `a` (source units) | 0.3977 |
| Exponent `n` | 0.3667 |
| Source rate unit | mm/s |
| Source flux unit | g/(cm²·s) |
| Propellant combination | N₂O / HTPB — **yes, specifically validated for this combination** |
| Provenance | `Provenance.SOURCED` |

No fitting, tuning or re-calibration of `a` or `n` was performed at any point.

**In addition**, three clearly labelled `Provenance.ILLUSTRATIVE` laws are
supplied for exponent-sensitivity work (§4.6). They are never presented as
experimental N₂O/HTPB data.

### 4.4 Exact unit handling

The one conversion in the model is performed once, at construction:

```
ṙ_SI = k_r · a_src · (G_SI / k_G)^n
     = [ k_r · a_src · k_G^(−n) ] · G_SI^n
     = a_SI · G_SI^n
```

where `k_r` is metres per second per source rate unit and `k_G` is
kilograms per square metre per second per source flux unit:

```
k_r = 1×10⁻³  (m/s per mm/s)
k_G = 10      (kg/(m²·s) per g/(cm²·s), since 1e−3 kg / 1e−4 m² = 10)
```

giving

```
a_SI = 1×10⁻³ · 0.3977 · 10^(−0.3667) = 1.709446806×10⁻⁴
ṙ [m/s] = 1.709446806×10⁻⁴ · (G_ox [kg/(m²·s)])^0.3667
```

This conversion is verified three independent ways (§7).

The source's reported data range, 3.5–12 g/(cm²·s), converts to
**35–120 kg/(m²·s)**. This range is stored on the law and every study row is
flagged `in-range` or `EXTRAPOLATED` against it.

### 4.5 Limitations of simple power-law correlations

Recorded honestly, from the sources themselves:

1. **No pressure dependence.** S2: classical non-metallised hybrids theoretically
   depend on `G^0.8` but not on chamber pressure — accurate for many systems, but
   contingent on neglecting radiation.
2. **Radiation breaks it at low flux.** S2: radiation "should not be neglected at
   low mass fluxes or for metallised or other heavily sooting fuels".
3. **Kinetics break it at high flux.** S2 reviews regimes where regression becomes
   nearly independent of `G` and instead varies with pressure, citing an
   explicitly pressure-dependent correlation `ṙ ∝ P^0.5 G^0.3 x^(−0.2)`.
4. **Axial variation is smeared out.** The `x^(−0.2)` axial term of the local law
   is absorbed into a space-averaged constant. S2 notes the axial exponent is
   close to zero for some motors and that injector effects commonly dominate it.
5. **Fits do not scale between motors.** S3 warns that regression laws written in
   terms of dimensional parameters "are accurate for interpolation purposes but
   potentially highly problematic when they are used to extrapolate to other
   scales", giving a worked example in which a 10× scale-up predicts a ~5×
   regression rate — "clearly ... highly unrealistic".
6. **Data reduction is not unique.** S3: methods of reducing firing data to a
   correlation differ, and different methods can yield significantly different
   results for the same tests.

### 4.6 Known discrepancies in the sourced coefficient — reported, not tuned away

Two are material and both are stated wherever the correlation is used.

**(a) The exponent is below the consensus band.** S1's `n = 0.3667` lies *below*
the 0.5–0.8 range S3 reports for most hybrid systems and well below the classical
0.8. Consequence: in this model `ṙ` responds only weakly to flux — over the study
sweep a 3.33× rise in `ṁ_ox` lifts `ṙ` by only 1.56×. This is not corrected. To
show what a steeper exponent would imply, three **illustrative** laws are provided
with `n = 0.5`, `0.62`, `0.8`. Their coefficients are fixed by one disclosed
condition only:

```
ṙ_variant(G_anchor) = ṙ_sourced(G_anchor),   G_anchor = 70 kg/(m²·s)
```

i.e. the curve is *pivoted* about a single stated anchor flux. This is a
transparent construction, **not** a calibration: no measurement informs the
resulting coefficient beyond the anchor value inherited from the sourced law. All
such laws are marked `ILLUSTRATIVE` in their label, reference string, printed
output and figure legends.

**(b) The correlation carries an unmodelled grain-length dependence.** S1 also
reports (its Eq. 11) a length-dependent fit from further tests in which fuel
length was varied from 150 mm to 310 mm:

```
ṙ = 0.07577 · G_o^0.364 · L^0.293      [ṙ mm/s, G_o g/(cm²·s), L mm]
```

Eq. (10) was fitted at a fixed `L ≈ 250 mm`. The representative grain here is
`L = 400 mm`. Taking Eq. (11) at face value, the length ratio alone implies
`(400/250)^0.293 ≈ 1.15`, i.e. the regression rates reported by this project may be
**~15 % low** for a 400 mm grain. This project does **not** apply that correction,
because 400 mm is itself outside the 150–310 mm range over which Eq. (11) was
fitted, and because S3 warns specifically against extrapolating dimensional fits
of this kind to other scales. The effect is documented rather than silently
absorbed.

### 4.7 What was *not* verified

* The Doran et al. AIAA 2006 N₂O/HTPB paper that S1 compares against could not be
  retrieved (the hosted copy presented an expired TLS certificate). S1's claim of
  agreement with it is therefore **unverified by this project**.
* The fit itself was not re-derived from S1's raw firing data; instead the SI
  implementation was checked against S1's 18 published measurements (§7).
* The HTPB density of 930 kg/m³ is a commonly used nominal value and is **not**
  traced here to a specific measurement. It is treated as an illustrative input.
* No coefficient was validated against any hardware by this project. Nothing here
  has been validated against a real motor.

---

## 5. Fuel-flow bookkeeping

```
ṁ_f = ρ_f · A_burn · ṙ = ρ_f · (π D_p L) · ṙ        [kg/s]
```

The instantaneous rate at which solid fuel is gasified by a surface receding at
`ṙ` over area `A_burn`. Zero regression gives exactly zero fuel flow. Negative
density, area or regression rate is rejected.

Note the two competing effects when the port grows at fixed `ṁ_ox`: `A_burn`
increases (`∝ D_p`) while `ṙ` decreases (`∝ D_p^(−2n)`). With `n = 0.3667` the
net exponent on `D_p` is `1 − 2n = +0.267`, so **`ṁ_f` rises slightly as the port
opens up** even though `ṙ` falls. The study output states this explicitly so the
non-monotonic-looking result is not mistaken for a bug.

---

## 6. O/F bookkeeping

```
O/F = ṁ_ox / ṁ_f        [–]
```

Undefined when `ṁ_f ≤ 0`; the implementation raises rather than returning
`inf`/`nan`.

**This is bookkeeping only.** At Milestone 1 `O/F` is *not* used to compute flame
temperature, `c*`, chamber pressure, combustion efficiency, nozzle performance,
specific impulse or thrust. It is reported so that the mixture the motor would be
running at is visible, and so that a later milestone has a verified input.

---

## 7. Verification strategy

Tests do **not** call a production function to build their own expectation. Every
expected value comes from one of:

1. **Independently arranged algebra** — e.g. `A_port` checked against `π r²`
   rather than `π D²/4`; `A_burn` against `2 π r L`; `G_ox` against
   `4 ṁ_ox / (π D_p²)`.
2. **Hand-computed decimal literals** — checked in to 12–16 significant figures.
3. **An independent unit route** — the reference implementation of the regression
   law converts SI flux *back* into the source's g/(cm²·s), evaluates the printed
   correlation in mm/s, and converts to m/s. It never touches `coefficient_si`.
4. **Published experimental data** — the strongest external check. All 18
   (`G_ox`, `ṙ`) pairs from S1's Table 4 are transcribed into the test suite, and
   the SI code path is required to reproduce each measurement to within 15 %
   (the largest residual observed is ~10 %, consistent with the scatter of the
   authors' own fit). Any silent factor-of-10 or factor-of-1000 slip anywhere in
   the unit chain fails these outright.
5. **Structural invariants** — e.g. `ṙ(kG)/ṙ(G) = kⁿ` exactly; doubling `D_p` at
   fixed `ṁ_ox` quarters `G_ox` exactly.

Required checks, and where they live:

| # | Check | File |
| --- | --- | --- |
| 1 | Circular port area | `tests/test_geometry.py` |
| 2 | Burning perimeter | `tests/test_geometry.py` |
| 3 | Burning area | `tests/test_geometry.py` |
| 4 | Remaining fuel volume | `tests/test_geometry.py` |
| 5 | Remaining fuel mass | `tests/test_geometry.py` |
| 6 | `G_ox = ṁ_ox / A_port` | `tests/test_operating_point.py` |
| 7 | Regression power law | `tests/test_regression.py` |
| 8 | Coefficient / unit conversion | `tests/test_regression.py` |
| 9 | `ṁ_f = ρ_f A_burn ṙ` | `tests/test_operating_point.py` |
| 10 | `O/F` | `tests/test_operating_point.py` |
| 11 | Zero oxidizer-flow behaviour | `tests/test_regression.py`, `tests/test_operating_point.py` |
| 12 | Negative-input rejection | all three |
| 13 | Invalid geometry rejection | `tests/test_geometry.py`, `tests/test_operating_point.py` |
| 14 | Scalar / array consistency | `tests/test_regression.py` |
| 15 | Monotonicity with mass flux | `tests/test_regression.py` |
| 16 | `ṙ` falls as port grows at fixed `ṁ_ox` | `tests/test_operating_point.py` |

Beyond the unit tests:

* `scripts/manual_check.py` re-derives the whole chain longhand from `math`
  primitives and compares against the package to a relative tolerance of 1×10⁻¹².
  It also asserts that no chamber-pressure, nozzle, `c*`, injector or thrust
  symbol is exported by the package — a scope guard against accidental creep.
* `scripts/regression_study.py` runs an explicit engineering sanity audit and
  exits non-zero if any check fails.
* Figures are regenerated twice and compared by SHA-256 to confirm byte-identical
  determinism. This holds for a fixed matplotlib/FreeType build (verified with
  matplotlib 3.10.9 / FreeType 2.6.1); PNG glyph rasterisation depends on the
  FreeType version, so hashes are not expected to be stable across environments
  with different builds, even though the figure content is identical.

Gate before commit: `pytest -W error -q` and `ruff check .` must both pass.

---

## 8. Assumptions

1. Single, central, circular port; stays circular and concentric.
2. Port is axially uniform; the port diameter is a single scalar.
3. Grain end faces are inhibited — only the port wall burns.
4. Fuel is homogeneous with uniform density.
5. Regression rate is spatially uniform along the port (the space–time averaged
   correlation is applied as a lumped value).
6. Regression is diffusion-limited and correlates with oxidizer flux alone.
7. No chamber-pressure dependence of regression (consistent with classical
   hybrid theory, S2).
8. Quasi-static: every result is an instantaneous snapshot; there is no time
   integration and no transient.
9. Oxidizer mass flow is prescribed and perfectly known.
10. Complete gasification of the receding fuel surface — no slag, no residue, no
    unburnt fuel accounting.

---

## 9. Limitations

* Coefficients are from **one** experimental campaign on **one** motor, with the
  discrepancies of §4.6.
* The study grain is longer (400 mm) than the grain the correlation was fitted on
  (250 mm), with the ~15 % implication quantified in §4.6(b).
* Extrapolation outside 35–120 kg/(m²·s) is flagged but not prevented; the 70 mm
  port case in the sensitivity study is genuinely outside the source data.
* No uncertainty quantification: the model returns point values, and no
  confidence interval on `a` or `n` is propagated.
* Nothing in this repository has been validated against hardware.

---

## 10. Deferred physics — explicitly NOT modelled at Milestone 1

Milestone 1 does **not** model, compute, approximate or report any of:

* N₂O tank thermodynamics (self-pressurisation, two-phase blowdown, saturation)
* Injector flow, discharge coefficient, cavitation, or injector sizing
* Feed-system pressure drop
* Chamber pressure
* Combustion efficiency
* Equilibrium chemistry / flame temperature
* Characteristic velocity `c*`
* Nozzle flow, expansion ratio, throat sizing, throat erosion
* Thrust, specific impulse, total impulse
* Structural design (case, bulkheads, joints, safety factors)
* Thermal design (liner, insulation, soak-back)
* Ignition
* Fabrication, assembly, handling or test procedures

`scripts/manual_check.py` enforces the code half of this list programmatically.

---

## 11. Milestone 1 status

Milestone 1 is complete and frozen. Its production physics, tests, scripts,
figures, coefficients and reported values are unchanged by Milestone 2. Part II
below builds on them without editing them; the only Milestone 1 files touched
are the package `__init__.py` (re-exports), `pyproject.toml` (a `scipy`
dependency and a version bump) and this document.

---

# Part II — Milestone 2

Transient port regression under prescribed oxidizer-flow histories.

---

## 12. Source audit for the geometric state equation

Milestone 2 introduces exactly one new physical statement: how the port geometry
evolves given a regression rate. No new regression coefficients were sought or
introduced; Part I's sourced N₂O/HTPB law is reused verbatim.

### 12.1 What was consulted

The same three sources audited in §4.1 were re-read for this purpose:

| # | Source | Used for |
| --- | --- | --- |
| S2 | Marquardt & Majdalani, **AIAA 2020-3758** | The definition of `ṙ` as a surface recession velocity |
| S3 | Karabeyoglu, Cantwell & Zilliac, **JPP** 23(4) (2007) | The port-diameter state equation, the quasi-steady assumption, the prescribed-flow boundary condition, an exact closed-form solution for `n = 0.5`, and a global mass balance as a numerical check |

### 12.2 What the literature directly supports

* **`ṙ` is a surface velocity, not a mass rate.** S2's Eq. (1) writes the wall
  heat flux as `Q̇_w = ṁ''_f h_v = ρ_f ṙ h_v`, i.e. `ρ_f ṙ` *is* the fuel mass
  flux leaving the surface. So `ṙ` has dimensions of length/time and is the speed
  at which the solid–gas interface recedes along its normal. This is what makes
  `ṁ_f = ρ_f A_burn ṙ` (§5) and the state equation below the *same* statement
  viewed two ways.
* **For a circular port, `2ṙ = ∂D/∂t`.** S3 states this explicitly, in those
  terms, immediately before its port-diameter equation: *"Based on the definition
  of regression rate, `2ṙ = ∂D/∂t`, the variation of the port diameter at an
  axial position x and time t can be written as"* — their Eq. (51),
  `∂D/∂t = C_D ṁⁿ x^m / D^(2n)`. With no length dependence (`m = 0`) this is
  precisely the equation integrated here. S3 derives it for *"the case of a
  single circular port"*, the same geometry as this project.
* **The regression law is applied quasi-steadily.** S3 introduces its companion
  port mass balance *"under the quasi-steady assumption"*. The empirical
  correlation is a steady-state result, so at each instant the model evaluates it
  at the instantaneous flux and treats the flow field as having already adjusted.
* **Oxidizer flow enters as an externally prescribed boundary condition.** S3
  solves its system with the initial condition `D(x, t=0) = D_i(x)` and the
  boundary condition `ṁ = ṁ_o(t)` — the oxidizer mass flow is *given*, not
  derived from tank or chamber state. That is exactly the framing used here.
* **A closed-form solution exists for `n = 0.5`.** S3 §VI.C derives an exact
  solution for a flux exponent of 0.5; with no length dependence and constant
  oxidizer mass flow it takes the form `D(t) = [D_i² + C t]^(1/2)`, i.e. the
  squared diameter is linear in time.
* **A global mass balance is the recommended numerical check.** S3's Eq. (63)
  states that the fuel mass added must equal the increase in port volume times
  the solid density, and says it *"can be used to check the accuracy of the
  numerical solutions"*. §17 below implements exactly that.

### 12.3 What was *not* verified

* No new source was consulted for regression coefficients, by design.
* S3's treatment is a partial differential equation in `(x, t)`; this project
  integrates only the `m = 0`, axially lumped ordinary differential equation. The
  axial structure of S3's solution is therefore **not** reproduced or verified
  here, and the sliver/non-uniformity behaviour it predicts is outside this model
  (see §22).

---

## 13. Derivation of the state equation

The fuel surface recedes along its outward normal at speed `ṙ`. For a circular
port the normal is radial and the surface is the cylinder `r = r_p`, so

```
dr_p/dt = ṙ                                                      (13.1)
```

and since `D_p = 2 r_p`,

```
dD_p/dt = 2 ṙ                                                    (13.2)
```

which is S3's stated identity. Everything else follows from Part I unchanged:

```
A_port(t) = π r_p(t)²
G_ox(t)   = ṁ_ox(t) / A_port(t)
ṙ(t)      = a_SI · G_ox(t)ⁿ
A_burn(t) = 2 π r_p(t) L
ṁ_f(t)    = ρ_f · A_burn(t) · ṙ(t)
O/F(t)    = ṁ_ox(t) / ṁ_f(t)
```

### 13.1 Assumptions specific to the transient model

1. The port stays circular, concentric and axially uniform at all times; `r_p` is
   a single scalar.
2. Regression is spatially uniform — the space–time averaged correlation is
   applied as a lumped instantaneous value (see §22 for what this costs).
3. The correlation is applied **quasi-steadily**: no thermal lag in the solid, no
   transient boundary-layer development, no ignition or shutdown dynamics.
4. Oxidizer mass flow is prescribed and perfectly known at every instant.
5. Grain end faces remain inhibited, so `A_burn = 2π r_p L` throughout.
6. Fuel is homogeneous, so `ρ_f` is constant in time and space.
7. There is no coupling back from a chamber, nozzle or feed system — by design.

---

## 14. Closed-form constant-flow solution

Substituting `G_ox = ṁ_ox/(π r²)` into `ṙ = a_SI G_ox ⁿ` gives

```
dr/dt = a_SI (ṁ_ox/π)ⁿ r^(−2n) = K r^(−2n),      K = a_SI (ṁ_ox/π)ⁿ    (14.1)
```

`K` carries SI units throughout because `a_SI` is the once-converted coefficient
from §4.4; no mm/s or g/(cm²·s) quantity ever enters. Separating,

```
r^(2n) dr = K dt   ⟹   r^(2n+1)/(2n+1) = K t + C
```

and applying `r(0) = r_0`,

```
r(t) = [ r_0^(2n+1) + (2n+1) K t ]^(1/(2n+1))                     (14.2)
```

### 14.1 Analytic burnout time

Setting `r = r_outer` and inverting:

```
t_burn = [ r_outer^(2n+1) − r_0^(2n+1) ] / [ (2n+1) K ]           (14.3)
```

For `ṁ_ox = 0` we get `K = 0` and the implementation returns `inf`: the port
never regresses, so burnout is never reached.

### 14.2 Structural check against the published exact solution

For `n = 0.5` the exponent `2n+1 = 2`, so (14.2) becomes `r² = r_0² + 2Kt` and
therefore `D² = D_i² + 8Kt` — the squared diameter is linear in time, matching
the form of S3's exact `n = 0.5` solution (§12.2). The test suite asserts both
the linearity and the slope `8K`.

### 14.3 Numerical values at the representative point

With `a_SI = 1.709446806 × 10⁻⁴`, `n = 0.3667`, `ṁ_ox = 0.100 kg/s`:

```
K       = 1.709446806e-4 × (0.100/π)^0.3667 = 4.828928 × 10⁻⁵  (SI)
2n + 1  = 1.7334
t_burn  = (0.045^1.7334 − 0.020^1.7334) / (1.7334 × 4.828928e-5) = 41.740618 s
```

Differentiating (14.2) at `t = 0` returns `K r_0^(−2n) = 8.508938 × 10⁻⁴ m/s`
= 0.850894 mm/s, identical to the Milestone 1 headline value — the transient
model starts exactly where the instantaneous model sits.

---

## 15. Numerical solver design

### 15.1 Augmented state

Cumulative masses are integrated as ODE states rather than post-processed, so
they inherit the integrator's own error control:

```
y = [ r_p , m_ox,cum , m_f,cum ]

dy/dt = [ ṙ , ṁ_ox(t) , ρ_f (2π r_p L) ṙ ]
```

### 15.2 Method

`scipy.integrate.solve_ivp` with `DOP853` by default (`RK45` is accepted and
tested), `rtol = 1e-10` / `atol = 1e-14` by default, and `dense_output=True` so
the reporting grid can be chosen independently of the steps actually taken.
`scipy` is therefore a new runtime dependency at Milestone 2; it is declared in
`pyproject.toml`.

### 15.3 Segment-wise integration

The integration is **restarted at every interior breakpoint** of the flow
history. A piecewise-constant schedule is therefore integrated exactly, segment
by segment, instead of having a step discontinuity smeared across an adaptive
step. Within a segment the prescribed flow is continuous, so the integrator sees
a smooth problem.

### 15.4 A defect found and fixed during development

The first implementation looked the flow up globally with half-open segments
`[start, end)`. An integrator legitimately evaluates the right-hand side *at* a
segment's right endpoint, and the global lookup there returned the **next**
segment's flow. The consequence was measurable and physical, not cosmetic: the
final step of each segment was integrated with the wrong derivative, leaving a
1.3 × 10⁻¹⁰ kg error in the accumulated oxidizer mass at the first breakpoint,
and — worse — a zero-flow interval was no longer exactly frozen, drifting by
3.6 × 10⁻¹¹ kg. The tests for "a zero-flow interval leaves the radius unchanged"
caught it.

The fix is `OxidizerFlowHistory.flow_on_segment(start, end)`, which returns a
flow law valid on the *closed* segment. `PiecewiseConstantOxidizerFlow` overrides
it to bind the single constant identified by the segment midpoint, so the law is
exactly constant across the whole segment including both endpoints. After the
fix, zero-flow intervals are frozen to *exactly* zero drift (spread 0.0 in
radius, cumulative oxidizer mass and cumulative fuel mass), and the oxidizer mass
at the breakpoint is exact to round-off.

### 15.5 Reporting grid and the right-limit convention

The reporting grid is `n_report` evenly spaced samples, unioned with the
termination time and with every interior breakpoint. At an interior breakpoint
the grid carries **two** nearly coincident samples — `nextafter(t_b, −∞)` and
`t_b` — so a discontinuous prescribed flow is reported by its left *and* right
limits and is never drawn as a spurious ramp. Derived quantities at `t_b` itself
use the right limit, matching the half-open `[start, end)` convention of the
schedule.

The grid affects sampling only: refining it from 11 to 4001 points changes
neither the burnout time nor the consumed fuel (§19).

---

## 16. Event handling

A terminal event `g(t, y) = r_p − r_outer` with `direction = +1` stops the
integration the instant the port reaches the outer grain radius. SciPy locates
the root by Brent's method on the dense output, so the burnout time is resolved
to solver precision rather than to a step boundary.

Consequences, all asserted by tests:

* the solution is **never continued into negative fuel thickness**;
* `r_p ≤ r_outer` and `D_p ≤ D_o` at every reported sample (the final sample's
  round-off is clipped to the boundary, a sub-nanometre correction);
* remaining fuel mass at burnout is zero to within 10⁻⁹ kg;
* the reported status is `FUEL_DEPLETED`, and the requested horizon is reported
  as *not* reached.

A run that reaches the end of its window with fuel remaining reports
`COMPLETED_TIME_WINDOW` and `burnout_time_s = None`. Nothing beyond burnout is
fabricated.

Milestone 1's `GrainGeometry.fuel_mass_kg` deliberately rejects `D_p ≥ D_o`,
which is correct for instantaneous bookkeeping but excludes the burnout point
itself. Rather than modify frozen Milestone 1 code, `transient.remaining_fuel_mass_kg`
evaluates the same geometry on the **closed** interval. Being independently
written, it doubles as a cross-check: a test asserts it agrees with the
Milestone 1 routine everywhere on the open interior.

---

## 17. Mass-conservation equations

Two independent bookkeeping paths must agree.

**Geometric fuel loss** — from the port radius alone:

```
m_f,consumed(t) = ρ_f L π ( r_p(t)² − r_p,0² )
m_f,remaining(t) = ρ_f L π ( r_outer² − r_p(t)² )
```

**Integrated fuel flow** — from the ODE accumulator:

```
m_f,cum(t) = ∫₀ᵗ ṁ_f dt = ∫₀ᵗ ρ_f (2π r_p L) ṙ dt
```

These are equal *identically*, because

```
ṁ_f = ρ_f (2π r_p L) dr_p/dt = ρ_f L π d(r_p²)/dt
```

so the residual `m_f,cum − ρ_f L π (r_p² − r_p,0²)` is **zero in exact
arithmetic**. It is therefore a pure measure of integration error and a sharp
check on the solver — not a modelling approximation with a physical tolerance.
This is the global mass balance S3 recommends (§12.2).

**Oxidizer bookkeeping** is checked the same way:

```
m_ox,used(t) = ∫₀ᵗ ṁ_ox dt
```

compared for constant flow against `ṁ_ox · t`, and for a schedule against a
hand-written sum of `flow × duration` over segments.

Measured residuals, worst case across all four study cases: **8.3 × 10⁻¹¹ kg
absolute, 7.9 × 10⁻¹¹ relative**, on a 1.899 kg grain. Constant-flow oxidizer
closure is below 10⁻¹¹ kg. Residuals are reported, never hidden.

---

## 18. Prescribed-flow schedule representation

Three history types, all validated at construction:

| Type | Meaning |
| --- | --- |
| `ConstantOxidizerFlow` | one constant flow over a window |
| `PiecewiseConstantOxidizerFlow` | a contiguous, gapless, strictly increasing sequence of `FlowSegment`s; `from_durations` builds one from `(duration, flow)` pairs |
| `CallableOxidizerFlow` | a user-supplied `ṁ_ox(t)`, sampled at construction and re-validated on every RHS evaluation |

Rejected inputs (all raise `ValueError`, or `TypeError` for a non-callable):

* negative or non-finite mass flow;
* negative time, or a window with `t_end ≤ t_start`;
* segments with non-positive duration;
* schedules with gaps, overlaps or out-of-order segments;
* an empty schedule;
* a flow query outside the schedule's window;
* an initial port diameter that is non-positive or not inside the grain
  (delegated to the Milestone 1 geometry);
* `n_report < 2`, or a `t_end_s` outside the flow history's window.

A callable is integrated as a single smooth segment; a history with genuine
discontinuities must use `PiecewiseConstantOxidizerFlow` so the solver restarts
at each jump. This is documented on the class.

### 18.1 Zero-flow convention

At `ṁ_ox = 0`:

```
G_ox = 0     ṙ = 0     ṁ_f = 0     dr_p/dt = 0     O/F = NaN
```

`ṙ = 0` is the exact continuous limit of `a G^n` as `G → 0` for `n > 0`, so no
special case is needed in the physics. The port radius, cumulative oxidizer mass
and cumulative fuel mass are all exactly frozen (verified to a spread of 0.0).

`O/F` is **undefined** when `ṁ_f = 0` — there is no fuel to form a ratio with —
so the history reports `NaN`. This is a deliberate choice over `inf` or a
sentinel: `NaN` propagates loudly, is trivially testable with `isnan`, and causes
plotting libraries to *break the line* rather than draw a fabricated value.
Figures M2-C and M2-E rely on exactly that: nothing is drawn across the coast.
The scalar helper `operating_point.mixture_ratio` (Milestone 1) still *raises*
for zero fuel flow; that difference is intentional — a scalar caller has made a
single meaningless request, whereas a history legitimately contains samples where
the quantity does not exist.

---

## 19. Validity and extrapolation classification

Milestone 1 recorded the source's measured flux range as 35–120 kg/(m²·s). Every
reported transient sample is classified into `FluxRangeStatus`:

| Status | Meaning |
| --- | --- |
| `WITHIN_SOURCE_FLUX_RANGE` | `G_ox` inside the measured range |
| `BELOW_SOURCE_FLUX_RANGE` | extrapolating below the data |
| `ABOVE_SOURCE_FLUX_RANGE` | extrapolating above the data |
| `ZERO_OXIDIZER_FLOW` | `ṁ_ox = 0`; the exact zero limit, not an extrapolation |
| `NO_SOURCE_RANGE_DECLARED` | the law carries no validated range (e.g. an illustrative law) |

`ZERO_OXIDIZER_FLOW` is kept distinct from `BELOW_SOURCE_FLUX_RANGE`
deliberately: at zero flow the model is not extrapolating an empirical
correlation at all, and lumping the two together would overstate how much of a
burn is extrapolated.

`TransientResult.time_fraction_with_status` weights by elapsed time
(trapezoidally over the reporting grid), not by sample count, so an uneven grid
cannot distort the answer.

**Nothing is clamped and nothing is discarded.** `G_ox` and `ṙ` are evaluated and
reported wherever the trajectory goes; the classification is metadata alongside
them. The study prints an explicit warning when a burnout prediction rests
mostly on extrapolation, and figures M2-B draws extrapolated segments dotted.

### 19.1 Result

| Case | `ṁ_ox` [kg/s] | Within range | Below range |
| --- | ---: | ---: | ---: |
| A | 0.100 | 33.7 % | 66.3 % |
| B | 0.050 | 3.8 % | 96.2 % |
| C | 0.150 | 61.6 % | 38.4 % |
| D | piecewise | 50.0 % | 35.0 % (plus 15.0 % zero-flow) |

Because port area grows at constant flow, `G_ox` falls monotonically and the
trajectory always ends up below the data. **The burnout predictions in Cases A
and B rest mostly on extrapolation**, and the study says so in those words.
Counter-intuitively, *lower* flow leaves the validated range sooner, not later.

---

## 20. Convergence study

The closed-form solution (§14) is the reference, so this is convergence against
an exact answer rather than against a finer numerical run.

| `rtol` | Method | `t_burn` [s] | Relative error | `r(10 s)` rel. error | Mass-closure residual [kg] |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1e-05 | DOP853 | 41.740617467 | 8.05e-09 | 9.41e-10 | 3.05e-07 |
| 1e-07 | DOP853 | 41.740617771 | 7.71e-10 | 1.30e-10 | 1.43e-08 |
| 1e-09 | DOP853 | 41.740617800 | 6.25e-11 | 2.83e-12 | 4.59e-10 |
| 1e-11 | DOP853 | 41.740617803 | 9.15e-13 | 4.81e-14 | 1.09e-11 |
| 1e-11 | RK45 | 41.740617803 | 6.11e-12 | 3.43e-13 | 3.92e-10 |

Closed-form reference: `t_burn = 41.740617803 s`, `r(10 s) = 0.027506697529 m`.

Every reported quantity converges monotonically towards the analytic value, and
two independent integrators agree. Reporting-grid refinement (11 → 4001 samples)
leaves the burnout time and consumed fuel unchanged to all printed digits, and a
`max_step = 0.25 s` restriction changes nothing. The results are not artefacts of
the tolerance or the grid.

---

## 21. Milestone 2 verification strategy

The independence rules of §7 continue to apply. Specifically, for Milestone 2:

* the closed-form solution is **re-derived inside the test file** from the
  coefficient as printed in the source paper (`0.3977`, mm/s, g/(cm²·s)) with its
  own unit arithmetic, and never read back from the production law object;
* piecewise trajectories are checked against the closed form **chained segment by
  segment by hand** — something the production solver never does, since it
  integrates numerically;
* the `n = 0.5` case is checked against the *published* exact solution's
  structure (`D²` linear in `t`, slope `8K`);
* mass closure is checked against geometry formed independently of the
  integrator's accumulator states;
* piecewise oxidizer mass is checked against a hand-written `Σ flow × duration`;
* the transient fuel-mass helper is cross-checked against Milestone 1's geometry
  routine on the open interior;
* several histories are re-checked sample-by-sample against Milestone 1's
  `OperatingPoint` at the corresponding port diameter.

`tests/test_transient.py` contains 119 tests covering all 30 required checks:
initial-state agreement with Milestone 1; `dr/dt = ṙ` and `dD/dt = 2ṙ` by central
differences; the analytic radius and burnout time; numerical-versus-analytic
agreement; monotonicity of radius, diameter, area and remaining fuel;
non-negative fuel; no overshoot; falling `G_ox` and `ṙ` at constant flow; the
fuel-flow and O/F equations; fuel and oxidizer mass conservation; piecewise mass
integration; all four zero-flow behaviours including the `NaN` convention; the
burnout event; a short non-burnout run; every rejection case; scalar/array
consistency; flux classification including the zero-flow category; and
convergence in tolerance, integrator, grid and `max_step`.

Whole-repository gate before commit: `pytest -W error -q` (286 tests) and
`ruff check .` must both pass.

---

## 22. Milestone 2 limitations

1. **Most of a constant-flow burn is extrapolated** (§19). The burnout times are
   extrapolated results. This is the single most important caveat on every
   Milestone 2 number.
2. **Spatially uniform regression.** The port is one scalar radius, so axial
   variation of `ṙ` is ignored. S3 shows that real single-port grains regress
   non-uniformly and leave a residual *sliver* of unburnt fuel at the point where
   the thinnest section breaks through. This model instead consumes the grain
   exactly and uniformly, so its burnout is an **idealised upper bound** on fuel
   utilisation and its burn time is correspondingly optimistic.
3. **Quasi-steady regression.** A steady-state correlation is applied
   instantaneously. There is no thermal lag in the solid, no transient
   boundary-layer response, no ignition transient and no tail-off.
4. **The grain-length dependence is still unmodelled.** The §4.6(b) caveat — that
   predictions may be ~15 % low for this 400 mm grain — applies to every
   transient result here, including the burnout times.
5. **Circular, concentric port for all time.** No coning, no non-circular ports,
   no multi-port geometry, no erosion of the port shape.
6. **No coupling back from the chamber.** The oxidizer flow is prescribed, so
   nothing represents how a real feed system and chamber would respond to the
   growing port. A real motor's oxidizer flow would itself depend on chamber
   pressure.
7. **No uncertainty quantification.** Point values only; no confidence interval
   on `a`, `n`, or the resulting burnout time is propagated.
8. **Nothing here has been validated against hardware.**

---

## 23. Deferred physics — still NOT modelled after Milestone 2

Everything listed in §10 remains out of scope and absent from production code:
N₂O tank thermodynamics, injector flow, feed-system pressure drop, chamber
pressure, combustion efficiency, equilibrium chemistry, `c*`, nozzle flow,
thrust, specific impulse, structural design, thermal design, ignition and
fabrication.

`scripts/manual_check.py` enforces the code half of this list programmatically,
and the Milestone 2 study prints its own scope notice.

**Milestone 3 (not started)** would couple the fuel and oxidizer flows computed
here to combustion properties, chamber pressure and nozzle flow, so that thrust
could eventually be predicted. No part of that exists in this repository.

---

# Part III — Milestone 3

Quasi-steady chamber/nozzle coupling and thrust prediction, still under a
prescribed oxidizer mass flow.

---

## 24. Source audit for the chamber and nozzle equations

### 24.1 What was consulted

| # | Source | Used for | Access |
| --- | --- | --- | --- |
| S4 | NASA Glenn Research Center, *Isentropic Flow Relations* (`isentrop.html`) | Static-to-total pressure and temperature ratios (their Eqs. 6 and 7) | Fetched and read |
| S5 | NASA Glenn Research Center, *Area Ratio as a Function of Mach Number* (`astar.html`) | The area-Mach relation | Fetched and read |
| S6 | NASA Glenn Research Center, *Mass Flow Rate at Choking Conditions* (`mflchk.html`) | Choked mass-flow relation and `V = M a = M sqrt(g R T)` | Fetched and read |
| S7 | NASA Glenn Research Center, *Rocket Thrust Equation* (`rockth.html`) | `F = m_dot V_e + (p_e - p_0) A_e` | Fetched and read |
| S8 | NASA Glenn Research Center, *Specific Impulse* (`specimp.html`) | `I_sp = F/(m_dot g_0)` | Fetched and read |
| S9 | J. Seitzman, *Thrust Coefficient, Characteristic Velocity and Ideal Nozzle Expansion*, Georgia Tech AE4451 | `c* = p_0 A_t / m_dot` (their Eq. IV.14/IV.16) and `C_F = F/(p_0 A_t)` | Fetched; equations partially recoverable from the PDF |
| S10 | M. C. Tarifa, L. Pizzuti, *Theoretical performance analysis of hybrid rocket propellants…*, EUCASS 2019-488, DOI 10.13009/EUCASS2019-488 | Ideal `c*` as a function of `k`, `T` and molar mass (their Eq. 1); adiabatic flame temperatures below 3000 K for N₂O-based hybrids | Fetched and read |
| S11 | R. Stark, *Flow Separation in Rocket Nozzles, a Simple Criteria*, AIAA 2005-3940 (DLR) | Summerfield's `p_sep/p_a ≈ 0.4` separation rule of thumb | Criterion confirmed via search summary |
| S1 | Rezaei, Soltani & Mohammadi, Scientia Iranica B 25(1) (2018) — the Milestone 1 regression source | Measured `c*` (1403–1587 m/s), chamber pressures (19.9–31.0 bar), thrust (15–29 kgf) and `I_sp` (167–223 s) for a comparable HTPB/N₂O motor, and the 94–98 % combustion efficiency band | Full text read in Milestone 1 |

### 24.2 Equation-by-equation verification

**A. Characteristic velocity.** S9 defines `c* = p_0 A_t / m_dot`. The same result
follows independently from S6: dividing `p_c A_t` by the choked mass flow

```
m_dot = (A p_t / sqrt(T_t)) sqrt(g/R) [(g+1)/2]^(-(g+1)/(2(g-1)))
```

gives `c* = sqrt(R T_c/g) [(g+1)/2]^((g+1)/(2(g-1)))`. S10's Eq. (1) states the
same dependence on `k`, `T` and molar mass. The test suite verifies both routes
agree.

**B. Quasi-steady chamber balance.** `m_dot_ox + m_dot_f = m_dot_nozzle`
rearranged with A gives `p_c = (m_dot_ox + m_dot_f) c* / A_t`. This is algebra on
the definition, not a separate physical claim.

**C, D. Area-Mach relation.** S5 prints it as

```
A/A* = {[(g+1)/2]^-k} / M * [1 + M^2 (g-1)/2]^k ,   k = (g+1)/(2(g-1))
```

which is algebraically identical to the `(2/(g+1))` form used in this project,
since `[(g+1)/2]^(-k) = (2/(g+1))^k`. The test suite checks the two arrangements
agree to 1e-12 across five values of `gamma`.

**E, F. Isentropic ratios.** S4 Eq. (6): `p/p_t = [1 + M^2 (g-1)/2]^(-g/(g-1))`.
S4 Eq. (7): `T/T_t = [1 + M^2 (g-1)/2]^(-1)`. Chamber conditions are treated as
stagnation conditions, the standard low-Mach-chamber assumption.

**G. Exit velocity.** S6: `V = M a = M sqrt(g R T)`.

**H. Thrust.** S7: `F = m_dot * Ve + (pe - p0) * Ae`, with S7 noting explicitly
that the pressure term is needed whenever the exit pressure differs from ambient.

**I. Thrust coefficient.** S9: `C_F = F/(p_0 A_t)`, so `F = C_F p_c A_t`.

**J. Specific impulse.** S8: `Isp = F / (mdot * g0)`. S8's schoolroom page rounds
`g_0` to 9.8 m/s²; this project uses the defined standard value
**9.80665 m/s²**.

### 24.3 What was *not* verified

* S9's slide deck extracts poorly from PDF (symbol fonts), so its equations were
  read structurally rather than quoted verbatim; the `c*` definition is
  independently corroborated by the S6 route above, which is the one the code and
  tests actually rely on.
* S11's separation criterion was confirmed from a search summary of the DLR/AIAA
  literature rather than from the primary PDF. It is used **only** as a warning
  flag, never in any calculation.
* **No N₂O/HTPB thermochemistry was computed or validated.** `gamma`, `T_c` and
  the molar mass are prescribed assumptions (§25), not results.

---

## 25. Prescribed combustion properties

### 25.1 Why they are prescribed

Milestone 3 deliberately contains **no equilibrium-chemistry solver and no CEA
interface**. Computing `c*`, `gamma` and `T_c` from first principles is a
separate, large piece of work; asserting values without computing them would be
worse. So they are declared as explicit, clearly labelled assumptions, and their
influence is quantified (§30).

### 25.2 The selection rule, fixed in advance

Four primitives were chosen **before any thrust was computed** and were not
adjusted afterwards:

| Primitive | Value | Basis |
| --- | ---: | --- |
| `gamma` | 1.20 | round value typical of hot rocket exhaust |
| `T_c` | 2600 K | round value below the sub-3000 K adiabatic flame temperatures S10 reports for N₂O-based hybrids, reduced because this model runs fuel-rich at O/F ≈ 2–2.5 |
| Molar mass `M` | 0.022 kg/mol | round value |
| `eta_c*` | 0.96 | midpoint of the 94–98 % combustion efficiency S1 measured for this propellant pair |

From these, `R = R_u/M = 377.930 J/(kg K)`,
`c*_ideal = 1528.486 m/s` and the delivered `c* = 0.96 x c*_ideal = 1467.347 m/s`.

### 25.3 Consistency, not validation

The delivered `c*` lands inside the 1403–1587 m/s that S1 measured, and the
implied chamber pressure (26.11 bar) lands inside the 19.9–31.0 bar S1 measured.
**Both are weak consistency checks, not validations.** The primitives are round
numbers and the efficiency is a stated midpoint; nothing was back-solved. Anyone
can re-derive the chain and confirm no fitting occurred.

`CombustionProperties` exposes `ideal_c_star_m_s` and `c_star_efficiency`
precisely so that this internal consistency is visible. A prescribed set implying
`eta_c* > 1` is physically impossible, and §30 shows exactly that happening when
`c*` is varied alone — which is the point of showing it.

---

## 26. Reference nozzle selection

Documented, and deliberately not optimised.

The throat was chosen by evaluating round candidates at the Milestone 1 total
mass flow, against the chamber-pressure band S1 measured for a comparable motor:

| `D_t` [mm] | `A_t` [m²] | `p_c` [bar] | Verdict |
| ---: | ---: | ---: | --- |
| 8 | 5.026548e-05 | 40.80 | above the 19.9–31.0 bar band |
| **10** | **7.853982e-05** | **26.11** | **inside the band — chosen** |
| 12 | 1.130973e-04 | 18.13 | below the band |
| 15 | 1.767146e-04 | 11.61 | well below the band |

A round 10 mm throat is the only round candidate landing mid-band. The expansion
ratio was set to a round `epsilon = 4`. Ideal sea-level expansion for this chamber
pressure would require `epsilon ≈ 4.34`, so the reference nozzle is *slightly*
underexpanded at sea level — a consequence of choosing a round number, not a
design decision. No thrust value influenced either choice.

---

## 27. Numerical method

### 27.1 Exit-Mach root solve

The area-Mach relation has a subsonic and a supersonic root for every
`A/A* > 1`. A choked converging-diverging nozzle running full sits on the
supersonic branch, where `A/A*` is monotonically increasing in `M`. The solver
therefore:

1. brackets the root on `[1, M_hi]`, doubling `M_hi` until the residual changes
   sign (capped at Mach 1e4);
2. solves with `scipy.optimize.brentq`, which is guaranteed to converge on a
   sign-changing bracket.

An unbracketed Newton iteration was deliberately **not** used: it can converge
onto the subsonic branch or diverge at large expansion ratios, and there is no
reason to accept that risk when a bracket is available for free. `epsilon = 1`
short-circuits to `M_e = 1` exactly.

Measured behaviour: the area-Mach residual falls to `0` at `xtol = 1e-15`, and
thrust is unchanged to 6 decimal places from `xtol = 1e-4` downwards (§31).

### 27.2 Transient coupling architecture

The coupling is **one-way and instantaneous**. Each Milestone 2 sample supplies
`m_dot_ox` and `m_dot_f`; the chamber and nozzle are solved at that state; nothing
flows back into the regression or port model. Milestone 2's `TransientResult` is
carried on the `ThrustHistory` unchanged, so mass bookkeeping remains the
Milestone 2 quantity and is not recomputed.

Because `epsilon` and `gamma` are constant, `M_e` and both isentropic ratios are
constant; because `T_c` is also constant, `T_e` and `V_e` are constant for an
entire burn. Only `p_e` moves, in proportion to `p_c`.

---

## 28. Idle / zero-flow behaviour

With `m_dot_ox + m_dot_f = 0` the chamber is reported as `IDLE_NO_FLOW` with
`p_c = 0`, and:

* thrust is exactly **0** — *not* `-p_a A_e`. Applying the pressure term to a
  nozzle merely full of still ambient gas would be meaningless; there is no jet.
* **no exit state is reported at all** (`exit_state is None`, and every exit array
  in a history is `NaN`), so no nozzle number can be mistaken for a firing
  condition.
* `C_F` and `I_sp` are `NaN`, because both are literally `0/0`.
* the mixture ratio is `NaN` whenever there is no fuel flow, matching the
  Milestone 2 convention (§18.1).

---

## 29. Thrust decomposition and integrated performance

```
momentum thrust = m_dot_total V_e
pressure thrust = (p_e - p_a) A_e
F               = momentum + pressure          (exactly, as floating-point addition)
C_F             = F / (p_c A_t)
I_sp            = F / (m_dot_total g_0),   g_0 = 9.80665 m/s^2
I_total         = integral F dt              (trapezoidal over the reporting grid)
I_sp,eq         = I_total / (m_propellant g_0)
mean thrust     = I_total / (active burn time)
```

The equivalent specific impulse is preferred over averaging instantaneous `I_sp`
because it has a clean integral definition and remains well defined when some
samples are idle. The active burn time is the trapezoidal integral of the firing
indicator, so a coast is excluded from the mean.

### 29.1 An exact structural identity

Since `p_c ∝ m_dot_total` while `V_e` and `p_e/p_c` are fixed,

```
F = m_dot_total [ V_e + c* epsilon (p_e/p_c) ] - p_a A_e
```

Thrust is exactly **affine** in the total mass flow. This is used as a third
independent verification route (§32) and explains why the thrust history is so
nearly a scaled copy of the mass-flow history.

---

## 30. Sensitivity methodology

Two deterministic, one-factor-at-a-time studies. Neither is uncertainty
propagation: no distributions are assumed and no intervals are produced.

**Ambient pressure** is swept at a fixed chamber state. The chamber and the entire
exit state are unchanged — the nozzle is choked — so only `(p_e - p_a) A_e` moves.
Sea level to vacuum gains exactly `p_a A_e = 31.832 N`.

**Combustion properties** are varied one at a time with mass flow and nozzle
geometry held fixed: `c*` by ±10 %, `gamma` over 1.15–1.25, `T_c` by ±10 %. The
implied `eta_c*` is reported for every variation, which makes visible that varying
one property alone renders the set internally inconsistent (the `c* +10 %` case
implies `eta_c* = 1.056`, physically impossible). That is a deliberate disclosure,
not an oversight.

Because `m_dot` is held fixed, the percentage change in `I_sp` is *identical* to
that in thrust in every row, so Figure M3-E plots the chamber-pressure response
instead — which genuinely differs.

---

## 31. Numerical convergence

Root solve (`epsilon = 4`, `gamma = 1.20`):

| `xtol` | `M_e` | Area-Mach residual | `F` [N] |
| ---: | ---: | ---: | ---: |
| 1e-04 | 2.619446781721 | 2.683e-08 | 310.047929 |
| 1e-08 | 2.619446776666 | -1.442e-12 | 310.047929 |
| 1e-12 | 2.619446776666 | 1.220e-12 | 310.047929 |
| 1e-15 | 2.619446776666 | 0.000e+00 | 310.047929 |

Transient integration and integrated performance:

| `rtol` | `n_report` | `t_burn` [s] | Peak `F` [N] | `I_total` [N s] | `I_sp,eq` [s] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1e-05 | 201 | 41.740617467 | 333.528161 | 13525.649028 | 227.103122 |
| 1e-07 | 401 | 41.740617771 | 333.528161 | 13525.651206 | 227.103159 |
| 1e-09 | 801 | 41.740617800 | 333.528161 | 13525.651738 | 227.103168 |
| 1e-11 | 1601 | 41.740617803 | 333.528161 | 13525.651869 | 227.103170 |

Total impulse converges to about 1e-7 relative and the burnout time to the
Milestone 2 closed-form value. Peak thrust is unchanged to 6 decimals throughout,
because it is an algebraic function of the burnout radius rather than an integral.

---

## 32. Milestone 3 verification strategy

`tests/test_nozzle.py` (113 tests), `tests/test_chamber.py` (51) and
`tests/test_performance.py` (51) cover all 40 required checks. Independence
comes from:

* **A different algebraic arrangement.** The area ratio is rebuilt in the tests in
  the `[(g+1)/2]^(-k)/M * [...]^k` form printed by S5, while the production code
  uses the `(2/(g+1))^k` form.
* **A plain-bisection area-Mach inversion** written in the test file, using no
  scipy and no production routine.
* **A full longhand reimplementation** of the chamber → nozzle → thrust chain from
  `math` primitives, compared field by field against the production result.
* **Three independent thrust routes** — exit-state, `C_F p_c A_t`, and the affine
  identity of §29.1 — none of which is used to generate another.
* **Bypassing the ODE integrator**: transient chamber pressures and thrusts are
  checked against the Milestone 2 *closed-form* radius solution re-derived in the
  test file.
* **Exact identities**: `p_c A_t - m_dot c*` and `F - (momentum + pressure)` must
  close to machine zero, not to a tolerance.
* **An independent route to ideal `c*`**: built from the S6 choked mass-flow
  relation and compared with the closed form.

Measured residuals at the reference point: `|A - B| = 0` N, `|A - C| = 5.7e-14` N,
`p_c A_t - m_dot c* = 0` N exactly.

Gate before commit: `pytest -W error -q` (501 tests) and `ruff check .`.

---

## 33. Milestone 3 limitations

1. **Oxidizer flow is prescribed and chamber pressure does not feed back into it.**
   The largest departure from a real motor. A real feed system would deliver less
   flow as `p_c` rises, damping the very thrust rise this model predicts. Every
   Milestone 3 thrust history should be read with that in mind.
2. **`c*`, `gamma` and `T_c` are prescribed constants** that do not vary with O/F,
   pressure or time. No equilibrium chemistry, no CEA.
3. **Combustion efficiency is represented solely by the prescribed `c*`.** No other
   loss mechanism exists anywhere in the model, which is why the predicted `I_sp`
   (226.2 s) sits above the 167–223 s S1 measured.
4. **The nozzle is 1-D, steady, adiabatic and isentropic** with a calorically
   perfect gas: no divergence loss, no boundary layer, no heat transfer, no
   two-phase flow, no shocks, no separation model, no side loads.
5. **The Summerfield flag warns but does not correct.** A flagged result is
   reported as doubtful, not fixed.
6. **No ignition, chamber-filling or tail-off transient.** The chamber has no gas
   storage term, so prescribed-flow steps produce instantaneous pressure and
   thrust steps that no real motor could follow.
7. **No combustion instability, structural sizing, thermal sizing, nozzle contour,
   materials or fabrication information.**
8. **Nothing has been validated against hardware**, and nothing here establishes a
   safe operating pressure or a flight-ready design.

---

## 34. Deferred physics — still NOT modelled after Milestone 3

Absent from production code: N₂O tank state and thermodynamics, saturation and
vapour-pressure models, tank blowdown, vapour-liquid equilibrium, injector CdA
sizing or mass-flow equations, feed-line pressure loss, chamber-to-feed coupling,
equilibrium chemistry / CEA, combustion instability, ignition transients, nozzle
contour generation, wall heat transfer, structural or casing stress analysis,
thermal sizing, and flight dynamics.

### 34.1 One documented change to a prior-milestone script

`scripts/manual_check.py` (Milestone 1) contained a scope guard asserting that no
chamber-pressure, nozzle or thrust symbol was exported by the package. That was
the correct boundary for Milestones 1 and 2, and it fired as designed the moment
Milestone 3 added `NozzleGeometry` and `ThrustHistory`.

This is **not** a defect in Milestone 1 — the guard did exactly its job — so it
was not "fixed"; it was **retargeted** at the physics that is still deferred
(§34). The change is confined to the `forbidden` tuple and its printed label, and
is annotated in place. Every numeric check in that script is untouched and every
Milestone 1 value it reports is unchanged:

```
A_port 1.256637061e-03   P_port 1.256637061e-01   A_burn 5.026548246e-02
m_f    1.899092759e+00   G_ox   7.957747155e+01   r_dot  8.508937503e-04
m_f    3.977664394e-02   O/F    2.514038141e+00
```

The boundary is additionally pinned by two permanent tests in
`tests/test_performance.py`: `test_package_exports_no_deferred_physics`, which
fails if any deferred-physics token reaches the public API, and
`test_oxidizer_flow_is_never_derived_from_chamber_pressure`, which fails if the
one-way coupling is ever reversed.

**Milestone 4 (not started)** would couple chamber pressure back to an N₂O
feed/injector/tank model, so that the oxidizer flow is no longer externally
prescribed. That needs properly sourced N₂O property data and two-phase /
feed-system modelling, with strong scope controls. No part of it exists here.
