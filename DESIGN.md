# DESIGN — Milestone 1

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

## 11. Next milestone

**Milestone 2 — transient port-radius evolution under prescribed oxidizer-flow
histories:** integrate `dD_p/dt = 2 ṙ` for a prescribed `ṁ_ox(t)`, handle fuel
depletion / burnout events when `D_p → D_o`, and produce time histories of `ṙ`,
`ṁ_f` and `O/F`. Still no chamber-pressure or nozzle coupling.
