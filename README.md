# Hybrid Rocket Motor — N₂O/HTPB Regression Rate Study

A rigorous, reduced-order engineering model of a **generic** nitrous-oxide /
HTPB hybrid rocket motor.

* **Milestone 1** — instantaneous prescribed-flow regression bookkeeping.
* **Milestone 2** — transient prescribed-flow port evolution, burnout and mass closure.
* **Milestone 3** — quasi-steady chamber pressure, 1-D nozzle and **thrust**.
* **Milestone 4** — coupled N₂O tank, injector and motor: **blowdown**.
* **Milestone 5** — O/F-dependent **equilibrium thermochemistry** from NASA CEA.

Milestone 5 removes the last prescribed combustion assumption. `c*`, `γ`, `T_c`
and the molar mass are interpolated from a frozen NASA CEA equilibrium table at
the *instantaneous* mixture ratio and chamber pressure, and fed back into the
coupled solve — so the chemistry and the feed system now move each other. Nothing
has been validated against hardware.

> **The Milestone 4 headline result.** Holding the oxidizer flow constant
> (Milestone 3) predicts thrust **rising** +7.6 % through the burn. With the tank
> coupled (Milestone 4) it **falls** −8.6 %. The constant-flow assumption gets the
> *sign* of the thrust-time slope wrong for this configuration.

> **The Milestone 5 headline result.** Letting the chemistry vary does **not**
> overturn that sign — it deepens it. The same burn now falls **−12.9 %**, about
> 1.5× the frozen-property decay, and the whole thrust curve drops with it
> (−10.7 % total impulse). The mechanism is new: as the port opens the mixture
> ratio falls, and on the fuel-rich side of the CEA table a falling O/F drags `c*`
> down with it. In Milestone 4 that channel did not exist, because `c*` was a
> constant. The tank and the chemistry now decay together.

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

## Milestone 2 scope

Milestone 2 makes the model time-dependent. The port radius becomes a state
variable integrated forward under a **prescribed** oxidizer mass-flow history:

* the geometric state equation `dr_p/dt = ṙ`, equivalently `dD_p/dt = 2ṙ`,
* constant, piecewise-constant and callable prescribed flow histories,
* a terminal fuel-depletion (burnout) event at `r_p = r_outer`,
* time histories of `G_ox`, `ṙ`, `ṁ_f`, `O/F`, remaining fuel and cumulative masses,
* a closed-form constant-flow solution used as an independent verification path,
* fuel- and oxidizer-mass conservation checks,
* per-sample tracking of whether the regression law is being **extrapolated**
  outside the flux range its source data covers.

**The oxidizer flow is still prescribed.** Milestone 2 adds no feed-system,
injector, tank or chamber physics, and still computes no chamber pressure, `c*`,
nozzle flow or thrust.

## Milestone 3 scope

Milestone 3 makes the model thrust-capable. Each instantaneous propellant-flow
state is coupled to a quasi-steady chamber mass balance and a 1-D isentropic
converging–diverging nozzle:

* chamber pressure from a **prescribed** characteristic velocity,
* supersonic exit Mach from the area–Mach relation (bracketed Brent solve),
* exit pressure, temperature and velocity from the isentropic relations,
* thrust, thrust coefficient and specific impulse,
* expansion-regime classification and a Summerfield separation flag,
* deterministic ambient-pressure and combustion-property sensitivity studies.

**Still prescribed, still not validated.** Milestone 3 adds **no** N₂O tank
thermodynamics, vapour–liquid equilibrium, blowdown, injector sizing or
pressure-drop, feed-line losses, ignition or chamber-filling transient,
combustion instability, equilibrium chemistry, nozzle contour, shock or
separation model, or structural/thermal sizing.

## Milestone 4 scope

Milestone 4 makes the oxidizer flow an **output**. Three new physical models are
solved together with the frozen Milestone 1–3 motor:

* a rigid, adiabatic, **saturated-equilibrium N₂O tank** whose pressure is
  `p_sat(T)` and whose temperature falls as liquid is withdrawn,
* an **injector** with three sourced models — incompressible SPI, homogeneous
  equilibrium HEM, and the Dyer/NHNE blend,
* a **coupled feed/chamber closure** solved as a bracketed scalar root,

integrated in time with terminal events for grain burnout, liquid depletion,
tank/chamber pressure equalisation and departure from the checked property band.

**Not modelled, deliberately:** equilibrium combustion chemistry, combustion
instability, ignition or chamber-filling transients, feed-line pressure losses,
nozzle contours, shocks or flow separation, structural or thermal sizing, and
the vapour-only discharge tail after liquid depletion.

### N₂O properties

Properties come from **CoolProp**, whose N₂O equation of state is the short
fundamental Helmholtz EOS of Lemmon & Span, *J. Chem. Eng. Data* 51 (2006)
785–850 — the same EOS the NIST Chemistry WebBook cites. No property correlation
is fitted, invented or transcribed anywhere in this project. The test suite
checks the module against an independently retrieved NIST WebBook table over
200–305 K; because NIST cites the same EOS, that verifies **units and call
structure**, not the EOS itself, and the tests say so.

### Tank model

A rigid, adiabatic, spatially uniform control volume in saturated equilibrium.
The open-system first law with `Q = 0`, `W = 0` and a saturated-liquid outflow:

```
dm/dt = −ṁ_out
dU/dt = −ṁ_out · h_l(T)

V = m_l/ρ_l(T) + m_v/ρ_v(T)     m = m_l + m_v
U = m_l·u_l(T) + m_v·u_v(T)     p = p_sat(T)
```

Given the conserved `(m, U)` the temperature follows from one bracketed root
solve. **Nothing prescribes a pressure-versus-time curve** — the pressure history
is an output of the energy balance.

### Injector models

All three are taken verbatim from Waxman, Zimmerman, Cantwell (Stanford) and
Zilliac (NASA Ames), NTRS 20190001326:

```
SPI   ṁ = C_d A √(2 ρ Δp)                                    their Eq. (2)
HEM   ṁ = C_d A ρ₂ √(2(h₁ − h₂)),  s₂ = s₁                   their Eqs. (3)–(4)
κ     = √((p₁ − p₂)/(p_v − p₂))                              their Eq. (8)
Dyer  ṁ = κ/(1+κ)·ṁ_SPI + 1/(1+κ)·ṁ_HEM                      their Eq. (9)
```

For a *saturated* tank `p₁ = p_v`, so **κ is exactly 1** and the Dyer model is an
equal blend of the two limits. Waxman et al. warn that the SPI equation "is often
applied inappropriately" to nitrous oxide because the liquid flashes in the
orifice; both models are provided and compared.

### Coupled closure

```
residual(ṁ_ox) = injector(tank, p_c(ṁ_ox, r_p)) − ṁ_ox = 0
```

solved by Brent's method on the closed-form bracket `[0, ṁ_SPI(p₂=0)]`. A
fixed-point iteration is deliberately not used.

### Chamber mass balance

Characteristic velocity is *defined* by `c* = p_c A_t / ṁ`. Combining it with a
quasi-steady chamber balance `ṁ_ox + ṁ_f = ṁ_nozzle` gives

```
p_c = (ṁ_ox + ṁ_f) · c* / A_t
```

with negligible chamber gas storage, no ignition transient, no finite-rate
combustion, and combustion efficiency absorbed into the prescribed `c*`.

### Nozzle relations

```
A_e/A_t = (1/M_e) [ (2/(γ+1)) (1 + (γ−1) M_e²/2) ]^((γ+1)/(2(γ−1)))
p_e/p_c = [1 + (γ−1) M_e²/2]^(−γ/(γ−1))
T_e/T_c = [1 + (γ−1) M_e²/2]^(−1)
V_e     = M_e √(γ R T_e)
F       = ṁ V_e + (p_e − p_a) A_e
C_F     = F / (p_c A_t)
I_sp    = F / (ṁ g₀)
```

### Prescribed combustion properties — ILLUSTRATIVE, not chemistry

| Property | Value | Status |
| --- | ---: | --- |
| `γ` | 1.20 | round value, prescribed |
| `T_c` | 2600 K | round value, prescribed |
| Molar mass `M` | 22.0 g/mol | round value, prescribed |
| `R` | 377.930 J/(kg·K) | derived, `R_u/M` |
| Ideal `c*` | 1528.486 m/s | derived from `γ`, `R`, `T_c` |
| `c*` efficiency | 0.96 | midpoint of the 94–98 % measured by Rezaei et al. (2018) |
| **Delivered `c*`** | **1467.347 m/s** | **derived from the four rows above** |

These four primitives were fixed **before** any thrust was computed and were not
adjusted afterwards. The delivered `c*` is a *consequence* of them, not a target.
That it lands inside the 1403–1587 m/s the same study measured for HTPB/N₂O is a
weak consistency check, **not** a validation. There is no equilibrium-chemistry
solver and no CEA interface anywhere in this repository.

The model exposes `ideal_c_star_m_s` and `c_star_efficiency` so that the internal
consistency of any prescribed set is visible rather than hidden.

### O/F-dependent combustion properties — SOURCED equilibrium chemistry

Milestone 5 replaces the four constants above with a lookup at the instantaneous
operating point:

```
c*_ideal, γ, T_c, M  =  CEA_equilibrium(O/F, p_c)      bilinear interpolation
c*_delivered         =  η_c* · c*_ideal                 η_c* = 0.96, applied once
R                    =  R_u / M                         recomputed, not interpolated
```

Provenance is **SOURCED**, not illustrative: the numbers come from NASA CEA
(Gordon & McBride, RP-1311) run through the `rocketcea` wrapper, with CEA's own
`N2O` and `HTPB` reactant cards quoted verbatim in the committed CSV header. No
HTPB formula was invented and no table was hand-transcribed. Sourced still does
**not** mean validated against a firing.

Two details carry most of the numerical weight:

* **`M, (1/n)`, not `MW`.** Where equilibrium produces condensed carbon, CEA
  prints a gas-phase mean molecular weight *and* a total-mass-per-mole-of-gas.
  A single-gas nozzle model needs the latter — the one that reproduces CEA's own
  sound speed, `M = γ R_u T / a²`. At the reference point the effective molar mass
  is 21.057 g/mol against a gas-phase 19.543 g/mol, 7.75 % higher; taking the
  gas-phase number instead would inflate `R`, `c*` and thrust. See
  [DESIGN.md](DESIGN.md) §50.
* **`R` is recomputed from `M`,** not interpolated alongside it, so the identity
  `R = R_u / M` holds exactly between nodes and not merely at them.

At the Milestone 1/3 reference flows (O/F = 2.514) the sourced properties are
markedly less favourable than the illustrative constants they replace:

| Property | M3 prescribed | M5 CEA table | Change |
| --- | --- | --- | --- |
| ideal `c*` [m/s] | 1528.49 | 1408.64 | **−7.84 %** |
| `γ` [-] | 1.20000 | 1.23201 | +2.67 % |
| `T_c` [K] | 2600.0 | 2172.4 | **−16.45 %** |
| molar mass [g/mol] | 22.000 | 21.057 | −4.29 % |
| `R` [J/(kg K)] | 377.93 | 394.86 | +4.48 % |

The Milestone 3 constants were **optimistic**. They were round illustrative
values, chosen before any chemistry was run and never tuned to the table; the
table has likewise not been retuned to reproduce them.

### Reference nozzle — illustrative, documented, not optimised

| Parameter | Value |
| --- | ---: |
| Throat diameter `D_t` | 10.000 mm |
| Exit diameter `D_e` | 20.000 mm |
| Throat area `A_t` | 7.853982 × 10⁻⁵ m² |
| Exit area `A_e` | 3.141593 × 10⁻⁴ m² |
| Expansion ratio `ε` | 4.000 |
| Reference ambient `p_a` | 101325 Pa |

Selection was by inspection of round candidates against the chamber pressure the
Milestone 1 flow implies, with the 19.9–31.0 bar measured by Rezaei et al. as the
plausibility band:

| `D_t` [mm] | `p_c` [bar] | Verdict |
| ---: | ---: | --- |
| 8 | 40.80 | outside the band |
| **10** | **26.11** | **inside — chosen** |
| 12 | 18.13 | outside the band |
| 15 | 11.61 | outside the band |

`ε = 4` is a round, modest choice. Ideal sea-level expansion for this chamber
pressure would need `ε ≈ 4.34`, so the reference nozzle is **slightly
underexpanded** at sea level. No thrust target influenced either choice.

### Transient state equation

```
dr_p/dt      = ṙ  =  a · G_ox(t)ⁿ          G_ox(t) = ṁ_ox(t) / (π r_p(t)²)
dm_ox,cum/dt = ṁ_ox(t)
dm_f,cum/dt  = ρ_f · (2π r_p L) · ṙ
```

integrated from `r_p(0) = D_p,0/2` and terminated at `r_p = D_o/2`.

### Closed-form constant-flow solution

Substituting `G_ox = ṁ_ox/(π r²)` gives `dr/dt = K r^(−2n)` with
`K = a_SI (ṁ_ox/π)ⁿ`, which separates exactly:

```
r(t)   = [ r_0^(2n+1) + (2n+1) K t ]^(1/(2n+1))
t_burn = [ r_outer^(2n+1) − r_0^(2n+1) ] / [ (2n+1) K ]
```

For `n = 0.5` this collapses to `D²` linear in `t`, matching the exact solution
published by Karabeyoglu, Cantwell & Zilliac (2007) for that special case. The
closed form is the primary reference the numerical solver is verified against —
the tests never use the solver to generate their own expected answer.

## Milestone 5 scope

Milestone 5 makes the combustion properties an **output** of the mixture ratio
rather than a prescribed constant:

* a **frozen NASA CEA equilibrium table** (Gordon & McBride, NASA RP-1311, run
  through the `rocketcea` wrapper) covering O/F 1.20–4.00 and 10–45 bar,
  committed to the repository so no chemistry solver is needed at run time,
* **bilinear interpolation** — deliberately not a spline, because the products
  step sharply at the soot boundary and a spline would overshoot there,
* an **implicit chamber closure**, because `c*` now depends on the pressure it
  helps set: `p_c = ṁ c*(O/F, p_c) / A_t`,
* a **nested coupled solve** — an outer bracketed root on the feed residual, an
  inner bracketed root on chamber pressure — with no fixed-point iteration
  anywhere,
* **explicit refusal to extrapolate**: outside the tabulated rectangle the solver
  fails with a named status and the integration stops. The mixture ratio is never
  clipped into the table and the chamber pressure is never clamped to an edge.

The Milestone 3 nozzle is reused **verbatim**; only the properties handed to it
have changed. Milestone 4 is untouched and is re-run as the controlled comparison.

**Not modelled, deliberately:** finite-rate combustion chemistry and kinetics
(equilibrium buys instantaneous, complete reaction — it says nothing about
*rates*), combustion instability, ignition or chamber-filling transients,
feed-line pressure losses, nozzle contours, shocks or flow separation, structural
or thermal sizing, and the vapour-only discharge tail after liquid depletion.

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

Milestone 2 closes the loop back on geometry — and only on geometry:

```
ṙ  ────────────────►  dr_p/dt = ṙ  ────►  r_p(t)  ────►  back into A_port, A_burn
```

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

## Milestone 1 representative results

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

## Milestone 2 representative results

Case A — constant prescribed `ṁ_ox = 0.100 kg/s`, integrated to burnout:

| Quantity | At `t = 0` | At burnout |
| --- | ---: | ---: |
| Port diameter `D_p` | 40.000 mm | 90.000 mm |
| Remaining web | 25.000 mm | 0.000 mm |
| `G_ox` | 79.577 kg/(m²·s) | 15.719 kg/(m²·s) |
| `ṙ` | 0.850894 mm/s | 0.469445 mm/s |
| `A_burn` | 0.050265 m² | 0.113097 m² |
| `ṁ_f` | 0.039777 kg/s | 0.049376 kg/s |
| `O/F` | 2.514 | 2.025 |
| Remaining fuel | 1.899093 kg | 0.000000 kg |

**Burnout time 41.740618 s**, against a closed-form value of 41.740618 s —
a relative difference of 3.9 × 10⁻¹². Oxidizer consumed 4.174 kg; fuel consumed
1.899093 kg, i.e. the whole grain.

Burnout times across the three constant-flow cases:

| Case | `ṁ_ox` [kg/s] | Initial `G_ox` [kg/(m²·s)] | `t_burn` numerical [s] | `t_burn` closed form [s] | Inside source flux range |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 0.100 | 79.58 | 41.740618 | 41.740618 | **33.7 %** |
| B | 0.050 | 39.79 | 53.819861 | 53.819861 | **3.8 %** |
| C | 0.150 | 119.37 | 35.973770 | 35.973770 | **61.6 %** |

### Burnout handling

Integration is stopped by a terminal event the moment `r_p` reaches the outer
grain radius. The solution is never continued into negative fuel thickness, the
port diameter never exceeds `D_o`, and the remaining fuel mass at burnout is zero
to within 10⁻⁹ kg. If the requested horizon extends past burnout the run
terminates early and reports `FUEL_DEPLETED`; a run that reaches the end of its
window without burning through reports `COMPLETED_TIME_WINDOW`.

### Extrapolation tracking — the headline honest finding

At constant oxidizer flow the port area grows, so `G_ox` falls monotonically and
eventually drops below the 35–120 kg/(m²·s) range over which the sourced
correlation has data. Every reported sample is classified, and **nothing is
clamped or discarded**:

> Only **33.7 %** of the representative burn (Case A) lies inside the measured
> flux range. For the lower-flow Case B it is **3.8 %**. The predicted burnout
> times therefore rest mostly on extrapolation *below* the data.

Counter-intuitively, lower oxidizer flow is **not** the safer choice for model
validity here: Case B starts at `G_ox = 39.8 kg/(m²·s)`, barely inside the band,
and leaves it after about 2 s. The highest-flow case stays inside longest.

A zero-flow interval is classified separately from "below range", because at
`ṁ_ox = 0` the model returns the exact `ṙ = 0` limit rather than extrapolating a
correlation.

### Mass conservation

Because `ṁ_f = ρ_f (2π r_p L) dr_p/dt = ρ_f L π d(r_p²)/dt`, the integrated fuel
mass and the fuel mass inferred from the port geometry must agree *exactly* in
exact arithmetic. Any difference is pure integration error, which makes it a
sharp check on the solver:

| Check | Worst residual |
| --- | ---: |
| `∫ṁ_f dt` vs `ρ_f L π (r_p² − r_p,0²)` | **8.3 × 10⁻¹¹ kg** absolute, 7.9 × 10⁻¹¹ relative, on a 1.899 kg grain |
| `∫ṁ_ox dt` vs `ṁ_ox · t` (constant flow) | below 10⁻¹¹ kg |

Convergence against the closed-form burnout time, driven purely by the solver
tolerance:

| `rtol` | Method | `t_burn` [s] | Relative error | Mass-closure residual [kg] |
| ---: | --- | ---: | ---: | ---: |
| 1 × 10⁻⁵ | DOP853 | 41.740617467 | 8.1 × 10⁻⁹ | 3.0 × 10⁻⁷ |
| 1 × 10⁻⁷ | DOP853 | 41.740617771 | 7.7 × 10⁻¹⁰ | 1.4 × 10⁻⁸ |
| 1 × 10⁻⁹ | DOP853 | 41.740617800 | 6.3 × 10⁻¹¹ | 4.6 × 10⁻¹⁰ |
| 1 × 10⁻¹¹ | DOP853 | 41.740617803 | 9.2 × 10⁻¹³ | 1.1 × 10⁻¹¹ |
| 1 × 10⁻¹¹ | RK45 | 41.740617803 | 6.1 × 10⁻¹² | 3.9 × 10⁻¹⁰ |

Refining the reporting grid from 11 to 4001 samples changes neither the burnout
time nor the consumed fuel: the grid only samples the dense solution.

### A result worth stating plainly

Over Case A the regression rate **falls** from 0.8509 to 0.4694 mm/s, yet the
fuel mass flow **rises** by a factor 1.24. The burning area grows as `r_p` while
`ṙ` falls only as `r_p^(−2n)` with `n = 0.3667`, so the net exponent on `r_p` is
`1 − 2n = +0.267`. O/F therefore drifts *downward*, from 2.514 to 2.025 — an
excursion of only about 19 % across an entire burn. That mildness is a direct
consequence of the low sourced flux exponent flagged in Milestone 1.

## Milestone 3 representative results

### Reference point — the Milestone 1 operating point, now with thrust

Prescribed `ṁ_ox = 0.100 kg/s`, `ṁ_f = 0.039776644 kg/s`, `p_a = 101325 Pa`:

| Quantity | Value |
| --- | ---: |
| `ṁ_total` | 0.139777 kg/s |
| `O/F` | 2.514 |
| **Chamber pressure `p_c`** | **2 611 424.85 Pa = 26.114 bar** |
| **Exit Mach `M_e`** | **2.619447** |
| **Exit pressure `p_e`** | **113 631.85 Pa** (`p_e/p_c` = 0.043513) |
| **Exit temperature `T_e`** | **1541.974 K** |
| Exit speed of sound `a_e` | 836.248 m/s |
| **Exit velocity `V_e`** | **2190.506 m/s** |
| **Momentum thrust** | **306.182 N** |
| **Pressure thrust** | **+3.866 N** |
| **Total thrust `F`** | **310.048 N** |
| **Thrust coefficient `C_F`** | **1.511685** |
| **Specific impulse `I_sp`** | **226.190 s** |
| Expansion regime | UNDEREXPANDED (attached flow expected) |

Three independent routes to the same thrust:

| Route | Expression | Result |
| --- | --- | ---: |
| A | `ṁ V_e + (p_e − p_a) A_e` | 310.047929470 N |
| B | `C_F p_c A_t` | 310.047929470 N |
| C | `ṁ [V_e + c* ε p_e/p_c] − p_a A_e` | 310.047929470 N |

Residuals: |A − B| = 0, |A − C| = 5.7 × 10⁻¹⁴ N. The defining identity
`p_c A_t − ṁ c*` closes to exactly 0.

### Transient thrust

Case A — constant prescribed `ṁ_ox = 0.100 kg/s`, integrated to burnout at
41.740618 s:

| Quantity | At `t = 0` | At burnout |
| --- | ---: | ---: |
| `ṁ_total` | 0.139777 kg/s | 0.149376 kg/s |
| `p_c` | 26.114 bar | 27.908 bar |
| Momentum thrust | 306.182 N | 327.210 N |
| Pressure thrust | +3.866 N | +6.318 N |
| **Thrust** | **310.048 N** | **333.528 N** |
| `C_F` | 1.511685 | 1.521660 |
| `I_sp` | 226.190 s | 227.683 s |

| Integrated | Value |
| --- | ---: |
| Peak thrust | 333.528 N |
| Mean thrust (active burn) | 324.041 N |
| **Total impulse** | **13 525.65 N·s** |
| Propellant consumed | 6.073155 kg |
| **Equivalent `I_sp`** | **227.103 s** |

Across the four cases:

| Case | `ṁ_ox` [kg/s] | Burn [s] | `p_c` [bar] | Thrust [N] | Total impulse [N·s] | `I_sp,eq` [s] | Regime at sea level |
| --- | ---: | ---: | --- | --- | ---: | ---: | --- |
| A | 0.100 | 41.741 | 26.11 → 27.91 | 310.0 → 333.5 | 13 525.65 | 227.10 | UNDEREXPANDED |
| B | 0.050 | 53.820 | 15.10 → 16.50 | 165.9 → 184.1 | 9 513.75 | 211.35 | **OVEREXPANDED** |
| C | 0.150 | 35.974 | 36.65 → 38.73 | 447.9 → 475.2 | 16 698.13 | 233.41 | UNDEREXPANDED |
| D | piecewise | 34.0 active | 26.11 → 38.01 | 310.0 → 465.8 | 12 367.03 | 229.34 | UNDEREXPANDED |

Chamber pressure and thrust both **rise** through every constant-flow burn,
because Milestone 2's fuel flow rises as the port opens while the throat is fixed.

**The low-flow case flips the sign of the pressure term.** Case B runs at only
15–16.5 bar, so its exit pressure (≈ 66 kPa) falls below ambient and the nozzle is
**overexpanded**: its pressure thrust is negative throughout, −11.18 to −9.28 N.
The same fixed `ε = 4` nozzle is therefore underexpanded at the nominal and high
flows and overexpanded at the low flow — a direct consequence of holding the
geometry fixed while the chamber pressure varies by a factor of 2.5 across cases.
Case B still sits above the Summerfield limit (`p_e/p_a ≈ 0.65` against 0.4), so
attached flow remains the expected condition; separation is not modelled.

### Ambient-pressure sensitivity

Deterministic nozzle sensitivity at the fixed reference chamber state — **not**
an altitude-performance envelope; no atmosphere model exists in this project.

| `p_a` [Pa] | Thrust [N] | Pressure thrust [N] | `C_F` | `I_sp` [s] | Regime |
| ---: | ---: | ---: | ---: | ---: | --- |
| 101 325 | 310.048 | +3.866 | 1.51169 | 226.190 | UNDEREXPANDED |
| 54 020 | 324.909 | +18.728 | 1.58414 | 237.032 | UNDEREXPANDED |
| 26 500 | 333.555 | +27.373 | 1.62630 | 243.339 | UNDEREXPANDED |
| 1 197 | 341.504 | +35.322 | 1.66505 | 249.138 | UNDEREXPANDED |
| 0 | 341.880 | +35.699 | 1.66689 | 249.413 | VACUUM |

Sea level to vacuum gains **31.832 N**, exactly `p_a A_e`. The momentum term and
the entire exit state are unchanged — the nozzle is choked, so only the pressure
term moves.

### Combustion-property sensitivity

Deterministic one-factor variation of the prescribed assumptions — **not**
uncertainty propagation. Mass flow and nozzle geometry are held fixed.

| Case | `p_c` [bar] | `M_e` | `V_e` [m/s] | Thrust [N] | ΔF | Implied `η_c*` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 26.114 | 2.6194 | 2190.51 | 310.048 | — | 0.960 |
| `c*` −10 % | 23.503 | 2.6194 | 2190.51 | 306.478 | −1.15 % | 0.864 |
| `c*` +10 % | 28.726 | 2.6194 | 2190.51 | 313.618 | +1.15 % | 1.056 |
| `γ` = 1.15 | 26.114 | 2.5458 | 2219.96 | 317.824 | +2.51 % | 0.945 |
| `γ` = 1.25 | 26.114 | 2.6956 | 2162.64 | 302.873 | −2.31 % | 0.974 |
| `T_c` −10 % | 26.114 | 2.6194 | 2078.10 | 294.336 | −5.07 % | 1.012 |
| `T_c` +10 % | 26.114 | 2.6194 | 2297.42 | 324.992 | +4.82 % | 0.915 |

The structure is worth stating plainly: a 10 % change in `c*` moves chamber
pressure by exactly 10 % but thrust by only 1.15 %, because `V_e` is untouched
and only the small pressure term scales. `T_c` moves `V_e` as `√T_c` and
dominates thrust while leaving `p_c` alone. `γ` shifts the area–Mach solution but
not `p_c`. Because `ṁ` is fixed, `ΔI_sp` is identical to `ΔF` in every row.

Note the `c*` +10 % row implies `η_c* = 1.056` — a physically impossible
efficiency. That is exactly the point of showing the column: varying one
prescribed property alone makes the set internally inconsistent, which is why
these are labelled assumptions rather than data.

### Numerical convergence

Exit-Mach root solve (`ε = 4`, `γ = 1.20`), and transient integration:

| `xtol` | `M_e` | Area–Mach residual |
| ---: | ---: | ---: |
| 1e−04 | 2.619446781721 | 2.7 × 10⁻⁸ |
| 1e−08 | 2.619446776666 | −1.4 × 10⁻¹² |
| 1e−12 | 2.619446776666 | 1.2 × 10⁻¹² |
| 1e−15 | 2.619446776666 | 0 |

| `rtol` | `n_report` | `t_burn` [s] | Peak `F` [N] | `I_total` [N·s] | `I_sp,eq` [s] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1e−05 | 201 | 41.740617467 | 333.528161 | 13 525.649028 | 227.103122 |
| 1e−07 | 401 | 41.740617771 | 333.528161 | 13 525.651206 | 227.103159 |
| 1e−09 | 801 | 41.740617800 | 333.528161 | 13 525.651738 | 227.103168 |
| 1e−11 | 1601 | 41.740617803 | 333.528161 | 13 525.651869 | 227.103170 |

### A result reported rather than corrected

The predicted `I_sp` of 226.2 s **exceeds** the 167–223 s measured by Rezaei et
al. for a comparable motor. That is expected: this is an ideal 1-D isentropic
nozzle with no divergence loss, no boundary layer, no heat loss and no two-phase
effects, and the `c*` efficiency is the *only* loss represented anywhere in the
model. The discrepancy is stated rather than absorbed into a fudge factor.

## Milestone 4 representative results

### Reference configuration (illustrative, generic)

| Parameter | Value | Status |
| --- | ---: | --- |
| Tank internal volume | 10.00 L | round generic value |
| Initial temperature | 293.15 K | round generic value |
| Initial liquid fill | 80 % by volume | round generic value |
| Initial tank pressure | 50.525 bar | **derived** from `p_sat(293.15 K)` |
| Initial N₂O mass | 6.5968 kg (6.2808 liquid) | derived |
| Effective injector area `A_eff` | 3.50 mm² | round, selected as below |
| Discharge coefficient `C_d` | 0.66 | Dyer et al.'s cross-injector average, quoted by Waxman et al. |
| Injector model | Dyer / NHNE | sourced |

The grain, regression law, combustion properties and nozzle are the **frozen**
Milestone 1–3 values, unchanged.

**Effective-area selection** — a comparison convenience, not an injector design.
`C_d` was fixed at the sourced 0.66 and one parameter, the area, was chosen from
a short list of round values so the coupled model starts near Milestone 3's
prescribed 0.100 kg/s. No thrust target was used and no second parameter moved:

| `A_eff` [mm²] | Initial `ṁ_ox` [kg/s] | Initial `p_c` [bar] | Verdict |
| ---: | ---: | ---: | --- |
| 2.5 | 0.075123 | 20.727 | outside the 0.08–0.12 band |
| 3.0 | 0.088644 | 23.671 | inside |
| **3.5** | **0.101316** | **26.396** | **inside — chosen** |
| 4.0 | 0.113049 | 28.894 | inside |
| 4.5 | 0.123791 | 31.164 | outside |

### Reference case (Case A)

| Quantity | At `t = 0` | At termination |
| --- | ---: | ---: |
| Tank pressure | 50.525 bar | 38.094 bar |
| Tank temperature | 293.150 K | 281.117 K |
| Tank liquid mass | 6.2808 kg | 1.7537 kg |
| **Oxidizer mass flow** | **0.101316 kg/s** | **0.083946 kg/s** |
| Fuel mass flow | 0.039968 kg/s | 0.046310 kg/s |
| Mixture ratio O/F | 2.535 | 1.813 |
| **Chamber pressure** | **26.396 bar** | **24.335 bar** |
| Pressure margin | 24.129 bar | 13.759 bar |
| Port diameter | 40.000 mm | 90.000 mm |
| **Thrust** | **313.733 N** | **286.756 N** |

Terminates on **grain burnout** at **42.917 s**. Total impulse **13 018.56 N·s**,
oxidizer consumed 3.982 kg, fuel consumed 1.899 kg, equivalent `I_sp` 225.73 s.
28 % of the liquid N₂O is still in the tank when the grain burns through.

### Case set

| Case | Variation | Termination | Burn [s] | `ṁ_ox` [kg/s] | Thrust [N] | `I_total` [N·s] | `I_sp,eq` [s] |
| --- | --- | --- | ---: | --- | --- | ---: | ---: |
| A | nominal | GRAIN_BURNOUT | 42.92 | 0.1013 → 0.0839 | 313.7 → 286.8 | 13 018.6 | 225.73 |
| B | tank at 283.15 K | GRAIN_BURNOUT | 44.96 | 0.0886 → 0.0748 | 277.9 → 259.8 | 12 198.7 | 223.22 |
| C | tank at 303.15 K | GRAIN_BURNOUT | 41.43 | 0.1115 → 0.0914 | 342.2 → 308.4 | 13 679.1 | 227.48 |
| D | `A_eff` = 2.5 mm² | GRAIN_BURNOUT | 47.41 | 0.0751 → 0.0661 | 239.5 → 233.5 | 11 334.5 | 220.11 |
| E | `A_eff` = 4.5 mm² | GRAIN_BURNOUT | 40.31 | 0.1238 → 0.0963 | 376.2 → 322.9 | 14 223.6 | 228.78 |
| F | SPI injector model | GRAIN_BURNOUT | 40.30 | 0.1259 → 0.0952 | 382.0 → 319.7 | 14 233.7 | 228.79 |
| G | 2.5 L tank | **LIQUID_DEPLETED** | 16.04 | 0.1013 → 0 | 313.7 → 0 | 4 580.4 | 224.40 |

**Thrust falls in every case.** Cases A–F all end on grain burnout — a 10 L tank
carries far more oxidizer than a 1.9 kg grain can consume, so the tank is never
the limiting factor. Case G, with a 2.5 L tank, terminates instead on liquid
depletion after 16.0 s with the port only 61.9 mm across and the grain unburnt.
Which event ends the burn is a model *outcome*, not an assumption.

### Milestone 3 versus Milestone 4

| Quantity | M3 (prescribed) | M4 (coupled) | Change |
| --- | ---: | ---: | ---: |
| Initial `ṁ_ox` [kg/s] | 0.1000 | 0.1013 | +1.32 % |
| Final `ṁ_ox` [kg/s] | 0.1000 | 0.0839 | **−16.05 %** |
| Initial `p_c` [bar] | 26.114 | 26.396 | +1.08 % |
| Final `p_c` [bar] | 27.908 | 24.335 | **−12.80 %** |
| Initial thrust [N] | 310.05 | 313.73 | +1.19 % |
| Final thrust [N] | 333.53 | 286.76 | **−14.02 %** |
| Peak thrust [N] | 333.53 | 313.73 | −5.93 % |
| Final O/F | 2.025 | 1.813 | −10.49 % |
| Burn duration [s] | 41.741 | 42.917 | +2.82 % |
| Total impulse [N·s] | 13 525.65 | 13 018.56 | −3.75 % |
| Equivalent `I_sp` [s] | 227.10 | 225.73 | −0.61 % |

The integrated quantities differ by only a few percent, but the **shape** of the
thrust curve is qualitatively different: prescribed flow rises +7.6 %, coupled
flow falls −8.6 %. Being more coupled does **not** make Milestone 4 validated —
it makes it a different set of assumptions, none of which has been tested against
hardware.

### Independent mass and energy audit

| Check | Residual |
| --- | ---: |
| Oxidizer drawn from tank vs integrated at injector | **2.2 × 10⁻¹⁵ kg** |
| Fuel from port geometry vs integrated `ṁ_f` | **7.0 × 10⁻⁸ kg** |
| Impulse state vs trapezoidal integral of thrust | 9.4 × 10⁻⁴ N·s (7.2 × 10⁻⁸ relative) |
| Adiabatic tank energy balance | 2.4 × 10⁻² J (**1.7 × 10⁻⁸ relative**) |
| Feed-coupling root residual | **1.4 × 10⁻¹³ kg/s** |

Each is reconstructed from the reported histories, not read back from the
solver's own accumulators.

### Convergence

| `rtol` | feed `xtol` | `n_report` | `t_end` [s] | Final `p_tank` [bar] | Final `p_c` [bar] | `I_total` [N·s] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1e−06 | 1e−09 | 101 | 42.917404 | 38.093597 | 24.335095 | 13 018.5817 |
| 1e−08 | 1e−12 | 201 | 42.917324 | 38.093619 | 24.335102 | 13 018.5583 |
| 1e−10 | 1e−12 | 401 | 42.917325 | 38.093618 | 24.335102 | 13 018.5588 |
| 1e−11 | 1e−14 | 801 | 42.917325 | 38.093618 | 24.335102 | 13 018.5588 |

Termination time converges to 8 significant figures and total impulse to
7 — the loosest setting differs by 2 × 10⁻⁶ relative in impulse.

## Milestone 5 representative results

### Static reference comparison — the chemistry alone

At the Milestone 1/3 prescribed flows (`ṁ_ox` = 0.100 kg/s, `ṁ_f` = 0.039776644
kg/s), with no feed coupling at all, so this isolates the chemistry:

| Quantity | M3 frozen | M5 CEA table | Change |
| --- | --- | --- | --- |
| `p_c` [bar] | 26.114 | 24.067 | −7.84 % |
| exit Mach [-] | 2.6194 | 2.6679 | +1.85 % |
| `p_e` [kPa] | 113.63 | 98.44 | −13.37 % |
| `V_e` [m/s] | 2190.5 | 2029.8 | −7.34 % |
| thrust [N] | 310.05 | 282.82 | −8.78 % |
| `C_F` [-] | 1.5117 | 1.4962 | −1.02 % |
| `I_sp` [s] | 226.19 | 206.32 | −8.78 % |

Two channels act, and they are not the same channel:

* the **lower `c*`** lowers chamber pressure, which lowers thrust;
* the **higher `γ`** raises the exit Mach number at the *same* area ratio, which
  lowers the exit pressure and changes the pressure-thrust term.

The chamber-solve residual is 4.7 × 10⁻¹⁰ Pa, and the frozen Milestone 3 nozzle
code is reused unmodified — only the properties handed to it changed.

### Reference case (Case A), coupled

| Quantity | M4 frozen | M5 variable | Change |
| --- | --- | --- | --- |
| burn time [s] | 42.917 | 42.522 | −0.92 % |
| total impulse [N s] | 13018.6 | 11627.0 | **−10.69 %** |
| mean thrust [N] | 303.34 | 273.44 | −9.86 % |
| peak thrust [N] | 313.73 | 290.76 | −7.32 % |
| final thrust [N] | 286.76 | 253.21 | −11.70 % |
| equivalent `I_sp` [s] | 225.73 | 199.45 | −11.64 % |
| oxidizer consumed [kg] | 3.9821 | 4.0454 | +1.59 % |
| fuel consumed [kg] | 1.8991 | 1.8991 | 0.00 % |
| **thrust change over the burn** | **−8.60 %** | **−12.92 %** | 1.50× |

Both runs terminate on grain burnout, so the fuel consumed is identical — it is
set by geometry, not by chemistry. The oxidizer consumed *rises*, because the
lower chamber pressure lets the injector pass more.

Over the burn the chemistry moves as the port opens:

| | start | end |
| --- | --- | --- |
| O/F [-] | 2.557 | 1.856 |
| `p_c` [bar] | 24.671 | 21.846 |
| delivered `c*` [m/s] | 1356.2 | 1279.8 |
| `γ` [-] | 1.23163 | 1.23490 |
| `T_c` [K] | 2186.0 | 1925.8 |

`γ` barely moves over this range, so the extra decay is almost entirely a `c*`
effect acting through chamber pressure.

### Case set

Every case is run **twice on identical conditions** — once with the frozen
Milestone 4 model, once with the Milestone 5 table — so each row differs only by
the combustion properties.

| Case | Configuration | M4 impulse [N s] | M5 impulse [N s] | Change | M4 decay | M5 decay |
| --- | --- | --- | --- | --- | --- | --- |
| A | nominal | 13018.6 | 11627.0 | −10.69 % | −8.60 % | −12.92 % |
| B | tank at 283.15 K | 12198.7 | 10772.8 | −11.69 % | −6.52 % | −10.57 % |
| C | tank at 303.15 K | 13679.1 | 12323.9 | −9.91 % | −9.86 % | −14.29 % |
| D | injector 2.5 mm² | 11334.5 | 9649.8 | −14.86 % | −2.49 % | −7.56 % |
| E | injector 4.5 mm² | 14223.6 | 13078.4 | −8.05 % | −14.15 % | −17.81 % |
| G | 2.5 L tank | 4580.4 | 4061.2 | −11.34 % | — | — |

Case F is the frozen Milestone 4 model itself — the M4 column of every row. It is
preserved unchanged as the controlled comparison, not re-derived or re-tuned.

Every case loses impulse, and every grain-limited case decays more steeply. The
effect is largest for the small injector (Case D, −14.9 %), which runs at the
lowest chamber pressure and hence the lowest O/F, furthest down the fuel-rich side
of the table. Cases A–E terminate on grain burnout; Case G, with a 2.5 L tank,
terminates on liquid depletion, so both terminal branches are exercised.

**Every case stayed strictly inside the tabulated rectangle.** No run approached
an O/F or pressure edge, so no result depends on the table's boundary behaviour.

### Does variable chemistry overturn the Milestone 4 result?

No — it deepens it. Milestone 4's finding was that thrust *falls* over the burn
even though the port is opening up, because the tank blows down faster than the
growing burning area can compensate. Milestone 5 keeps the sign and multiplies the
magnitude by 1.50 in Case A.

The mechanism is genuinely new. In Milestone 4 the only decay channel was the
tank, because `c*` was a constant. Here the falling mixture ratio drags `c*` down
with it, and the two effects add. Separately, the whole curve sits **lower**,
because the equilibrium `c*` at this motor's mixture ratio is below the
illustrative constant Milestone 3 prescribed. Level and slope are separate
findings and are reported separately.

### Sensitivity

**`c*` efficiency** — swept 0.92 → 1.00:

| `η_c*` | total impulse [N s] | mean thrust [N] | burn time [s] |
| --- | --- | --- | --- |
| 0.92 | 11632.4 | 274.34 | 42.402 |
| 0.94 | 11630.1 | 273.90 | 42.461 |
| 0.96 | 11627.0 | 273.44 | 42.522 |
| 0.98 | 11623.2 | 272.95 | 42.584 |
| 1.00 | 11618.6 | 272.43 | 42.648 |

Total impulse barely moves — about 0.1 % across an 8.7 % change in `η_c*` — and it
moves *downwards* as `η` rises. That is a property of the coupling, not of the
chemistry: a higher `c*` raises chamber pressure (23.78 → 25.54 bar), which cuts
the pressure drop the injector sees so the oxidizer flow falls (0.10334 → 0.10203
kg/s), while the grain-limited burn fixes the fuel mass. The pressure-fed system
self-compensates and what is left is a slightly longer, slightly softer burn.

So `η_c*` is far less consequential for total impulse than its size suggests — but
it still moves chamber pressure by ~7 %, which is what any pressure-driven
conclusion would rest on.

**Interpolation grid resolution** — coarser tables built by taking every n-th node
of the committed table, so the coarse nodes are an exact subset of the same CEA
solutions and the difference is pure interpolation error:

| O/F step | `p_c` step [bar] | nodes | total impulse [N s] | vs full |
| --- | --- | --- | --- | --- |
| 0.025 | 2.5 | 1695 | 11627.03 | — |
| 0.050 | 2.5 | 855 | 11626.97 | −0.0005 % |
| 0.100 | 5.0 | 232 | 11626.74 | −0.0025 % |
| 0.200 | 10.0 | 75 | 11625.77 | −0.0109 % |
| 0.400 | 17.5 | 24 | 11622.61 | −0.0380 % |

The error falls monotonically and is already below 0.04 % on the coarsest grid
tried. The committed resolution is far finer than the answer needs, deliberately:
the table is generated once and costs nothing at run time.

**2-D table vs O/F-only interpolation** — freezing every property at 25 bar and
keeping only the mixture-ratio dependence changes total impulse by +0.003 %. Across
the whole tabulated span (10 → 45 bar) the ideal `c*` at O/F 2.5 moves only
+0.029 %, so this is expected rather than surprising. The second dimension was
tabulated anyway: it costs nothing at run time and it lets the chamber closure be
solved implicitly in `p_c` rather than assumed independent of it.

**Equilibrium vs frozen nozzle chemistry** — CEA's frozen option converged at 19
of 29 sampled O/F points. It **failed to converge over O/F 2.10–3.00**, where
condensed carbon is present, and that band straddles this motor's operating range.
That non-convergence is reported rather than worked around: no frozen number is
invented where the solver could not produce one. Where it did converge, frozen
vacuum `I_sp` is −1.40 % to −0.29 % relative to equilibrium (mean −0.80 %).
Bracketing the operating band from both sides, the expansion-chemistry assumption
is worth roughly 1 % in `I_sp` — an order of magnitude smaller than the 8–12 % the
property change is worth, but not negligible, and it is **not** modelled here.

### Independent verification

Every residual below is zero in exact arithmetic and is re-formed from the
reported histories rather than read back from the solver:

| Case | ox mass [kg] | fuel mass [kg] | `c*` identity [Pa] | O/F identity | feed [kg/s] | impulse |
| --- | --- | --- | --- | --- | --- | --- |
| A | 4.7e−15 | 5.1e−07 | 9.3e−10 | 0.0 | 3.1e−13 | 4.3e−08 |
| B | 1.0e−14 | 4.0e−07 | 9.3e−10 | 0.0 | 3.1e−13 | 1.4e−08 |
| C | 4.4e−15 | 5.7e−07 | 9.3e−10 | 0.0 | 2.8e−13 | 7.6e−08 |
| D | 4.4e−15 | 3.2e−07 | 7.0e−10 | 0.0 | 1.6e−13 | 2.5e−08 |
| E | 4.9e−15 | 1.5e−08 | 9.3e−10 | 0.0 | 2.7e−13 | 3.1e−08 |
| G | 1.0e−15 | 4.0e−08 | 9.3e−10 | 0.0 | 3.0e−13 | 1.1e−03 |

The mixture-ratio identity holds **exactly** — the reported O/F is the same float
as `ṁ_ox / ṁ_f`. Case G's impulse residual is large because that run ends on
liquid depletion, where thrust steps discontinuously to zero and a trapezoidal
reconstruction cannot represent the step; the integrator's own impulse state is
the accurate quantity, as documented for Milestone 4.

## Installation

```bash
pip install -e ".[dev]"
```

Requires Python ≥ 3.11. Runtime dependencies are `numpy`, `scipy` (ODE solvers
and root finders), `CoolProp` (the N₂O equation of state used by Milestone 4)
and `matplotlib`. The `dev` extra installs `pytest` and `ruff`, which is
everything the verification workflow below needs.

The Milestone 5 thermochemistry table is **committed**, so using it needs nothing
extra. Only *regenerating* it needs NASA CEA:

```bash
pip install -e ".[cea]"
```

`rocketcea` wraps the CEA FORTRAN source and needs a Fortran compiler, which is
exactly why it is an optional extra rather than a runtime dependency.

## Running the study

```bash
python scripts/regression_study.py
```

Prints the configuration, the oxidizer-flow sweep, the port-size sensitivity, the
illustrative exponent sensitivity, the observations, and an engineering sanity
audit. Exits non-zero if any sanity check fails.

Longhand cross-check of the whole Milestone 1 model chain against hand
arithmetic:

```bash
python scripts/manual_check.py
```

Run the Milestone 2 transient study (four prescribed-flow cases, the convergence
study and the engineering sanity audit):

```bash
python scripts/transient_regression_study.py
```

Run the Milestone 3 thrust study (reference point with three independent thrust
routes, four transient cases, ambient-pressure and combustion-property
sensitivity, convergence and the sanity audit):

```bash
python scripts/thrust_prediction_study.py
```

Run the Milestone 4 coupled blowdown study (property and injector audit,
effective-area calibration, seven cases, the M3-vs-M4 comparison, an independent
mass/energy audit, convergence and the sanity audit):

```bash
python scripts/blowdown_thrust_study.py
```

Run the Milestone 5 O/F-dependent thermochemistry study (table provenance, the
M3-vs-M5 property and static comparison, six cases each run under both models, the
four sensitivities and the closure audit):

```bash
python scripts/variable_thermochemistry_study.py
```

Regenerate the committed CEA table and prove the committed copy is what the
current generator produces (needs the optional `cea` extra):

```bash
python scripts/generate_thermochemistry_table.py --check
```

Regenerate the figures (deterministic — byte-identical on repeated runs within a
given environment; see note below):

```bash
python scripts/generate_m1_figures.py
```

```bash
python scripts/generate_m2_figures.py
```

```bash
python scripts/generate_m3_figures.py
```

```bash
python scripts/generate_m4_figures.py
```

```bash
python scripts/generate_m5_figures.py
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

854 tests (167 Milestone 1, 119 Milestone 2, 215 Milestone 3, 211 Milestone 4,
142 Milestone 5).
Expected values are written independently of the production code —
independently arranged algebra, hand-computed literals, an independent
unit-conversion route, all 18 published `(G_ox, ṙ)` measurements from the source
paper's own data table, a closed-form transient solution re-derived inside the
test file, a longhand reimplementation of the chamber/nozzle chain with a
plain-bisection area–Mach inversion, and for Milestone 4 an independently
retrieved NIST WebBook property table, CoolProp's own density–internal-energy
flash as a second route to the tank temperature, hand-written SPI/HEM/Dyer
formulas, and a plain-bisection re-solve of the coupled feed closure.

For Milestone 5 the committed CSV is re-parsed with the standard library rather
than the package loader; `R = R_u / M` and the ideal-`c*` relation are re-derived
from the raw columns; the implicit chamber closure is re-solved by an independent
scan-plus-secant method; all six coupled feed relations are re-formed from the
frozen Milestone 1–4 components; the closed-form O/F inversion is round-tripped
through the forward regression chain; blowdown closures are rebuilt from the
reported histories; sampled states are re-solved from scratch with the static
solver; and the frozen Milestone 4 model is re-run in the same session and
asserted to give its own unchanged answer.

## Figures

Milestone 1:

| Figure | File |
| --- | --- |
| **A** — Generic grain geometry schematic | [`figures/fig_a_grain_geometry.png`](figures/fig_a_grain_geometry.png) |
| **B** — Regression rate vs oxidizer mass flux | [`figures/fig_b_regression_vs_flux.png`](figures/fig_b_regression_vs_flux.png) |
| **C** — Port-size sensitivity at fixed `ṁ_ox` | [`figures/fig_c_port_sensitivity.png`](figures/fig_c_port_sensitivity.png) |

Milestone 2:

| Figure | File |
| --- | --- |
| **M2-A** — Port diameter and remaining web vs time | [`figures/fig_m2_a_port_growth.png`](figures/fig_m2_a_port_growth.png) |
| **M2-B** — `G_ox` and `ṙ` histories, with the source flux band | [`figures/fig_m2_b_flux_and_rate.png`](figures/fig_m2_b_flux_and_rate.png) |
| **M2-C** — `ṁ_f` and O/F histories | [`figures/fig_m2_c_fuel_flow_and_of.png`](figures/fig_m2_c_fuel_flow_and_of.png) |
| **M2-D** — Fuel mass bookkeeping and conservation residual | [`figures/fig_m2_d_mass_closure.png`](figures/fig_m2_d_mass_closure.png) |
| **M2-E** — Response to a piecewise prescribed schedule | [`figures/fig_m2_e_piecewise.png`](figures/fig_m2_e_piecewise.png) |

Milestone 3:

| Figure | File |
| --- | --- |
| **M3-A** — Chamber pressure vs time | [`figures/fig_m3_a_chamber_pressure.png`](figures/fig_m3_a_chamber_pressure.png) |
| **M3-B** — **Predicted thrust vs time** | [`figures/fig_m3_b_thrust.png`](figures/fig_m3_b_thrust.png) |
| **M3-C** — Nozzle pressure ratios and expansion regime | [`figures/fig_m3_c_nozzle_state.png`](figures/fig_m3_c_nozzle_state.png) |
| **M3-D** — Ambient-pressure sensitivity | [`figures/fig_m3_d_ambient_sensitivity.png`](figures/fig_m3_d_ambient_sensitivity.png) |
| **M3-E** — Prescribed-property sensitivity | [`figures/fig_m3_e_property_sensitivity.png`](figures/fig_m3_e_property_sensitivity.png) |
| **M3-F** — Piecewise-schedule thrust response | [`figures/fig_m3_f_piecewise_thrust.png`](figures/fig_m3_f_piecewise_thrust.png) |

Milestone 4:

| Figure | File |
| --- | --- |
| **M4-A** — Tank blowdown: pressure, temperature, mass | [`figures/fig_m4_a_tank_blowdown.png`](figures/fig_m4_a_tank_blowdown.png) |
| **M4-B** — Feed coupling: flow and chamber pressure, M3 vs M4 | [`figures/fig_m4_b_feed_coupling.png`](figures/fig_m4_b_feed_coupling.png) |
| **M4-C** — **Thrust: prescribed vs coupled** | [`figures/fig_m4_c_thrust_comparison.png`](figures/fig_m4_c_thrust_comparison.png) |
| **M4-D** — Injector characteristics: SPI, HEM, Dyer | [`figures/fig_m4_d_injector_models.png`](figures/fig_m4_d_injector_models.png) |
| **M4-E** — Parametric cases and conservation closure | [`figures/fig_m4_e_case_comparison.png`](figures/fig_m4_e_case_comparison.png) |

![Generic hybrid grain geometry](figures/fig_a_grain_geometry.png)

![Regression rate versus oxidizer mass flux](figures/fig_b_regression_vs_flux.png)

![Port-size sensitivity](figures/fig_c_port_sensitivity.png)

![Port growth under constant prescribed oxidizer flow](figures/fig_m2_a_port_growth.png)

![Oxidizer mass flux and regression rate histories](figures/fig_m2_b_flux_and_rate.png)

![Fuel mass flow and mixture-ratio histories](figures/fig_m2_c_fuel_flow_and_of.png)

![Fuel mass bookkeeping and conservation](figures/fig_m2_d_mass_closure.png)

![Response to a piecewise prescribed oxidizer-flow schedule](figures/fig_m2_e_piecewise.png)

![Chamber pressure history](figures/fig_m3_a_chamber_pressure.png)

![Predicted thrust history](figures/fig_m3_b_thrust.png)

![Nozzle state](figures/fig_m3_c_nozzle_state.png)

![Ambient-pressure sensitivity](figures/fig_m3_d_ambient_sensitivity.png)

![Prescribed-property sensitivity](figures/fig_m3_e_property_sensitivity.png)

![Piecewise-schedule thrust response](figures/fig_m3_f_piecewise_thrust.png)

![Tank blowdown](figures/fig_m4_a_tank_blowdown.png)

![Feed coupling](figures/fig_m4_b_feed_coupling.png)

![Thrust: prescribed versus coupled](figures/fig_m4_c_thrust_comparison.png)

![Injector characteristics](figures/fig_m4_d_injector_models.png)

![Parametric cases and closure](figures/fig_m4_e_case_comparison.png)

Milestone 5:

| Figure | File |
| --- | --- |
| **M5-A** — CEA equilibrium properties vs O/F, with the M3 constants | [`figures/fig_m5_a_thermochemistry_table.png`](figures/fig_m5_a_thermochemistry_table.png) |
| **M5-B** — Coupled thrust and impulse, M4 vs M5 | [`figures/fig_m5_b_thrust_comparison.png`](figures/fig_m5_b_thrust_comparison.png) |
| **M5-C** — What the chemistry does during a burn | [`figures/fig_m5_c_property_histories.png`](figures/fig_m5_c_property_histories.png) |
| **M5-D** — The case set and its closure residuals | [`figures/fig_m5_d_case_comparison.png`](figures/fig_m5_d_case_comparison.png) |
| **M5-E** — Efficiency, grid resolution and expansion chemistry | [`figures/fig_m5_e_sensitivities.png`](figures/fig_m5_e_sensitivities.png) |

## Repository layout

```
src/hybrid_rocket_motor/
    geometry.py          grain geometry, areas, volumes, masses, validation      [M1]
    operating_point.py   prescribed-ṁ_ox snapshot: G_ox, ṁ_f, O/F bookkeeping    [M1]
    regression.py        ṙ = a G_ox ⁿ with explicit source-unit handling         [M1]
    transient.py         flow histories, ODE solver, burnout event, closed form  [M2]
    nozzle.py            geometry, area-Mach solve, isentropic exit state        [M3]
    chamber.py           prescribed combustion properties, quasi-steady p_c      [M3]
    performance.py       chamber + nozzle coupling, thrust, total impulse        [M3]
    nitrous_properties.py  saturated N2O properties (CoolProp / Lemmon-Span)     [M4]
    tank.py              rigid adiabatic saturated-equilibrium N2O tank          [M4]
    injector.py          SPI / HEM / Dyer-NHNE injector flow models              [M4]
    feed_system.py       coupled injector-chamber-regression root solve          [M4]
    blowdown.py          coupled time integration with terminal events           [M4]
    thermochemistry.py   frozen CEA table, bilinear interpolation, refusal       [M5]
    variable_chamber.py  implicit p_c = ṁ c*(O/F, p_c) / A_t closure             [M5]
    variable_feed_system.py  nested feed / chemistry root solve                  [M5]
    variable_nozzle.py   adapter driving the frozen M3 nozzle from the table     [M5]
    variable_blowdown.py coupled integration with table-validity events          [M5]
data/
    n2o_htpb_equilibrium.csv           1695 CEA points, committed, no timestamp  [M5]
    n2o_htpb_equilibrium_metadata.json solver and version provenance             [M5]
    n2o_htpb_nozzle_chemistry.csv      equilibrium vs frozen expansion           [M5]
tests/                   independent verification (854 tests)
scripts/
    manual_check.py                longhand arithmetic cross-check + scope guard [M1]
    regression_study.py            the Milestone 1 study and sanity audit        [M1]
    generate_m1_figures.py         deterministic figure generation               [M1]
    transient_regression_study.py  the Milestone 2 study, convergence and audit  [M2]
    generate_m2_figures.py         deterministic figure generation               [M2]
    thrust_prediction_study.py     the Milestone 3 study, sensitivity and audit  [M3]
    generate_m3_figures.py         deterministic figure generation               [M3]
    blowdown_thrust_study.py       the Milestone 4 study, comparison and audit   [M4]
    generate_m4_figures.py         deterministic figure generation               [M4]
    generate_thermochemistry_table.py  one-time CEA table generation (needs CEA) [M5]
    variable_thermochemistry_study.py  the Milestone 5 study, sensitivity, audit [M5]
    generate_m5_figures.py         deterministic figure generation               [M5]
figures/                 twenty-four PNG figures
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
* **Nothing here has been validated against hardware.**

Additional Milestone 2 limitations:

* **Most of a constant-flow burn is extrapolated.** Only 33.7 % / 3.8 % / 61.6 %
  of Cases A / B / C lie inside the correlation's measured flux range, so the
  burnout times are extrapolated results, flagged as such rather than hidden.
* **Spatially uniform regression.** The port is a single scalar radius, so axial
  variation of `ṙ` is ignored entirely. Real single-port grains regress
  non-uniformly and leave a residual sliver at burnout; this model instead
  consumes the grain exactly and uniformly, so its burnout is an idealised upper
  bound on fuel utilisation.
* **Quasi-steady regression.** The steady-state correlation is applied
  instantaneously at each time, with no thermal lag in the solid, no transient
  boundary-layer response and no combustion-delay dynamics.
* **The grain length dependence is still unmodelled**, and the Milestone 1
  caveat that predictions may be ~15 % low for this 400 mm grain applies to every
  transient result here as well.
* **No coupling back from the chamber.** Oxidizer flow is prescribed, so nothing
  represents how a real feed system and chamber would interact with the growing
  port.

Additional Milestone 3 limitations:

* **`ṁ_ox` is prescribed externally and chamber pressure does NOT feed back into
  it.** This is the single largest departure from a real motor: a real oxidizer
  feed would deliver less flow as chamber pressure rises, damping exactly the
  thrust rise this model predicts.
* **`c*`, `γ` and `T_c` are prescribed constants** and do not vary with O/F,
  chamber pressure or time. Real combustion properties vary with all three. There
  is no equilibrium chemistry and no CEA interface.
* **Combustion efficiency is represented only by the prescribed `c*`.** No other
  loss mechanism exists anywhere in the model.
* **The nozzle is 1-D, steady, adiabatic and isentropic** with a calorically
  perfect gas. No divergence loss, no boundary layer, no heat transfer, no
  two-phase flow, no shocks, no separation model, no side loads and no
  altitude-compensating behaviour. The Summerfield flag only *marks* results that
  would be doubtful; it does not correct them.
* **No ignition, chamber-filling or tail-off transient.** The chamber has no gas
  storage term, so every step in prescribed flow produces an instantaneous step in
  pressure and thrust. A real motor cannot do that.
* **No combustion instability, no structural sizing, no thermal sizing, no nozzle
  contour, no materials and no fabrication information.**
* **No experimental thrust validation.** The reference geometry is illustrative
  and not taken from, nor validated against, any real hardware. Nothing here
  establishes a safe operating pressure.

Additional Milestone 4 limitations:

* **The equilibrium tank model is the simplest defensible closure, and it is
  known to be imperfect.** Zimmerman (Stanford PhD, 2015) shows experimentally
  that self-pressurising tanks exhibit an initial transient with rapid pressure
  fluctuations and bubble nucleation that equilibrium models do not reproduce.
  The equilibrium model is also reported to over-predict tank pressure and is
  outperformed by the Zilliac–Karabeyoglu and Casalino–Pastrone models. It is
  used here because it is a genuine thermodynamic closure rather than a fitted
  decay curve, and needs saturation properties only.
* **The tank is adiabatic, uniform and in instantaneous phase equilibrium.** No
  wall thermal mass, no heat leak, no thermal stratification, no ullage dynamics,
  no metastable superheat.
* **The vapour-only discharge tail is not modelled.** A run that empties its
  liquid terminates there and says so, rather than inventing a gas-blowdown model.
* **The injector is an effective area only.** `C_d A` is a lumped generic
  parameter. No hole count, hole diameter, plate geometry or manufacturing
  information is produced, and `C_d = 0.66` is a literature average for *other*
  injectors, not a measurement of anything modelled here.
* **The Dyer/NHNE blend is an empirical correlation**, reported by its authors as
  accurate to about ±15 % against a limited hot-fire set. This project uses the
  plain HEM branch at the actual backpressure, not the critical (maximised) HEM
  value; Waxman et al. note HEM "is not necessarily suitable" beyond the critical
  pressure drop.
* **No feed-line pressure loss, no valve dynamics, no injector transient.** The
  tank connects to the chamber through a single lumped area with no plumbing.
* **Combustion properties were still prescribed constants at Milestone 4** and did
  not respond to the large O/F excursion the coupling produces (2.53 → 1.81 in
  Case A). **Milestone 5 removes this limitation** and quantifies what it was
  worth: −10.7 % total impulse and a 1.5× deeper thrust decay.
* **Being more coupled is not being more validated.** Milestone 4 replaces one
  assumption with three new ones. None of them has been checked against hardware.

Milestone 5 adds its own:

* **Equilibrium is not kinetics.** The table assumes reaction is instantaneous and
  complete at every instant. Real hybrids have finite-rate chemistry, incomplete
  mixing and a diffusion flame that does not reach equilibrium everywhere. Nothing
  here models reaction rates.
* **One `η_c*` stands in for every combustion loss** — mixing, incomplete
  combustion, heat loss, residence time — and it is not measured. It is the
  largest unvalidated number left in the chain. Notably it barely moves total
  impulse (§ Sensitivity), but it moves chamber pressure by ~7 %.
* **The nozzle expands with a single chamber `γ`.** That sits between CEA's
  shifting and frozen limits, and CEA cannot even run the frozen case across
  O/F 2.10–3.00 because of condensed carbon — a band that straddles this motor's
  operating range. The assumption is worth roughly 1 % in `I_sp` where it can be
  measured at all.
* **Condensed carbon is carried as mass, not as a phase.** Soot contributes to the
  effective molar mass but there is no two-phase flow, no particle lag and no
  particulate drag in the nozzle.
* **The table is a rectangle, and the model refuses to leave it.** That is a
  deliberate safety property, not a physical statement: a configuration whose O/F
  or chamber pressure leaves the tabulated span simply stops, and would need the
  table regenerated over a wider range.

## Next milestone

**Milestone 6 (not started).** With the chemistry now responding to the operating
point, the largest remaining gaps are no longer about *what* the products are but
about *how fast* and *how completely* they form, and about what the model still
refuses to look at:

* **combustion efficiency as physics rather than a constant** — `η_c*` is one
  unmeasured number standing in for mixing, residence time and heat loss, and it
  is the single largest unvalidated quantity left;
* **the vapour-only discharge tail**, still unmodelled: runs that deplete their
  liquid stop rather than continuing on vapour blowdown;
* **ignition and chamber-filling transients**, still absent at both ends of every
  burn.

Nothing in the current repository models reaction rates, ignition, instability,
nozzle contours, structures or thermal response.

## License

MIT
