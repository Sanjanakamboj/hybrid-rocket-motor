# Verification

Every number on this page comes from a run made during the Milestone 6 audit
session, not copied from an earlier checkpoint report. Reproduce it with:

```bash
python scripts/independent_audit.py
```

```bash
python scripts/final_robustness_study.py
```

**Independence strategy.** The "expected" side of every check below is built from
an algebraic identity written out longhand, a closed-form solution derived
independently of the solver, an alternate numerical route sharing no code with the
production path, published experimental data, or raw CSV columns parsed with the
standard library. Nothing is verified by calling the same function twice.

**Two tolerance standards, kept separate.** `identity` rows are algebraic
relations that must hold to floating-point precision. `discretisation` rows are
produced by ODE integration or quadrature, where the honest expectation is the
integrator's own tolerance — **these are never reported as if they should reach
machine epsilon**. `empirical` rows compare against measured data, where the
tolerance is the scatter of the experiment.

**Headline: 59 of 59 independent checks pass. The worst algebraic-identity
residual across the whole chain is 9.313 × 10⁻¹⁰ Pa**, on a chamber pressure of
about 2.4 × 10⁶ Pa — a relative residual near 4 × 10⁻¹⁶.

---

## M1 regression and geometry

| Check | Abs. residual | Tolerance | Kind |
| --- | --- | --- | --- |
| `A_port = π D²/4` | 0 | 1e-15 | identity |
| `A_burn = π D L` | 0 | 1e-15 | identity |
| `G_ox = ṁ_ox / A_port` | 0 | 1e-15 | identity |
| `ṙ` rebuilt from the **source-unit** correlation | 0 | 1e-12 | identity |
| `ṁ_f = ρ_f A_burn ṙ` | 0 | 1e-12 | identity |
| `O/F = ṁ_ox / ṁ_f` | 0 | 1e-12 | identity |
| Worst of 18 published firings | 1.085e-01 | 2.0e-01 | empirical |
| Mean deviation over those firings | 5.022e-02 | 1.0e-01 | empirical |

The regression law is rebuilt here from the paper's own units — `a` = 0.3977
mm/s, `n` = 0.3667 with flux in g/(cm²·s) — and converted inside the audit, so the
whole unit chain is exercised rather than the stored SI coefficient being trusted.

All 18 `(G_ox, ṙ)` pairs from Table 4 of Rezaei et al. (2018) are reproduced within
10.85 % worst-case and 5.02 % mean. That is agreement with the **scatter of the
source experiment**, which is the only meaningful standard for a fitted
correlation — not a claim of predictive accuracy for any other motor.

Reference point, `ṁ_ox` = 0.100 kg/s:

| Quantity | Value |
| --- | --- |
| `G_ox` | 79.577472 kg/(m²·s) |
| `ṙ` | 0.850894 mm/s |
| `ṁ_f` | 0.039776644 kg/s |
| `O/F` | 2.514038 |

## M2 transient regression

| Check | Abs. residual | Tolerance | Kind |
| --- | --- | --- | --- |
| Analytic burnout time | 1.647e-10 s | 1e-9 rel | identity |
| Analytic radius `r(t)`, worst of 5 samples | 1.397e-11 | 1e-7 | discretisation |
| Oxidizer mass integral `= ṁ t` | 8.882e-16 kg | 1e-9 rel | identity |
| Fuel mass: integral vs swept geometry | 4.288e-12 kg | 1e-8 | discretisation |
| Loaded fuel mass `= ρ π (R_o² − R_i²) L` | 4.441e-16 kg | 1e-14 rel | identity |

The closed form is re-derived inside the audit rather than imported. With
`ṙ = a (ṁ/(π r²))ⁿ`,

```
r(t) = [ r₀^(2n+1) + (2n+1) a (ṁ/π)ⁿ t ]^(1/(2n+1))
```

giving a burnout time of **41.740618 s** for the constant 0.100 kg/s case, matched
by the integrator to 1.6 × 10⁻¹⁰ s. Port diameter runs 40 → 90 mm.

## M3 chamber and nozzle

| Check | Rel. residual | Tolerance | Kind |
| --- | --- | --- | --- |
| `p_c A_t = ṁ c*` | 0 | 1e-14 | identity |
| Exit Mach by independent bisection | 1.695e-16 | 1e-11 | identity |
| `p_e / p_c` isentropic | 0 | 1e-11 | identity |
| `T_e / T_c` isentropic | 0 | 1e-11 | identity |
| `V_e = M_e √(γ R T_e)` | 0 | 1e-11 | identity |
| `F = ṁ V_e + (p_e − p_a) A_e` | 0 | 1e-11 | identity |
| `C_F = F / (p_c A_t)` | 0 | 1e-13 | identity |
| Thrust via the `C_F p_c A_t` route | 0 | 1e-13 | identity |
| `I_sp = F / (ṁ g₀)` | 0 | 1e-13 | identity |
| `I_sp = c* C_F / g₀` | 0 | 1e-13 | identity |
| Total impulse vs independent trapezoid | 0 | 1e-6 | discretisation |

The area–Mach relation is written out longhand and solved by **plain bisection**
using no production routine; it agrees with the production solve to 1.7 × 10⁻¹⁶
relative. Thrust is confirmed through two independent routes — the momentum plus
pressure decomposition, and the thrust-coefficient route — which agree exactly.

Reference point (`ṁ_ox` = 0.100, `ṁ_f` = 0.039776644 kg/s, sea level):

| Quantity | Value |
| --- | --- |
| `p_c` | 26.114 bar |
| Thrust | 310.048 N |
| `I_sp` | 226.190 s |
| Total impulse (constant-flow burn) | 13 525.65 N·s |
| Equivalent `I_sp` | 227.103 s |

## M4 tank, feed and blowdown

| Check | Residual | Tolerance | Kind |
| --- | --- | --- | --- |
| Tank phase mass closure | 0 | 1e-14 | identity |
| Tank volume closure | 0 | 1e-14 | identity |
| Vapour quality `= m_v / m_total` | 0 | 1e-14 | identity |
| `p_sat(T)` via an independent CoolProp call | 0 | 1e-12 | identity |
| Liquid density via an independent CoolProp call | 0 | 1e-12 | identity |
| `ṁ_SPI = C_d A √(2 ρ Δp)` longhand | 0 | 1e-13 | identity |
| Dyer `κ = 1` for a saturated tank | 0 | 1e-12 | identity |
| Dyer `= (SPI + HEM)/2` when `κ = 1` | 0 | 1e-13 | identity |
| Feed residual re-formed at 3 samples | 1.006e-13 kg/s | 1e-11 | identity |
| Oxidizer mass: tank drawn vs integrated | 1.776e-15 kg | 1e-12 | discretisation |
| Fuel mass: integral vs swept geometry | 6.964e-08 kg | 1e-6 | discretisation |
| Tank energy balance, relative | 6.777e-08 | 1e-6 | discretisation |
| Impulse vs trapezoid, relative | 2.911e-07 | 1e-6 | discretisation |

The saturation pressure and liquid density are re-fetched through a **different
CoolProp entry point** than the production path uses. The honesty caveat recorded
since Milestone 4 still stands: NIST cites the same Lemmon & Span (2006) equation
of state that CoolProp implements, so this agreement verifies units and call
structure, **not the equation of state itself**.

Nominal case:

| Quantity | Value |
| --- | --- |
| Initial tank pressure | 50.525 bar |
| Final tank pressure | 38.094 bar |
| Final tank temperature | 281.117 K |
| Thrust | 313.733 → 286.756 N (**−8.599 %**) |
| Burn duration | 42.917 s |
| `O/F` | 2.535 → 1.813 |
| Terminal event | `GRAIN_BURNOUT` |

A 2.5 L tank instead terminates on `LIQUID_DEPLETED` at 16.035 s. Substituting the
bare SPI injector model raises initial thrust from 313.733 N to **381.967 N**,
which is **+21.75 %**.

> **Correction to the record.** The Milestone 4 commit message stated this as
> "+18.6 %" alongside the thrust pair "313.7 → 382.0 N". Those two statements are
> mutually inconsistent: 382.0 / 313.7 − 1 = +21.75 %. The newton values are
> correct and match live computation and the README's Case F row; the percentage
> in that commit message is wrong. Commit messages are immutable, so the value is
> corrected here rather than by rewriting history. No code and no published table
> was affected.

## M5 thermochemistry

| Check | Residual | Tolerance | Kind |
| --- | --- | --- | --- |
| Grid is complete and rectangular (1695 rows) | 0 | exact | identity |
| `O/F` axis strictly increasing | 0 | exact | identity |
| Pressure axis strictly increasing | 0 | exact | identity |
| Grid-point retrieval, worst of 9 nodes | 0 | 1e-14 | identity |
| Bilinear cell centre `=` mean of 4 raw corners | 0 | 1e-14 | identity |
| `R = R_u / M` off-node | 0 | 1e-15 | identity |
| `c*_delivered = η c*_ideal` | 0 | 1e-15 | identity |
| `η` does not touch `γ` | 0 | 1e-15 | identity |
| Implicit `p_c` by independent bisection | 0 | 1e-12 | identity |
| Chamber residual `ṁ c*/A_t − p_c` | 4.657e-10 Pa | 1e-6 | identity |
| Coupled: `G_ox` re-formed | 1.739e-16 | 1e-14 | identity |
| Coupled: `ṁ_f` re-formed | 0 | 1e-13 | identity |
| Coupled: `O/F` identity | 0 | 1e-15 | identity |
| Coupled: `p_c = ṁ c*/A_t` | 1.887e-16 | 1e-11 | identity |
| Coupled: injector agrees with `p_c` | 1.216e-15 | 1e-12 | identity |
| Table edge classification (5 cases) | 0 | exact | identity |
| Extrapolation refused (4 cases) | 0 | exact | identity |
| Blowdown: oxidizer mass closure | 4.441e-15 kg | 1e-12 | discretisation |
| Blowdown: `c*` identity over history | 9.313e-10 Pa | 1e-6 | identity |
| Blowdown: `O/F` identity over history | **0 exactly** | 1e-14 | identity |
| Blowdown: feed residual over history | 2.909e-13 kg/s | 1e-11 | identity |
| Blowdown: stayed inside the table | 0 | exact | identity |

**Table provenance.** 1695 NASA CEA equilibrium points (RP-1311, via `rocketcea`
1.2.3) over `O/F` 1.20–4.00 and 10–45 bar, committed to the repository. The CSV
carries no timestamp, so `python scripts/generate_thermochemistry_table.py --check`
proves byte-identical regeneration; both committed data files pass. The runtime
never imports `rocketcea`.

**Interpolation.** Bilinear, verified two ways: exact retrieval at grid points read
straight from the raw CSV, and a cell-centre value confirmed to equal the mean of
its four raw corners. The no-overshoot bound is what makes bilinear the right
choice at the soot boundary, where `γ` steps ~4.5 % in one interval.

**Refusal to extrapolate.** All five edge classifications are correct and all four
out-of-range evaluations raise. Nothing is clamped.

Nominal case:

| Quantity | Value |
| --- | --- |
| `O/F` | 2.5569 → 1.8560 |
| Delivered `c*` | 1356.21 → 1279.81 m/s (−5.633 %) |
| `γ` | 1.23163 → 1.23490 (+0.265 %) |
| `T_c` | 2186.0 → 1925.8 K |
| `p_c` | 24.671 → 21.846 bar |
| Thrust | 290.760 → 253.208 N (**−12.915 %**) |
| Total impulse vs M4 | **−10.689 %** |

At the M1/M3 reference flows the tabulated properties differ from the Milestone 3
prescribed constants by: ideal `c*` **−7.84 %**, `T_c` **−16.45 %**, `γ` +2.67 %,
molar mass −4.29 %, `R` +4.48 %. (Both evaluated at the converged variable-property
chamber pressure, which is the convention used throughout; evaluating at Milestone
3's frozen `p_c` instead gives −16.43 % for `T_c`.)

## M6 final robustness

**Sensitivity methodology.** Deterministic, one-factor-at-a-time. Each family
changes exactly one assumption to stated alternative values and reruns the entire
coupled model; everything else stays at the reference case. No probability
distribution is assumed for any input, and no sampling is performed.

**Ranking metric**, disclosed and computed from the runs rather than assigned:

```
swing(family, metric) = max | 100 (variant_metric / baseline_metric − 1) |
                        over the VALID variants of that family
```

Families with no valid variant are dropped rather than ranked at zero, which would
understate them as insensitive.

**Ranking, total impulse:** injector model 13.589 % · regression coefficient `a`
9.028 % · initial tank temperature 7.347 % · injector effective area 5.381 % ·
regression exponent `n` 2.606 % · `c*` efficiency 0.073 % · chemistry grid
resolution 0.038 % · chemistry pressure dependence 0.003 %.

**Ranking, thrust decay:** injector model 52.770 % · regression coefficient `a`
18.310 % · initial tank temperature 18.189 % · injector effective area 14.071 % ·
regression exponent `n` 9.308 % · `c*` efficiency 1.857 % · chemistry pressure
dependence 0.063 % · chemistry grid resolution 0.042 %.

**Ranking, burn duration:** regression coefficient `a` 11.401 % · injector model
6.845 % · initial tank temperature 4.320 % · regression exponent `n` 3.271 % ·
injector effective area 3.115 % · `c*` efficiency 0.297 % · chemistry grid
resolution 0.004 % · chemistry pressure dependence 0.000 %.

The dominant assumption is **not** the same for all three outputs — burn duration
is led by the regression coefficient while impulse and slope are led by the
injector model — so no single universal ranking is claimed.

**Invalid case, reported not repaired.** Driving the coupling with the bare **HEM**
limit produces an oxidizer flow too low for the chamber to reach the table's 10 bar
floor, so no firing state is ever established. It is recorded with its terminal
status, excluded from every ranking, and **not clamped back into validity**. All
other 14 perturbed cases stayed inside both the tabulated rectangle and the checked
N₂O property band.

**Qualitative conclusion support counts:**

| # | Conclusion | Support | Verdict |
| --- | --- | --- | --- |
| 1 | Coupling the tank reverses the M3 rising-thrust trend | 5/5 | robust within tested deterministic family |
| 2 | Variable thermochemistry preserves the falling-thrust sign | 15/15 | robust within tested deterministic family |
| 3 | Variable thermochemistry deepens the decay relative to M4 | 5/5 | robust within tested deterministic family |
| 4 | `c*` variation dominates the M5 chemistry effect over `γ` | 15/15 | robust within tested deterministic family |
| 5 | Feed coupling changes thrust shape more than the chemistry does | 5/5 | robust within tested deterministic family |
| 6 | `O/F` drifts downward during the burn | 15/15 | robust within tested deterministic family |
| 7 | The 10 L reference case burns through before liquid runs out | 15/15 | robust within tested deterministic family |
| 8 | A small tank can instead terminate on liquid depletion | 1/1 | robust within tested deterministic family |

"Robust within tested deterministic family" means every case in **this** tested set
agreed. It is not a claim about untested regions of the assumption space, and it is
not a probability.

Conclusion 4 in numbers: over the nominal burn delivered `c*` moves 5.633 % while
`γ` moves 0.265 % — a factor of **21.2**.

### One defect found and fixed this milestone

The audit drove the coupled solver into a region no published case reaches, and
exposed a **diagnostic** defect in `variable_feed_system.py`.

* **Demonstrated independently.** With the bare HEM injector the feed root sits at
  `ṁ_ox` ≈ 0.0346 kg/s, giving `O/F` ≈ 1.2845 — comfortably inside the tabulated
  span [1.20, 4.00]. The binding constraint is the chamber-pressure floor: the
  lower search edge set by the O/F bound has a pressure margin of −148 174 Pa, so
  the pressure floor raises that edge to `ṁ_ox` = 0.037579 kg/s.
* **Consequence.** The solver reported `OF_OUTSIDE_TABLE` for **every** root below
  the window regardless of which bound trimmed it, pointing a user at the wrong
  assumption. Diagnostic only: every Milestone 4 and Milestone 5 study case is
  `FLOWING`, so no published number was ever affected.
* **Fixed minimally.** The trimming loop now records which constraint set each
  edge and reports that reason. The complementary upper-edge case is handled the
  same way, and `NO_ROOT` is kept for the genuinely different situation where the
  injector's free-discharge flow — not the table — is what bounds the window.
* **Regression tests added.** Three, in `tests/test_variable_feed_system.py`: the
  HEM case must report `PRESSURE_OUTSIDE_TABLE`; a starved injector must still
  report `OF_OUTSIDE_TABLE` (so the fix did not simply relabel everything); and the
  nominal operating point must be unchanged to 1e-12 relative.
* **Verified no published output moved.** `scripts/variable_thermochemistry_study.py`
  produces **byte-identical** output before and after the fix, and the nominal
  coupled solve still gives `ṁ_ox` = 0.10270411484637884 kg/s and
  `p_c` = 2 467 088.5060716257 Pa.

## Regression — all frozen headline values

Reproduced this session, 40 of 40 checks passing.

| Milestone | Quantity | Frozen value | This session |
| --- | --- | --- | --- |
| M1 | `G_ox` | 79.57747 | 79.577472 |
| M1 | `ṙ` | 0.850894 mm/s | 0.850894 |
| M1 | `ṁ_f` | 0.039776644 kg/s | 0.039776644 |
| M1 | `O/F` | 2.514038 | 2.514038 |
| M2 | Burnout | 41.740618 s | 41.740618 |
| M2 | Port diameter | 40 → 90 mm | 40 → 90 |
| M2 | Fuel closure | ≈ 8.3e-11 kg | 7.943e-11 |
| M3 | `p_c` | 26.114 bar | 26.114248 |
| M3 | Thrust | 310.048 N | 310.047929 |
| M3 | `I_sp` | 226.190 s | 226.190078 |
| M3 | Total impulse | 13 525.65 N·s | 13 525.651216 |
| M3 | Equivalent `I_sp` | 227.103 s | 227.103159 |
| M4 | Initial tank pressure | 50.53 bar | 50.525093 |
| M4 | Thrust | 313.7 → 286.8 N | 313.733 → 286.756 |
| M4 | Thrust slope | −8.6 % | −8.598864 |
| M4 | Burnout | 42.917 s | 42.917324 |
| M4 | Final tank pressure | 38.09 bar | 38.093619 |
| M4 | Final tank temperature | 281.12 K | 281.117427 |
| M4 | `O/F` | 2.53 → 1.81 | 2.534935 → 1.812801 |
| M4 | 2.5 L liquid depletion | ≈ 16.0 s | 16.035367 |
| M4 | SPI initial thrust | 382.0 N | 381.966849 (**+21.75 %**, see correction above) |
| M5 | Thrust decay | −12.915 % | −12.914915 |
| M5 | Decay ratio vs M4 | 1.50 | 1.501933 |
| M5 | Impulse vs M4 | −10.689 % | −10.688801 |
| M5 | `O/F` | 2.557 → 1.856 | 2.556881 → 1.855960 |
| M5 | Reference ideal `c*` vs M3 | −7.84 % | −7.840636 |
| M5 | Reference `T_c` vs M3 | −16.45 % | −16.445027 |

One rounding correction was made to the documentation this milestone: the M5
thrust decay is −12.914915 %, which rounds to **−12.91 %** at two decimals. The
Milestone 5 README reported −12.92 %. The underlying value never changed.

## Final test suite

```
pytest -W error -q      857 passed, 0 failed, 0 warnings
ruff check .            All checks passed!
```

857 tests: 167 Milestone 1, 119 Milestone 2, 215 Milestone 3, 211 Milestone 4,
142 Milestone 5, 3 Milestone 6 regression tests.

`-W error` means any warning fails the suite; none is raised. All 11 primary
scripts exit 0, and all 29 figures regenerate byte-identically within one
environment.
