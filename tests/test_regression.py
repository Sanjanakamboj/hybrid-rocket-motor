"""Independent verification of the power-law regression module.

The strongest checks here are:

* the SI unit conversion is re-derived by an *independent route* (convert the SI
  flux back into the source's g/(cm^2 s), evaluate the printed correlation in
  mm/s, then convert the answer to m/s) rather than by reusing the production
  code's pre-collapsed SI coefficient; and
* the production code is required to reproduce the *published experimental data*
  of the source paper, in the source's own units, to within the fit scatter the
  paper itself displays.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from hybrid_rocket_motor.regression import (
    ILLUSTRATIVE_ANCHOR_FLUX_SI,
    ILLUSTRATIVE_SENSITIVITY_LAWS,
    REZAEI_2018_N2O_HTPB,
    MassFluxUnit,
    PowerLawRegressionLaw,
    Provenance,
    RegressionRateUnit,
    illustrative_exponent_variant,
)

# Coefficients exactly as printed in Rezaei et al. (2018), Eq. (10):
#     r_dot [mm/s] = 0.3977 * (G_ox [g/(cm^2 s)]) ** 0.3667
SOURCE_A_MM_S = 0.3977
SOURCE_N = 0.3667


def expected_rate_si_via_source_units(flux_si: float) -> float:
    """Independent reference implementation.

    Deliberately takes the long way round: SI flux -> g/(cm^2 s) -> the printed
    correlation -> mm/s -> m/s.  It never touches ``coefficient_si``.
    """
    flux_g_cm2_s = flux_si / 10.0  # 1 g/(cm^2 s) == 10 kg/(m^2 s)
    rate_mm_s = SOURCE_A_MM_S * flux_g_cm2_s**SOURCE_N
    return rate_mm_s / 1000.0


# --- 7. regression power law ----------------------------------------------


@pytest.mark.parametrize("flux_si", [1.0, 10.0, 35.0, 79.577471545947, 120.0, 250.0])
def test_regression_rate_matches_independent_unit_route(flux_si):
    assert REZAEI_2018_N2O_HTPB.regression_rate_si(flux_si) == pytest.approx(
        expected_rate_si_via_source_units(flux_si), rel=1e-12
    )


def test_regression_rate_matches_hardcoded_literal_at_representative_flux():
    """At G_ox = 79.57747154594766 kg/(m^2 s) the law gives 0.8508937503295065 mm/s.

    Cross-checked by hand: G = 7.957747154594766 g/(cm^2 s);
    0.3977 * 7.957747154594766^0.3667 = 0.85089375... mm/s.
    """
    rate = REZAEI_2018_N2O_HTPB.regression_rate_si(79.57747154594766)
    assert rate == pytest.approx(8.508937503295065e-4, rel=1e-12)


def test_power_law_obeys_its_own_scaling_identity():
    """r(k G) / r(G) must equal k^n exactly, for any k -- a law-shape check."""
    law = REZAEI_2018_N2O_HTPB
    base = law.regression_rate_si(40.0)
    for k in (0.5, 2.0, 3.0):
        assert law.regression_rate_si(40.0 * k) / base == pytest.approx(k**law.exponent, rel=1e-12)


# --- 8. coefficient / unit conversion -------------------------------------


def test_si_coefficient_equals_independently_derived_value():
    """a_SI = 1e-3 * 0.3977 * 10^(-0.3667) = 1.7094470...e-4 (m/s)/(kg/m^2 s)^n."""
    independent = 1.0e-3 * SOURCE_A_MM_S * math.pow(10.0, -SOURCE_N)
    assert REZAEI_2018_N2O_HTPB.coefficient_si == pytest.approx(independent, rel=1e-14)
    assert REZAEI_2018_N2O_HTPB.coefficient_si == pytest.approx(1.709447e-4, rel=1e-6)


def test_unit_conversion_factors_are_correct():
    assert MassFluxUnit.G_PER_CM2_S.si_per_unit == 10.0
    assert MassFluxUnit.KG_PER_M2_S.si_per_unit == 1.0
    assert RegressionRateUnit.MM_PER_S.si_per_unit == 1.0e-3
    assert RegressionRateUnit.M_PER_S.si_per_unit == 1.0


def test_declaring_the_law_directly_in_si_gives_the_same_curve():
    """Same physics declared in SI units must produce an identical curve."""
    in_si = PowerLawRegressionLaw(
        coefficient_source_units=1.0e-3 * SOURCE_A_MM_S * math.pow(10.0, -SOURCE_N),
        exponent=SOURCE_N,
        flux_unit=MassFluxUnit.KG_PER_M2_S,
        rate_unit=RegressionRateUnit.M_PER_S,
        label="same law, SI declaration",
        provenance=Provenance.SOURCED,
        reference="unit-conversion round-trip check",
    )
    for flux in (5.0, 50.0, 500.0):
        assert in_si.regression_rate_si(flux) == pytest.approx(
            REZAEI_2018_N2O_HTPB.regression_rate_si(flux), rel=1e-12
        )


def test_mm_per_s_helper_is_exactly_1000x_the_si_value():
    for flux in (0.0, 12.5, 80.0):
        si = REZAEI_2018_N2O_HTPB.regression_rate_si(flux)
        assert REZAEI_2018_N2O_HTPB.regression_rate_mm_s(flux) == pytest.approx(
            si * 1000.0, rel=1e-12
        )


# --- source-data consistency (strongest external check) --------------------

# (mean oxidizer mass flux [g/(cm^2 s)], measured regression rate [mm/s])
# transcribed from Table 4 of Rezaei et al. (2018) -- the 18 firings used to fit
# Eq. (10).  These are experimental measurements, not model output.
REZAEI_TABLE_4 = [
    (6.88, 0.779),
    (11.95, 0.891),
    (11.61, 1.044),
    (7.70, 0.882),
    (6.36, 0.775),
    (5.93, 0.810),
    (4.65, 0.653),
    (4.89, 0.740),
    (11.61, 1.006),
    (8.79, 0.940),
    (5.47, 0.710),
    (12.02, 0.899),
    (5.12, 0.693),
    (6.19, 0.739),
    (8.64, 0.953),
    (10.20, 0.948),
    (4.66, 0.732),
    (3.24, 0.613),
]


@pytest.mark.parametrize(("flux_g_cm2_s", "measured_mm_s"), REZAEI_TABLE_4)
def test_si_implementation_reproduces_published_measurements(flux_g_cm2_s, measured_mm_s):
    """SI code path must match the paper's own measured data to within its scatter.

    This exercises the whole unit chain end to end: a flux quoted in the source's
    g/(cm^2 s) is converted to SI, pushed through the SI coefficient, and the
    resulting m/s is compared against a mm/s laboratory measurement.  A silent
    factor-of-10 or factor-of-1000 slip anywhere would fail this outright.
    """
    flux_si = flux_g_cm2_s * 10.0
    predicted_mm_s = REZAEI_2018_N2O_HTPB.regression_rate_si(flux_si) * 1000.0
    assert predicted_mm_s == pytest.approx(measured_mm_s, rel=0.15)


def test_published_measurements_bracket_the_expected_millimetre_per_second_scale():
    """Sanity: HTPB/N2O regression is order 1 mm/s, not 1 m/s or 1 um/s."""
    for flux_g_cm2_s, _ in REZAEI_TABLE_4:
        rate_mm_s = REZAEI_2018_N2O_HTPB.regression_rate_si(flux_g_cm2_s * 10.0) * 1000.0
        assert 0.1 < rate_mm_s < 10.0


# --- 11. zero flux, and 12. negative-input rejection ----------------------


def test_zero_flux_gives_exactly_zero_regression_rate():
    assert REZAEI_2018_N2O_HTPB.regression_rate_si(0.0) == 0.0
    assert REZAEI_2018_N2O_HTPB.regression_rate_mm_s(0.0) == 0.0


def test_zero_flux_inside_an_array_gives_zero_without_nan():
    rates = REZAEI_2018_N2O_HTPB.regression_rate_si(np.array([0.0, 1.0, 50.0]))
    assert rates[0] == 0.0
    assert np.all(np.isfinite(rates))


def test_regression_rate_is_continuous_approaching_zero_flux():
    tiny = REZAEI_2018_N2O_HTPB.regression_rate_si(1e-12)
    assert 0.0 < tiny < 1e-6


@pytest.mark.parametrize("bad_flux", [-1e-9, -1.0, -100.0])
def test_negative_flux_is_rejected(bad_flux):
    with pytest.raises(ValueError):
        REZAEI_2018_N2O_HTPB.regression_rate_si(bad_flux)


def test_negative_flux_inside_an_array_is_rejected():
    with pytest.raises(ValueError):
        REZAEI_2018_N2O_HTPB.regression_rate_si(np.array([10.0, -10.0]))


@pytest.mark.parametrize("bad_flux", [math.nan, math.inf, -math.inf])
def test_non_finite_flux_is_rejected(bad_flux):
    with pytest.raises(ValueError):
        REZAEI_2018_N2O_HTPB.regression_rate_si(bad_flux)


@pytest.mark.parametrize(("a", "n"), [(0.0, 0.5), (-1.0, 0.5), (1.0, 0.0), (1.0, -0.5)])
def test_invalid_coefficients_are_rejected(a, n):
    with pytest.raises(ValueError):
        PowerLawRegressionLaw(
            coefficient_source_units=a,
            exponent=n,
            flux_unit=MassFluxUnit.KG_PER_M2_S,
            rate_unit=RegressionRateUnit.M_PER_S,
            label="bad",
            provenance=Provenance.ILLUSTRATIVE,
            reference="invalid-input test",
        )


@pytest.mark.parametrize("bad_range", [(0.0, 10.0), (-1.0, 10.0), (10.0, 10.0), (20.0, 10.0)])
def test_invalid_validity_ranges_are_rejected(bad_range):
    with pytest.raises(ValueError):
        PowerLawRegressionLaw(
            coefficient_source_units=1.0,
            exponent=0.5,
            flux_unit=MassFluxUnit.KG_PER_M2_S,
            rate_unit=RegressionRateUnit.M_PER_S,
            label="bad range",
            provenance=Provenance.ILLUSTRATIVE,
            reference="invalid-range test",
            valid_flux_range_si=bad_range,
        )


# --- 14. scalar / array consistency ---------------------------------------


def test_scalar_and_array_evaluation_agree_elementwise():
    fluxes = np.linspace(0.0, 200.0, 41)
    vectorised = REZAEI_2018_N2O_HTPB.regression_rate_si(fluxes)
    assert isinstance(vectorised, np.ndarray)
    assert vectorised.shape == fluxes.shape
    for flux, value in zip(fluxes.tolist(), vectorised.tolist(), strict=True):
        assert value == pytest.approx(REZAEI_2018_N2O_HTPB.regression_rate_si(flux), rel=1e-12)


def test_scalar_input_returns_a_plain_float():
    assert isinstance(REZAEI_2018_N2O_HTPB.regression_rate_si(50.0), float)
    assert isinstance(REZAEI_2018_N2O_HTPB.regression_rate_si(np.float64(50.0)), float)


def test_list_input_is_accepted_and_returns_an_array():
    result = REZAEI_2018_N2O_HTPB.regression_rate_si([10.0, 20.0])
    assert isinstance(result, np.ndarray)
    assert result.shape == (2,)


def test_array_evaluation_produces_no_nan_or_inf_over_the_study_domain():
    fluxes = np.linspace(0.0, 300.0, 601)
    rates = REZAEI_2018_N2O_HTPB.regression_rate_si(fluxes)
    assert np.all(np.isfinite(rates))
    assert np.all(rates >= 0.0)


# --- 15. monotonicity with mass flux --------------------------------------


def test_regression_rate_is_strictly_increasing_with_flux():
    fluxes = np.linspace(1.0, 200.0, 200)
    rates = REZAEI_2018_N2O_HTPB.regression_rate_si(fluxes)
    assert np.all(np.diff(rates) > 0.0)


def test_all_laws_are_strictly_increasing_with_flux():
    fluxes = np.linspace(1.0, 200.0, 100)
    for law in (REZAEI_2018_N2O_HTPB, *ILLUSTRATIVE_SENSITIVITY_LAWS):
        assert np.all(np.diff(law.regression_rate_si(fluxes)) > 0.0), law.label


def test_regression_rate_is_positive_for_any_positive_flux():
    for flux in (1e-6, 1.0, 42.0, 1000.0):
        assert REZAEI_2018_N2O_HTPB.regression_rate_si(flux) > 0.0


# --- validated-range reporting --------------------------------------------


def test_sourced_law_records_the_papers_reported_flux_range():
    """Table 4 covers 3.5-12 g/(cm^2 s), i.e. 35-120 kg/(m^2 s)."""
    assert REZAEI_2018_N2O_HTPB.valid_flux_range_si == (35.0, 120.0)
    assert REZAEI_2018_N2O_HTPB.is_within_validated_flux_range(80.0)
    assert not REZAEI_2018_N2O_HTPB.is_within_validated_flux_range(25.0)
    assert not REZAEI_2018_N2O_HTPB.is_within_validated_flux_range(150.0)


def test_a_law_without_a_recorded_range_never_claims_validity():
    for law in ILLUSTRATIVE_SENSITIVITY_LAWS:
        assert law.valid_flux_range_si is None
        assert not law.is_within_validated_flux_range(70.0)


# --- provenance labelling --------------------------------------------------


def test_sourced_law_is_marked_sourced_and_cites_its_paper():
    assert REZAEI_2018_N2O_HTPB.provenance is Provenance.SOURCED
    assert not REZAEI_2018_N2O_HTPB.is_illustrative
    assert "Scientia Iranica" in REZAEI_2018_N2O_HTPB.reference
    assert "10.24200/sci.2017.4317" in REZAEI_2018_N2O_HTPB.reference


def test_illustrative_laws_are_unmistakably_labelled():
    for law in ILLUSTRATIVE_SENSITIVITY_LAWS:
        assert law.provenance is Provenance.ILLUSTRATIVE
        assert law.is_illustrative
        assert "ILLUSTRATIVE" in law.label
        assert "ILLUSTRATIVE SENSITIVITY CASE" in law.reference
        assert "not experimental data" in law.reference


def test_illustrative_sweep_covers_the_reported_hybrid_exponent_band():
    """Karabeyoglu et al. (2007) report n in 0.5-0.8 for most hybrid systems."""
    exponents = sorted(law.exponent for law in ILLUSTRATIVE_SENSITIVITY_LAWS)
    assert exponents == [0.5, 0.62, 0.8]


def test_source_and_si_equation_strings_state_their_units():
    assert "mm/s" in REZAEI_2018_N2O_HTPB.source_equation
    assert "g/(cm^2 s)" in REZAEI_2018_N2O_HTPB.source_equation
    assert "m/s" in REZAEI_2018_N2O_HTPB.si_equation
    assert "kg/(m^2 s)" in REZAEI_2018_N2O_HTPB.si_equation


# --- illustrative variant construction ------------------------------------


def test_illustrative_variant_matches_base_law_exactly_at_the_anchor_flux():
    base = REZAEI_2018_N2O_HTPB
    for law in ILLUSTRATIVE_SENSITIVITY_LAWS:
        assert law.regression_rate_si(ILLUSTRATIVE_ANCHOR_FLUX_SI) == pytest.approx(
            base.regression_rate_si(ILLUSTRATIVE_ANCHOR_FLUX_SI), rel=1e-12
        )


def test_illustrative_variant_coefficient_matches_independent_pivot_algebra():
    """a_var = r_base(G_a) / G_a^n_var, computed here from the printed correlation."""
    anchor = ILLUSTRATIVE_ANCHOR_FLUX_SI
    rate_at_anchor = expected_rate_si_via_source_units(anchor)
    for law in ILLUSTRATIVE_SENSITIVITY_LAWS:
        expected_a = rate_at_anchor / anchor**law.exponent
        assert law.coefficient_si == pytest.approx(expected_a, rel=1e-12)


def test_higher_exponent_is_steeper_about_the_anchor():
    """A larger n must give a lower rate below the anchor and a higher rate above it."""
    low, high = ILLUSTRATIVE_SENSITIVITY_LAWS[0], ILLUSTRATIVE_SENSITIVITY_LAWS[-1]
    assert high.exponent > low.exponent
    below = 0.5 * ILLUSTRATIVE_ANCHOR_FLUX_SI
    above = 2.0 * ILLUSTRATIVE_ANCHOR_FLUX_SI
    assert high.regression_rate_si(below) < low.regression_rate_si(below)
    assert high.regression_rate_si(above) > low.regression_rate_si(above)


@pytest.mark.parametrize("bad_anchor", [0.0, -1.0, math.nan, math.inf])
def test_illustrative_variant_rejects_invalid_anchor(bad_anchor):
    with pytest.raises(ValueError):
        illustrative_exponent_variant(REZAEI_2018_N2O_HTPB, 0.6, bad_anchor)


@pytest.mark.parametrize("bad_exponent", [0.0, -0.5, math.nan, math.inf])
def test_illustrative_variant_rejects_invalid_exponent(bad_exponent):
    with pytest.raises(ValueError):
        illustrative_exponent_variant(REZAEI_2018_N2O_HTPB, bad_exponent)


def test_law_is_immutable():
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        REZAEI_2018_N2O_HTPB.exponent = 0.8
