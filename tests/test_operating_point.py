"""Independent verification of prescribed-oxidizer operating-point bookkeeping.

Expected values are formed from independent algebra, e.g. ``G_ox = 4 m_dot / (pi D^2)``
rather than ``m_dot / A_port``, so that a mistake in the production port-area
routine cannot cancel itself out.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from hybrid_rocket_motor.geometry import REPRESENTATIVE_GRAIN, GrainGeometry
from hybrid_rocket_motor.operating_point import (
    OperatingPoint,
    fuel_mass_flow,
    mixture_ratio,
    oxidizer_mass_flux,
)
from hybrid_rocket_motor.regression import ILLUSTRATIVE_SENSITIVITY_LAWS, REZAEI_2018_N2O_HTPB

LAW = REZAEI_2018_N2O_HTPB

#: Representative prescribed oxidizer mass flow for the Milestone 1 study [kg/s].
REPRESENTATIVE_M_DOT_OX = 0.100


# --- 6. G_ox = m_dot_ox / A_port ------------------------------------------


@pytest.mark.parametrize(("m_dot", "diameter"), [(0.05, 0.040), (0.10, 0.040), (0.10, 0.070)])
def test_flux_matches_independent_four_over_pi_d_squared_formula(m_dot, diameter):
    """G_ox = 4 m_dot / (pi D^2), formed without calling port_area_m2."""
    point = OperatingPoint(REPRESENTATIVE_GRAIN, m_dot, diameter)
    expected = 4.0 * m_dot / (math.pi * diameter**2)
    assert point.oxidizer_mass_flux_si == pytest.approx(expected, rel=1e-12)


def test_flux_matches_hardcoded_literal_at_the_representative_point():
    """0.100 kg/s through pi(0.02)^2 m^2 is 79.57747154594766 kg/(m^2 s)."""
    point = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX)
    assert point.oxidizer_mass_flux_si == pytest.approx(79.57747154594766, rel=1e-12)


def test_free_function_and_operating_point_agree():
    point = OperatingPoint(REPRESENTATIVE_GRAIN, 0.08, 0.05)
    assert oxidizer_mass_flux(0.08, point.port_area_m2) == pytest.approx(
        point.oxidizer_mass_flux_si, rel=1e-12
    )


def test_flux_is_proportional_to_oxidizer_mass_flow():
    single = OperatingPoint(REPRESENTATIVE_GRAIN, 0.05).oxidizer_mass_flux_si
    double = OperatingPoint(REPRESENTATIVE_GRAIN, 0.10).oxidizer_mass_flux_si
    assert double == pytest.approx(2.0 * single, rel=1e-12)


def test_zero_oxidizer_flow_gives_zero_flux():
    assert OperatingPoint(REPRESENTATIVE_GRAIN, 0.0).oxidizer_mass_flux_si == 0.0


# --- 9. m_dot_f = rho_f * A_burn * r_dot ----------------------------------


def test_fuel_mass_flow_matches_independent_formula():
    """m_dot_f = rho * pi D L * r_dot, assembled from scratch."""
    point = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX)
    r_dot = point.regression_rate_si(LAW)
    expected = (
        REPRESENTATIVE_GRAIN.fuel_density_kg_m3
        * math.pi
        * REPRESENTATIVE_GRAIN.initial_port_diameter_m
        * REPRESENTATIVE_GRAIN.length_m
        * r_dot
    )
    assert point.fuel_mass_flow_kg_s(LAW) == pytest.approx(expected, rel=1e-12)


def test_fuel_mass_flow_matches_hardcoded_literal():
    """930 * 0.05026548245743669 * 8.508937503295065e-4 = 0.039776643938707 kg/s."""
    point = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX)
    assert point.fuel_mass_flow_kg_s(LAW) == pytest.approx(0.039776643938707, rel=1e-12)


def test_fuel_mass_flow_free_function_against_exact_integers():
    """1000 kg/m^3 * 2 m^2 * 0.003 m/s = 6 kg/s exactly."""
    assert fuel_mass_flow(1000.0, 2.0, 0.003) == pytest.approx(6.0, rel=1e-14)


def test_fuel_mass_flow_is_linear_in_each_factor():
    assert fuel_mass_flow(2000.0, 2.0, 0.003) == pytest.approx(12.0, rel=1e-14)
    assert fuel_mass_flow(1000.0, 4.0, 0.003) == pytest.approx(12.0, rel=1e-14)
    assert fuel_mass_flow(1000.0, 2.0, 0.006) == pytest.approx(12.0, rel=1e-14)


def test_fuel_mass_flow_is_positive_for_positive_regression():
    point = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX)
    assert point.regression_rate_si(LAW) > 0.0
    assert point.fuel_mass_flow_kg_s(LAW) > 0.0


def test_zero_regression_rate_gives_zero_fuel_flow():
    assert fuel_mass_flow(930.0, 0.05, 0.0) == 0.0


@pytest.mark.parametrize(
    "args",
    [
        (0.0, 2.0, 0.003),
        (-930.0, 2.0, 0.003),
        (1000.0, 0.0, 0.003),
        (1000.0, -2.0, 0.003),
        (1000.0, 2.0, -0.003),
        (math.nan, 2.0, 0.003),
        (1000.0, math.inf, 0.003),
        (1000.0, 2.0, math.nan),
    ],
)
def test_fuel_mass_flow_rejects_invalid_inputs(args):
    with pytest.raises(ValueError):
        fuel_mass_flow(*args)


# --- 10. O/F ---------------------------------------------------------------


def test_mixture_ratio_against_exact_integers():
    assert mixture_ratio(6.0, 2.0) == pytest.approx(3.0, rel=1e-14)
    assert mixture_ratio(1.0, 4.0) == pytest.approx(0.25, rel=1e-14)


def test_mixture_ratio_matches_independent_ratio_at_representative_point():
    point = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX)
    expected = REPRESENTATIVE_M_DOT_OX / point.fuel_mass_flow_kg_s(LAW)
    assert point.mixture_ratio(LAW) == pytest.approx(expected, rel=1e-12)


def test_mixture_ratio_matches_hardcoded_literal():
    """0.100 / 0.039776643938707 = 2.5140381414302 (dimensionless)."""
    point = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX)
    assert point.mixture_ratio(LAW) == pytest.approx(2.5140381414302, rel=1e-12)


def test_mixture_ratio_is_finite_and_positive_across_the_study_range():
    for m_dot in np.linspace(0.045, 0.150, 8):
        point = OperatingPoint(REPRESENTATIVE_GRAIN, float(m_dot))
        of = point.mixture_ratio(LAW)
        assert math.isfinite(of)
        assert of > 0.0


def test_mixture_ratio_is_undefined_without_fuel_flow():
    with pytest.raises(ValueError):
        mixture_ratio(0.1, 0.0)
    with pytest.raises(ValueError):
        mixture_ratio(0.1, -1e-9)


def test_mixture_ratio_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        mixture_ratio(-0.1, 0.04)
    with pytest.raises(ValueError):
        mixture_ratio(math.nan, 0.04)
    with pytest.raises(ValueError):
        mixture_ratio(0.1, math.nan)


# --- 11. zero oxidizer-flow behaviour -------------------------------------


def test_zero_oxidizer_flow_gives_zero_regression_and_zero_fuel_flow():
    point = OperatingPoint(REPRESENTATIVE_GRAIN, 0.0)
    assert point.oxidizer_mass_flux_si == 0.0
    assert point.regression_rate_si(LAW) == 0.0
    assert point.fuel_mass_flow_kg_s(LAW) == 0.0


def test_mixture_ratio_is_undefined_at_zero_oxidizer_flow():
    point = OperatingPoint(REPRESENTATIVE_GRAIN, 0.0)
    with pytest.raises(ValueError):
        point.mixture_ratio(LAW)


# --- 12. negative-input rejection -----------------------------------------


@pytest.mark.parametrize("bad_flow", [-1e-9, -0.1, math.nan, math.inf])
def test_negative_or_nonfinite_oxidizer_flow_is_rejected(bad_flow):
    with pytest.raises(ValueError):
        OperatingPoint(REPRESENTATIVE_GRAIN, bad_flow)


@pytest.mark.parametrize("bad_area", [0.0, -1.0, math.nan])
def test_flux_rejects_invalid_port_area(bad_area):
    with pytest.raises(ValueError):
        oxidizer_mass_flux(0.1, bad_area)


@pytest.mark.parametrize("bad_flow", [-0.1, math.nan])
def test_flux_rejects_invalid_mass_flow(bad_flow):
    with pytest.raises(ValueError):
        oxidizer_mass_flux(bad_flow, 1e-3)


# --- 13. invalid geometry rejection at the operating point ----------------


@pytest.mark.parametrize("bad_port", [0.0, -0.01, 0.090, 0.20, math.nan])
def test_operating_point_rejects_invalid_port_diameter(bad_port):
    with pytest.raises(ValueError):
        OperatingPoint(REPRESENTATIVE_GRAIN, 0.1, bad_port)


def test_operating_point_defaults_to_the_initial_port_diameter():
    point = OperatingPoint(REPRESENTATIVE_GRAIN, 0.1)
    assert point.port_diameter_m == REPRESENTATIVE_GRAIN.initial_port_diameter_m


def test_operating_point_is_immutable():
    point = OperatingPoint(REPRESENTATIVE_GRAIN, 0.1)
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        point.oxidizer_mass_flow_kg_s = 0.2


# --- 16. growing the port at fixed m_dot_ox lowers G_ox and r_dot ---------

PORT_DIAMETERS_M = (0.040, 0.050, 0.060, 0.070)


def test_larger_port_reduces_flux_at_fixed_oxidizer_flow():
    fluxes = [
        OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX, d).oxidizer_mass_flux_si
        for d in PORT_DIAMETERS_M
    ]
    assert all(b < a for a, b in itertools.pairwise(fluxes))


def test_larger_port_reduces_regression_rate_at_fixed_oxidizer_flow():
    rates = [
        OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX, d).regression_rate_si(LAW)
        for d in PORT_DIAMETERS_M
    ]
    assert all(b < a for a, b in itertools.pairwise(rates))
    assert all(r > 0.0 for r in rates)


def test_port_size_coupling_holds_for_every_law_including_illustrative_ones():
    for law in (LAW, *ILLUSTRATIVE_SENSITIVITY_LAWS):
        rates = [
            OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX, d).regression_rate_si(law)
            for d in PORT_DIAMETERS_M
        ]
        assert all(b < a for a, b in itertools.pairwise(rates)), law.label


def test_flux_falls_exactly_as_the_inverse_square_of_port_diameter():
    """Doubling D_p at fixed m_dot_ox must quarter G_ox."""
    small = OperatingPoint(REPRESENTATIVE_GRAIN, 0.1, 0.040).oxidizer_mass_flux_si
    large = OperatingPoint(REPRESENTATIVE_GRAIN, 0.1, 0.080).oxidizer_mass_flux_si
    assert large == pytest.approx(small / 4.0, rel=1e-12)


def test_larger_port_increases_burning_area_but_still_lowers_regression_rate():
    """The two competing effects on m_dot_f are both present and correctly signed."""
    small = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX, 0.040)
    large = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX, 0.070)
    assert large.burning_area_m2 > small.burning_area_m2
    assert large.regression_rate_si(LAW) < small.regression_rate_si(LAW)


# --- engineering sanity ----------------------------------------------------


def test_remaining_fuel_mass_is_positive_across_the_port_sensitivity_range():
    for d in PORT_DIAMETERS_M:
        point = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX, d)
        assert point.remaining_fuel_mass_kg > 0.0
        assert point.port_area_m2 > 0.0


def test_no_nan_or_inf_anywhere_in_the_study_domain():
    for m_dot in np.linspace(0.045, 0.150, 8):
        for d in PORT_DIAMETERS_M:
            point = OperatingPoint(REPRESENTATIVE_GRAIN, float(m_dot), d)
            values = (
                point.port_area_m2,
                point.port_perimeter_m,
                point.burning_area_m2,
                point.oxidizer_mass_flux_si,
                point.regression_rate_si(LAW),
                point.fuel_mass_flow_kg_s(LAW),
                point.mixture_ratio(LAW),
                point.remaining_fuel_mass_kg,
            )
            assert all(math.isfinite(v) for v in values)
            assert all(v > 0.0 for v in values)


def test_representative_flux_sits_inside_the_sources_validated_range():
    point = OperatingPoint(REPRESENTATIVE_GRAIN, REPRESENTATIVE_M_DOT_OX)
    assert LAW.is_within_validated_flux_range(point.oxidizer_mass_flux_si)


def test_geometry_reference_is_carried_by_the_operating_point():
    other = GrainGeometry(
        length_m=0.25,
        initial_port_diameter_m=0.030,
        outer_diameter_m=0.080,
        fuel_density_kg_m3=930.0,
    )
    point = OperatingPoint(other, 0.1)
    assert point.geometry is other
    assert point.burning_area_m2 == pytest.approx(math.pi * 0.030 * 0.25, rel=1e-12)
