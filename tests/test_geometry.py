"""Independent verification of the grain-geometry module.

Expected values are written from first principles using algebra deliberately
different from the production code (radius-based rather than diameter-based
forms) and, where practical, hard-coded decimal literals.  No test calls the
function under test to build its own expectation.
"""

from __future__ import annotations

import itertools
import math

import pytest

from hybrid_rocket_motor.geometry import REPRESENTATIVE_GRAIN, GrainGeometry

# A deliberately simple grain whose exact answers are clean multiples of pi.
SIMPLE = GrainGeometry(
    length_m=2.0,
    initial_port_diameter_m=1.0,
    outer_diameter_m=3.0,
    fuel_density_kg_m3=1000.0,
)


# --- 1. circular port area -------------------------------------------------


def test_port_area_matches_independent_radius_formula():
    """A_port from pi r^2, not pi d^2 / 4."""
    for d in (0.5, 1.0, 2.5):
        radius = d / 2.0
        assert SIMPLE.port_area_m2(d) == pytest.approx(math.pi * radius**2, rel=1e-12)


def test_port_area_matches_hardcoded_literal():
    """A_port for D_p = 1 m is pi/4 = 0.78539816339744830961... m^2."""
    assert SIMPLE.port_area_m2(1.0) == pytest.approx(0.7853981633974483, rel=1e-12)


def test_port_area_scales_as_diameter_squared():
    """Doubling the port diameter must quadruple the area."""
    assert SIMPLE.port_area_m2(2.0) == pytest.approx(4.0 * SIMPLE.port_area_m2(1.0), rel=1e-12)


# --- 2. burning perimeter --------------------------------------------------


def test_port_perimeter_matches_independent_radius_formula():
    """P_port from 2 pi r, not pi d."""
    for d in (0.5, 1.0, 2.5):
        assert SIMPLE.port_perimeter_m(d) == pytest.approx(2.0 * math.pi * (d / 2.0), rel=1e-12)


def test_port_perimeter_matches_hardcoded_literal():
    """P_port for D_p = 1 m is pi = 3.14159265358979311600 m."""
    assert SIMPLE.port_perimeter_m(1.0) == pytest.approx(3.141592653589793, rel=1e-12)


# --- 3. burning area -------------------------------------------------------


def test_burning_area_matches_independent_formula():
    """A_burn = 2 pi r L, built without using port_perimeter_m."""
    for d in (0.5, 1.0, 2.5):
        expected = 2.0 * math.pi * (d / 2.0) * SIMPLE.length_m
        assert SIMPLE.burning_area_m2(d) == pytest.approx(expected, rel=1e-12)


def test_burning_area_matches_hardcoded_literal():
    """A_burn for D_p = 1 m, L = 2 m is 2 pi = 6.28318530717958623200 m^2."""
    assert SIMPLE.burning_area_m2(1.0) == pytest.approx(6.283185307179586, rel=1e-12)


def test_burning_area_is_linear_in_length():
    longer = GrainGeometry(
        length_m=4.0,
        initial_port_diameter_m=1.0,
        outer_diameter_m=3.0,
        fuel_density_kg_m3=1000.0,
    )
    assert longer.burning_area_m2(1.0) == pytest.approx(
        2.0 * SIMPLE.burning_area_m2(1.0), rel=1e-12
    )


# --- 4. remaining fuel volume ---------------------------------------------


def test_fuel_cross_section_is_outer_disc_minus_port_disc():
    """Annulus area built as a difference of two independently formed discs."""
    for d in (0.5, 1.0, 2.5):
        outer_disc = math.pi * (SIMPLE.outer_diameter_m / 2.0) ** 2
        port_disc = math.pi * (d / 2.0) ** 2
        assert SIMPLE.fuel_cross_section_area_m2(d) == pytest.approx(
            outer_disc - port_disc, rel=1e-12
        )


def test_fuel_volume_matches_independent_formula():
    """V = (pi/4)(D_o^2 - D_p^2) L, arranged differently from production code."""
    for d in (0.5, 1.0, 2.5):
        expected = 0.25 * math.pi * (SIMPLE.outer_diameter_m**2 - d**2) * SIMPLE.length_m
        assert SIMPLE.fuel_volume_m3(d) == pytest.approx(expected, rel=1e-12)


def test_fuel_volume_matches_hardcoded_literal():
    """V for D_o = 3, D_p = 1, L = 2 is 4 pi = 12.5663706143591724640 m^3."""
    assert SIMPLE.fuel_volume_m3(1.0) == pytest.approx(12.566370614359172, rel=1e-12)


def test_fuel_volume_decreases_as_port_grows():
    volumes = [SIMPLE.fuel_volume_m3(d) for d in (1.0, 1.5, 2.0, 2.5)]
    assert all(later < earlier for earlier, later in itertools.pairwise(volumes))


# --- 5. remaining fuel mass -----------------------------------------------


def test_fuel_mass_matches_independent_formula():
    """m = rho * (pi/4)(D_o^2 - D_p^2) L."""
    for d in (0.5, 1.0, 2.5):
        expected = (
            SIMPLE.fuel_density_kg_m3
            * 0.25
            * math.pi
            * (SIMPLE.outer_diameter_m**2 - d**2)
            * SIMPLE.length_m
        )
        assert SIMPLE.fuel_mass_kg(d) == pytest.approx(expected, rel=1e-12)


def test_fuel_mass_matches_hardcoded_literal():
    """m = 1000 * 4 pi = 12566.3706143591724640 kg for the simple grain."""
    assert SIMPLE.fuel_mass_kg(1.0) == pytest.approx(12566.370614359172, rel=1e-12)


def test_fuel_mass_is_proportional_to_density():
    denser = GrainGeometry(
        length_m=2.0,
        initial_port_diameter_m=1.0,
        outer_diameter_m=3.0,
        fuel_density_kg_m3=2000.0,
    )
    assert denser.fuel_mass_kg(1.0) == pytest.approx(2.0 * SIMPLE.fuel_mass_kg(1.0), rel=1e-12)


def test_fuel_mass_tends_to_zero_as_port_approaches_outer_diameter():
    nearly_burnt = SIMPLE.fuel_mass_kg(2.999999)
    assert 0.0 < nearly_burnt < 1e-2 * SIMPLE.fuel_mass_kg(1.0)


# --- 13. invalid geometry rejection ---------------------------------------


@pytest.mark.parametrize(
    "kwargs",
    [
        {"length_m": 0.0},
        {"length_m": -0.4},
        {"initial_port_diameter_m": 0.0},
        {"initial_port_diameter_m": -0.04},
        {"outer_diameter_m": 0.0},
        {"outer_diameter_m": -0.09},
        {"fuel_density_kg_m3": 0.0},
        {"fuel_density_kg_m3": -930.0},
        {"length_m": math.nan},
        {"fuel_density_kg_m3": math.inf},
    ],
)
def test_nonpositive_or_nonfinite_dimensions_are_rejected(kwargs):
    base = {
        "length_m": 0.40,
        "initial_port_diameter_m": 0.040,
        "outer_diameter_m": 0.090,
        "fuel_density_kg_m3": 930.0,
    }
    base.update(kwargs)
    with pytest.raises(ValueError):
        GrainGeometry(**base)


@pytest.mark.parametrize("port", [0.090, 0.100])
def test_port_not_smaller_than_outer_diameter_is_rejected(port):
    with pytest.raises(ValueError):
        GrainGeometry(
            length_m=0.40,
            initial_port_diameter_m=port,
            outer_diameter_m=0.090,
            fuel_density_kg_m3=930.0,
        )


@pytest.mark.parametrize("port", [0.0, -0.01, 0.090, 0.15, math.nan, math.inf])
def test_invalid_current_port_diameter_is_rejected_by_every_query(port):
    for method in (
        REPRESENTATIVE_GRAIN.port_area_m2,
        REPRESENTATIVE_GRAIN.port_perimeter_m,
        REPRESENTATIVE_GRAIN.burning_area_m2,
        REPRESENTATIVE_GRAIN.fuel_cross_section_area_m2,
        REPRESENTATIVE_GRAIN.fuel_volume_m3,
        REPRESENTATIVE_GRAIN.fuel_mass_kg,
    ):
        with pytest.raises(ValueError):
            method(port)


def test_geometry_is_immutable():
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        REPRESENTATIVE_GRAIN.length_m = 1.0


# --- representative grain sanity ------------------------------------------


def test_representative_grain_has_expected_illustrative_dimensions():
    assert REPRESENTATIVE_GRAIN.length_m == 0.40
    assert REPRESENTATIVE_GRAIN.initial_port_diameter_m == 0.040
    assert REPRESENTATIVE_GRAIN.outer_diameter_m == 0.090
    assert REPRESENTATIVE_GRAIN.fuel_density_kg_m3 == 930.0


def test_representative_grain_initial_values_against_hand_computed_literals():
    # A_port = pi (0.02)^2                      = 1.256637061435917e-3 m^2
    # A_burn = pi (0.040)(0.40)                 = 5.026548245743669e-2 m^2
    # V      = pi/4 (0.090^2 - 0.040^2)(0.40)   = 2.042035224833366e-3 m^3
    # m      = 930 * V                          = 1.899092759095030    kg
    assert REPRESENTATIVE_GRAIN.initial_port_area_m2 == pytest.approx(
        1.256637061435918e-3, rel=1e-12
    )
    assert REPRESENTATIVE_GRAIN.initial_burning_area_m2 == pytest.approx(
        5.026548245743669e-2, rel=1e-12
    )
    assert REPRESENTATIVE_GRAIN.fuel_volume_m3(0.040) == pytest.approx(
        2.042035224833366e-3, rel=1e-12
    )
    assert REPRESENTATIVE_GRAIN.initial_fuel_mass_kg == pytest.approx(1.899092759095030, rel=1e-12)


def test_representative_initial_web_thickness():
    """(0.090 - 0.040)/2 = 0.025 m."""
    assert REPRESENTATIVE_GRAIN.initial_web_thickness_m == pytest.approx(0.025, rel=1e-12)


def test_representative_grain_has_positive_area_and_mass():
    assert REPRESENTATIVE_GRAIN.initial_port_area_m2 > 0.0
    assert REPRESENTATIVE_GRAIN.initial_fuel_mass_kg > 0.0
