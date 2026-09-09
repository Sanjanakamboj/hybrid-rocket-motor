"""Independent verification of the 1-D isentropic nozzle model.

Independence strategy
---------------------
Expected values are formed from the NASA Glenn relations written out here in a
different algebraic arrangement from the production code (for example the area
ratio is rebuilt from ``[(g+1)/2]^(-k)/M * [...]^k``, the form printed on
``astar.html``, rather than the ``(2/(g+1))`` form used in ``nozzle.py``), plus
hand-computed literals and closed-form inversions.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from hybrid_rocket_motor.nozzle import (
    NEAR_IDEAL_PRESSURE_TOLERANCE,
    REFERENCE_NOZZLE,
    SUMMERFIELD_SEPARATION_RATIO,
    ExpansionRegime,
    NozzleGeometry,
    area_mach_ratio,
    classify_expansion_regime,
    isentropic_pressure_ratio,
    isentropic_temperature_ratio,
    solve_nozzle_exit,
    solve_supersonic_exit_mach,
    speed_of_sound,
)

GAMMAS = (1.15, 1.20, 1.25, 1.30, 1.40)


def nasa_area_ratio(mach: float, gamma: float) -> float:
    """NASA Glenn ``astar.html`` form, deliberately arranged differently.

    ``A/A* = {[(g+1)/2]^-k} / M * [1 + M^2 (g-1)/2]^k``  with ``k=(g+1)/(2(g-1))``
    """
    k = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (((gamma + 1.0) / 2.0) ** (-k)) / mach * (1.0 + mach * mach * (gamma - 1.0) / 2.0) ** k


# --- 9. area / diameter conversion ----------------------------------------


def test_areas_follow_from_diameters():
    nozzle = NozzleGeometry(throat_diameter_m=0.010, exit_diameter_m=0.020)
    assert nozzle.throat_area_m2 == pytest.approx(math.pi * 0.005**2, rel=1e-14)
    assert nozzle.exit_area_m2 == pytest.approx(math.pi * 0.010**2, rel=1e-14)


def test_area_literals_for_the_reference_nozzle():
    """A_t = pi(0.005)^2 = 7.853981633974483e-5, A_e = 4x that."""
    assert REFERENCE_NOZZLE.throat_area_m2 == pytest.approx(7.853981633974483e-5, rel=1e-14)
    assert REFERENCE_NOZZLE.exit_area_m2 == pytest.approx(3.141592653589793e-4, rel=1e-14)


def test_from_areas_round_trips():
    nozzle = NozzleGeometry.from_areas(7.853981633974483e-5, 3.141592653589793e-4)
    assert nozzle.throat_diameter_m == pytest.approx(0.010, rel=1e-14)
    assert nozzle.exit_diameter_m == pytest.approx(0.020, rel=1e-14)


def test_from_throat_and_expansion_ratio_round_trips():
    for eps in (1.0, 2.5, 4.0, 12.0):
        nozzle = NozzleGeometry.from_throat_and_expansion_ratio(0.012, eps)
        assert nozzle.throat_diameter_m == pytest.approx(0.012, rel=1e-14)
        assert nozzle.expansion_ratio == pytest.approx(eps, rel=1e-13)


# --- 10. expansion ratio identity -----------------------------------------


def test_expansion_ratio_is_the_squared_diameter_ratio():
    nozzle = NozzleGeometry(throat_diameter_m=0.010, exit_diameter_m=0.025)
    assert nozzle.expansion_ratio == pytest.approx((0.025 / 0.010) ** 2, rel=1e-14)


def test_reference_nozzle_expansion_ratio_is_exactly_four():
    assert REFERENCE_NOZZLE.expansion_ratio == pytest.approx(4.0, rel=1e-13)
    assert REFERENCE_NOZZLE.throat_diameter_m == 0.010


# --- 30. invalid geometry rejection ---------------------------------------


@pytest.mark.parametrize("bad", [0.0, -0.01, math.nan, math.inf])
def test_invalid_diameters_are_rejected(bad):
    with pytest.raises(ValueError):
        NozzleGeometry(throat_diameter_m=bad, exit_diameter_m=0.02)
    with pytest.raises(ValueError):
        NozzleGeometry(throat_diameter_m=0.01, exit_diameter_m=bad)


def test_exit_smaller_than_throat_is_rejected():
    with pytest.raises(ValueError):
        NozzleGeometry(throat_diameter_m=0.020, exit_diameter_m=0.010)
    with pytest.raises(ValueError):
        NozzleGeometry.from_areas(3.0e-4, 1.0e-4)


@pytest.mark.parametrize("bad_eps", [0.0, 0.99, -1.0, math.nan])
def test_expansion_ratio_below_one_is_rejected(bad_eps):
    with pytest.raises(ValueError):
        NozzleGeometry.from_throat_and_expansion_ratio(0.01, bad_eps)
    with pytest.raises(ValueError):
        solve_supersonic_exit_mach(bad_eps, 1.2)


def test_geometry_is_immutable():
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        REFERENCE_NOZZLE.throat_diameter_m = 0.02


# --- 11. area-Mach at M = 1 gives epsilon = 1 ------------------------------


@pytest.mark.parametrize("gamma", GAMMAS)
def test_area_ratio_is_exactly_one_at_the_throat(gamma):
    assert area_mach_ratio(1.0, gamma) == pytest.approx(1.0, rel=1e-14)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_area_ratio_matches_the_nasa_arrangement(gamma):
    for mach in (0.2, 0.5, 1.0, 1.5, 2.6194467766658, 4.0, 8.0):
        assert area_mach_ratio(mach, gamma) == pytest.approx(
            nasa_area_ratio(mach, gamma), rel=1e-12
        )


def test_area_ratio_is_monotonic_on_the_supersonic_branch():
    machs = np.linspace(1.0, 8.0, 200)
    ratios = area_mach_ratio(machs, 1.2)
    assert np.all(np.diff(ratios) > 0.0)


def test_area_ratio_is_monotonically_decreasing_on_the_subsonic_branch():
    machs = np.linspace(0.05, 1.0, 200)
    ratios = area_mach_ratio(machs, 1.2)
    assert np.all(np.diff(ratios) < 0.0)


def test_area_ratio_scalar_and_array_agree():
    machs = np.array([1.2, 2.0, 3.5])
    vectorised = area_mach_ratio(machs, 1.2)
    assert isinstance(vectorised, np.ndarray)
    for mach, value in zip(machs.tolist(), vectorised.tolist(), strict=True):
        assert value == pytest.approx(area_mach_ratio(mach, 1.2), rel=1e-15)


# --- 12, 13. supersonic root solve ----------------------------------------


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("eps", [1.5, 2.0, 4.0, 10.0, 50.0, 200.0])
def test_solved_exit_mach_is_supersonic_and_reconstructs_the_area_ratio(gamma, eps):
    mach = solve_supersonic_exit_mach(eps, gamma)
    assert mach > 1.0
    assert nasa_area_ratio(mach, gamma) == pytest.approx(eps, rel=1e-12)


def test_expansion_ratio_of_one_returns_sonic_throat():
    for gamma in GAMMAS:
        assert solve_supersonic_exit_mach(1.0, gamma) == 1.0


def test_reference_exit_mach_hardcoded_literal():
    """epsilon = 4, gamma = 1.2 -> M_e = 2.6194467766658 (area residual < 1e-13)."""
    mach = solve_supersonic_exit_mach(4.0, 1.20)
    assert mach == pytest.approx(2.6194467766658, rel=1e-11)
    assert abs(nasa_area_ratio(mach, 1.20) - 4.0) < 1e-13


def test_exit_mach_increases_with_expansion_ratio():
    machs = [solve_supersonic_exit_mach(eps, 1.2) for eps in (2.0, 4.0, 8.0, 16.0)]
    assert all(later > earlier for earlier, later in itertools.pairwise(machs))


def test_exit_mach_is_stable_across_root_tolerances():
    reference = solve_supersonic_exit_mach(4.0, 1.2, xtol=1e-15)
    for xtol in (1e-6, 1e-9, 1e-12):
        assert solve_supersonic_exit_mach(4.0, 1.2, xtol=xtol) == pytest.approx(
            reference, rel=1e-6
        )


# --- 14, 15. isentropic relations -----------------------------------------


@pytest.mark.parametrize("gamma", GAMMAS)
def test_pressure_ratio_matches_the_hand_written_relation(gamma):
    for mach in (0.0, 0.5, 1.0, 2.6194467766658, 5.0):
        expected = math.pow(1.0 + (gamma - 1.0) * mach * mach / 2.0, -gamma / (gamma - 1.0))
        assert isentropic_pressure_ratio(mach, gamma) == pytest.approx(expected, rel=1e-13)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_temperature_ratio_matches_the_hand_written_relation(gamma):
    for mach in (0.0, 0.5, 1.0, 2.6194467766658, 5.0):
        expected = 1.0 / (1.0 + (gamma - 1.0) * mach * mach / 2.0)
        assert isentropic_temperature_ratio(mach, gamma) == pytest.approx(expected, rel=1e-14)


def test_isentropic_ratios_are_one_at_rest():
    for gamma in GAMMAS:
        assert isentropic_pressure_ratio(0.0, gamma) == pytest.approx(1.0, rel=1e-15)
        assert isentropic_temperature_ratio(0.0, gamma) == pytest.approx(1.0, rel=1e-15)


def test_isentropic_ratios_lie_strictly_between_zero_and_one_when_supersonic():
    for gamma in GAMMAS:
        for mach in (1.5, 2.6194467766658, 6.0):
            p_ratio = isentropic_pressure_ratio(mach, gamma)
            t_ratio = isentropic_temperature_ratio(mach, gamma)
            assert 0.0 < p_ratio < 1.0
            assert 0.0 < t_ratio < 1.0


def test_pressure_and_temperature_ratios_are_linked_by_the_isentropic_exponent():
    """p/p_t = (T/T_t)^(g/(g-1)) is an independent consistency identity."""
    for gamma in GAMMAS:
        for mach in (0.7, 1.0, 3.0):
            p_ratio = isentropic_pressure_ratio(mach, gamma)
            t_ratio = isentropic_temperature_ratio(mach, gamma)
            assert p_ratio == pytest.approx(t_ratio ** (gamma / (gamma - 1.0)), rel=1e-13)


def test_choked_throat_pressure_ratio_matches_the_classical_value():
    """At M = 1, p/p_t = [2/(g+1)]^(g/(g-1))."""
    for gamma in GAMMAS:
        expected = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
        assert isentropic_pressure_ratio(1.0, gamma) == pytest.approx(expected, rel=1e-13)


# --- 16, 17. speed of sound and velocity ----------------------------------


def test_speed_of_sound_matches_sqrt_gamma_r_t():
    assert speed_of_sound(1.4, 287.0, 288.15) == pytest.approx(
        math.sqrt(1.4 * 287.0 * 288.15), rel=1e-14
    )


def test_speed_of_sound_hardcoded_literal_for_air():
    """sqrt(1.4 * 287 * 288.15) = 340.2626485525556 m/s."""
    assert speed_of_sound(1.4, 287.0, 288.15) == pytest.approx(340.2626485525556, rel=1e-12)


def test_speed_of_sound_scales_as_sqrt_temperature():
    base = speed_of_sound(1.2, 377.93, 1500.0)
    assert speed_of_sound(1.2, 377.93, 6000.0) == pytest.approx(2.0 * base, rel=1e-14)


@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf])
def test_speed_of_sound_rejects_invalid_inputs(bad):
    with pytest.raises(ValueError):
        speed_of_sound(1.2, 377.93, bad)
    with pytest.raises(ValueError):
        speed_of_sound(1.2, bad, 1500.0)


@pytest.mark.parametrize("bad_gamma", [1.0, 0.9, 0.0, -1.2, math.nan, math.inf])
def test_invalid_gamma_is_rejected_everywhere(bad_gamma):
    with pytest.raises(ValueError):
        area_mach_ratio(2.0, bad_gamma)
    with pytest.raises(ValueError):
        isentropic_pressure_ratio(2.0, bad_gamma)
    with pytest.raises(ValueError):
        isentropic_temperature_ratio(2.0, bad_gamma)
    with pytest.raises(ValueError):
        solve_supersonic_exit_mach(4.0, bad_gamma)


# --- full exit-state solve -------------------------------------------------

REFERENCE_PC = 2611424.847596018
REFERENCE_TC = 2600.0
REFERENCE_GAMMA = 1.20
REFERENCE_R = 377.930119
SEA_LEVEL_PA = 101325.0


@pytest.fixture(scope="module")
def reference_exit():
    return solve_nozzle_exit(
        REFERENCE_NOZZLE, REFERENCE_PC, REFERENCE_TC, REFERENCE_GAMMA, REFERENCE_R, SEA_LEVEL_PA
    )


def test_reference_exit_state_against_hand_computed_values(reference_exit):
    mach = solve_supersonic_exit_mach(4.0, REFERENCE_GAMMA)
    expected_p = REFERENCE_PC * math.pow(
        1.0 + 0.1 * mach * mach, -REFERENCE_GAMMA / (REFERENCE_GAMMA - 1.0)
    )
    expected_t = REFERENCE_TC / (1.0 + 0.1 * mach * mach)
    expected_a = math.sqrt(REFERENCE_GAMMA * REFERENCE_R * expected_t)
    assert reference_exit.mach == pytest.approx(mach, rel=1e-14)
    assert reference_exit.pressure_pa == pytest.approx(expected_p, rel=1e-13)
    assert reference_exit.temperature_k == pytest.approx(expected_t, rel=1e-13)
    assert reference_exit.speed_of_sound_m_s == pytest.approx(expected_a, rel=1e-13)
    assert reference_exit.velocity_m_s == pytest.approx(mach * expected_a, rel=1e-13)


def test_exit_velocity_is_mach_times_speed_of_sound(reference_exit):
    assert reference_exit.velocity_m_s == pytest.approx(
        reference_exit.mach * reference_exit.speed_of_sound_m_s, rel=1e-15
    )


def test_exit_state_ratios_are_consistent_with_the_absolute_values(reference_exit):
    assert reference_exit.pressure_pa == pytest.approx(
        REFERENCE_PC * reference_exit.pressure_ratio, rel=1e-14
    )
    assert reference_exit.temperature_k == pytest.approx(
        REFERENCE_TC * reference_exit.temperature_ratio, rel=1e-14
    )
    assert 0.0 < reference_exit.pressure_ratio < 1.0
    assert 0.0 < reference_exit.temperature_ratio < 1.0
    assert reference_exit.velocity_m_s > 0.0


def test_area_mach_residual_is_negligible(reference_exit):
    assert abs(reference_exit.area_mach_residual) < 1.0e-12


def test_exit_state_is_independent_of_chamber_pressure_except_through_scaling():
    """M_e, T_e and V_e depend only on epsilon, gamma, R and T_c -- not on p_c."""
    low = solve_nozzle_exit(REFERENCE_NOZZLE, 5.0e5, REFERENCE_TC, REFERENCE_GAMMA,
                            REFERENCE_R, SEA_LEVEL_PA)
    high = solve_nozzle_exit(REFERENCE_NOZZLE, 5.0e6, REFERENCE_TC, REFERENCE_GAMMA,
                             REFERENCE_R, SEA_LEVEL_PA)
    assert low.mach == pytest.approx(high.mach, rel=1e-15)
    assert low.temperature_k == pytest.approx(high.temperature_k, rel=1e-15)
    assert low.velocity_m_s == pytest.approx(high.velocity_m_s, rel=1e-15)
    assert high.pressure_pa == pytest.approx(10.0 * low.pressure_pa, rel=1e-13)


@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf])
def test_exit_solver_rejects_invalid_chamber_conditions(bad):
    with pytest.raises(ValueError):
        solve_nozzle_exit(REFERENCE_NOZZLE, bad, REFERENCE_TC, REFERENCE_GAMMA,
                          REFERENCE_R, SEA_LEVEL_PA)
    with pytest.raises(ValueError):
        solve_nozzle_exit(REFERENCE_NOZZLE, REFERENCE_PC, bad, REFERENCE_GAMMA,
                          REFERENCE_R, SEA_LEVEL_PA)
    with pytest.raises(ValueError):
        solve_nozzle_exit(REFERENCE_NOZZLE, REFERENCE_PC, REFERENCE_TC, REFERENCE_GAMMA,
                          bad, SEA_LEVEL_PA)


@pytest.mark.parametrize("bad_pa", [-1.0, -101325.0, math.nan, math.inf])
def test_exit_solver_rejects_invalid_ambient_pressure(bad_pa):
    with pytest.raises(ValueError):
        solve_nozzle_exit(REFERENCE_NOZZLE, REFERENCE_PC, REFERENCE_TC, REFERENCE_GAMMA,
                          REFERENCE_R, bad_pa)


# --- 23, 24, 25. expansion-regime classification --------------------------


def test_underexpanded_classification():
    assert classify_expansion_regime(2.0e5, 1.0e5) is ExpansionRegime.UNDEREXPANDED


def test_overexpanded_classification():
    assert classify_expansion_regime(5.0e4, 1.0e5) is ExpansionRegime.OVEREXPANDED


def test_near_ideal_classification():
    assert classify_expansion_regime(1.0e5, 1.0e5) is ExpansionRegime.NEAR_IDEAL_EXPANSION
    edge = 1.0e5 * (1.0 + 0.5 * NEAR_IDEAL_PRESSURE_TOLERANCE)
    assert classify_expansion_regime(edge, 1.0e5) is ExpansionRegime.NEAR_IDEAL_EXPANSION


def test_near_ideal_tolerance_boundary_is_respected():
    just_outside = 1.0e5 * (1.0 + 2.0 * NEAR_IDEAL_PRESSURE_TOLERANCE)
    assert classify_expansion_regime(just_outside, 1.0e5) is ExpansionRegime.UNDEREXPANDED


def test_vacuum_is_reported_as_its_own_regime():
    assert classify_expansion_regime(1.0e4, 0.0) is ExpansionRegime.VACUUM


def test_reference_nozzle_is_slightly_underexpanded_at_sea_level(reference_exit):
    assert reference_exit.expansion_regime is ExpansionRegime.UNDEREXPANDED
    assert reference_exit.pressure_pa > SEA_LEVEL_PA


def test_negative_ambient_pressure_is_rejected_by_the_classifier():
    with pytest.raises(ValueError):
        classify_expansion_regime(1.0e5, -1.0)


# --- Summerfield separation flag ------------------------------------------


def test_reference_case_does_not_expect_separation(reference_exit):
    assert not reference_exit.is_separation_expected
    assert reference_exit.pressure_pa > SUMMERFIELD_SEPARATION_RATIO * SEA_LEVEL_PA


def test_a_strongly_overexpanded_case_is_flagged():
    """A large expansion ratio at sea level drives p_e below 0.4 p_a."""
    wide = NozzleGeometry.from_throat_and_expansion_ratio(0.010, 40.0)
    state = solve_nozzle_exit(wide, REFERENCE_PC, REFERENCE_TC, REFERENCE_GAMMA,
                              REFERENCE_R, SEA_LEVEL_PA)
    assert state.expansion_regime is ExpansionRegime.OVEREXPANDED
    assert state.is_separation_expected
    assert state.pressure_pa < SUMMERFIELD_SEPARATION_RATIO * SEA_LEVEL_PA


def test_vacuum_never_flags_separation():
    state = solve_nozzle_exit(REFERENCE_NOZZLE, REFERENCE_PC, REFERENCE_TC, REFERENCE_GAMMA,
                              REFERENCE_R, 0.0)
    assert not state.is_separation_expected
    assert state.expansion_regime is ExpansionRegime.VACUUM
