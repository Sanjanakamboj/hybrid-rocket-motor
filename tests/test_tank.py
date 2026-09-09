"""Independent verification of the saturated-equilibrium N2O tank.

Independence strategy
---------------------
The temperature recovered by this project's explicit bracketed root solve is
checked against **CoolProp's own density-internal-energy flash**, which is a
completely different numerical route to the same thermodynamic state.  Volume,
mass and energy closures are re-formed here from the saturation properties
rather than read back from the tank object.
"""

from __future__ import annotations

import itertools
import math

import CoolProp.CoolProp as CP
import pytest

from hybrid_rocket_motor.nitrous_properties import (
    MIN_VERIFIED_TEMPERATURE_K,
    PropertyRangeStatus,
    saturated_state,
)
from hybrid_rocket_motor.tank import (
    REFERENCE_TANK,
    REFERENCE_TANK_FILL_FRACTION,
    REFERENCE_TANK_TEMPERATURE_K,
    NitrousTank,
    TankPhase,
)

TANK = REFERENCE_TANK
T0 = REFERENCE_TANK_TEMPERATURE_K
FILL0 = REFERENCE_TANK_FILL_FRACTION


@pytest.fixture(scope="module")
def initial():
    return TANK.initial_state(T0, FILL0)


# --- 10. initial state closure --------------------------------------------


def test_initial_state_matches_an_independent_construction(initial):
    """Rebuild the initial state longhand from the saturation properties."""
    sat = saturated_state(T0)
    liquid_volume = FILL0 * TANK.volume_m3
    vapour_volume = (1.0 - FILL0) * TANK.volume_m3
    liquid_mass = liquid_volume * sat.liquid_density_kg_m3
    vapour_mass = vapour_volume * sat.vapour_density_kg_m3
    assert initial.liquid_mass_kg == pytest.approx(liquid_mass, rel=1e-14)
    assert initial.vapour_mass_kg == pytest.approx(vapour_mass, rel=1e-14)
    assert initial.total_mass_kg == pytest.approx(liquid_mass + vapour_mass, rel=1e-14)
    assert initial.total_internal_energy_j == pytest.approx(
        liquid_mass * sat.liquid_internal_energy_j_kg
        + vapour_mass * sat.vapour_internal_energy_j_kg,
        rel=1e-14,
    )
    assert initial.pressure_pa == pytest.approx(sat.pressure_pa, rel=1e-14)


def test_reference_tank_initial_values_hardcoded():
    """10 L, 80 % fill, 293.15 K -> 50.525 bar, 6.597 kg total, 6.281 kg liquid."""
    state = TANK.initial_state(T0, FILL0)
    assert TANK.volume_m3 == 0.010
    assert state.pressure_pa / 1e5 == pytest.approx(50.525, rel=1e-4)
    assert state.total_mass_kg == pytest.approx(6.5968, rel=1e-4)
    assert state.liquid_mass_kg == pytest.approx(6.2808, rel=1e-4)
    assert state.vapour_mass_kg == pytest.approx(0.3160, rel=1e-3)


# --- 8, 9. mass and volume closure ----------------------------------------


def test_phase_masses_sum_to_the_total(initial):
    assert initial.liquid_mass_kg + initial.vapour_mass_kg == pytest.approx(
        initial.total_mass_kg, rel=1e-14
    )


def test_phase_volumes_sum_to_the_tank_volume(initial):
    assert initial.liquid_volume_m3 + initial.vapour_volume_m3 == pytest.approx(
        TANK.volume_m3, rel=1e-14
    )


def test_closures_hold_after_discharge(initial):
    for removed in (0.5, 2.0, 4.0):
        state = TANK.state_from_mass_and_energy(
            initial.total_mass_kg - removed,
            initial.total_internal_energy_j - removed * initial.liquid_enthalpy_j_kg,
        )
        assert state.liquid_mass_kg + state.vapour_mass_kg == pytest.approx(
            state.total_mass_kg, rel=1e-12
        )
        assert state.liquid_volume_m3 + state.vapour_volume_m3 == pytest.approx(
            TANK.volume_m3, rel=1e-12
        )
        # Volume closure re-formed independently from the saturation densities.
        sat = saturated_state(state.temperature_k)
        rebuilt = (
            state.liquid_mass_kg / sat.liquid_density_kg_m3
            + state.vapour_mass_kg / sat.vapour_density_kg_m3
        )
        assert rebuilt == pytest.approx(TANK.volume_m3, rel=1e-10)


# --- 12. the flash against an independent route ---------------------------


def test_temperature_solve_matches_coolprops_own_density_energy_flash(initial):
    """Our bracketed root solve vs CoolProp's internal D-U flash."""
    for removed in (0.0, 1.0, 3.0, 5.0):
        mass = initial.total_mass_kg - removed
        energy = initial.total_internal_energy_j - removed * initial.liquid_enthalpy_j_kg
        state = TANK.state_from_mass_and_energy(mass, energy)
        reference = CP.PropsSI(
            "T", "D", mass / TANK.volume_m3, "U", energy / mass, "NitrousOxide"
        )
        assert state.temperature_k == pytest.approx(reference, abs=1e-6)


def test_flash_round_trips_the_initial_state(initial):
    state = TANK.state_from_mass_and_energy(
        initial.total_mass_kg, initial.total_internal_energy_j
    )
    assert state.temperature_k == pytest.approx(T0, abs=1e-7)
    assert state.pressure_pa == pytest.approx(initial.pressure_pa, rel=1e-9)
    assert state.liquid_mass_kg == pytest.approx(initial.liquid_mass_kg, rel=1e-7)


def test_energy_closure_of_the_solved_state(initial):
    """U must equal m_l u_l + m_v u_v at the solved temperature."""
    mass = initial.total_mass_kg - 2.5
    energy = initial.total_internal_energy_j - 2.5 * initial.liquid_enthalpy_j_kg
    state = TANK.state_from_mass_and_energy(mass, energy)
    sat = saturated_state(state.temperature_k)
    rebuilt = (
        state.liquid_mass_kg * sat.liquid_internal_energy_j_kg
        + state.vapour_mass_kg * sat.vapour_internal_energy_j_kg
    )
    assert rebuilt == pytest.approx(energy, rel=1e-9)


# --- 13. saturation pressure ----------------------------------------------


def test_pressure_is_the_saturation_pressure_of_the_solved_temperature(initial):
    for removed in (0.0, 1.5, 3.5):
        state = TANK.state_from_mass_and_energy(
            initial.total_mass_kg - removed,
            initial.total_internal_energy_j - removed * initial.liquid_enthalpy_j_kg,
        )
        assert state.pressure_pa == pytest.approx(
            saturated_state(state.temperature_k).pressure_pa, rel=1e-14
        )


# --- 11, 14, 15. discharge behaviour --------------------------------------


def test_liquid_mass_and_pressure_fall_monotonically_during_discharge(initial):
    removed_masses = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
    states = [
        TANK.state_from_mass_and_energy(
            initial.total_mass_kg - m,
            initial.total_internal_energy_j - m * initial.liquid_enthalpy_j_kg,
        )
        for m in removed_masses
    ]
    liquids = [s.liquid_mass_kg for s in states]
    pressures = [s.pressure_pa for s in states]
    temperatures = [s.temperature_k for s in states]
    assert all(later < earlier for earlier, later in itertools.pairwise(liquids))
    assert all(later < earlier for earlier, later in itertools.pairwise(pressures))
    assert all(later < earlier for earlier, later in itertools.pairwise(temperatures))


def test_self_pressurising_blowdown_cools_the_tank(initial):
    """Withdrawing liquid must cool the remaining propellant, not warm it."""
    state = TANK.state_from_mass_and_energy(
        initial.total_mass_kg - 2.0,
        initial.total_internal_energy_j - 2.0 * initial.liquid_enthalpy_j_kg,
    )
    assert state.temperature_k < T0
    assert state.temperature_k == pytest.approx(286.21, abs=0.05)
    assert state.pressure_pa < initial.pressure_pa


def test_vapour_quality_rises_as_liquid_is_withdrawn(initial):
    qualities = [
        TANK.state_from_mass_and_energy(
            initial.total_mass_kg - m,
            initial.total_internal_energy_j - m * initial.liquid_enthalpy_j_kg,
        ).vapour_quality
        for m in (0.0, 1.0, 2.0, 3.0, 4.0)
    ]
    assert all(later > earlier for earlier, later in itertools.pairwise(qualities))


def test_phase_masses_are_never_negative(initial):
    for removed in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5):
        state = TANK.state_from_mass_and_energy(
            initial.total_mass_kg - removed,
            initial.total_internal_energy_j - removed * initial.liquid_enthalpy_j_kg,
        )
        assert state.liquid_mass_kg >= 0.0
        assert state.vapour_mass_kg >= 0.0
        assert 0.0 <= state.vapour_quality <= 1.0
        assert 0.0 <= state.liquid_fill_fraction <= 1.0


# --- 16. phase / status transitions ---------------------------------------


def test_two_phase_state_is_reported_while_liquid_remains(initial):
    state = TANK.state_from_mass_and_energy(
        initial.total_mass_kg - 3.0,
        initial.total_internal_energy_j - 3.0 * initial.liquid_enthalpy_j_kg,
    )
    assert state.phase is TankPhase.TWO_PHASE
    assert state.has_liquid


def test_vapour_only_states_are_rejected_rather_than_silently_extrapolated():
    """A charge far past liquid depletion must be refused, not given a fake T and p."""
    tank = NitrousTank(volume_m3=0.010)
    vapour_state = saturated_state(250.0)
    mass = 0.010 * vapour_state.vapour_density_kg_m3 * 0.5
    energy = mass * vapour_state.vapour_internal_energy_j_kg * 1.2
    with pytest.raises(ValueError, match="single-phase vapour"):
        tank.state_from_mass_and_energy(mass, energy)


def test_compressed_liquid_states_are_rejected():
    """A deeply over-filled tank has no two-phase root either.

    The quality tolerance is far less sensitive on the liquid side than on the
    vapour side, because ``v_v - v_l`` is three orders of magnitude larger than
    the liquid specific volume itself: a modest mass overfill barely moves the
    quality.  A genuinely compressed charge is needed to trip the guard, and in
    practice :meth:`NitrousTank.initial_state` prevents over-filling by
    construction.
    """
    tank = NitrousTank(volume_m3=0.010)
    liquid_state = saturated_state(250.0)
    mass = 0.010 * liquid_state.liquid_density_kg_m3 * 3.0
    energy = mass * liquid_state.liquid_internal_energy_j_kg
    with pytest.raises(ValueError, match="compressed liquid"):
        tank.state_from_mass_and_energy(mass, energy)


def test_the_unclamped_quality_crosses_one_smoothly_at_depletion(initial):
    """Event detection needs a signed quantity; physical masses stay clamped."""
    fractions = []
    for removed in (5.9, 6.0, 6.1):
        state = TANK.state_from_mass_and_energy(
            initial.total_mass_kg - removed,
            initial.total_internal_energy_j - removed * initial.liquid_enthalpy_j_kg,
        )
        fractions.append(state.liquid_fraction_remaining)
        assert state.liquid_mass_kg >= 0.0
    assert all(later < earlier for earlier, later in itertools.pairwise(fractions))


def test_reference_tank_state_is_inside_the_checked_property_band(initial):
    assert initial.range_status is PropertyRangeStatus.VERIFIED
    assert initial.is_within_verified_property_range
    assert initial.temperature_k > MIN_VERIFIED_TEMPERATURE_K


# --- input validation ------------------------------------------------------


@pytest.mark.parametrize("bad", [0.0, -0.01, math.nan, math.inf])
def test_invalid_tank_volume_is_rejected(bad):
    with pytest.raises(ValueError):
        NitrousTank(volume_m3=bad)


@pytest.mark.parametrize("bad_fill", [0.0, 1.0, -0.1, 1.5, math.nan, math.inf])
def test_invalid_fill_fraction_is_rejected(bad_fill):
    with pytest.raises(ValueError):
        TANK.initial_state(T0, bad_fill)


@pytest.mark.parametrize("bad_temperature", [180.0, 310.0, math.nan, math.inf])
def test_invalid_initial_temperature_is_rejected(bad_temperature):
    with pytest.raises(ValueError):
        TANK.initial_state(bad_temperature, 0.8)


@pytest.mark.parametrize("bad_mass", [0.0, -1.0, math.nan, math.inf])
def test_invalid_mass_is_rejected(bad_mass, initial):
    with pytest.raises(ValueError):
        TANK.state_from_mass_and_energy(bad_mass, initial.total_internal_energy_j)


@pytest.mark.parametrize("bad_energy", [math.nan, math.inf, -math.inf])
def test_non_finite_energy_is_rejected(bad_energy, initial):
    with pytest.raises(ValueError):
        TANK.state_from_mass_and_energy(initial.total_mass_kg, bad_energy)


def test_tank_and_state_are_immutable(initial):
    with pytest.raises(Exception):  # noqa: B017 - dataclass raises FrozenInstanceError
        TANK.volume_m3 = 0.02
    with pytest.raises(Exception):  # noqa: B017
        initial.temperature_k = 300.0


def test_a_larger_tank_at_the_same_fill_holds_proportionally_more():
    small = NitrousTank(volume_m3=0.005).initial_state(T0, 0.8)
    large = NitrousTank(volume_m3=0.010).initial_state(T0, 0.8)
    assert large.total_mass_kg == pytest.approx(2.0 * small.total_mass_kg, rel=1e-12)
    assert large.pressure_pa == pytest.approx(small.pressure_pa, rel=1e-14)
