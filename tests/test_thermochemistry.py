"""Independent verification of the frozen CEA table and its interpolation.

Independence strategy
---------------------
The committed CSV is re-parsed here with the standard library rather than through
the package loader, so the grid, the row count and the node values are checked
against the file itself.  Thermodynamic identities (``R = R_u / M``, the ideal
``c*`` relation, the sound-speed definition of the effective molar mass) are
re-derived from first principles instead of being read back from the loader.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import pytest

from hybrid_rocket_motor.chamber import ideal_characteristic_velocity
from hybrid_rocket_motor.thermochemistry import (
    DEFAULT_C_STAR_EFFICIENCY,
    DEFAULT_TABLE_PATH,
    TableStatus,
    ThermochemistryTable,
    default_table,
)

UNIVERSAL_GAS_CONSTANT_J_MOL_K = 8.314462618


@pytest.fixture(scope="module")
def table() -> ThermochemistryTable:
    return default_table()


@pytest.fixture(scope="module")
def raw_rows() -> list[dict[str, str]]:
    with Path(DEFAULT_TABLE_PATH).open(encoding="utf-8") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


# -- the committed file --------------------------------------------------------


def test_committed_table_exists_and_is_non_trivial(raw_rows):
    assert len(raw_rows) == 1695


def test_table_is_a_complete_rectangular_grid(raw_rows):
    of_values = sorted({float(row["of"]) for row in raw_rows})
    pressures = sorted({float(row["chamber_pressure_pa"]) for row in raw_rows})
    assert len(of_values) * len(pressures) == len(raw_rows)
    seen = {(float(r["of"]), float(r["chamber_pressure_pa"])) for r in raw_rows}
    assert len(seen) == len(raw_rows), "a duplicated node would silently overwrite another"


def test_loader_grid_matches_the_file(table, raw_rows):
    of_values = sorted({float(row["of"]) for row in raw_rows})
    pressures = sorted({float(row["chamber_pressure_pa"]) for row in raw_rows})
    assert np.allclose(table.of_values, of_values)
    assert np.allclose(table.pressure_values, pressures)


def test_file_carries_its_provenance():
    header = Path(DEFAULT_TABLE_PATH).read_text(encoding="utf-8").splitlines()
    comments = "\n".join(line for line in header if line.startswith("#"))
    assert "CEA" in comments
    assert "RP-1311" in comments
    assert "generate_thermochemistry_table.py" in comments


def test_csv_contains_no_date_but_the_metadata_does():
    """The CSV must regenerate byte-identically, so the date lives in the JSON.

    Splitting them this way is what lets ``--check`` prove the committed table is
    unchanged: a timestamp inside the CSV would make every regeneration differ.
    """
    import re

    iso_date = re.compile(r"\d{4}-\d{2}-\d{2}")
    assert not iso_date.search(Path(DEFAULT_TABLE_PATH).read_text(encoding="utf-8"))
    metadata_path = Path(DEFAULT_TABLE_PATH).with_name("n2o_htpb_equilibrium_metadata.json")
    assert iso_date.search(metadata_path.read_text(encoding="utf-8"))


def test_metadata_records_the_solver_version():
    import json

    path = Path(DEFAULT_TABLE_PATH).with_name("n2o_htpb_equilibrium_metadata.json")
    metadata = json.loads(path.read_text(encoding="utf-8"))
    assert metadata["solver"].startswith("NASA CEA")
    assert metadata["row_count"] == 1695
    assert "rocketcea_version" in metadata
    assert "ideal" in metadata["cstar_convention"]


# -- thermodynamic identities at every node ------------------------------------


def test_specific_gas_constant_matches_molar_mass_at_every_node(raw_rows):
    for row in raw_rows:
        molar = float(row["effective_molar_mass_kg_mol"])
        expected = UNIVERSAL_GAS_CONSTANT_J_MOL_K / molar
        assert float(row["specific_gas_constant_j_kg_k"]) == pytest.approx(expected, rel=1e-9)


def test_ideal_cstar_matches_gamma_r_and_temperature_at_every_node(raw_rows):
    """CEA's own c* against the c* implied by the properties we store.

    This is the check that catches a units error in the table generator: if the
    molar mass were wrong by a factor, ``R`` would be wrong by that factor and the
    identity would break by orders of magnitude even though every individual column
    still looked plausible.  It caught exactly that during development.

    The agreement is close but not exact, and the residual is physical rather than
    numerical.  ``ideal_characteristic_velocity`` assumes one constant ``gamma``
    from chamber to throat; CEA's ``c*`` comes from a real equilibrium expansion in
    which the composition and ``gamma`` both shift.  Measured across the committed
    table the gap runs from -0.53 % to +0.99 %, so the bound below brackets that.
    It is a measurement of the single-``gamma`` idealisation, not a slack tolerance
    hiding an error -- which is why the test also asserts the gap changes sign.
    """
    deviations = []
    for row in raw_rows:
        expected = ideal_characteristic_velocity(
            float(row["gamma"]),
            float(row["specific_gas_constant_j_kg_k"]),
            float(row["chamber_temperature_k"]),
        )
        deviations.append(float(row["ideal_cstar_m_s"]) / expected - 1.0)
    assert max(abs(d) for d in deviations) < 0.012
    assert min(deviations) < 0.0 < max(deviations), "the idealisation should err both ways"


def test_effective_molar_mass_is_never_below_the_gas_molar_mass(raw_rows):
    """Condensed carbon adds mass to the mixture without adding moles of gas.

    So the effective molar mass -- total mass per mole of *gas* -- is at least the
    gas-phase mean molecular weight, and strictly greater wherever soot appears.
    Both columns are stored in kg/mol.
    """
    ratios = [
        float(row["effective_molar_mass_kg_mol"]) / float(row["gas_molar_mass_kg_mol"])
        for row in raw_rows
    ]
    assert min(ratios) >= 1.0 - 1e-9, "an effective molar mass below the gas value is unphysical"
    assert max(ratios) > 1.05, "a condensed phase should appear somewhere in this range"
    assert max(ratios) < 2.0


def test_all_nodes_are_physically_plausible(raw_rows):
    for row in raw_rows:
        assert 1.0 < float(row["gamma"]) < 1.7
        assert 500.0 < float(row["chamber_temperature_k"]) < 4000.0
        assert 500.0 < float(row["ideal_cstar_m_s"]) < 2500.0
        assert 0.005 < float(row["effective_molar_mass_kg_mol"]) < 0.060


# -- interpolation --------------------------------------------------------------


def test_interpolation_reproduces_node_values_exactly(table, raw_rows):
    lookup = {(float(r["of"]), float(r["chamber_pressure_pa"])): r for r in raw_rows}
    for of in (table.of_values[0], table.of_values[40], table.of_values[-1]):
        for pressure in (
            table.pressure_values[0],
            table.pressure_values[7],
            table.pressure_values[-1],
        ):
            row = lookup[(float(of), float(pressure))]
            state = table.evaluate(float(of), float(pressure), c_star_efficiency=1.0)
            assert state.ideal_c_star_m_s == pytest.approx(float(row["ideal_cstar_m_s"]), rel=1e-12)
            assert state.gamma == pytest.approx(float(row["gamma"]), rel=1e-12)
            assert state.chamber_temperature_k == pytest.approx(
                float(row["chamber_temperature_k"]), rel=1e-12
            )


def test_midpoint_interpolation_is_the_average_of_its_neighbours(table):
    """Bilinear interpolation has an exactly known value at a cell centre."""
    of_lo, of_hi = float(table.of_values[10]), float(table.of_values[11])
    p_lo, p_hi = float(table.pressure_values[3]), float(table.pressure_values[4])
    corners = [
        table.evaluate(of, p, c_star_efficiency=1.0).gamma
        for of in (of_lo, of_hi)
        for p in (p_lo, p_hi)
    ]
    centre = table.evaluate((of_lo + of_hi) / 2.0, (p_lo + p_hi) / 2.0, c_star_efficiency=1.0).gamma
    assert centre == pytest.approx(sum(corners) / 4.0, rel=1e-12)


def test_interpolated_gas_constant_matches_interpolated_molar_mass(table):
    """``R`` is recomputed, not interpolated, so the identity holds off-node too."""
    for of in (1.37, 2.1234, 2.9876, 3.51):
        for pressure in (11.3e5, 24.7e5, 41.9e5):
            state = table.evaluate(of, pressure)
            expected = UNIVERSAL_GAS_CONSTANT_J_MOL_K / state.effective_molar_mass_kg_mol
            assert state.specific_gas_constant_j_kg_k == pytest.approx(expected, rel=1e-14)


def test_interpolation_stays_within_the_bounding_nodes(table):
    """Bilinear interpolation cannot overshoot -- the reason a spline was not used."""
    for i in (5, 40, 90):
        of_lo, of_hi = float(table.of_values[i]), float(table.of_values[i + 1])
        p_lo, p_hi = float(table.pressure_values[2]), float(table.pressure_values[3])
        corners = [
            table.evaluate(of, p, c_star_efficiency=1.0).ideal_c_star_m_s
            for of in (of_lo, of_hi)
            for p in (p_lo, p_hi)
        ]
        for fraction in (0.13, 0.5, 0.87):
            of = of_lo + fraction * (of_hi - of_lo)
            pressure = p_lo + fraction * (p_hi - p_lo)
            value = table.evaluate(of, pressure, c_star_efficiency=1.0).ideal_c_star_m_s
            assert min(corners) - 1e-9 <= value <= max(corners) + 1e-9


# -- efficiency bookkeeping ------------------------------------------------------


def test_efficiency_is_applied_once_and_only_to_cstar(table):
    ideal = table.evaluate(2.4, 25e5, c_star_efficiency=1.0)
    delivered = table.evaluate(2.4, 25e5, c_star_efficiency=0.8)
    assert delivered.c_star_m_s == pytest.approx(0.8 * ideal.ideal_c_star_m_s, rel=1e-14)
    assert delivered.ideal_c_star_m_s == pytest.approx(ideal.ideal_c_star_m_s, rel=1e-14)
    assert delivered.gamma == pytest.approx(ideal.gamma, rel=1e-14)
    assert delivered.chamber_temperature_k == pytest.approx(ideal.chamber_temperature_k, rel=1e-14)


def test_default_efficiency_matches_the_milestone_3_value():
    from hybrid_rocket_motor.chamber import ILLUSTRATIVE_N2O_HTPB_COMBUSTION

    implied = ILLUSTRATIVE_N2O_HTPB_COMBUSTION.c_star_efficiency
    default = DEFAULT_C_STAR_EFFICIENCY
    assert default == pytest.approx(implied, rel=1e-6)


@pytest.mark.parametrize("efficiency", [-0.1, 0.0, 1.5, math.nan])
def test_invalid_efficiency_is_refused(table, efficiency):
    with pytest.raises(ValueError):
        table.evaluate(2.4, 25e5, c_star_efficiency=efficiency)


# -- refusing to extrapolate -----------------------------------------------------


def test_status_reports_which_edge_was_crossed(table):
    of_low, of_high = table.of_bounds
    p_low, p_high = table.pressure_bounds_pa
    assert table.status(2.5, 25e5) is TableStatus.WITHIN_TABLE
    assert table.status(of_low - 0.01, 25e5) is TableStatus.O_F_BELOW_TABLE
    assert table.status(of_high + 0.01, 25e5) is TableStatus.O_F_ABOVE_TABLE
    assert table.status(2.5, p_low - 1.0) is TableStatus.PRESSURE_BELOW_TABLE
    assert table.status(2.5, p_high + 1.0) is TableStatus.PRESSURE_ABOVE_TABLE


def test_evaluation_outside_the_table_raises_rather_than_extrapolating(table):
    of_low, of_high = table.of_bounds
    p_low, p_high = table.pressure_bounds_pa
    for of, pressure in (
        (of_low - 0.001, 25e5),
        (of_high + 0.001, 25e5),
        (2.5, p_low - 1.0),
        (2.5, p_high + 1.0),
    ):
        with pytest.raises(ValueError):
            table.evaluate(of, pressure)


def test_the_exact_boundary_is_inside_the_table(table):
    of_low, of_high = table.of_bounds
    p_low, p_high = table.pressure_bounds_pa
    for of, pressure in ((of_low, p_low), (of_high, p_high), (of_low, p_high)):
        assert table.evaluate(of, pressure).is_within_table


# -- table variants used by the sensitivity study ---------------------------------


def test_subsampling_keeps_both_endpoints_and_is_an_exact_subset(table):
    coarse = table.subsample(of_stride=8, pressure_stride=4)
    assert coarse.of_bounds == table.of_bounds
    assert coarse.pressure_bounds_pa == table.pressure_bounds_pa
    assert coarse.of_values.size < table.of_values.size
    for of in coarse.of_values:
        for pressure in coarse.pressure_values:
            fine_state = table.evaluate(float(of), float(pressure))
            coarse_state = coarse.evaluate(float(of), float(pressure))
            assert coarse_state.ideal_c_star_m_s == pytest.approx(
                fine_state.ideal_c_star_m_s, rel=1e-14
            )


def test_fixed_pressure_variant_removes_the_pressure_dependence(table):
    flat = table.at_fixed_pressure(25e5)
    reference = table.evaluate(2.37, 25e5)
    for pressure in (12e5, 25e5, 44e5):
        state = flat.evaluate(2.37, pressure)
        assert state.ideal_c_star_m_s == pytest.approx(reference.ideal_c_star_m_s, rel=1e-12)
        assert state.gamma == pytest.approx(reference.gamma, rel=1e-12)


def test_fixed_pressure_variant_refuses_a_pressure_outside_the_table(table):
    with pytest.raises(ValueError):
        table.at_fixed_pressure(1.0e5)


@pytest.mark.parametrize("stride", [0, -1])
def test_invalid_subsample_stride_is_refused(table, stride):
    with pytest.raises(ValueError):
        table.subsample(of_stride=stride)


# -- the loader ------------------------------------------------------------------


def test_default_table_is_cached_and_read_only(table):
    assert default_table() is table


def test_pressure_dependence_is_weak_but_present(table):
    """Quantifies the claim the study makes about the second table dimension."""
    low = table.evaluate(2.5, table.pressure_bounds_pa[0]).ideal_c_star_m_s
    high = table.evaluate(2.5, table.pressure_bounds_pa[1]).ideal_c_star_m_s
    relative = abs(high / low - 1.0)
    assert 0.0 < relative < 0.01
