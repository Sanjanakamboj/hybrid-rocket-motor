"""Generate the frozen N2O/HTPB equilibrium thermochemistry table with NASA CEA.

This script is a **generation-time tool only**.  It requires ``rocketcea``, which
wraps the actual NASA CEA FORTRAN code (Gordon & McBride, NASA RP-1311) and needs
a Fortran compiler to build.  The runtime package never imports it: everything
downstream reads the committed CSV, so a normal install has no CEA dependency.

Install the optional extra to regenerate::

    pip install -e ".[cea]"
    python scripts/generate_thermochemistry_table.py

Outputs
-------
``data/n2o_htpb_equilibrium.csv``
    The table itself, with a ``#``-prefixed provenance header.  Deliberately
    contains **no timestamp**, so regeneration is byte-identical.
``data/n2o_htpb_equilibrium_metadata.json``
    Generation date and tool versions.  This file *does* change on regeneration.

Reactant representation
-----------------------
Taken verbatim from RocketCEA's propellant library, which uses NASA CEA's own
reactant cards::

    fuel R-45(HTPB FROM_RPL_DATA) C 7.3165 H 10.3360 O 0.1063  wt%=100.00
    h,cal= 1200.0  t(k)=298.15  rho=0.9220

    oxid NitrousOxide  N 2.0 O 1.0  wt%=100.00
    h,cal= 19467.0  t(k)=298.15

So HTPB is the R-45 hydroxyl-terminated polybutadiene surrogate with empirical
formula C(7.3165) H(10.3360) O(0.1063) and an assigned enthalpy of +1200 cal/mol
at 298.15 K.  **This is a surrogate composition, not a unique molecule**; real
HTPB varies with formulation, cure system and any additives, none of which are
represented.

The oxidizer is N2O at its 298.15 K standard state.  The Milestone 4 tank
actually delivers *saturated liquid* N2O at 281-293 K, whose enthalpy is lower by
roughly the latent heat.  That difference is **not** accounted for in the
baseline table; the ``N2O_nbp`` card (liquid at the 184.4 K normal boiling point)
bounds the effect and is used as a sensitivity case in the study script.

Effective molecular weight
--------------------------
At the fuel-rich mixture ratios this motor runs at, CEA predicts **condensed
carbon** in the chamber.  CEA then reports two different molecular weights:

* ``MW, MOL WT`` -- the gas-phase molecular weight;
* ``M, (1/n)``   -- the mass of the *total* mixture per mole of gas.

Only the second reproduces CEA's own sound speed, and it is the one a
single-species nozzle model needs, because that model carries the whole mass flow
as if it were gas.  RocketCEA's ``get_Chamber_MolWt_gamma`` returns the first.
This script therefore derives the effective value from CEA's own consistent
outputs::

    M_eff = gamma R_u T_c / a_chamber^2

which reproduces the printed ``M, (1/n)`` to four decimal places.  Both are
stored; the ratio flags where a condensed phase is present.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import sys
from pathlib import Path

UNIVERSAL_GAS_CONSTANT_J_MOL_K = 8.314462618

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
TABLE_PATH = DATA_DIR / "n2o_htpb_equilibrium.csv"
METADATA_PATH = DATA_DIR / "n2o_htpb_equilibrium_metadata.json"
NOZZLE_CHEMISTRY_PATH = DATA_DIR / "n2o_htpb_nozzle_chemistry.csv"

OX_NAME = "N2O"
FUEL_NAME = "HTPB"

#: Mixture-ratio grid.  The Milestone 4 blowdown envelope is O/F 1.558-2.909
#: across all its cases; this spans it with generous margin on both sides,
#: because the Milestone 5 trajectories are not the Milestone 4 ones.
OF_MIN, OF_MAX, OF_STEP = 1.20, 4.00, 0.025

#: Chamber-pressure grid [Pa].  The Milestone 4 envelope is 20.3-31.6 bar.
PC_MIN_BAR, PC_MAX_BAR, PC_STEP_BAR = 10.0, 45.0, 2.5

#: Passed to RocketCEA calls that require it.  Chamber properties do not depend
#: on the expansion ratio; the script asserts that.
REFERENCE_EPS = 4.0

#: Chamber pressure at which the equilibrium/frozen expansion comparison is
#: taken [Pa].  Near the middle of the Milestone 5 blowdown envelope.
NOZZLE_CHEMISTRY_PC_PA = 25.0e5

#: Coarser O/F grid for that comparison -- it is a sensitivity check, not a
#: runtime input, so it does not need the resolution of the main table.
NOZZLE_CHEMISTRY_OF_STEP = 0.10

NOZZLE_CHEMISTRY_COLUMNS = (
    "of",
    "chamber_pressure_pa",
    "expansion_ratio",
    "vacuum_isp_equilibrium_s",
    "vacuum_isp_frozen_s",
    "frozen_converged",
)

COLUMNS = (
    "of",
    "chamber_pressure_pa",
    "ideal_cstar_m_s",
    "gamma",
    "chamber_temperature_k",
    "effective_molar_mass_kg_mol",
    "gas_molar_mass_kg_mol",
    "specific_gas_constant_j_kg_k",
)


def _grid(minimum: float, maximum: float, step: float) -> list[float]:
    count = round((maximum - minimum) / step) + 1
    return [round(minimum + i * step, 10) for i in range(count)]


def _evaluate(cea, of: float, pressure_pa: float) -> dict[str, float]:
    cstar = float(cea.get_Cstar(Pc=pressure_pa, MR=of))
    temperature = float(cea.get_Tcomb(Pc=pressure_pa, MR=of))
    gas_molar_mass, gamma = cea.get_Chamber_MolWt_gamma(Pc=pressure_pa, MR=of, eps=REFERENCE_EPS)
    sonic = float(cea.get_SonicVelocities(Pc=pressure_pa, MR=of, eps=REFERENCE_EPS)[0])
    # a^2 = gamma R T with R = R_u / M, so M = gamma R_u T / a^2 comes out
    # directly in kg/mol when R_u is in J/(mol K). RocketCEA's gas molar mass, by
    # contrast, is reported in g/mol and does need converting.
    effective_molar_mass = float(gamma) * UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature
    effective_molar_mass /= sonic * sonic
    return {
        "of": of,
        "chamber_pressure_pa": pressure_pa,
        "ideal_cstar_m_s": cstar,
        "gamma": float(gamma),
        "chamber_temperature_k": temperature,
        "effective_molar_mass_kg_mol": effective_molar_mass,
        "gas_molar_mass_kg_mol": float(gas_molar_mass) * 1.0e-3,
        "specific_gas_constant_j_kg_k": UNIVERSAL_GAS_CONSTANT_J_MOL_K / effective_molar_mass,
    }


def _sanity_check(cea) -> None:
    """Assert the assumptions this generator relies on, before writing anything."""
    wide = cea.get_Chamber_MolWt_gamma(Pc=26.11e5, MR=2.514, eps=4.0)
    narrow = cea.get_Chamber_MolWt_gamma(Pc=26.11e5, MR=2.514, eps=12.0)
    if abs(wide[1] - narrow[1]) > 1e-9 or abs(wide[0] - narrow[0]) > 1e-9:
        raise RuntimeError("chamber properties unexpectedly depend on the expansion ratio")

    row = _evaluate(cea, 2.514, 26.11e5)
    sonic = math.sqrt(
        row["gamma"] * row["specific_gas_constant_j_kg_k"] * row["chamber_temperature_k"]
    )
    reference = float(cea.get_SonicVelocities(Pc=26.11e5, MR=2.514, eps=REFERENCE_EPS)[0])
    if abs(sonic - reference) > 1e-6 * reference:
        raise RuntimeError(
            f"effective molar mass does not reproduce the CEA sound speed: {sonic} vs {reference}"
        )


def generate(verbose: bool = True) -> list[dict[str, float]]:
    try:
        from rocketcea.cea_obj_w_units import CEA_Obj
    except ImportError as error:  # pragma: no cover - generation-time only
        raise SystemExit(
            "rocketcea is required to regenerate the table. Install the optional extra:\n"
            '    pip install -e ".[cea]"\n'
            "It wraps the NASA CEA FORTRAN source and needs a Fortran compiler."
        ) from error

    cea = CEA_Obj(
        oxName=OX_NAME,
        fuelName=FUEL_NAME,
        cstar_units="m/s",
        pressure_units="Pa",
        temperature_units="K",
        sonic_velocity_units="m/s",
    )
    _sanity_check(cea)

    of_values = _grid(OF_MIN, OF_MAX, OF_STEP)
    pressure_values = [round(p * 1e5, 6) for p in _grid(PC_MIN_BAR, PC_MAX_BAR, PC_STEP_BAR)]

    rows = []
    for of in of_values:
        for pressure in pressure_values:
            rows.append(_evaluate(cea, of, pressure))
    if verbose:
        print(
            f"  evaluated {len(of_values)} O/F x {len(pressure_values)} pressures "
            f"= {len(rows)} CEA points"
        )
    return rows


def generate_nozzle_chemistry(verbose: bool = True) -> list[dict[str, float]]:
    """Equilibrium vs frozen expansion at the reference pressure and area ratio.

    Our nozzle expands with a single ``gamma`` taken from the chamber, which sits
    between the two limits CEA can compute: *shifting* (equilibrium) composition,
    which keeps recombining and releases more energy, and *frozen* composition,
    which does not.  The gap between them bounds what the single-gamma assumption
    can be worth.

    CEA's frozen option does **not** converge across part of this range.  Where it
    does not, the row is written with a zero Isp and ``frozen_converged = 0``
    rather than being dropped or filled in, so the gap is visible in the data.
    """
    from rocketcea.cea_obj import CEA_Obj

    cea = CEA_Obj(oxName=OX_NAME, fuelName=FUEL_NAME)
    pressure_psia = NOZZLE_CHEMISTRY_PC_PA / 6894.757293168361

    rows = []
    for of in _grid(OF_MIN, OF_MAX, NOZZLE_CHEMISTRY_OF_STEP):
        equilibrium = float(cea.get_IvacCstrTc(Pc=pressure_psia, MR=of, eps=REFERENCE_EPS)[0])
        frozen = float(cea.getFrozen_IvacCstrTc(Pc=pressure_psia, MR=of, eps=REFERENCE_EPS)[0])
        rows.append(
            {
                "of": of,
                "chamber_pressure_pa": NOZZLE_CHEMISTRY_PC_PA,
                "expansion_ratio": REFERENCE_EPS,
                "vacuum_isp_equilibrium_s": equilibrium,
                "vacuum_isp_frozen_s": frozen,
                "frozen_converged": 1 if frozen > 0.0 else 0,
            }
        )
    if verbose:
        converged = sum(int(row["frozen_converged"]) for row in rows)
        print(f"  frozen expansion converged at {converged} of {len(rows)} O/F points")
    return rows


def write_nozzle_chemistry(rows: list[dict[str, float]]) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    failures = [row["of"] for row in rows if not row["frozen_converged"]]
    span = f"O/F {min(failures):.2f} to {max(failures):.2f}" if failures else "nowhere"
    header = [
        "# Equilibrium vs frozen nozzle expansion, N2O/HTPB.",
        "# Generated by scripts/generate_thermochemistry_table.py -- do not edit by hand.",
        "#",
        "# Solver     : NASA CEA via rocketcea (Gordon & McBride, NASA RP-1311).",
        f"# Condition  : p_c = {NOZZLE_CHEMISTRY_PC_PA / 1e5:.1f} bar, eps = {REFERENCE_EPS:.1f}.",
        "#",
        "# Vacuum Isp is quoted because it isolates the expansion assumption from",
        "# any ambient-pressure bookkeeping.  CEA's frozen option fails to converge",
        f"# over {span}, where condensed carbon is present; those rows carry",
        "# vacuum_isp_frozen_s = 0 and frozen_converged = 0.  They are kept rather",
        "# than dropped so the gap is visible in the committed data.",
        "#",
    ]
    with NOZZLE_CHEMISTRY_PATH.open("w", encoding="utf-8", newline="\n") as handle:
        for line in header:
            handle.write(line + "\n")
        writer = csv.DictWriter(handle, fieldnames=NOZZLE_CHEMISTRY_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "of": f"{row['of']:.6f}",
                    "chamber_pressure_pa": f"{row['chamber_pressure_pa']:.6f}",
                    "expansion_ratio": f"{row['expansion_ratio']:.6f}",
                    "vacuum_isp_equilibrium_s": f"{row['vacuum_isp_equilibrium_s']:.6f}",
                    "vacuum_isp_frozen_s": f"{row['vacuum_isp_frozen_s']:.6f}",
                    "frozen_converged": int(row["frozen_converged"]),
                }
            )


def write_table(rows: list[dict[str, float]]) -> None:
    from rocketcea.cea_obj import fuelCards, oxCards

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    header = [
        "# N2O/HTPB equilibrium combustion properties for the chamber.",
        "# Generated by scripts/generate_thermochemistry_table.py -- do not edit by hand.",
        "#",
        "# Solver     : NASA CEA via rocketcea (Gordon & McBride, NASA RP-1311).",
        "# Oxidizer   : " + " ".join(line.strip() for line in oxCards[OX_NAME]),
        "# Fuel       : " + " ".join(line.strip() for line in fuelCards[FUEL_NAME]),
        "# Assumption : equilibrium composition; chamber (stagnation) properties.",
        "# O/F        : oxidizer mass flow / fuel mass flow.",
        "# c*         : IDEAL equilibrium characteristic velocity. A combustion",
        "#              efficiency is applied downstream and is NOT included here.",
        "# gamma      : CEA chamber GAMMAs (isentropic exponent).",
        "# molar mass : effective_molar_mass is CEA's M,(1/n), the total mixture mass",
        "#              per mole of gas -- the quantity that reproduces CEA's own sound",
        "#              speed when a condensed phase is present. gas_molar_mass is the",
        "#              gas-phase MW, kept only as a diagnostic; their ratio flags soot.",
        "# Units      : O/F [-], pressure [Pa], c* [m/s], gamma [-], T [K],",
        "#              molar mass [kg/mol], R [J/(kg K)].",
        "# NOT experimentally validated. Generic reduced-order educational study.",
    ]
    with TABLE_PATH.open("w", newline="", encoding="utf-8") as handle:
        for line in header:
            handle.write(line + "\n")
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: repr(row[key]) for key in COLUMNS})


def write_metadata(rows: list[dict[str, float]]) -> None:
    import rocketcea

    metadata = {
        "generator": "scripts/generate_thermochemistry_table.py",
        "solver": "NASA CEA via rocketcea",
        "rocketcea_version": rocketcea.__version__,
        "python_version": sys.version.split()[0],
        "generated_utc": datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d"),
        "oxidizer_card": OX_NAME,
        "fuel_card": FUEL_NAME,
        "of_grid": {"min": OF_MIN, "max": OF_MAX, "step": OF_STEP},
        "chamber_pressure_grid_bar": {"min": PC_MIN_BAR, "max": PC_MAX_BAR, "step": PC_STEP_BAR},
        "row_count": len(rows),
        "cstar_convention": "ideal equilibrium; combustion efficiency applied downstream",
        "note": (
            "This metadata file carries the generation date and tool versions and "
            "therefore changes on regeneration. The CSV itself contains no timestamp "
            "and regenerates byte-identically."
        ),
    }
    METADATA_PATH.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="regenerate in memory and confirm the committed CSV is unchanged",
    )
    args = parser.parse_args()

    print("Generating the N2O/HTPB equilibrium table with NASA CEA...")
    rows = generate()
    print("Generating the equilibrium/frozen expansion comparison...")
    nozzle_rows = generate_nozzle_chemistry()

    if args.check:
        stale = 0
        for path, writer, payload in (
            (TABLE_PATH, write_table, rows),
            (NOZZLE_CHEMISTRY_PATH, write_nozzle_chemistry, nozzle_rows),
        ):
            existing = path.read_text(encoding="utf-8")
            writer(payload)
            regenerated = path.read_text(encoding="utf-8")
            path.write_text(existing, encoding="utf-8")
            if existing == regenerated:
                print(f"  {path.name} is byte-identical to a fresh generation")
            else:
                print(f"  *** {path.name} DIFFERS from a fresh generation ***")
                stale = 1
        return stale

    write_table(rows)
    write_nozzle_chemistry(nozzle_rows)
    write_metadata(rows)
    print(f"  wrote {TABLE_PATH.relative_to(DATA_DIR.parent)} ({TABLE_PATH.stat().st_size} bytes)")
    print(f"  wrote {NOZZLE_CHEMISTRY_PATH.relative_to(DATA_DIR.parent)}")
    print(f"  wrote {METADATA_PATH.relative_to(DATA_DIR.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
