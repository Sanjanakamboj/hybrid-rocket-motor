"""Quasi-steady chamber mass balance with prescribed combustion properties.

Scope (Milestone 3)
-------------------
An *algebraic* chamber closure.  Given the oxidizer and fuel mass flows and a
**prescribed** characteristic velocity, the chamber pressure follows from the
definition of ``c*`` applied to a choked throat.

Nothing here models N2O tank state, vapour pressure, blowdown, injector flow,
feed-line losses, ignition, chamber filling, finite-rate combustion or
combustion instability.  The oxidizer mass flow remains an external input
exactly as in Milestones 1 and 2, so **chamber pressure does not feed back into
oxidizer flow**.

Governing relations
-------------------
Characteristic velocity is *defined* by ``c* = p_c A_t / m_dot`` (Seitzman,
Georgia Tech AE4451; the same definition follows from the NASA Glenn choked
mass-flow relation, see ``DESIGN.md`` §24).  Combining it with a quasi-steady
chamber mass balance ``m_dot_ox + m_dot_f = m_dot_nozzle`` gives

    p_c = (m_dot_ox + m_dot_f) c* / A_t

Assumptions, all of them consequential
--------------------------------------
* **Negligible chamber gas storage** — whatever is produced leaves through the
  throat in the same instant.  There is no chamber filling transient and no
  ``d(rho V)/dt`` term.
* **No ignition transient and no tail-off.**
* **No finite-rate combustion**; conversion is instantaneous and complete.
* **``c*`` is prescribed**, so combustion efficiency is *absorbed into* the
  prescribed value rather than modelled.
* **The feed system is unaffected by chamber pressure**, because ``m_dot_ox``
  is prescribed externally.  A real motor's oxidizer flow would depend on the
  chamber pressure it is feeding.

Units
-----
SI throughout: kg/s, m^2, Pa, K, J/(kg K), m/s.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from .regression import Provenance

__all__ = [
    "ILLUSTRATIVE_N2O_HTPB_COMBUSTION",
    "UNIVERSAL_GAS_CONSTANT_J_MOL_K",
    "ChamberState",
    "ChamberStatus",
    "CombustionProperties",
    "ideal_characteristic_velocity",
    "solve_chamber",
]

#: CODATA/SI defined molar gas constant [J/(mol K)].
UNIVERSAL_GAS_CONSTANT_J_MOL_K = 8.314462618


class ChamberStatus(Enum):
    """Whether the chamber is producing gas at this operating point."""

    FIRING = "FIRING"
    IDLE_NO_FLOW = "IDLE_NO_FLOW"


def _check_positive(name: str, value: float) -> float:
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if v <= 0.0:
        raise ValueError(f"{name} must be strictly positive, got {v!r}")
    return v


def _check_non_negative(name: str, value: float) -> float:
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if v < 0.0:
        raise ValueError(f"{name} must be non-negative, got {v!r}")
    return v


def ideal_characteristic_velocity(
    gamma: float, gas_constant_j_kg_k: float, chamber_temperature_k: float
) -> float:
    """Ideal ``c*`` implied by a set of thermochemical properties [m/s].

    Dividing the NASA Glenn choked mass-flow relation

        m_dot = (A_t p_c / sqrt(T_c)) sqrt(g/R) [(g+1)/2]^(-(g+1)/(2(g-1)))

    into ``p_c A_t`` gives

        c*_ideal = sqrt(R T_c / g) [(g+1)/2]^((g+1)/(2(g-1)))

    This is a *consistency* quantity: comparing it against a separately
    prescribed ``c*`` exposes what combustion efficiency the prescription
    implies, instead of leaving the two assumptions silently unrelated.
    """
    g = float(gamma)
    if not math.isfinite(g) or g <= 1.0:
        raise ValueError(f"gamma must be finite and greater than 1, got {gamma!r}")
    r_gas = _check_positive("gas_constant_j_kg_k", gas_constant_j_kg_k)
    t_c = _check_positive("chamber_temperature_k", chamber_temperature_k)
    exponent = (g + 1.0) / (2.0 * (g - 1.0))
    return math.sqrt(r_gas * t_c / g) * ((g + 1.0) / 2.0) ** exponent


@dataclass(frozen=True)
class CombustionProperties:
    """PRESCRIBED combustion gas properties.

    .. warning::

       These are **prescribed, illustrative assumptions**, not computed
       chemistry and not experimentally validated values for any real motor.
       Milestone 3 deliberately contains **no equilibrium-chemistry solver and
       no CEA interface**.  Every value here is an input the user may change,
       and the accompanying study quantifies how much the answers depend on it.

    Parameters
    ----------
    c_star_m_s:
        Delivered characteristic velocity ``c*`` [m/s].  Combustion efficiency
        is absorbed into this number.
    gamma:
        Ratio of specific heats of the combustion products [-].
    chamber_temperature_k:
        Stagnation (adiabatic flame) temperature ``T_c`` [K].
    molar_mass_kg_mol:
        Mean molar mass of the combustion products [kg/mol].
    label, reference:
        Human-readable identification and provenance note.
    """

    c_star_m_s: float
    gamma: float
    chamber_temperature_k: float
    molar_mass_kg_mol: float
    label: str
    reference: str
    provenance: Provenance = Provenance.ILLUSTRATIVE
    gas_constant_j_kg_k: float = field(init=False)

    def __post_init__(self) -> None:
        _check_positive("c_star_m_s", self.c_star_m_s)
        g = float(self.gamma)
        if not math.isfinite(g) or g <= 1.0:
            raise ValueError(f"gamma must be finite and greater than 1, got {self.gamma!r}")
        _check_positive("chamber_temperature_k", self.chamber_temperature_k)
        molar = _check_positive("molar_mass_kg_mol", self.molar_mass_kg_mol)
        object.__setattr__(
            self, "gas_constant_j_kg_k", UNIVERSAL_GAS_CONSTANT_J_MOL_K / molar
        )

    @classmethod
    def from_ideal_with_efficiency(
        cls,
        gamma: float,
        chamber_temperature_k: float,
        molar_mass_kg_mol: float,
        c_star_efficiency: float,
        label: str,
        reference: str,
    ) -> CombustionProperties:
        """Build with ``c*`` derived from the other properties and an efficiency.

        This is the preferred constructor: it makes ``c*`` a *consequence* of the
        thermochemical assumptions plus a stated combustion efficiency, rather
        than a fourth independent number that might quietly contradict them.
        """
        eta = float(c_star_efficiency)
        if not math.isfinite(eta) or not (0.0 < eta <= 1.0):
            raise ValueError(f"c_star_efficiency must lie in (0, 1], got {c_star_efficiency!r}")
        r_gas = UNIVERSAL_GAS_CONSTANT_J_MOL_K / _check_positive(
            "molar_mass_kg_mol", molar_mass_kg_mol
        )
        ideal = ideal_characteristic_velocity(gamma, r_gas, chamber_temperature_k)
        return cls(
            c_star_m_s=eta * ideal,
            gamma=gamma,
            chamber_temperature_k=chamber_temperature_k,
            molar_mass_kg_mol=molar_mass_kg_mol,
            label=label,
            reference=reference,
        )

    @property
    def ideal_c_star_m_s(self) -> float:
        """``c*`` implied by ``gamma``, ``R`` and ``T_c`` alone [m/s]."""
        return ideal_characteristic_velocity(
            self.gamma, self.gas_constant_j_kg_k, self.chamber_temperature_k
        )

    @property
    def c_star_efficiency(self) -> float:
        """Prescribed ``c*`` divided by the ideal ``c*`` implied by the rest [-].

        A value above 1 would mean the prescribed properties are mutually
        inconsistent (a real motor cannot beat its own ideal ``c*``).
        """
        return self.c_star_m_s / self.ideal_c_star_m_s

    @property
    def is_illustrative(self) -> bool:
        return self.provenance is Provenance.ILLUSTRATIVE

    def replace(self, **changes) -> CombustionProperties:
        """Return a copy with selected fields changed (used by sensitivity studies)."""
        allowed = {"c_star_m_s", "gamma", "chamber_temperature_k", "molar_mass_kg_mol",
                   "label", "reference"}
        unknown = set(changes) - allowed
        if unknown:
            raise ValueError(f"unknown field(s) for replace: {sorted(unknown)}")
        return CombustionProperties(
            c_star_m_s=changes.get("c_star_m_s", self.c_star_m_s),
            gamma=changes.get("gamma", self.gamma),
            chamber_temperature_k=changes.get(
                "chamber_temperature_k", self.chamber_temperature_k
            ),
            molar_mass_kg_mol=changes.get("molar_mass_kg_mol", self.molar_mass_kg_mol),
            label=changes.get("label", self.label),
            reference=changes.get("reference", self.reference),
        )


#: PRESCRIBED, ILLUSTRATIVE combustion properties for the generic N2O/HTPB study.
#:
#: Selection rule, fixed before any thrust was computed and not adjusted after
#: (see ``DESIGN.md`` §25):
#:
#: * ``gamma = 1.20``            - round value typical of hot rocket exhaust;
#: * ``T_c = 2600 K``            - round value below the sub-3000 K peak that
#:   Tarifa & Pizzuti (EUCASS 2019-488) report for N2O-based hybrids, reduced
#:   because this model runs fuel-rich at O/F ~ 2-2.5;
#: * ``M = 0.022 kg/mol``        - round value;
#: * ``eta_c* = 0.96``           - midpoint of the 94-98 % combustion efficiency
#:   measured by Rezaei, Soltani & Mohammadi (2018) for HTPB/N2O.
#:
#: The resulting delivered ``c*`` is a *consequence* of those four choices, not a
#: target. That it lands inside the 1403-1587 m/s the same study measured is a
#: weak consistency check, **not** a validation.
ILLUSTRATIVE_N2O_HTPB_COMBUSTION = CombustionProperties.from_ideal_with_efficiency(
    gamma=1.20,
    chamber_temperature_k=2600.0,
    molar_mass_kg_mol=0.022,
    c_star_efficiency=0.96,
    label="ILLUSTRATIVE prescribed N$_2$O/HTPB combustion properties",
    reference=(
        "PRESCRIBED / ILLUSTRATIVE - NOT experimentally validated and NOT computed "
        "chemistry. gamma = 1.20, T_c = 2600 K, M = 0.022 kg/mol chosen as round "
        "values; c* efficiency 0.96 is the midpoint of the 94-98% measured by "
        "Rezaei, Soltani & Mohammadi, Scientia Iranica B 25(1) (2018) 253-265. "
        "Chamber temperature bounded by Tarifa & Pizzuti, EUCASS 2019-488, who "
        "report sub-3000 K adiabatic flame temperatures for N2O-based hybrids."
    ),
)


@dataclass(frozen=True)
class ChamberState:
    """Quasi-steady chamber solution at one instant."""

    oxidizer_mass_flow_kg_s: float
    fuel_mass_flow_kg_s: float
    total_mass_flow_kg_s: float
    mixture_ratio: float
    chamber_pressure_pa: float
    c_star_m_s: float
    throat_area_m2: float
    status: ChamberStatus

    @property
    def is_firing(self) -> bool:
        return self.status is ChamberStatus.FIRING


def solve_chamber(
    combustion: CombustionProperties,
    throat_area_m2: float,
    oxidizer_mass_flow_kg_s: float,
    fuel_mass_flow_kg_s: float,
) -> ChamberState:
    """Quasi-steady chamber pressure from the prescribed ``c*`` closure.

    ``p_c = (m_dot_ox + m_dot_f) c* / A_t``.

    With no propellant flow at all the chamber is reported as
    :attr:`ChamberStatus.IDLE_NO_FLOW` with ``p_c = 0``: the quasi-steady balance
    has no gas to sustain any pressure, and no firing state is fabricated.  The
    mixture ratio is ``NaN`` whenever there is no fuel flow to form a ratio with,
    matching the Milestone 2 convention.
    """
    a_t = _check_positive("throat_area_m2", throat_area_m2)
    c_star = _check_positive("c_star_m_s", combustion.c_star_m_s)
    m_ox = _check_non_negative("oxidizer_mass_flow_kg_s", oxidizer_mass_flow_kg_s)
    m_f = _check_non_negative("fuel_mass_flow_kg_s", fuel_mass_flow_kg_s)

    m_total = m_ox + m_f
    of_ratio = m_ox / m_f if m_f > 0.0 else math.nan

    if m_total == 0.0:
        return ChamberState(
            oxidizer_mass_flow_kg_s=m_ox,
            fuel_mass_flow_kg_s=m_f,
            total_mass_flow_kg_s=0.0,
            mixture_ratio=math.nan,
            chamber_pressure_pa=0.0,
            c_star_m_s=c_star,
            throat_area_m2=a_t,
            status=ChamberStatus.IDLE_NO_FLOW,
        )

    return ChamberState(
        oxidizer_mass_flow_kg_s=m_ox,
        fuel_mass_flow_kg_s=m_f,
        total_mass_flow_kg_s=m_total,
        mixture_ratio=of_ratio,
        chamber_pressure_pa=m_total * c_star / a_t,
        c_star_m_s=c_star,
        throat_area_m2=a_t,
        status=ChamberStatus.FIRING,
    )
