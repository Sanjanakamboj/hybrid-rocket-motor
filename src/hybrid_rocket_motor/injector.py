"""Reduced-order nitrous oxide injector flow models.

Scope
-----
This module models an injector as an **effective flow area** only.  It produces
no hole count, no hole diameter, no plate geometry and no manufacturing
information of any kind; ``C_d`` and ``A_eff`` are lumped generic parameters.

Sourced equations
-----------------
All three models are taken from

    B. S. Waxman, J. E. Zimmerman, B. J. Cantwell (Stanford University) and
    G. G. Zilliac (NASA Ames), "Mass Flow Rate and Isolation Characteristics of
    Injectors for Use with Self-Pressurizing Oxidizers in Hybrid Rockets",
    NASA NTRS 20190001326.

* **SPI** -- single-phase incompressible, their Eq. (2), the familiar "C_d A"
  form they attribute to Sutton & Biblarz::

      m_dot_SPI = C_d A sqrt(2 rho dP),      dP = P1 - P2

  Their Eq. (1) carries an extra ``1/(1 - (A/A_1)^2)`` velocity-of-approach term;
  Eq. (2) is the ``A << A_1`` limit, and they note the correction is "often
  wrapped into C_d".  This project uses Eq. (2).

* **HEM** -- homogeneous equilibrium, their Eqs. (3)-(4)::

      m_dot_HEM = C_d A rho_2 sqrt(2 (h_1 - h_2)),     s_2 = s_1

  i.e. the downstream state is found by isentropic expansion to the chamber
  pressure.  For nitrous oxide that downstream state is generally two-phase.

* **Dyer / NHNE** -- their Eqs. (8)-(9), the non-equilibrium blend first proposed
  by Dyer et al. and corrected by Solomon::

      kappa      = sqrt( (P1 - P2) / (Pv - P2) )
      m_dot_DYER = kappa/(1+kappa) m_dot_SPI + 1/(1+kappa) m_dot_HEM

  ``kappa`` compares the characteristic bubble-growth time with the liquid
  residence time in the orifice.

  A transcription note, stated rather than silently corrected: the printed
  Eq. (9) in the retrieved PDF renders as ``m_dot_DYER = A [ ... ]``, with a
  leading ``A`` that is an artefact of the large bracket in the original
  typesetting.  It cannot be a real factor -- ``m_dot_SPI`` and ``m_dot_HEM``
  already contain ``C_d A``, so multiplying by an area again would give
  kg m^2/s.  The dimensionally consistent weighted mean above is what is
  implemented, and it is the form universally quoted for the NHNE model.

Why the SPI model alone is not enough
-------------------------------------
Waxman et al. note that the "C_d A" equation "is often applied inappropriately to
injectors of high vapor pressure propellants such as nitrous oxide", because the
liquid flashes inside the orifice and the flow chokes; below a critical
downstream pressure the real mass flow becomes insensitive to backpressure, which
SPI cannot reproduce.  Both models are therefore provided and compared.

Units
-----
SI throughout: m^2, Pa, kg/m^3, J/kg, kg/s.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from .nitrous_properties import isentropic_downstream_enthalpy_j_kg
from .tank import TankState

__all__ = [
    "DYER_REFERENCE_DISCHARGE_COEFFICIENT",
    "Injector",
    "InjectorFlowStatus",
    "InjectorModel",
    "InjectorResult",
]

#: Discharge coefficient used by Dyer et al. when comparing their model against
#: Stanford hot-fire data, quoted by Waxman et al. as "a constant C_d of 0.66 ...
#: based on the average value for all of the injectors that were included".  A
#: generic literature value, not a measurement of any injector modelled here.
DYER_REFERENCE_DISCHARGE_COEFFICIENT = 0.66


class InjectorModel(Enum):
    """Which flow model to evaluate."""

    SPI = "SPI"
    HEM = "HEM"
    DYER_NHNE = "DYER_NHNE"


class InjectorFlowStatus(Enum):
    """Why an injector result is what it is."""

    FLOWING = "FLOWING"
    NO_FLOW_NO_PRESSURE_DROP = "NO_FLOW_NO_PRESSURE_DROP"
    NO_FLOW_NO_LIQUID = "NO_FLOW_NO_LIQUID"


@dataclass(frozen=True)
class InjectorResult:
    """One injector evaluation."""

    mass_flow_kg_s: float
    model: InjectorModel
    status: InjectorFlowStatus
    upstream_pressure_pa: float
    downstream_pressure_pa: float
    pressure_drop_pa: float
    spi_mass_flow_kg_s: float
    hem_mass_flow_kg_s: float
    non_equilibrium_parameter: float

    @property
    def is_flowing(self) -> bool:
        return self.status is InjectorFlowStatus.FLOWING


def _check_positive(name: str, value: float) -> float:
    v = float(value)
    if not math.isfinite(v):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if v <= 0.0:
        raise ValueError(f"{name} must be strictly positive, got {v!r}")
    return v


@dataclass(frozen=True)
class Injector:
    """A generic injector described only by an effective area and a ``C_d``.

    Parameters
    ----------
    effective_area_m2:
        Total geometric flow area ``A`` [m^2].  A lumped effective value; this
        project does not model or recommend any physical hole arrangement.
    discharge_coefficient:
        Dimensionless ``C_d``.  Must lie in ``(0, 1]``.
    model:
        Which of the three sourced models to evaluate.
    """

    effective_area_m2: float
    discharge_coefficient: float = DYER_REFERENCE_DISCHARGE_COEFFICIENT
    model: InjectorModel = InjectorModel.DYER_NHNE

    def __post_init__(self) -> None:
        _check_positive("effective_area_m2", self.effective_area_m2)
        cd = float(self.discharge_coefficient)
        if not math.isfinite(cd):
            raise ValueError(f"discharge_coefficient must be finite, got {cd!r}")
        if not (0.0 < cd <= 1.0):
            raise ValueError(f"discharge_coefficient must lie in (0, 1], got {cd!r}")
        if not isinstance(self.model, InjectorModel):
            raise TypeError(f"model must be an InjectorModel, got {self.model!r}")

    @property
    def effective_cd_area_m2(self) -> float:
        """The lumped product ``C_d A`` [m^2], which is what actually sets flow."""
        return self.discharge_coefficient * self.effective_area_m2

    # -- individual models -------------------------------------------------

    def spi_mass_flow_kg_s(self, liquid_density_kg_m3: float, pressure_drop_pa: float) -> float:
        """``m_dot_SPI = C_d A sqrt(2 rho dP)`` [kg/s]; zero for ``dP <= 0``."""
        rho = _check_positive("liquid_density_kg_m3", liquid_density_kg_m3)
        drop = float(pressure_drop_pa)
        if not math.isfinite(drop):
            raise ValueError(f"pressure_drop_pa must be finite, got {pressure_drop_pa!r}")
        if drop <= 0.0:
            return 0.0
        return self.effective_cd_area_m2 * math.sqrt(2.0 * rho * drop)

    def hem_mass_flow_kg_s(
        self,
        upstream_entropy_j_kg_k: float,
        upstream_enthalpy_j_kg: float,
        downstream_pressure_pa: float,
    ) -> float:
        """``m_dot_HEM = C_d A rho_2 sqrt(2 (h_1 - h_2))`` at ``s_2 = s_1`` [kg/s].

        Returns zero when the isentropic expansion would not lower the enthalpy,
        which is the no-flow limit.
        """
        h2, rho2 = isentropic_downstream_enthalpy_j_kg(
            upstream_entropy_j_kg_k, downstream_pressure_pa
        )
        delta_h = float(upstream_enthalpy_j_kg) - h2
        if delta_h <= 0.0:
            return 0.0
        return self.effective_cd_area_m2 * rho2 * math.sqrt(2.0 * delta_h)

    @staticmethod
    def non_equilibrium_parameter(
        upstream_pressure_pa: float, downstream_pressure_pa: float, vapour_pressure_pa: float
    ) -> float:
        """``kappa = sqrt((P1 - P2)/(Pv - P2))`` [-].

        For a *saturated* tank the upstream pressure equals the vapour pressure,
        so ``kappa`` is exactly 1 and the Dyer model reduces to an equal blend of
        SPI and HEM.  When the downstream pressure reaches or exceeds the vapour
        pressure there is no flashing potential at all; ``inf`` is returned, which
        drives the Dyer weighting to pure SPI -- the physically correct limit.
        """
        p1 = float(upstream_pressure_pa)
        p2 = float(downstream_pressure_pa)
        pv = float(vapour_pressure_pa)
        if not all(math.isfinite(v) for v in (p1, p2, pv)):
            raise ValueError("pressures must be finite")
        flash_margin = pv - p2
        if flash_margin <= 0.0:
            return math.inf
        drop = p1 - p2
        if drop <= 0.0:
            return 0.0
        return math.sqrt(drop / flash_margin)

    # -- combined evaluation -----------------------------------------------

    def mass_flow(self, tank_state: TankState, chamber_pressure_pa: float) -> InjectorResult:
        """Evaluate the selected model for liquid discharge into a chamber.

        No flow is produced when the tank holds no liquid, or when the chamber
        pressure has risen to meet the tank pressure.  Nothing is clipped or
        fabricated: a non-flowing condition is reported as such.
        """
        p2 = float(chamber_pressure_pa)
        if not math.isfinite(p2) or p2 < 0.0:
            raise ValueError(
                f"chamber_pressure_pa must be finite and non-negative, got {chamber_pressure_pa!r}"
            )
        if not isinstance(tank_state, TankState):
            raise TypeError(f"tank_state must be a TankState, got {type(tank_state).__name__}")

        p1 = tank_state.pressure_pa
        drop = p1 - p2

        def empty(status: InjectorFlowStatus) -> InjectorResult:
            return InjectorResult(
                mass_flow_kg_s=0.0,
                model=self.model,
                status=status,
                upstream_pressure_pa=p1,
                downstream_pressure_pa=p2,
                pressure_drop_pa=drop,
                spi_mass_flow_kg_s=0.0,
                hem_mass_flow_kg_s=0.0,
                non_equilibrium_parameter=math.nan,
            )

        if not tank_state.has_liquid:
            return empty(InjectorFlowStatus.NO_FLOW_NO_LIQUID)
        if drop <= 0.0:
            return empty(InjectorFlowStatus.NO_FLOW_NO_PRESSURE_DROP)

        spi = self.spi_mass_flow_kg_s(tank_state.liquid_density_kg_m3, drop)
        hem = self.hem_mass_flow_kg_s(
            tank_state.liquid_entropy_j_kg_k, tank_state.liquid_enthalpy_j_kg, p2
        )
        kappa = self.non_equilibrium_parameter(p1, p2, tank_state.pressure_pa)

        if self.model is InjectorModel.SPI:
            flow = spi
        elif self.model is InjectorModel.HEM:
            flow = hem
        elif math.isinf(kappa):
            flow = spi
        else:
            flow = (kappa / (1.0 + kappa)) * spi + (1.0 / (1.0 + kappa)) * hem

        return InjectorResult(
            mass_flow_kg_s=flow,
            model=self.model,
            status=InjectorFlowStatus.FLOWING,
            upstream_pressure_pa=p1,
            downstream_pressure_pa=p2,
            pressure_drop_pa=drop,
            spi_mass_flow_kg_s=spi,
            hem_mass_flow_kg_s=hem,
            non_equilibrium_parameter=kappa,
        )
