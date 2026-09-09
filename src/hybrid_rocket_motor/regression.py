"""Empirical hybrid fuel regression law ``r_dot = a * G_ox ** n``.

Physical background
-------------------
For a classical (non-liquefying, non-metallised) hybrid, Marxman's
diffusion-limited theory gives a local regression rate that scales as
``G ** 0.8`` with essentially no chamber-pressure dependence.  Experimenters
almost universally report a *space-time averaged* fit of the same power-law
form, ``r_dot = a * G_ox ** n``, with the flux exponent ``n`` treated as a free
parameter; values in the range 0.5-0.8 are reported for most hybrid systems.
See ``DESIGN.md`` for the full source audit.

Unit handling
-------------
Published correlations are almost never quoted in SI.  A :class:`PowerLawRegressionLaw`
therefore stores the coefficient *in the units of its source* together with an
explicit declaration of those units, and converts **once** to SI at construction
time.  The conversion is

    r_dot_SI = k_r * a_src * (G_SI / k_G) ** n
             = [k_r * a_src * k_G ** (-n)] * G_SI ** n
             = a_SI * G_SI ** n

where ``k_r`` is the number of m/s per source rate unit and ``k_G`` the number of
kg/(m^2 s) per source flux unit.  ``mm/s`` and ``m/s`` are never implicitly mixed,
and neither are ``g/(cm^2 s)`` and ``kg/(m^2 s)``: every public evaluation
function takes SI flux and returns SI regression rate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

__all__ = [
    "ILLUSTRATIVE_ANCHOR_FLUX_SI",
    "ILLUSTRATIVE_SENSITIVITY_LAWS",
    "REZAEI_2018_N2O_HTPB",
    "MassFluxUnit",
    "PowerLawRegressionLaw",
    "Provenance",
    "RegressionRateUnit",
    "illustrative_exponent_variant",
]


class MassFluxUnit(Enum):
    """Mass-flux unit, with its value in SI kg/(m^2 s) per unit."""

    KG_PER_M2_S = ("kg/(m^2 s)", 1.0)
    G_PER_CM2_S = ("g/(cm^2 s)", 10.0)  # 1e-3 kg / 1e-4 m^2 = 10 kg/(m^2 s)

    def __init__(self, label: str, si_per_unit: float) -> None:
        self.label = label
        self.si_per_unit = si_per_unit


class RegressionRateUnit(Enum):
    """Regression-rate unit, with its value in SI m/s per unit."""

    M_PER_S = ("m/s", 1.0)
    MM_PER_S = ("mm/s", 1.0e-3)

    def __init__(self, label: str, si_per_unit: float) -> None:
        self.label = label
        self.si_per_unit = si_per_unit


class Provenance(Enum):
    """Where a coefficient pair came from.

    ``SOURCED``
        Taken verbatim from a published, peer-reviewed correlation for this
        propellant combination.
    ``ILLUSTRATIVE``
        Constructed by this project purely to explore sensitivity.  Such
        coefficients are **not** experimental data and must never be presented
        as validated N2O/HTPB measurements.
    """

    SOURCED = "sourced"
    ILLUSTRATIVE = "ILLUSTRATIVE (not experimentally validated)"


@dataclass(frozen=True)
class PowerLawRegressionLaw:
    """Power-law regression correlation ``r_dot = a * G_ox ** n``.

    The coefficient is supplied in the units of the source publication and is
    converted once, at construction, into SI.  Evaluation is SI-in / SI-out.

    Parameters
    ----------
    coefficient_source_units:
        ``a`` exactly as printed by the source.
    exponent:
        ``n``, dimensionless.
    flux_unit:
        Unit in which the source's ``G_ox`` is expressed.
    rate_unit:
        Unit in which the source's ``r_dot`` is expressed.
    label:
        Short human-readable name, used in tables, legends and figures.
    provenance:
        :class:`Provenance` of the coefficient pair.
    reference:
        Bibliographic reference, or a note describing how an illustrative pair
        was constructed.
    valid_flux_range_si:
        ``(G_min, G_max)`` in kg/(m^2 s) over which the source reports data, or
        ``None`` if unknown.  Purely informational: evaluation is never clipped,
        but :meth:`is_within_validated_flux_range` lets callers flag
        extrapolation honestly.
    """

    coefficient_source_units: float
    exponent: float
    flux_unit: MassFluxUnit
    rate_unit: RegressionRateUnit
    label: str
    provenance: Provenance
    reference: str
    valid_flux_range_si: tuple[float, float] | None = None
    coefficient_si: float = field(init=False)

    def __post_init__(self) -> None:
        a = float(self.coefficient_source_units)
        n = float(self.exponent)
        if not math.isfinite(a) or a <= 0.0:
            raise ValueError(f"coefficient must be finite and positive, got {a!r}")
        if not math.isfinite(n) or n <= 0.0:
            raise ValueError(f"exponent must be finite and positive, got {n!r}")
        if self.valid_flux_range_si is not None:
            lo, hi = self.valid_flux_range_si
            if not (math.isfinite(lo) and math.isfinite(hi)):
                raise ValueError("valid_flux_range_si entries must be finite")
            if lo <= 0.0 or hi <= lo:
                raise ValueError(f"valid_flux_range_si must satisfy 0 < lo < hi, got {(lo, hi)!r}")
        # Single, explicit unit conversion (see module docstring).
        a_si = self.rate_unit.si_per_unit * a * self.flux_unit.si_per_unit ** (-n)
        object.__setattr__(self, "coefficient_si", a_si)

    # -- description -------------------------------------------------------

    @property
    def is_illustrative(self) -> bool:
        """True if these coefficients are illustrative rather than sourced."""
        return self.provenance is Provenance.ILLUSTRATIVE

    @property
    def source_equation(self) -> str:
        """The correlation as printed by its source, with units."""
        return (
            f"r_dot [{self.rate_unit.label}] = {self.coefficient_source_units:g} "
            f"* (G_ox [{self.flux_unit.label}])^{self.exponent:g}"
        )

    @property
    def si_equation(self) -> str:
        """The correlation after conversion to SI."""
        return (
            f"r_dot [m/s] = {self.coefficient_si:.6e} "
            f"* (G_ox [kg/(m^2 s)])^{self.exponent:g}"
        )

    # -- evaluation --------------------------------------------------------

    def regression_rate_si(self, oxidizer_mass_flux_si):
        """Regression rate [m/s] for oxidizer mass flux [kg/(m^2 s)].

        Accepts a scalar or any array-like; returns a float for scalar input and
        a :class:`numpy.ndarray` for array input.  ``G_ox = 0`` yields exactly
        ``0.0`` (no oxidizer flux, no diffusion-limited regression), which is the
        continuous limit of the power law for ``n > 0``.  Negative flux is
        rejected: it is not physically meaningful here and ``G ** n`` for
        non-integer ``n`` would otherwise produce NaN.
        """
        g = np.asarray(oxidizer_mass_flux_si, dtype=float)
        if np.any(~np.isfinite(g)):
            raise ValueError("oxidizer mass flux must be finite")
        if np.any(g < 0.0):
            raise ValueError("oxidizer mass flux must be non-negative")
        # np.power(0.0, n) is already 0.0 for n > 0; np.where keeps it explicit
        # and avoids any 0**n edge behaviour creeping in.
        rate = np.where(g > 0.0, self.coefficient_si * np.power(g, self.exponent), 0.0)
        if np.isscalar(oxidizer_mass_flux_si) or np.ndim(oxidizer_mass_flux_si) == 0:
            return float(rate)
        return rate

    def regression_rate_mm_s(self, oxidizer_mass_flux_si):
        """Regression rate in mm/s (readability only) for SI flux input."""
        return (
            self.regression_rate_si(oxidizer_mass_flux_si)
            / RegressionRateUnit.MM_PER_S.si_per_unit
        )

    def is_within_validated_flux_range(self, oxidizer_mass_flux_si) -> bool:
        """True if the flux lies inside the source's reported data range.

        Always ``False`` when no validated range is recorded, so that callers
        default to declaring extrapolation rather than assuming validity.
        """
        if self.valid_flux_range_si is None:
            return False
        lo, hi = self.valid_flux_range_si
        return bool(lo <= float(oxidizer_mass_flux_si) <= hi)


# ---------------------------------------------------------------------------
# Sourced correlation
# ---------------------------------------------------------------------------

#: Peer-reviewed N2O/HTPB space-time averaged correlation, Eq. (10) of
#: Rezaei, Soltani & Mohammadi (2018).  Printed there as
#: ``r_dot = 0.3977 * G_o ** 0.3667`` with r_dot in mm/s and G_o in g/(cm^2 s),
#: fitted to 18 firings of a 250 mm HTPB grain with 25-49 mm initial port
#: diameter over a mean oxidizer mass flux of 3.5-12 g/(cm^2 s).
REZAEI_2018_N2O_HTPB = PowerLawRegressionLaw(
    coefficient_source_units=0.3977,
    exponent=0.3667,
    flux_unit=MassFluxUnit.G_PER_CM2_S,
    rate_unit=RegressionRateUnit.MM_PER_S,
    label="Rezaei et al. (2018) N$_2$O/HTPB",
    provenance=Provenance.SOURCED,
    reference=(
        "H. Rezaei, M. R. Soltani, A. R. Mohammadi, 'Experimental study of fuel "
        "regression rate in an HTPB/N2O hybrid rocket motor', Scientia Iranica B "
        "25(1) (2018) 253-265, Eq. (10). DOI 10.24200/sci.2017.4317"
    ),
    # 3.5-12 g/(cm^2 s) -> 35-120 kg/(m^2 s)
    valid_flux_range_si=(35.0, 120.0),
)


# ---------------------------------------------------------------------------
# Illustrative sensitivity coefficients -- NOT experimental data
# ---------------------------------------------------------------------------

#: Flux at which illustrative exponent variants are pinned to the sourced law
#: [kg/(m^2 s)].  Chosen near the centre of the Milestone 1 study range.
ILLUSTRATIVE_ANCHOR_FLUX_SI = 70.0


def illustrative_exponent_variant(
    base_law: PowerLawRegressionLaw,
    exponent: float,
    anchor_flux_si: float = ILLUSTRATIVE_ANCHOR_FLUX_SI,
) -> PowerLawRegressionLaw:
    """Build an ILLUSTRATIVE law with a different flux exponent.

    The exponent reported by :data:`REZAEI_2018_N2O_HTPB` (``n = 0.3667``) sits
    below the 0.5-0.8 band that Karabeyoglu, Cantwell & Zilliac (2007) report for
    most hybrid systems.  To explore what that discrepancy implies, this helper
    returns a law with the requested exponent whose coefficient is fixed by the
    single stated condition

        r_dot_variant(anchor_flux) == r_dot_base(anchor_flux)

    i.e. the curve is *pivoted* about one disclosed anchor flux.  This is a
    transparent construction, **not** a calibration to N2O/HTPB data: no
    experimental measurement informs the resulting coefficient beyond the anchor
    value inherited from the base law.  The returned law is marked
    :attr:`Provenance.ILLUSTRATIVE` and carries no validated flux range.
    """
    anchor = float(anchor_flux_si)
    if not math.isfinite(anchor) or anchor <= 0.0:
        raise ValueError(f"anchor_flux_si must be finite and positive, got {anchor_flux_si!r}")
    n = float(exponent)
    if not math.isfinite(n) or n <= 0.0:
        raise ValueError(f"exponent must be finite and positive, got {exponent!r}")

    rate_at_anchor_si = base_law.regression_rate_si(anchor)
    a_si = rate_at_anchor_si / anchor**n

    return PowerLawRegressionLaw(
        coefficient_source_units=a_si,
        exponent=n,
        flux_unit=MassFluxUnit.KG_PER_M2_S,
        rate_unit=RegressionRateUnit.M_PER_S,
        label=f"ILLUSTRATIVE $n={n:g}$",
        provenance=Provenance.ILLUSTRATIVE,
        reference=(
            "ILLUSTRATIVE SENSITIVITY CASE - not experimental data. Exponent set "
            f"to {n:g} (classical/reported hybrid range 0.5-0.8); coefficient fixed "
            f"only by pivoting about G_ox = {anchor:g} kg/(m^2 s) on "
            f"'{base_law.label}'."
        ),
        valid_flux_range_si=None,
    )


#: Clearly-labelled illustrative exponent sweep.  Present in every study and
#: figure so that the sensitivity of the reduced-order model to the flux
#: exponent is visible, and never presented as N2O/HTPB measurement.
ILLUSTRATIVE_SENSITIVITY_LAWS: tuple[PowerLawRegressionLaw, ...] = tuple(
    illustrative_exponent_variant(REZAEI_2018_N2O_HTPB, n) for n in (0.5, 0.62, 0.8)
)
