"""Deterministic one-factor robustness analysis of the coupled model.

This is the Milestone 6 analysis layer. It adds **no new physics**. It runs the
existing frozen Milestone 4 and Milestone 5 models over a disclosed set of
one-at-a-time perturbations to the assumptions that are actually uncertain, and
records what changes.

What this is, and is not
------------------------
This is a **deterministic** sensitivity study. Each case changes exactly one
assumption to a stated alternative value and reruns the whole coupled model. No
probability distribution is assumed for any input, no sampling is performed, and
nothing here is a confidence interval. A "±10 %" below means two extra runs at
those two values -- not a claim that the true value lies in that band.

Invalid cases are reported, never repaired
------------------------------------------
If a perturbation drives the coupled state off the thermochemistry table or out
of the checked property range, the case is recorded with its terminal status and
excluded from the ranking with an explicit note. It is never clamped back into
validity, because doing so would silently replace the case with a different one.

Units
-----
SI throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from .blowdown import BlowdownResult, simulate_blowdown
from .chamber import CombustionProperties
from .geometry import GrainGeometry
from .injector import Injector, InjectorModel
from .nozzle import NozzleGeometry
from .regression import PowerLawRegressionLaw, Provenance
from .tank import NitrousTank, TankState
from .thermochemistry import DEFAULT_C_STAR_EFFICIENCY, ThermochemistryTable
from .variable_blowdown import VariableBlowdownResult, simulate_variable_blowdown

__all__ = [
    "CaseOutcome",
    "ModelKind",
    "SensitivityFamily",
    "SensitivityRanking",
    "conclusion_support",
    "rank_families",
    "summarise",
]


class ModelKind(Enum):
    """Which model produced a case."""

    FROZEN_PROPERTIES = "M4 frozen properties"
    VARIABLE_PROPERTIES = "M5 variable thermochemistry"


@dataclass(frozen=True)
class CaseOutcome:
    """Everything the robustness study records about one coupled run.

    ``valid`` is False when the run could not be completed on the stated
    assumptions -- for example because the coupled state left the tabulated
    thermochemistry rectangle. Such a case keeps its ``termination`` for the
    record and is excluded from every quantitative comparison.
    """

    label: str
    family: str
    model: ModelKind
    termination: str
    valid: bool
    burn_time_s: float
    initial_thrust_n: float
    final_thrust_n: float
    thrust_slope_percent: float
    peak_thrust_n: float
    mean_thrust_n: float
    total_impulse_n_s: float
    equivalent_specific_impulse_s: float
    initial_mixture_ratio: float
    final_mixture_ratio: float
    initial_chamber_pressure_pa: float
    final_chamber_pressure_pa: float
    oxidizer_consumed_kg: float
    fuel_consumed_kg: float
    #: Fractional change in delivered ``c*`` and in ``gamma`` across the firing
    #: part of the burn. ``nan`` for the frozen-property model, where both are
    #: constant by assumption rather than by physics.
    c_star_swing: float = math.nan
    gamma_swing: float = math.nan
    note: str = ""


def summarise(
    label: str,
    family: str,
    result: BlowdownResult | VariableBlowdownResult,
    *,
    model: ModelKind,
    note: str = "",
) -> CaseOutcome:
    """Reduce a coupled run to the recorded robustness quantities.

    The mixture ratio and chamber pressure are read at the first and last
    *firing* samples. A run that ends on liquid depletion has a final sample at
    which nothing is flowing, and reporting that sample's ``nan`` O/F as the
    "final" mixture ratio would be misleading.
    """
    firing = result.firing
    if not firing.any():
        # A perturbation can push the coupled state off the tabulated rectangle or
        # out of the checked property band immediately, so that nothing ever
        # fires. That is a real, reportable outcome of that assumption -- not an
        # error to be repaired -- so it is recorded as an invalid case carrying its
        # terminal status, and excluded from every quantitative comparison.
        statuses = {s.value if hasattr(s, "value") else str(s) for s in result.feed_status}
        return CaseOutcome(
            label=label,
            family=family,
            model=model,
            termination=result.termination.value,
            valid=False,
            burn_time_s=float(result.burn_time_s),
            initial_thrust_n=math.nan,
            final_thrust_n=math.nan,
            thrust_slope_percent=math.nan,
            peak_thrust_n=math.nan,
            mean_thrust_n=math.nan,
            total_impulse_n_s=math.nan,
            equivalent_specific_impulse_s=math.nan,
            initial_mixture_ratio=math.nan,
            final_mixture_ratio=math.nan,
            initial_chamber_pressure_pa=math.nan,
            final_chamber_pressure_pa=math.nan,
            oxidizer_consumed_kg=float(result.oxidizer_consumed_kg),
            fuel_consumed_kg=float(result.fuel_consumed_kg),
            note=note
            or ("never established a firing state; feed status " + "/".join(sorted(statuses))),
        )
    mixture = result.mixture_ratio[firing]
    pressure = result.chamber_pressure_pa[firing]
    thrust = result.thrust_n[firing]

    within_table = getattr(result, "stayed_within_table", True)
    termination = result.termination.value
    valid = bool(within_table)
    if not valid and not note:
        note = "left the tabulated thermochemistry rectangle"

    # Only the variable-property model has property histories to swing.
    c_star_swing, gamma_swing = math.nan, math.nan
    if hasattr(result, "c_star_m_s"):
        c_star = result.c_star_m_s[firing]
        gamma = result.gamma[firing]
        c_star_swing = float(c_star[-1] / c_star[0] - 1.0)
        gamma_swing = float(gamma[-1] / gamma[0] - 1.0)

    return CaseOutcome(
        label=label,
        family=family,
        model=model,
        termination=termination,
        valid=valid,
        burn_time_s=float(result.burn_time_s),
        initial_thrust_n=float(thrust[0]),
        final_thrust_n=float(thrust[-1]),
        thrust_slope_percent=100.0 * (float(thrust[-1]) / float(thrust[0]) - 1.0),
        peak_thrust_n=float(result.peak_thrust_n),
        mean_thrust_n=float(result.mean_thrust_n),
        total_impulse_n_s=float(result.total_impulse_n_s),
        equivalent_specific_impulse_s=float(result.equivalent_specific_impulse_s),
        initial_mixture_ratio=float(mixture[0]),
        final_mixture_ratio=float(mixture[-1]),
        initial_chamber_pressure_pa=float(pressure[0]),
        final_chamber_pressure_pa=float(pressure[-1]),
        oxidizer_consumed_kg=float(result.oxidizer_consumed_kg),
        fuel_consumed_kg=float(result.fuel_consumed_kg),
        c_star_swing=c_star_swing,
        gamma_swing=gamma_swing,
        note=note,
    )


@dataclass(frozen=True)
class SensitivityFamily:
    """One assumption, varied to one or more alternatives about a baseline."""

    name: str
    description: str
    baseline: CaseOutcome
    variants: tuple[CaseOutcome, ...]

    @property
    def valid_variants(self) -> tuple[CaseOutcome, ...]:
        return tuple(v for v in self.variants if v.valid)

    @property
    def invalid_variants(self) -> tuple[CaseOutcome, ...]:
        return tuple(v for v in self.variants if not v.valid)

    def swing_percent(self, metric: str) -> float:
        """Largest absolute % change from baseline across the valid variants.

        This is the disclosed ranking metric. It is a *range* over the tested
        alternatives, not a standard deviation and not an uncertainty.
        """
        reference = getattr(self.baseline, metric)
        if reference == 0.0 or not math.isfinite(reference):
            return math.nan
        swings = [abs(100.0 * (getattr(v, metric) / reference - 1.0)) for v in self.valid_variants]
        return max(swings) if swings else math.nan

    def signed_extremes(self, metric: str) -> tuple[float, float]:
        """Most negative and most positive % change from baseline [%, %]."""
        reference = getattr(self.baseline, metric)
        if reference == 0.0 or not math.isfinite(reference):
            return (math.nan, math.nan)
        changes = [100.0 * (getattr(v, metric) / reference - 1.0) for v in self.valid_variants]
        if not changes:
            return (math.nan, math.nan)
        return (min(changes), max(changes))


@dataclass(frozen=True)
class SensitivityRanking:
    """Families ordered by their effect on one metric."""

    metric: str
    entries: tuple[tuple[str, float], ...]

    @property
    def dominant(self) -> str:
        return self.entries[0][0] if self.entries else "none"


def rank_families(families: list[SensitivityFamily], metric: str) -> SensitivityRanking:
    """Order families by ``swing_percent`` on ``metric``, largest first.

    Families whose swing is undefined (no valid variant) are dropped rather than
    ranked at zero, which would understate them as "insensitive".
    """
    scored = [
        (family.name, family.swing_percent(metric))
        for family in families
        if math.isfinite(family.swing_percent(metric))
    ]
    scored.sort(key=lambda item: item[1], reverse=True)
    return SensitivityRanking(metric=metric, entries=tuple(scored))


@dataclass(frozen=True)
class ConclusionSupport:
    """How a qualitative claim fared across the tested deterministic family."""

    claim: str
    supporting: int
    total: int
    exceptions: tuple[str, ...]

    @property
    def verdict(self) -> str:
        """Classification. Deliberately avoids universal language.

        ``robust within tested deterministic family`` -- every tested case agrees.
        ``conditionally robust`` -- at least four fifths agree, exceptions named.
        ``model-sensitive`` -- otherwise.
        """
        if self.total == 0:
            return "untested"
        if self.supporting == self.total:
            return "robust within tested deterministic family"
        if self.supporting >= 0.8 * self.total:
            return "conditionally robust"
        return "model-sensitive"


def conclusion_support(claim: str, outcomes: list[CaseOutcome], predicate) -> ConclusionSupport:
    """Evaluate ``predicate`` over the valid outcomes and count support."""
    valid = [o for o in outcomes if o.valid]
    exceptions = tuple(o.label for o in valid if not predicate(o))
    return ConclusionSupport(
        claim=claim,
        supporting=len(valid) - len(exceptions),
        total=len(valid),
        exceptions=exceptions,
    )


# ---------------------------------------------------------------------------
# perturbation helpers -- each returns a modified INPUT, never a modified module
# ---------------------------------------------------------------------------


def scaled_regression_law(law: PowerLawRegressionLaw, factor: float) -> PowerLawRegressionLaw:
    """The same correlation with its coefficient scaled by ``factor``.

    Marked ILLUSTRATIVE, because a scaled coefficient is no longer the published
    fit: the provenance must not claim a source it no longer has.
    """
    return PowerLawRegressionLaw(
        coefficient_source_units=law.coefficient_source_units * factor,
        exponent=law.exponent,
        flux_unit=law.flux_unit,
        rate_unit=law.rate_unit,
        label=f"{law.label} x{factor:g} coefficient",
        reference=f"{law.reference} -- coefficient scaled by {factor:g} for sensitivity only",
        valid_flux_range_si=law.valid_flux_range_si,
        provenance=Provenance.ILLUSTRATIVE,
    )


def law_with_exponent(
    law: PowerLawRegressionLaw, exponent: float, *, anchor_flux_si: float
) -> PowerLawRegressionLaw:
    """The correlation re-fitted to a different exponent through a fixed anchor.

    Changing ``n`` alone would change the regression rate everywhere, so the
    comparison would confound "different exponent" with "different overall rate".
    The coefficient is therefore re-chosen so that the perturbed law reproduces
    the baseline law **exactly at ``anchor_flux_si``**, which is the reference
    operating flux. The comparison is then specifically about the *slope* of the
    correlation, and that convention is stated wherever the result is reported.
    """
    baseline_rate = law.regression_rate_si(anchor_flux_si)
    coefficient_si = baseline_rate / anchor_flux_si**exponent
    # Convert the SI coefficient back to the source units the constructor takes.
    # The source-to-SI factor depends on the EXPONENT (the flux unit conversion is
    # raised to the power n), so it has to be measured at the *new* exponent, not
    # carried over from the baseline law. A probe with a unit source coefficient
    # reports exactly that factor.
    probe = PowerLawRegressionLaw(
        coefficient_source_units=1.0,
        exponent=exponent,
        flux_unit=law.flux_unit,
        rate_unit=law.rate_unit,
        label="unit probe",
        reference="internal unit-conversion probe",
        provenance=Provenance.ILLUSTRATIVE,
    )
    return PowerLawRegressionLaw(
        coefficient_source_units=coefficient_si / probe.coefficient_si,
        exponent=exponent,
        flux_unit=law.flux_unit,
        rate_unit=law.rate_unit,
        label=f"{law.label} n={exponent:g} (anchored)",
        reference=(
            f"{law.reference} -- exponent varied to {exponent:g} with the coefficient "
            f"re-anchored at G_ox = {anchor_flux_si:g} kg/(m^2 s) for sensitivity only"
        ),
        valid_flux_range_si=law.valid_flux_range_si,
        provenance=Provenance.ILLUSTRATIVE,
    )


def run_frozen(
    tank: NitrousTank,
    state: TankState,
    injector: Injector,
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    combustion: CombustionProperties,
    nozzle: NozzleGeometry,
    ambient_pressure_pa: float,
    **kwargs,
) -> BlowdownResult:
    """Run the frozen Milestone 4 model, unmodified."""
    return simulate_blowdown(
        tank, state, injector, geometry, law, combustion, nozzle, ambient_pressure_pa, **kwargs
    )


def run_variable(
    tank: NitrousTank,
    state: TankState,
    injector: Injector,
    geometry: GrainGeometry,
    law: PowerLawRegressionLaw,
    table: ThermochemistryTable,
    nozzle: NozzleGeometry,
    ambient_pressure_pa: float,
    *,
    c_star_efficiency: float = DEFAULT_C_STAR_EFFICIENCY,
    **kwargs,
) -> VariableBlowdownResult:
    """Run the frozen Milestone 5 model, unmodified."""
    return simulate_variable_blowdown(
        tank,
        state,
        injector,
        geometry,
        law,
        table,
        nozzle,
        ambient_pressure_pa,
        c_star_efficiency=c_star_efficiency,
        **kwargs,
    )


def injector_with(area_m2: float, model: InjectorModel | None = None) -> Injector:
    """An injector with a different effective area and/or flow model."""
    if model is None:
        return Injector(area_m2)
    return Injector(area_m2, model=model)
