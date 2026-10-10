"""Statistical helpers for strand-aware PON evaluation."""

from dataclasses import dataclass
from enum import Enum
import math
import sys

from .evidence import AggregatedEvidence


_CHANNELS = (
    "alt_forward",
    "alt_reverse",
    "non_alt_forward",
    "non_alt_reverse",
)

DEFAULT_TRUNCATE = 0.1


def _validate_positive_finite(value: float, *, name: str) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and > 0")


def _validate_probability_bounds(
    lower: float, upper: float, *, lower_name: str, upper_name: str,
) -> None:
    if not (
        math.isfinite(lower) and math.isfinite(upper) and 0.0 < lower <= upper < 1.0
    ):
        raise ValueError(
            f"{lower_name} and {upper_name} must be finite with "
            f"0 < {lower_name} <= {upper_name} < 1"
        )


def _validate_model_parameters(
    *,
    truncate: float | None = None,
    pseudocount: float | None = None,
    prior_artifact_probability: float | None = None,
) -> None:
    """Share statistical parameter contracts between annotation and scoring."""
    if truncate is not None and (
        not math.isfinite(truncate) or not 0.0 < truncate <= 1.0
    ):
        raise ValueError("truncate must be greater than 0 and no greater than 1, and finite")
    if pseudocount is not None:
        _validate_positive_finite(pseudocount, name="pseudocount")
    if prior_artifact_probability is not None and (
        not math.isfinite(prior_artifact_probability)
        or not 0.0 < prior_artifact_probability < 1.0
    ):
        raise ValueError("prior_artifact_probability must be finite and between 0 and 1")


class AssessmentStatus(str, Enum):
    """Whether the available evidence meets the configured assessment policy."""

    ASSESSED = "ASSESSED"
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"


@dataclass(frozen=True)
class AssessmentThresholds:
    """Inclusive usable-depth minima; zero disables sample/strand requirements.

    Defaults only exclude absent case or normal evidence. They are not
    assay-validated coverage requirements. Normal limits apply after truncation.
    """

    min_case_depth: int = 1
    min_normal_depth: int = 1
    min_normal_samples: int = 0
    min_case_strand_depth: int = 0
    min_normal_strand_depth: int = 0

    def __post_init__(self) -> None:
        for name, minimum in (
            ("min_case_depth", 1),
            ("min_normal_depth", 1),
            ("min_normal_samples", 0),
            ("min_case_strand_depth", 0),
            ("min_normal_strand_depth", 0),
        ):
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")


def _assessment_reasons(
    case: AggregatedEvidence,
    normal: AggregatedEvidence,
    normal_sample_count: int | None,
    thresholds: AssessmentThresholds,
) -> tuple[str, ...]:
    case_forward = case.alt_forward + case.non_alt_forward
    case_reverse = case.alt_reverse + case.non_alt_reverse
    normal_forward = normal.alt_forward + normal.non_alt_forward
    normal_reverse = normal.alt_reverse + normal.non_alt_reverse
    reasons = []
    if case_forward + case_reverse < thresholds.min_case_depth:
        reasons.append("CASE_DEPTH")
    if normal_forward + normal_reverse < thresholds.min_normal_depth:
        reasons.append("NORMAL_DEPTH")
    if thresholds.min_normal_samples > 0:
        if normal_sample_count is None:
            reasons.append("NORMAL_SAMPLE_COUNT_UNAVAILABLE")
        elif normal_sample_count < thresholds.min_normal_samples:
            reasons.append("NORMAL_SAMPLE_COUNT")
    if min(case_forward, case_reverse) < thresholds.min_case_strand_depth:
        reasons.append("CASE_STRAND_DEPTH")
    if min(normal_forward, normal_reverse) < thresholds.min_normal_strand_depth:
        reasons.append("NORMAL_STRAND_DEPTH")
    return tuple(reasons)


@dataclass(frozen=True)
class Stats:
    """Case/PON summary with model scores only for assessed evidence.

    ``assessment_reasons`` lists all unmet requirements, or is empty when
    ``assessment_status`` is ASSESSED. The log Bayes factor and artifact
    posterior are None when evidence is insufficient; counts remain available.
    """

    case_counts: dict[str, int]
    normal_counts: dict[str, int]
    background_rate_by_channel: dict[str, float]
    expected_case_counts: dict[str, float]
    log_bayes_factor_artifact_vs_variant: float | None
    artifact_posterior: float | None
    dispersion_rho: float
    pseudocount: float
    assessment_status: AssessmentStatus
    assessment_reasons: tuple[str, ...]


def _bound(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _log_beta(a: float, b: float) -> float:
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def _logbb(x: int, n: int, mu_scaled: float, disp: float) -> float:
    """Log beta-binomial term (without binomial coefficient), following deepSNV."""
    return _log_beta(x + mu_scaled, n - x - mu_scaled + disp) - _log_beta(mu_scaled, disp - mu_scaled)


def truncated_normal_evidences(
    per_sample_evidences: list[AggregatedEvidence],
    *,
    truncate: float = DEFAULT_TRUNCATE,
    epsilon: float = sys.float_info.epsilon,
) -> list[AggregatedEvidence]:
    """Return per-sample normal evidences retained by the truncation rule."""
    _validate_model_parameters(truncate=truncate)
    _validate_positive_finite(epsilon, name="epsilon")
    return [
        sample
        for sample in per_sample_evidences
        if (
            (
                sample.alt_forward
                + sample.alt_reverse
                + epsilon
            )
            /
            (
                sample.alt_forward
                + sample.alt_reverse
                + sample.non_alt_forward
                + sample.non_alt_reverse
                + epsilon
            )
        )
        < truncate
    ]


def aggregate_evidence(evidences: list[AggregatedEvidence]) -> AggregatedEvidence:
    """Aggregate a list of evidence objects into one strand-aware summary."""
    unusable_by_reason: dict = {}
    for evidence in evidences:
        for reason, count in evidence.unusable_by_reason.items():
            unusable_by_reason[reason] = unusable_by_reason.get(reason, 0) + count

    return AggregatedEvidence(
        alt_forward=sum(evidence.alt_forward for evidence in evidences),
        alt_reverse=sum(evidence.alt_reverse for evidence in evidences),
        non_alt_forward=sum(evidence.non_alt_forward for evidence in evidences),
        non_alt_reverse=sum(evidence.non_alt_reverse for evidence in evidences),
        usable=sum(evidence.usable for evidence in evidences),
        unusable=sum(evidence.unusable for evidence in evidences),
        unusable_by_reason=unusable_by_reason,
    )


def estimate_rho(
    per_sample_evidences: list[AggregatedEvidence],
    *,
    truncate: float = DEFAULT_TRUNCATE,
    rho_min: float = 1e-4,
    rho_max: float = 0.1,
    pseudo: float = sys.float_info.epsilon,
) -> float:
    """Estimate beta-binomial overdispersion (rho) from per-sample PON evidence.

    Adapts the method-of-moments estimator from Shearwater's estimateRho(),
    using a two-channel tensor-like representation of the available evidence:
    alt and non-alt, each combined across strands for rho estimation.
    Each channel's mean uses counts and depths from the same retained normals.
    Returns the alt-channel rho bounded to [rho_min, rho_max].
    """
    _validate_model_parameters(truncate=truncate)
    _validate_positive_finite(pseudo, name="pseudo")
    _validate_probability_bounds(
        rho_min, rho_max, lower_name="rho_min", upper_name="rho_max",
    )
    if len(per_sample_evidences) < 2:
        return rho_min

    ncol = 2
    total_depth_by_sample = [
        sample.alt_forward
        + sample.alt_reverse
        + sample.non_alt_forward
        + sample.non_alt_reverse
        for sample in per_sample_evidences
    ]
    x_by_channel = [
        [sample.alt_forward + sample.alt_reverse for sample in per_sample_evidences],
        [sample.non_alt_forward + sample.non_alt_reverse for sample in per_sample_evidences],
    ]
    rho_by_channel: list[float] = []
    for channel_index in range(ncol):
        mu_values = [
            (x_by_channel[channel_index][sample_index] + pseudo)
            / (total_depth_by_sample[sample_index] + ncol * pseudo)
            for sample_index in range(len(per_sample_evidences))
        ]
        included = [mu_value < truncate for mu_value in mu_values]
        included_count = sum(included)
        if included_count < 2:
            rho_by_channel.append(rho_min)
            continue

        xix = sum(
            x_by_channel[channel_index][sample_index]
            for sample_index in range(len(per_sample_evidences))
            if included[sample_index]
        )
        retained_depth = sum(
            total_depth_by_sample[sample_index]
            for sample_index in range(len(per_sample_evidences))
            if included[sample_index]
        )
        nu = (xix + pseudo) / (retained_depth + ncol * pseudo)

        valid_depths = [
            total_depth_by_sample[sample_index]
            for sample_index in range(len(per_sample_evidences))
            if included[sample_index] and total_depth_by_sample[sample_index] > 0
        ]
        valid_mu = [
            mu_values[sample_index]
            for sample_index in range(len(per_sample_evidences))
            if included[sample_index] and total_depth_by_sample[sample_index] > 0
        ]
        if included_count < 2 or not valid_depths:
            rho_by_channel.append(rho_min)
            continue

        sum_valid_depths = sum(valid_depths)
        s2 = (
            included_count
            * sum(
                valid_depths[value_index] * (valid_mu[value_index] - nu) ** 2
                for value_index in range(len(valid_depths))
            )
            / ((included_count - 1) * sum_valid_depths)
        )

        sum_inv_nix = sum(1.0 / depth for depth in valid_depths)
        denom = included_count - sum_inv_nix
        if denom <= 0 or nu <= 0.0 or nu >= 1.0:
            rho_by_channel.append(rho_min)
            continue

        rho_hat = (
            included_count * (s2 / nu / (1.0 - nu)) - sum_inv_nix
        ) / denom
        if not math.isfinite(rho_hat):
            rho_by_channel.append(rho_min)
            continue

        rho_hat = _bound(rho_hat, 0.0, 1.0)
        rho_hat = _bound(rho_hat, rho_min, rho_max)
        rho_by_channel.append(rho_hat)

    return rho_by_channel[0]


def compute_stats(
    case_evidence: AggregatedEvidence,
    normal_evidence: AggregatedEvidence,
    *,
    rho: float = 1e-4,
    per_sample_evidences: list[AggregatedEvidence] | None = None,
    truncate: float = DEFAULT_TRUNCATE,
    pseudocount: float = sys.float_info.epsilon,
    prior_artifact_probability: float = 0.5,
    mu_min: float = 1e-6,
    mu_max: float = 1 - 1e-6,
    assessment_thresholds: AssessmentThresholds = AssessmentThresholds(),
) -> Stats:
    """Compute a Shearwater-style beta-binomial Bayes-factor summary.

    The Bayes factor is oriented as artifact-vs-variant (null/alternative),
    consistent with the original deepSNV Shearwater code path. The reported
    posterior probability matches Shearwater's posterior for the null/artifact
    model M0, so lower values indicate stronger evidence for a true variant.

    When ``per_sample_evidences`` is supplied, rho is estimated from the
    per-sample PON evidence using the Shearwater method-of-moments estimator
    (``estimate_rho``), replacing the fixed ``rho`` default.
    The retained per-sample pool also replaces ``normal_evidence`` for counts,
    background summaries, scoring, and eligibility, including for an empty list.

    ``assessment_thresholds`` controls whether model scores are available.
    Ineligible evidence returns None for the log Bayes factor and posterior;
    counts and diagnostics remain available. Eligible scores are unchanged.
    Defaults exclude zero case/normal depth without imposing assay-specific
    sample-count or strand requirements. An aggregate-only pool cannot satisfy
    a positive minimum sample count: its reason is NORMAL_SAMPLE_COUNT_UNAVAILABLE.

    All parameters must be finite, with pseudocount > 0, 0 < truncate <= 1,
    0 < rho < 1, 0 < prior_artifact_probability < 1, and
    0 < mu_min <= mu_max < 1. Validation also applies at zero depth and before
    replacing rho with a per-sample estimate.
    """
    _validate_model_parameters(
        truncate=truncate,
        pseudocount=pseudocount,
        prior_artifact_probability=prior_artifact_probability,
    )
    if not math.isfinite(rho) or not 0.0 < rho < 1.0:
        raise ValueError("rho must be finite and between 0 and 1")
    _validate_probability_bounds(mu_min, mu_max, lower_name="mu_min", upper_name="mu_max")
    normal_sample_count = None
    if per_sample_evidences is not None:
        rho = estimate_rho(per_sample_evidences, truncate=truncate)
        retained_normals = truncated_normal_evidences(per_sample_evidences, truncate=truncate)
        normal_sample_count = len(retained_normals)
        normal_evidence = aggregate_evidence(retained_normals)
    case_counts = {
        "alt_forward": case_evidence.alt_forward,
        "alt_reverse": case_evidence.alt_reverse,
        "non_alt_forward": case_evidence.non_alt_forward,
        "non_alt_reverse": case_evidence.non_alt_reverse,
    }
    normal_counts = {
        "alt_forward": normal_evidence.alt_forward,
        "alt_reverse": normal_evidence.alt_reverse,
        "non_alt_forward": normal_evidence.non_alt_forward,
        "non_alt_reverse": normal_evidence.non_alt_reverse,
    }

    case_total = sum(case_counts.values())
    normal_total = sum(normal_counts.values())

    background_rate_by_channel = {
        channel: (normal_counts[channel] / normal_total) if normal_total > 0 else 0.0
        for channel in _CHANNELS
    }
    expected_case_counts = {
        channel: case_total * background_rate_by_channel[channel]
        for channel in _CHANNELS
    }

    x_fw = case_counts["alt_forward"]
    x_bw = case_counts["alt_reverse"]
    n_fw = x_fw + case_counts["non_alt_forward"]
    n_bw = x_bw + case_counts["non_alt_reverse"]

    X_fw = normal_counts["alt_forward"]
    X_bw = normal_counts["alt_reverse"]
    N_fw = X_fw + normal_counts["non_alt_forward"]
    N_bw = X_bw + normal_counts["non_alt_reverse"]

    assessment_reasons = _assessment_reasons(
        case_evidence, normal_evidence, normal_sample_count, assessment_thresholds,
    )
    log_bayes_factor: float | None = None
    artifact_posterior: float | None = None
    if case_total > 0:
        rho = _bound(rho, 1e-6, 1 - 1e-6)
    if not assessment_reasons:
        disp = (1.0 - rho) / rho

        mu = _bound(
            (x_fw + x_bw + pseudocount) / (n_fw + n_bw + 2.0 * pseudocount),
            mu_min,
            mu_max,
        )
        nu0_fw = _bound(
            (X_fw + x_fw + pseudocount) / (N_fw + n_fw + 2.0 * pseudocount),
            mu_min,
            mu_max,
        )
        nu0_bw = _bound(
            (X_bw + x_bw + pseudocount) / (N_bw + n_bw + 2.0 * pseudocount),
            mu_min,
            mu_max,
        )
        nu_fw = _bound((X_fw + pseudocount) / (N_fw + 2.0 * pseudocount), mu_min, mu_max)
        nu_bw = _bound((X_bw + pseudocount) / (N_bw + 2.0 * pseudocount), mu_min, mu_max)

        # Shearwater floor: prevent variant-rate mu from dropping below
        # strand-specific null-rate estimates.
        mu = max(mu, nu0_fw, nu0_bw)

        mu_scaled = mu * disp
        nu0_fw_scaled = nu0_fw * disp
        nu0_bw_scaled = nu0_bw * disp
        nu_fw_scaled = nu_fw * disp
        nu_bw_scaled = nu_bw * disp

        # AND-model style Bayes factor terms from deepSNV Shearwater formulation.
        log_bayes_factor = (
            _logbb(x_fw, n_fw, nu0_fw_scaled, disp)
            + _logbb(X_fw, N_fw, nu0_fw_scaled, disp)
            + _logbb(x_bw, n_bw, nu0_bw_scaled, disp)
            + _logbb(X_bw, N_bw, nu0_bw_scaled, disp)
            - _logbb(x_fw, n_fw, mu_scaled, disp)
            - _logbb(X_fw, N_fw, nu_fw_scaled, disp)
            - _logbb(x_bw, n_bw, mu_scaled, disp)
            - _logbb(X_bw, N_bw, nu_bw_scaled, disp)
        )

        prior_artifact_probability = _bound(
            prior_artifact_probability,
            1e-12,
            1 - 1e-12,
        )
        odds_artifact = prior_artifact_probability / (1.0 - prior_artifact_probability)
        log_posterior_odds_artifact = log_bayes_factor + math.log(odds_artifact)
        if log_posterior_odds_artifact >= 0:
            exp_neg_delta = math.exp(-log_posterior_odds_artifact)
            artifact_posterior = 1.0 / (1.0 + exp_neg_delta)
        else:
            exp_delta = math.exp(log_posterior_odds_artifact)
            artifact_posterior = exp_delta / (1.0 + exp_delta)

    return Stats(
        case_counts=case_counts,
        normal_counts=normal_counts,
        background_rate_by_channel=background_rate_by_channel,
        expected_case_counts=expected_case_counts,
        log_bayes_factor_artifact_vs_variant=log_bayes_factor,
        artifact_posterior=artifact_posterior,
        dispersion_rho=rho,
        pseudocount=pseudocount,
        assessment_status=(
            AssessmentStatus.INSUFFICIENT_EVIDENCE
            if assessment_reasons else AssessmentStatus.ASSESSED
        ),
        assessment_reasons=assessment_reasons,
    )
