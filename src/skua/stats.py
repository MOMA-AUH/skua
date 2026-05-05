"""Statistical helpers for strand-aware PON evaluation."""

from dataclasses import dataclass
import math
import sys

from .evidence import AggregatedEvidence


_CHANNELS = (
    "alt_forward",
    "alt_reverse",
    "non_alt_forward",
    "non_alt_reverse",
)


@dataclass(frozen=True)
class Stats:
    """Typed strand-aware summary for case vs panel-of-normals background."""

    case_counts: dict[str, int]
    normal_counts: dict[str, int]
    background_rate_by_channel: dict[str, float]
    expected_case_counts: dict[str, float]
    bayes_factor: float
    log_bayes_factor_artifact_vs_variant: float
    artifact_posterior: float
    dispersion_rho: float
    pseudocount: float

    def to_dict(self) -> dict[str, dict[str, float | int] | float]:
        """Return a JSON-serializable representation."""
        return {
            "bayes_factor": self.bayes_factor,
            "artifact_posterior": self.artifact_posterior,
        }


def _bound(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _log_beta(a: float, b: float) -> float:
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def _logbb(x: int, n: int, mu_scaled: float, disp: float) -> float:
    """Log beta-binomial term (without binomial coefficient), following deepSNV."""
    return _log_beta(x + mu_scaled, n - x - mu_scaled + disp) - _log_beta(mu_scaled, disp - mu_scaled)


def estimate_rho(
    per_sample_evidences: list[AggregatedEvidence],
    *,
    truncate: float = 0.1,
    rho_min: float = 1e-4,
    rho_max: float = 0.1,
    pseudo: float = sys.float_info.epsilon,
) -> float:
    """Estimate beta-binomial overdispersion (rho) from per-sample PON evidence.

    Implements the method-of-moments estimator from Shearwater's estimateRho(),
    using a two-channel tensor-like representation of the available evidence:
    alt and non-alt, each combined across strands for rho estimation.
    Returns the alt-channel rho bounded to [rho_min, rho_max].
    """
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
    total_depth_all = sum(total_depth_by_sample)
    total_by_channel = [sum(channel_counts) for channel_counts in x_by_channel]

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
        nu = (xix + pseudo) / (total_depth_all + ncol * pseudo)

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
    truncate: float = 0.1,
    pseudocount: float = sys.float_info.epsilon,
    prior_variant_probability: float = 0.5,
    mu_min: float = 1e-6,
    mu_max: float = 1 - 1e-6,
) -> Stats:
    """Compute a Shearwater-style beta-binomial Bayes-factor summary.

    The Bayes factor is oriented as artifact-vs-variant (null/alternative),
    consistent with the original deepSNV Shearwater code path. The reported
    posterior probability matches Shearwater's posterior for the null/artifact
    model M0, so lower values indicate stronger evidence for a true variant.

    When ``per_sample_evidences`` is supplied, rho is estimated from the
    per-sample PON evidence using the Shearwater method-of-moments estimator
    (``estimate_rho``), replacing the fixed ``rho`` default.
    """
    if per_sample_evidences is not None:
        rho = estimate_rho(per_sample_evidences, truncate=truncate)
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

    if per_sample_evidences:
        # Shearwater-style truncation: exclude high-background PON samples
        # from the background pool used in Bayes-factor terms.
        masked = []
        for sample in per_sample_evidences:
            sample_alt = sample.alt_forward + sample.alt_reverse
            sample_depth = (
                sample.alt_forward
                + sample.alt_reverse
                + sample.non_alt_forward
                + sample.non_alt_reverse
            )
            sample_mu = (sample_alt + sys.float_info.epsilon) / (sample_depth + sys.float_info.epsilon)
            if sample_mu < truncate:
                masked.append(sample)

        X_fw = sum(sample.alt_forward for sample in masked)
        X_bw = sum(sample.alt_reverse for sample in masked)
        N_fw = sum(sample.alt_forward + sample.non_alt_forward for sample in masked)
        N_bw = sum(sample.alt_reverse + sample.non_alt_reverse for sample in masked)

    if case_total == 0:
        log_bayes_factor = 0.0
    else:
        rho = _bound(rho, 1e-6, 1 - 1e-6)
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

    if log_bayes_factor > 700:
        bayes_factor = float("inf")
    elif log_bayes_factor < -700:
        bayes_factor = 0.0
    else:
        bayes_factor = math.exp(log_bayes_factor)

    prior_variant_probability = _bound(prior_variant_probability, 1e-12, 1 - 1e-12)
    odds_variant = prior_variant_probability / (1.0 - prior_variant_probability)
    artifact_posterior = bayes_factor / (bayes_factor + odds_variant)

    return Stats(
        case_counts=case_counts,
        normal_counts=normal_counts,
        background_rate_by_channel=background_rate_by_channel,
        expected_case_counts=expected_case_counts,
        bayes_factor=bayes_factor,
        log_bayes_factor_artifact_vs_variant=log_bayes_factor,
        artifact_posterior=artifact_posterior,
        dispersion_rho=rho,
        pseudocount=pseudocount,
    )
