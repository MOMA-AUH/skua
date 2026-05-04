"""Statistical helpers for strand-aware PON evaluation."""

from dataclasses import dataclass
import math

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


def compute_stats(
    case_evidence: AggregatedEvidence,
    normal_evidence: AggregatedEvidence,
    *,
    rho: float = 1e-4,
    pseudocount: float = 0.5,
    prior_variant_probability: float = 0.5,
    mu_min: float = 1e-6,
    mu_max: float = 1 - 1e-6,
) -> Stats:
    """Compute a Shearwater-style beta-binomial Bayes-factor summary.

    The Bayes factor is oriented as artifact-vs-variant (null/alternative),
    consistent with the original deepSNV Shearwater code path. The reported
    posterior probability matches Shearwater's posterior for the null/artifact
    model M0, so lower values indicate stronger evidence for a true variant.
    """
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
