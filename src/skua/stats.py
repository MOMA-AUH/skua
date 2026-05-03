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
class StrandAwarePonStats:
    """Typed strand-aware summary for case vs panel-of-normals background."""

    case_counts: dict[str, int]
    normal_counts: dict[str, int]
    background_rate_by_channel: dict[str, float]
    expected_case_counts: dict[str, float]
    channel_score_contribution: dict[str, float]
    combined_score: float
    p_value: float
    method: str
    degrees_of_freedom: int
    pseudocount: float
    min_expected_count: float
    approximation_warning: bool

    def to_dict(self) -> dict[str, dict[str, float | int] | float]:
        """Return a JSON-serializable representation."""
        return {
            "case_counts": dict(self.case_counts),
            "normal_counts": dict(self.normal_counts),
            "background_rate_by_channel": dict(self.background_rate_by_channel),
            "expected_case_counts": dict(self.expected_case_counts),
            "channel_score_contribution": dict(self.channel_score_contribution),
            "combined_score": self.combined_score,
            "p_value": self.p_value,
            "method": self.method,
            "degrees_of_freedom": self.degrees_of_freedom,
            "pseudocount": self.pseudocount,
            "min_expected_count": self.min_expected_count,
            "approximation_warning": self.approximation_warning,
        }


def compute_strand_aware_pon_stats(
    case_evidence: AggregatedEvidence,
    normal_evidence: AggregatedEvidence,
    *,
    pseudocount: float = 0.5,
) -> StrandAwarePonStats:
    """Compute basic strand-aware background estimates and score contributions."""
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
    min_expected_count = min(expected_case_counts.values()) if expected_case_counts else 0.0
    # Rule of thumb: chi-square approximation is less reliable when expected counts are small.
    approximation_warning = min_expected_count < 5.0
    channel_score_contribution = {
        channel: (
            (case_counts[channel] - expected_case_counts[channel])
            / math.sqrt(expected_case_counts[channel] + pseudocount)
            if case_total > 0
            else 0.0
        )
        for channel in _CHANNELS
    }
    combined_score = math.sqrt(
        sum(contribution * contribution for contribution in channel_score_contribution.values())
    )
    chi2_stat = combined_score * combined_score
    # Survival function for Chi-square with k=4 degrees of freedom:
    # sf(x) = exp(-x/2) * (1 + x/2)
    p_value = math.exp(-chi2_stat / 2.0) * (1.0 + chi2_stat / 2.0)

    return StrandAwarePonStats(
        case_counts=case_counts,
        normal_counts=normal_counts,
        background_rate_by_channel=background_rate_by_channel,
        expected_case_counts=expected_case_counts,
        channel_score_contribution=channel_score_contribution,
        combined_score=combined_score,
        p_value=p_value,
        method="chi_square_4channel_approx",
        degrees_of_freedom=4,
        pseudocount=pseudocount,
        min_expected_count=min_expected_count,
        approximation_warning=approximation_warning,
    )