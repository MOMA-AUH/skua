from skua.evidence import AggregatedEvidence
from skua.stats import Stats, compute_stats


def test_compute_stats_returns_typed_background_and_score() -> None:
    case_evidence = AggregatedEvidence(
        alt_forward=8,
        alt_reverse=0,
        non_alt_forward=2,
        non_alt_reverse=0,
        usable=10,
        unusable=0,
        unusable_by_reason={},
    )
    normal_evidence = AggregatedEvidence(
        alt_forward=1,
        alt_reverse=1,
        non_alt_forward=9,
        non_alt_reverse=9,
        usable=20,
        unusable=0,
        unusable_by_reason={},
    )

    stats = compute_stats(case_evidence, normal_evidence)

    assert isinstance(stats, Stats)
    assert stats.case_counts["alt_forward"] == 8
    assert stats.normal_counts["non_alt_forward"] == 9
    assert stats.background_rate_by_channel["alt_forward"] == 0.05
    assert stats.expected_case_counts["alt_forward"] == 0.5
    assert stats.bayes_factor >= 0.0
    assert 0.0 <= stats.artifact_posterior <= 1.0
    assert stats.dispersion_rho == 1e-4
    assert stats.pseudocount == 0.5


def test_compute_stats_is_stable_for_zero_depth() -> None:
    zero_evidence = AggregatedEvidence(
        alt_forward=0,
        alt_reverse=0,
        non_alt_forward=0,
        non_alt_reverse=0,
        usable=0,
        unusable=0,
        unusable_by_reason={},
    )

    stats = compute_stats(zero_evidence, zero_evidence)

    assert stats.background_rate_by_channel == {
        "alt_forward": 0.0,
        "alt_reverse": 0.0,
        "non_alt_forward": 0.0,
        "non_alt_reverse": 0.0,
    }
    assert stats.expected_case_counts == {
        "alt_forward": 0.0,
        "alt_reverse": 0.0,
        "non_alt_forward": 0.0,
        "non_alt_reverse": 0.0,
    }
    assert stats.bayes_factor == 1.0
    assert stats.log_bayes_factor_artifact_vs_variant == 0.0
    assert stats.artifact_posterior == 0.5
    assert stats.dispersion_rho == 1e-4
    assert stats.pseudocount == 0.5


def test_compute_stats_null_posterior_decreases_with_stronger_signal() -> None:
    normal_evidence = AggregatedEvidence(
        alt_forward=1,
        alt_reverse=1,
        non_alt_forward=9,
        non_alt_reverse=9,
        usable=20,
        unusable=0,
        unusable_by_reason={},
    )

    weaker_case = AggregatedEvidence(
        alt_forward=3,
        alt_reverse=0,
        non_alt_forward=7,
        non_alt_reverse=0,
        usable=10,
        unusable=0,
        unusable_by_reason={},
    )
    stronger_case = AggregatedEvidence(
        alt_forward=8,
        alt_reverse=0,
        non_alt_forward=2,
        non_alt_reverse=0,
        usable=10,
        unusable=0,
        unusable_by_reason={},
    )

    weaker_stats = compute_stats(weaker_case, normal_evidence)
    stronger_stats = compute_stats(stronger_case, normal_evidence)

    assert stronger_stats.bayes_factor < weaker_stats.bayes_factor
    assert stronger_stats.artifact_posterior < weaker_stats.artifact_posterior
