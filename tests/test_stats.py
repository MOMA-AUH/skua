from skua.evidence import AggregatedEvidence
from skua.stats import StrandAwarePonStats, compute_strand_aware_pon_stats


def test_compute_strand_aware_pon_stats_returns_typed_background_and_score() -> None:
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

    stats = compute_strand_aware_pon_stats(case_evidence, normal_evidence)

    assert isinstance(stats, StrandAwarePonStats)
    assert stats.case_counts["alt_forward"] == 8
    assert stats.normal_counts["non_alt_forward"] == 9
    assert stats.background_rate_by_channel["alt_forward"] == 0.05
    assert stats.expected_case_counts["alt_forward"] == 0.5
    assert stats.channel_score_contribution["alt_forward"] > 0.0
    assert stats.combined_score > 0.0
    assert 0.0 <= stats.p_value <= 1.0
    assert stats.method == "chi_square_4channel_approx"
    assert stats.degrees_of_freedom == 4
    assert stats.pseudocount == 0.5
    assert stats.min_expected_count == 0.5
    assert stats.approximation_warning is True


def test_compute_strand_aware_pon_stats_is_stable_for_zero_depth() -> None:
    zero_evidence = AggregatedEvidence(
        alt_forward=0,
        alt_reverse=0,
        non_alt_forward=0,
        non_alt_reverse=0,
        usable=0,
        unusable=0,
        unusable_by_reason={},
    )

    stats = compute_strand_aware_pon_stats(zero_evidence, zero_evidence)

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
    assert stats.channel_score_contribution == {
        "alt_forward": 0.0,
        "alt_reverse": 0.0,
        "non_alt_forward": 0.0,
        "non_alt_reverse": 0.0,
    }
    assert stats.combined_score == 0.0
    assert stats.p_value == 1.0
    assert stats.method == "chi_square_4channel_approx"
    assert stats.degrees_of_freedom == 4
    assert stats.pseudocount == 0.5
    assert stats.min_expected_count == 0.0
    assert stats.approximation_warning is True


def test_compute_strand_aware_pon_stats_p_value_decreases_with_stronger_signal() -> None:
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

    weaker_stats = compute_strand_aware_pon_stats(weaker_case, normal_evidence)
    stronger_stats = compute_strand_aware_pon_stats(stronger_case, normal_evidence)

    assert stronger_stats.combined_score > weaker_stats.combined_score
    assert stronger_stats.p_value < weaker_stats.p_value


def test_compute_strand_aware_pon_stats_can_disable_warning_with_large_expected_counts() -> None:
    case_evidence = AggregatedEvidence(
        alt_forward=50,
        alt_reverse=50,
        non_alt_forward=50,
        non_alt_reverse=50,
        usable=200,
        unusable=0,
        unusable_by_reason={},
    )
    normal_evidence = AggregatedEvidence(
        alt_forward=250,
        alt_reverse=250,
        non_alt_forward=250,
        non_alt_reverse=250,
        usable=1000,
        unusable=0,
        unusable_by_reason={},
    )

    stats = compute_strand_aware_pon_stats(case_evidence, normal_evidence)

    assert stats.min_expected_count == 50.0
    assert stats.approximation_warning is False