from skua.evidence import AggregatedEvidence
from skua.stats import Stats, compute_stats, estimate_rho
import sys


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
    assert stats.pseudocount == sys.float_info.epsilon


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
    assert stats.pseudocount == sys.float_info.epsilon


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


def _make_normal(alt_fw: int, alt_bw: int, depth: int) -> AggregatedEvidence:
    non = depth - alt_fw - alt_bw
    return AggregatedEvidence(
        alt_forward=alt_fw,
        alt_reverse=alt_bw,
        non_alt_forward=non // 2,
        non_alt_reverse=non - non // 2,
        usable=depth,
        unusable=0,
        unusable_by_reason={},
    )


def _estimate_rho_reference(
    per_sample_evidences: list[AggregatedEvidence],
    *,
    truncate: float = 0.1,
    rho_min: float = 1e-4,
    rho_max: float = 0.1,
    pseudo: float = sys.float_info.epsilon,
) -> float:
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
        if not valid_depths:
            rho_by_channel.append(rho_min)
            continue

        s2 = (
            included_count
            * sum(
                valid_depths[value_index] * (valid_mu[value_index] - nu) ** 2
                for value_index in range(len(valid_depths))
            )
            / ((included_count - 1) * sum(valid_depths))
        )
        sum_inv_nix = sum(1.0 / depth for depth in valid_depths)
        denom = included_count - sum_inv_nix
        if denom <= 0 or nu <= 0.0 or nu >= 1.0:
            rho_by_channel.append(rho_min)
            continue

        rho_hat = (included_count * (s2 / nu / (1.0 - nu)) - sum_inv_nix) / denom
        if not sys.float_info.max > abs(rho_hat) or rho_hat != rho_hat:
            rho_by_channel.append(rho_min)
            continue

        rho_by_channel.append(min(max(min(max(rho_hat, 0.0), 1.0), rho_min), rho_max))

    return rho_by_channel[0]


def test_estimate_rho_empty_returns_rho_min() -> None:
    assert estimate_rho([]) == 1e-4


def test_estimate_rho_single_sample_returns_rho_min() -> None:
    assert estimate_rho([_make_normal(1, 1, 100)]) == 1e-4


def test_estimate_rho_uniform_low_background_returns_rho_min() -> None:
    # All samples have the same very low error rate -> no overdispersion -> rho_min
    samples = [_make_normal(1, 1, 1000) for _ in range(10)]
    rho = estimate_rho(samples)
    assert rho == 1e-4


def test_estimate_rho_overdispersed_samples_returns_higher_rho() -> None:
    # Half samples with ~1% error, half with ~5% error -> measurable overdispersion
    low = [_make_normal(10, 10, 1000) for _ in range(5)]   # ~2%
    high = [_make_normal(50, 50, 1000) for _ in range(5)]  # ~10% -> truncated
    rho = estimate_rho(low + high)
    assert rho > 1e-4


def test_estimate_rho_result_is_within_bounds() -> None:
    samples = [_make_normal(i, i, 500) for i in range(1, 21)]
    rho = estimate_rho(samples)
    assert 1e-4 <= rho <= 0.1


def test_estimate_rho_matches_reference_two_channel_formula() -> None:
    samples = [
        _make_normal(1, 0, 120),
        _make_normal(2, 1, 140),
        _make_normal(3, 0, 160),
        _make_normal(4, 1, 180),
        _make_normal(12, 8, 150),
    ]

    assert estimate_rho(samples) == _estimate_rho_reference(samples)


def test_compute_stats_uses_estimated_rho_from_per_sample_evidences() -> None:
    # When per_sample_evidences is supplied, dispersion_rho should differ from default
    case_evidence = _make_normal(8, 0, 10)
    # Build a set of moderately overdispersed normals
    per_sample = [_make_normal(i % 3, (i + 1) % 3, 200) for i in range(20)]
    normal_aggregate = AggregatedEvidence(
        alt_forward=sum(s.alt_forward for s in per_sample),
        alt_reverse=sum(s.alt_reverse for s in per_sample),
        non_alt_forward=sum(s.non_alt_forward for s in per_sample),
        non_alt_reverse=sum(s.non_alt_reverse for s in per_sample),
        usable=sum(s.usable for s in per_sample),
        unusable=0,
        unusable_by_reason={},
    )
    stats_fixed = compute_stats(case_evidence, normal_aggregate)
    stats_estimated = compute_stats(
        case_evidence, normal_aggregate, per_sample_evidences=per_sample
    )
    # When per_sample_evidences is provided, rho is estimated (may equal rho_min,
    # but dispersion_rho reflects the actual estimated value not the kwarg default)
    assert stats_estimated.dispersion_rho == estimate_rho(per_sample)
    assert stats_fixed.dispersion_rho == 1e-4


def test_compute_stats_applies_truncation_to_background_pool() -> None:
    case_evidence = _make_normal(4, 0, 10)

    low_background = [_make_normal(1, 0, 200) for _ in range(10)]
    high_background_outlier = _make_normal(40, 10, 100)
    per_sample = low_background + [high_background_outlier]

    normal_aggregate = AggregatedEvidence(
        alt_forward=sum(s.alt_forward for s in per_sample),
        alt_reverse=sum(s.alt_reverse for s in per_sample),
        non_alt_forward=sum(s.non_alt_forward for s in per_sample),
        non_alt_reverse=sum(s.non_alt_reverse for s in per_sample),
        usable=sum(s.usable for s in per_sample),
        unusable=0,
        unusable_by_reason={},
    )

    untruncated = compute_stats(
        case_evidence,
        normal_aggregate,
        per_sample_evidences=per_sample,
        truncate=1.0,
    )
    truncated = compute_stats(
        case_evidence,
        normal_aggregate,
        per_sample_evidences=per_sample,
        truncate=0.1,
    )

    assert truncated.bayes_factor < untruncated.bayes_factor
    assert truncated.artifact_posterior < untruncated.artifact_posterior
