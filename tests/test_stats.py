from skua.evidence import AggregatedEvidence
from skua import compute_stats
from skua.stats import Stats, aggregate_evidence, estimate_rho, truncated_normal_evidences
import sys
import pytest


@pytest.mark.parametrize("pseudocount", [float("nan"), float("inf"), -float("inf"), 0.0, -1.0])
@pytest.mark.parametrize("depth", [0, 200])
def test_compute_stats_rejects_invalid_pseudocount(pseudocount, depth) -> None:
    case = _make_normal(10 if depth else 0, 10 if depth else 0, depth)
    normal = _make_normal(0, 0, depth)
    with pytest.raises(ValueError, match="pseudocount must be finite and > 0"):
        compute_stats(case, normal, pseudocount=pseudocount)


@pytest.mark.parametrize(
    ("parameters", "log_bayes_factor", "posterior"),
    [
        ({}, -14.294928638148122, 6.191431955750486e-7),
        ({"pseudocount": 0.5}, -13.327051319298334, 1.629800490504908e-6),
        ({"rho": 0.02, "pseudocount": 1, "truncate": 1, "prior_artifact_probability": 0.2},
         -6.393538514971425, 0.00041790732073216204),
    ],
)
def test_valid_parameters_preserve_numerical_results(parameters, log_bayes_factor, posterior) -> None:
    # Values captured before parameter validation changed (base eca1d8c).
    stats = compute_stats(_make_normal(10, 10, 200), _make_normal(0, 0, 200), **parameters)
    assert stats.log_bayes_factor_artifact_vs_variant == pytest.approx(log_bayes_factor)
    assert stats.artifact_posterior == pytest.approx(posterior)


def test_statistical_parameter_boundaries_accept_equal_interior_limits() -> None:
    empty = _make_normal(0, 0, 0)
    stats = compute_stats(empty, empty, truncate=1, mu_min=0.5, mu_max=0.5)
    assert stats.artifact_posterior == 0.5
    assert estimate_rho([], rho_min=0.01, rho_max=0.01, truncate=1) == 0.01


@pytest.mark.parametrize("parameter", ["rho", "mu_min", "mu_max", "prior_artifact_probability", "truncate"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf"), 0.0, -0.1, 1.1])
@pytest.mark.parametrize("normals", [None, []])
def test_compute_stats_rejects_invalid_probability_parameters(parameter, value, normals) -> None:
    empty = _make_normal(0, 0, 0)
    with pytest.raises(ValueError, match=parameter):
        compute_stats(empty, empty, per_sample_evidences=normals, **{parameter: value})


@pytest.mark.parametrize(
    "parameters",
    [{"rho": 1}, {"mu_min": 1}, {"mu_max": 1}, {"mu_min": 0.8, "mu_max": 0.2},
     {"prior_artifact_probability": 1}],
)
def test_compute_stats_rejects_invalid_probability_bounds(parameters) -> None:
    empty = _make_normal(0, 0, 0)
    with pytest.raises(ValueError, match="rho|mu_min|mu_max|prior_artifact_probability"):
        compute_stats(empty, empty, **parameters)


@pytest.mark.parametrize(
    "parameters",
    [
        {"pseudo": float("nan")}, {"pseudo": float("inf")}, {"pseudo": 0},
        {"truncate": float("nan")}, {"truncate": 0}, {"truncate": 1.1},
        {"rho_min": float("nan")}, {"rho_max": float("inf")},
        {"rho_min": 0}, {"rho_max": 1}, {"rho_min": 0.2, "rho_max": 0.1},
    ],
)
def test_estimate_rho_validates_parameters_even_without_samples(parameters) -> None:
    with pytest.raises(ValueError, match="pseudo|truncate|rho_min|rho_max"):
        estimate_rho([], **parameters)


@pytest.mark.parametrize(
    "parameters",
    [{"epsilon": float("nan")}, {"epsilon": float("inf")}, {"epsilon": 0},
     {"truncate": float("nan")}, {"truncate": 0}, {"truncate": 1.1}],
)
def test_truncated_normals_validates_parameters_even_without_samples(parameters) -> None:
    with pytest.raises(ValueError, match="epsilon|truncate"):
        truncated_normal_evidences([], **parameters)


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
    assert isinstance(stats.log_bayes_factor_artifact_vs_variant, float)
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

    assert stronger_stats.log_bayes_factor_artifact_vs_variant < weaker_stats.log_bayes_factor_artifact_vs_variant
    assert stronger_stats.artifact_posterior < weaker_stats.artifact_posterior


def test_compute_stats_artifact_prior_changes_only_the_posterior() -> None:
    case_evidence = _make_normal(4, 0, 10)
    normal_evidence = _make_normal(1, 0, 100)

    low_prior = compute_stats(
        case_evidence,
        normal_evidence,
        prior_artifact_probability=0.2,
    )
    high_prior = compute_stats(
        case_evidence,
        normal_evidence,
        prior_artifact_probability=0.8,
    )

    assert low_prior.log_bayes_factor_artifact_vs_variant == (
        high_prior.log_bayes_factor_artifact_vs_variant
    )
    assert low_prior.artifact_posterior < high_prior.artifact_posterior


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


def test_estimate_rho_empty_returns_rho_min() -> None:
    assert estimate_rho([]) == 1e-4


def test_estimate_rho_single_sample_returns_rho_min() -> None:
    assert estimate_rho([_make_normal(1, 1, 100)]) == 1e-4


def test_estimate_rho_uniform_low_background_returns_rho_min() -> None:
    # All samples have the same very low error rate -> no overdispersion -> rho_min
    samples = [_make_normal(1, 1, 1000) for _ in range(10)]
    rho = estimate_rho(samples)
    assert rho == 1e-4


@pytest.mark.parametrize("variable_background", [False, True])
@pytest.mark.parametrize(
    "excluded_counts",
    [
        [(50, 50, 1000)] * 5,  # Exactly at the truncation threshold.
        [(500, 500, 10000)] * 5,  # Same ALT fraction, greater depth.
        [(100, 50, 1000), (2000, 1000, 10000)],  # Above the threshold.
    ],
)
def test_excluded_normals_do_not_change_dispersion_or_scores(
    variable_background, excluded_counts,
) -> None:
    retained = [_make_normal(10, 10, 1000) for _ in range(5)]
    if variable_background:
        retained += [_make_normal(40, 40, 1000) for _ in range(5)]
    case = _make_normal(10, 10, 200)
    baseline = compute_stats(
        case, aggregate_evidence(retained), per_sample_evidences=retained,
    )
    normals = retained + [_make_normal(*counts) for counts in excluded_counts]
    stats = compute_stats(
        case, aggregate_evidence(normals), per_sample_evidences=normals,
    )

    assert truncated_normal_evidences(normals) == retained
    assert estimate_rho(normals) == pytest.approx(estimate_rho(retained))
    assert stats.dispersion_rho == pytest.approx(baseline.dispersion_rho)
    assert stats.normal_counts == baseline.normal_counts
    assert stats.log_bayes_factor_artifact_vs_variant == pytest.approx(
        baseline.log_bayes_factor_artifact_vs_variant,
    )
    assert stats.artifact_posterior == pytest.approx(baseline.artifact_posterior)


def test_estimate_rho_overdispersed_samples_returns_higher_rho() -> None:
    # Both the 2% and 8% groups are retained and contribute real variation.
    low = [_make_normal(10, 10, 1000) for _ in range(5)]
    high = [_make_normal(40, 40, 1000) for _ in range(5)]
    rho = estimate_rho(low + high)
    assert rho > 1e-4


def test_estimate_rho_result_is_within_bounds() -> None:
    samples = [_make_normal(i, i, 500) for i in range(1, 21)]
    rho = estimate_rho(samples)
    assert 1e-4 <= rho <= 0.1


@pytest.mark.parametrize(
    ("counts", "rho", "log_bayes_factor", "posterior"),
    [
        ([(10, 10, 1000)] * 5,
         1e-4, -15.199404710030649, 2.5060071016458127e-7),
        ([(10, 10, 1000)] * 5 + [(40, 40, 1000)] * 5,
         0.020072704283230606, -1.5856301824479715, 0.16999958817075578),
        ([(1, 0, 120), (2, 1, 140), (3, 0, 160), (4, 1, 180)],
         1e-4, -12.199449784340686, 5.033198870934455e-6),
    ],
)
def test_scores_are_preserved_when_all_normals_are_retained(
    counts, rho, log_bayes_factor, posterior,
) -> None:
    # Baselines captured at v0.7.4 (0ab3cef), with no excluded normals.
    samples = [_make_normal(*values) for values in counts]
    stats = compute_stats(
        _make_normal(10, 10, 200), aggregate_evidence(samples),
        per_sample_evidences=samples,
    )

    assert truncated_normal_evidences(samples) == samples
    assert estimate_rho(samples) == pytest.approx(rho)
    assert stats.dispersion_rho == pytest.approx(rho)
    assert stats.log_bayes_factor_artifact_vs_variant == pytest.approx(log_bayes_factor)
    assert stats.artifact_posterior == pytest.approx(posterior)


def test_estimate_rho_preserves_pseudocount_with_excluded_normals() -> None:
    retained = [_make_normal(10, 10, 1000) for _ in range(5)]
    retained += [_make_normal(40, 40, 1000) for _ in range(5)]
    excluded = [_make_normal(50, 50, 1000), _make_normal(2000, 1000, 10000)]
    # v0.7.4 baseline with pseudo=0.5 and only the retained panel.
    expected = 0.02001653472902874
    assert estimate_rho(retained, pseudo=0.5) == pytest.approx(expected)
    assert estimate_rho(retained + excluded, pseudo=0.5) == pytest.approx(expected)


@pytest.mark.parametrize("retained_count", [0, 1])
@pytest.mark.parametrize("truncate", [0.1, 0.2])
def test_estimate_rho_falls_back_with_fewer_than_two_retained_normals(
    retained_count, truncate,
) -> None:
    retained = [_make_normal(10, 10, 1000) for _ in range(retained_count)]
    excluded = [_make_normal(100, 100, 1000), _make_normal(2000, 1000, 10000)]
    assert estimate_rho(retained + excluded, truncate=truncate, rho_min=0.002) == 0.002


@pytest.mark.parametrize(("rho_min", "rho_max", "expected"), [(0.03, 0.1, 0.03), (1e-4, 0.01, 0.01)])
def test_estimate_rho_clips_retained_dispersion_to_configured_bounds(
    rho_min, rho_max, expected,
) -> None:
    retained = [_make_normal(10, 10, 1000) for _ in range(5)]
    retained += [_make_normal(40, 40, 1000) for _ in range(5)]
    excluded = [_make_normal(500, 500, 10000)]
    assert estimate_rho(retained + excluded, rho_min=rho_min, rho_max=rho_max) == expected


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

    assert truncated.log_bayes_factor_artifact_vs_variant < untruncated.log_bayes_factor_artifact_vs_variant
    assert truncated.artifact_posterior < untruncated.artifact_posterior
