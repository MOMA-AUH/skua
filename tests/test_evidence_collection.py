from skua.evidence import UnusableReason, collect_snv_evidence
from tests.helpers import FakeRead, build_linear_pairs



def test_collect_snv_evidence_combines_classification_and_aggregation() -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]

    counts = collect_snv_evidence(
        reads,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert counts.alt_forward == 1
    assert counts.alt_reverse == 1
    assert counts.non_alt_forward == 1
    assert counts.non_alt_reverse == 0
    assert counts.usable == 3
    assert counts.unusable == 0



def test_collect_snv_evidence_propagates_unusable_reason_counts() -> None:
    reads = [
        FakeRead(
            mapping_quality=5,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAANAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35, 35, 35, 35, 35, 10, 35, 35, 35, 35],
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]

    counts = collect_snv_evidence(
        reads,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert counts.alt_forward == 0
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 0
    assert counts.usable == 0
    assert counts.unusable == 3
    assert counts.unusable_by_reason[UnusableReason.LOW_MAPQ] == 1
    assert counts.unusable_by_reason[UnusableReason.INVALID_BASE] == 1
    assert counts.unusable_by_reason[UnusableReason.LOW_BASEQ] == 1



def test_collect_snv_evidence_handles_empty_reads() -> None:
    counts = collect_snv_evidence(
        [],
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
    )

    assert counts.alt_forward == 0
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 0
    assert counts.usable == 0
    assert counts.unusable == 0
    assert counts.unusable_by_reason == {}
