import gc
import tracemalloc

import pysam
import pytest

from skua import annotate_variant, annotate_variants
from skua.evidence import (
    AggregatedEvidence,
    AlleleSupport,
    UnusableReason,
    classify_variant_read,
    collect_evidence_from_alignment_batch,
)
from skua.variants import Variant
from tests.helpers import FakeAlignmentFile
from tests.test_evidence_pysam_integration import build_aligned_segment


@pytest.mark.parametrize("operation", [2, 3], ids=["deletion", "skip"])
@pytest.mark.parametrize("mode", ["snv", "deletion", "batch"])
def test_irrelevant_cigar_gap_length_does_not_increase_allocation(operation, mode):
    peaks = []
    for length in (100, 200_000):
        read = build_aligned_segment(
            query_name="gap", query_sequence="TG", reference_start=100,
            cigar=((0, 1), (operation, length), (0, 1)),
        )
        variants = [
            Variant("chr1", 100, "A", "T"),
            Variant("chr1", 101 + length, "A", "G"),
        ]
        alignment = FakeAlignmentFile([read])
        gc.collect()
        tracemalloc.start()
        try:
            if mode == "batch":
                calls = collect_evidence_from_alignment_batch(alignment, variants)
            else:
                call = classify_variant_read(
                    read, ref_pos0=100,
                    ref_base="TC" if mode == "deletion" else "A",
                    alt_base="T",
                )
            peaks.append(tracemalloc.get_traced_memory()[1])
        finally:
            tracemalloc.stop()

        if mode == "batch":
            assert [evidence.alt_forward for evidence in calls] == [1, 1]
        else:
            expected = (
                AlleleSupport.ALT if mode == "snv" else
                AlleleSupport.NON_ALT if operation == 2 else AlleleSupport.UNUSABLE
            )
            assert call.support == expected

    # Generous budgets for a two-base read; catch expansion in either the parser
    # or deletion classification without making CI depend on wall-clock timing.
    assert max(peaks) < 256 * 1024, peaks
    assert peaks[1] - peaks[0] < 64 * 1024, peaks


@pytest.mark.parametrize("operation", [2, 3], ids=["deletion", "skip"])
def test_long_gap_boundaries_and_indels_agree_in_singleton_and_batch_bam_apis(
    tmp_path, operation,
):
    length = 200_000
    path = tmp_path / "gap.bam"
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1_000_000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for reverse in (False, True):
            bam.write(build_aligned_segment(
                query_name=f"gap_{reverse}", query_sequence="TG", reference_start=100,
                cigar=((0, 1), (operation, length), (0, 1)), is_reverse=reverse,
            ))
    pysam.index(str(path))

    variants = [
        Variant("chr1", 100, "A", "T"),
        Variant("chr1", 101, "A", "T"),          # First gap position.
        Variant("chr1", 100 + length // 2, "A", "T"),
        Variant("chr1", 100 + length, "A", "T"),  # Last gap position.
        Variant("chr1", 101 + length, "A", "G"),
        Variant("chr1", 100, "AA", "TG"),        # MNV crossing the gap.
        Variant("chr1", 100, "T", "TA"),         # Insertion next to a gap.
        Variant("chr1", 100, "T" + "A" * length, "T"),
        Variant("chr1", 100, "TA", "T"),         # Shorter deletion.
    ]
    alt = AggregatedEvidence(1, 1, 0, 0, 2, 0, {})
    non_alt = AggregatedEvidence(0, 0, 1, 1, 2, 0, {})
    unusable = AggregatedEvidence(0, 0, 0, 0, 0, 2, {UnusableReason.NO_BASE_AT_SITE: 2})
    expected = [
        alt, unusable, unusable, unusable, alt, unusable, unusable,
        alt if operation == 2 else unusable,
        non_alt if operation == 2 else unusable,
    ]
    with pysam.AlignmentFile(path, "rb") as alignment:
        assert [annotate_variant(alignment, variant) for variant in variants] == expected
        assert list(annotate_variants(alignment, variants)) == list(zip(variants, expected))


@pytest.mark.parametrize(
    ("cigar", "sequence", "support"),
    [
        ("1M199999D1D1M", "TG", AlleleSupport.ALT),
        ("1M200000D", "T", AlleleSupport.UNUSABLE),
        ("1M200000D1I1M", "TCG", AlleleSupport.UNUSABLE),
        ("1M200000D1N1M", "TG", AlleleSupport.UNUSABLE),
        ("1M1N199999D1M", "TG", AlleleSupport.UNUSABLE),
        ("1M200000D1S1M", "TCG", AlleleSupport.UNUSABLE),
        ("1M200000D1H1M", "TG", AlleleSupport.UNUSABLE),
        ("1M200000D1P1M", "TG", AlleleSupport.UNUSABLE),
        ("1S1M200000D1M1S", "CTGC", AlleleSupport.ALT),
    ],
)
def test_compact_deletion_preserves_flank_and_complex_event_rules(cigar, sequence, support):
    read = build_aligned_segment(
        query_name="deletion", query_sequence=sequence, reference_start=100,
    )
    read.cigarstring = cigar
    call = classify_variant_read(read, ref_pos0=100, ref_base="T" + "A" * 200_000, alt_base="T")
    assert call.support == support
    assert call.reason == (UnusableReason.NO_BASE_AT_SITE if support == AlleleSupport.UNUSABLE else None)


@pytest.mark.parametrize("operation", [2, 3], ids=["deletion", "skip"])
@pytest.mark.parametrize("indel", ["insertion", "deletion"])
@pytest.mark.parametrize("quality", [40, 5, 255])
def test_indel_after_long_gap_preserves_anchor_and_quality_rules(operation, indel, quality):
    length = 200_000
    suffix = ((0, 1), (1, 1), (0, 1)) if indel == "insertion" else ((0, 1), (2, 1), (0, 1))
    read = build_aligned_segment(
        query_name="downstream", query_sequence="ATCG" if indel == "insertion" else "ATG",
        reference_start=100, cigar=((0, 1), (operation, length), *suffix),
    )
    qualities = list(read.query_qualities)
    qualities[1] = quality
    read.query_qualities = qualities
    call = classify_variant_read(
        read, ref_pos0=101 + length,
        ref_base="T" if indel == "insertion" else "TC",
        alt_base="TC" if indel == "insertion" else "T",
    )
    assert call.support == (AlleleSupport.ALT if quality == 40 else AlleleSupport.UNUSABLE)
    assert call.reason == {
        40: None, 5: UnusableReason.LOW_BASEQ, 255: UnusableReason.MISSING_BASEQ,
    }[quality]
