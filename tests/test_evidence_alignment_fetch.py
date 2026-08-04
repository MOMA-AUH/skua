from skua.evidence import (
    UnusableReason,
    collect_evidence_from_alignment,
    collect_evidence_from_alignment_batch,
)
from tests.helpers import FakeAlignmentFile, FakeRead, build_linear_pairs
from skua.variants import Variant



def test_collect_evidence_from_alignment_fetches_one_locus_window() -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        )
    ]
    alignment_file = FakeAlignmentFile(reads)

    counts = collect_evidence_from_alignment(
        alignment_file,
        contig="chr7",
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
    )

    assert alignment_file.fetch_calls == [("chr7", 105, 106)]
    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 0



def test_collect_evidence_from_alignment_propagates_mixed_counts() -> None:
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
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=5,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    counts = collect_evidence_from_alignment(
        alignment_file,
        contig="chr1",
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 1
    assert counts.usable == 2
    assert counts.unusable == 1
    assert counts.unusable_by_reason[UnusableReason.LOW_MAPQ] == 1


def test_collect_evidence_from_alignment_limits_reads_to_allowed_read_groups() -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "case-rg"},
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "other-rg"},
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    counts = collect_evidence_from_alignment(
        alignment_file,
        contig="chr1",
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        allowed_read_group_ids=frozenset({"case-rg"}),
    )

    assert counts.alt_forward == 1
    assert counts.non_alt_forward == 0
    assert counts.usable == 1


def test_batch_collection_fetches_dense_variants_once_with_site_parity() -> None:
    variants = (
        Variant(contig="chr1", ref_pos0=105, ref="A", alt="T"),
        Variant(contig="chr1", ref_pos0=108, ref="A", alt="C"),
    )
    reads = [
        FakeRead(
            query_name="alt-both",
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAACA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            query_name="ref-both",
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    batch_alignment = FakeAlignmentFile(reads)
    site_alignment = FakeAlignmentFile(reads)

    batch_evidences = collect_evidence_from_alignment_batch(
        batch_alignment,
        variants,
    )
    site_evidences = tuple(
        collect_evidence_from_alignment(
            site_alignment,
            contig=variant.contig,
            ref_pos0=variant.ref_pos0,
            ref_base=variant.ref,
            alt_base=variant.alt,
        )
        for variant in variants
    )

    assert batch_alignment.fetch_calls == [("chr1", 105, 109)]
    assert batch_evidences == site_evidences


def test_batch_collection_preserves_fragment_and_read_group_semantics() -> None:
    variants = (
        Variant(contig="chr1", ref_pos0=105, ref="A", alt="T"),
        Variant(contig="chr1", ref_pos0=108, ref="A", alt="C"),
    )
    reads = [
        FakeRead(
            query_name="conflicting-pair",
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAACA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            flag=0x43,
            tags={"RG": "case-rg"},
        ),
        FakeRead(
            query_name="conflicting-pair",
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            flag=0x83,
            tags={"RG": "case-rg"},
        ),
        FakeRead(
            query_name="excluded",
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAACA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "other-rg"},
        ),
    ]
    batch_alignment = FakeAlignmentFile(reads)
    site_alignment = FakeAlignmentFile(reads)

    batch_evidences = collect_evidence_from_alignment_batch(
        batch_alignment,
        variants,
        allowed_read_group_ids=frozenset({"case-rg"}),
    )
    site_evidences = tuple(
        collect_evidence_from_alignment(
            site_alignment,
            contig=variant.contig,
            ref_pos0=variant.ref_pos0,
            ref_base=variant.ref,
            alt_base=variant.alt,
            allowed_read_group_ids=frozenset({"case-rg"}),
        )
        for variant in variants
    )

    assert batch_evidences == site_evidences
    assert all(evidence.unusable == 1 for evidence in batch_evidences)
    assert all(
        evidence.unusable_by_reason == {UnusableReason.CONFLICTING_MATES: 1}
        for evidence in batch_evidences
    )


def test_batch_collection_matches_site_collection_for_insertion_and_snv() -> None:
    variants = (
        Variant(contig="chr1", ref_pos0=100, ref="A", alt="AT"),
        Variant(contig="chr1", ref_pos0=105, ref="A", alt="T"),
    )
    reads = [
        FakeRead(
            query_name="alt",
            mapping_quality=60,
            is_reverse=False,
            query_sequence="ATAAAATAAA",
            query_qualities=[35] * 10,
            aligned_pairs=[
                (0, 100),
                (1, None),
                (2, 101),
                (3, 102),
                (4, 103),
                (5, 104),
                (6, 105),
                (7, 106),
                (8, 107),
                (9, 108),
            ],
        ),
        FakeRead(
            query_name="ref",
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    batch_alignment = FakeAlignmentFile(reads)
    site_alignment = FakeAlignmentFile(reads)

    batch_evidences = collect_evidence_from_alignment_batch(
        batch_alignment,
        variants,
    )
    site_evidences = tuple(
        collect_evidence_from_alignment(
            site_alignment,
            contig=variant.contig,
            ref_pos0=variant.ref_pos0,
            ref_base=variant.ref,
            alt_base=variant.alt,
        )
        for variant in variants
    )

    assert batch_evidences == site_evidences
