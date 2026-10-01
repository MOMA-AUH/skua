from pathlib import Path

import pysam
import pytest

from skua import (
    annotate_variant,
    annotate_variants,
    annotate_vcf_with_normals,
    annotate_vcf_with_pon,
    build_pon,
    read_pon_evidence,
)
from skua.evidence import (
    AggregatedEvidence,
    UnusableReason,
    collect_evidence_from_alignment,
    collect_evidence_from_alignment_batch,
)
from skua.variants import Variant


HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000}],
    "RG": [
        {"ID": "rg1", "SM": "sample"},
        {"ID": "rg2", "SM": "sample"},
    ],
}

def build_aligned_segment(
    *,
    query_name: str,
    query_sequence: str,
    reference_start: int,
    mapping_quality: int = 60,
    is_reverse: bool = False,
    flag: int | None = None,
    cigar: tuple[tuple[int, int], ...] | None = None,
    read_group: str | None = None,
) -> pysam.AlignedSegment:
    segment = pysam.AlignedSegment()
    segment.query_name = query_name
    segment.query_sequence = query_sequence
    segment.flag = flag if flag is not None else (147 if is_reverse else 99)
    segment.reference_id = 0
    segment.reference_start = reference_start
    segment.mapping_quality = mapping_quality
    segment.cigar = cigar if cigar is not None else ((0, len(query_sequence)),)
    segment.next_reference_id = 0
    segment.next_reference_start = reference_start
    segment.template_length = -len(query_sequence) if is_reverse else len(query_sequence)
    segment.query_qualities = pysam.qualitystring_to_array("I" * len(query_sequence))
    if read_group is not None:
        segment.set_tag("RG", read_group)
    return segment



def create_test_bam(tmp_path: Path, reads: list[pysam.AlignedSegment]) -> Path:
    unsorted_bam = tmp_path / "reads.unsorted.bam"
    sorted_bam = tmp_path / "reads.bam"

    with pysam.AlignmentFile(unsorted_bam, "wb", header=HEADER) as bam_file:
        for read in reads:
            bam_file.write(read)

    pysam.sort("-o", str(sorted_bam), str(unsorted_bam))
    pysam.index(str(sorted_bam))
    return sorted_bam


@pytest.mark.parametrize("is_reverse", [False, True], ids=["forward", "reverse"])
@pytest.mark.parametrize("inserted_baseq", [40, 5], ids=["high-baseq", "low-baseq"])
@pytest.mark.parametrize("batch", [False, True], ids=["singleton", "batch"])
def test_mnv_internal_insertion_is_unusable_in_singleton_and_batch_apis(
    tmp_path: Path, is_reverse: bool, inserted_baseq: int, batch: bool,
) -> None:
    reads = []
    for index in range(5):
        read = build_aligned_segment(
            query_name=f"internal_insertion_{index}",
            query_sequence="TGC",
            reference_start=100,
            is_reverse=is_reverse,
            cigar=((0, 1), (1, 1), (0, 1)),  # 1M1I1M.
        )
        read.query_qualities = [40, inserted_baseq, 40]
        reads.append(read)
    bam_path = create_test_bam(tmp_path, reads)
    variants = [
        Variant(contig="chr1", ref_pos0=100, ref="AA", alt="TC"),
        Variant(contig="chr1", ref_pos0=100, ref="AA", alt="GG"),
    ]
    expected = AggregatedEvidence(
        alt_forward=0,
        alt_reverse=0,
        non_alt_forward=0,
        non_alt_reverse=0,
        usable=0,
        unusable=5,
        unusable_by_reason={UnusableReason.NO_BASE_AT_SITE: 5},
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment:
        if batch:
            assert list(annotate_variants(alignment, variants)) == [
                (variants[0], expected), (variants[1], expected),
            ]
        else:
            for variant in variants:
                assert annotate_variant(alignment, variant) == expected


@pytest.mark.parametrize(
    ("sequence", "cigar"),
    [
        ("TC", ((0, 2),)),  # 2M.
        ("TC", ((8, 1), (8, 1))),  # Consecutive X operations are aligned.
        ("GTC", ((1, 1), (0, 2))),  # Insertion before the MNV interval.
        ("TCG", ((0, 2), (1, 1))),  # Insertion after the MNV interval.
        ("GTCG", ((4, 1), (0, 2), (4, 1))),  # Flanking soft clips.
    ],
)
def test_contiguous_mnv_retains_strand_counts_in_singleton_and_batch_apis(
    tmp_path: Path, sequence: str, cigar: tuple[tuple[int, int], ...],
) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name=f"contiguous_{reverse}",
                query_sequence=sequence,
                reference_start=100,
                is_reverse=reverse,
                cigar=cigar,
            )
            for reverse in (False, True)
        ],
    )
    variants = [
        Variant(contig="chr1", ref_pos0=100, ref="AA", alt="TC"),
        Variant(contig="chr1", ref_pos0=100, ref="AA", alt="GG"),
    ]
    expected = [
        AggregatedEvidence(1, 1, 0, 0, 2, 0, {}),
        AggregatedEvidence(0, 0, 1, 1, 2, 0, {}),
    ]
    with pysam.AlignmentFile(bam_path, "rb") as alignment:
        assert [annotate_variant(alignment, variant) for variant in variants] == expected
        assert list(annotate_variants(alignment, variants)) == list(zip(variants, expected))


def test_mnv_direct_normals_and_cached_pon_agree_with_internal_insertions(
    tmp_path: Path,
) -> None:
    bam_paths = []
    for sample, sequence, count in (("case", "TC", 2), ("normal", "AA", 20)):
        sample_path = tmp_path / sample
        sample_path.mkdir()
        reads = [
            build_aligned_segment(
                query_name=f"{sample}_aligned_{index}",
                query_sequence=sequence,
                reference_start=100,
                is_reverse=bool(index % 2),
                read_group="rg1",
            )
            for index in range(count)
        ]
        for reverse in (False, True):
            for inserted_baseq in (40, 5):
                read = build_aligned_segment(
                    query_name=f"{sample}_complex_{reverse}_{inserted_baseq}",
                    query_sequence="TGC",
                    reference_start=100,
                    is_reverse=reverse,
                    cigar=((0, 1), (1, 1), (0, 1)),
                    read_group="rg1",
                )
                read.query_qualities = [40, inserted_baseq, 40]
                reads.append(read)
        bam_paths.append(create_test_bam(sample_path, reads))

    targets = tmp_path / "targets.vcf"
    targets.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t101\t.\tAA\tTC\t.\tPASS\t.\n"
        "chr1\t101\t.\tAA\tGG\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    pon_path = tmp_path / "normals.pon.bcf"
    direct_output = tmp_path / "direct.vcf"
    cached_output = tmp_path / "cached.vcf"
    with (
        pysam.AlignmentFile(bam_paths[0], "rb") as case,
        pysam.AlignmentFile(bam_paths[1], "rb") as normal,
    ):
        annotate_vcf_with_normals(
            case, targets, normal_alignments=[normal], output_path=direct_output,
        )
        build_pon(targets, normal_alignments=[normal], output_path=pon_path)
        annotate_vcf_with_pon(case, pon_path, output_path=cached_output)

    cached_evidence = list(read_pon_evidence(pon_path))
    assert len(cached_evidence) == 2
    for _variant, normal_evidences in cached_evidence:
        assert normal_evidences == (AggregatedEvidence(0, 0, 10, 10, 20, 4, {}),)

    with (
        pysam.VariantFile(str(direct_output)) as direct_vcf,
        pysam.VariantFile(str(cached_output)) as cached_vcf,
    ):
        direct_records = list(direct_vcf)
        cached_records = list(cached_vcf)
        assert len(direct_records) == len(cached_records) == 2
        for direct, cached in zip(direct_records, cached_records, strict=True):
            # Compare all evidence and model scores, then check literal counts
            # so two equally incorrect workflows cannot satisfy this regression.
            assert cached.alleles == direct.alleles
            assert dict(cached.info) == dict(direct.info)
            assert dict(cached.samples["sample"]) == dict(direct.samples["sample"])
            assert direct.info["SKUA_PON_SAMPLE_COUNT"] == 1
            assert direct.info["SKUA_PON_ALT_FWD"] == 0
            assert direct.info["SKUA_PON_ALT_REV"] == 0
            assert direct.info["SKUA_PON_NON_ALT_FWD"] == 10
            assert direct.info["SKUA_PON_NON_ALT_REV"] == 10
            assert direct.info["SKUA_PON_USABLE"] == 20
            assert direct.info["SKUA_PON_UNUSABLE"] == 4
            sample = direct.samples["sample"]
            assert sample["SKUA_USABLE"] == 2
            assert sample["SKUA_UNUSABLE"] == 4
            if direct.alts == ("TC",):
                assert sample["SKUA_ALT_FWD"] == sample["SKUA_ALT_REV"] == 1
                assert sample["SKUA_NON_ALT_FWD"] == sample["SKUA_NON_ALT_REV"] == 0
            else:
                assert sample["SKUA_ALT_FWD"] == sample["SKUA_ALT_REV"] == 0
                assert sample["SKUA_NON_ALT_FWD"] == sample["SKUA_NON_ALT_REV"] == 1



def test_collect_evidence_from_alignment_with_real_bam(tmp_path: Path) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="alt_forward",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                is_reverse=False,
            ),
            build_aligned_segment(
                query_name="alt_reverse",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                is_reverse=True,
            ),
            build_aligned_segment(
                query_name="ref_reverse",
                query_sequence="AAAAAAAAAA",
                reference_start=100,
                is_reverse=True,
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
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
    assert counts.alt_reverse == 1
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 1
    assert counts.usable == 3
    assert counts.unusable == 0


def test_batch_collection_matches_site_collection_with_real_bam(tmp_path: Path) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="alt-both",
                query_sequence="AAAAATAACA",
                reference_start=100,
            ),
            build_aligned_segment(
                query_name="ref-both",
                query_sequence="AAAAAAAAAA",
                reference_start=100,
                is_reverse=True,
            ),
        ],
    )
    variants = (
        Variant(contig="chr1", ref_pos0=105, ref="A", alt="T"),
        Variant(contig="chr1", ref_pos0=108, ref="A", alt="C"),
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        batch_evidences = collect_evidence_from_alignment_batch(
            alignment_file,
            variants,
        )
        site_evidences = tuple(
            collect_evidence_from_alignment(
                alignment_file,
                contig=variant.contig,
                ref_pos0=variant.ref_pos0,
                ref_base=variant.ref,
                alt_base=variant.alt,
            )
            for variant in variants
        )

    assert batch_evidences == site_evidences



def test_collect_evidence_from_alignment_tracks_real_bam_unusable_reads(tmp_path: Path) -> None:
    low_mapq = build_aligned_segment(
        query_name="low_mapq",
        query_sequence="AAAAATAAAA",
        reference_start=100,
        mapping_quality=5,
    )
    invalid_base = build_aligned_segment(
        query_name="invalid_base",
        query_sequence="AAAAANAAAA",
        reference_start=100,
    )
    low_baseq = build_aligned_segment(
        query_name="low_baseq",
        query_sequence="AAAAATAAAA",
        reference_start=100,
    )
    low_baseq.query_qualities = pysam.qualitystring_to_array("IIIII+IIII")

    bam_path = create_test_bam(tmp_path, [low_mapq, invalid_base, low_baseq])

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
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


def test_collect_evidence_from_alignment_excludes_rejected_sam_flags(tmp_path: Path) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="accepted",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=99,  # 0x63: primary, mapped, first mate in a proper pair.
            ),
            build_aligned_segment(
                query_name="unpaired",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=0,  # No SAM flags: mapped but unpaired.
            ),
            build_aligned_segment(
                query_name="improper_pair",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=65,  # 0x41: paired first mate, but not a proper pair.
            ),
            build_aligned_segment(
                query_name="secondary",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=99 | 0x100,  # Add SECONDARY.
            ),
            build_aligned_segment(
                query_name="qc_fail",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=99 | 0x200,  # Add failed quality-control checks.
            ),
            build_aligned_segment(
                query_name="duplicate",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=99 | 0x400,  # Add PCR/optical DUPLICATE.
            ),
            build_aligned_segment(
                query_name="supplementary",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=99 | 0x800,  # Add SUPPLEMENTARY.
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=105,
            ref_base="A",
            alt_base="T",
        )

    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 0
    assert counts.usable == 1
    assert counts.unusable == 0


def test_collect_evidence_from_alignment_counts_overlapping_mates_once(tmp_path: Path) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="same_fragment",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=99,  # 0x63: primary, mapped, first mate in a proper pair.
                read_group="rg1",
            ),
            build_aligned_segment(
                query_name="same_fragment",
                query_sequence="AAAAAATAAA",
                reference_start=99,
                flag=147,  # 0x93: primary, mapped, second mate in a proper pair.
                read_group="rg1",
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=105,
            ref_base="A",
            alt_base="T",
        )

    # The second mate is leftmost and fetched first, but the first mate defines
    # the fragment strand when both calls agree.
    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.usable == 1
    assert counts.unusable == 0


def test_collect_evidence_from_alignment_marks_conflicting_mates_unusable(
    tmp_path: Path,
) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="conflicting_fragment",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=99,  # 0x63: first mate supports ALT.
            ),
            build_aligned_segment(
                query_name="conflicting_fragment",
                query_sequence="AAAAAAAAAA",
                reference_start=99,
                flag=147,  # 0x93: second mate supports the reference allele.
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=105,
            ref_base="A",
            alt_base="T",
        )

    assert counts.alt_forward == 0
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 0
    assert counts.usable == 0
    assert counts.unusable == 1
    assert counts.unusable_by_reason[UnusableReason.CONFLICTING_MATES] == 1


def test_collect_evidence_from_alignment_uses_usable_mate(tmp_path: Path) -> None:
    usable_first_mate = build_aligned_segment(
        query_name="partly_usable_fragment",
        query_sequence="AAAAATAAAA",
        reference_start=100,
        flag=99,  # 0x63: first mate supports ALT.
    )
    unusable_second_mate = build_aligned_segment(
        query_name="partly_usable_fragment",
        query_sequence="AAAAAATAAA",
        reference_start=99,
        flag=147,  # 0x93: second mate has low base quality at the variant.
    )
    second_mate_qualities = unusable_second_mate.query_qualities
    second_mate_qualities[6] = 10
    unusable_second_mate.query_qualities = second_mate_qualities
    bam_path = create_test_bam(tmp_path, [usable_first_mate, unusable_second_mate])

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=105,
            ref_base="A",
            alt_base="T",
        )

    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.usable == 1
    assert counts.unusable == 0


def test_collect_evidence_from_alignment_counts_two_unusable_mates_once(
    tmp_path: Path,
) -> None:
    first_mate = build_aligned_segment(
        query_name="unusable_fragment",
        query_sequence="AAAAATAAAA",
        reference_start=100,
        mapping_quality=5,
        flag=99,  # 0x63: first mate has low mapping quality.
    )
    second_mate = build_aligned_segment(
        query_name="unusable_fragment",
        query_sequence="AAAAAATAAA",
        reference_start=99,
        flag=147,  # 0x93: second mate has low base quality at the variant.
    )
    second_mate_qualities = second_mate.query_qualities
    second_mate_qualities[6] = 10
    second_mate.query_qualities = second_mate_qualities
    bam_path = create_test_bam(tmp_path, [first_mate, second_mate])

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=105,
            ref_base="A",
            alt_base="T",
        )

    assert counts.usable == 0
    assert counts.unusable == 1


def test_collect_evidence_from_alignment_requires_exact_deletion_length(
    tmp_path: Path,
) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="exact_deletion",
                query_sequence="AAAAAAAAAA",
                reference_start=100,
                flag=99,  # 0x63: primary, mapped, first mate in a proper pair.
                cigar=((0, 5), (2, 1), (0, 5)),  # 5M1D5M: exact 1-base deletion.
            ),
            build_aligned_segment(
                query_name="larger_deletion",
                query_sequence="AAAAAAAAAA",
                reference_start=100,
                flag=99,  # 0x63: primary, mapped, first mate in a proper pair.
                cigar=((0, 5), (2, 2), (0, 5)),  # 5M2D5M: deletion is too long.
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=104,
            ref_base="AT",
            alt_base="A",
        )

    assert counts.alt_forward == 1
    assert counts.non_alt_forward == 1
    assert counts.usable == 2


def test_collect_evidence_does_not_treat_soft_clip_as_insertion(tmp_path: Path) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="soft_clipped",
                query_sequence="AT",
                reference_start=100,
                cigar=((0, 1), (4, 1)),  # 1M1S, not 1M1I.
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=100,
            ref_base="A",
            alt_base="AT",
        )

    assert counts.alt_forward == 0
    assert counts.usable == 0
    assert counts.unusable == 1
    assert counts.unusable_by_reason[UnusableReason.NO_BASE_AT_SITE] == 1


def test_collect_evidence_counts_flanked_cigar_insertion(tmp_path: Path) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="insertion",
                query_sequence="ATA",
                reference_start=100,
                cigar=((0, 1), (1, 1), (0, 1)),  # 1M1I1M.
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=100,
            ref_base="A",
            alt_base="AT",
        )

    assert counts.alt_forward == 1
    assert counts.usable == 1
    assert counts.unusable == 0


def test_collect_evidence_does_not_treat_reference_skip_as_deletion(
    tmp_path: Path,
) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="reference_skip",
                query_sequence="AA",
                reference_start=100,
                cigar=((0, 1), (3, 1), (0, 1)),  # 1M1N1M, not 1M1D1M.
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=100,
            ref_base="AT",
            alt_base="A",
        )

    assert counts.alt_forward == 0
    assert counts.usable == 0
    assert counts.unusable == 1
    assert counts.unusable_by_reason[UnusableReason.NO_BASE_AT_SITE] == 1


@pytest.mark.parametrize(
    ("query_sequence", "cigar", "ref_base", "alt_base"),
    [
        ("A", ((0, 1),), "A", "AT"),  # No right flank after the anchor.
        ("AT", ((0, 1), (1, 1)), "A", "AT"),  # Terminal insertion.
        ("ATA", ((0, 1), (1, 1), (2, 1), (0, 1)), "A", "AT"),
        ("A", ((0, 1), (2, 1)), "AT", "A"),  # Terminal deletion.
        ("ATA", ((0, 1), (2, 1), (1, 1), (0, 1)), "AT", "A"),
        ("ATAA", ((0, 1), (1, 1), (0, 2)), "AT", "A"),
    ],
)
def test_collect_evidence_requires_clean_right_flank_for_indels(
    tmp_path: Path,
    query_sequence: str,
    cigar: tuple[tuple[int, int], ...],
    ref_base: str,
    alt_base: str,
) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="unflanked_or_complex",
                query_sequence=query_sequence,
                reference_start=100,
                cigar=cigar,
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=100,
            ref_base=ref_base,
            alt_base=alt_base,
        )

    assert counts.usable == 0
    assert counts.unusable == 1
    assert counts.unusable_by_reason == {UnusableReason.NO_BASE_AT_SITE: 1}


def test_collect_evidence_handles_missing_base_qualities(tmp_path: Path) -> None:
    read = build_aligned_segment(
        query_name="missing_qualities",
        query_sequence="T",
        reference_start=100,
    )
    read.query_qualities = None
    bam_path = create_test_bam(tmp_path, [read])

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=100,
            ref_base="A",
            alt_base="T",
        )

    assert counts.alt_forward == 0
    assert counts.usable == 0
    assert counts.unusable == 1
    assert counts.unusable_by_reason[UnusableReason.MISSING_BASEQ] == 1


def test_collect_evidence_handles_missing_query_sequence(tmp_path: Path) -> None:
    read = build_aligned_segment(
        query_name="missing_sequence",
        query_sequence="A",
        reference_start=100,
    )
    read.query_sequence = None
    bam_path = create_test_bam(tmp_path, [read])

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=100,
            ref_base="A",
            alt_base="T",
        )

    assert counts.alt_forward == 0
    assert counts.usable == 0
    assert counts.unusable == 1
    assert counts.unusable_by_reason[UnusableReason.NO_BASE_AT_SITE] == 1


def test_collect_evidence_from_alignment_scopes_query_names_to_read_group(
    tmp_path: Path,
) -> None:
    bam_path = create_test_bam(
        tmp_path,
        [
            build_aligned_segment(
                query_name="reused_name",
                query_sequence="AAAAATAAAA",
                reference_start=100,
                flag=99,  # 0x63: primary, mapped, first mate in a proper pair.
                read_group="rg1",
            ),
            build_aligned_segment(
                query_name="reused_name",
                query_sequence="AAAAAAAAAA",
                reference_start=100,
                flag=99,  # 0x63: a different fragment in another read group.
                read_group="rg2",
            ),
        ],
    )

    with pysam.AlignmentFile(bam_path, "rb") as alignment_file:
        counts = collect_evidence_from_alignment(
            alignment_file,
            contig="chr1",
            ref_pos0=105,
            ref_base="A",
            alt_base="T",
        )

    assert counts.alt_forward == 1
    assert counts.non_alt_forward == 1
    assert counts.usable == 2
    assert counts.unusable == 0
