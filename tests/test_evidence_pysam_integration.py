from pathlib import Path

import pysam

from skua.evidence import UnusableReason, collect_evidence_from_alignment


HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000}],
}

def build_aligned_segment(
    *,
    query_name: str,
    query_sequence: str,
    reference_start: int,
    mapping_quality: int = 60,
    is_reverse: bool = False,
    flag: int | None = None,
) -> pysam.AlignedSegment:
    segment = pysam.AlignedSegment()
    segment.query_name = query_name
    segment.query_sequence = query_sequence
    segment.flag = flag if flag is not None else (147 if is_reverse else 99)
    segment.reference_id = 0
    segment.reference_start = reference_start
    segment.mapping_quality = mapping_quality
    segment.cigar = ((0, len(query_sequence)),)
    segment.next_reference_id = 0
    segment.next_reference_start = reference_start
    segment.template_length = -len(query_sequence) if is_reverse else len(query_sequence)
    segment.query_qualities = pysam.qualitystring_to_array("I" * len(query_sequence))
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
