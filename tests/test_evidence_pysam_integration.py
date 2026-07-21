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
) -> pysam.AlignedSegment:
    segment = pysam.AlignedSegment()
    segment.query_name = query_name
    segment.query_sequence = query_sequence
    segment.flag = 16 if is_reverse else 0
    segment.reference_id = 0
    segment.reference_start = reference_start
    segment.mapping_quality = mapping_quality
    segment.cigar = ((0, len(query_sequence)),)
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
