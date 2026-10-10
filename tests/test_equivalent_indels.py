from pathlib import Path

import pysam
import pytest

from skua import (
    annotate_variant, annotate_variants, annotate_vcf_with_normals,
    annotate_vcf_with_pon, build_pon,
)
from skua.variants import Variant
from skua.evidence import UnusableReason
from tests.test_evidence_pysam_integration import build_aligned_segment, create_test_bam


def reference_file(tmp_path: Path, sequence: str) -> Path:
    path = tmp_path / "reference.fa"
    path.write_text(">chr1\n" + sequence + "G" * (1000 - len(sequence)) + "\n")
    pysam.faidx(str(path))
    return path


def test_shifted_homopolymer_insertion_supports_alt(tmp_path: Path) -> None:
    reference = reference_file(tmp_path, "CAAAAT")
    read = build_aligned_segment(
        query_name="insertion", query_sequence="CAAAAAT", reference_start=0,
        cigar=((0, 3), (1, 1), (0, 3)),
    )
    bam = create_test_bam(tmp_path, [read])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        evidence = annotate_variant(
            alignment, Variant("chr1", 0, "C", "CA"), reference_path=reference,
        )
    assert evidence.alt_forward == 1
    assert evidence.non_alt_forward == evidence.unusable == 0


@pytest.mark.parametrize(
    "sequence,read_sequence,cigar,variant",
    [
        ("CAAAAT", "CAAAAAT", "3M1I3M", Variant("chr1", 0, "C", "CA")),
        ("CAAAAT", "CAAAAAT", "1M1I5M", Variant("chr1", 3, "A", "AA")),
        ("CAAAAT", "CAAAT", "3M1D2M", Variant("chr1", 0, "CA", "C")),
        ("CAAAAT", "CAAT", "1M2D3M", Variant("chr1", 1, "AAA", "A")),
        ("CATATATG", "CATATATATG", "2M2I6M", Variant("chr1", 0, "C", "CAT")),
        ("CATATATG", "CATATATATG", "1M2I7M", Variant("chr1", 3, "A", "ATA")),
        ("CATATATG", "CATATG", "2M2D4M", Variant("chr1", 0, "CAT", "C")),
        ("CATATATG", "CATATG", "1M2D5M", Variant("chr1", 3, "ATA", "A")),
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
def test_shifted_repeat_indels_agree_in_singleton_and_batch(
    tmp_path: Path, sequence: str, read_sequence: str, cigar: str, variant: Variant,
    reverse: bool,
) -> None:
    reference = reference_file(tmp_path, sequence)
    read = build_aligned_segment(
        query_name="repeat", query_sequence=read_sequence, reference_start=0,
        is_reverse=reverse,
    )
    read.cigarstring = cigar
    bam = create_test_bam(tmp_path, [read])
    variants = [Variant("chr1", 0, "C", "G"), variant]
    with pysam.AlignmentFile(bam, "rb") as alignment:
        single = annotate_variant(alignment, variant, reference_path=reference)
        batch = list(annotate_variants(alignment, variants, reference_path=reference))
    assert batch[1] == (variant, single)
    assert single.alt_reverse == int(reverse)
    assert single.alt_forward == int(not reverse)
    assert single.usable == 1
    assert single.unusable == 0


def test_shifted_homopolymer_deletion_supports_alt(tmp_path: Path) -> None:
    reference = reference_file(tmp_path, "CAAAAT")
    read = build_aligned_segment(
        query_name="deletion", query_sequence="CAAAT", reference_start=0,
        cigar=((0, 3), (2, 1), (0, 2)),
    )
    bam = create_test_bam(tmp_path, [read])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        evidence = annotate_variant(
            alignment, Variant("chr1", 0, "CA", "C"), reference_path=reference,
        )
    assert evidence.alt_forward == 1
    assert evidence.non_alt_forward == evidence.unusable == 0


@pytest.mark.parametrize("separate_vcf", [False, True])
def test_shifted_indels_direct_and_cached_counts_agree(tmp_path: Path, separate_vcf: bool) -> None:
    reference = reference_file(tmp_path, "CAAAAT")
    vcf = tmp_path / "targets.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1\t.\tC\tCA\t.\tPASS\t.\n"
        "chr1\t1\t.\tCA\tC\t.\tPASS\t.\n"
    )
    paths = []
    for sample in ("case", "normal"):
        directory = tmp_path / sample
        directory.mkdir()
        reads = []
        for name, sequence, cigar in (
            ("insertion", "CAAAAAT", "3M1I3M"),
            ("deletion", "CAAAT", "3M1D2M"),
            ("reference", "CAAAAT", "6M"),
        ):
            read = build_aligned_segment(
                query_name=name, query_sequence=sequence, reference_start=0,
                read_group="rg1",
            )
            read.cigarstring = cigar
            reads.append(read)
        paths.append(create_test_bam(directory, reads, sample_name=sample))
    pon, direct, cached = (tmp_path / name for name in ("pon.bcf", "direct.vcf", "cached.vcf"))
    with pysam.AlignmentFile(paths[0], "rb") as case, pysam.AlignmentFile(paths[1], "rb") as normal:
        build_pon(vcf, normal_alignments=[normal], output_path=pon, reference_path=reference)
        annotate_vcf_with_normals(
            case, vcf, normal_alignments=[normal], output_path=direct, reference_path=reference,
            truncate=1.0,
        )
        annotate_vcf_with_pon(
            case, pon, vcf_path=vcf if separate_vcf else None,
            output_path=cached, reference_path=reference,
            truncate=1.0,
        )
        with pytest.raises(ValueError, match="requires.*reference"):
            annotate_vcf_with_pon(case, pon, output_path=tmp_path / "missing-reference.vcf")
    with pysam.VariantFile(direct) as live, pysam.VariantFile(cached) as stored:
        for observed, expected in zip(stored, live, strict=True):
            assert observed.samples["case"]["SKUA_ALT_FWD"] == 1
            assert observed.info["SKUA_PON_ALT_FWD"] == 1
            assert dict(observed.samples["case"]) == dict(expected.samples["case"])
            assert dict(observed.info) == dict(expected.info)
    with pysam.AlignmentFile(paths[0], "rb") as case, pysam.AlignmentFile(paths[1], "rb") as normal:
        exact_pon, exact_output = tmp_path / "exact.bcf", tmp_path / "exact.vcf"
        build_pon(vcf, normal_alignments=[normal], output_path=exact_pon)
        annotate_vcf_with_pon(
            case, exact_pon, vcf_path=vcf if separate_vcf else None,
            output_path=exact_output, reference_path=reference,
        )
    with pysam.VariantFile(exact_output) as stored:
        assert 'IndelMatching="exact_anchor"' in str(stored.header)
        assert [record.samples["case"]["SKUA_ALT_FWD"] for record in stored] == [0, 0]


@pytest.mark.parametrize("sequence,cigar,ref,alt", [
    ("CAAA", "3M1I", "C", "CA"),
    ("CAA", "3M1D", "CA", "C"),
    ("CAAAAT", "3M1I1D2M", "C", "CA"),
])
def test_shifted_indel_requires_clean_aligned_right_flank(
    tmp_path: Path, sequence: str, cigar: str, ref: str, alt: str,
) -> None:
    reference = reference_file(tmp_path, "CAAAAT")
    read = build_aligned_segment(query_name="complex", query_sequence=sequence, reference_start=0)
    read.cigarstring = cigar
    bam = create_test_bam(tmp_path, [read])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        evidence = annotate_variant(alignment, Variant("chr1", 0, ref, alt), reference_path=reference)
    assert evidence.usable == 0
    assert evidence.unusable_by_reason == {UnusableReason.NO_BASE_AT_SITE: 1}


@pytest.mark.parametrize("kind,position", [
    ("insertion", position) for position in range(5)
] + [("deletion", position) for position in range(4)])
@pytest.mark.parametrize("quality,reason", [(5, UnusableReason.LOW_BASEQ), (255, UnusableReason.MISSING_BASEQ)])
def test_shifted_indel_checks_every_comparison_base_quality(
    tmp_path: Path, kind: str, position: int, quality: int, reason: UnusableReason,
) -> None:
    reference = reference_file(tmp_path, "CAAAAT")
    insertion = kind == "insertion"
    variant = Variant("chr1", 0, "C" if insertion else "CA", "CA" if insertion else "C")
    read = build_aligned_segment(
        query_name="quality", query_sequence="CAAAAAT" if insertion else "CAAAT", reference_start=0,
        cigar=((0, 3), (1, 1), (0, 3)) if insertion else ((0, 3), (2, 1), (0, 2)),
    )
    qualities = list(read.query_qualities)
    qualities[position] = quality
    read.query_qualities = qualities
    bam = create_test_bam(tmp_path, [read])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        single = annotate_variant(alignment, variant, reference_path=reference)
        batch = list(annotate_variants(alignment, [variant, variant], reference_path=reference))
    assert batch == [(variant, single), (variant, single)]
    assert single.usable == 0
    assert single.unusable_by_reason == {reason: 1}


@pytest.mark.parametrize("sequence,read_sequence,cigar", [
    ("CAAAAT", "CAACAAT", "3M1I3M"),  # Wrong inserted allele.
    ("CAGGGT", "CAAAAAT", "3M1I3M"),  # Read mismatches cannot prove reference equivalence.
    ("CAAAAT", "CAAAAT", "6M"),  # An explicit event is necessary.
    ("CAAAAT", "CAAAAT", "3M1S2M"),  # Soft clipping is not insertion.
    ("CANAAT", "CAAAAAT", "3M1I3M"),  # Ambiguous reference context.
])
def test_shifted_insertion_rejects_non_equivalent_or_implicit_events(
    tmp_path: Path, sequence: str, read_sequence: str, cigar: str,
) -> None:
    reference = reference_file(tmp_path, sequence)
    read = build_aligned_segment(query_name="other", query_sequence=read_sequence, reference_start=0)
    read.cigarstring = cigar
    bam = create_test_bam(tmp_path, [read])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        evidence = annotate_variant(alignment, Variant("chr1", 0, "C", "CA"), reference_path=reference)
    assert evidence.alt_forward == evidence.alt_reverse == 0


@pytest.mark.parametrize("distance,expected_alt", [(100, 1), (101, 0)])
def test_shifted_matching_has_a_bounded_search(tmp_path: Path, distance: int, expected_alt: int) -> None:
    reference = reference_file(tmp_path, "C" + "A" * 110 + "T")
    read = build_aligned_segment(
        query_name="long-repeat", query_sequence="C" + "A" * 111 + "T", reference_start=0,
        cigar=((0, distance + 1), (1, 1), (0, 111 - distance)),
    )
    bam = create_test_bam(tmp_path, [read])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        variant = Variant("chr1", 0, "C", "CA")
        assert annotate_variant(alignment, variant, reference_path=reference).alt_forward == expected_alt
        assert annotate_variant(alignment, variant).alt_forward == 0


@pytest.mark.parametrize("mate_sequence,mate_cigar,mapq,expected_alt,reason", [
    ("CAAAAAT", "1M1I5M", 60, 1, None),
    ("CAAAAT", "6M", 60, 0, UnusableReason.CONFLICTING_MATES),
    ("CAAAAT", "6M", 10, 1, None),
])
def test_shifted_evidence_preserves_fragment_resolution(
    tmp_path: Path, mate_sequence: str, mate_cigar: str, mapq: int,
    expected_alt: int, reason: UnusableReason | None,
) -> None:
    reference = reference_file(tmp_path, "CAAAAT")
    first = build_aligned_segment(
        query_name="pair", query_sequence="CAAAAAT", reference_start=0,
        cigar=((0, 3), (1, 1), (0, 3)),
    )
    second = build_aligned_segment(
        query_name="pair", query_sequence=mate_sequence, reference_start=0,
        mapping_quality=mapq, is_reverse=True,
    )
    second.cigarstring = mate_cigar
    bam = create_test_bam(tmp_path, [first, second])
    variant = Variant("chr1", 0, "C", "CA")
    with pysam.AlignmentFile(bam, "rb") as alignment:
        single = annotate_variant(alignment, variant, reference_path=reference)
        batch = list(annotate_variants(alignment, [variant, variant], reference_path=reference))
    assert batch == [(variant, single), (variant, single)]
    assert single.alt_forward == expected_alt
    assert single.alt_reverse == 0
    assert single.usable + single.unusable == 1
    assert single.unusable_by_reason == ({} if reason is None else {reason: 1})


@pytest.mark.parametrize("base,quality,reason", [
    ("C", 5, UnusableReason.LOW_BASEQ),
    ("C", 255, UnusableReason.MISSING_BASEQ),
    ("N", 40, UnusableReason.INVALID_BASE),
])
def test_uncertain_shifted_insertion_is_not_reference_evidence(
    tmp_path: Path, base: str, quality: int, reason: UnusableReason,
) -> None:
    reference = reference_file(tmp_path, "CAAAAT")
    read = build_aligned_segment(
        query_name="uncertain", query_sequence="CAA" + base + "AAT", reference_start=0,
        cigar=((0, 3), (1, 1), (0, 3)),
    )
    read.query_qualities = [40, 40, 40, quality, 40, 40, 40]
    bam = create_test_bam(tmp_path, [read])
    variant = Variant("chr1", 0, "C", "CA")
    with pysam.AlignmentFile(bam, "rb") as alignment:
        single = annotate_variant(alignment, variant, reference_path=reference)
        batch = list(annotate_variants(alignment, [variant, variant], reference_path=reference))
    assert batch == [(variant, single), (variant, single)]
    assert single.usable == 0
    assert single.unusable_by_reason == {reason: 1}
