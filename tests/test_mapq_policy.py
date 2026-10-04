import pysam
import pytest

from skua import (
    Variant, annotate_variant, annotate_variants, annotate_vcf_with_normals,
    annotate_vcf_with_pon, build_pon,
)
from tests.test_evidence_pysam_integration import build_aligned_segment, create_test_bam


VARIANTS = [Variant.from_vcf_fields(contig="chr1", pos1=pos, ref="A", alt="T") for pos in (101, 102)]


@pytest.mark.parametrize("mapq,minimum,usable,reason", [
    (0, 0, 1, None),
    (0, 20, 0, "low_mapq"),
    (19, 20, 0, "low_mapq"),
    (20, 20, 1, None),
    (21, 20, 1, None),
    (254, 20, 1, None),
    (255, 20, 1, None),
    (255, 0, 1, None),
    (255, 255, 1, None),
    (255, 256, 0, "low_mapq"),
])
def test_mapq_uses_numeric_threshold(tmp_path, mapq, minimum, usable, reason):
    bam = create_test_bam(tmp_path, [build_aligned_segment(
        query_name="fragment", query_sequence="T" * 30, reference_start=100,
        mapping_quality=mapq, read_group="rg1",
    )])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        single = [annotate_variant(alignment, variant, min_mapq=minimum) for variant in VARIANTS]
        batch = [evidence for _, evidence in annotate_variants(alignment, VARIANTS, min_mapq=minimum)]
    assert single == batch
    for evidence in batch:
        assert evidence.alt_forward == evidence.usable == usable
        assert evidence.unusable == 1 - usable
        assert {key.value: value for key, value in evidence.unusable_by_reason.items()} == (
            {} if reason is None else {reason: 1}
        )


@pytest.mark.parametrize("mate_base,mate_mapq,usable,reason", [
    ("T", 60, 1, None),
    ("A", 60, 0, "conflicting_mates"),
    ("A", 19, 1, None),
])
def test_mapq_255_follows_normal_fragment_collapsing(tmp_path, mate_base, mate_mapq, usable, reason):
    bam = create_test_bam(tmp_path, [
        build_aligned_segment(query_name="pair", query_sequence="T" * 30,
                              reference_start=100, mapping_quality=255, read_group="rg1"),
        build_aligned_segment(query_name="pair", query_sequence=mate_base * 30,
                              reference_start=100, mapping_quality=mate_mapq, is_reverse=True, read_group="rg1"),
    ])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        single = [annotate_variant(alignment, variant) for variant in VARIANTS]
        batch = [evidence for _, evidence in annotate_variants(alignment, VARIANTS)]
    assert single == batch
    for evidence in batch:
        assert evidence.usable == evidence.alt_forward == usable
        assert evidence.non_alt_reverse == evidence.alt_reverse == evidence.non_alt_forward == 0
        assert evidence.unusable == 1 - usable
        assert {key.value: value for key, value in evidence.unusable_by_reason.items()} == (
            {} if reason is None else {reason: 1}
        )


def test_live_case_and_normals_match_cached_mapq_policy(tmp_path):
    paths = []
    for sample, base in (("case", "T"), ("normal", "A")):
        directory = tmp_path / sample
        directory.mkdir()
        paths.append(create_test_bam(directory, [build_aligned_segment(
            query_name=f"fragment-{mapq}", query_sequence=base * 30, reference_start=100,
            mapping_quality=mapq, read_group="rg1",
        ) for mapq in (19, 20, 21, 255)], sample_name=sample))
    targets = tmp_path / "targets.vcf"
    targets.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t101\t.\tA\tT\t.\tPASS\t.\n"
        "chr1\t102\t.\tA\tT\t.\tPASS\t.\n"
    )
    panel = tmp_path / "panel.bcf"
    with pysam.AlignmentFile(paths[0], "rb") as case, pysam.AlignmentFile(paths[1], "rb") as normal:
        annotate_vcf_with_normals(case, targets, normal_alignments=[normal], output_path=tmp_path / "live.vcf")
        build_pon(targets, normal_alignments=[normal], output_path=panel)
        for mode, vcf in (("cached", None), ("cached-vcf", targets)):
            annotate_vcf_with_pon(case, panel, vcf_path=vcf, output_path=tmp_path / f"{mode}.vcf")
    for mode in ("live", "cached", "cached-vcf"):
        with pysam.VariantFile(str(tmp_path / f"{mode}.vcf")) as output:
            records = list(output)
            assert len(records) == 2
            for record in records:
                assert record.samples["case"]["SKUA_ALT_FWD"] == 3
                assert record.samples["case"]["SKUA_USABLE"] == 3
                assert record.samples["case"]["SKUA_UNUSABLE"] == 1
                assert record.info["SKUA_PON_SAMPLE_COUNT"] == 1
                assert record.info["SKUA_PON_NON_ALT_FWD"] == 3
                assert record.info["SKUA_PON_USABLE"] == 3
                assert record.info["SKUA_PON_UNUSABLE"] == 1
