import pysam
import pytest

from skua import (
    Variant, annotate_variant, annotate_variants, annotate_vcf_with_normals,
    annotate_vcf_with_pon, build_pon, read_pon_metadata,
)
from skua.pon import inspect_pon, validate_pon
from pysam import bcftools
from tests.test_evidence_pysam_integration import build_aligned_segment, create_test_bam


VARIANTS = [Variant.from_vcf_fields(contig="chr1", pos1=pos, ref="A", alt="T") for pos in (101, 102)]


@pytest.mark.parametrize("mapq,minimum,usable,reason", [
    (19, 20, 0, "low_mapq"),
    (20, 20, 1, None),
    (21, 20, 1, None),
    (254, 20, 1, None),
    (255, 20, 0, "unavailable_mapq"),
    (255, 0, 0, "unavailable_mapq"),
    (255, 256, 0, "unavailable_mapq"),
])
def test_mapq_unavailable_is_distinct_from_numeric_thresholds(tmp_path, mapq, minimum, usable, reason):
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


def test_unavailable_mapq_mate_cannot_override_a_passing_mate(tmp_path):
    bam = create_test_bam(tmp_path, [
        build_aligned_segment(query_name="pair", query_sequence="T" * 30,
                              reference_start=100, mapping_quality=255, read_group="rg1"),
        build_aligned_segment(query_name="pair", query_sequence="A" * 30,
                              reference_start=100, mapping_quality=60, is_reverse=True, read_group="rg1"),
    ])
    with pysam.AlignmentFile(bam, "rb") as alignment:
        single = [annotate_variant(alignment, variant) for variant in VARIANTS]
        batch = [evidence for _, evidence in annotate_variants(alignment, VARIANTS)]
    assert single == batch
    for evidence in batch:
        assert evidence.usable == evidence.non_alt_reverse == 1
        assert evidence.alt_forward == evidence.unusable == 0


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
                assert record.samples["case"]["SKUA_ALT_FWD"] == 2
                assert record.samples["case"]["SKUA_USABLE"] == 2
                assert record.samples["case"]["SKUA_UNUSABLE"] == 2
                assert record.info["SKUA_PON_SAMPLE_COUNT"] == 1
                assert record.info["SKUA_PON_NON_ALT_FWD"] == 2
                assert record.info["SKUA_PON_USABLE"] == 2
                assert record.info["SKUA_PON_UNUSABLE"] == 2
    assert read_pon_metadata(panel).evidence_policy_version == 6
    serialized = bcftools.view("-Ov", str(panel))
    old = serialized.replace('EvidencePolicyVersion="6"', 'EvidencePolicyVersion="5"')
    assert old != serialized
    (tmp_path / "old.vcf").write_text(old)
    incompatible = tmp_path / "old.bcf"
    bcftools.view("-Ob", "-o", str(incompatible), str(tmp_path / "old.vcf"), catch_stdout=False)
    bcftools.index(str(incompatible))
    assert inspect_pon(incompatible).evidence_policy_version == "5"
    assert not validate_pon(incompatible).valid
    with pytest.raises(ValueError, match="evidence policy.*rebuild"):
        read_pon_metadata(incompatible)
