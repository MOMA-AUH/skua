from contextlib import ExitStack
from pathlib import Path

import pysam
import pytest
from pysam import bcftools

from skua import (
    Variant, annotate_variant_with_normals, annotate_variants_with_normals,
    annotate_vcf_with_normals, annotate_vcf_with_pon, annotate_variants_from_pon,
    build_pon, read_pon_evidence, read_pon_metadata,
)
from skua.pon import inspect_pon, validate_pon
from tests.helpers import FakeAlignmentFile


def targets(path: Path) -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t21\t.\tA\tT\t.\tPASS\t.\n", encoding="utf-8",
    )
    return path


def alignment(path: Path, sample: str, *, read_groups=None, tags=None):
    groups = [{"ID": sample, "SM": sample}] if read_groups is None else read_groups
    with pysam.AlignmentFile(str(path), "wb", header={
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 200}], "RG": groups,
    }) as output:
        for index, tag in enumerate([sample] if tags is None else tags):
            read = pysam.AlignedSegment(output.header)
            read.query_name = f"fragment-{index}"
            read.flag = 99
            read.reference_id = 0
            read.reference_start = 10
            read.next_reference_id = 0
            read.next_reference_start = 50
            read.template_length = 70
            read.mapping_quality = 60
            read.cigarstring = "30M"
            read.query_sequence = "A" * 30
            read.query_qualities = [35] * 30
            if tag is not None:
                read.set_tag("RG", tag)
            output.write(read)
    pysam.index(str(path))
    return pysam.AlignmentFile(str(path), "rb")


@pytest.mark.parametrize("mode", ["single", "batch", "vcf", "build"])
@pytest.mark.parametrize("duplicate", ["handle", "path", "symlink", "hardlink", "sample"])
def test_repeated_normal_is_rejected_before_it_can_duplicate_evidence(tmp_path, mode, duplicate):
    vcf = targets(tmp_path / "targets.vcf")
    variant = Variant.from_vcf_fields(contig="chr1", pos1=21, ref="A", alt="T")
    with ExitStack() as stack:
        case = stack.enter_context(alignment(tmp_path / "case.bam", "CASE"))
        normal = stack.enter_context(alignment(tmp_path / "normal.bam", "NORMAL"))
        repeated = normal
        if duplicate == "sample":
            repeated = stack.enter_context(alignment(tmp_path / "other.bam", "NORMAL"))
        elif duplicate != "handle":
            path = tmp_path / "normal.bam"
            if duplicate in {"symlink", "hardlink"}:
                path = tmp_path / "alias.bam"
                getattr(path, "symlink_to" if duplicate == "symlink" else "hardlink_to")(tmp_path / "normal.bam")
                Path(f"{path}.bai").symlink_to(tmp_path / "normal.bam.bai")
            repeated = stack.enter_context(pysam.AlignmentFile(str(path), "rb"))
        options = {"normal_alignments": [normal, repeated]}
        with pytest.raises(ValueError, match="[Dd]uplicate|must be unique"):
            if mode == "single":
                annotate_variant_with_normals(case, variant, **options)
            elif mode == "batch":
                list(annotate_variants_with_normals(case, [variant], **options))
            elif mode == "vcf":
                annotate_vcf_with_normals(case, vcf, output_path=tmp_path / "out.vcf", **options)
            else:
                build_pon(vcf, output_path=tmp_path / "panel.bcf", **options)
    assert not (tmp_path / "out.vcf").exists()
    assert not (tmp_path / "panel.bcf").exists()


@pytest.mark.parametrize("mode", ["single", "batch", "vcf", "cached", "cached-vcf", "cached-iterator"])
def test_case_cannot_be_a_member_of_its_own_panel(tmp_path, mode):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    variant = Variant.from_vcf_fields(contig="chr1", pos1=21, ref="A", alt="T")
    with (
        alignment(tmp_path / "case.bam", "CASE") as case,
        alignment(tmp_path / "normal.bam", "CASE") as normal,
    ):
        build_pon(vcf, normal_alignments=[normal], output_path=panel)
        with pytest.raises(ValueError, match="[Cc]ase.*(panel|normal)"):
            if mode == "single":
                annotate_variant_with_normals(case, variant, normal_alignments=[normal])
            elif mode == "batch":
                list(annotate_variants_with_normals(case, [variant], normal_alignments=[normal]))
            elif mode == "vcf":
                annotate_vcf_with_normals(case, vcf, normal_alignments=[normal], output_path=tmp_path / "out.vcf")
            elif mode == "cached-iterator":
                list(annotate_variants_from_pon(case, panel))
            else:
                annotate_vcf_with_pon(
                    case, panel, vcf_path=vcf if mode == "cached-vcf" else None,
                    output_path=tmp_path / "out.vcf",
                )
    assert not (tmp_path / "out.vcf").exists()


@pytest.mark.parametrize("cached", [False, True])
def test_membership_checks_only_the_selected_case_sample(tmp_path, cached):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    with (
        alignment(tmp_path / "case.bam", "CASE", read_groups=[
            {"ID": "CASE", "SM": "CASE"}, {"ID": "NORMAL", "SM": "NORMAL"},
        ]) as case,
        alignment(tmp_path / "normal.bam", "NORMAL") as normal,
    ):
        if cached:
            build_pon(vcf, normal_alignments=[normal], output_path=panel)
            annotate_vcf_with_pon(case, panel, sample_name="CASE", output_path=tmp_path / "out.vcf")
        else:
            annotate_vcf_with_normals(case, vcf, normal_alignments=[normal], sample_name="CASE", output_path=tmp_path / "out.vcf")
    with pysam.VariantFile(str(tmp_path / "out.vcf")) as output:
        record = next(output)
        assert record.info["SKUA_PON_SAMPLE_COUNT"] == 1
        assert record.info["SKUA_PON_USABLE"] == 1


@pytest.mark.parametrize("mode", ["single", "batch", "vcf", "cached", "cached-vcf"])
def test_distinct_normals_count_only_reads_assigned_to_their_sample(tmp_path, mode):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    variant = Variant.from_vcf_fields(contig="chr1", pos1=21, ref="A", alt="T")
    with (
        alignment(tmp_path / "case.bam", "CASE") as case,
        alignment(tmp_path / "normal.bam", "NORMAL", read_groups=[
            {"ID": "lib1", "SM": "NORMAL"}, {"ID": "lib2", "SM": "NORMAL"},
            {"ID": "unassigned"},
        ], tags=["lib1", "lib2", "unassigned", "unknown", None]) as normal,
        alignment(tmp_path / "normal2.bam", "NORMAL2") as normal2,
    ):
        normals = [normal, normal2]
        if mode == "single":
            result = annotate_variant_with_normals(case, variant, normal_alignments=normals)
        elif mode == "batch":
            [(_, result)] = annotate_variants_with_normals(case, [variant], normal_alignments=normals)
        else:
            if mode == "vcf":
                annotate_vcf_with_normals(case, vcf, normal_alignments=normals, output_path=tmp_path / "out.vcf")
            else:
                build_pon(vcf, normal_alignments=normals, output_path=panel)
                [(_, evidence)] = read_pon_evidence(panel)
                assert [item.usable for item in evidence] == [2, 1]
                annotate_vcf_with_pon(case, panel, vcf_path=vcf if mode == "cached-vcf" else None, output_path=tmp_path / "out.vcf")
            with pysam.VariantFile(str(tmp_path / "out.vcf")) as output:
                record = next(output)
                assert record.info["SKUA_PON_SAMPLE_COUNT"] == 2
                assert record.info["SKUA_PON_USABLE"] == 3
                assert record.info["SKUA_PON_UNUSABLE"] == 0
            return
        assert [item.usable for item in result.normal_evidences] == [2, 1]
        assert result.normal_aggregate_evidence.usable == 3
        assert result.normal_aggregate_evidence.unusable == 0


@pytest.mark.parametrize("groups", [
    [{"SM": "NORMAL"}],
    [{"ID": "duplicate", "SM": "NORMAL"}, {"ID": "duplicate"}],
])
def test_normal_without_unambiguous_read_group_ids_is_rejected(tmp_path, groups):
    vcf = targets(tmp_path / "targets.vcf")
    with alignment(tmp_path / "normal.bam", "NORMAL", read_groups=groups) as normal:
        with pytest.raises(ValueError, match="read-group ID|could not parse header"):
            build_pon(vcf, normal_alignments=[normal], output_path=tmp_path / "panel.bcf")


@pytest.mark.parametrize("duplicate", ["handle", "path", "symlink", "hardlink"])
def test_duplicates_are_rejected_even_without_sample_metadata(tmp_path, duplicate):
    normal = FakeAlignmentFile([])
    original = tmp_path / "normal.bam"
    original.touch()
    normal.filename = bytes(original)
    repeated = normal
    if duplicate != "handle":
        repeated = FakeAlignmentFile([])
        path = original
        if duplicate in {"symlink", "hardlink"}:
            path = tmp_path / "alias.bam"
            getattr(path, "symlink_to" if duplicate == "symlink" else "hardlink_to")(original)
        repeated.filename = bytes(path)
    variant = Variant.from_vcf_fields(contig="chr1", pos1=21, ref="A", alt="T")
    with pytest.raises(ValueError, match="Duplicate normal alignment"):
        annotate_variant_with_normals(FakeAlignmentFile([]), variant, normal_alignments=[normal, repeated])
    with pytest.raises(ValueError, match="Case alignment.*normal"):
        annotate_variant_with_normals(normal, variant, normal_alignments=[repeated])


@pytest.mark.parametrize("old_policy", ["4", "7"])
def test_schema_two_panel_with_old_normal_policy_requires_rebuild(tmp_path, old_policy):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    incompatible = tmp_path / "old-policy.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL") as normal:
        build_pon(vcf, normal_alignments=[normal], output_path=panel)
    serialized = bcftools.view("-Ov", str(panel))
    old = serialized.replace('EvidencePolicyVersion="8"', f'EvidencePolicyVersion="{old_policy}"')
    assert old != serialized
    (tmp_path / "old.vcf").write_text(old)
    bcftools.view("-Ob", "-o", str(incompatible), str(tmp_path / "old.vcf"), catch_stdout=False)
    bcftools.index(str(incompatible))
    assert inspect_pon(incompatible).evidence_policy_version == old_policy
    assert not validate_pon(incompatible).valid
    with pytest.raises(ValueError, match="evidence policy.*rebuild"):
        read_pon_metadata(incompatible)
    with alignment(tmp_path / "case.bam", "CASE") as case:
        with pytest.raises(ValueError, match="evidence policy.*rebuild"):
            annotate_vcf_with_pon(case, incompatible, output_path=tmp_path / "out.vcf")
    assert not (tmp_path / "out.vcf").exists()
