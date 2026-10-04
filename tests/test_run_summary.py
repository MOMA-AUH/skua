import builtins
import hashlib
import json
import re
from pathlib import Path

import pysam
import pytest

import skua
from skua.pon import EVIDENCE_POLICY_VERSION, validate_pon
from tests.helpers import FakeAlignmentFile, FakeRead, build_linear_pairs
from tests.test_normal_identity import alignment, targets


def read_summary(path):
    with pysam.VariantFile(str(path)) as source:
        assert {r.key for r in source.header.records if r.key.startswith("SKUA_")} == {
            "SKUA_RUN", "SKUA_REFERENCE", "SKUA_REFERENCE_STATUS",
        }
        [record] = [r for r in source.header.records if r.key == "SKUA_RUN"]
        return {key: value.strip('"') for key, value in record.items() if key != "IDX"}


def assert_no_input_details(value):
    text = json.dumps(value)
    for omitted in ("inputs", "sha256", "size_bytes", "normal_samples", "normal_selections",
                    "panel_build", "sample_name", "read_group_ids", "path"):
        assert json.dumps(omitted) not in text
    for private in ("PRIVATE_NORMAL", "private-normal-rg", "private-case-rg", ".bam", ".bcf", ".vcf"):
        assert private not in text


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_live_annotation_records_only_effective_settings(tmp_path, suffix):
    vcf = targets(tmp_path / "targets.vcf")
    output = tmp_path / f"output{suffix}"
    with (
        alignment(tmp_path / "case.bam", "CASE", read_groups=[{"ID": "private-case-rg", "SM": "CASE"}],
                  tags=["private-case-rg"]) as case,
        alignment(tmp_path / "normal.bam", "PRIVATE_NORMAL",
                  read_groups=[{"ID": "private-normal-rg", "SM": "PRIVATE_NORMAL"}],
                  tags=["private-normal-rg"]) as normal,
    ):
        skua.annotate_vcf_with_normals(
            case, vcf, normal_alignments=[normal], output_path=output,
            min_baseq=23, min_mapq=31, truncate=0.2, pseudocount=0.01,
            prior_artifact_probability=0.3,
            assessment_thresholds=skua.AssessmentThresholds(
                min_case_depth=3, min_normal_depth=4, min_normal_samples=2,
                min_case_strand_depth=1, min_normal_strand_depth=2,
            ),
        )
    summary = read_summary(output)
    assert summary == {
        "SchemaVersion": "1", "SkuaVersion": skua.__version__, "Mode": "live_normals",
        "EvidencePolicyVersion": str(EVIDENCE_POLICY_VERSION), "MinBaseQ": "23", "MinMapQ": "31",
        "CaseReadGroups": "assigned_to_sample",
        "NormalReadGroups": "assigned_to_sample", "Truncate": "0.2", "Pseudocount": "0.01",
        "PriorPolicy": "record_info_then_fallback", "PriorFallback": "0.3",
        "MinCaseDepth": "3", "MinNormalDepth": "4", "MinNormalSamples": "2",
        "MinCaseStrandDepth": "1", "MinNormalStrandDepth": "2",
    }
    assert_no_input_details(summary)
    with pysam.VariantFile(str(output)) as annotated:
        assert "PRIVATE_NORMAL" not in str(annotated.header)
        assert "private-normal-rg" not in str(annotated.header)
        assert "private-case-rg" not in str(annotated.header)
        assert "##SKUA_REFERENCE_STATUS=INSUFFICIENT_METADATA" in str(annotated.header)
        assert next(annotated).info["SKUA_ARTIFACT_PRIOR"] == pytest.approx((0.3,))


@pytest.mark.parametrize("separate_vcf", [False, True])
def test_cached_summary_uses_panel_thresholds_without_exporting_membership_or_history(tmp_path, separate_vcf):
    vcf = targets(tmp_path / "targets.vcf")
    vcf.write_text(vcf.read_text().replace(
        "#CHROM", '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior">\n#CHROM',
    ).replace("PASS\t.", "PASS\tSKUA_ARTIFACT_PRIOR=0.8"))
    panel = tmp_path / "panel.bcf"
    with alignment(tmp_path / "normal.bam", "PRIVATE_NORMAL") as normal:
        skua.build_pon(vcf, normal_alignments=[normal], output_path=panel, min_baseq=25, min_mapq=35)
    assert skua.read_pon_metadata(panel).sample_names == ("PRIVATE_NORMAL",)
    assert validate_pon(panel).valid
    output = tmp_path / "out.vcf"
    with alignment(tmp_path / "case.bam", "CASE") as case:
        skua.annotate_vcf_with_pon(
            case, panel, vcf_path=vcf if separate_vcf else None, output_path=output,
            truncate=0.25, pseudocount=0.02, prior_artifact_probability=0.4,
        )
    summary = read_summary(output)
    assert summary["Mode"] == "cached_pon"
    assert summary["MinBaseQ"] == "25"
    assert summary["MinMapQ"] == "35"
    assert summary["Truncate"] == "0.25"
    assert summary["Pseudocount"] == "0.02"
    assert summary["PriorFallback"] == "0.4"
    assert_no_input_details(summary)
    with pysam.VariantFile(str(output)) as result:
        assert tuple(result.header.samples) == ("CASE",)
        assert "PRIVATE_NORMAL" not in str(result.header)
        assert "private-normal-rg" not in str(result.header)
        assert next(result).info["SKUA_ARTIFACT_PRIOR"] == pytest.approx((0.8,))


@pytest.mark.parametrize("unknown_metadata", [False, True])
def test_force_replaces_metadata_and_preserves_unrelated_annotations(tmp_path, unknown_metadata):
    vcf = targets(tmp_path / "targets.vcf")
    vcf.write_text(vcf.read_text().replace(
        "#CHROM", '##external=/upstream/private/input.vcf\n##INFO=<ID=KEEP,Number=1,Type=Integer,Description="Keep">\n#CHROM',
    ).replace("PASS\t.", "PASS\tKEEP=7"))
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf.gz"
    with alignment(tmp_path / 'case, "ø".bam', "CASE") as case:
        skua.annotate_vcf_with_normals(case, vcf, output_path=first, min_mapq=25)
        if unknown_metadata:
            first.write_text(re.sub(r"^##SKUA_RUN=.*$", "##SKUA_UNUSED=discard_me",
                                   first.read_text(), flags=re.MULTILINE))
        with pytest.raises(ValueError, match="already contains Skua annotations"):
            skua.annotate_vcf_with_normals(case, first, output_path=second)
        skua.annotate_vcf_with_normals(case, first, output_path=second, min_mapq=40, force=True)
    assert read_summary(second)["MinMapQ"] == "40"
    with pysam.VariantFile(str(second)) as result:
        assert "##external=/upstream/private/input.vcf" in str(result.header)
        assert next(result).info["KEEP"] == 7


def test_case_only_summary_keeps_reference_metadata_without_reference_filename(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    reference = tmp_path / "private-reference.fa"
    reference.write_text(">chr1\n" + "A" * 200 + "\n")
    pysam.faidx(str(reference))
    output = tmp_path / "case-only.vcf"
    with alignment(tmp_path / "case.bam", "CASE") as case:
        skua.annotate_vcf(case, vcf, output_path=output, reference_path=reference)
    summary = read_summary(output)
    assert summary["Mode"] == "case_only"
    assert summary["NormalReadGroups"] == "not_applicable"
    assert "Truncate" not in summary
    assert "MinCaseDepth" not in summary
    with pysam.VariantFile(str(output)) as result:
        assert "##SKUA_REFERENCE=<ID=chr1" in str(result.header)
        assert "private-reference.fa" not in str(result.header)


@pytest.mark.parametrize("with_normals", [False, True])
def test_json_document_retains_minimal_summary_and_actual_python_selection_policy(tmp_path, with_normals):
    vcf = targets(tmp_path / "targets.vcf")
    output = tmp_path / "result.json"
    with alignment(tmp_path / "case.bam", "CASE") as case, alignment(tmp_path / "normal.bam", "PRIVATE_NORMAL") as normal:
        if with_normals:
            payload = skua.annotate_vcf_to_json_with_normals(
                case, vcf, normal_alignments=[normal], output_path=output,
                min_baseq=24, min_mapq=32, truncate=0.2,
                pseudocount=0.01, prior_artifact_probability=0.3,
            )
        else:
            payload = skua.annotate_vcf_to_json(case, vcf, output_path=output, min_baseq=24, min_mapq=32)
    document = json.loads(payload)
    assert json.loads(output.read_text()) == document
    assert set(document) == {"run_summary", "records"}
    assert len(document["records"]) == 1
    summary = document["run_summary"]
    assert set(summary) == {"schema_version", "skua_version", "mode", "evidence", "model", "case_read_groups", "reference"}
    assert summary["schema_version"] == 1
    assert summary["skua_version"] == skua.__version__
    assert summary["evidence"] == {
        "policy_version": EVIDENCE_POLICY_VERSION, "min_baseq": 24, "min_mapq": 32,
        "normal_read_groups": "assigned_to_sample" if with_normals else "not_applicable",
    }
    assert summary["case_read_groups"] == "all_alignment_reads"
    assert_no_input_details(document)
    if with_normals:
        assert summary["model"]["prior"] == {"policy": "constant", "fallback": 0.3}
        assert document["records"][0]["artifact_prior"] == 0.3
    else:
        assert summary["model"] is None


def test_python_objects_without_backing_files_keep_summary_without_input_details(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    document = json.loads(skua.annotate_vcf_to_json(FakeAlignmentFile([]), vcf))
    assert document["run_summary"]["case_read_groups"] == "all_alignment_reads"
    assert_no_input_details(document)


def test_mixed_headerless_and_named_normals_record_selection_policy_without_identities(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    headerless = FakeAlignmentFile([FakeRead(
        mapping_quality=60, is_reverse=False, query_sequence="A" * 30,
        query_qualities=[35] * 30, aligned_pairs=build_linear_pairs(30, 10),
    )])
    with alignment(tmp_path / "normal.bam", "PRIVATE_NORMAL") as normal:
        document = json.loads(skua.annotate_vcf_to_json_with_normals(
            FakeAlignmentFile([]), vcf, normal_alignments=[headerless, normal],
        ))
    assert document["records"][0]["counts"]["normal"]["usable"] == 2
    assert document["run_summary"]["evidence"]["normal_read_groups"] == "per_normal_selection"
    assert_no_input_details(document)


@pytest.mark.parametrize("mode", ["case", "live", "build", "cached", "cached_targets", "json", "json_normals"])
def test_metadata_does_not_hash_or_read_raw_input_files(tmp_path, monkeypatch, mode):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 200 + "\n")
    pysam.faidx(str(reference))
    with alignment(tmp_path / "case.bam", "CASE") as case, alignment(tmp_path / "normal.bam", "NORMAL") as normal:
        if mode.startswith("cached"):
            skua.build_pon(vcf, normal_alignments=[normal], output_path=panel)
        original_open = builtins.open

        def no_hash(*args, **kwargs):
            pytest.fail("Metadata must not hash input files")

        def no_raw_input_read(file, mode="r", *args, **kwargs):
            if mode == "rb" and Path(file).suffix in {".bam", ".bcf", ".fa", ".vcf"}:
                pytest.fail("Metadata must not make a raw input-file read")
            return original_open(file, mode, *args, **kwargs)

        monkeypatch.setattr(hashlib, "sha256", no_hash)
        monkeypatch.setattr(builtins, "open", no_raw_input_read)
        options = dict(output_path=tmp_path / "out.vcf", reference_path=reference)
        if mode == "case":
            skua.annotate_vcf(case, vcf, **options)
        elif mode == "live":
            skua.annotate_vcf_with_normals(case, vcf, normal_alignments=[normal], **options)
        elif mode == "build":
            skua.build_pon(vcf, normal_alignments=[normal], output_path=panel, reference_path=reference)
        elif mode.startswith("cached"):
            skua.annotate_vcf_with_pon(case, panel, vcf_path=vcf if mode == "cached_targets" else None, **options)
        elif mode == "json":
            skua.annotate_vcf_to_json(case, vcf)
        else:
            skua.annotate_vcf_to_json_with_normals(case, vcf, normal_alignments=[normal])


@pytest.mark.parametrize("force", [False, True])
def test_pon_build_rejects_or_replaces_unknown_skua_metadata(tmp_path, force):
    vcf = targets(tmp_path / "targets.vcf")
    vcf.write_text(vcf.read_text().replace("#CHROM", "##SKUA_UNUSED=discard_me\n#CHROM"))
    panel = tmp_path / "panel.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL") as normal:
        if not force:
            with pytest.raises(ValueError, match="already contains Skua annotations: SKUA_UNUSED"):
                skua.build_pon(vcf, normal_alignments=[normal], output_path=panel)
            assert not panel.exists()
            return
        skua.build_pon(vcf, normal_alignments=[normal], output_path=panel, force=True)
    assert validate_pon(panel).valid
    with pysam.VariantFile(str(panel)) as source:
        assert {r.key for r in source.header.records if r.key.startswith("SKUA_")} == {
            "SKUA_PON", "SKUA_REFERENCE", "SKUA_REFERENCE_STATUS",
        }
