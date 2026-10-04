import hashlib
import json
import base64
import re

import pysam
import pytest
from pysam import bcftools

import skua
from tests.test_normal_identity import alignment, targets
from skua.pon import inspect_pon, validate_pon
from tests.helpers import FakeAlignmentFile


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_live_annotation_records_effective_settings_and_input_content_identities(tmp_path, suffix):
    vcf = targets(tmp_path / "targets.vcf")
    case_path, normal_path = tmp_path / "case.bam", tmp_path / "normal.bam"
    output = tmp_path / f"output{suffix}"
    with alignment(case_path, "CASE") as case, alignment(normal_path, "NORMAL") as normal:
        skua.annotate_vcf_with_normals(
            case, vcf, normal_alignments=[normal], output_path=output,
            min_baseq=23, min_mapq=31, truncate=0.2, pseudocount=0.01,
            prior_artifact_probability=0.3,
            assessment_thresholds=skua.AssessmentThresholds(min_case_depth=3),
        )
    provenance = skua.read_provenance(output)
    assert provenance["schema_version"] == 1
    assert provenance["skua_version"] == skua.__version__
    assert provenance["mode"] == "live_normals"
    assert provenance["evidence"] == {
        "policy_version": 6, "min_baseq": 23, "min_mapq": 31,
        "mapq_255": "exclude", "normal_read_groups": "assigned_to_sample",
    }
    assert provenance["model"]["truncate"] == 0.2
    assert provenance["model"]["pseudocount"] == 0.01
    assert provenance["model"]["prior"] == {"policy": "record_info_then_fallback", "fallback": 0.3}
    assert provenance["model"]["assessment_thresholds"]["min_case_depth"] == 3
    assert provenance["case"]["sample_name"] == "CASE"
    assert provenance["case"]["read_group_ids"] == ["CASE"]
    assert provenance["normal_samples"] == ["NORMAL"]
    for identity, path in (
        (provenance["inputs"]["case"], case_path),
        (provenance["inputs"]["normals"][0], normal_path),
        (provenance["inputs"]["targets"], vcf),
    ):
        assert identity["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        assert identity["size_bytes"] == path.stat().st_size
        assert identity["identity_method"] == "sha256_file_bytes"
    assert provenance["reference"]["status"] == "INSUFFICIENT_METADATA"
    with pysam.VariantFile(str(output)) as annotated:
        assert next(annotated).info["SKUA_ARTIFACT_PRIOR"] == pytest.approx((0.3,))


@pytest.mark.parametrize("separate_vcf", [False, True])
def test_cached_run_identifies_exact_panel_membership_and_build_inputs(tmp_path, separate_vcf):
    vcf = targets(tmp_path / "targets.vcf")
    vcf.write_text(vcf.read_text().replace(
        "#CHROM", '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior">\n#CHROM',
    ).replace("PASS\t.", "PASS\tSKUA_ARTIFACT_PRIOR=0.8"))
    panel = tmp_path / "panel.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL") as normal:
        skua.build_pon(vcf, normal_alignments=[normal], output_path=panel, min_baseq=25, min_mapq=35)
    build = skua.read_provenance(panel)
    assert build["mode"] == "pon_build"
    assert build["normal_samples"] == ["NORMAL"]
    assert build["inputs"]["normals"][0]["sha256"] == hashlib.sha256((tmp_path / "normal.bam").read_bytes()).hexdigest()
    assert skua.read_pon_metadata(panel).provenance == build
    output = tmp_path / "out.vcf"
    with alignment(tmp_path / "case.bam", "CASE") as case:
        skua.annotate_vcf_with_pon(
            case, panel, vcf_path=vcf if separate_vcf else None, output_path=output,
            truncate=0.25, pseudocount=0.02, prior_artifact_probability=0.4,
        )
    run = skua.read_provenance(output)
    assert run["mode"] == "cached_pon"
    assert run["normal_samples"] == ["NORMAL"]
    assert run["panel_build"] == build
    assert run["inputs"]["pon"]["sha256"] == hashlib.sha256(panel.read_bytes()).hexdigest()
    assert run["evidence"]["min_baseq"] == 25
    assert run["evidence"]["min_mapq"] == 35
    assert run["model"]["truncate"] == 0.25
    assert run["model"]["pseudocount"] == 0.02
    assert run["model"]["prior"]["fallback"] == 0.4
    assert run["case"]["sample_name"] == "CASE"
    inputs, model = run["inputs"], run["model"]
    replay = tmp_path / "replayed.vcf.gz"
    with pysam.AlignmentFile(inputs["case"]["path"], "rb") as case:
        skua.annotate_vcf_with_pon(
            case, inputs["pon"]["path"], output_path=replay,
            vcf_path=None if inputs["targets"]["sha256"] == inputs["pon"]["sha256"] else inputs["targets"]["path"],
            sample_name=run["case"]["sample_name"],
            strict=run["options"]["strict"], force=run["options"]["force"],
            truncate=model["truncate"], pseudocount=model["pseudocount"],
            prior_artifact_probability=model["prior"]["fallback"],
            assessment_thresholds=skua.AssessmentThresholds(**model["assessment_thresholds"]),
        )
    with pysam.VariantFile(str(output)) as original, pysam.VariantFile(str(replay)) as repeated:
        expected, observed = next(original), next(repeated)
        assert dict(observed.info) == dict(expected.info)
        assert dict(observed.samples["CASE"]) == dict(expected.samples["CASE"])
        assert observed.info["SKUA_ARTIFACT_PRIOR"] == pytest.approx((0.8,))


def test_force_replaces_run_metadata_and_preserves_unrelated_annotations(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    vcf.write_text(vcf.read_text().replace(
        "#CHROM", '##external=preserve_me\n##INFO=<ID=KEEP,Number=1,Type=Integer,Description="Keep">\n#CHROM',
    ).replace("PASS\t.", "PASS\tKEEP=7"))
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf.gz"
    with alignment(tmp_path / 'case, "ø".bam', "CASE") as case:
        skua.annotate_vcf_with_normals(case, vcf, output_path=first, min_mapq=25)
        skua.annotate_vcf_with_normals(case, first, output_path=second, min_mapq=40, force=True)
    provenance = skua.read_provenance(second)
    assert provenance["evidence"]["min_mapq"] == 40
    assert provenance["options"] == {"force": True, "strict": False}
    assert provenance["inputs"]["case"]["path"].endswith('case, "ø".bam')
    assert provenance["inputs"]["targets"]["sha256"] == hashlib.sha256(first.read_bytes()).hexdigest()
    with pysam.VariantFile(str(second)) as result:
        assert len([r for r in result.header.records if r.key == "SKUA_PROVENANCE"]) == 1
        assert "##external=preserve_me" in str(result.header)
        assert next(result).info["KEEP"] == 7


def test_case_only_annotation_records_reference_input_identity(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 200 + "\n")
    pysam.faidx(str(reference))
    output = tmp_path / "case-only.vcf"
    with alignment(tmp_path / "case.bam", "CASE") as case:
        skua.annotate_vcf(case, vcf, output_path=output, reference_path=reference)
    provenance = skua.read_provenance(output)
    assert provenance["mode"] == "case_only"
    assert provenance["model"] is None
    assert provenance["inputs"]["reference"]["sha256"] == hashlib.sha256(reference.read_bytes()).hexdigest()


@pytest.mark.parametrize("with_normals", [False, True])
def test_json_document_retains_provenance_and_actual_python_selection_policy(tmp_path, with_normals):
    vcf = targets(tmp_path / "targets.vcf")
    output = tmp_path / "result.json"
    with alignment(tmp_path / "case.bam", "CASE") as case, alignment(tmp_path / "normal.bam", "NORMAL") as normal:
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
    assert len(document["records"]) == 1
    assert document["provenance"]["evidence"]["min_baseq"] == 24
    assert document["provenance"]["evidence"]["min_mapq"] == 32
    assert document["provenance"]["case"]["sample_name"] == "CASE"
    assert document["provenance"]["case"]["selection"] == "all_alignment_reads"
    assert document["provenance"]["case"]["read_group_ids"] is None
    if with_normals:
        assert document["provenance"]["model"]["prior"] == {"policy": "constant", "fallback": 0.3}
        assert document["records"][0]["artifact_prior"] == 0.3


def rewrite_panel_provenance(panel, output, provenance):
    serialized = bcftools.view("-Ov", str(panel))
    replacement = "" if provenance is None else "##SKUA_PROVENANCE=" + base64.b64encode(json.dumps(provenance).encode()).decode() + "\n"
    serialized = re.sub(r"^##SKUA_PROVENANCE=.*\n", replacement, serialized, flags=re.MULTILINE)
    source = output.with_suffix(".vcf")
    source.write_text(serialized)
    bcftools.view("-Ob", "-o", str(output), str(source), catch_stdout=False)
    bcftools.index(str(output))


@pytest.mark.parametrize("conflict", ["sample", "threshold", "reference", "mode", "version"])
def test_contradictory_panel_provenance_is_rejected_before_output(tmp_path, conflict):
    vcf = targets(tmp_path / "targets.vcf")
    panel, damaged = tmp_path / "panel.bcf", tmp_path / "damaged.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL") as normal:
        skua.build_pon(vcf, normal_alignments=[normal], output_path=panel)
    provenance = skua.read_provenance(panel)
    if conflict == "sample":
        provenance["normal_samples"] = ["OTHER"]
    elif conflict == "threshold":
        provenance["evidence"]["min_mapq"] = 99
    elif conflict == "reference":
        provenance["reference"]["contigs"][0]["length"] = 999
    elif conflict == "mode":
        provenance["mode"] = "live_normals"
    else:
        provenance["schema_version"] = 999
    rewrite_panel_provenance(panel, damaged, provenance)
    with pytest.raises(ValueError, match="provenance"):
        skua.read_pon_metadata(damaged)
    assert not validate_pon(damaged).valid
    output = tmp_path / "out.vcf"
    output.write_text("preserve existing output\n")
    with alignment(tmp_path / "case.bam", "CASE") as case:
        with pytest.raises(ValueError, match="provenance"):
            skua.annotate_vcf_with_pon(case, damaged, output_path=output, force=True)
    assert output.read_text() == "preserve existing output\n"


def test_legacy_build_metadata_is_explicitly_unavailable_but_panel_identity_is_exact(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    panel, legacy = tmp_path / "panel.bcf", tmp_path / "legacy.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL") as normal:
        skua.build_pon(vcf, normal_alignments=[normal], output_path=panel)
    rewrite_panel_provenance(panel, legacy, None)
    assert skua.read_pon_metadata(legacy).provenance is None
    assert inspect_pon(legacy).as_dict()["provenance"] is None
    assert validate_pon(legacy).valid
    with alignment(tmp_path / "case.bam", "CASE") as case:
        skua.annotate_vcf_with_pon(case, legacy, output_path=tmp_path / "out.vcf")
    run = skua.read_provenance(tmp_path / "out.vcf")
    assert run["panel_build"] is None
    assert run["normal_samples"] == ["NORMAL"]
    assert run["inputs"]["pon"]["sha256"] == hashlib.sha256(legacy.read_bytes()).hexdigest()


def test_python_objects_without_backing_files_do_not_claim_content_identity(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    document = json.loads(skua.annotate_vcf_to_json(FakeAlignmentFile([]), vcf))
    assert document["provenance"]["inputs"]["case"] == {
        "path": None, "identity_method": "unavailable", "sha256": None, "size_bytes": None,
    }
    assert document["provenance"]["case"]["sample_name"] is None


@pytest.mark.parametrize("metadata", [
    "##SKUA_PROVENANCE=not-base64!\n",
    "##SKUA_PROVENANCE=e30=\n",
    "##SKUA_PROVENANCE=e30=\n##SKUA_PROVENANCE=e30=\n",
])
def test_public_reader_rejects_malformed_or_duplicate_provenance(tmp_path, metadata):
    vcf = targets(tmp_path / "targets.vcf")
    vcf.write_text(vcf.read_text().replace("#CHROM", metadata + "#CHROM"))
    with pytest.raises(ValueError, match="(?i)provenance"):
        skua.read_provenance(vcf)
