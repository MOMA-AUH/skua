"""Assessment eligibility through the supported Python and CLI interfaces."""

from contextlib import ExitStack
import json

import pysam
import pytest

import skua
from skua.cli import main


def _evidence(alt_forward=0, alt_reverse=0, non_alt_forward=0, non_alt_reverse=0):
    return skua.AggregatedEvidence(
        alt_forward=alt_forward,
        alt_reverse=alt_reverse,
        non_alt_forward=non_alt_forward,
        non_alt_reverse=non_alt_reverse,
        usable=alt_forward + alt_reverse + non_alt_forward + non_alt_reverse,
        unusable=0,
        unusable_by_reason={},
    )


def test_zero_evidence_has_no_model_scores() -> None:
    stats = skua.compute_stats(
        _evidence(), _evidence(), per_sample_evidences=[],
        prior_artifact_probability=0.001,
    )

    assert stats.assessment_status == "INSUFFICIENT_EVIDENCE"
    assert stats.assessment_reasons == ("CASE_DEPTH", "NORMAL_DEPTH")
    assert stats.artifact_posterior is None
    assert stats.log_bayes_factor_artifact_vs_variant is None
    assert stats.case_counts["alt_forward"] == 0
    assert stats.normal_counts["non_alt_forward"] == 0


@pytest.mark.parametrize(
    ("overrides", "reasons"),
    [
        ({}, ()),
        ({"min_case_depth": 5}, ("CASE_DEPTH",)),
        ({"min_normal_depth": 9}, ("NORMAL_DEPTH",)),
        ({"min_normal_samples": 3}, ("NORMAL_SAMPLE_COUNT",)),
        ({"min_case_strand_depth": 3}, ("CASE_STRAND_DEPTH",)),
        ({"min_normal_strand_depth": 5}, ("NORMAL_STRAND_DEPTH",)),
    ],
)
def test_assessment_thresholds_gate_scores_and_preserve_counts(overrides, reasons) -> None:
    limits = dict(
        min_case_depth=4, min_normal_depth=8, min_normal_samples=2,
        min_case_strand_depth=2, min_normal_strand_depth=4,
    )
    limits.update(overrides)
    case = _evidence(1, 1, 1, 1)
    normal = _evidence(0, 0, 4, 4)
    samples = [_evidence(0, 0, 2, 2), _evidence(0, 0, 2, 2)]
    baseline = skua.compute_stats(case, normal, per_sample_evidences=samples)
    stats = skua.compute_stats(
        case, normal, per_sample_evidences=samples,
        assessment_thresholds=skua.AssessmentThresholds(**limits),
    )
    assert stats.assessment_status == ("INSUFFICIENT_EVIDENCE" if reasons else "ASSESSED")
    assert stats.assessment_reasons == reasons
    assert stats.case_counts == baseline.case_counts
    assert stats.normal_counts == baseline.normal_counts
    if reasons:
        assert stats.artifact_posterior is None
        assert stats.log_bayes_factor_artifact_vs_variant is None
    else:
        assert stats.artifact_posterior == baseline.artifact_posterior
        assert stats.log_bayes_factor_artifact_vs_variant == baseline.log_bayes_factor_artifact_vs_variant


@pytest.mark.parametrize(("rho", "expected"), [(1e-8, 1e-6), (1 - 1e-8, 1 - 1e-6)])
def test_assessment_thresholds_preserve_bounded_dispersion(rho, expected) -> None:
    case, normal = _evidence(1, 1, 1, 1), _evidence(0, 0, 10, 10)
    assessed = skua.compute_stats(case, normal, rho=rho)
    insufficient = skua.compute_stats(
        case, normal, rho=rho,
        assessment_thresholds=skua.AssessmentThresholds(min_case_depth=5),
    )
    assert assessed.assessment_status == "ASSESSED"
    assert insufficient.assessment_status == "INSUFFICIENT_EVIDENCE"
    assert assessed.dispersion_rho == insufficient.dispersion_rho == expected
    assert insufficient.artifact_posterior is None
    assert insufficient.log_bayes_factor_artifact_vs_variant is None


@pytest.mark.parametrize("name", [
    "min_case_depth", "min_normal_depth", "min_normal_samples",
    "min_case_strand_depth", "min_normal_strand_depth",
])
@pytest.mark.parametrize("value", [-1, 1.5, float("nan"), float("inf"), True])
def test_invalid_assessment_thresholds_are_rejected(name, value) -> None:
    with pytest.raises(ValueError, match=name):
        skua.compute_stats(
            _evidence(), _evidence(),
            assessment_thresholds=skua.AssessmentThresholds(**{name: value}),
        )


@pytest.mark.parametrize("name", ["min_case_depth", "min_normal_depth"])
def test_zero_evidence_cannot_be_enabled_by_lowering_thresholds(name) -> None:
    with pytest.raises(ValueError, match=name):
        skua.compute_stats(
            _evidence(), _evidence(),
            assessment_thresholds=skua.AssessmentThresholds(**{name: 0}),
        )


def _write_bam(path, sample, counts) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "RG": [{"ID": "rg1", "SM": sample}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for channel, count in enumerate((*counts, 1)):
            for index in range(count):
                read = pysam.AlignedSegment()
                read.query_name = f"{channel}-{index}"
                read.query_sequence = "AAAAATAAAA" if channel < 2 else "AAAAAAAAAA"
                read.flag = 147 if channel % 2 else 99
                read.reference_id = 0
                read.reference_start = 100
                read.mapping_quality = 0 if channel == 4 else 60
                read.cigartuples = [(0, 10)]
                read.next_reference_id = 0
                read.next_reference_start = 100
                read.query_qualities = pysam.qualitystring_to_array("I" * 10)
                read.set_tag("RG", "rg1")
                bam.write(read)
    pysam.index(str(path))


def _write_inputs(tmp_path, case, normals):
    targets = tmp_path / "targets.vcf"
    targets.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t106\t.\tA\tT\t.\tPASS\t.\n",
    )
    case_path = tmp_path / "case.bam"
    _write_bam(case_path, "CASE", case)
    normal_paths = []
    for index, counts in enumerate(normals):
        path = tmp_path / f"normal{index}.bam"
        _write_bam(path, f"NORMAL{index}", counts)
        normal_paths.append(path)
    normal_list = tmp_path / "normals.lst"
    normal_list.write_text("".join(f"{path}\n" for path in normal_paths))
    return targets, case_path, normal_paths, normal_list


@pytest.mark.parametrize("cached", [False, True], ids=["direct", "cached"])
@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_reasons_format_is_present_on_all_model_annotated_records(tmp_path, cached, suffix):
    targets, case_path, normal_paths, _ = _write_inputs(
        tmp_path, (1, 1, 1, 1), [(0, 0, 2, 2)],
    )
    targets.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE\tOTHER\n"
        "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0\n"
        "chr1\t206\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0\n",
    )
    output = tmp_path / f"calls{suffix}"
    with (
        pysam.AlignmentFile(str(case_path)) as case,
        pysam.AlignmentFile(str(normal_paths[0])) as normal,
    ):
        if cached:
            panel = tmp_path / "panel.bcf"
            skua.build_pon(targets, normal_alignments=[normal], output_path=panel)
            skua.annotate_vcf_with_pon(case, panel, vcf_path=targets, output_path=output)
        else:
            skua.annotate_vcf_with_normals(
                case, targets, normal_alignments=[normal], output_path=output,
            )

    with pysam.VariantFile(str(output)) as vcf:
        assert "SKUA_ASSESSMENT_REASONS" in vcf.header.formats
        assessed, insufficient = list(vcf)
        assert assessed.samples["CASE"]["SKUA_ASSESSMENT_STATUS"] == "ASSESSED"
        assert assessed.samples["CASE"]["SKUA_ASSESSMENT_REASONS"] == (".",)
        assert tuple(assessed.format) == tuple(insufficient.format)
        assert insufficient.samples["CASE"]["SKUA_ASSESSMENT_STATUS"] == "INSUFFICIENT_EVIDENCE"
        assert insufficient.samples["CASE"]["SKUA_ASSESSMENT_REASONS"] == ("CASE_DEPTH", "NORMAL_DEPTH")
        assert insufficient.samples["OTHER"]["SKUA_ASSESSMENT_REASONS"] == (".",)
        for field in ("SKUA_LOG_BAYES_FACTOR", "SKUA_ARTIFACT_POSTERIOR"):
            assert isinstance(assessed.samples["CASE"][field], float)
            assert insufficient.samples["CASE"][field] is None
            assert assessed.samples["OTHER"][field] is None
            assert insufficient.samples["OTHER"][field] is None
        for record in (assessed, insufficient):
            assert record.samples["OTHER"]["SKUA_ASSESSMENT_STATUS"] == "."
            assert record.samples["OTHER"]["SKUA_ASSESSMENT_REASONS"] == (".",)
            assert record.samples["CASE"]["GT"] == (0, 1)
            assert record.samples["OTHER"]["GT"] == (0, 0)


@pytest.mark.parametrize("mode", ["direct", "pon", "vcf-pon"])
@pytest.mark.parametrize(
    ("case", "normals", "thresholds", "reasons", "case_depth", "normal_depth", "normal_samples"),
    [
        pytest.param((0, 0, 0, 0), [(0, 0, 0, 0)], {},
                     ("CASE_DEPTH", "NORMAL_DEPTH"), 0, 0, 0, id="zero-evidence"),
        pytest.param((0, 0, 0, 0), [(0, 0, 2, 2)], {},
                     ("CASE_DEPTH",), 0, 4, 1, id="zero-case"),
        pytest.param((1, 1, 1, 1), [(0, 0, 0, 0)], {},
                     ("NORMAL_DEPTH",), 4, 0, 0, id="zero-normal"),
        pytest.param((1, 1, 1, 1), [(2, 2, 0, 0)], {},
                     ("NORMAL_DEPTH",), 4, 0, 0, id="all-truncated"),
        pytest.param((2, 0, 0, 0), [(0, 0, 4, 0)], {},
                     (), 2, 4, 1, id="one-strand-default"),
        pytest.param((2, 0, 0, 0), [(0, 0, 4, 0)],
                     {"min_case_strand_depth": 1, "min_normal_strand_depth": 1},
                     ("CASE_STRAND_DEPTH", "NORMAL_STRAND_DEPTH"), 2, 4, 1, id="one-strand-required"),
        pytest.param((1, 1, 1, 1), [(0, 0, 2, 2), (0, 0, 2, 2)],
                     dict(min_case_depth=4, min_normal_depth=8, min_normal_samples=2,
                          min_case_strand_depth=2, min_normal_strand_depth=4),
                     (), 4, 8, 2, id="exact-boundaries"),
        pytest.param((1, 1, 1, 1), [(0, 0, 2, 2), (2, 2, 0, 0)],
                     dict(min_case_depth=5, min_normal_depth=5, min_normal_samples=2,
                          min_case_strand_depth=3, min_normal_strand_depth=3),
                     ("CASE_DEPTH", "NORMAL_DEPTH", "NORMAL_SAMPLE_COUNT",
                      "CASE_STRAND_DEPTH", "NORMAL_STRAND_DEPTH"),
                     4, 4, 1, id="below-all-boundaries-after-truncation"),
    ],
)
def test_cli_and_json_assessment_preserve_counts_and_gate_scores(
    tmp_path, mode, case, normals, thresholds, reasons, case_depth, normal_depth, normal_samples,
) -> None:
    targets, case_path, normal_paths, normal_list = _write_inputs(tmp_path, case, normals)
    output = tmp_path / "calls.vcf"
    args = [
        "annotate", "--alignment", str(case_path), "--output", str(output),
        "--prior-artifact-probability", "0.001",
    ]
    if mode == "direct":
        args += ["--vcf", str(targets), "--normal-list", str(normal_list)]
    else:
        panel = tmp_path / "panel.bcf"
        assert main([
            "pon", "build", "--vcf", str(targets), "--normal-list", str(normal_list),
            "--output", str(panel),
        ]) == 0
        args += ["--pon", str(panel)]
        if mode == "vcf-pon":
            args += ["--vcf", str(targets)]
    for name, value in thresholds.items():
        args += [f"--{name.replace('_', '-')}", str(value)]
    assert main(args) == 0
    with pysam.VariantFile(str(output)) as vcf:
        record = next(vcf)
        sample = record.samples["CASE"]
        assert record.info["SKUA_STATUS"] == "ANNOTATED"
        assert sample["SKUA_ASSESSMENT_STATUS"] == ("INSUFFICIENT_EVIDENCE" if reasons else "ASSESSED")
        assert "SKUA_ASSESSMENT_REASONS" in vcf.header.formats
        if reasons:
            assert sample["SKUA_ASSESSMENT_REASONS"] == reasons
        else:
            assert sample["SKUA_ASSESSMENT_REASONS"] == (".",)
        assert sample["SKUA_USABLE"] == case_depth
        assert sample["SKUA_UNUSABLE"] == 1
        assert record.info["SKUA_PON_SAMPLE_COUNT"] == normal_samples
        assert record.info["SKUA_PON_USABLE"] == normal_depth
        assert record.info["SKUA_PON_UNUSABLE"] == normal_samples
        assert record.info["SKUA_ARTIFACT_PRIOR"] == pytest.approx((0.001,))
        if reasons:
            assert sample["SKUA_ARTIFACT_POSTERIOR"] is None
            assert sample["SKUA_LOG_BAYES_FACTOR"] is None
            fields = str(record).strip().split("\t")
            serialized_sample = dict(zip(fields[8].split(":"), fields[9].split(":"), strict=True))
            assert serialized_sample["SKUA_ARTIFACT_POSTERIOR"] == "."
            assert serialized_sample["SKUA_LOG_BAYES_FACTOR"] == "."
        else:
            assert isinstance(sample["SKUA_ARTIFACT_POSTERIOR"], float)
            assert isinstance(sample["SKUA_LOG_BAYES_FACTOR"], float)

    with ExitStack() as stack:
        alignment = stack.enter_context(pysam.AlignmentFile(str(case_path)))
        normal_alignments = [
            stack.enter_context(pysam.AlignmentFile(str(path))) for path in normal_paths
        ]
        [row] = json.loads(skua.annotate_vcf_to_json_with_normals(
            alignment, targets, normal_alignments=normal_alignments,
            prior_artifact_probability=0.001,
            assessment_thresholds=skua.AssessmentThresholds(**thresholds),
        ))["records"]
    assert row["stats"]["assessment_status"] == ("INSUFFICIENT_EVIDENCE" if reasons else "ASSESSED")
    assert row["stats"]["assessment_reasons"] == list(reasons)
    assert row["counts"]["case"]["usable"] == case_depth
    assert row["counts"]["case"]["unusable"] == 1
    assert row["counts"]["normal"]["usable"] == normal_depth
    assert row["stats"]["pon_sample_count"] == normal_samples
    assert row["artifact_prior"] == 0.001
    if reasons:
        assert row["stats"]["artifact_posterior"] is None
        assert row["stats"]["log_bayes_factor_artifact_vs_variant"] is None
    else:
        # Text VCF serialization retains fewer significant digits than JSON.
        assert row["stats"]["artifact_posterior"] == pytest.approx(
            sample["SKUA_ARTIFACT_POSTERIOR"], rel=1e-5,
        )
        assert row["stats"]["log_bayes_factor_artifact_vs_variant"] == pytest.approx(
            sample["SKUA_LOG_BAYES_FACTOR"], rel=1e-5,
        )


@pytest.mark.parametrize("samples", [[], [_evidence(2, 2)]])
def test_stats_checks_retained_normals_even_when_supplied_aggregate_has_depth(samples) -> None:
    stats = skua.compute_stats(
        _evidence(1, 1), _evidence(0, 0, 20, 20), per_sample_evidences=samples,
        assessment_thresholds=skua.AssessmentThresholds(min_normal_samples=1),
    )
    assert stats.assessment_status == "INSUFFICIENT_EVIDENCE"
    assert stats.assessment_reasons == ("NORMAL_DEPTH", "NORMAL_SAMPLE_COUNT")
    assert stats.normal_counts == {
        "alt_forward": 0, "alt_reverse": 0, "non_alt_forward": 0, "non_alt_reverse": 0,
    }
    assert stats.artifact_posterior is None
    assert stats.log_bayes_factor_artifact_vs_variant is None


def test_aggregate_only_stats_cannot_claim_to_meet_a_sample_count_requirement() -> None:
    case, normal = _evidence(1, 1), _evidence(0, 0, 20, 20)
    assert skua.compute_stats(case, normal).assessment_status == "ASSESSED"
    stats = skua.compute_stats(
        case, normal, assessment_thresholds=skua.AssessmentThresholds(min_normal_samples=2),
    )
    assert stats.assessment_status == "INSUFFICIENT_EVIDENCE"
    assert stats.assessment_reasons == ("NORMAL_SAMPLE_COUNT_UNAVAILABLE",)
    assert stats.artifact_posterior is None
    assert stats.log_bayes_factor_artifact_vs_variant is None


@pytest.mark.parametrize("source", ["--normal-list", "--pon"])
@pytest.mark.parametrize(("option", "value"), [
    ("min-case-depth", "0"), ("min-normal-depth", "0"), ("min-normal-samples", "-1"),
    ("min-case-strand-depth", "-1"), ("min-normal-strand-depth", "-1"),
    ("min-case-depth", "1.5"),
])
def test_invalid_cli_thresholds_preserve_existing_output_and_indexes(tmp_path, capsys, source, option, value):
    output = tmp_path / "calls.vcf.gz"
    protected = [output, tmp_path / "calls.vcf.gz.tbi", tmp_path / "calls.vcf.gz.csi"]
    for path in protected:
        path.write_bytes(b"existing output or index")
    with pytest.raises(SystemExit) as error:
        main([
            "annotate", "--vcf", "not-opened.vcf", "--alignment", "not-opened.bam",
            source, "not-opened-panel", "--output", str(output), "--force",
            f"--{option}={value}",
        ])
    assert error.value.code == 2
    assert f"--{option}" in capsys.readouterr().err
    assert all(path.read_bytes() == b"existing output or index" for path in protected)


@pytest.mark.parametrize("cached", [False, True], ids=["direct", "vcf-pon"])
@pytest.mark.parametrize("min_case_depth", [1, 5])
def test_force_replaces_assessment_and_scores_without_assessing_unsupported_or_other_samples(
    tmp_path, cached, min_case_depth,
):
    targets, case_path, _, normal_list = _write_inputs(
        tmp_path, (1, 1, 1, 1), [(0, 0, 2, 2)],
    )
    source_args = ["--normal-list", str(normal_list)]
    if cached:
        panel = tmp_path / "panel.bcf"
        assert main([
            "pon", "build", "--vcf", str(targets), "--normal-list", str(normal_list),
            "--output", str(panel),
        ]) == 0
        source_args = ["--pon", str(panel)]
    candidates = tmp_path / "candidates.vcf"
    candidates.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=SKUA_ASSESSMENT_STATUS,Number=1,Type=String,Description="Old">\n'
        '##FORMAT=<ID=SKUA_ASSESSMENT_REASONS,Number=.,Type=String,Description="Old">\n'
        '##FORMAT=<ID=SKUA_ARTIFACT_POSTERIOR,Number=1,Type=Float,Description="Old">\n'
        '##FORMAT=<ID=SKUA_LOG_BAYES_FACTOR,Number=1,Type=Float,Description="Old">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE\tOTHER\n"
        "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT:SKUA_ASSESSMENT_STATUS:SKUA_ASSESSMENT_REASONS"
        ":SKUA_ARTIFACT_POSTERIOR:SKUA_LOG_BAYES_FACTOR"
        "\t0/1:INSUFFICIENT_EVIDENCE:CASE_DEPTH:0.001:-8\t0/0:ASSESSED:.:0.2:-1\n"
        "chr1\t206\t.\tA\t<DEL>\t.\tPASS\t.\tGT:SKUA_ASSESSMENT_STATUS:SKUA_ASSESSMENT_REASONS"
        ":SKUA_ARTIFACT_POSTERIOR:SKUA_LOG_BAYES_FACTOR"
        "\t0/1:ASSESSED:.:0.001:-8\t0/0:ASSESSED:.:0.2:-1\n",
    )
    output = tmp_path / "calls.vcf"
    assert main([
        "annotate", "--vcf", str(candidates), "--alignment", str(case_path),
        "--output", str(output), "--force", "--min-case-depth", str(min_case_depth), *source_args,
    ]) == 0
    with pysam.VariantFile(str(output)) as vcf:
        supported, unsupported = list(vcf)
        sample = supported.samples["CASE"]
        assert sample["SKUA_USABLE"] == 4
        if min_case_depth == 1:
            assert sample["SKUA_ASSESSMENT_STATUS"] == "ASSESSED"
            assert sample["SKUA_ASSESSMENT_REASONS"] == (".",)
            assert isinstance(sample["SKUA_ARTIFACT_POSTERIOR"], float)
            assert isinstance(sample["SKUA_LOG_BAYES_FACTOR"], float)
        else:
            assert sample["SKUA_ASSESSMENT_STATUS"] == "INSUFFICIENT_EVIDENCE"
            assert sample["SKUA_ASSESSMENT_REASONS"] == ("CASE_DEPTH",)
            assert sample["SKUA_ARTIFACT_POSTERIOR"] is None
            assert sample["SKUA_LOG_BAYES_FACTOR"] is None
        assert supported.samples["OTHER"]["SKUA_ARTIFACT_POSTERIOR"] is None
        assert supported.samples["OTHER"]["SKUA_LOG_BAYES_FACTOR"] is None
        assert supported.samples["OTHER"].get("SKUA_ASSESSMENT_STATUS") in (None, ".")
        assert supported.samples["OTHER"]["SKUA_ASSESSMENT_REASONS"] == (".",)
        assert unsupported.info["SKUA_STATUS"] == "UNSUPPORTED_SYMBOLIC_ALLELE"
        assert "SKUA_ASSESSMENT_STATUS" not in unsupported.format
        assert "SKUA_ASSESSMENT_REASONS" not in unsupported.format
        assert "SKUA_ARTIFACT_POSTERIOR" not in unsupported.format
        assert "SKUA_LOG_BAYES_FACTOR" not in unsupported.format
        assert unsupported.samples["CASE"]["GT"] == (0, 1)
        assert unsupported.samples["OTHER"]["GT"] == (0, 0)
