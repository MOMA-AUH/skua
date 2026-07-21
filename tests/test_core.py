import json

import pytest

from skua.core import (
    PonAnnotation,
    annotate_variant,
    annotate_variant_with_normals,
    annotate_variants_from_vcf,
    annotate_vcf,
    annotate_vcf_to_json,
    annotate_vcf_with_normals,
    format_annotation_results,
    render_annotation_results_json,
    write_annotation_results_json,
)
from skua.evidence import AggregatedEvidence, UnusableReason
from tests.helpers import FakeAlignmentFile, FakeAlignmentHeader, FakeRead, build_linear_pairs
from skua.variants import Variant


def test_annotate_variant_collects_evidence_for_single_variant() -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)
    variant = Variant(contig="chr1", ref_pos0=105, ref="A", alt="T")

    counts = annotate_variant(
        alignment_file,
        variant,
        min_baseq=20,
        min_mapq=20,
    )

    assert alignment_file.fetch_calls == [("chr1", 105, 106)]
    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 1
    assert counts.usable == 2
    assert counts.unusable == 0


def test_annotate_variants_from_vcf_processes_simple_records_only(tmp_path) -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.",
                "chr1\t200\t.\tA\tAT\t.\tPASS\t.",
                "chr1\t300\t.\tAT\tA\t.\tPASS\t.",
                "chr1\t400\t.\tC\tG,T\t.\tPASS\t.",
            ]
        )
        + "\n"
    )

    results = list(
        annotate_variants_from_vcf(
            alignment_file,
            vcf_path,
            min_baseq=20,
            min_mapq=20,
        )
    )

    assert [variant for variant, _counts in results] == [
        Variant(contig="chr1", ref_pos0=105, ref="A", alt="T"),
        Variant(contig="chr1", ref_pos0=199, ref="A", alt="AT"),
        Variant(contig="chr1", ref_pos0=299, ref="AT", alt="A"),
    ]
    assert alignment_file.fetch_calls == [("chr1", 105, 106), ("chr1", 199, 200), ("chr1", 299, 300)]


def test_format_annotation_results_returns_json_ready_records() -> None:
    results = [
        (
            Variant(contig="chr1", ref_pos0=105, ref="A", alt="T"),
            AggregatedEvidence(
                alt_forward=1,
                alt_reverse=2,
                non_alt_forward=3,
                non_alt_reverse=4,
                usable=10,
                unusable=2,
                unusable_by_reason={UnusableReason.LOW_MAPQ: 2},
            ),
        )
    ]

    records = format_annotation_results(results)

    assert records == [
        {
            "contig": "chr1",
            "pos1": 106,
            "ref": "A",
            "alt": "T",
            "counts": {
                "case": {
                    "alt_forward": 1,
                    "alt_reverse": 2,
                    "non_alt_forward": 3,
                    "non_alt_reverse": 4,
                    "usable": 10,
                    "unusable": 2,
                    "unusable_by_reason": {"low_mapq": 2},
                },
            },
        }
    ]


def test_verify_and_format_from_vcf_end_to_end(tmp_path) -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=5,
            is_reverse=True,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.",
            ]
        )
        + "\n"
    )

    rows = format_annotation_results(
        annotate_variants_from_vcf(
            alignment_file,
            vcf_path,
            min_baseq=20,
            min_mapq=20,
        )
    )

    assert rows == [
        {
            "contig": "chr1",
            "pos1": 106,
            "ref": "A",
            "alt": "T",
            "counts": {
                "case": {
                    "alt_forward": 1,
                    "alt_reverse": 0,
                    "non_alt_forward": 0,
                    "non_alt_reverse": 1,
                    "usable": 2,
                    "unusable": 1,
                    "unusable_by_reason": {"low_mapq": 1},
                },
            },
        }
    ]


def test_render_annotation_results_json_returns_json_text() -> None:
    rows = [
        {
            "contig": "chr1",
            "pos1": 106,
            "ref": "A",
            "alt": "T",
            "counts": {
                "case": {
                    "alt_forward": 1,
                    "alt_reverse": 0,
                    "non_alt_forward": 0,
                    "non_alt_reverse": 1,
                    "usable": 2,
                    "unusable": 1,
                    "unusable_by_reason": {"low_mapq": 1},
                },
            },
        }
    ]

    payload = render_annotation_results_json(rows)

    assert json.loads(payload) == rows


def test_write_annotation_results_json_writes_payload_to_file(tmp_path) -> None:
    rows = [
        {
            "contig": "chr1",
            "pos1": 106,
            "ref": "A",
            "alt": "T",
            "case": {
                "alt_forward": 1,
                "alt_reverse": 0,
                "non_alt_forward": 0,
                "non_alt_reverse": 1,
                "usable": 2,
                "unusable": 1,
                "unusable_by_reason": {"low_mapq": 1},
            },
        }
    ]
    output_path = tmp_path / "verification.json"

    write_annotation_results_json(rows, output_path)

    assert json.loads(output_path.read_text(encoding="utf-8")) == rows


def test_annotate_vcf_to_json_returns_payload_and_writes_file(tmp_path) -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.",
            ]
        )
        + "\n"
    )
    output_path = tmp_path / "verification.json"

    payload = annotate_vcf_to_json(
        alignment_file,
        vcf_path,
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

    expected_rows = [
        {
            "contig": "chr1",
            "pos1": 106,
            "ref": "A",
            "alt": "T",
            "counts": {
                "case": {
                    "alt_forward": 1,
                    "alt_reverse": 0,
                    "non_alt_forward": 0,
                    "non_alt_reverse": 1,
                    "usable": 2,
                    "unusable": 0,
                    "unusable_by_reason": {},
                },
            },
        }
    ]
    assert json.loads(payload) == expected_rows
    assert json.loads(output_path.read_text(encoding="utf-8")) == expected_rows


def test_annotate_vcf_writes_case_format_fields(tmp_path) -> None:
    import pysam

    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=5,
            is_reverse=True,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"

    payload = annotate_vcf(
        alignment_file,
        vcf_path,
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

    assert "SKUA_ALT_FWD" in payload
    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        sample = record.samples["CASE"]
        assert sample["SKUA_ALT_FWD"] == 1
        assert sample["SKUA_ALT_REV"] == 0
        assert sample["SKUA_NON_ALT_FWD"] == 0
        assert sample["SKUA_NON_ALT_REV"] == 1
        assert sample["SKUA_USABLE"] == 2
        assert sample["SKUA_UNUSABLE"] == 1


def test_annotate_vcf_supports_simple_insertion(tmp_path) -> None:
    import pysam

    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="ATAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=[
                (0, 100),
                (1, None),
                (2, 101),
                (3, 102),
                (4, 103),
                (5, 104),
                (6, 105),
                (7, 106),
                (8, 107),
                (9, 108),
            ],
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t101\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated_insertion.vcf"

    payload = annotate_vcf(
        alignment_file,
        vcf_path,
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

    assert "SKUA_ALT_FWD" in payload
    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        sample = record.samples["CASE"]
        assert sample["SKUA_ALT_FWD"] == 1
        assert sample["SKUA_NON_ALT_FWD"] == 0
        assert sample["SKUA_USABLE"] == 1
        assert sample["SKUA_UNUSABLE"] == 0


def test_annotate_vcf_supports_bgzipped_output(tmp_path) -> None:
    import pysam

    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf.gz"

    payload = annotate_vcf(
        alignment_file,
        vcf_path,
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

    assert output_path.read_bytes()[:2] == b"\x1f\x8b"
    assert "#CHROM" in payload
    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        sample = record.samples["CASE"]
        assert sample["SKUA_ALT_FWD"] == 1


def test_annotate_vcf_with_normals_adds_sample_for_site_only_vcf(tmp_path) -> None:
    import pysam

    case_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    case_alignment = FakeAlignmentFile(
        case_reads,
        header=FakeAlignmentHeader([{"SM": "CASE"}]),
    )

    vcf_path = tmp_path / "site_only.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated_site_only.vcf"

    payload = annotate_vcf_with_normals(
        case_alignment,
        vcf_path,
        normal_alignments=[],
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

    assert "SKUA_ARTIFACT_POSTERIOR" in payload
    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        assert list(annotated_vcf.header.samples) == ["CASE"]
        record = next(iter(annotated_vcf))
        sample = record.samples["CASE"]
        assert sample["SKUA_ALT_FWD"] == 1
        assert sample["SKUA_ALT_REV"] == 0
        assert sample["SKUA_NON_ALT_FWD"] == 0
        assert sample["SKUA_NON_ALT_REV"] == 0
        assert sample["SKUA_USABLE"] == 1
        assert sample["SKUA_UNUSABLE"] == 0
        assert 0.0 <= sample["SKUA_ARTIFACT_POSTERIOR"] <= 1.0
        assert isinstance(sample["SKUA_LOG_BAYES_FACTOR"], float)


def test_annotate_vcf_with_normals_requires_single_alignment_sample_name(tmp_path) -> None:
    case_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    vcf_path = tmp_path / "site_only.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    no_sm_alignment = FakeAlignmentFile(case_reads, header=FakeAlignmentHeader([]))
    with pytest.raises(ValueError, match="usable read-group SM"):
        annotate_vcf_with_normals(
            no_sm_alignment,
            vcf_path,
            normal_alignments=[],
            min_baseq=20,
            min_mapq=20,
        )

    multi_sm_alignment = FakeAlignmentFile(
        case_reads,
        header=FakeAlignmentHeader([{"SM": "CASE"}, {"SM": "TUMOR"}]),
    )
    with pytest.raises(ValueError, match="multiple distinct read-group SM tags"):
        annotate_vcf_with_normals(
            multi_sm_alignment,
            vcf_path,
            normal_alignments=[],
            min_baseq=20,
            min_mapq=20,
        )


def test_annotate_vcf_with_normals_writes_info_and_format(tmp_path) -> None:
    import pysam

    case_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    case_alignment = FakeAlignmentFile(case_reads)

    normal_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=5,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    normal_alignment = FakeAlignmentFile(normal_reads)

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated_pon.vcf"

    payload = annotate_vcf_with_normals(
        case_alignment,
        vcf_path,
        normal_alignments=[normal_alignment],
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

    assert "SKUA_LOG_BAYES_FACTOR" in payload
    assert "SKUA_ARTIFACT_POSTERIOR" in payload
    
    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        sample = record.samples["CASE"]
        assert sample["SKUA_ALT_FWD"] == 1
        assert 0.0 <= sample["SKUA_ARTIFACT_POSTERIOR"] <= 1.0
        assert isinstance(sample["SKUA_LOG_BAYES_FACTOR"], float)
        assert record.info["SKUA_PON_SAMPLE_COUNT"] == 1
        assert record.info["SKUA_PON_ALT_FWD"] == 0
        assert record.info["SKUA_PON_ALT_REV"] == 0
        assert record.info["SKUA_PON_NON_ALT_FWD"] == 1
        assert record.info["SKUA_PON_NON_ALT_REV"] == 0
        assert record.info["SKUA_PON_USABLE"] == 1
        assert record.info["SKUA_PON_UNUSABLE"] == 1
        assert record.info["SKUA_PON_DISPERSION_FACTOR"] == pytest.approx(1e-4)


def test_annotate_variant_with_normals_returns_case_and_normal_evidence() -> None:
    case_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    case_alignment = FakeAlignmentFile(case_reads)

    normal1_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    normal1_alignment = FakeAlignmentFile(normal1_reads)

    normal2_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    normal2_alignment = FakeAlignmentFile(normal2_reads)

    variant = Variant(contig="chr1", ref_pos0=105, ref="A", alt="T")

    result = annotate_variant_with_normals(
        case_alignment,
        variant,
        normal_alignments=[normal1_alignment, normal2_alignment],
        min_baseq=20,
        min_mapq=20,
    )

    assert result.case_evidence.alt_forward == 1
    assert result.case_evidence.non_alt_forward == 0
    assert len(result.normal_evidences) == 2
    assert result.normal_aggregate_evidence.alt_forward == 1
    assert result.normal_aggregate_evidence.non_alt_forward == 1
    assert result.normal_aggregate_evidence.usable == 2
    assert result.normal_aggregate_evidence.unusable == 0


def test_annotate_vcf_to_json_with_normals_returns_pon_payload(tmp_path) -> None:
    from skua.core import annotate_vcf_to_json_with_normals

    case_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    case_alignment = FakeAlignmentFile(case_reads)

    normal_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    normal_alignment = FakeAlignmentFile(normal_reads)

    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.",
            ]
        )
        + "\n"
    )

    payload = annotate_vcf_to_json_with_normals(
        case_alignment,
        vcf_path,
        normal_alignments=[normal_alignment],
        min_baseq=20,
        min_mapq=20,
    )

    import json
    result = json.loads(payload)
    assert len(result) == 1
    assert result[0]["contig"] == "chr1"
    assert result[0]["pos1"] == 106
    assert result[0]["counts"]["normal"]["alt_forward"] == 0
    assert result[0]["counts"]["normal"]["alt_reverse"] == 0
    assert result[0]["counts"]["normal"]["non_alt_forward"] == 1
    assert result[0]["counts"]["normal"]["non_alt_reverse"] == 0
    assert result[0]["counts"]["normal"]["usable"] == 1
    assert result[0]["counts"]["normal"]["unusable"] == 0
    assert result[0]["counts"]["normal"]["unusable_by_reason"] == {}
    assert result[0]["counts"]["case"]["alt_forward"] == 1
    assert result[0]["counts"]["case"]["alt_reverse"] == 0
    assert result[0]["counts"]["case"]["non_alt_forward"] == 0
    assert result[0]["counts"]["case"]["non_alt_reverse"] == 0
    assert result[0]["counts"]["case"]["usable"] == 1
    assert result[0]["counts"]["case"]["unusable"] == 0
    assert result[0]["counts"]["case"]["unusable_by_reason"] == {}
    assert 0.0 <= result[0]["stats"]["artifact_posterior"] <= 1.0
    assert isinstance(result[0]["stats"]["log_bayes_factor_artifact_vs_variant"], float)
    assert result[0]["stats"]["dispersion_factor"] == 1e-4
    assert result[0]["stats"]["pon_sample_count"] == 1
    assert list(result[0].keys()) == [
        "contig",
        "pos1",
        "ref",
        "alt",
        "stats",
        "counts",
    ]
    assert list(result[0]["stats"].keys()) == [
        "artifact_posterior",
        "log_bayes_factor_artifact_vs_variant",
        "dispersion_factor",
        "pon_sample_count",
    ]


def test_format_annotation_results_with_normals_excludes_truncated_normals() -> None:
    from skua.core import format_annotation_results_with_normals

    variant = Variant(contig="chr1", ref_pos0=105, ref="A", alt="T")
    case_evidence = AggregatedEvidence(
        alt_forward=2,
        alt_reverse=0,
        non_alt_forward=8,
        non_alt_reverse=0,
        usable=10,
        unusable=0,
        unusable_by_reason={},
    )

    low_background = AggregatedEvidence(
        alt_forward=1,
        alt_reverse=0,
        non_alt_forward=99,
        non_alt_reverse=0,
        usable=100,
        unusable=0,
        unusable_by_reason={},
    )
    high_background_outlier = AggregatedEvidence(
        alt_forward=20,
        alt_reverse=0,
        non_alt_forward=80,
        non_alt_reverse=0,
        usable=100,
        unusable=0,
        unusable_by_reason={},
    )

    rows = format_annotation_results_with_normals(
        [
            (
                variant,
                PonAnnotation(
                    case_evidence=case_evidence,
                    normal_evidences=(low_background, high_background_outlier),
                    normal_aggregate_evidence=AggregatedEvidence(
                        alt_forward=21,
                        alt_reverse=0,
                        non_alt_forward=179,
                        non_alt_reverse=0,
                        usable=200,
                        unusable=0,
                        unusable_by_reason={},
                    ),
                ),
            )
        ]
    )

    assert rows[0]["stats"]["pon_sample_count"] == 1
    assert rows[0]["counts"]["normal"]["alt_forward"] == 1
    assert rows[0]["counts"]["normal"]["non_alt_forward"] == 99
    assert rows[0]["counts"]["normal"]["usable"] == 100
