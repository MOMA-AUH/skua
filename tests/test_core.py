import json

from skua.core import (
    format_verification_results,
    render_verification_results_json,
    verify_snv_vcf_to_json,
    verify_snv_variant,
    verify_snv_variant_with_normals,
    verify_snv_variants_from_vcf,
    write_verification_results_json,
)
from skua.evidence import AggregatedEvidence, UnusableReason
from tests.helpers import FakeAlignmentFile, FakeRead, build_linear_pairs
from skua.variants import Variant


def test_verify_snv_variant_collects_evidence_for_single_variant() -> None:
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

    counts = verify_snv_variant(
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


def test_verify_snv_variants_from_vcf_processes_snv_records_only(tmp_path) -> None:
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
                "chr1\t300\t.\tC\tG,T\t.\tPASS\t.",
            ]
        )
        + "\n"
    )

    results = list(
        verify_snv_variants_from_vcf(
            alignment_file,
            vcf_path,
            min_baseq=20,
            min_mapq=20,
        )
    )

    assert len(results) == 1
    variant, counts = results[0]
    assert variant == Variant(contig="chr1", ref_pos0=105, ref="A", alt="T")
    assert alignment_file.fetch_calls == [("chr1", 105, 106)]
    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 1
    assert counts.usable == 2
    assert counts.unusable == 0


def test_format_verification_results_returns_json_ready_records() -> None:
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

    records = format_verification_results(results)

    assert records == [
        {
            "contig": "chr1",
            "pos1": 106,
            "ref": "A",
            "alt": "T",
            "case": {
                "alt_forward": 1,
                "alt_reverse": 2,
                "non_alt_forward": 3,
                "non_alt_reverse": 4,
                "usable": 10,
                "unusable": 2,
                "unusable_by_reason": {"low_mapq": 2},
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
                "chr1\t200\t.\tA\tAT\t.\tPASS\t.",
            ]
        )
        + "\n"
    )

    rows = format_verification_results(
        verify_snv_variants_from_vcf(
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


def test_render_verification_results_json_returns_json_text() -> None:
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

    payload = render_verification_results_json(rows)

    assert json.loads(payload) == rows


def test_write_verification_results_json_writes_payload_to_file(tmp_path) -> None:
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

    write_verification_results_json(rows, output_path)

    assert json.loads(output_path.read_text(encoding="utf-8")) == rows


def test_verify_snv_vcf_to_json_returns_payload_and_writes_file(tmp_path) -> None:
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
            ]
        )
        + "\n"
    )
    output_path = tmp_path / "verification.json"

    payload = verify_snv_vcf_to_json(
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
            "case": {
                "alt_forward": 1,
                "alt_reverse": 0,
                "non_alt_forward": 0,
                "non_alt_reverse": 1,
                "usable": 2,
                "unusable": 0,
                "unusable_by_reason": {},
            },
        }
    ]
    assert json.loads(payload) == expected_rows
    assert json.loads(output_path.read_text(encoding="utf-8")) == expected_rows


def test_verify_snv_variant_with_normals_returns_case_and_normal_evidence() -> None:
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

    result = verify_snv_variant_with_normals(
        case_alignment,
        variant,
        normal_alignments=[normal1_alignment, normal2_alignment],
        min_baseq=20,
        min_mapq=20,
    )

    assert result["case_evidence"].alt_forward == 1
    assert result["case_evidence"].non_alt_forward == 0
    assert len(result["normal_evidences"]) == 2
    assert result["normal_aggregate_evidence"].alt_forward == 1
    assert result["normal_aggregate_evidence"].non_alt_forward == 1
    assert result["normal_aggregate_evidence"].usable == 2
    assert result["normal_aggregate_evidence"].unusable == 0
    assert "normals_with_alt" not in result
    assert "normals_with_ref_only" not in result


def test_verify_snv_vcf_to_json_with_normals_returns_pon_payload(tmp_path) -> None:
    from skua.core import verify_snv_vcf_to_json_with_normals

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

    payload = verify_snv_vcf_to_json_with_normals(
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
    assert result[0]["normal"]["alt_forward"] == 0
    assert result[0]["normal"]["alt_reverse"] == 0
    assert result[0]["normal"]["non_alt_forward"] == 1
    assert result[0]["normal"]["non_alt_reverse"] == 0
    assert result[0]["normal"]["usable"] == 1
    assert result[0]["normal"]["unusable"] == 0
    assert result[0]["normal"]["unusable_by_reason"] == {}
    assert result[0]["case"]["alt_forward"] == 1
    assert result[0]["case"]["alt_reverse"] == 0
    assert result[0]["case"]["non_alt_forward"] == 0
    assert result[0]["case"]["non_alt_reverse"] == 0
    assert result[0]["case"]["usable"] == 1
    assert result[0]["case"]["unusable"] == 0
    assert result[0]["case"]["unusable_by_reason"] == {}
    assert result[0]["bayes_factor"] >= 0.0
    assert 0.0 <= result[0]["artifact_posterior"] <= 1.0
    assert set(
        key
        for key in result[0].keys()
        if key in {"bayes_factor", "artifact_posterior"}
    ) == {
        "bayes_factor",
        "artifact_posterior",
    }
    assert "stats" not in result[0]


