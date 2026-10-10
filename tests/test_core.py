import json
import skua

import pytest

from skua.core import (
    ANNOTATION_STATUS_INFO_FIELD_DEFINITION,
    MODEL_SCORE_FORMAT_FIELD_DEFINITIONS,
    PON_INFO_FIELD_DEFINITIONS,
    READ_COUNT_FORMAT_FIELD_DEFINITIONS,
    AnnotationStatus,
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
from skua.pon import PON_EVIDENCE_FORMAT_FIELDS
from tests.helpers import FakeAlignmentFile, FakeAlignmentHeader, FakeRead, build_linear_pairs
from skua.variants import Variant


@pytest.mark.parametrize(
    "api",
    ["annotate_vcf_with_normals", "annotate_vcf_with_pon", "annotate_vcf_to_json_with_normals"],
)
@pytest.mark.parametrize("pseudocount", [float("nan"), float("inf"), -float("inf"), 0, -1])
def test_annotation_apis_reject_invalid_model_before_io(tmp_path, api, pseudocount) -> None:
    output_path = tmp_path / "existing-output"
    output_path.write_bytes(b"keep existing output")
    options = {} if api == "annotate_vcf_to_json_with_normals" else {"force": True}
    with pytest.raises(ValueError, match="pseudocount must be finite and > 0"):
        getattr(skua, api)(
            FakeAlignmentFile([]), tmp_path / "not-opened-input",
            output_path=output_path, pseudocount=pseudocount, **options,
        )
    assert output_path.read_bytes() == b"keep existing output"


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
    assert alignment_file.fetch_calls == [("chr1", 105, 300)]


def test_annotate_variants_from_vcf_keeps_sparse_variants_as_site_fetches(tmp_path) -> None:
    alignment_file = FakeAlignmentFile([])
    vcf_path = tmp_path / "sparse.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.",
                "chr1\t1001\t.\tA\tC\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    list(annotate_variants_from_vcf(alignment_file, vcf_path))

    assert alignment_file.fetch_calls == [
        ("chr1", 105, 106),
        ("chr1", 1000, 1001),
    ]


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
    assert json.loads(payload)["records"] == expected_rows
    assert json.loads(output_path.read_text(encoding="utf-8"))["records"] == expected_rows


def test_annotate_vcf_writes_case_format_fields(tmp_path) -> None:
    import pysam

    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "case-rg"},
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "case-rg"},
        ),
        FakeRead(
            mapping_quality=5,
            is_reverse=True,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "case-rg"},
        ),
    ]
    alignment_file = FakeAlignmentFile(
        reads,
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )

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

    assert annotate_vcf(
        alignment_file,
        vcf_path,
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    ) is None
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
            reference_start=100,
            cigartuples=((0, 1), (1, 1), (0, 8)),
            tags={"RG": "case-rg"},
        ),
    ]
    alignment_file = FakeAlignmentFile(
        reads,
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )

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

    annotate_vcf(
        alignment_file,
        vcf_path,
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        sample = record.samples["CASE"]
        assert sample["SKUA_ALT_FWD"] == 1
        assert sample["SKUA_NON_ALT_FWD"] == 0
        assert sample["SKUA_USABLE"] == 1
        assert sample["SKUA_UNUSABLE"] == 0


def test_annotate_vcf_supports_uppercase_bgzipped_output(tmp_path) -> None:
    import pysam

    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "case-rg"},
        ),
    ]
    alignment_file = FakeAlignmentFile(
        reads,
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )

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
    output_path = tmp_path / "annotated.VCF.GZ"

    annotate_vcf(
        alignment_file,
        vcf_path,
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

    assert output_path.read_bytes()[:2] == b"\x1f\x8b"
    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        sample = record.samples["CASE"]
        assert sample["SKUA_ALT_FWD"] == 1


def test_annotate_vcf_does_not_clobber_an_existing_output_by_default(tmp_path) -> None:
    alignment_file = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"
    output_path.write_text("do not replace\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="already exists"):
        annotate_vcf(alignment_file, vcf_path, output_path=output_path)

    assert output_path.read_text(encoding="utf-8") == "do not replace\n"


def test_annotate_vcf_force_atomically_replaces_an_existing_output(tmp_path) -> None:
    import pysam

    alignment_file = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"
    output_path.write_text("replace me\n", encoding="utf-8")

    annotate_vcf(alignment_file, vcf_path, output_path=output_path, force=True)

    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        assert next(iter(annotated_vcf)).info["SKUA_STATUS"] == "ANNOTATED"


def test_annotate_vcf_late_failure_does_not_publish_a_partial_output(
    monkeypatch,
    tmp_path,
) -> None:
    import skua.core as core

    alignment_file = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
                "chr1\t107\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"
    original = core.annotate_variants

    def stop_after_one(*args, **kwargs):
        yield next(original(*args, **kwargs))

    monkeypatch.setattr(core, "annotate_variants", stop_after_one)

    with pytest.raises(RuntimeError, match="ended before"):
        annotate_vcf(alignment_file, vcf_path, output_path=output_path)

    assert not output_path.exists()
    assert list(tmp_path.glob(".*annotated.vcf*")) == []


@pytest.mark.parametrize(
    "definition",
    [
        f'##FORMAT=<ID={field_id},Number=1,Type=Integer,Description="Old">'
        for field_id, _description in READ_COUNT_FORMAT_FIELD_DEFINITIONS
    ]
    + [
        f'##FORMAT=<ID={field_id},Number=1,Type={field_type},Description="Old">'
        for field_id, field_type, _description in MODEL_SCORE_FORMAT_FIELD_DEFINITIONS
    ]
    + [
        f'##INFO=<ID={field_id},Number=1,Type={field_type},Description="Old">'
        for field_id, field_type, _description in PON_INFO_FIELD_DEFINITIONS
    ]
    + [
        '##INFO=<ID={},Number=1,Type={},Description="Old">'.format(
            ANNOTATION_STATUS_INFO_FIELD_DEFINITION[0],
            ANNOTATION_STATUS_INFO_FIELD_DEFINITION[1],
        )
    ]
    + [
        f'##FORMAT=<ID={field_id},Number=1,Type=Integer,Description="Old">'
        for field_id, _description in PON_EVIDENCE_FORMAT_FIELDS
    ]
    + [
        "##SKUA_PON=<SchemaVersion=1,EvidencePolicyVersion=2,MinBaseQ=20,"
        'MinMapQ=20,SkuaVersion="0.7.0">'
    ]
    + [
        '##INFO=<ID=SKUA_UNUSED,Number=1,Type=Integer,Description="Unused">',
        '##FORMAT=<ID=SKUA_UNUSED_FMT,Number=1,Type=Integer,Description="Unused">',
    ],
)
def test_annotate_vcf_rejects_preannotated_input(definition, tmp_path) -> None:
    alignment_file = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                definition,
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="already contains Skua annotations"):
        annotate_vcf(alignment_file, vcf_path, output_path=tmp_path / "output.vcf")


def test_annotate_vcf_force_replaces_stale_annotations_on_every_record_and_sample(
    tmp_path,
) -> None:
    import pysam

    from skua.core import _ensure_skua_vcf_header_fields

    input_path = tmp_path / "preannotated.vcf"
    output_path = tmp_path / "reannotated.vcf"
    header = pysam.VariantHeader()
    header.contigs.add("chr1")
    header.add_line(
        '##INFO=<ID=CALLER_SCORE,Number=1,Type=Integer,Description="Caller score">'
    )
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    _ensure_skua_vcf_header_fields(header, include_pon_info=True)
    header.add_line(
        '##INFO=<ID=SKUA_UNUSED,Number=1,Type=Integer,Description="Unused">'
    )
    header.add_line(
        '##FORMAT=<ID=SKUA_UNUSED_FMT,Number=1,Type=Integer,Description="Unused">'
    )
    header.add_sample("CASE")
    header.add_sample("CONTROL")

    with pysam.VariantFile(str(input_path), "w", header=header) as preannotated:
        for position, alts in ((106, ("T",)), (107, ("C", "G"))):
            record = preannotated.new_record(
                contig="chr1",
                start=position - 1,
                stop=position,
                alleles=("A", *alts),
            )
            record.info["CALLER_SCORE"] = position
            record.info["SKUA_STATUS"] = "ANNOTATED"
            record.info["SKUA_UNUSED"] = 99
            record.info["SKUA_ARTIFACT_PRIOR"] = tuple(0.2 for _alt in alts)
            for field_id, field_type, _description in PON_INFO_FIELD_DEFINITIONS:
                record.info[field_id] = 9.0 if field_type == "Float" else 99
            for sample_name in ("CASE", "CONTROL"):
                sample = record.samples[sample_name]
                sample["GT"] = (0, 1)
                sample["SKUA_UNUSED_FMT"] = 99
                for field_id, _description in READ_COUNT_FORMAT_FIELD_DEFINITIONS:
                    sample[field_id] = 99
                for field_id, _field_type, _description in (
                    MODEL_SCORE_FORMAT_FIELD_DEFINITIONS
                ):
                    sample[field_id] = 0.99
            preannotated.write(record)

    alignment = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    annotate_vcf_with_normals(
        alignment,
        input_path,
        normal_alignments=[],
        output_path=output_path,
        force=True,
    )

    with pysam.VariantFile(str(output_path)) as reannotated:
        assert "SKUA_UNUSED" not in reannotated.header.info
        assert "SKUA_UNUSED_FMT" not in reannotated.header.formats
        records = list(reannotated)

    supported, unsupported = records
    assert supported.info["SKUA_STATUS"] == "ANNOTATED"
    assert supported.info["SKUA_ARTIFACT_PRIOR"][0] == pytest.approx(0.2)
    assert supported.info["CALLER_SCORE"] == 106
    assert supported.samples["CASE"]["SKUA_ALT_FWD"] == 0
    assert supported.samples["CASE"]["SKUA_USABLE"] == 0
    assert supported.samples["CONTROL"]["SKUA_ALT_FWD"] is None
    assert supported.samples["CONTROL"]["SKUA_ARTIFACT_POSTERIOR"] is None

    assert unsupported.info["SKUA_STATUS"] == "UNSUPPORTED_MULTIALLELIC"
    assert unsupported.info["SKUA_ARTIFACT_PRIOR"] == pytest.approx((0.2, 0.2))
    assert unsupported.info["CALLER_SCORE"] == 107
    assert all(
        field_id not in unsupported.info
        for field_id, *_rest in PON_INFO_FIELD_DEFINITIONS
    )
    assert all(
        field_id not in unsupported.format
        for field_id, _description in READ_COUNT_FORMAT_FIELD_DEFINITIONS
    )
    assert all(
        field_id not in unsupported.format
        for field_id, _field_type, _description in MODEL_SCORE_FORMAT_FIELD_DEFINITIONS
    )


@pytest.mark.parametrize("input_format", ["vcf", "bcf"])
@pytest.mark.parametrize(
    ("alt", "case_gt", "control_gt"),
    [
        ("T", "0|1/1", ".|1/0"),
        ("T,G", "0|1/2", "0/1|2"),
        ("T", "0|1", "1/0"),
        ("T", "./.", ".|."),
        ("T", "1", "."),
        ("T,G", ".|1/2|.", "0/.|2/1"),
    ],
)
def test_forced_annotation_preserves_genotype_phase_encoding(
    tmp_path, input_format, alt, case_gt, control_gt,
) -> None:
    import pysam

    input_path = tmp_path / "input.vcf"
    output_path = tmp_path / "output.vcf"
    input_path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
        '##INFO=<ID=SKUA_STATUS,Number=1,Type=String,Description="Old status">\n'
        '##FORMAT=<ID=SKUA_OLD,Number=1,Type=Integer,Description="Old count">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=CALLER_NOTE,Number=1,Type=String,Description="Caller note">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE\tCONTROL\n"
        f"chr1\t106\ths1\tA\t{alt}\t30\tPASS\tSKUA_STATUS=ANNOTATED\t"
        f"GT:SKUA_OLD:DP:CALLER_NOTE\t{case_gt}:99:15:case-note\t"
        f"{control_gt}:99:20:control-note\n"
    )
    if input_format == "bcf":
        bcf_path = tmp_path / "input.bcf"
        with pysam.VariantFile(str(input_path)) as source:
            with pysam.VariantFile(str(bcf_path), "wb", header=source.header) as output:
                for record in source:
                    output.write(record)
        input_path = bcf_path
    alignment = FakeAlignmentFile(
        [], header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )
    annotate_vcf(alignment, input_path, output_path=output_path, force=True)

    with pysam.VariantFile(str(output_path)) as output:
        record = next(iter(output))
        assert "SKUA_OLD" not in output.header.formats
        assert "SKUA_OLD" not in record.format
        expected_status = "UNSUPPORTED_MULTIALLELIC" if "," in alt else "ANNOTATED"
        assert record.info["SKUA_STATUS"] == expected_status
        columns = str(record).rstrip("\n").split("\t")
        fields = columns[8].split(":")
        for sample_name, values, gt, depth, note in (
            ("CASE", columns[9], case_gt, 15, "case-note"),
            ("CONTROL", columns[10], control_gt, 20, "control-note"),
        ):
            assert dict(zip(fields, values.split(":")))["GT"] == gt
            assert record.samples[sample_name]["DP"] == depth
            assert record.samples[sample_name]["CALLER_NOTE"] == note


def test_annotate_vcf_force_replaces_incompatible_owned_definition(tmp_path) -> None:
    import pysam

    input_path = tmp_path / "preannotated.vcf"
    input_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=SKUA_STATUS,Number=2,Type=Integer,Description="Wrong">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "reannotated.vcf"
    alignment = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )

    annotate_vcf(alignment, input_path, output_path=output_path, force=True)

    with pysam.VariantFile(str(output_path)) as reannotated:
        status = reannotated.header.info["SKUA_STATUS"]
        assert status.number == 1
        assert status.type == "String"
        assert next(iter(reannotated)).info["SKUA_STATUS"] == "ANNOTATED"


@pytest.mark.parametrize("use_symlink", [False, True], ids=["same-path", "resolved-same-path"])
@pytest.mark.parametrize("with_normals", [False, True], ids=["case", "case-with-normals"])
def test_annotate_vcf_rejects_output_path_that_would_overwrite_input(
    tmp_path, use_symlink, with_normals
) -> None:
    alignment_file = FakeAlignmentFile([])
    vcf_path = tmp_path / "input.vcf"
    original_contents = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1>",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
            "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
        ]
    ) + "\n"
    vcf_path.write_text(original_contents, encoding="utf-8")
    output_path = vcf_path
    if use_symlink:
        output_path = tmp_path / "output-link.vcf"
        output_path.symlink_to(vcf_path)

    with pytest.raises(ValueError, match="input VCF"):
        if with_normals:
            annotate_vcf_with_normals(
                alignment_file,
                vcf_path,
                normal_alignments=[],
                output_path=output_path,
            )
        else:
            annotate_vcf(alignment_file, vcf_path, output_path=output_path)

    assert vcf_path.read_text(encoding="utf-8") == original_contents


def test_annotate_vcf_counts_lowercase_alleles_as_alt_evidence(tmp_path) -> None:
    import pysam

    alignment_file = FakeAlignmentFile(
        [
            FakeRead(
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAATAAAA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
                tags={"RG": "case-rg"},
            ),
        ],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\ta\tt\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"

    annotate_vcf(alignment_file, vcf_path, output_path=output_path)

    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        sample = record.samples["CASE"]
        assert sample["SKUA_ALT_FWD"] == 1
        assert sample["SKUA_NON_ALT_FWD"] == 0


def test_annotate_vcf_with_normals_adds_sample_for_site_only_vcf(tmp_path) -> None:
    import pysam

    case_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "case-rg"},
        ),
    ]
    case_alignment = FakeAlignmentFile(
        case_reads,
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
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

    assert annotate_vcf_with_normals(
        case_alignment,
        vcf_path,
        normal_alignments=[],
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    ) is None
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
            output_path=tmp_path / "unused.vcf",
            min_baseq=20,
            min_mapq=20,
        )

    multi_sm_alignment = FakeAlignmentFile(
        case_reads,
        header=FakeAlignmentHeader([{"SM": "CASE"}, {"SM": "TUMOR"}]),
    )
    with pytest.raises(ValueError, match="multiple samples; specify --sample"):
        annotate_vcf_with_normals(
            multi_sm_alignment,
            vcf_path,
            normal_alignments=[],
            output_path=tmp_path / "unused.vcf",
            min_baseq=20,
            min_mapq=20,
        )


def test_annotate_vcf_selects_the_only_matching_case_sample_and_filters_reads(tmp_path) -> None:
    import pysam

    case_alignment = FakeAlignmentFile(
        [
            FakeRead(
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAATAAAA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
                tags={"RG": "case-rg"},
            ),
            FakeRead(
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAAAAAAA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
                tags={"RG": "other-rg"},
            ),
        ],
        header=FakeAlignmentHeader(
            [
                {"ID": "case-rg", "SM": "CASE"},
                {"ID": "other-rg", "SM": "UNRELATED"},
            ]
        ),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE\tCONTROL",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"

    annotate_vcf(case_alignment, vcf_path, output_path=output_path)

    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        assert record.samples["CASE"]["SKUA_ALT_FWD"] == 1
        assert record.samples["CASE"]["SKUA_NON_ALT_FWD"] == 0
        assert record.samples["CONTROL"]["SKUA_ALT_FWD"] is None


def test_annotate_vcf_filters_untagged_and_sm_less_read_groups_for_single_sample_alignment(
    tmp_path,
) -> None:
    import pysam

    case_alignment = FakeAlignmentFile(
        [
            FakeRead(
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAATAAAA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
                tags={"RG": "case-rg"},
            ),
            FakeRead(
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAAAAAAA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
                tags={"RG": "unassigned-rg"},
            ),
            FakeRead(
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAAAAAAA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
            ),
        ],
        header=FakeAlignmentHeader(
            [
                {"ID": "case-rg", "SM": "CASE"},
                {"ID": "unassigned-rg"},
            ]
        ),
    )
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

    annotate_vcf(case_alignment, vcf_path, output_path=output_path)

    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        record = next(iter(annotated_vcf))
        assert record.samples["CASE"]["SKUA_ALT_FWD"] == 1
        assert record.samples["CASE"]["SKUA_NON_ALT_FWD"] == 0
        assert record.samples["CASE"]["SKUA_USABLE"] == 1


def test_annotate_vcf_requires_explicit_sample_when_multiple_samples_match(tmp_path) -> None:
    case_alignment = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader(
            [
                {"ID": "case-rg", "SM": "CASE"},
                {"ID": "control-rg", "SM": "CONTROL"},
            ]
        ),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE\tCONTROL",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="--sample"):
        annotate_vcf(case_alignment, vcf_path, output_path=tmp_path / "unused.vcf")


def test_annotate_vcf_rejects_multi_sample_normal_alignment(tmp_path) -> None:
    case_alignment = FakeAlignmentFile([])
    normal_alignment = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader(
            [
                {"ID": "normal-a", "SM": "NORMAL_A"},
                {"ID": "normal-b", "SM": "NORMAL_B"},
            ]
        ),
    )
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

    with pytest.raises(ValueError, match="Normal alignment 1"):
        annotate_vcf_with_normals(
            case_alignment,
            vcf_path,
            normal_alignments=[normal_alignment],
            output_path=tmp_path / "unused.vcf",
        )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"min_baseq": -1}, "min_baseq"),
        ({"min_mapq": -1}, "min_mapq"),
        ({"truncate": 0}, "truncate"),
        ({"pseudocount": 0}, "pseudocount"),
        ({"prior_artifact_probability": 1}, "prior_artifact_probability"),
    ],
)
def test_annotate_vcf_with_normals_rejects_invalid_parameters(tmp_path, kwargs, message) -> None:
    alignment_file = FakeAlignmentFile([])
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

    with pytest.raises(ValueError, match=message):
        annotate_vcf_with_normals(
            alignment_file,
            vcf_path,
            normal_alignments=[],
            output_path=tmp_path / "unused.vcf",
            **kwargs,
        )


def test_annotate_vcf_rejects_reference_mismatch_before_writing_output(tmp_path) -> None:
    import pysam

    alignment_file = FakeAlignmentFile([], references=("chr1",))
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tG\tT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(">chr1\n" + "A" * 200 + "\n", encoding="utf-8")
    pysam.faidx(str(reference_path))
    output_path = tmp_path / "annotated.vcf"

    with pytest.raises(ValueError, match="REF allele"):
        annotate_vcf(
            alignment_file,
            vcf_path,
            output_path=output_path,
            reference_path=reference_path,
        )

    assert not output_path.exists()


def test_annotate_vcf_rejects_missing_case_contig_before_writing_output(tmp_path) -> None:
    alignment_file = FakeAlignmentFile([], references=("chr2",))
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

    with pytest.raises(ValueError, match="Case alignment does not contain contig 'chr1'"):
        annotate_vcf(alignment_file, vcf_path, output_path=output_path)

    assert not output_path.exists()


def test_annotate_vcf_rejects_alignment_without_an_index(tmp_path) -> None:
    alignment_file = FakeAlignmentFile([], indexed=False)
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

    with pytest.raises(ValueError, match="Case alignment must be indexed"):
        annotate_vcf(alignment_file, vcf_path, output_path=tmp_path / "unused.vcf")


def test_annotate_vcf_with_normals_reports_record_statuses_and_summary(tmp_path) -> None:
    import pysam

    alignment_file = FakeAlignmentFile(
        [
            FakeRead(
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAATAAAA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
                tags={"RG": "case-rg"},
            )
        ],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
                "chr1\t107\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/1",
                "chr1\t108\t.\tA\t<DEL>\t.\tPASS\t.\tGT\t0/1",
                "chr1\t109\t.\tA\tA]chr2:42]\t.\tPASS\t.\tGT\t0/1",
                "chr1\t110\t.\tAT\tGCA\t.\tPASS\t.\tGT\t0/1",
                "chr1\t111\t.\tA\t.\t.\tPASS\t.\tGT\t0/0",
                "chr1\t112\t.\tA\tA\t.\tPASS\t.\tGT\t0/1",
                "chr1\t113\t.\tA\tCT\t.\tPASS\t.\tGT\t0/1",
                "chr1\t114\t.\tAT\tC\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"

    annotate_vcf_with_normals(
        alignment_file,
        vcf_path,
        normal_alignments=[],
        output_path=output_path,
    )

    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        records = list(annotated_vcf)
        assert [record.info["SKUA_STATUS"] for record in records] == [
            "ANNOTATED",
            "UNSUPPORTED_MULTIALLELIC",
            "UNSUPPORTED_SYMBOLIC_ALLELE",
            "UNSUPPORTED_BREAKEND",
            "UNSUPPORTED_COMPLEX_ALLELE",
            "UNSUPPORTED_RECORD",
            "UNSUPPORTED_COMPLEX_ALLELE",
            "UNSUPPORTED_COMPLEX_ALLELE",
            "UNSUPPORTED_COMPLEX_ALLELE",
        ]
        assert records[0].samples["CASE"]["SKUA_ALT_FWD"] == 1
        assert "SKUA_ALT_FWD" not in records[1].format
        assert all("SKUA_ALT_FWD" not in record.format for record in records[6:])


def test_annotate_vcf_with_normals_strict_mode_rejects_unsupported_records_before_output(
    tmp_path,
) -> None:
    alignment_file = FakeAlignmentFile([])
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tCT\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"

    with pytest.raises(ValueError, match="UNSUPPORTED_COMPLEX_ALLELE"):
        annotate_vcf_with_normals(
            alignment_file,
            vcf_path,
            normal_alignments=[],
            output_path=output_path,
            strict=True,
        )

    assert not output_path.exists()


def test_annotate_vcf_with_normals_writes_info_and_format(tmp_path) -> None:
    import pysam

    case_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "case-rg"},
        ),
    ]
    case_alignment = FakeAlignmentFile(
        case_reads,
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )

    normal_reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "normal-rg"},
        ),
        FakeRead(
            mapping_quality=5,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
            tags={"RG": "normal-rg"},
        ),
    ]
    normal_alignment = FakeAlignmentFile(
        normal_reads,
        header=FakeAlignmentHeader([{"ID": "normal-rg", "SM": "NORMAL"}]),
    )

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

    annotate_vcf_with_normals(
        case_alignment,
        vcf_path,
        normal_alignments=[normal_alignment],
        output_path=output_path,
        min_baseq=20,
        min_mapq=20,
    )

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


def test_annotate_vcf_with_normals_uses_and_writes_effective_artifact_priors(
    tmp_path,
) -> None:
    import pysam

    alignment_file = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior probability that the ALT allele is an artifact before Skua evidence">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\texplicit\tA\tT\t.\tPASS\tSKUA_ARTIFACT_PRIOR=0.2\tGT\t0/1",
                "chr1\t107\tabsent\tA\tC\t.\tPASS\t.\tGT\t0/1",
                "chr1\t108\tmissing\tA\tG\t.\tPASS\tSKUA_ARTIFACT_PRIOR=.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    compressed_vcf_path = tmp_path / "input.vcf.gz"
    pysam.tabix_compress(str(vcf_path), str(compressed_vcf_path), force=True)
    output_path = tmp_path / "annotated.vcf"

    annotate_vcf_with_normals(
        alignment_file,
        compressed_vcf_path,
        normal_alignments=[],
        output_path=output_path,
        prior_artifact_probability=0.25,
    )

    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        prior_header = annotated_vcf.header.info["SKUA_ARTIFACT_PRIOR"]
        assert prior_header.number == "A"
        assert prior_header.type == "Float"
        records = list(annotated_vcf)

    assert [record.info["SKUA_ARTIFACT_PRIOR"][0] for record in records] == pytest.approx(
        [0.2, 0.25, 0.25]
    )
    assert [
        record.samples["CASE"]["SKUA_LOG_BAYES_FACTOR"] for record in records
    ] == pytest.approx([0.0, 0.0, 0.0])
    assert [
        record.samples["CASE"]["SKUA_ARTIFACT_POSTERIOR"] for record in records
    ] == pytest.approx([0.2, 0.25, 0.25])


@pytest.mark.parametrize(
    ("raw_prior", "message"),
    [
        ("abc", "malformed"),
        ("", "malformed"),
        ("0.2,0.3", "exactly one value"),
        ("nan", "finite and between 0 and 1"),
        ("inf", "finite and between 0 and 1"),
        ("-inf", "finite and between 0 and 1"),
        ("0", "finite and between 0 and 1"),
        ("1", "finite and between 0 and 1"),
        ("-0.1", "finite and between 0 and 1"),
        ("1.1", "finite and between 0 and 1"),
    ],
)
def test_annotate_vcf_with_normals_rejects_invalid_artifact_prior_before_output(
    tmp_path,
    raw_prior,
    message,
) -> None:
    alignment_file = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior probability that the ALT allele is an artifact before Skua evidence">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                f"chr1\t106\t.\tA\tT\t.\tPASS\tSKUA_ARTIFACT_PRIOR={raw_prior}\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"

    with pytest.raises(
        ValueError,
        match=rf"SKUA_ARTIFACT_PRIOR.*chr1:106.*{message}",
    ):
        annotate_vcf_with_normals(
            alignment_file,
            vcf_path,
            normal_alignments=[],
            output_path=output_path,
        )

    assert not output_path.exists()


@pytest.mark.parametrize(
    "definition",
    [
        '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=1,Type=Float,Description="Wrong number">',
        '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=String,Description="Wrong type">',
    ],
)
def test_annotate_vcf_with_normals_rejects_incompatible_artifact_prior_header(
    tmp_path,
    definition,
) -> None:
    alignment_file = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    vcf_path = tmp_path / "input.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                definition,
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\tSKUA_ARTIFACT_PRIOR=0.2\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "annotated.vcf"

    with pytest.raises(ValueError, match="incompatible SKUA_ARTIFACT_PRIOR"):
        annotate_vcf_with_normals(
            alignment_file,
            vcf_path,
            normal_alignments=[],
            output_path=output_path,
        )

    assert not output_path.exists()


def test_annotate_vcf_with_normals_batches_dense_records_per_alignment(tmp_path) -> None:
    import pysam

    case_alignment = FakeAlignmentFile(
        [
            FakeRead(
                query_name="case-alt",
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAATAACA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
                tags={"RG": "case-rg"},
            )
        ],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
    )
    normal_alignment = FakeAlignmentFile(
        [
            FakeRead(
                query_name="normal-ref",
                mapping_quality=60,
                is_reverse=False,
                query_sequence="AAAAAAAAAA",
                query_qualities=[35] * 10,
                aligned_pairs=build_linear_pairs(10, 100),
                tags={"RG": "normal-rg"},
            )
        ],
        header=FakeAlignmentHeader([{"ID": "normal-rg", "SM": "NORMAL"}]),
    )
    vcf_path = tmp_path / "dense.vcf"
    vcf_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
                "chr1\t109\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "dense.annotated.vcf"

    annotate_vcf_with_normals(
        case_alignment,
        vcf_path,
        normal_alignments=[normal_alignment],
        output_path=output_path,
    )

    assert case_alignment.fetch_calls == [("chr1", 105, 109)]
    assert normal_alignment.fetch_calls == [("chr1", 105, 109)]
    with pysam.VariantFile(str(output_path)) as annotated_vcf:
        records = list(annotated_vcf)
        assert [record.samples["CASE"]["SKUA_ALT_FWD"] for record in records] == [1, 1]
        assert [record.info["SKUA_PON_NON_ALT_FWD"] for record in records] == [1, 1]


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
    result = json.loads(payload)["records"]
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
        "artifact_prior",
    ]
    assert list(result[0]["stats"].keys()) == [
        "artifact_posterior",
        "log_bayes_factor_artifact_vs_variant",
        "dispersion_factor",
        "pon_sample_count",
        "assessment_status",
        "assessment_reasons",
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
