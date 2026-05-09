import gzip

import pytest

from skua.variants import Variant, parse_vcf_variant_line, read_vcf_variant_file


def test_variant_from_vcf_fields_converts_pos1_to_pos0() -> None:
    variant = Variant.from_vcf_fields(contig="chr7", pos1=106, ref="A", alt="T")

    assert variant.contig == "chr7"
    assert variant.ref_pos0 == 105
    assert variant.ref == "A"
    assert variant.alt == "T"


def test_variant_from_vcf_fields_parses_simple_deletion() -> None:
    variant = Variant.from_vcf_fields(contig="chr1", pos1=10, ref="AT", alt="A")

    assert variant.contig == "chr1"
    assert variant.ref_pos0 == 9
    assert variant.ref == "AT"
    assert variant.alt == "A"
    assert variant.kind.value == "deletion"


def test_variant_from_vcf_fields_parses_simple_insertion() -> None:
    variant = Variant.from_vcf_fields(contig="chr1", pos1=10, ref="A", alt="AT")

    assert variant.contig == "chr1"
    assert variant.ref_pos0 == 9
    assert variant.ref == "A"
    assert variant.alt == "AT"
    assert variant.kind.value == "insertion"


def test_parse_vcf_variant_line_parses_data_line() -> None:
    line = "chr1\t106\t.\tA\tT\t.\tPASS\t." 

    variant = parse_vcf_variant_line(line)

    assert variant is not None
    assert variant.contig == "chr1"
    assert variant.ref_pos0 == 105
    assert variant.ref == "A"
    assert variant.alt == "T"


def test_parse_vcf_variant_line_skips_header_lines() -> None:
    assert parse_vcf_variant_line("##fileformat=VCFv4.2") is None
    assert parse_vcf_variant_line("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") is None


def test_parse_vcf_variant_line_parses_simple_indels() -> None:
    insertion = parse_vcf_variant_line("chr1\t106\t.\tA\tAT\t.\tPASS\t.")
    deletion = parse_vcf_variant_line("chr1\t106\t.\tAT\tA\t.\tPASS\t.")

    assert insertion == Variant(contig="chr1", ref_pos0=105, ref="A", alt="AT")
    assert deletion == Variant(contig="chr1", ref_pos0=105, ref="AT", alt="A")


def test_parse_vcf_variant_line_skips_multiallelic_records() -> None:
    assert parse_vcf_variant_line("chr1\t106\t.\tA\tT,C\t.\tPASS\t.") is None


def test_read_vcf_variant_file_yields_simple_records_only(tmp_path) -> None:
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

    variants = list(read_vcf_variant_file(vcf_path))

    assert variants == [
        Variant(contig="chr1", ref_pos0=105, ref="A", alt="T"),
        Variant(contig="chr1", ref_pos0=199, ref="A", alt="AT"),
        Variant(contig="chr1", ref_pos0=299, ref="AT", alt="A"),
    ]


def test_read_vcf_variant_file_yields_simple_records_from_gzipped_input(tmp_path) -> None:
    vcf_path = tmp_path / "input.vcf.gz"
    with gzip.open(vcf_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "##fileformat=VCFv4.2",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                    "chr1\t106\t.\tA\tT\t.\tPASS\t.",
                    "chr1\t200\t.\tA\tAT\t.\tPASS\t.",
                    "chr1\t300\t.\tAT\tA\t.\tPASS\t.",
                ]
            )
            + "\n"
        )

    variants = list(read_vcf_variant_file(vcf_path))

    assert variants == [
        Variant(contig="chr1", ref_pos0=105, ref="A", alt="T"),
        Variant(contig="chr1", ref_pos0=199, ref="A", alt="AT"),
        Variant(contig="chr1", ref_pos0=299, ref="AT", alt="A"),
    ]
