import gzip

import pytest

from skua.variants import Variant, parse_vcf_snv_line, read_vcf_snv_file


def test_variant_from_vcf_fields_converts_pos1_to_pos0() -> None:
    variant = Variant.from_vcf_fields(contig="chr7", pos1=106, ref="A", alt="T")

    assert variant.contig == "chr7"
    assert variant.ref_pos0 == 105
    assert variant.ref == "A"
    assert variant.alt == "T"


def test_variant_from_vcf_fields_rejects_non_snv_ref() -> None:
    with pytest.raises(ValueError, match="SNV"):
        Variant.from_vcf_fields(contig="chr1", pos1=10, ref="AT", alt="A")


def test_variant_from_vcf_fields_rejects_non_snv_alt() -> None:
    with pytest.raises(ValueError, match="SNV"):
        Variant.from_vcf_fields(contig="chr1", pos1=10, ref="A", alt="AT")


def test_parse_vcf_snv_line_parses_data_line() -> None:
    line = "chr1\t106\t.\tA\tT\t.\tPASS\t." 

    variant = parse_vcf_snv_line(line)

    assert variant is not None
    assert variant.contig == "chr1"
    assert variant.ref_pos0 == 105
    assert variant.ref == "A"
    assert variant.alt == "T"


def test_parse_vcf_snv_line_skips_header_lines() -> None:
    assert parse_vcf_snv_line("##fileformat=VCFv4.2") is None
    assert parse_vcf_snv_line("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") is None


def test_parse_vcf_snv_line_skips_non_snv_records() -> None:
    assert parse_vcf_snv_line("chr1\t106\t.\tA\tAT\t.\tPASS\t.") is None
    assert parse_vcf_snv_line("chr1\t106\t.\tAT\tA\t.\tPASS\t.") is None


def test_parse_vcf_snv_line_skips_multiallelic_records() -> None:
    assert parse_vcf_snv_line("chr1\t106\t.\tA\tT,C\t.\tPASS\t.") is None


def test_read_vcf_snv_file_yields_snv_records_only(tmp_path) -> None:
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

    variants = list(read_vcf_snv_file(vcf_path))

    assert variants == [Variant(contig="chr1", ref_pos0=105, ref="A", alt="T")]


def test_read_vcf_snv_file_yields_snv_records_from_gzipped_input(tmp_path) -> None:
    vcf_path = tmp_path / "input.vcf.gz"
    with gzip.open(vcf_path, "wt", encoding="utf-8") as handle:
        handle.write(
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

    variants = list(read_vcf_snv_file(vcf_path))

    assert variants == [Variant(contig="chr1", ref_pos0=105, ref="A", alt="T")]
