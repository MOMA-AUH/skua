import gzip

import pytest

from skua.variants import Variant, parse_vcf_variant_line, read_vcf_variant_file


def test_variant_from_vcf_fields_converts_pos1_to_pos0() -> None:
    variant = Variant.from_vcf_fields(contig="chr7", pos1=106, ref="A", alt="T")

    assert variant.contig == "chr7"
    assert variant.ref_pos0 == 105
    assert variant.ref == "A"
    assert variant.alt == "T"


def test_variant_from_vcf_fields_normalizes_alleles_to_uppercase() -> None:
    variant = Variant.from_vcf_fields(contig="chr7", pos1=106, ref="a", alt="t")

    assert variant.ref == "A"
    assert variant.alt == "T"


def test_direct_variant_construction_normalizes_alleles_to_uppercase() -> None:
    variant = Variant(contig="chr7", ref_pos0=105, ref="a", alt="t")

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


@pytest.mark.parametrize(
    ("ref", "alt"),
    [
        ("A", "A"),
        ("A", "CT"),
        ("AT", "C"),
    ],
)
def test_variant_rejects_identity_and_non_anchored_alleles(ref: str, alt: str) -> None:
    with pytest.raises(ValueError, match="simple|different"):
        Variant.from_vcf_fields(contig="chr1", pos1=10, ref=ref, alt=alt)

    with pytest.raises(ValueError, match="simple|different"):
        Variant(contig="chr1", ref_pos0=9, ref=ref, alt=alt)


@pytest.mark.parametrize(
    ("contig", "ref_pos0", "ref", "alt", "message"),
    [
        ("", 9, "A", "T", "contig"),
        ("chr1", -1, "A", "T", "ref_pos0"),
        ("chr1", 9, "", "T", "non-empty"),
        ("chr1", 9, "N", "T", "only A"),
        ("chr1", 9, "AT", "GCA", "simple"),
    ],
)
def test_direct_variant_construction_rejects_invalid_values(
    contig: str,
    ref_pos0: int,
    ref: str,
    alt: str,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        Variant(contig=contig, ref_pos0=ref_pos0, ref=ref, alt=alt)


def test_variant_from_vcf_fields_rejects_non_positive_position() -> None:
    with pytest.raises(ValueError, match="VCF POS"):
        Variant.from_vcf_fields(contig="chr1", pos1=0, ref="A", alt="T")


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


@pytest.mark.parametrize(
    ("ref", "alt"),
    [("A", "A"), ("A", "CT"), ("AT", "C")],
)
def test_parse_vcf_variant_line_skips_identity_and_non_anchored_alleles(
    ref: str,
    alt: str,
) -> None:
    assert parse_vcf_variant_line(f"chr1\t106\t.\t{ref}\t{alt}\t.\tPASS\t.") is None


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
