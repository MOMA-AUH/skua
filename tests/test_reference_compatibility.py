from pathlib import Path

import pysam
import pytest
from pysam import bcftools

from skua import annotate_vcf_with_normals, annotate_vcf_with_pon, build_pon, read_pon_metadata
from skua import annotate_variant_with_normals, annotate_variants_with_normals, annotate_variants_from_pon, Variant
from skua.pon import inspect_pon, validate_pon


def targets(path: Path) -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t21\t.\tA\tT\t.\tPASS\t.\n", encoding="utf-8",
    )
    return path


def alignment(path: Path, sample: str, length: int, md5: str | None = None, extra=None):
    sequence = {"SN": "chr1", "LN": length}
    if md5 is not None:
        sequence["M5"] = md5
    with pysam.AlignmentFile(str(path), "wb", header={
        "HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [sequence] + ([] if extra is None else [extra]),
        "RG": [{"ID": sample, "SM": sample}],
    }):
        pass
    pysam.index(str(path))
    return pysam.AlignmentFile(str(path), "rb")


def test_live_annotation_rejects_conflicting_reference_lengths_before_output(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    output = tmp_path / "result.vcf"
    with (
        alignment(tmp_path / "case.bam", "CASE", 200) as case,
        alignment(tmp_path / "normal.bam", "NORMAL", 300) as normal,
    ):
        with pytest.raises(ValueError, match="reference length.*chr1"):
            annotate_vcf_with_normals(case, vcf, normal_alignments=[normal], output_path=output)
    assert not output.exists()


def test_live_annotation_rejects_same_length_different_reference_checksums(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    output = tmp_path / "result.vcf"
    with (
        alignment(tmp_path / "case.bam", "CASE", 200, "a" * 32) as case,
        alignment(tmp_path / "normal.bam", "NORMAL", 200, "b" * 32) as normal,
    ):
        with pytest.raises(ValueError, match="reference checksum.*chr1"):
            annotate_vcf_with_normals(case, vcf, normal_alignments=[normal], output_path=output)
    assert not output.exists()


@pytest.mark.parametrize("mismatch", ["length", "checksum"])
def test_reference_fasta_is_compared_beyond_target_ref_bases(tmp_path, mismatch):
    vcf = targets(tmp_path / "targets.vcf")
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "A" * (300 if mismatch == "length" else 200) + "\n")
    pysam.faidx(str(fasta))
    with alignment(tmp_path / "case.bam", "CASE", 200, "a" * 32) as case:
        with pytest.raises(ValueError, match=f"reference {mismatch}.*chr1"):
            annotate_vcf_with_normals(
                case, vcf, output_path=tmp_path / "result.vcf", reference_path=fasta,
            )
    assert not (tmp_path / "result.vcf").exists()


@pytest.mark.parametrize("checksum,expected", [(None, "INSUFFICIENT_METADATA"), ("a" * 32, "VERIFIED")])
def test_output_distinguishes_matching_identity_from_missing_checksums(tmp_path, checksum, expected):
    vcf = targets(tmp_path / "targets.vcf")
    output = tmp_path / "result.vcf"
    with (
        alignment(tmp_path / "case.bam", "CASE", 200, checksum) as case,
        alignment(tmp_path / "normal.bam", "NORMAL", 200, checksum) as normal,
    ):
        annotate_vcf_with_normals(case, vcf, normal_alignments=[normal], output_path=output)
    with pysam.VariantFile(str(output)) as result:
        [status] = [r for r in result.header.records if r.key == "SKUA_REFERENCE_STATUS"]
        assert status.value == expected


def test_pon_round_trips_reference_identity_and_rejects_incompatible_case(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL", 200, "a" * 32) as normal:
        build_pon(vcf, normal_alignments=[normal], output_path=panel)
    metadata = read_pon_metadata(panel)
    assert metadata.schema_version == 2
    assert metadata.reference_identity.status == "VERIFIED"
    [contig] = metadata.reference_identity.contigs
    assert (contig.name, contig.length, contig.md5) == ("chr1", 200, "a" * 32)
    with alignment(tmp_path / "case.bam", "CASE", 200, "b" * 32) as case:
        for source in (None, vcf):
            with pytest.raises(ValueError, match="reference checksum.*chr1"):
                annotate_vcf_with_pon(
                    case, panel, vcf_path=source, output_path=tmp_path / "result.vcf",
                )
    assert not (tmp_path / "result.vcf").exists()


def test_validate_pon_checks_full_reference_identity_and_reports_missing_metadata(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL", 200, "a" * 32) as normal:
        build_pon(vcf, normal_alignments=[normal], output_path=panel)
    assert inspect_pon(panel).as_dict()["reference_status"] == "VERIFIED"
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "A" * 200 + "\n")
    pysam.faidx(str(fasta))
    validation = validate_pon(panel, reference_path=fasta)
    assert not validation.valid
    assert any("reference checksum" in error for error in validation.errors)


def test_evidence_only_apis_reject_conflicting_alignment_and_panel_identities(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    variant = Variant.from_vcf_fields(contig="chr1", pos1=21, ref="A", alt="T")
    with (
        alignment(tmp_path / "case.bam", "CASE", 200, "a" * 32) as case,
        alignment(tmp_path / "normal.bam", "NORMAL", 200, "b" * 32) as normal,
    ):
        build_pon(vcf, normal_alignments=[normal], output_path=panel)
        with pytest.raises(ValueError, match="reference checksum"):
            annotate_variant_with_normals(case, variant, normal_alignments=[normal])
        with pytest.raises(ValueError, match="reference checksum"):
            list(annotate_variants_with_normals(case, [variant], normal_alignments=[normal]))
        with pytest.raises(ValueError, match="reference checksum"):
            list(annotate_variants_from_pon(case, panel))


@pytest.mark.parametrize("normal_checksum,expected", [
    (None, "INSUFFICIENT_METADATA"), ("16bf06b3717d1f238252870e699c2a2e", "VERIFIED"),
])
def test_matching_fasta_and_cached_identity_preserve_verification_strength(tmp_path, normal_checksum, expected):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "a" * 200 + "\n")
    pysam.faidx(str(fasta))
    with alignment(tmp_path / "normal.bam", "NORMAL", 200, normal_checksum) as normal:
        build_pon(vcf, normal_alignments=[normal], output_path=panel, reference_path=fasta)
    assert read_pon_metadata(panel).reference_identity.status == expected
    assert validate_pon(panel, reference_path=fasta).valid
    with alignment(tmp_path / "case.bam", "CASE", 200, "16bf06b3717d1f238252870e699c2a2e") as case:
        output = tmp_path / "result.vcf"
        annotate_vcf_with_pon(case, panel, reference_path=fasta, output_path=output)
        with pysam.VariantFile(str(output)) as result:
            [status] = [r.value for r in result.header.records if r.key == "SKUA_REFERENCE_STATUS"]
            assert status == expected


def test_cached_subset_only_requires_reference_identity_on_used_contigs(tmp_path):
    subset = targets(tmp_path / "subset.vcf")
    all_targets = tmp_path / "all.vcf"
    all_targets.write_text(subset.read_text().replace(
        "##contig=<ID=chr1>", "##contig=<ID=chr1>\n##contig=<ID=chr2>",
    ) + "chr2\t21\t.\tA\tT\t.\tPASS\t.\n")
    panel = tmp_path / "panel.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL", 200, "a" * 32,
                   extra={"SN": "chr2", "LN": 300, "M5": "b" * 32}) as normal:
        build_pon(all_targets, normal_alignments=[normal], output_path=panel)
    with alignment(tmp_path / "case.bam", "CASE", 200, "a" * 32) as case:
        output = tmp_path / "result.vcf"
        annotate_vcf_with_pon(case, panel, vcf_path=subset, output_path=output)
        with pysam.VariantFile(str(output)) as result:
            assert [r.contig for r in result] == ["chr1"]


def test_conflicting_unused_contigs_do_not_reject_live_annotation(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    with (
        alignment(tmp_path / "case.bam", "CASE", 200, "a" * 32,
                  extra={"SN": "unused", "LN": 300}) as case,
        alignment(tmp_path / "normal.bam", "NORMAL", 200, "a" * 32,
                  extra={"SN": "unused", "LN": 400}) as normal,
    ):
        annotate_vcf_with_normals(case, vcf, normal_alignments=[normal], output_path=tmp_path / "out.vcf")


@pytest.mark.parametrize("attribute,value", [("length", "300"), ("md5", "b" * 32)])
def test_target_vcf_reference_metadata_is_checked(tmp_path, attribute, value):
    vcf = targets(tmp_path / "targets.vcf")
    vcf.write_text(vcf.read_text().replace("ID=chr1>", f"ID=chr1,{attribute}={value}>"))
    with alignment(tmp_path / "case.bam", "CASE", 200, "a" * 32) as case:
        with pytest.raises(ValueError, match="Conflicting reference"):
            annotate_vcf_with_normals(case, vcf, output_path=tmp_path / "out.vcf")


def test_pon_validation_rejects_conflicting_target_dictionary(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL", 200, "a" * 32) as normal:
        build_pon(vcf, normal_alignments=[normal], output_path=panel)
    vcf.write_text(vcf.read_text().replace("ID=chr1>", "ID=chr1,length=300>"))
    validation = validate_pon(panel, target_vcf_path=vcf)
    assert not validation.valid
    assert any("reference length" in error for error in validation.errors)


@pytest.mark.parametrize("old,new,error", [
    ('Length="200"', 'Length="0"', "Invalid reference length"),
    ('MD5="' + "a" * 32 + '"', 'MD5="invalid"', "Invalid reference checksum"),
    (',MD5="' + "a" * 32 + '"', "", "Invalid reference verification"),
    ("SKUA_REFERENCE_STATUS=VERIFIED", "SKUA_REFERENCE_STATUS=INSUFFICIENT_METADATA", "inconsistent"),
    ('EvidencePolicyVersion="4"', 'EvidencePolicyVersion="999"', "Unsupported PON evidence policy"),
    ('ID=chr1,Verified="1"', 'ID=other,Verified="1"', "missing reference identity"),
])
def test_malformed_panel_reference_metadata_cannot_publish_output(tmp_path, old, new, error):
    vcf = targets(tmp_path / "targets.vcf")
    panel = tmp_path / "panel.bcf"
    with alignment(tmp_path / "normal.bam", "NORMAL", 200, "a" * 32) as normal:
        build_pon(vcf, normal_alignments=[normal], output_path=panel)
    with pysam.VariantFile(str(panel)) as source:
        text = str(source.header) + "".join(str(record) for record in source)
    assert old in text
    damaged_vcf = tmp_path / "damaged.vcf"
    damaged_vcf.write_text(text.replace(old, new))
    damaged_panel = tmp_path / "damaged.bcf"
    with pysam.VariantFile(str(damaged_vcf)) as source:
        with pysam.VariantFile(str(damaged_panel), "wb", header=source.header) as output:
            for record in source:
                output.write(record)
    bcftools.index(str(damaged_panel))
    assert any(error in message for message in validate_pon(damaged_panel).errors)
    destination = tmp_path / "result.vcf"
    destination.write_bytes(b"original output")
    with alignment(tmp_path / "case.bam", "CASE", 200, "a" * 32) as case:
        for target in (None, vcf):
            with pytest.raises(ValueError, match=error):
                annotate_vcf_with_pon(case, damaged_panel, vcf_path=target,
                                      output_path=destination, force=True)
    assert destination.read_bytes() == b"original output"


def test_forced_annotation_refreshes_reference_metadata_preserving_other_headers(tmp_path):
    vcf = targets(tmp_path / "targets.vcf")
    vcf.write_text(vcf.read_text().replace("##contig", "##external=keep\n##contig"))
    first = tmp_path / "first.vcf"
    second = tmp_path / "second.vcf"
    with alignment(tmp_path / "case.bam", "CASE", 200, "a" * 32) as case:
        annotate_vcf_with_normals(case, vcf, output_path=first)
    with alignment(tmp_path / "new-case.bam", "CASE", 200) as case:
        annotate_vcf_with_normals(case, first, output_path=second, force=True)
    with pysam.VariantFile(str(second)) as result:
        statuses = [r.value for r in result.header.records if r.key == "SKUA_REFERENCE_STATUS"]
        assert statuses == ["INSUFFICIENT_METADATA"]
        assert "##external=keep" in str(result.header)
        assert "MD5=" not in str(result.header)
