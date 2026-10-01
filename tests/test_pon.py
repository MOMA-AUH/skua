import pysam
import pytest
from pysam import bcftools

from skua import (
    annotate_vcf_with_normals,
    annotate_vcf_with_pon,
    build_pon,
    read_pon_evidence,
    read_pon_metadata,
)
from skua.pon import inspect_pon, validate_pon
from tests.helpers import (
    FakeAlignmentFile,
    FakeAlignmentHeader,
    FakeRead,
    build_linear_pairs,
)


def _write_targets(path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=HOTSPOT,Number=0,Type=Flag,Description="Known hotspot">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\ths1\tA\tT\t.\tPASS\tHOTSPOT",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _read(sequence: str, *, reverse: bool = False, read_group: str | None = None) -> FakeRead:
    tags = {} if read_group is None else {"RG": read_group}
    return FakeRead(
        mapping_quality=60,
        is_reverse=reverse,
        query_sequence=sequence,
        query_qualities=[35] * len(sequence),
        aligned_pairs=build_linear_pairs(len(sequence), 100),
        tags=tags,
    )


def _normal(sample_name: str, reads: list[FakeRead]) -> FakeAlignmentFile:
    return FakeAlignmentFile(
        reads,
        header=FakeAlignmentHeader([{"ID": f"{sample_name}-rg", "SM": sample_name}]),
        references=("chr1",),
    )


def _write_legacy_pon(path, policy_version: int) -> None:
    """Write a structurally valid PON produced under an older evidence policy."""
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000)
    header.add_meta(
        "SKUA_PON",
        items=[
            ("SchemaVersion", "1"),
            ("EvidencePolicyVersion", str(policy_version)),
            ("MinBaseQ", "20"),
            ("MinMapQ", "20"),
            ("SkuaVersion", "0.6.0" if policy_version == 1 else "0.7.1"),
        ],
    )
    for field_id in (
        "SKUA_PON_AF",
        "SKUA_PON_AR",
        "SKUA_PON_NF",
        "SKUA_PON_NR",
        "SKUA_PON_U",
        "SKUA_PON_X",
    ):
        header.add_line(
            f'##FORMAT=<ID={field_id},Number=1,Type=Integer,Description="Legacy evidence">'
        )
    header.add_sample("N1")

    with pysam.VariantFile(str(path), "wb", header=header) as pon_file:
        record = pon_file.new_record(
            contig="chr1",
            start=105,
            stop=106,
            alleles=("A", "T"),
        )
        sample = record.samples["N1"]
        sample["SKUA_PON_AF"] = 0
        sample["SKUA_PON_AR"] = 0
        sample["SKUA_PON_NF"] = 1
        sample["SKUA_PON_NR"] = 0
        sample["SKUA_PON_U"] = 1
        sample["SKUA_PON_X"] = 0
        pon_file.write(record)

    bcftools.index("--force", str(path))


def test_build_pon_round_trips_per_sample_evidence_and_metadata(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    output_path = tmp_path / "hotspots.pon.bcf"
    _write_targets(target_path)
    normals = [
        _normal("N1", [_read("AAAAATAAAA")]),
        _normal("N2", [_read("AAAAAAAAAA", reverse=True)]),
    ]

    build_pon(
        target_path,
        normal_alignments=normals,
        output_path=output_path,
        min_baseq=25,
        min_mapq=30,
    )

    metadata = read_pon_metadata(output_path)
    assert metadata.schema_version == 1
    assert metadata.evidence_policy_version == 3
    assert metadata.min_baseq == 25
    assert metadata.min_mapq == 30
    assert metadata.sample_names == ("N1", "N2")

    [(variant, evidences)] = list(read_pon_evidence(output_path))
    assert (variant.contig, variant.ref_pos0, variant.ref, variant.alt) == (
        "chr1",
        105,
        "A",
        "T",
    )
    assert evidences[0].alt_forward == 1
    assert evidences[0].usable == 1
    assert evidences[1].non_alt_reverse == 1
    assert evidences[1].usable == 1

    with pysam.VariantFile(str(output_path)) as artifact:
        record = next(iter(artifact))
        assert record.id == "hs1"
        assert record.info["HOTSPOT"]


def test_build_pon_indexes_bcf_so_opening_it_does_not_log_an_index_error(
    tmp_path,
    capfd,
) -> None:
    target_path = tmp_path / "hotspots.vcf"
    output_path = tmp_path / "hotspots.pon.bcf"
    _write_targets(target_path)

    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [_read("AAAAAAAAAA")])],
        output_path=output_path,
    )
    capfd.readouterr()

    read_pon_metadata(output_path)
    list(read_pon_evidence(output_path))

    assert (tmp_path / "hotspots.pon.bcf.csi").exists()
    assert "Could not retrieve index file" not in capfd.readouterr().err


def test_build_pon_does_not_clobber_an_existing_artifact_by_default(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    output_path = tmp_path / "hotspots.pon.bcf"
    _write_targets(target_path)
    output_path.write_bytes(b"existing bcf")
    index_path = tmp_path / "hotspots.pon.bcf.csi"
    index_path.write_bytes(b"existing index")

    with pytest.raises(FileExistsError, match="already exists"):
        build_pon(
            target_path,
            normal_alignments=[_normal("N1", [])],
            output_path=output_path,
        )

    assert output_path.read_bytes() == b"existing bcf"
    assert index_path.read_bytes() == b"existing index"


def test_build_pon_force_strips_target_annotations_but_preserves_prior(tmp_path) -> None:
    target_path = tmp_path / "preannotated.vcf"
    output_path = tmp_path / "hotspots.pon.bcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=CALLER_SCORE,Number=1,Type=Integer,Description="Caller score">',
                '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior">',
                '##INFO=<ID=SKUA_STATUS,Number=1,Type=String,Description="Old">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=SKUA_ALT_FWD,Number=1,Type=Integer,Description="Old">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASE",
                "chr1\t106\ths1\tA\tT\t.\tPASS\t"
                "CALLER_SCORE=7;SKUA_ARTIFACT_PRIOR=0.2;SKUA_STATUS=ANNOTATED\t"
                "GT:SKUA_ALT_FWD\t0/1:99",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [])],
        output_path=output_path,
        force=True,
    )

    with pysam.VariantFile(str(output_path)) as artifact:
        assert "SKUA_STATUS" not in artifact.header.info
        assert "SKUA_ALT_FWD" not in artifact.header.formats
        record = next(iter(artifact))
        assert record.info["CALLER_SCORE"] == 7
        assert record.info["SKUA_ARTIFACT_PRIOR"][0] == pytest.approx(0.2)


def test_build_pon_force_preserves_existing_pair_when_indexing_fails(
    monkeypatch,
    tmp_path,
) -> None:
    import skua.pon as pon

    target_path = tmp_path / "hotspots.vcf"
    output_path = tmp_path / "hotspots.pon.bcf"
    _write_targets(target_path)
    normal = _normal("N1", [])
    build_pon(
        target_path,
        normal_alignments=[normal],
        output_path=output_path,
    )
    index_path = tmp_path / "hotspots.pon.bcf.csi"
    old_bcf = output_path.read_bytes()
    old_csi = index_path.read_bytes()

    def fail_index(*args) -> None:
        raise pysam.SamtoolsError("forced index failure")

    monkeypatch.setattr(pon.bcftools, "index", fail_index)

    with pytest.raises(ValueError, match="coordinate-sorted"):
        build_pon(
            target_path,
            normal_alignments=[normal],
            output_path=output_path,
            min_baseq=20,
            min_mapq=20,
            force=True,
        )

    assert output_path.read_bytes() == old_bcf
    assert index_path.read_bytes() == old_csi
    assert list(tmp_path.glob(".*hotspots.pon.bcf*")) == []


def test_inspect_pon_reports_header_metadata_without_scanning_targets(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    output_path = tmp_path / "hotspots.pon.bcf"
    _write_targets(target_path)

    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [_read("AAAAAAAAAA")])],
        output_path=output_path,
        min_baseq=25,
        min_mapq=30,
    )

    inspection = inspect_pon(output_path)

    assert inspection.format == "BCF"
    assert inspection.index_present is True
    assert inspection.metadata_record_count == 1
    assert inspection.schema_version == "1"
    assert inspection.evidence_policy_version == "3"
    assert inspection.min_baseq == "25"
    assert inspection.min_mapq == "30"
    assert inspection.sample_names == ("N1",)


@pytest.mark.parametrize("policy_version", [1, 2])
def test_legacy_pon_can_be_inspected_but_not_read(tmp_path, policy_version) -> None:
    pon_path = tmp_path / "legacy.pon.bcf"
    _write_legacy_pon(pon_path, policy_version)

    inspection = inspect_pon(pon_path)

    assert inspection.evidence_policy_version == str(policy_version)
    with pytest.raises(
        ValueError,
        match=f"Unsupported PON evidence policy version {policy_version}; expected 3",
    ):
        read_pon_metadata(pon_path)
    with pytest.raises(
        ValueError,
        match=f"Unsupported PON evidence policy version {policy_version}; expected 3",
    ):
        list(read_pon_evidence(pon_path))


@pytest.mark.parametrize("policy_version", [1, 2])
def test_validate_pon_rejects_legacy_artifact(tmp_path, policy_version) -> None:
    pon_path = tmp_path / "legacy.pon.bcf"
    _write_legacy_pon(pon_path, policy_version)

    result = validate_pon(pon_path)

    assert result.valid is False
    assert result.errors == (
        f"Unsupported PON evidence policy version {policy_version}; expected 3",
    )


@pytest.mark.parametrize("policy_version", [1, 2])
def test_cached_annotation_rejects_legacy_pon_before_output(tmp_path, policy_version) -> None:
    pon_path = tmp_path / "legacy.pon.bcf"
    output_path = tmp_path / "annotated.vcf"
    _write_legacy_pon(pon_path, policy_version)
    case = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )

    with pytest.raises(
        ValueError,
        match=f"Unsupported PON evidence policy version {policy_version}; expected 3",
    ):
        annotate_vcf_with_pon(case, pon_path, output_path=output_path)

    assert not output_path.exists()


def test_validate_pon_checks_records_and_requires_a_csi_index(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    output_path = tmp_path / "hotspots.pon.bcf"
    _write_targets(target_path)

    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [_read("AAAAAAAAAA")])],
        output_path=output_path,
    )

    assert validate_pon(output_path).valid is True

    (tmp_path / "hotspots.pon.bcf.csi").unlink()
    result = validate_pon(output_path)

    assert result.valid is False
    assert result.inspection is not None
    assert result.errors == ("PON artifact is missing its .csi index",)


def test_validate_pon_optionally_checks_reference_and_target_vcf(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    output_path = tmp_path / "hotspots.pon.bcf"
    reference_path = tmp_path / "reference.fa"
    mismatched_targets = tmp_path / "different-targets.vcf"
    _write_targets(target_path)
    reference_path.write_text(">chr1\n" + "A" * 200 + "\n", encoding="utf-8")
    pysam.faidx(str(reference_path))
    mismatched_targets.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\ths1\tA\tC\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [_read("AAAAAAAAAA")])],
        output_path=output_path,
    )

    assert validate_pon(
        output_path,
        reference_path=reference_path,
        target_vcf_path=target_path,
    ).valid is True

    result = validate_pon(output_path, target_vcf_path=mismatched_targets)

    assert result.valid is False
    assert result.errors == ("PON artifact targets do not match the supplied target VCF",)


def test_validate_pon_accepts_a_fresh_coordinate_sorted_multicontig_artifact(tmp_path) -> None:
    target_path = tmp_path / "multicontig-targets.vcf"
    output_path = tmp_path / "multicontig.pon.bcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##contig=<ID=chr2>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\t.",
                "chr2\t106\t.\tA\tC\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    normal = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "N1-rg", "SM": "N1"}]),
        references=("chr1", "chr2"),
    )

    build_pon(target_path, normal_alignments=[normal], output_path=output_path)

    assert validate_pon(output_path).valid is True


def test_build_pon_rejects_an_unsorted_target_vcf_before_writing_an_artifact(tmp_path) -> None:
    target_path = tmp_path / "unsorted-targets.vcf"
    output_path = tmp_path / "unsorted.pon.bcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t109\t.\tA\tT\t.\tPASS\t.",
                "chr1\t106\t.\tA\tC\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    normal = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "N1-rg", "SM": "N1"}]),
        references=("chr1",),
    )

    with pytest.raises(ValueError, match="coordinate-sorted"):
        build_pon(target_path, normal_alignments=[normal], output_path=output_path)

    assert not output_path.exists()


def test_build_pon_rejects_target_records_out_of_header_contig_order(tmp_path) -> None:
    target_path = tmp_path / "lexically-sorted-targets.vcf"
    output_path = tmp_path / "lexically-sorted.pon.bcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "##contig=<ID=chr2>",
                "##contig=<ID=chr19>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr19\t33302346\t.\tA\tT\t.\tPASS\t.",
                "chr2\t25234307\t.\tA\tC\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    normal = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "N1-rg", "SM": "N1"}]),
        references=("chr1", "chr2", "chr19"),
    )

    with pytest.raises(ValueError, match="coordinate-sorted"):
        build_pon(target_path, normal_alignments=[normal], output_path=output_path)

    assert not output_path.exists()


def test_read_pon_metadata_rejects_vcf_artifacts(tmp_path) -> None:
    vcf_path = tmp_path / "not-a-pon.vcf"
    _write_targets(vcf_path)

    with pytest.raises(ValueError, match="PON artifact must be BCF"):
        read_pon_metadata(vcf_path)


def test_annotate_vcf_with_pon_counts_only_case_and_preserves_targets(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    pon_path = tmp_path / "hotspots.pon.bcf"
    output_path = tmp_path / "calls.vcf"
    _write_targets(target_path)
    normal = _normal("N1", [_read("AAAAAAAAAA")])
    build_pon(
        target_path,
        normal_alignments=[normal],
        output_path=pon_path,
    )
    normal.fetch_calls.clear()

    case = FakeAlignmentFile(
        [_read("AAAAATAAAA", read_group="case-rg")],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )
    output_path.write_text("replace me\n", encoding="utf-8")
    annotate_vcf_with_pon(case, pon_path, output_path=output_path, force=True)

    assert normal.fetch_calls == []
    with pysam.VariantFile(str(output_path)) as calls:
        assert tuple(calls.header.samples) == ("CASE",)
        assert "SKUA_PON_AF" not in calls.header.formats
        assert not any(record.key == "SKUA_PON" for record in calls.header.records)
        record = next(iter(calls))
        assert record.id == "hs1"
        assert record.info["HOTSPOT"]
        assert record.info["SKUA_STATUS"] == "ANNOTATED"
        assert record.info["SKUA_PON_SAMPLE_COUNT"] == 1
        assert record.info["SKUA_PON_NON_ALT_FWD"] == 1
        assert record.samples["CASE"]["SKUA_ALT_FWD"] == 1
        assert record.samples["CASE"]["SKUA_USABLE"] == 1
        assert record.samples["CASE"]["SKUA_ARTIFACT_POSTERIOR"] is not None


def test_pon_preserves_and_uses_target_artifact_prior(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    pon_path = tmp_path / "hotspots.pon.bcf"
    output_path = tmp_path / "calls.vcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior probability that the ALT allele is an artifact before Skua evidence">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\ths1\tA\tT\t.\tPASS\tSKUA_ARTIFACT_PRIOR=0.8",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [])],
        output_path=pon_path,
    )

    with pysam.VariantFile(str(pon_path)) as artifact:
        pon_record = next(iter(artifact))
        assert pon_record.info["SKUA_ARTIFACT_PRIOR"][0] == pytest.approx(0.8)

    case = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )
    annotate_vcf_with_pon(case, pon_path, output_path=output_path)

    with pysam.VariantFile(str(output_path)) as calls:
        record = next(iter(calls))
        assert record.info["SKUA_ARTIFACT_PRIOR"][0] == pytest.approx(0.8)
        assert record.samples["CASE"]["SKUA_ARTIFACT_POSTERIOR"] == pytest.approx(0.8)


def test_annotate_vcf_with_pon_uses_input_vcf_records_and_matching_cached_evidence(
    tmp_path,
) -> None:
    target_path = tmp_path / "hotspots.vcf"
    pon_path = tmp_path / "hotspots.pon.bcf"
    input_path = tmp_path / "case-candidates.vcf"
    output_path = tmp_path / "calls.vcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\tpanel-1\tA\tT\t.\tPASS\t.",
                "chr1\t109\tpanel-2\tA\tC\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    input_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=CALLER_SCORE,Number=1,Type=Integer,Description="Caller score">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t109\tcase-call\tA\tC\t42\tPASS\tCALLER_SCORE=7",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [_read("AAAAAAAAAA")])],
        output_path=pon_path,
    )
    case = FakeAlignmentFile(
        [_read("AAAAAAAACA", read_group="case-rg")],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )

    annotate_vcf_with_pon(
        case,
        pon_path,
        vcf_path=input_path,
        output_path=output_path,
    )

    assert case.fetch_calls == [("chr1", 108, 109)]
    with pysam.VariantFile(str(output_path)) as calls:
        records = list(calls)
        assert len(records) == 1
        record = records[0]
        assert record.id == "case-call"
        assert record.qual == 42
        assert record.info["CALLER_SCORE"] == 7
        assert record.info["SKUA_PON_SAMPLE_COUNT"] == 1
        assert record.info["SKUA_PON_NON_ALT_FWD"] == 1
        assert record.samples["CASE"]["SKUA_ALT_FWD"] == 1


def test_input_vcf_artifact_prior_owns_cached_annotation_precedence(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    pon_path = tmp_path / "hotspots.pon.bcf"
    input_path = tmp_path / "case-candidates.vcf"
    output_path = tmp_path / "calls.vcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior probability that the ALT allele is an artifact before Skua evidence">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\tpanel\tA\tT\t.\tPASS\tSKUA_ARTIFACT_PRIOR=0.8",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    input_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior probability that the ALT allele is an artifact before Skua evidence">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\texplicit\tA\tT\t.\tPASS\tSKUA_ARTIFACT_PRIOR=0.2",
                "chr1\t106\tmissing\tA\tT\t.\tPASS\tSKUA_ARTIFACT_PRIOR=.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [])],
        output_path=pon_path,
    )
    case = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )

    annotate_vcf_with_pon(
        case,
        pon_path,
        vcf_path=input_path,
        output_path=output_path,
        prior_artifact_probability=0.4,
    )

    with pysam.VariantFile(str(output_path)) as calls:
        records = list(calls)
    assert [record.info["SKUA_ARTIFACT_PRIOR"][0] for record in records] == pytest.approx(
        [0.2, 0.4]
    )
    assert [
        record.samples["CASE"]["SKUA_ARTIFACT_POSTERIOR"] for record in records
    ] == pytest.approx([0.2, 0.4])


def test_annotate_vcf_with_pon_rejects_input_variant_missing_from_pon_before_output(
    tmp_path,
) -> None:
    target_path = tmp_path / "hotspots.vcf"
    pon_path = tmp_path / "hotspots.pon.bcf"
    input_path = tmp_path / "case-candidates.vcf"
    output_path = tmp_path / "calls.vcf"
    _write_targets(target_path)
    input_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t107\tcase-call\tA\tG\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [])],
        output_path=pon_path,
    )
    case = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )

    with pytest.raises(ValueError, match=r"chr1:107 A>G.*not present in the PON"):
        annotate_vcf_with_pon(
            case,
            pon_path,
            vcf_path=input_path,
            output_path=output_path,
        )

    assert not output_path.exists()


def test_precomputed_pon_matches_live_normal_annotation(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    pon_path = tmp_path / "hotspots.pon.bcf"
    live_output = tmp_path / "live.vcf"
    cached_output = tmp_path / "cached.vcf"
    _write_targets(target_path)
    normals = [
        _normal("N1", [_read("AAAAAAAAAA")]),
        _normal("N2", [_read("AAAAATAAAA", reverse=True)]),
    ]
    case = FakeAlignmentFile(
        [_read("AAAAATAAAA", read_group="case-rg")],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )

    annotate_vcf_with_normals(
        case,
        target_path,
        normal_alignments=normals,
        output_path=live_output,
    )
    build_pon(target_path, normal_alignments=normals, output_path=pon_path)
    annotate_vcf_with_pon(case, pon_path, output_path=cached_output)

    with (
        pysam.VariantFile(str(live_output)) as live_vcf,
        pysam.VariantFile(str(cached_output)) as cached_vcf,
    ):
        live = next(iter(live_vcf))
        cached = next(iter(cached_vcf))
        for field in (
            "SKUA_PON_SAMPLE_COUNT",
            "SKUA_PON_ALT_FWD",
            "SKUA_PON_ALT_REV",
            "SKUA_PON_NON_ALT_FWD",
            "SKUA_PON_NON_ALT_REV",
            "SKUA_PON_USABLE",
            "SKUA_PON_UNUSABLE",
            "SKUA_PON_DISPERSION_FACTOR",
        ):
            assert cached.info[field] == live.info[field]
        for field in (
            "SKUA_ALT_FWD",
            "SKUA_ALT_REV",
            "SKUA_NON_ALT_FWD",
            "SKUA_NON_ALT_REV",
            "SKUA_USABLE",
            "SKUA_UNUSABLE",
            "SKUA_LOG_BAYES_FACTOR",
            "SKUA_ARTIFACT_POSTERIOR",
        ):
            assert cached.samples["CASE"][field] == live.samples["CASE"][field]


def test_annotate_vcf_with_pon_uses_artifact_evidence_thresholds(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    pon_path = tmp_path / "hotspots.pon.bcf"
    output_path = tmp_path / "calls.vcf"
    _write_targets(target_path)
    build_pon(
        target_path,
        normal_alignments=[_normal("N1", [])],
        output_path=pon_path,
        min_baseq=25,
        min_mapq=30,
    )
    case_read = _read("AAAAATAAAA", read_group="case-rg")
    case_read.mapping_quality = 29
    case = FakeAlignmentFile(
        [case_read],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )

    annotate_vcf_with_pon(case, pon_path, output_path=output_path)

    with pysam.VariantFile(str(output_path)) as calls:
        record = next(iter(calls))
        assert record.samples["CASE"]["SKUA_USABLE"] == 0
        assert record.samples["CASE"]["SKUA_UNUSABLE"] == 1


def test_build_pon_rejects_duplicate_normal_sample_names(tmp_path) -> None:
    target_path = tmp_path / "hotspots.vcf"
    _write_targets(target_path)

    with pytest.raises(ValueError, match="sample names must be unique"):
        build_pon(
            target_path,
            normal_alignments=[_normal("N1", []), _normal("N1", [])],
            output_path=tmp_path / "unused.bcf",
        )


def test_build_pon_rejects_duplicate_target_alleles_before_writing(tmp_path) -> None:
    target_path = tmp_path / "duplicate-targets.vcf"
    output_path = tmp_path / "unused.bcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\tfirst\tA\tT\t.\tPASS\t.",
                "chr1\t106\tsecond\tA\tT\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"duplicate target allele.*chr1:106 A>T"):
        build_pon(
            target_path,
            normal_alignments=[_normal("N1", [])],
            output_path=output_path,
        )

    assert not output_path.exists()


def test_build_pon_rejects_invalid_artifact_prior_before_writing(tmp_path) -> None:
    target_path = tmp_path / "invalid-prior.vcf"
    output_path = tmp_path / "unused.bcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                '##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior probability that the ALT allele is an artifact before Skua evidence">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT\t.\tPASS\tSKUA_ARTIFACT_PRIOR=abc",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"SKUA_ARTIFACT_PRIOR.*malformed"):
        build_pon(
            target_path,
            normal_alignments=[_normal("N1", [])],
            output_path=output_path,
        )

    assert not output_path.exists()
    assert not (tmp_path / "unused.bcf.csi").exists()


@pytest.mark.parametrize(
    ("ref", "alt", "status"),
    [
        ("A", "T,C", "UNSUPPORTED_MULTIALLELIC"),
        ("A", "A", "UNSUPPORTED_COMPLEX_ALLELE"),
        ("A", "CT", "UNSUPPORTED_COMPLEX_ALLELE"),
        ("AT", "C", "UNSUPPORTED_COMPLEX_ALLELE"),
    ],
)
def test_build_pon_rejects_unsupported_target_before_writing(
    tmp_path,
    ref: str,
    alt: str,
    status: str,
) -> None:
    target_path = tmp_path / "unsupported.vcf"
    output_path = tmp_path / "unused.bcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                f"chr1\t106\t.\t{ref}\t{alt}\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=status):
        build_pon(
            target_path,
            normal_alignments=[_normal("N1", [])],
            output_path=output_path,
        )
    assert not output_path.exists()
    assert not (tmp_path / "unused.bcf.csi").exists()
