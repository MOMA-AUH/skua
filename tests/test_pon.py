import pysam
import pytest

from skua import (
    annotate_vcf_with_normals,
    annotate_vcf_with_pon,
    build_pon,
    read_pon_evidence,
    read_pon_metadata,
)
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
    assert metadata.evidence_policy_version == 1
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
    annotate_vcf_with_pon(case, pon_path, output_path=output_path)

    assert normal.fetch_calls == []
    with pysam.VariantFile(str(output_path)) as calls:
        assert tuple(calls.header.samples) == ("CASE",)
        record = next(iter(calls))
        assert record.id == "hs1"
        assert record.info["HOTSPOT"]
        assert record.info["SKUA_STATUS"] == "ANNOTATED"
        assert record.info["SKUA_PON_SAMPLE_COUNT"] == 1
        assert record.info["SKUA_PON_NON_ALT_FWD"] == 1
        assert record.samples["CASE"]["SKUA_ALT_FWD"] == 1
        assert record.samples["CASE"]["SKUA_USABLE"] == 1
        assert record.samples["CASE"]["SKUA_ARTIFACT_POSTERIOR"] is not None


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


def test_build_pon_rejects_unsupported_target_before_writing(tmp_path) -> None:
    target_path = tmp_path / "unsupported.vcf"
    output_path = tmp_path / "unused.bcf"
    target_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t106\t.\tA\tT,C\t.\tPASS\t.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="UNSUPPORTED_MULTIALLELIC"):
        build_pon(
            target_path,
            normal_alignments=[_normal("N1", [])],
            output_path=output_path,
        )
    assert not output_path.exists()
