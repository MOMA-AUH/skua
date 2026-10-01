import os
from pathlib import Path

import pysam
import pytest

from skua import annotate_vcf, annotate_vcf_with_normals, annotate_vcf_with_pon, build_pon
from tests.helpers import FakeAlignmentFile, FakeAlignmentHeader


def _write_vcf(path: Path, pos1: int) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{pos1}\t.\tA\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )


@pytest.fixture(params=[("tbi",), ("csi",), ("tbi", "csi")])
def indexed_output(tmp_path, request):
    old_vcf = tmp_path / "old.vcf"
    output_path = tmp_path / "annotated.vcf.gz"
    _write_vcf(old_vcf, 101)
    pysam.tabix_compress(str(old_vcf), str(output_path))
    indexes = []
    for suffix in request.param:
        pysam.tabix_index(str(output_path), preset="vcf", csi=suffix == "csi")
        indexes.append(Path(f"{output_path}.{suffix}"))

    input_path = tmp_path / "input.vcf"
    _write_vcf(input_path, 900006)
    alignment = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )
    return alignment, input_path, output_path, indexes


@pytest.mark.parametrize("mode", ["case-only", "live-pon", "cached-pon"])
def test_forced_vcf_replacement_retires_companion_indexes(indexed_output, mode) -> None:
    alignment, input_path, output_path, _indexes = indexed_output

    if mode == "case-only":
        annotate_vcf(alignment, input_path, output_path=output_path, force=True)
    else:
        normal = FakeAlignmentFile(
            [],
            header=FakeAlignmentHeader([{"ID": "normal-rg", "SM": "NORMAL"}]),
            references=("chr1",),
        )
        if mode == "live-pon":
            annotate_vcf_with_normals(
                alignment,
                input_path,
                normal_alignments=[normal],
                output_path=output_path,
                force=True,
            )
        else:
            pon_path = output_path.parent / "panel.bcf"
            build_pon(input_path, normal_alignments=[normal], output_path=pon_path)
            annotate_vcf_with_pon(
                alignment, pon_path, output_path=output_path, force=True
            )

    assert not Path(f"{output_path}.tbi").exists()
    assert not Path(f"{output_path}.csi").exists()
    with pysam.VariantFile(str(output_path)) as output:
        records = list(output)
        assert [record.pos for record in records] == [900006]
        assert records[0].info["SKUA_STATUS"] == "ANNOTATED"
        with pytest.raises(ValueError, match="fetch requires an index"):
            list(output.fetch("chr1", 900000, 900010))
    assert list(output_path.parent.glob(".*annotated.vcf.gz*")) == []


@pytest.mark.parametrize("failure", ["backup", "remove-index", "replace-output"])
def test_failed_vcf_publication_restores_output_and_indexes(
    indexed_output, monkeypatch, failure
) -> None:
    alignment, input_path, output_path, indexes = indexed_output
    original_bytes = {path: path.read_bytes() for path in (output_path, *indexes)}
    original_link = os.link
    original_unlink = Path.unlink
    original_replace = os.replace

    def fail_backup(source, target, **kwargs):
        if Path(source) == indexes[-1]:
            raise OSError("forced backup failure")
        return original_link(source, target, **kwargs)

    def fail_index_removal(path, *args, **kwargs):
        if path == indexes[-1]:
            raise OSError("forced index removal failure")
        return original_unlink(path, *args, **kwargs)

    def fail_output_replacement(source, target):
        if Path(target) == output_path:
            assert all(not path.exists() for path in indexes)
            raise OSError("forced output replacement failure")
        return original_replace(source, target)

    if failure == "backup":
        monkeypatch.setattr(os, "link", fail_backup)
    elif failure == "remove-index":
        monkeypatch.setattr(Path, "unlink", fail_index_removal)
    else:
        monkeypatch.setattr(os, "replace", fail_output_replacement)

    with pytest.raises(OSError, match="forced .* failure"):
        annotate_vcf(alignment, input_path, output_path=output_path, force=True)

    assert {path: path.read_bytes() for path in original_bytes} == original_bytes
    for index in indexes:
        with pysam.VariantFile(str(output_path), index_filename=str(index)) as output:
            assert [record.pos for record in output] == [101]
            assert [record.pos for record in output.fetch("chr1", 100, 110)] == [101]
    assert list(output_path.parent.glob(".*annotated.vcf.gz*")) == []


@pytest.mark.parametrize("output_present", [True, False])
def test_vcf_publication_requires_force_for_existing_output_or_indexes(
    indexed_output, output_present
) -> None:
    alignment, input_path, output_path, indexes = indexed_output
    if not output_present:
        output_path.unlink()
    paths = (output_path, *indexes) if output_present else indexes
    original_bytes = {path: path.read_bytes() for path in paths}

    with pytest.raises(FileExistsError, match="already exists"):
        annotate_vcf(alignment, input_path, output_path=output_path)

    assert output_path.exists() == output_present
    assert {path: path.read_bytes() for path in original_bytes} == original_bytes
    assert list(output_path.parent.glob(".*annotated.vcf.gz*")) == []


def test_annotation_failure_preserves_existing_output_and_indexes(
    indexed_output, monkeypatch
) -> None:
    alignment, input_path, output_path, indexes = indexed_output
    original_bytes = {path: path.read_bytes() for path in (output_path, *indexes)}

    def fail_fetch(*args, **kwargs):
        raise OSError("forced alignment read failure")

    monkeypatch.setattr(alignment, "fetch", fail_fetch)

    with pytest.raises(OSError, match="forced alignment read failure"):
        annotate_vcf(alignment, input_path, output_path=output_path, force=True)

    assert {path: path.read_bytes() for path in original_bytes} == original_bytes
    assert list(output_path.parent.glob(".*annotated.vcf.gz*")) == []


@pytest.mark.parametrize("force", [False, True])
def test_new_vcf_output_is_not_automatically_indexed(tmp_path, force) -> None:
    input_path = tmp_path / "input.vcf"
    output_path = tmp_path / "annotated.vcf.gz"
    _write_vcf(input_path, 900006)
    alignment = FakeAlignmentFile(
        [],
        header=FakeAlignmentHeader([{"ID": "case-rg", "SM": "CASE"}]),
        references=("chr1",),
    )

    annotate_vcf(alignment, input_path, output_path=output_path, force=force)

    with pysam.VariantFile(str(output_path)) as output:
        assert [record.pos for record in output] == [900006]
    assert not Path(f"{output_path}.tbi").exists()
    assert not Path(f"{output_path}.csi").exists()
