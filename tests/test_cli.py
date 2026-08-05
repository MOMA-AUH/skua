from pathlib import Path

import skua.cli as cli
from skua import __version__


def test_main_version_prints_version_and_exits_successfully(capsys) -> None:
    try:
        cli.main(["--version"])
    except SystemExit as exc:
        assert exc.code == 0
    else:
        raise AssertionError("Expected SystemExit for --version")

    assert capsys.readouterr().out == f"{__version__}\n"


def test_main_annotate_requires_normal_list_or_precomputed_pon(capsys) -> None:
    try:
        cli.main(
            [
                "annotate",
                "--vcf",
                "input.vcf",
                "--alignment",
                "reads.bam",
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected SystemExit for missing PON source")

    assert "one of the arguments --normal-list --pon is required" in capsys.readouterr().err


def test_main_annotate_with_precomputed_pon_counts_only_case(monkeypatch, tmp_path) -> None:
    opened_paths: list[str] = []
    calls: list[dict[str, object]] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str, **kwargs) -> None:
            self.path = path
            opened_paths.append(path)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

    def fake_annotate(alignment_file, pon_path, **kwargs):
        calls.append(
            {
                "case_path": alignment_file.path,
                "pon_path": str(pon_path),
                **kwargs,
            }
        )

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(cli, "annotate_vcf_with_pon", fake_annotate)

    assert cli.main(
        [
            "annotate",
            "--alignment",
            "case.bam",
            "--pon",
            "hotspots.pon.bcf",
            "--output",
            "calls.vcf.gz",
        ]
    ) == 0

    assert opened_paths == ["case.bam"]
    assert calls == [
        {
            "case_path": "case.bam",
            "pon_path": "hotspots.pon.bcf",
            "output_path": "calls.vcf.gz",
            "sample_name": None,
            "reference_path": None,
            "min_baseq": None,
            "min_mapq": None,
            "truncate": 0.1,
            "prior_variant_probability": 0.5,
        }
    ]


def test_main_pon_build_opens_normals_and_writes_bcf(monkeypatch, tmp_path) -> None:
    calls: list[dict[str, object]] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str, **kwargs) -> None:
            self.path = path

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

    def fake_build(vcf_path, **kwargs):
        calls.append(
            {
                "vcf_path": str(vcf_path),
                "normal_paths": [normal.path for normal in kwargs["normal_alignments"]],
                **{key: value for key, value in kwargs.items() if key != "normal_alignments"},
            }
        )

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(cli, "build_pon", fake_build)
    normal_list = tmp_path / "normals.lst"
    normal_list.write_text("normal1.bam\nnormal2.bam\n", encoding="utf-8")

    assert cli.main(
        [
            "pon",
            "build",
            "--vcf",
            "hotspots.vcf.gz",
            "--normal-list",
            str(normal_list),
            "--output",
            "hotspots.pon.bcf",
            "--min-baseq",
            "25",
        ]
    ) == 0

    assert calls == [
        {
            "vcf_path": "hotspots.vcf.gz",
            "normal_paths": ["normal1.bam", "normal2.bam"],
            "output_path": Path("hotspots.pon.bcf"),
            "reference_path": None,
            "min_baseq": 25,
            "min_mapq": 20,
        }
    ]


def test_main_annotate_with_normal_uses_pon_functions(monkeypatch, capsys, tmp_path) -> None:
    calls: list[dict[str, object]] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str, **kwargs) -> None:
            self.path = path
            self.mode = mode
            self.kwargs = kwargs

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def close(self) -> None:
            pass

    def fake_verify_with_normals(alignment_file, vcf_path, **kwargs):
        calls.append(
            {
                "case_path": alignment_file.path,
                "vcf_path": str(vcf_path),
                "normal_count": len(kwargs.get("normal_alignments", [])),
                **{k: v for k, v in kwargs.items() if k != "normal_alignments"},
            }
        )
        return None

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(
        cli,
        "annotate_vcf_with_normals",
        fake_verify_with_normals,
    )

    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.bam\nnormal2.bam\n", encoding="utf-8")

    exit_code = cli.main(
        [
            "annotate",
            "--vcf",
            "input.vcf",
            "--alignment",
            "case.bam",
            "--normal-list",
            str(normal_list_path),
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "case_path": "case.bam",
            "vcf_path": "input.vcf",
            "normal_count": 2,
            "output_path": "-",
            "sample_name": None,
            "reference_path": None,
            "strict": False,
            "min_baseq": 20,
            "min_mapq": 20,
            "truncate": 0.1,
            "prior_variant_probability": 0.5,
        }
    ]
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_main_annotate_with_normal_uses_output_path_and_does_not_print(
    monkeypatch,
    capsys,
    tmp_path,
) -> None:
    calls: list[dict[str, object]] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str, **kwargs) -> None:
            self.path = path
            self.mode = mode
            self.kwargs = kwargs

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def close(self) -> None:
            pass

    def fake_verify_with_normals(alignment_file, vcf_path, **kwargs):
        calls.append(
            {
                "case_path": alignment_file.path,
                "vcf_path": str(vcf_path),
                "normal_count": len(kwargs.get("normal_alignments", [])),
                **{k: v for k, v in kwargs.items() if k != "normal_alignments"},
            }
        )
        return None

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(
        cli,
        "annotate_vcf_with_normals",
        fake_verify_with_normals,
    )

    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.bam\n", encoding="utf-8")

    exit_code = cli.main(
        [
            "annotate",
            "--vcf",
            "input.vcf",
            "--alignment",
            "case.bam",
            "--normal-list",
            str(normal_list_path),
            "--output",
            "out.vcf.gz",
            "--min-baseq",
            "15",
            "--min-mapq",
            "12",
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "case_path": "case.bam",
            "vcf_path": "input.vcf",
            "normal_count": 1,
            "output_path": "out.vcf.gz",
            "sample_name": None,
            "reference_path": None,
            "strict": False,
            "min_baseq": 15,
            "min_mapq": 12,
            "truncate": 0.1,
            "prior_variant_probability": 0.5,
        }
    ]
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_main_annotate_forwards_requested_sample(monkeypatch, tmp_path) -> None:
    calls: list[dict[str, object]] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str, **kwargs) -> None:
            self.path = path

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

    def fake_annotate(alignment_file, vcf_path, **kwargs):
        calls.append(kwargs)
        return None

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(cli, "annotate_vcf_with_normals", fake_annotate)
    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.bam\n", encoding="utf-8")

    assert cli.main(
        [
            "annotate",
            "--vcf",
            "input.vcf",
            "--alignment",
            "case.bam",
            "--normal-list",
            str(normal_list_path),
            "--sample",
            "CASE",
        ]
    ) == 0

    assert calls[0]["sample_name"] == "CASE"


def test_main_annotate_accepts_alignment_path_for_cram(monkeypatch, capsys, tmp_path) -> None:
    calls: list[dict[str, object]] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str, **kwargs) -> None:
            self.path = path
            self.mode = mode
            self.kwargs = kwargs

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

        def close(self) -> None:
            pass

    def fake_verify_with_normals(alignment_file, vcf_path, **kwargs):
        calls.append(
            {
                "alignment_path": alignment_file.path,
                "alignment_mode": alignment_file.mode,
                "alignment_kwargs": alignment_file.kwargs,
                "vcf_path": str(vcf_path),
                **{k: v for k, v in kwargs.items() if k != "normal_alignments"},
                "normal_count": len(kwargs.get("normal_alignments", [])),
            }
        )
        return None

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(
        cli,
        "annotate_vcf_with_normals",
        fake_verify_with_normals,
    )

    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.cram\n", encoding="utf-8")

    exit_code = cli.main(
        [
            "annotate",
            "--vcf",
            "input.vcf",
            "--alignment",
            "reads.cram",
            "--reference",
            "ref.fa",
            "--normal-list",
            str(normal_list_path),
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "alignment_path": "reads.cram",
            "alignment_mode": "rb",
            "alignment_kwargs": {"reference_filename": "ref.fa"},
            "vcf_path": "input.vcf",
            "output_path": "-",
            "sample_name": None,
            "strict": False,
            "min_baseq": 20,
            "min_mapq": 20,
            "truncate": 0.1,
            "prior_variant_probability": 0.5,
            "reference_path": "ref.fa",
            "normal_count": 1,
        }
    ]
    assert capsys.readouterr().out == ""


def test_main_annotate_requires_reference_for_cram(capsys, tmp_path) -> None:
    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.bam\n", encoding="utf-8")

    try:
        cli.main(
            [
                "annotate",
                "--vcf",
                "input.vcf",
                "--alignment",
                "reads.cram",
                "--normal-list",
                str(normal_list_path),
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected SystemExit for missing CRAM reference")

    assert "--reference is required for CRAM input" in capsys.readouterr().err


def test_main_annotate_requires_reference_for_cram_in_normal_list(capsys, tmp_path) -> None:
    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.cram\n", encoding="utf-8")

    try:
        cli.main(
            [
                "annotate",
                "--vcf",
                "input.vcf",
                "--alignment",
                "reads.bam",
                "--normal-list",
                str(normal_list_path),
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected SystemExit for missing CRAM reference")

    assert "--reference is required for CRAM input" in capsys.readouterr().err


def test_main_annotate_rejects_empty_normal_list(capsys, tmp_path) -> None:
    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("# comment only\n", encoding="utf-8")

    try:
        cli.main(
            [
                "annotate",
                "--vcf",
                "input.vcf",
                "--alignment",
                "reads.bam",
                "--normal-list",
                str(normal_list_path),
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected SystemExit for empty normal list")

    assert "--normal-list must include at least one normal alignment path" in capsys.readouterr().err


def test_main_annotate_rejects_invalid_parameter_before_opening_files(capsys, tmp_path) -> None:
    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.bam\n", encoding="utf-8")

    try:
        cli.main(
            [
                "annotate",
                "--vcf",
                "input.vcf",
                "--alignment",
                "reads.bam",
                "--normal-list",
                str(normal_list_path),
                "--truncate",
                "0",
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected SystemExit for invalid --truncate")

    assert "--truncate must be greater than 0 and no greater than 1" in capsys.readouterr().err


def test_main_annotate_closes_opened_normals_when_later_open_fails(monkeypatch, tmp_path) -> None:
    opened_files: list[FakeAlignmentFile] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str, **kwargs) -> None:
            if path == "bad-normal.bam":
                raise OSError("failed to open normal")
            self.path = path
            self.mode = mode
            self.kwargs = kwargs
            self.closed = False
            opened_files.append(self)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            self.close()
            return None

        def close(self) -> None:
            self.closed = True

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)

    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("good-normal.bam\nbad-normal.bam\n", encoding="utf-8")

    try:
        cli.main(
            [
                "annotate",
                "--vcf",
                "input.vcf",
                "--alignment",
                "case.bam",
                "--normal-list",
                str(normal_list_path),
            ]
        )
    except OSError as exc:
        assert str(exc) == "failed to open normal"
    else:
        raise AssertionError("Expected OSError when a later normal alignment open fails")

    assert [alignment.path for alignment in opened_files] == ["case.bam", "good-normal.bam"]
    assert all(alignment.closed for alignment in opened_files)


def test_main_annotate_rejects_output_path_with_non_vcf_suffix(capsys, tmp_path) -> None:
    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.bam\n", encoding="utf-8")

    try:
        cli.main(
            [
                "annotate",
                "--vcf",
                "input.vcf",
                "--alignment",
                "reads.bam",
                "--normal-list",
                str(normal_list_path),
                "--output",
                "out.txt",
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected SystemExit for invalid output extension")

    assert "--output must end with .vcf or .vcf.gz" in capsys.readouterr().err
