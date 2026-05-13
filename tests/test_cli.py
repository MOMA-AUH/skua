import skua.cli as cli


def test_main_annotate_requires_normal_list(capsys) -> None:
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
        raise AssertionError("Expected SystemExit for missing --normal-list")

    assert "the following arguments are required: --normal-list" in capsys.readouterr().err


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
        return "##fileformat=VCFv4.2\n"

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(
        cli,
        "annotate_snv_vcf_with_normals",
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
            "output_path": None,
            "min_baseq": 20,
            "min_mapq": 20,
            "truncate": 0.1,
            "prior_variant_probability": 0.5,
        }
    ]
    assert capsys.readouterr().out == "##fileformat=VCFv4.2\n"


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
        return "##fileformat=VCFv4.2\n"

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(
        cli,
        "annotate_snv_vcf_with_normals",
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
            "min_baseq": 15,
            "min_mapq": 12,
            "truncate": 0.1,
            "prior_variant_probability": 0.5,
        }
    ]
    assert capsys.readouterr().out == ""


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
        return "##fileformat=VCFv4.2\n"

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(
        cli,
        "annotate_snv_vcf_with_normals",
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
            "output_path": None,
            "min_baseq": 20,
            "min_mapq": 20,
            "truncate": 0.1,
            "prior_variant_probability": 0.5,
            "normal_count": 1,
        }
    ]
    assert capsys.readouterr().out == "##fileformat=VCFv4.2\n"


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


def test_main_annotate_rejects_removed_output_format_flag(capsys, tmp_path) -> None:
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
                "--output-format",
                "vcf",
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected SystemExit for removed --output-format")

    assert "unrecognized arguments: --output-format vcf" in capsys.readouterr().err


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
