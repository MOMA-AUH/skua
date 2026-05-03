import skua.cli as cli


def test_main_verify_prints_payload_without_output_path(monkeypatch, capsys) -> None:
    calls: list[dict[str, object]] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str) -> None:
            self.path = path
            self.mode = mode

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

    def fake_verify(alignment_file, vcf_path, **kwargs):
        calls.append(
            {
                "alignment_path": alignment_file.path,
                "alignment_mode": alignment_file.mode,
                "vcf_path": str(vcf_path),
                **kwargs,
            }
        )
        return "{\"ok\": true}"

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(cli, "verify_snv_vcf_to_json", fake_verify)

    exit_code = cli.main(["verify", "--vcf", "input.vcf", "--alignment", "reads.bam"])

    assert exit_code == 0
    assert calls == [
        {
            "alignment_path": "reads.bam",
            "alignment_mode": "rb",
            "vcf_path": "input.vcf",
            "output_path": None,
            "min_baseq": 20,
            "min_mapq": 20,
        }
    ]
    assert capsys.readouterr().out == "{\"ok\": true}\n"


def test_main_verify_uses_output_path_and_does_not_print(monkeypatch, capsys) -> None:
    calls: list[dict[str, object]] = []

    class FakeAlignmentFile:
        def __init__(self, path: str, mode: str) -> None:
            self.path = path
            self.mode = mode

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

    def fake_verify(alignment_file, vcf_path, **kwargs):
        calls.append(
            {
                "alignment_path": alignment_file.path,
                "alignment_mode": alignment_file.mode,
                "vcf_path": str(vcf_path),
                **kwargs,
            }
        )
        return "{\"ok\": true}"

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(cli, "verify_snv_vcf_to_json", fake_verify)

    exit_code = cli.main(
        [
            "verify",
            "--vcf",
            "input.vcf",
            "--alignment",
            "reads.bam",
            "--output",
            "out.json",
            "--min-baseq",
            "15",
            "--min-mapq",
            "12",
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "alignment_path": "reads.bam",
            "alignment_mode": "rb",
            "vcf_path": "input.vcf",
            "output_path": "out.json",
            "min_baseq": 15,
            "min_mapq": 12,
        }
    ]
    assert capsys.readouterr().out == ""


def test_main_verify_accepts_alignment_path_for_cram(monkeypatch, capsys) -> None:
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

    def fake_verify(alignment_file, vcf_path, **kwargs):
        calls.append(
            {
                "alignment_path": alignment_file.path,
                "alignment_mode": alignment_file.mode,
                "alignment_kwargs": alignment_file.kwargs,
                "vcf_path": str(vcf_path),
                **kwargs,
            }
        )
        return "{\"ok\": true}"

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(cli, "verify_snv_vcf_to_json", fake_verify)

    exit_code = cli.main(
        [
            "verify",
            "--vcf",
            "input.vcf",
            "--alignment",
            "reads.cram",
            "--reference",
            "ref.fa",
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
        }
    ]
    assert capsys.readouterr().out == "{\"ok\": true}\n"


def test_main_verify_requires_reference_for_cram(capsys) -> None:
    try:
        cli.main(
            [
                "verify",
                "--vcf",
                "input.vcf",
                "--alignment",
                "reads.cram",
            ]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected SystemExit for missing CRAM reference")

    assert "--reference is required for CRAM input" in capsys.readouterr().err


def test_main_verify_prints_tsv_when_format_tsv(monkeypatch, capsys) -> None:
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

    def fake_verify_tsv(alignment_file, vcf_path, **kwargs):
        calls.append(
            {
                "alignment_path": alignment_file.path,
                "alignment_mode": alignment_file.mode,
                "alignment_kwargs": alignment_file.kwargs,
                "vcf_path": str(vcf_path),
                **kwargs,
            }
        )
        return "col1\tcol2\nA\tB\n"

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(cli, "verify_snv_vcf_to_tsv", fake_verify_tsv)

    exit_code = cli.main(
        [
            "verify",
            "--vcf",
            "input.vcf",
            "--alignment",
            "reads.bam",
            "--format",
            "tsv",
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "alignment_path": "reads.bam",
            "alignment_mode": "rb",
            "alignment_kwargs": {},
            "vcf_path": "input.vcf",
            "output_path": None,
            "min_baseq": 20,
            "min_mapq": 20,
        }
    ]
    assert capsys.readouterr().out == "col1\tcol2\nA\tB\n\n"


def test_main_verify_with_normal_uses_pon_functions(monkeypatch, capsys, tmp_path) -> None:
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
        return "{\"ok\": true}"

    monkeypatch.setattr(cli.pysam, "AlignmentFile", FakeAlignmentFile)
    monkeypatch.setattr(cli, "verify_snv_vcf_to_json_with_normals", fake_verify_with_normals)

    normal_list_path = tmp_path / "normals.txt"
    normal_list_path.write_text("normal1.bam\nnormal2.bam\n", encoding="utf-8")

    exit_code = cli.main(
        [
            "verify",
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
        }
    ]
    assert capsys.readouterr().out == "{\"ok\": true}\n"
