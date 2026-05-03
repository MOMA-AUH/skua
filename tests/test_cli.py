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

    exit_code = cli.main(["verify", "--vcf", "input.vcf", "--bam", "reads.bam"])

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
            "--bam",
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
