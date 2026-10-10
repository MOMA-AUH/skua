"""Exercise an installed distribution, independently of the source checkout.

Run with ``python -I tests/artifact_smoke.py --expected-version VERSION``.
All inputs are small synthetic, coordinate-sorted, indexed files. The case has
four paired fragments (two ALT, two reference); the normal has ten reference
fragments. Both mates overlap the first two targets and must count once per
fragment. A third target has no coverage and must have missing model scores.
"""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from typing import Any

import pysam
import skua


def cli(*arguments: str) -> str:
    executable = Path(sys.executable).parent / "skua"
    result = subprocess.run([str(executable), *arguments], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"skua {arguments} failed:\n{result.stdout}\n{result.stderr}")
    return result.stdout


def alignment(root: Path, sample: str, reference: Path, suffix: str) -> Path:
    path = root / f"{sample}.{suffix}"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 200}],
        "RG": [{"ID": sample, "SM": sample}],
    }
    with pysam.AlignmentFile(
        str(path), "wc" if suffix == "cram" else "wb", header=header,
        reference_filename=str(reference),
    ) as output:
        for index in range(4 if sample == "CASE" else 10):
            for mate in range(2):
                read = pysam.AlignedSegment(output.header)
                read.query_name = f"{sample}-{index}"
                sequence = list("A" * 40)
                if sample == "CASE" and index < 2:
                    sequence[10] = "T"
                read.query_sequence = "".join(sequence)
                read.query_qualities = pysam.qualitystring_to_array("I" * 40)
                read.flag = ((99, 147) if index % 2 == 0 else (83, 163))[mate]
                read.reference_id = read.next_reference_id = 0
                read.reference_start = read.next_reference_start = 10
                read.template_length = 40 if mate == 0 else -40
                read.mapping_quality = 60
                read.cigarstring = "40M"
                read.set_tag("RG", sample)
                output.write(read)
    pysam.index(str(path))
    return path


def annotation_values(path: Path) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    with pysam.VariantFile(str(path)) as source:
        return [(dict(record.info), dict(record.samples["CASE"])) for record in source]


def exercise(root: Path, reference: Path, targets: Path, suffix: str) -> list:
    case = alignment(root, "CASE", reference, suffix)
    normal = alignment(root, "NORMAL", reference, suffix)
    normal_list = root / f"normals-{suffix}.txt"
    normal_list.write_text(f"{normal}\n", encoding="utf-8")
    panel = root / f"panel-{suffix}.bcf"
    cli("pon", "build", "--vcf", str(targets), "--normal-list", str(normal_list),
        "--reference", str(reference), "--output", str(panel))
    validation = json.loads(cli("pon", "validate", str(panel), "--reference", str(reference),
                                "--targets", str(targets), "--json"))
    assert validation["valid"], validation
    results = []
    for mode, source_args in (
        ("direct", ["--vcf", str(targets), "--normal-list", str(normal_list)]),
        ("cached", ["--pon", str(panel)]),
        ("subset", ["--vcf", str(targets), "--pon", str(panel)]),
    ):
        output = root / f"{suffix}-{mode}.vcf"
        cli("annotate", "--alignment", str(case), "--reference", str(reference),
            "--output", str(output), *source_args)
        values = annotation_values(output)
        assert len(values) == 3
        info, sample = values[0]
        assert sample["SKUA_ALT_FWD"] == sample["SKUA_ALT_REV"] == 1
        assert sample["SKUA_NON_ALT_FWD"] == sample["SKUA_NON_ALT_REV"] == 1
        assert sample["SKUA_USABLE"] == 4
        assert sample["SKUA_ASSESSMENT_STATUS"] == "ASSESSED"
        assert info["SKUA_PON_USABLE"] == 10
        assert info["SKUA_PON_SAMPLE_COUNT"] == 1
        assert sample["SKUA_ARTIFACT_POSTERIOR"] is not None
        assert values[1][1]["SKUA_ALT_FWD"] == values[1][1]["SKUA_ALT_REV"] == 0
        uncovered_info, uncovered_sample = values[2]
        assert uncovered_info["SKUA_STATUS"] == "ANNOTATED"
        assert uncovered_info["SKUA_PON_USABLE"] == uncovered_sample["SKUA_USABLE"] == 0
        assert uncovered_sample["SKUA_ASSESSMENT_STATUS"] == "INSUFFICIENT_EVIDENCE"
        assert uncovered_sample["SKUA_ASSESSMENT_REASONS"] == ("CASE_DEPTH", "NORMAL_DEPTH")
        assert uncovered_sample["SKUA_ARTIFACT_POSTERIOR"] is None
        assert uncovered_sample["SKUA_LOG_BAYES_FACTOR"] is None
        results.append(values)
    assert results[0] == results[1] == results[2], "Direct/cached evidence or scores differ"
    return results[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-version", required=True)
    args = parser.parse_args()
    package_path = Path(skua.__file__).resolve()
    assert package_path.is_relative_to(Path(sys.prefix).resolve()), (
        f"Skua imported outside the installed environment: {package_path}"
    )
    assert skua.__version__ == version("skua") == cli("--version").strip() == args.expected_version
    with TemporaryDirectory(prefix="skua-artifact-") as directory:
        root = Path(directory)
        reference = root / "reference.fa"
        reference.write_text(">chr1\n" + "A" * 200 + "\n", encoding="utf-8")
        pysam.faidx(str(reference))
        targets = root / "targets.vcf"
        targets.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t21\t.\tA\tT\t.\tPASS\t.\n"
            "chr1\t25\t.\tA\tG\t.\tPASS\t.\n"
            "chr1\t151\t.\tA\tT\t.\tPASS\t.\n", encoding="utf-8",
        )
        bam_result = exercise(root, reference, targets, "bam")
        cram_result = exercise(root, reference, targets, "cram")
        assert bam_result == cram_result, "BAM/CRAM evidence or scores differ"
    print(f"Installed Skua {args.expected_version}: BAM/CRAM and direct/cached parity passed")


if __name__ == "__main__":
    main()
