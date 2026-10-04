"""Command-line interface for skua."""

import argparse
from contextlib import ExitStack
import json
from pathlib import Path
from typing import Any

import pysam

from . import __version__
from .core import (
    _validate_annotation_parameters,
    annotate_vcf_with_normals,
    annotate_vcf_with_pon,
    build_pon,
)
from .pon import PonInspection, inspect_pon, validate_pon
from .stats import AssessmentThresholds


class OptionalDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Show defaults only for non-required options with concrete defaults."""

    def _get_help_string(self, action: argparse.Action) -> str:
        help_text = action.help
        if help_text is None:
            help_text = ""
        if (
            "%(default)" not in help_text
            and action.default is not argparse.SUPPRESS
            and action.default is not None
            and not action.required
        ):
            help_text += " (default: %(default)s)"
        return help_text


def _add_evidence_arguments(
    parser: argparse.ArgumentParser,
    *,
    default: int | None,
) -> None:
    inherited_help = (
        "; defaults to 20 with --normal-list and cannot be used with --pon"
        if default is None
        else ""
    )
    parser.add_argument(
        "--min-baseq",
        type=int,
        default=default,
        help="Minimum base quality" + inherited_help,
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=default,
        help="Minimum mapping quality" + inherited_help,
    )


def _add_model_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--truncate",
        type=float,
        default=0.1,
        help="Truncation threshold for PON sample inclusion",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=None,
        help="Optional pseudocount for beta-binomial rate estimates",
    )
    parser.add_argument(
        "--prior-artifact-probability",
        type=float,
        default=0.5,
        help="Fallback artifact prior when SKUA_ARTIFACT_PRIOR is missing",
    )
    defaults = AssessmentThresholds()
    for name, description in (
        ("min_case_depth", "Minimum usable case depth for assessment"),
        ("min_normal_depth", "Minimum pooled normal usable depth after truncation"),
        ("min_normal_samples", "Minimum retained normals; 0 disables this requirement"),
        ("min_case_strand_depth", "Minimum usable case depth on each strand; 0 disables"),
        ("min_normal_strand_depth", "Minimum pooled normal depth on each strand after truncation; 0 disables"),
    ):
        parser.add_argument(
            "--" + name.replace("_", "-"), type=int,
            default=getattr(defaults, name), help=description,
        )


def _validate_parameters(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    *,
    include_model: bool,
) -> None:
    """Reject invalid parameters before opening input files."""
    try:
        _validate_annotation_parameters(
            min_baseq=args.min_baseq,
            min_mapq=args.min_mapq,
            truncate=args.truncate if include_model else None,
            pseudocount=args.pseudocount if include_model else None,
            prior_artifact_probability=(
                args.prior_artifact_probability if include_model else None
            ),
        )
        if include_model:
            args.assessment_thresholds = AssessmentThresholds(
                min_case_depth=args.min_case_depth,
                min_normal_depth=args.min_normal_depth,
                min_normal_samples=args.min_normal_samples,
                min_case_strand_depth=args.min_case_strand_depth,
                min_normal_strand_depth=args.min_normal_strand_depth,
            )
    except ValueError as exc:
        parser.error(f"--{str(exc).replace('_', '-')}")


def _read_normal_paths(parser: argparse.ArgumentParser, list_path: str) -> list[str]:
    normal_paths = []
    for line in Path(list_path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            normal_paths.append(stripped)
    if not normal_paths:
        parser.error("--normal-list must include at least one normal alignment path")
    return normal_paths


def _require_reference_for_crams(
    parser: argparse.ArgumentParser,
    paths: list[str],
    *,
    reference: str | None,
) -> None:
    if reference is None and any(Path(path).suffix.lower() == ".cram" for path in paths):
        parser.error("--reference is required for CRAM input")


def _alignment_kwargs(reference: str | None) -> dict[str, str]:
    if reference is None:
        return {}
    return {"reference_filename": reference}


def _pon_model_kwargs(args: argparse.Namespace) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "truncate": args.truncate,
        "prior_artifact_probability": args.prior_artifact_probability,
        "assessment_thresholds": args.assessment_thresholds,
    }
    if args.pseudocount is not None:
        kwargs["pseudocount"] = args.pseudocount
    return kwargs


def _validate_vcf_output(parser: argparse.ArgumentParser, output: str | None) -> None:
    if output is None:
        return
    output_path = Path(output)
    if not (
        output_path.suffix.lower() == ".vcf"
        or output_path.name.lower().endswith(".vcf.gz")
    ):
        parser.error("--output must end with .vcf or .vcf.gz")


def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="skua",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
        help="Show the skua version and exit",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate targets with read counts and PON statistics",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    annotate_parser.add_argument(
        "--vcf",
        help="Input target VCF path (required with --normal-list; optional with --pon)",
    )
    annotate_parser.add_argument(
        "--alignment",
        required=True,
        help="Input case BAM/CRAM path (required)",
    )
    annotate_parser.add_argument(
        "--sample",
        help="Case sample name (required when BAM sample matching is ambiguous)",
    )
    annotate_parser.add_argument("--reference", help="Reference FASTA path (required for CRAM)")
    annotate_parser.add_argument("--output", help="Optional output VCF path (.vcf or .vcf.gz)")
    annotate_parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing output and recompute existing Skua annotations",
    )
    pon_source_group = annotate_parser.add_mutually_exclusive_group(required=True)
    pon_source_group.add_argument(
        "--normal-list",
        help="File listing normal BAM/CRAM paths; requires --vcf",
    )
    pon_source_group.add_argument(
        "--pon",
        help="Precomputed PON BCF; its records define targets when --vcf is omitted",
    )
    _add_evidence_arguments(annotate_parser, default=None)
    _add_model_arguments(annotate_parser)
    annotate_parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail when an input VCF record cannot be annotated",
    )

    pon_parser = subparsers.add_parser(
        "pon",
        help="Build and inspect precomputed panels of normals",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    pon_subparsers = pon_parser.add_subparsers(dest="pon_command", required=True)
    pon_build_parser = pon_subparsers.add_parser(
        "build",
        help="Precompute normal evidence for target variants",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    pon_build_parser.add_argument("--vcf", required=True, help="Target VCF path (required)")
    pon_build_parser.add_argument(
        "--normal-list",
        required=True,
        help="File listing normal BAM/CRAM paths, one per line (required)",
    )
    pon_build_parser.add_argument(
        "--output",
        required=True,
        help="Output precomputed PON path ending in .bcf (required)",
    )
    pon_build_parser.add_argument(
        "--reference",
        help="Reference FASTA path (required for CRAM)",
    )
    pon_build_parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing output pair and target Skua annotations",
    )
    _add_evidence_arguments(pon_build_parser, default=20)

    pon_inspect_parser = pon_subparsers.add_parser(
        "inspect",
        help="Show PON header metadata without validating compatibility",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    pon_inspect_parser.add_argument("pon", help="Precomputed PON BCF path")
    pon_inspect_parser.add_argument(
        "--json",
        action="store_true",
        help="Write inspection metadata as JSON",
    )

    pon_validate_parser = pon_subparsers.add_parser(
        "validate",
        help="Validate PON structure and optional reference/target compatibility",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    pon_validate_parser.add_argument("pon", help="Precomputed PON BCF path")
    pon_validate_parser.add_argument(
        "--reference",
        help="Reference FASTA used to validate PON REF alleles",
    )
    pon_validate_parser.add_argument(
        "--targets",
        help="Target VCF expected to define exactly the PON targets",
    )
    pon_validate_parser.add_argument(
        "--json",
        action="store_true",
        help="Write validation result as JSON",
    )

    return parser


def _run_annotate(parser: argparse.ArgumentParser, args: argparse.Namespace) -> int:
    if args.pon is not None and (
        args.min_baseq is not None or args.min_mapq is not None
    ):
        parser.error(
            "--min-baseq and --min-mapq cannot be used with --pon; "
            "the PON defines both thresholds"
        )
    _validate_parameters(parser, args, include_model=True)
    _validate_vcf_output(parser, args.output)
    if args.normal_list is not None and args.vcf is None:
        parser.error("--vcf is required with --normal-list")

    _require_reference_for_crams(
        parser,
        [args.alignment],
        reference=args.reference,
    )
    alignment_kwargs = _alignment_kwargs(args.reference)
    normal_paths: list[str] = []
    if args.normal_list is not None:
        normal_paths = _read_normal_paths(parser, args.normal_list)
        _require_reference_for_crams(parser, normal_paths, reference=args.reference)

    with ExitStack() as stack:
        alignment_file = stack.enter_context(
            pysam.AlignmentFile(args.alignment, "rb", **alignment_kwargs)
        )
        try:
            if args.pon is not None:
                annotate_vcf_with_pon(
                    alignment_file,
                    Path(args.pon),
                    vcf_path=Path(args.vcf) if args.vcf is not None else None,
                    output_path=args.output if args.output is not None else "-",
                    sample_name=args.sample,
                    reference_path=args.reference,
                    strict=args.strict,
                    force=args.force,
                    **_pon_model_kwargs(args),
                )
            else:
                min_baseq = args.min_baseq if args.min_baseq is not None else 20
                min_mapq = args.min_mapq if args.min_mapq is not None else 20
                normal_alignments = [
                    stack.enter_context(
                        pysam.AlignmentFile(path, "rb", **alignment_kwargs)
                    )
                    for path in normal_paths
                ]
                annotate_vcf_with_normals(
                    alignment_file,
                    Path(args.vcf),
                    normal_alignments=normal_alignments,
                    output_path=args.output if args.output is not None else "-",
                    sample_name=args.sample,
                    reference_path=args.reference,
                    strict=args.strict,
                    min_baseq=min_baseq,
                    min_mapq=min_mapq,
                    force=args.force,
                    **_pon_model_kwargs(args),
                )
        except (FileExistsError, ValueError) as exc:
            parser.error(str(exc))
    return 0


def _run_pon_build(parser: argparse.ArgumentParser, args: argparse.Namespace) -> int:
    _validate_parameters(parser, args, include_model=False)
    if Path(args.output).suffix.lower() != ".bcf":
        parser.error("--output must end with .bcf")

    normal_paths = _read_normal_paths(parser, args.normal_list)
    _require_reference_for_crams(parser, normal_paths, reference=args.reference)
    alignment_kwargs = _alignment_kwargs(args.reference)
    with ExitStack() as stack:
        normal_alignments = [
            stack.enter_context(pysam.AlignmentFile(path, "rb", **alignment_kwargs))
            for path in normal_paths
        ]
        try:
            build_pon(
                Path(args.vcf),
                normal_alignments=normal_alignments,
                output_path=Path(args.output),
                reference_path=args.reference,
                min_baseq=args.min_baseq,
                min_mapq=args.min_mapq,
                **({"force": True} if args.force else {}),
            )
        except (FileExistsError, ValueError) as exc:
            parser.error(str(exc))
    return 0


def _format_pon_inspection(inspection: PonInspection) -> str:
    """Render header-level PON metadata for a human terminal."""
    def value_or_missing(value: str | None) -> str:
        return value if value is not None else "<missing>"

    sample_names = ", ".join(inspection.sample_names) or "<none>"
    return "\n".join(
        (
            f"PON: {inspection.path}",
            f"Format: {inspection.format}",
            f"CSI index: {'present' if inspection.index_present else 'missing'}",
            f"PON metadata records: {inspection.metadata_record_count}",
            f"Schema version: {value_or_missing(inspection.schema_version)}",
            "Evidence policy version: "
            f"{value_or_missing(inspection.evidence_policy_version)}",
            f"Minimum base quality: {value_or_missing(inspection.min_baseq)}",
            f"Minimum mapping quality: {value_or_missing(inspection.min_mapq)}",
            f"Skua version: {value_or_missing(inspection.skua_version)}",
            "Reference status: " + (
                inspection.reference_identity.status if inspection.reference_identity else "<missing or invalid>"
            ),
            "Legacy build provenance: " + ("available" if inspection.provenance is not None else "not present or invalid (optional)"),
            f"Normal samples ({len(inspection.sample_names)}): {sample_names}",
        )
    )


def _run_pon_inspect(parser: argparse.ArgumentParser, args: argparse.Namespace) -> int:
    try:
        inspection = inspect_pon(Path(args.pon))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    if args.json:
        print(json.dumps(inspection.as_dict(), indent=2, sort_keys=True))
    else:
        print(_format_pon_inspection(inspection))
    return 0


def _run_pon_validate(_parser: argparse.ArgumentParser, args: argparse.Namespace) -> int:
    result = validate_pon(
        Path(args.pon),
        reference_path=Path(args.reference) if args.reference is not None else None,
        target_vcf_path=Path(args.targets) if args.targets is not None else None,
    )
    if args.json:
        print(json.dumps(result.as_dict(), indent=2, sort_keys=True))
    elif result.valid:
        print(f"PON validation passed: {args.pon}")
    else:
        for error in result.errors:
            print(f"PON validation failed: {error}")
    return 0 if result.valid else 1


def main(argv: list[str] | None = None) -> int:
    """Run the skua CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "annotate":
        return _run_annotate(parser, args)
    if args.command == "pon":
        if args.pon_command == "build":
            return _run_pon_build(parser, args)
        if args.pon_command == "inspect":
            return _run_pon_inspect(parser, args)
        if args.pon_command == "validate":
            return _run_pon_validate(parser, args)

    parser.error(f"Unknown command: {args.command}")
    return 2
