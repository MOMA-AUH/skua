"""Command-line interface for skua."""

import argparse
from contextlib import ExitStack
from pathlib import Path

import pysam

from . import __version__
from .core import (
    annotate_snv_vcf_with_normals,
)


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
        help="Annotate VCF with read counts and PON statistics",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    annotate_parser.add_argument("--vcf", required=True, help="Input VCF path (required)")
    annotate_parser.add_argument("--alignment", required=True, help="Input BAM/CRAM path (required)")
    annotate_parser.add_argument("--reference", help="Reference FASTA path (required for CRAM)")
    annotate_parser.add_argument(
        "--output",
        help="Optional output VCF path (.vcf or .vcf.gz)",
    )
    annotate_parser.add_argument(
        "--normal-list",
        required=True,
        help="Path to file listing normal sample BAM/CRAM paths, one per line (required)",
    )
    annotate_parser.add_argument("--min-baseq", type=int, default=20, help="Minimum base quality")
    annotate_parser.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality")
    annotate_parser.add_argument(
        "--truncate",
        type=float,
        default=0.1,
        help="Truncation threshold for PON sample inclusion",
    )
    annotate_parser.add_argument(
        "--pseudocount",
        type=float,
        default=None,
        help="Optional pseudocount for beta-binomial rate estimates",
    )
    annotate_parser.add_argument(
        "--prior-variant-probability",
        type=float,
        default=0.5,
        help="Prior probability for the variant model",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the skua CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "annotate":
        if args.output is not None:
            output_path = Path(args.output)
            if not (
                output_path.suffix.lower() == ".vcf"
                or output_path.name.lower().endswith(".vcf.gz")
            ):
                parser.error("--output must end with .vcf or .vcf.gz")

        alignment_path = Path(args.alignment)
        if alignment_path.suffix.lower() == ".cram" and args.reference is None:
            parser.error("--reference is required for CRAM input")

        alignment_kwargs: dict[str, str] = {}
        if args.reference is not None:
            alignment_kwargs["reference_filename"] = args.reference

        normal_paths: list[str] = []
        for line in Path(args.normal_list).read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                normal_paths.append(stripped)

        if not normal_paths:
            parser.error("--normal-list must include at least one normal alignment path")

        for normal_path in normal_paths:
            normal_path_obj = Path(normal_path)
            if normal_path_obj.suffix.lower() == ".cram" and args.reference is None:
                parser.error("--reference is required for CRAM input")

        with ExitStack() as stack:
            alignment_file = stack.enter_context(
                pysam.AlignmentFile(args.alignment, "rb", **alignment_kwargs)
            )

            normal_alignments = []
            for normal_path in normal_paths:
                normal_alignments.append(
                    stack.enter_context(
                        pysam.AlignmentFile(normal_path, "rb", **alignment_kwargs)
                    )
                )

            pon_kwargs = {
                "truncate": args.truncate,
                "prior_variant_probability": args.prior_variant_probability,
            }
            if args.pseudocount is not None:
                pon_kwargs["pseudocount"] = args.pseudocount
            payload = annotate_snv_vcf_with_normals(
                alignment_file,
                Path(args.vcf),
                normal_alignments=normal_alignments,
                output_path=args.output,
                min_baseq=args.min_baseq,
                min_mapq=args.min_mapq,
                **pon_kwargs,
            )

        if args.output is None:
            print(payload, end="")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2
