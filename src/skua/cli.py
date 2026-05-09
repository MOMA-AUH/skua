"""Command-line interface for skua."""

import argparse
from pathlib import Path

import pysam

from .core import (
    verify_snv_vcf_to_annotated_vcf,
    verify_snv_vcf_to_annotated_vcf_with_normals,
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
    subparsers = parser.add_subparsers(dest="command", required=True)

    verify_parser = subparsers.add_parser(
        "verify",
        help="Verify SNV evidence and emit annotated VCF output",
        formatter_class=OptionalDefaultsHelpFormatter,
    )
    verify_parser.add_argument("--vcf", required=True, help="Input VCF path (required)")
    verify_parser.add_argument("--alignment", required=True, help="Input BAM/CRAM path (required)")
    verify_parser.add_argument("--reference", help="Reference FASTA path (required for CRAM)")
    verify_parser.add_argument(
        "--output",
        help="Optional output VCF path (.vcf or .vcf.gz)",
    )
    verify_parser.add_argument(
        "--normal-list",
        help="Path to file listing normal sample BAM/CRAM paths, one per line",
    )
    verify_parser.add_argument("--min-baseq", type=int, default=20, help="Minimum base quality")
    verify_parser.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality")
    verify_parser.add_argument(
        "--truncate",
        type=float,
        default=0.1,
        help="Truncation threshold for PON sample inclusion",
    )
    verify_parser.add_argument(
        "--pseudocount",
        type=float,
        default=None,
        help="Optional pseudocount for beta-binomial rate estimates",
    )
    verify_parser.add_argument(
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

    if args.command == "verify":
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

        # Open normal alignments if provided
        normal_alignments = []
        normal_paths: list[str] = []
        if args.normal_list is not None:
            for line in Path(args.normal_list).read_text(encoding="utf-8").splitlines():
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    normal_paths.append(stripped)

        if normal_paths:
            for normal_path in normal_paths:
                normal_path_obj = Path(normal_path)
                if normal_path_obj.suffix.lower() == ".cram" and args.reference is None:
                    parser.error("--reference is required for CRAM input")
                normal_alignments.append(
                    pysam.AlignmentFile(normal_path, "rb", **alignment_kwargs)
                )

        try:
            with pysam.AlignmentFile(args.alignment, "rb", **alignment_kwargs) as alignment_file:
                if normal_paths:
                    pon_kwargs = {
                        "truncate": args.truncate,
                        "prior_variant_probability": args.prior_variant_probability,
                    }
                    if args.pseudocount is not None:
                        pon_kwargs["pseudocount"] = args.pseudocount
                    payload = verify_snv_vcf_to_annotated_vcf_with_normals(
                        alignment_file,
                        Path(args.vcf),
                        normal_alignments=normal_alignments,
                        output_path=args.output,
                        min_baseq=args.min_baseq,
                        min_mapq=args.min_mapq,
                        **pon_kwargs,
                    )
                else:
                    payload = verify_snv_vcf_to_annotated_vcf(
                        alignment_file,
                        Path(args.vcf),
                        output_path=args.output,
                        min_baseq=args.min_baseq,
                        min_mapq=args.min_mapq,
                    )
        finally:
            for normal_alignment in normal_alignments:
                normal_alignment.close()

        if args.output is None:
            print(payload, end="")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2
