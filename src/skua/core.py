"""Core public API for skua."""

import json
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Iterable
from typing import Iterator

from .evidence import AggregatedEvidence, collect_snv_evidence_from_alignment
from .stats import compute_stats
from .variants import Variant, read_vcf_snv_file


def verify_snv_variant(
    alignment_file: Any,
    variant: Variant,
    *,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> AggregatedEvidence:
    """Collect strand-aware evidence for one SNV variant from one alignment."""
    return collect_snv_evidence_from_alignment(
        alignment_file,
        contig=variant.contig,
        ref_pos0=variant.ref_pos0,
        ref_base=variant.ref,
        alt_base=variant.alt,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
    )


def verify_snv_variant_with_normals(
    alignment_file: Any,
    variant: Variant,
    *,
    normal_alignments: list[Any] | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> dict[str, Any]:
    """Collect case and normal evidence for one SNV variant."""
    if normal_alignments is None:
        normal_alignments = []

    case_evidence = verify_snv_variant(
        alignment_file,
        variant,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
    )

    normal_evidences: list[AggregatedEvidence] = []
    normal_aggregate_evidence = AggregatedEvidence(
        alt_forward=0,
        alt_reverse=0,
        non_alt_forward=0,
        non_alt_reverse=0,
        usable=0,
        unusable=0,
        unusable_by_reason={},
    )

    for normal_alignment in normal_alignments:
        normal_evidence = verify_snv_variant(
            normal_alignment,
            variant,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        )
        normal_evidences.append(normal_evidence)

        normal_unusable_by_reason = dict(normal_aggregate_evidence.unusable_by_reason)
        for reason, count in normal_evidence.unusable_by_reason.items():
            normal_unusable_by_reason[reason] = normal_unusable_by_reason.get(reason, 0) + count

        normal_aggregate_evidence = AggregatedEvidence(
            alt_forward=normal_aggregate_evidence.alt_forward + normal_evidence.alt_forward,
            alt_reverse=normal_aggregate_evidence.alt_reverse + normal_evidence.alt_reverse,
            non_alt_forward=normal_aggregate_evidence.non_alt_forward + normal_evidence.non_alt_forward,
            non_alt_reverse=normal_aggregate_evidence.non_alt_reverse + normal_evidence.non_alt_reverse,
            usable=normal_aggregate_evidence.usable + normal_evidence.usable,
            unusable=normal_aggregate_evidence.unusable + normal_evidence.unusable,
            unusable_by_reason=normal_unusable_by_reason,
        )

    return {
        "case_evidence": case_evidence,
        "normal_evidences": normal_evidences,
        "normal_aggregate_evidence": normal_aggregate_evidence,
    }


def verify_snv_variants_from_vcf(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> Iterator[tuple[Variant, AggregatedEvidence]]:
    """Yield per-variant evidence for SNV records from a VCF file."""
    for variant in read_vcf_snv_file(vcf_path):
        yield (
            variant,
            verify_snv_variant(
                alignment_file,
                variant,
                min_baseq=min_baseq,
                min_mapq=min_mapq,
            ),
        )


def format_verification_results(
    results: Iterable[tuple[Variant, AggregatedEvidence]],
) -> list[dict[str, Any]]:
    """Convert verification results to JSON/tabular-ready row dictionaries."""
    rows: list[dict[str, Any]] = []
    for variant, evidence in results:
        rows.append(
            {
                "contig": variant.contig,
                "pos1": variant.ref_pos0 + 1,
                "ref": variant.ref,
                "alt": variant.alt,
                "case": {
                    "alt_forward": evidence.alt_forward,
                    "alt_reverse": evidence.alt_reverse,
                    "non_alt_forward": evidence.non_alt_forward,
                    "non_alt_reverse": evidence.non_alt_reverse,
                    "usable": evidence.usable,
                    "unusable": evidence.unusable,
                    "unusable_by_reason": {
                        reason.value: count
                        for reason, count in evidence.unusable_by_reason.items()
                    },
                },
            }
        )
    return rows


def render_verification_results_json(rows: Iterable[dict[str, Any]]) -> str:
    """Render formatted verification rows as JSON text."""
    return json.dumps(list(rows), indent=2)


def write_verification_results_json(
    rows: Iterable[dict[str, Any]],
    output_path: str | Path,
) -> None:
    """Write formatted verification rows to a JSON file."""
    Path(output_path).write_text(
        render_verification_results_json(rows),
        encoding="utf-8",
    )


def _build_verification_rows(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    min_baseq: int,
    min_mapq: int,
) -> list[dict[str, Any]]:
    """Build formatted verification rows from one alignment and one VCF."""
    return format_verification_results(
        verify_snv_variants_from_vcf(
            alignment_file,
            vcf_path,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        )
    )


def _render_and_optionally_write(
    rows: Iterable[dict[str, Any]],
    *,
    renderer: Callable[[Iterable[dict[str, Any]]], str],
    output_path: str | Path | None,
) -> str:
    """Render rows and optionally persist the payload to disk."""
    payload = renderer(rows)
    if output_path is not None:
        Path(output_path).write_text(payload, encoding="utf-8")
    return payload


def verify_snv_vcf_to_json(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    output_path: str | Path | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> str:
    """Run SNV verification from VCF and return JSON output, optionally writing to file."""
    rows = _build_verification_rows(
        alignment_file,
        vcf_path,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
    )
    return _render_and_optionally_write(
        rows,
        renderer=render_verification_results_json,
        output_path=output_path,
    )


def verify_snv_variants_from_vcf_with_normals(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    normal_alignments: list[Any] | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> Iterator[tuple[Variant, dict[str, Any]]]:
    """Yield per-variant case+normal evidence for SNV records from a VCF file."""
    if normal_alignments is None:
        normal_alignments = []

    for variant in read_vcf_snv_file(vcf_path):
        yield (
            variant,
            verify_snv_variant_with_normals(
                alignment_file,
                variant,
                normal_alignments=normal_alignments,
                min_baseq=min_baseq,
                min_mapq=min_mapq,
            ),
        )


def format_verification_results_with_normals(
    results: Iterable[tuple[Variant, dict[str, Any]]],
) -> list[dict[str, Any]]:
    """Convert PON verification results to JSON/tabular-ready row dictionaries."""
    rows: list[dict[str, Any]] = []
    for variant, pon_result in results:
        evidence = pon_result["case_evidence"]
        normal_aggregate_evidence = pon_result["normal_aggregate_evidence"]
        stats = compute_stats(
            evidence,
            normal_aggregate_evidence,
            per_sample_evidences=pon_result["normal_evidences"],
        )
        rows.append(
            {
                "contig": variant.contig,
                "pos1": variant.ref_pos0 + 1,
                "ref": variant.ref,
                "alt": variant.alt,
                "bayes_factor": stats.bayes_factor,
                "artifact_posterior": stats.artifact_posterior,
                "case": {
                    "alt_forward": evidence.alt_forward,
                    "alt_reverse": evidence.alt_reverse,
                    "non_alt_forward": evidence.non_alt_forward,
                    "non_alt_reverse": evidence.non_alt_reverse,
                    "usable": evidence.usable,
                    "unusable": evidence.unusable,
                    "unusable_by_reason": {
                        reason.value: count
                        for reason, count in evidence.unusable_by_reason.items()
                    },
                },
                "normal": {
                    "alt_forward": normal_aggregate_evidence.alt_forward,
                    "alt_reverse": normal_aggregate_evidence.alt_reverse,
                    "non_alt_forward": normal_aggregate_evidence.non_alt_forward,
                    "non_alt_reverse": normal_aggregate_evidence.non_alt_reverse,
                    "usable": normal_aggregate_evidence.usable,
                    "unusable": normal_aggregate_evidence.unusable,
                    "unusable_by_reason": {
                        reason.value: count
                        for reason, count in normal_aggregate_evidence.unusable_by_reason.items()
                    },
                },
            }
        )
    return rows


def _build_verification_rows_with_normals(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    normal_alignments: list[Any] | None,
    min_baseq: int,
    min_mapq: int,
) -> list[dict[str, Any]]:
    """Build formatted PON verification rows from case + normal alignments and one VCF."""
    return format_verification_results_with_normals(
        verify_snv_variants_from_vcf_with_normals(
            alignment_file,
            vcf_path,
            normal_alignments=normal_alignments,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        )
    )


def verify_snv_vcf_to_json_with_normals(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    normal_alignments: list[Any] | None = None,
    output_path: str | Path | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> str:
    """Run PON SNV verification from VCF and return JSON output, optionally writing to file."""
    if normal_alignments is None:
        normal_alignments = []

    rows = _build_verification_rows_with_normals(
        alignment_file,
        vcf_path,
        normal_alignments=normal_alignments,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
    )
    return _render_and_optionally_write(
        rows,
        renderer=render_verification_results_json,
        output_path=output_path,
    )


