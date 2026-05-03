"""Core public API for skua."""

from pathlib import Path
from typing import Any
from typing import Iterable
from typing import Iterator

from .evidence import AggregatedEvidence, collect_snv_evidence_from_alignment
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
            }
        )
    return rows
