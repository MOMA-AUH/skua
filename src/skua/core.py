"""Core public API for skua."""

from typing import Any

from .evidence import AggregatedEvidence, collect_snv_evidence_from_alignment
from .variants import Variant


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
