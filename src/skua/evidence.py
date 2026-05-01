"""Read-level evidence classification primitives for variant verification."""

from dataclasses import dataclass
from enum import Enum
from typing import Any


class AlleleSupport(str, Enum):
    """High-level read support classification at a variant locus."""

    ALT = "alt"
    NON_ALT = "non_alt"
    UNUSABLE = "unusable"


class UnusableReason(str, Enum):
    """Reason for excluding a read from evidence counting."""

    LOW_MAPQ = "low_mapq"
    LOW_BASEQ = "low_baseq"
    NO_BASE_AT_SITE = "no_base_at_site"
    INVALID_BASE = "invalid_base"


@dataclass(frozen=True)
class ReadAlleleCall:
    """Result of classifying one read at one SNV locus."""

    support: AlleleSupport
    is_reverse: bool
    reason: UnusableReason | None = None
    observed_base: str | None = None
    base_quality: int | None = None


def classify_snv_read(
    read: Any,
    *,
    ref_pos0: int,
    ref_base: str,
    alt_base: str,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> ReadAlleleCall:
    """Classify one read as ALT, NON_ALT, or UNUSABLE for a single SNV."""
    if read.mapping_quality < min_mapq:
        return ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=read.is_reverse,
            reason=UnusableReason.LOW_MAPQ,
        )

    query_pos: int | None = None
    for qpos, rpos in read.aligned_pairs:
        if rpos == ref_pos0:
            query_pos = qpos
            break

    if query_pos is None:
        return ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=read.is_reverse,
            reason=UnusableReason.NO_BASE_AT_SITE,
        )

    observed_base = read.query_sequence[query_pos]
    base_quality = read.query_qualities[query_pos]

    if observed_base not in {"A", "C", "G", "T"}:
        return ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=read.is_reverse,
            reason=UnusableReason.INVALID_BASE,
            observed_base=observed_base,
            base_quality=base_quality,
        )

    if base_quality < min_baseq:
        return ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=read.is_reverse,
            reason=UnusableReason.LOW_BASEQ,
            observed_base=observed_base,
            base_quality=base_quality,
        )

    support = AlleleSupport.ALT if observed_base == alt_base else AlleleSupport.NON_ALT
    return ReadAlleleCall(
        support=support,
        is_reverse=read.is_reverse,
        observed_base=observed_base,
        base_quality=base_quality,
    )
