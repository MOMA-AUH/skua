"""Read-level evidence classification primitives for variant verification."""

from collections import Counter
from collections.abc import Iterable
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
    """Result of classifying one read at one variant locus."""

    support: AlleleSupport
    is_reverse: bool
    reason: UnusableReason | None = None
    observed_base: str | None = None
    base_quality: int | None = None


@dataclass(frozen=True)
class AggregatedEvidence:
    """Strand-aware summary of read-level allele calls at one locus."""

    alt_forward: int
    alt_reverse: int
    non_alt_forward: int
    non_alt_reverse: int
    usable: int
    unusable: int
    unusable_by_reason: dict[UnusableReason, int]


def _query_position_bases_and_qualities(read: Any, query_positions: list[int], *, min_baseq: int) -> tuple[str | None, UnusableReason | None, str | None]:
    """Return the observed read sequence across query positions or an unusable reason."""
    observed_bases: list[str] = []
    sequence = read.query_sequence
    qualities = read.query_qualities

    for query_pos in query_positions:
        if query_pos is None or query_pos < 0 or query_pos >= len(sequence):
            return None, UnusableReason.NO_BASE_AT_SITE, None

        observed_base = sequence[query_pos]
        if observed_base not in {"A", "C", "G", "T"}:
            return None, UnusableReason.INVALID_BASE, observed_base

        if qualities[query_pos] < min_baseq:
            return None, UnusableReason.LOW_BASEQ, observed_base

        observed_bases.append(observed_base)

    return "".join(observed_bases), None, None


def _ref_position_map(read: Any) -> dict[int, int | None]:
    """Map each reference position in the alignment to its query position or None."""
    ref_to_query: dict[int, int | None] = {}
    for query_pos, ref_pos in read.aligned_pairs:
        if ref_pos is None or ref_pos in ref_to_query:
            continue
        ref_to_query[ref_pos] = query_pos
    return ref_to_query


def _query_positions_for_ref_span(read: Any, *, ref_pos0: int, ref_span_len: int) -> list[int] | None:
    """Return query positions covering a contiguous reference span, if fully aligned."""
    ref_to_query = _ref_position_map(read)
    query_positions: list[int] = []
    for target_ref_pos in range(ref_pos0, ref_pos0 + ref_span_len):
        if target_ref_pos not in ref_to_query or ref_to_query[target_ref_pos] is None:
            return None
        query_positions.append(ref_to_query[target_ref_pos])
    return query_positions


def _query_positions_for_insertion(read: Any, *, ref_pos0: int) -> list[int] | None:
    """Return inserted query positions immediately after the anchor base, if any."""
    insertion_query_positions: list[int] = []
    seen_anchor = False

    for query_pos, ref_pos in read.aligned_pairs:
        if ref_pos == ref_pos0 and query_pos is not None:
            seen_anchor = True
            continue

        if not seen_anchor:
            continue

        if ref_pos is None and query_pos is not None:
            insertion_query_positions.append(query_pos)
            continue

        if ref_pos is not None:
            break

    return insertion_query_positions


def classify_variant_read(
    read: Any,
    *,
    ref_pos0: int,
    ref_base: str,
    alt_base: str,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> ReadAlleleCall:
    """Classify one read as ALT, NON_ALT, or UNUSABLE for a variant."""
    if read.mapping_quality < min_mapq:
        return ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=read.is_reverse,
            reason=UnusableReason.LOW_MAPQ,
        )

    ref_len = len(ref_base)
    alt_len = len(alt_base)

    # Simple substitutions, including MNVs.
    if ref_len == alt_len:
        query_positions = _query_positions_for_ref_span(read, ref_pos0=ref_pos0, ref_span_len=ref_len)
        if query_positions is None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=UnusableReason.NO_BASE_AT_SITE,
            )

        observed_sequence, unusable_reason, observed_base = _query_position_bases_and_qualities(
            read,
            query_positions,
            min_baseq=min_baseq,
        )
        if unusable_reason is not None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=unusable_reason,
                observed_base=observed_base,
            )

        support = AlleleSupport.ALT if observed_sequence == alt_base else AlleleSupport.NON_ALT
        return ReadAlleleCall(
            support=support,
            is_reverse=read.is_reverse,
            observed_base=observed_sequence,
            base_quality=min(read.query_qualities[qpos] for qpos in query_positions),
        )

    # Simple insertion.
    if ref_len == 1 and alt_len > 1:
        query_positions = _query_positions_for_ref_span(read, ref_pos0=ref_pos0, ref_span_len=1)
        if query_positions is None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=UnusableReason.NO_BASE_AT_SITE,
            )

        anchor_bases, unusable_reason, observed_base = _query_position_bases_and_qualities(
            read,
            query_positions,
            min_baseq=min_baseq,
        )
        if unusable_reason is not None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=unusable_reason,
                observed_base=observed_base,
            )

        inserted_query_positions = _query_positions_for_insertion(read, ref_pos0=ref_pos0)
        inserted_sequence, unusable_reason, observed_base = _query_position_bases_and_qualities(
            read,
            inserted_query_positions,
            min_baseq=min_baseq,
        )
        if unusable_reason is not None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=unusable_reason,
                observed_base=observed_base,
            )

        support = AlleleSupport.ALT if inserted_sequence == alt_base[1:] and anchor_bases == ref_base else AlleleSupport.NON_ALT
        return ReadAlleleCall(
            support=support,
            is_reverse=read.is_reverse,
            observed_base=inserted_sequence,
            base_quality=min(read.query_qualities[qpos] for qpos in query_positions),
        )

    # Simple deletion.
    if ref_len > 1 and alt_len == 1:
        ref_to_query = _ref_position_map(read)
        if ref_pos0 not in ref_to_query or ref_to_query[ref_pos0] is None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=UnusableReason.NO_BASE_AT_SITE,
            )

        query_positions = [ref_to_query[ref_pos0]]
        anchor_bases, unusable_reason, observed_base = _query_position_bases_and_qualities(
            read,
            query_positions,
            min_baseq=min_baseq,
        )
        if unusable_reason is not None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=unusable_reason,
                observed_base=observed_base,
            )

        deletion_query_positions: list[int | None] = []
        for target_ref_pos in range(ref_pos0 + 1, ref_pos0 + ref_len):
            if target_ref_pos not in ref_to_query:
                return ReadAlleleCall(
                    support=AlleleSupport.UNUSABLE,
                    is_reverse=read.is_reverse,
                    reason=UnusableReason.NO_BASE_AT_SITE,
                )
            deletion_query_positions.append(ref_to_query[target_ref_pos])

        support = AlleleSupport.ALT if all(query_pos is None for query_pos in deletion_query_positions) and anchor_bases == ref_base[:1] else AlleleSupport.NON_ALT
        return ReadAlleleCall(
            support=support,
            is_reverse=read.is_reverse,
            observed_base=ref_base,
            base_quality=min(read.query_qualities[qpos] for qpos in query_positions if qpos is not None),
        )

    raise ValueError("Only simple substitutions and simple indels are supported")


def aggregate_read_calls(calls: Iterable[ReadAlleleCall]) -> AggregatedEvidence:
    """Aggregate read-level calls into strand-aware evidence counts."""
    alt_forward = 0
    alt_reverse = 0
    non_alt_forward = 0
    non_alt_reverse = 0
    usable = 0
    unusable = 0
    unusable_by_reason: Counter[UnusableReason] = Counter()

    for call in calls:
        if call.support == AlleleSupport.ALT:
            usable += 1
            if call.is_reverse:
                alt_reverse += 1
            else:
                alt_forward += 1
            continue

        if call.support == AlleleSupport.NON_ALT:
            usable += 1
            if call.is_reverse:
                non_alt_reverse += 1
            else:
                non_alt_forward += 1
            continue

        unusable += 1
        if call.reason is not None:
            unusable_by_reason[call.reason] += 1

    return AggregatedEvidence(
        alt_forward=alt_forward,
        alt_reverse=alt_reverse,
        non_alt_forward=non_alt_forward,
        non_alt_reverse=non_alt_reverse,
        usable=usable,
        unusable=unusable,
        unusable_by_reason=dict(unusable_by_reason),
    )


def collect_evidence(
    reads: Iterable[Any],
    *,
    ref_pos0: int,
    ref_base: str,
    alt_base: str,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> AggregatedEvidence:
    """Collect strand-aware evidence for one variant from an iterable of reads."""
    calls = [
            classify_variant_read(
            read,
            ref_pos0=ref_pos0,
            ref_base=ref_base,
            alt_base=alt_base,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        )
        for read in reads
    ]
    return aggregate_read_calls(calls)


def collect_evidence_from_alignment(
    alignment_file: Any,
    *,
    contig: str,
    ref_pos0: int,
    ref_base: str,
    alt_base: str,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> AggregatedEvidence:
    """Fetch overlapping reads for one variant and collect strand-aware evidence."""
    reads = alignment_file.fetch(contig, ref_pos0, ref_pos0 + 1)
    return collect_evidence(
        reads,
        ref_pos0=ref_pos0,
        ref_base=ref_base,
        alt_base=alt_base,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
    )
