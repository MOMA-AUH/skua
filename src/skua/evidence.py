"""Read-level evidence classification primitives for variant verification."""

from bisect import bisect_left
from collections import Counter
from collections.abc import Iterable
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum, IntFlag
from typing import Any

from .variants import Variant, _normalize_simple_alleles


class AlleleSupport(str, Enum):
    """High-level read support classification at a variant locus."""

    ALT = "alt"
    NON_ALT = "non_alt"
    UNUSABLE = "unusable"


class UnusableReason(str, Enum):
    """Reason for excluding a read from evidence counting."""

    LOW_MAPQ = "low_mapq"
    LOW_BASEQ = "low_baseq"
    MISSING_BASEQ = "missing_baseq"
    NO_BASE_AT_SITE = "no_base_at_site"
    INVALID_BASE = "invalid_base"
    CONFLICTING_MATES = "conflicting_mates"


class SamFlag(IntFlag):
    """SAM flag bits used by the alignment-record acceptance policy."""

    PAIRED = 0x1
    PROPER_PAIR = 0x2
    UNMAPPED = 0x4
    MATE_UNMAPPED = 0x8
    FIRST_IN_PAIR = 0x40
    SECOND_IN_PAIR = 0x80
    SECONDARY = 0x100
    QC_FAIL = 0x200
    DUPLICATE = 0x400
    SUPPLEMENTARY = 0x800


REQUIRED_SAM_FLAGS = SamFlag.PAIRED | SamFlag.PROPER_PAIR
REJECTED_SAM_FLAGS = (
    SamFlag.UNMAPPED
    | SamFlag.MATE_UNMAPPED
    | SamFlag.SECONDARY
    | SamFlag.QC_FAIL
    | SamFlag.DUPLICATE
    | SamFlag.SUPPLEMENTARY
)


def is_accepted_sam_flag(flag: int) -> bool:
    """Return whether a SAM flag identifies a primary, usable proper-pair record."""
    sam_flag = SamFlag(flag)
    return (
        sam_flag & REQUIRED_SAM_FLAGS == REQUIRED_SAM_FLAGS
        and not sam_flag & REJECTED_SAM_FLAGS
    )


def _read_group_id(read: Any) -> str | None:
    """Return the alignment record's read-group ID when present."""
    has_tag = getattr(read, "has_tag", None)
    if has_tag is None or not has_tag("RG"):
        return None
    return str(read.get_tag("RG"))


def _group_reads_by_fragment(reads: Iterable[Any]) -> Iterable[list[Any]]:
    """Group overlapping alignment records by read group and query name."""
    reads_by_fragment: dict[tuple[str | None, str], list[Any]] = {}
    unnamed_reads: list[Any] = []
    for read in reads:
        query_name = getattr(read, "query_name", None)
        if query_name is None:
            unnamed_reads.append(read)
            continue
        fragment_key = (_read_group_id(read), query_name)
        reads_by_fragment.setdefault(fragment_key, []).append(read)

    yield from reads_by_fragment.values()
    for read in unnamed_reads:
        yield [read]


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


_CIGAR_ALIGNED_OPS = frozenset({0, 7, 8})  # M, =, X
_CIGAR_INSERTION = 1
_CIGAR_DELETION = 2
_CIGAR_REFERENCE_SKIP = 3
_CIGAR_SOFT_CLIP = 4
_CIGAR_HARD_CLIP = 5
_CIGAR_PADDING = 6


@dataclass(frozen=True)
class _AlignmentPositions:
    """Operation-aware alignment positions reused across variant calls."""

    pairs: tuple[tuple[int | None, int | None, int | None], ...]
    ref_to_query: dict[int, int | None]
    ref_to_pair_index: dict[int, int]


def _preferred_mate(
    read_calls: list[tuple[Any, ReadAlleleCall]],
) -> tuple[Any, ReadAlleleCall]:
    """Choose the first mate as the stable fragment representative when present."""
    return min(
        read_calls,
        key=lambda read_call: not (
            SamFlag(getattr(read_call[0], "flag", 0)) & SamFlag.FIRST_IN_PAIR
        ),
    )


def _resolve_fragment_call(
    read_calls: list[tuple[Any, ReadAlleleCall]],
) -> ReadAlleleCall:
    """Resolve all overlapping mate calls into one fragment-level observation."""
    usable_read_calls = [
        read_call
        for read_call in read_calls
        if read_call[1].support != AlleleSupport.UNUSABLE
    ]

    if not usable_read_calls:
        return _preferred_mate(read_calls)[1]

    if len(usable_read_calls) == 1:
        return usable_read_calls[0][1]

    supports = {call.support for _read, call in usable_read_calls}
    if len(supports) > 1:
        _read, representative = _preferred_mate(usable_read_calls)
        return ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=representative.is_reverse,
            reason=UnusableReason.CONFLICTING_MATES,
        )

    return _preferred_mate(usable_read_calls)[1]


def _query_position_bases_and_qualities(
    read: Any,
    query_positions: Sequence[int],
    *,
    min_baseq: int,
) -> tuple[str | None, UnusableReason | None, str | None]:
    """Return the observed read sequence across query positions or an unusable reason."""
    observed_bases: list[str] = []
    sequence = read.query_sequence
    qualities = read.query_qualities

    if sequence is None:
        return None, UnusableReason.NO_BASE_AT_SITE, None
    if qualities is None:
        return None, UnusableReason.MISSING_BASEQ, None

    for query_pos in query_positions:
        if query_pos < 0 or query_pos >= len(sequence):
            return None, UnusableReason.NO_BASE_AT_SITE, None

        if query_pos >= len(qualities):
            return None, UnusableReason.MISSING_BASEQ, None

        observed_base = sequence[query_pos]
        if observed_base not in {"A", "C", "G", "T"}:
            return None, UnusableReason.INVALID_BASE, observed_base

        base_quality = qualities[query_pos]
        if base_quality == 255:
            return None, UnusableReason.MISSING_BASEQ, observed_base
        if base_quality < min_baseq:
            return None, UnusableReason.LOW_BASEQ, observed_base

        observed_bases.append(observed_base)

    return "".join(observed_bases), None, None


def _alignment_positions(read: Any) -> _AlignmentPositions:
    """Expand a read's CIGAR while retaining the operation behind each gap."""
    cigar = getattr(read, "cigartuples", None)
    reference_start = getattr(read, "reference_start", None)
    pairs: list[tuple[int | None, int | None, int | None]] = []

    if cigar is not None and reference_start is not None:
        query_pos = 0
        ref_pos = int(reference_start)
        for cigar_op, length in cigar:
            if length <= 0:
                pairs.append((None, None, cigar_op))
                break

            if cigar_op in _CIGAR_ALIGNED_OPS:
                pairs.extend(
                    (query_pos + offset, ref_pos + offset, cigar_op)
                    for offset in range(length)
                )
                query_pos += length
                ref_pos += length
                continue

            if cigar_op in {_CIGAR_INSERTION, _CIGAR_SOFT_CLIP}:
                pairs.extend(
                    (query_pos + offset, None, cigar_op)
                    for offset in range(length)
                )
                query_pos += length
                continue

            if cigar_op in {_CIGAR_DELETION, _CIGAR_REFERENCE_SKIP}:
                pairs.extend(
                    (None, ref_pos + offset, cigar_op)
                    for offset in range(length)
                )
                ref_pos += length
                continue

            if cigar_op in {_CIGAR_HARD_CLIP, _CIGAR_PADDING}:
                pairs.append((None, None, cigar_op))
                continue

            # Unsupported or malformed operations make the rest of the CIGAR
            # unsafe to interpret at a variant boundary.
            pairs.append((None, None, cigar_op))
            break
    else:
        # Two-column aligned pairs are sufficient for aligned reference bases,
        # but their gaps cannot safely be interpreted as I/D rather than S/N.
        pairs.extend(
            (
                query_pos,
                ref_pos,
                0 if query_pos is not None and ref_pos is not None else None,
            )
            for query_pos, ref_pos in read.aligned_pairs
        )

    ref_to_query: dict[int, int | None] = {}
    ref_to_pair_index: dict[int, int] = {}
    for pair_index, (query_pos, ref_pos, _cigar_op) in enumerate(pairs):
        if ref_pos is None or ref_pos in ref_to_query:
            continue
        ref_to_query[ref_pos] = query_pos
        ref_to_pair_index[ref_pos] = pair_index

    return _AlignmentPositions(
        pairs=tuple(pairs),
        ref_to_query=ref_to_query,
        ref_to_pair_index=ref_to_pair_index,
    )


def _ref_position_map(read: Any) -> dict[int, int | None]:
    """Map each reference position in the alignment to its query position or None."""
    return _alignment_positions(read).ref_to_query


def _query_positions_for_ref_span(
    read: Any,
    *,
    ref_pos0: int,
    ref_span_len: int,
    ref_to_query: dict[int, int | None] | None = None,
) -> list[int] | None:
    """Return consecutive query positions spanning every requested reference base."""
    if ref_to_query is None:
        ref_to_query = _ref_position_map(read)
    query_positions: list[int] = []
    for target_ref_pos in range(ref_pos0, ref_pos0 + ref_span_len):
        query_pos = ref_to_query.get(target_ref_pos)
        if query_pos is None:
            return None
        # An internal insertion is complex evidence, not a simple MNV allele.
        if query_positions and query_pos != query_positions[-1] + 1:
            return None
        query_positions.append(query_pos)
    return query_positions


def _query_positions_for_contiguous_ref_span(
    alignment_positions: _AlignmentPositions,
    *,
    ref_pos0: int,
    ref_span_len: int,
) -> list[int] | None:
    """Return query positions for an uninterrupted aligned CIGAR span."""
    first_pair_index = alignment_positions.ref_to_pair_index.get(ref_pos0)
    if first_pair_index is None:
        return None

    query_positions: list[int] = []
    for offset in range(ref_span_len):
        pair_index = first_pair_index + offset
        if pair_index >= len(alignment_positions.pairs):
            return None
        query_pos, ref_pos, cigar_op = alignment_positions.pairs[pair_index]
        if (
            query_pos is None
            or ref_pos != ref_pos0 + offset
            or cigar_op not in _CIGAR_ALIGNED_OPS
        ):
            return None
        query_positions.append(query_pos)
    return query_positions


def _query_positions_for_insertion(
    alignment_positions: _AlignmentPositions,
    *,
    ref_pos0: int,
) -> list[int] | None:
    """Return a flanked CIGAR insertion after the anchor, or None if unobservable."""
    anchor_pair_index = alignment_positions.ref_to_pair_index.get(ref_pos0)
    if anchor_pair_index is None:
        return None

    anchor_query_pos, anchor_ref_pos, anchor_op = alignment_positions.pairs[
        anchor_pair_index
    ]
    if (
        anchor_query_pos is None
        or anchor_ref_pos != ref_pos0
        or anchor_op not in _CIGAR_ALIGNED_OPS
    ):
        return None

    insertion_query_positions: list[int] = []
    pair_index = anchor_pair_index + 1
    while pair_index < len(alignment_positions.pairs):
        query_pos, ref_pos, cigar_op = alignment_positions.pairs[pair_index]
        if cigar_op != _CIGAR_INSERTION:
            break
        if query_pos is None or ref_pos is not None:
            return None
        insertion_query_positions.append(query_pos)
        pair_index += 1

    if pair_index >= len(alignment_positions.pairs):
        return None

    flank_query_pos, flank_ref_pos, flank_op = alignment_positions.pairs[pair_index]
    if (
        flank_query_pos is None
        or flank_ref_pos != ref_pos0 + 1
        or flank_op not in _CIGAR_ALIGNED_OPS
    ):
        return None
    return insertion_query_positions


def _deleted_reference_positions(
    alignment_positions: _AlignmentPositions,
    *,
    ref_pos0: int,
) -> list[int] | None:
    """Return a flanked CIGAR deletion after the anchor, or None if absent/unsafe."""
    anchor_pair_index = alignment_positions.ref_to_pair_index.get(ref_pos0)
    if anchor_pair_index is None:
        return None

    anchor_query_pos, anchor_ref_pos, anchor_op = alignment_positions.pairs[
        anchor_pair_index
    ]
    if (
        anchor_query_pos is None
        or anchor_ref_pos != ref_pos0
        or anchor_op not in _CIGAR_ALIGNED_OPS
    ):
        return None

    deleted_ref_positions: list[int] = []
    pair_index = anchor_pair_index + 1
    while pair_index < len(alignment_positions.pairs):
        query_pos, ref_pos, cigar_op = alignment_positions.pairs[pair_index]
        if cigar_op != _CIGAR_DELETION:
            break
        if query_pos is not None or ref_pos is None:
            return None
        deleted_ref_positions.append(ref_pos)
        pair_index += 1

    if not deleted_ref_positions or pair_index >= len(alignment_positions.pairs):
        return None

    flank_query_pos, flank_ref_pos, flank_op = alignment_positions.pairs[pair_index]
    if (
        flank_query_pos is None
        or flank_ref_pos != deleted_ref_positions[-1] + 1
        or flank_op not in _CIGAR_ALIGNED_OPS
    ):
        return None
    return deleted_ref_positions


def classify_variant_read(
    read: Any,
    *,
    ref_pos0: int,
    ref_base: str,
    alt_base: str,
    min_baseq: int = 20,
    min_mapq: int = 20,
    ref_to_query: dict[int, int | None] | None = None,
    alignment_positions: _AlignmentPositions | None = None,
) -> ReadAlleleCall:
    """Classify one read as ALT, NON_ALT, or UNUSABLE for a variant.

    The alignment caches can be supplied by a batch caller to reuse a read's
    parsed positions across every variant that the read overlaps.
    """
    ref_base, alt_base = _normalize_simple_alleles(ref_base, alt_base)

    if read.mapping_quality < min_mapq:
        return ReadAlleleCall(
            support=AlleleSupport.UNUSABLE,
            is_reverse=read.is_reverse,
            reason=UnusableReason.LOW_MAPQ,
        )

    ref_len = len(ref_base)
    alt_len = len(alt_base)

    if alignment_positions is not None:
        ref_to_query = alignment_positions.ref_to_query
    elif ref_to_query is None or ref_len != alt_len:
        alignment_positions = _alignment_positions(read)
        ref_to_query = alignment_positions.ref_to_query

    # Simple substitutions, including MNVs.
    if ref_len == alt_len:
        query_positions = _query_positions_for_ref_span(
            read,
            ref_pos0=ref_pos0,
            ref_span_len=ref_len,
            ref_to_query=ref_to_query,
        )
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
        assert alignment_positions is not None
        query_positions = _query_positions_for_ref_span(
            read,
            ref_pos0=ref_pos0,
            ref_span_len=1,
            ref_to_query=ref_to_query,
        )
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

        inserted_query_positions = _query_positions_for_insertion(
            alignment_positions,
            ref_pos0=ref_pos0,
        )
        if inserted_query_positions is None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=UnusableReason.NO_BASE_AT_SITE,
            )

        allele_query_positions = query_positions + inserted_query_positions
        observed_allele, unusable_reason, observed_base = (
            _query_position_bases_and_qualities(
                read,
                allele_query_positions,
                min_baseq=min_baseq,
            )
        )
        if unusable_reason is not None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=unusable_reason,
                observed_base=observed_base,
            )

        inserted_sequence = observed_allele[1:]
        support = (
            AlleleSupport.ALT
            if inserted_sequence == alt_base[1:] and anchor_bases == ref_base
            else AlleleSupport.NON_ALT
        )
        return ReadAlleleCall(
            support=support,
            is_reverse=read.is_reverse,
            observed_base=inserted_sequence,
            base_quality=min(
                read.query_qualities[qpos] for qpos in allele_query_positions
            ),
        )

    # Simple deletion.
    if ref_len > 1 and alt_len == 1:
        assert alignment_positions is not None
        if ref_pos0 not in ref_to_query or ref_to_query[ref_pos0] is None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=UnusableReason.NO_BASE_AT_SITE,
            )

        query_positions = [ref_to_query[ref_pos0]]
        anchor_bases, unusable_reason, observed_base = (
            _query_position_bases_and_qualities(
                read,
                query_positions,
                min_baseq=min_baseq,
            )
        )
        if unusable_reason is not None:
            return ReadAlleleCall(
                support=AlleleSupport.UNUSABLE,
                is_reverse=read.is_reverse,
                reason=unusable_reason,
                observed_base=observed_base,
            )

        reference_query_positions = _query_positions_for_contiguous_ref_span(
            alignment_positions,
            ref_pos0=ref_pos0,
            ref_span_len=ref_len + 1,
        )
        if reference_query_positions is not None:
            support = AlleleSupport.NON_ALT
        else:
            deleted_ref_positions = _deleted_reference_positions(
                alignment_positions,
                ref_pos0=ref_pos0,
            )
            if deleted_ref_positions is None:
                return ReadAlleleCall(
                    support=AlleleSupport.UNUSABLE,
                    is_reverse=read.is_reverse,
                    reason=UnusableReason.NO_BASE_AT_SITE,
                )

            expected_deleted_positions = list(
                range(ref_pos0 + 1, ref_pos0 + ref_len)
            )
            support = (
                AlleleSupport.ALT
                if deleted_ref_positions == expected_deleted_positions
                and anchor_bases == ref_base[:1]
                else AlleleSupport.NON_ALT
            )

        return ReadAlleleCall(
            support=support,
            is_reverse=read.is_reverse,
            observed_base=ref_base,
            base_quality=min(read.query_qualities[qpos] for qpos in query_positions),
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
    ref_base, alt_base = _normalize_simple_alleles(ref_base, alt_base)
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
    allowed_read_group_ids: frozenset[str] | None = None,
) -> AggregatedEvidence:
    """Fetch overlapping reads for one variant and collect strand-aware evidence.

    When ``allowed_read_group_ids`` is supplied, only reads assigned to one of
    those read groups contribute evidence. This is used to isolate one sample
    from a multi-sample alignment.
    """
    ref_base, alt_base = _normalize_simple_alleles(ref_base, alt_base)
    reads = (
        read
        for read in alignment_file.fetch(contig, ref_pos0, ref_pos0 + 1)
        if is_accepted_sam_flag(read.flag)
        and (
            allowed_read_group_ids is None
            or _read_group_id(read) in allowed_read_group_ids
        )
    )
    fragment_calls: list[ReadAlleleCall] = []
    for fragment_reads in _group_reads_by_fragment(reads):
        read_calls = [
            (
                read,
                classify_variant_read(
                    read,
                    ref_pos0=ref_pos0,
                    ref_base=ref_base,
                    alt_base=alt_base,
                    min_baseq=min_baseq,
                    min_mapq=min_mapq,
                ),
            )
            for read in fragment_reads
        ]
        fragment_calls.append(_resolve_fragment_call(read_calls))

    return aggregate_read_calls(fragment_calls)


def _read_reference_span(read: Any) -> tuple[int, int] | None:
    """Return the half-open reference span covered by an alignment record."""
    reference_start = getattr(read, "reference_start", None)
    reference_end = getattr(read, "reference_end", None)
    if reference_start is not None and reference_end is not None:
        return int(reference_start), int(reference_end)

    reference_positions = [
        ref_pos
        for _query_pos, ref_pos in read.aligned_pairs
        if ref_pos is not None
    ]
    if not reference_positions:
        return None
    return min(reference_positions), max(reference_positions) + 1


def collect_evidence_from_alignment_batch(
    alignment_file: Any,
    variants: Sequence[Variant],
    *,
    min_baseq: int = 20,
    min_mapq: int = 20,
    allowed_read_group_ids: frozenset[str] | None = None,
) -> tuple[AggregatedEvidence, ...]:
    """Collect evidence for same-contig variants with one alignment fetch.

    Results correspond positionally to ``variants``. Each variant retains the
    same read filtering and fragment-resolution semantics as
    :func:`collect_evidence_from_alignment`.
    """
    if not variants:
        return ()

    contig = variants[0].contig
    if any(variant.contig != contig for variant in variants):
        raise ValueError("A variant batch must contain exactly one contig")

    sorted_variants = sorted(
        enumerate(variants),
        key=lambda indexed_variant: indexed_variant[1].ref_pos0,
    )
    sorted_positions = [variant.ref_pos0 for _index, variant in sorted_variants]
    fragments_by_variant: list[dict[Any, list[tuple[Any, ReadAlleleCall]]]] = [
        {} for _variant in variants
    ]
    unnamed_read_index = 0

    for read in alignment_file.fetch(
        contig,
        sorted_positions[0],
        sorted_positions[-1] + 1,
    ):
        if not is_accepted_sam_flag(read.flag):
            continue
        read_group_id = _read_group_id(read)
        if (
            allowed_read_group_ids is not None
            and read_group_id not in allowed_read_group_ids
        ):
            continue

        reference_span = _read_reference_span(read)
        if reference_span is None:
            continue
        reference_start, reference_end = reference_span
        first_variant = bisect_left(sorted_positions, reference_start)
        after_last_variant = bisect_left(sorted_positions, reference_end)
        if first_variant == after_last_variant:
            continue

        query_name = getattr(read, "query_name", None)
        if query_name is None:
            fragment_key: Any = (None, unnamed_read_index)
            unnamed_read_index += 1
        else:
            fragment_key = (read_group_id, query_name)

        alignment_positions = (
            _alignment_positions(read) if read.mapping_quality >= min_mapq else None
        )
        ref_to_query = (
            alignment_positions.ref_to_query
            if alignment_positions is not None
            else None
        )

        for sorted_variant_index in range(first_variant, after_last_variant):
            original_index, variant = sorted_variants[sorted_variant_index]
            read_call = classify_variant_read(
                read,
                ref_pos0=variant.ref_pos0,
                ref_base=variant.ref,
                alt_base=variant.alt,
                min_baseq=min_baseq,
                min_mapq=min_mapq,
                ref_to_query=ref_to_query,
                alignment_positions=alignment_positions,
            )
            fragments_by_variant[original_index].setdefault(fragment_key, []).append(
                (read, read_call)
            )

    return tuple(
        aggregate_read_calls(
            _resolve_fragment_call(read_calls)
            for read_calls in fragment_calls.values()
        )
        for fragment_calls in fragments_by_variant
    )
