from types import SimpleNamespace

import pytest

from skua.evidence import AlleleSupport, UnusableReason, classify_variant_read
from tests.helpers import FakeRead, build_linear_pairs


def test_classify_alt_on_forward_strand() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="AAAAATAAAA",
        query_qualities=[35] * 10,
        aligned_pairs=build_linear_pairs(10, 100),
    )

    call = classify_variant_read(
        read,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.ALT
    assert call.is_reverse is False
    assert call.reason is None
    assert call.observed_base == "T"


def test_classify_substitution_reuses_supplied_position_map() -> None:
    read = SimpleNamespace(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="T",
        query_qualities=[35],
    )

    call = classify_variant_read(
        read,
        ref_pos0=100,
        ref_base="A",
        alt_base="T",
        ref_to_query={100: 0},
    )

    assert call.support == AlleleSupport.ALT


@pytest.mark.parametrize(
    ("ref_base", "alt_base"),
    [("A", "A"), ("A", "CT"), ("AT", "C"), ("AT", "GCA")],
)
def test_classify_rejects_invalid_simple_alleles(
    ref_base: str,
    alt_base: str,
) -> None:
    read = SimpleNamespace(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="ATG",
        query_qualities=[35] * 3,
        reference_start=100,
        cigartuples=((0, 1), (1, 1), (0, 1)),
    )

    with pytest.raises(ValueError, match="simple|different"):
        classify_variant_read(
            read,
            ref_pos0=100,
            ref_base=ref_base,
            alt_base=alt_base,
        )


def test_classify_alt_on_reverse_strand() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=True,
        query_sequence="AAAAATAAAA",
        query_qualities=[35] * 10,
        aligned_pairs=build_linear_pairs(10, 100),
    )

    call = classify_variant_read(
        read,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.ALT
    assert call.is_reverse is True
    assert call.reason is None
    assert call.observed_base == "T"


def test_classify_non_alt_read() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="AAAAAAAAAA",
        query_qualities=[35] * 10,
        aligned_pairs=build_linear_pairs(10, 100),
    )

    call = classify_variant_read(
        read,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.NON_ALT
    assert call.reason is None
    assert call.observed_base == "A"


def test_classify_unusable_for_deletion_at_locus() -> None:
    aligned_pairs: list[tuple[int | None, int | None]] = [
        (0, 100),
        (1, 101),
        (2, 102),
        (3, 103),
        (4, 104),
        (None, 105),
        (5, 106),
        (6, 107),
        (7, 108),
        (8, 109),
    ]
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="AAAAAAAAA",
        query_qualities=[35] * 9,
        aligned_pairs=aligned_pairs,
    )

    call = classify_variant_read(
        read,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.UNUSABLE
    assert call.reason == UnusableReason.NO_BASE_AT_SITE
    assert call.observed_base is None


def test_classify_unusable_for_low_base_quality() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="AAAAATAAAA",
        query_qualities=[35, 35, 35, 35, 35, 10, 35, 35, 35, 35],
        aligned_pairs=build_linear_pairs(10, 100),
    )

    call = classify_variant_read(
        read,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.UNUSABLE
    assert call.reason == UnusableReason.LOW_BASEQ
    assert call.observed_base == "T"


def test_classify_unusable_for_missing_base_quality_sentinel() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="T",
        query_qualities=[255],
        aligned_pairs=build_linear_pairs(1, 100),
    )

    call = classify_variant_read(
        read,
        ref_pos0=100,
        ref_base="A",
        alt_base="T",
        min_baseq=0,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.UNUSABLE
    assert call.reason == UnusableReason.MISSING_BASEQ


def test_classify_unusable_for_low_mapping_quality() -> None:
    read = FakeRead(
        mapping_quality=5,
        is_reverse=False,
        query_sequence="AAAAATAAAA",
        query_qualities=[35] * 10,
        aligned_pairs=build_linear_pairs(10, 100),
    )

    call = classify_variant_read(
        read,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.UNUSABLE
    assert call.reason == UnusableReason.LOW_MAPQ
    assert call.observed_base is None


def test_classify_unusable_for_invalid_base() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="AAAAANAAAA",
        query_qualities=[35] * 10,
        aligned_pairs=build_linear_pairs(10, 100),
    )

    call = classify_variant_read(
        read,
        ref_pos0=105,
        ref_base="A",
        alt_base="T",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.UNUSABLE
    assert call.reason == UnusableReason.INVALID_BASE
    assert call.observed_base == "N"


def test_classify_simple_insertion_as_alt() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="ATAAAAAAAA",
        query_qualities=[35] * 10,
        aligned_pairs=[
            (0, 100),
            (1, None),
            (2, 101),
            (3, 102),
            (4, 103),
            (5, 104),
            (6, 105),
            (7, 106),
            (8, 107),
            (9, 108),
        ],
        reference_start=100,
        cigartuples=((0, 1), (1, 1), (0, 8)),
    )

    call = classify_variant_read(
        read,
        ref_pos0=100,
        ref_base="A",
        alt_base="AT",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.ALT
    assert call.reason is None
    assert call.observed_base == "T"


def test_classify_insertion_reports_minimum_anchor_and_inserted_quality() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="ATGAAAAAAA",
        query_qualities=[35, 30, 21, 35, 35, 35, 35, 35, 35, 35],
        aligned_pairs=[
            (0, 100),
            (1, None),
            (2, None),
            (3, 101),
            (4, 102),
            (5, 103),
            (6, 104),
            (7, 105),
            (8, 106),
            (9, 107),
        ],
        reference_start=100,
        cigartuples=((0, 1), (1, 2), (0, 7)),
    )

    call = classify_variant_read(
        read,
        ref_pos0=100,
        ref_base="A",
        alt_base="ATG",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.ALT
    assert call.base_quality == 21


def test_classify_simple_deletion_as_alt() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="AAAAAAA",
        query_qualities=[35] * 7,
        aligned_pairs=[
            (0, 100),
            (None, 101),
            (1, 102),
            (2, 103),
            (3, 104),
            (4, 105),
            (5, 106),
            (6, 107),
        ],
        reference_start=100,
        cigartuples=((0, 1), (2, 1), (0, 6)),
    )

    call = classify_variant_read(
        read,
        ref_pos0=100,
        ref_base="AT",
        alt_base="A",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.ALT
    assert call.reason is None
    assert call.observed_base == "AT"


def test_classify_reference_spanning_read_as_non_alt_for_deletion() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="AAA",
        query_qualities=[35] * 3,
        aligned_pairs=build_linear_pairs(3, 100),
    )

    call = classify_variant_read(
        read,
        ref_pos0=100,
        ref_base="AT",
        alt_base="A",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.NON_ALT
    assert call.reason is None


def test_classify_larger_deletion_as_non_alt() -> None:
    read = FakeRead(
        mapping_quality=60,
        is_reverse=False,
        query_sequence="AAAAAAA",
        query_qualities=[35] * 7,
        aligned_pairs=[
            (0, 100),
            (None, 101),
            (None, 102),
            (1, 103),
            (2, 104),
            (3, 105),
            (4, 106),
            (5, 107),
            (6, 108),
        ],
        reference_start=100,
        cigartuples=((0, 1), (2, 2), (0, 6)),
    )

    call = classify_variant_read(
        read,
        ref_pos0=100,
        ref_base="AT",
        alt_base="A",
        min_baseq=20,
        min_mapq=20,
    )

    assert call.support == AlleleSupport.NON_ALT
    assert call.reason is None
