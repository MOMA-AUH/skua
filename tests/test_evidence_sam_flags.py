import pytest

from skua.evidence import is_accepted_sam_flag


@pytest.mark.parametrize(
    "flag",
    [
        99,   # 0x63: paired, proper pair, mate reverse, first mate
        147,  # 0x93: paired, proper pair, read reverse, second mate
    ],
)
def test_accepts_primary_mapped_proper_pair_records(flag: int) -> None:
    assert is_accepted_sam_flag(flag)


@pytest.mark.parametrize(
    "flag",
    [
        0,             # No SAM flags: unpaired record.
        0x1,           # PAIRED, but not PROPER_PAIR.
        99 | 0x4,     # UNMAPPED.
        99 | 0x8,     # MATE_UNMAPPED.
        99 | 0x100,   # SECONDARY alignment.
        99 | 0x200,   # Failed vendor/platform quality checks.
        99 | 0x400,   # PCR or optical DUPLICATE.
        99 | 0x800,   # SUPPLEMENTARY alignment.
    ],
)
def test_rejects_records_outside_the_sam_flag_policy(flag: int) -> None:
    assert not is_accepted_sam_flag(flag)
