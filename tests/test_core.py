from dataclasses import dataclass

from skua.core import verify_snv_variant
from skua.variants import Variant


@dataclass
class FakeRead:
    mapping_quality: int
    is_reverse: bool
    query_sequence: str
    query_qualities: list[int]
    aligned_pairs: list[tuple[int | None, int | None]]


class FakeAlignmentFile:
    def __init__(self, reads: list[FakeRead]) -> None:
        self._reads = reads
        self.fetch_calls: list[tuple[str, int, int]] = []

    def fetch(self, contig: str, start: int, stop: int):
        self.fetch_calls.append((contig, start, stop))
        return iter(self._reads)


def build_linear_pairs(read_len: int, ref_start: int) -> list[tuple[int, int]]:
    return [(qpos, ref_start + qpos) for qpos in range(read_len)]


def test_verify_snv_variant_collects_evidence_for_single_variant() -> None:
    reads = [
        FakeRead(
            mapping_quality=60,
            is_reverse=False,
            query_sequence="AAAAATAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
        FakeRead(
            mapping_quality=60,
            is_reverse=True,
            query_sequence="AAAAAAAAAA",
            query_qualities=[35] * 10,
            aligned_pairs=build_linear_pairs(10, 100),
        ),
    ]
    alignment_file = FakeAlignmentFile(reads)
    variant = Variant(contig="chr1", ref_pos0=105, ref="A", alt="T")

    counts = verify_snv_variant(
        alignment_file,
        variant,
        min_baseq=20,
        min_mapq=20,
    )

    assert alignment_file.fetch_calls == [("chr1", 105, 106)]
    assert counts.alt_forward == 1
    assert counts.alt_reverse == 0
    assert counts.non_alt_forward == 0
    assert counts.non_alt_reverse == 1
    assert counts.usable == 2
    assert counts.unusable == 0
