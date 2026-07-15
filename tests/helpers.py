from dataclasses import dataclass


@dataclass
class FakeAlignmentHeader:
    read_groups: list[dict[str, str]]

    def to_dict(self) -> dict[str, list[dict[str, str]]]:
        return {"RG": self.read_groups}


@dataclass
class FakeRead:
    mapping_quality: int
    is_reverse: bool
    query_sequence: str
    query_qualities: list[int]
    aligned_pairs: list[tuple[int | None, int | None]]


class FakeAlignmentFile:
    def __init__(self, reads: list[FakeRead], header: FakeAlignmentHeader | None = None) -> None:
        self._reads = reads
        self.header = header
        self.fetch_calls: list[tuple[str, int, int]] = []

    def fetch(self, contig: str, start: int, stop: int):
        self.fetch_calls.append((contig, start, stop))
        return iter(self._reads)


def build_linear_pairs(read_len: int, ref_start: int) -> list[tuple[int, int]]:
    return [(qpos, ref_start + qpos) for qpos in range(read_len)]
