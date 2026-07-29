from dataclasses import dataclass, field


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
    flag: int = 0x3
    tags: dict[str, str] = field(default_factory=dict)

    def has_tag(self, name: str) -> bool:
        return name in self.tags

    def get_tag(self, name: str) -> str:
        return self.tags[name]


class FakeAlignmentFile:
    def __init__(
        self,
        reads: list[FakeRead],
        header: FakeAlignmentHeader | None = None,
        *,
        references: tuple[str, ...] | None = None,
        indexed: bool = True,
    ) -> None:
        self._reads = reads
        self.header = header
        self.references = references
        self._indexed = indexed
        self.fetch_calls: list[tuple[str, int, int]] = []

    def fetch(self, contig: str, start: int, stop: int):
        self.fetch_calls.append((contig, start, stop))
        return iter(self._reads)

    def has_index(self) -> bool:
        return self._indexed


def build_linear_pairs(read_len: int, ref_start: int) -> list[tuple[int, int]]:
    return [(qpos, ref_start + qpos) for qpos in range(read_len)]
