"""Safe removal of htslib header records whose proxies can be invalidated."""

from typing import Any, Iterable


def remove_header_records(header: Any, keys: Iterable[str]) -> None:
    """Refresh record proxies after each deletion, including duplicate keys."""
    for key in set(keys):
        while True:
            record = next((r for r in header.records if r.key == key), None)
            if record is None:
                break
            record.remove()
