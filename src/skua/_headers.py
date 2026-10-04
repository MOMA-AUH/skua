"""Safe removal of htslib header records whose proxies can be invalidated."""

from typing import Any, Iterable


def unquote_header_value(value: Any) -> str:
    """Decode the optional surrounding quotes retained by pysam header proxies."""
    text = str(value)
    return text[1:-1] if len(text) >= 2 and text[0] == text[-1] == '"' else text


def remove_header_records(header: Any, keys: Iterable[str]) -> None:
    """Refresh record proxies after each deletion, including duplicate keys."""
    for key in set(keys):
        while True:
            record = next((r for r in header.records if r.key == key), None)
            if record is None:
                break
            record.remove()
