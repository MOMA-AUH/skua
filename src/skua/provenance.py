"""Portable run provenance and streaming identities for local input files."""

import base64
import binascii
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import pysam

from ._headers import remove_header_records


PROVENANCE_HEADER_KEY = "SKUA_PROVENANCE"


def input_identity(path: str | bytes | Path | None) -> dict[str, Any]:
    """Hash local file bytes in bounded memory; explicitly label unavailable identity."""
    identity: dict[str, Any] = {
        "path": None if path is None else os.fsdecode(path),
        "identity_method": "unavailable", "sha256": None, "size_bytes": None,
    }
    if path is None or not Path(os.fsdecode(path)).is_file():
        return identity
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        before = os.fstat(source.fileno())
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
        after = os.fstat(source.fileno())
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Input changed while computing provenance: {os.fsdecode(path)}")
    identity.update(
        path=str(Path(os.fsdecode(path)).resolve()), identity_method="sha256_file_bytes",
        sha256=digest.hexdigest(), size_bytes=after.st_size,
    )
    return identity


def write_provenance_header(header: Any, provenance: dict[str, Any]) -> None:
    """Replace Skua's record; base64 keeps arbitrary paths safe in VCF headers."""
    payload = json.dumps(provenance, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    remove_header_records(header, (PROVENANCE_HEADER_KEY,))
    header.add_meta(PROVENANCE_HEADER_KEY, value=base64.b64encode(payload).decode("ascii"))


def read_provenance_header(header: Any) -> dict[str, Any] | None:
    """Read a versioned record, distinguishing old outputs from malformed metadata."""
    records = [record for record in header.records if record.key == PROVENANCE_HEADER_KEY]
    if not records:
        return None
    if len(records) != 1:
        raise ValueError("Expected exactly one SKUA_PROVENANCE record")
    try:
        payload = json.loads(base64.b64decode(records[0].value, validate=True))
    except (ValueError, TypeError, binascii.Error, UnicodeDecodeError) as exc:
        raise ValueError("Invalid SKUA_PROVENANCE payload") from exc
    if not isinstance(payload, dict) or type(payload.get("schema_version")) is not int or payload["schema_version"] != 1:
        raise ValueError("Unsupported Skua provenance schema version")
    return payload


def read_provenance(path: str | Path) -> dict[str, Any] | None:
    """Return run/build provenance from VCF, bgzip VCF, or BCF; None for legacy files."""
    with pysam.VariantFile(str(path)) as source:
        return read_provenance_header(source.header)
