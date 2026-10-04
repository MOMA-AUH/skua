"""Compact run summaries and read compatibility for legacy provenance blobs."""

import base64
import binascii
import json
from pathlib import Path
from typing import Any

import pysam

from ._headers import remove_header_records


PROVENANCE_HEADER_KEY = "SKUA_PROVENANCE"
RUN_HEADER_KEY = "SKUA_RUN"


def evidence_provenance(
    policy_version: int, min_baseq: int, min_mapq: int, *,
    normal_read_groups: str = "assigned_to_sample",
) -> dict[str, Any]:
    """Keep serialization and artifact validation on one evidence-policy description."""
    return {
        "policy_version": policy_version, "min_baseq": min_baseq, "min_mapq": min_mapq,
        "mapq_255": "exclude", "normal_read_groups": normal_read_groups,
    }


def write_run_summary_header(header: Any, summary: dict[str, Any]) -> None:
    """Write effective settings as readable VCF fields, replacing old run metadata."""
    evidence = summary["evidence"]
    items = [
        ("SchemaVersion", summary["schema_version"]),
        ("SkuaVersion", summary["skua_version"]),
        ("Mode", summary["mode"]),
        ("EvidencePolicyVersion", evidence["policy_version"]),
        ("MinBaseQ", evidence["min_baseq"]),
        ("MinMapQ", evidence["min_mapq"]),
        ("MapQ255", evidence["mapq_255"]),
        ("CaseReadGroups", summary["case_read_groups"]),
        ("NormalReadGroups", evidence["normal_read_groups"]),
    ]
    model = summary["model"]
    if model is not None:
        items.extend([
            ("Truncate", model["truncate"]),
            ("Pseudocount", model["pseudocount"]),
            ("PriorPolicy", model["prior"]["policy"]),
            ("PriorFallback", model["prior"]["fallback"]),
        ])
        items.extend(
            ("".join(part.capitalize() for part in name.split("_")), value)
            for name, value in model["assessment_thresholds"].items()
        )
    remove_header_records(header, (PROVENANCE_HEADER_KEY, RUN_HEADER_KEY))
    header.add_meta(RUN_HEADER_KEY, items=[(key, str(value)) for key, value in items])


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
    """Read a legacy provenance blob; return None for files without one, including new outputs."""
    with pysam.VariantFile(str(path)) as source:
        return read_provenance_header(source.header)
