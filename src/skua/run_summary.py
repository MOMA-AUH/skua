"""Compact, human-readable summaries of effective annotation settings."""

from typing import Any

from ._headers import remove_header_records


RUN_SUMMARY_SCHEMA_VERSION = 1
RUN_HEADER_KEY = "SKUA_RUN"


def evidence_summary(
    policy_version: int, min_baseq: int, min_mapq: int, *,
    normal_read_groups: str = "assigned_to_sample",
) -> dict[str, Any]:
    """Describe quality filters and normal read-selection policy."""
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
    remove_header_records(header, (RUN_HEADER_KEY,))
    header.add_meta(RUN_HEADER_KEY, items=[(key, str(value)) for key, value in items])
