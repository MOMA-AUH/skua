"""Versioned storage for allele-targeted panel-of-normals evidence."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator

import pysam
from pysam import bcftools

from ._version import __version__
from .evidence import AggregatedEvidence
from .variants import Variant


PON_SCHEMA_VERSION = 1
EVIDENCE_POLICY_VERSION = 1
PON_HEADER_KEY = "SKUA_PON"

PON_EVIDENCE_FORMAT_FIELDS: tuple[tuple[str, str], ...] = (
    ("SKUA_PON_AF", "PON sample ALT-supporting forward reads"),
    ("SKUA_PON_AR", "PON sample ALT-supporting reverse reads"),
    ("SKUA_PON_NF", "PON sample non-ALT forward reads"),
    ("SKUA_PON_NR", "PON sample non-ALT reverse reads"),
    ("SKUA_PON_U", "PON sample usable reads"),
    ("SKUA_PON_X", "PON sample unusable reads"),
)

_EVIDENCE_ATTRIBUTES_BY_FIELD = {
    "SKUA_PON_AF": "alt_forward",
    "SKUA_PON_AR": "alt_reverse",
    "SKUA_PON_NF": "non_alt_forward",
    "SKUA_PON_NR": "non_alt_reverse",
    "SKUA_PON_U": "usable",
    "SKUA_PON_X": "unusable",
}


@dataclass(frozen=True)
class PonArtifactMetadata:
    """Provenance required to interpret a precomputed PON artifact."""

    schema_version: int
    evidence_policy_version: int
    min_baseq: int
    min_mapq: int
    skua_version: str
    sample_names: tuple[str, ...]


def _unquote_header_value(value: Any) -> str:
    text = str(value)
    if len(text) >= 2 and text[0] == text[-1] == '"':
        return text[1:-1]
    return text


def _metadata_items(header: Any) -> dict[str, str]:
    records = [record for record in header.records if record.key == PON_HEADER_KEY]
    if len(records) != 1:
        raise ValueError(
            f"PON artifact must contain exactly one {PON_HEADER_KEY} metadata record"
        )
    return {
        key: _unquote_header_value(value)
        for key, value in records[0].items()
        if key != "IDX"
    }


def _parse_metadata(header: Any) -> PonArtifactMetadata:
    items = _metadata_items(header)
    required = {
        "SchemaVersion",
        "EvidencePolicyVersion",
        "MinBaseQ",
        "MinMapQ",
        "SkuaVersion",
    }
    missing = sorted(required - items.keys())
    if missing:
        raise ValueError("PON artifact metadata is missing: " + ", ".join(missing))

    try:
        schema_version = int(items["SchemaVersion"])
        evidence_policy_version = int(items["EvidencePolicyVersion"])
        min_baseq = int(items["MinBaseQ"])
        min_mapq = int(items["MinMapQ"])
    except ValueError as exc:
        raise ValueError("PON artifact contains invalid integer metadata") from exc

    if schema_version != PON_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported PON schema version {schema_version}; expected {PON_SCHEMA_VERSION}"
        )
    if evidence_policy_version != EVIDENCE_POLICY_VERSION:
        raise ValueError(
            "Unsupported PON evidence policy version "
            f"{evidence_policy_version}; expected {EVIDENCE_POLICY_VERSION}"
        )

    sample_names = tuple(header.samples)
    if not sample_names:
        raise ValueError("PON artifact must contain at least one normal sample")

    missing_fields = [
        field_id for field_id, _description in PON_EVIDENCE_FORMAT_FIELDS
        if field_id not in header.formats
    ]
    if missing_fields:
        raise ValueError("PON artifact is missing FORMAT fields: " + ", ".join(missing_fields))

    return PonArtifactMetadata(
        schema_version=schema_version,
        evidence_policy_version=evidence_policy_version,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
        skua_version=items["SkuaVersion"],
        sample_names=sample_names,
    )


def read_pon_metadata(path: str | Path) -> PonArtifactMetadata:
    """Read and validate PON artifact metadata without loading its evidence."""
    with pysam.VariantFile(str(path)) as pon_file:
        return _parse_metadata(pon_file.header)


def _add_pon_header_fields(
    header: Any,
    *,
    min_baseq: int,
    min_mapq: int,
) -> None:
    if any(record.key == PON_HEADER_KEY for record in header.records):
        raise ValueError("Target VCF already contains SKUA_PON metadata")
    header.add_meta(
        PON_HEADER_KEY,
        items=[
            ("SchemaVersion", str(PON_SCHEMA_VERSION)),
            ("EvidencePolicyVersion", str(EVIDENCE_POLICY_VERSION)),
            ("MinBaseQ", str(min_baseq)),
            ("MinMapQ", str(min_mapq)),
            ("SkuaVersion", __version__),
        ],
    )
    for field_id, description in PON_EVIDENCE_FORMAT_FIELDS:
        if field_id in header.formats:
            field = header.formats[field_id]
            if field.number != 1 or field.type != "Integer":
                raise ValueError(
                    f"Target VCF contains an incompatible {field_id} FORMAT definition"
                )
            continue
        definition = (
            f'##FORMAT=<ID={field_id},Number=1,Type=Integer,'
            f'Description="{description}">'
        )
        header.add_line(definition)


def _copy_target_record(record: Any, output_file: Any) -> Any:
    copied = output_file.new_record(
        contig=record.contig,
        start=record.start,
        stop=record.stop,
        id=record.id,
        alleles=record.alleles,
        qual=record.qual,
    )
    for filter_id in record.filter.keys():
        copied.filter.add(filter_id)
    for key, value in record.info.items():
        copied.info[key] = value
    return copied


def write_pon_artifact(
    target_vcf_path: str | Path,
    output_path: str | Path,
    *,
    sample_names: tuple[str, ...],
    evidence_records: Iterable[tuple[Variant, tuple[AggregatedEvidence, ...]]],
    min_baseq: int,
    min_mapq: int,
) -> None:
    """Write per-normal, per-allele evidence to an immutable BCF artifact."""
    if not sample_names:
        raise ValueError("PON artifact requires at least one normal sample")
    if len(set(sample_names)) != len(sample_names):
        raise ValueError("PON normal sample names must be unique")

    with pysam.VariantFile(str(target_vcf_path)) as target_vcf:
        target_vcf.subset_samples([])
        header = target_vcf.header.copy()
        _add_pon_header_fields(header, min_baseq=min_baseq, min_mapq=min_mapq)
        for sample_name in sample_names:
            header.add_sample(sample_name)

        evidence_iterator = iter(evidence_records)
        with pysam.VariantFile(str(output_path), "wb", header=header) as output_file:
            for target_record in target_vcf:
                try:
                    variant, normal_evidences = next(evidence_iterator)
                except StopIteration as exc:
                    raise RuntimeError("PON evidence ended before the target VCF") from exc

                record_variant = Variant.from_vcf_fields(
                    contig=target_record.contig,
                    pos1=target_record.pos,
                    ref=target_record.ref,
                    alt=target_record.alts[0],
                )
                if record_variant != variant:
                    raise RuntimeError("PON evidence variants are out of target VCF order")
                if len(normal_evidences) != len(sample_names):
                    raise RuntimeError("PON evidence sample count does not match its header")

                output_record = _copy_target_record(target_record, output_file)
                for sample_name, evidence in zip(
                    sample_names,
                    normal_evidences,
                    strict=True,
                ):
                    sample = output_record.samples[sample_name]
                    for field_id, attribute in _EVIDENCE_ATTRIBUTES_BY_FIELD.items():
                        sample[field_id] = getattr(evidence, attribute)
                output_file.write(output_record)

            try:
                next(evidence_iterator)
            except StopIteration:
                pass
            else:
                raise RuntimeError("PON evidence contains more variants than the target VCF")

    try:
        bcftools.index("--force", str(output_path))
    except pysam.SamtoolsError as exc:
        raise ValueError(
            "PON targets must be coordinate-sorted so the BCF can be indexed"
        ) from exc


def _evidence_from_sample(sample: Any, *, sample_name: str, variant: Variant) -> AggregatedEvidence:
    values: dict[str, int] = {}
    for field_id, attribute in _EVIDENCE_ATTRIBUTES_BY_FIELD.items():
        value = sample[field_id]
        if value is None:
            raise ValueError(
                f"PON sample {sample_name!r} has missing {field_id} at "
                f"{variant.contig}:{variant.ref_pos0 + 1}"
            )
        integer_value = int(value)
        if integer_value < 0:
            raise ValueError(
                f"PON sample {sample_name!r} has negative {field_id} at "
                f"{variant.contig}:{variant.ref_pos0 + 1}"
            )
        values[attribute] = integer_value

    if values["usable"] != sum(
        values[key]
        for key in ("alt_forward", "alt_reverse", "non_alt_forward", "non_alt_reverse")
    ):
        raise ValueError(
            f"PON sample {sample_name!r} has inconsistent usable counts at "
            f"{variant.contig}:{variant.ref_pos0 + 1}"
        )

    return AggregatedEvidence(**values, unusable_by_reason={})


def read_pon_evidence(
    path: str | Path,
) -> Iterator[tuple[Variant, tuple[AggregatedEvidence, ...]]]:
    """Yield target variants and their per-normal evidence in artifact order."""
    with pysam.VariantFile(str(path)) as pon_file:
        metadata = _parse_metadata(pon_file.header)
        for record in pon_file:
            alts = record.alts or ()
            if len(alts) != 1:
                raise ValueError(
                    f"PON artifact contains a non-biallelic record at {record.contig}:{record.pos}"
                )
            try:
                variant = Variant.from_vcf_fields(
                    contig=record.contig,
                    pos1=record.pos,
                    ref=record.ref,
                    alt=alts[0],
                )
            except ValueError as exc:
                raise ValueError(
                    f"PON artifact contains an unsupported allele at {record.contig}:{record.pos}"
                ) from exc

            yield variant, tuple(
                _evidence_from_sample(
                    record.samples[sample_name],
                    sample_name=sample_name,
                    variant=variant,
                )
                for sample_name in metadata.sample_names
            )
