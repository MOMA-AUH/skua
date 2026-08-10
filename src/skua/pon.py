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
EVIDENCE_POLICY_VERSION = 2
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


@dataclass(frozen=True)
class PonInspection:
    """Header-level facts reported by ``skua pon inspect``."""

    path: str
    format: str
    index_present: bool
    metadata_record_count: int
    schema_version: str | None
    evidence_policy_version: str | None
    min_baseq: str | None
    min_mapq: str | None
    skua_version: str | None
    sample_names: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        """Return a JSON-ready representation of this inspection."""
        return {
            "path": self.path,
            "format": self.format,
            "index_present": self.index_present,
            "metadata_record_count": self.metadata_record_count,
            "schema_version": self.schema_version,
            "evidence_policy_version": self.evidence_policy_version,
            "min_baseq": self.min_baseq,
            "min_mapq": self.min_mapq,
            "skua_version": self.skua_version,
            "sample_count": len(self.sample_names),
            "sample_names": list(self.sample_names),
        }


@dataclass(frozen=True)
class PonValidationResult:
    """Outcome reported by ``skua pon validate``."""

    inspection: PonInspection | None
    errors: tuple[str, ...]

    @property
    def valid(self) -> bool:
        """Return whether the artifact passed every requested validation."""
        return not self.errors

    def as_dict(self) -> dict[str, object]:
        """Return a JSON-ready representation of this validation result."""
        return {
            "valid": self.valid,
            "errors": list(self.errors),
            "inspection": None if self.inspection is None else self.inspection.as_dict(),
        }


def _unquote_header_value(value: Any) -> str:
    text = str(value)
    if len(text) >= 2 and text[0] == text[-1] == '"':
        return text[1:-1]
    return text


def _metadata_records(header: Any) -> list[Any]:
    """Return every PON provenance record found in a VCF header."""
    return [record for record in header.records if record.key == PON_HEADER_KEY]


def _metadata_items_from_record(record: Any) -> dict[str, str]:
    """Return normalized metadata values from one PON provenance record."""
    return {
        key: _unquote_header_value(value)
        for key, value in record.items()
        if key != "IDX"
    }


def _metadata_items(header: Any) -> dict[str, str]:
    records = _metadata_records(header)
    if len(records) != 1:
        raise ValueError(
            f"PON artifact must contain exactly one {PON_HEADER_KEY} metadata record"
        )
    return _metadata_items_from_record(records[0])


def inspect_pon(path: str | Path) -> PonInspection:
    """Read PON header metadata without requiring a supported artifact schema."""
    path_obj = Path(path)
    with pysam.VariantFile(str(path_obj)) as pon_file:
        metadata_records = _metadata_records(pon_file.header)
        metadata = (
            _metadata_items_from_record(metadata_records[0])
            if len(metadata_records) == 1
            else {}
        )
        return PonInspection(
            path=str(path_obj),
            format=pon_file.format,
            index_present=Path(f"{path_obj}.csi").is_file(),
            metadata_record_count=len(metadata_records),
            schema_version=metadata.get("SchemaVersion"),
            evidence_policy_version=metadata.get("EvidencePolicyVersion"),
            min_baseq=metadata.get("MinBaseQ"),
            min_mapq=metadata.get("MinMapQ"),
            skua_version=metadata.get("SkuaVersion"),
            sample_names=tuple(pon_file.header.samples),
        )


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
    if min_baseq < 0 or min_mapq < 0:
        raise ValueError("PON artifact contains negative evidence thresholds")

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
        if pon_file.format != "BCF":
            raise ValueError("PON artifact must be BCF")
        return _parse_metadata(pon_file.header)


def _variant_from_record(record: Any, *, source_label: str) -> Variant:
    """Build a supported PON target variant or raise a descriptive error."""
    alts = record.alts or ()
    if len(alts) != 1:
        raise ValueError(
            f"{source_label} contains a non-biallelic record at {record.contig}:{record.pos}"
        )
    alt = alts[0]
    if any(base not in {"A", "C", "G", "T"} for base in record.ref.upper() + alt.upper()):
        raise ValueError(
            f"{source_label} contains a non-standard allele at {record.contig}:{record.pos}"
        )
    try:
        return Variant.from_vcf_fields(
            contig=record.contig,
            pos1=record.pos,
            ref=record.ref,
            alt=alt,
        )
    except ValueError as exc:
        raise ValueError(
            f"{source_label} contains an unsupported allele at {record.contig}:{record.pos}"
        ) from exc


def _validate_index(path: Path, errors: list[str]) -> None:
    """Require a readable CSI index for a PON artifact."""
    if not Path(f"{path}.csi").is_file():
        errors.append("PON artifact is missing its .csi index")
        return
    try:
        with pysam.VariantFile(str(path)) as pon_file:
            next(pon_file.fetch(), None)
    except (OSError, ValueError, pysam.SamtoolsError) as exc:
        errors.append(f"PON artifact has an unreadable .csi index: {exc}")


def _validate_header(header: Any, errors: list[str]) -> PonArtifactMetadata | None:
    """Validate all header invariants needed to interpret PON evidence."""
    try:
        metadata = _parse_metadata(header)
    except ValueError as exc:
        errors.append(str(exc))
        return None

    if len(set(metadata.sample_names)) != len(metadata.sample_names):
        errors.append("PON artifact contains duplicate normal sample names")

    for field_id, _description in PON_EVIDENCE_FORMAT_FIELDS:
        field = header.formats[field_id]
        if field.number != 1 or field.type != "Integer":
            errors.append(
                f"PON artifact has an incompatible {field_id} FORMAT definition; "
                "expected Number=1,Type=Integer"
            )

    return metadata


def _validate_reference(
    reference_path: str | Path,
    variants: tuple[Variant, ...],
    errors: list[str],
) -> None:
    """Check every PON target REF allele against an optional reference FASTA."""
    try:
        with pysam.FastaFile(str(reference_path)) as reference_file:
            reference_contigs = frozenset(reference_file.references)
            for variant in variants:
                if variant.contig not in reference_contigs:
                    errors.append(
                        f"Reference FASTA does not contain contig {variant.contig!r}"
                    )
                    continue
                observed_ref = reference_file.fetch(
                    variant.contig,
                    variant.ref_pos0,
                    variant.ref_pos0 + len(variant.ref),
                ).upper()
                if observed_ref != variant.ref:
                    errors.append(
                        "PON REF allele at "
                        f"{variant.contig}:{variant.ref_pos0 + 1} is {variant.ref!r}, "
                        f"but the reference FASTA contains {observed_ref!r}"
                    )
    except (OSError, ValueError) as exc:
        errors.append(f"Could not read reference FASTA: {exc}")


def _target_variants(path: str | Path) -> tuple[Variant, ...]:
    """Read exactly the supported target alleles expected from a target VCF."""
    variants: list[Variant] = []
    seen_variants: set[Variant] = set()
    with pysam.VariantFile(str(path)) as target_vcf:
        for record in target_vcf:
            variant = _variant_from_record(record, source_label="Target VCF")
            if variant in seen_variants:
                raise ValueError(
                    "Target VCF contains duplicate target allele "
                    f"{variant.contig}:{variant.ref_pos0 + 1} {variant.ref}>{variant.alt}"
                )
            seen_variants.add(variant)
            variants.append(variant)
    return tuple(variants)


def validate_pon(
    path: str | Path,
    *,
    reference_path: str | Path | None = None,
    target_vcf_path: str | Path | None = None,
) -> PonValidationResult:
    """Validate a PON artifact and optional reference and target-VCF compatibility."""
    try:
        inspection = inspect_pon(path)
    except (OSError, ValueError) as exc:
        return PonValidationResult(inspection=None, errors=(str(exc),))

    errors: list[str] = []
    path_obj = Path(path)
    variants: list[Variant] = []

    if inspection.format != "BCF":
        errors.append("PON artifact must be BCF")
        return PonValidationResult(inspection=inspection, errors=tuple(errors))

    _validate_index(path_obj, errors)
    try:
        with pysam.VariantFile(str(path_obj)) as pon_file:
            metadata = _validate_header(pon_file.header, errors)
            sample_names = () if metadata is None else metadata.sample_names
            seen_variants: set[Variant] = set()
            previous_location: tuple[int, int] | None = None
            contig_order = {
                contig: index for index, contig in enumerate(pon_file.header.contigs)
            }

            for record in pon_file:
                location = (contig_order.get(record.contig, len(contig_order)), record.pos)
                if previous_location is not None and location < previous_location:
                    errors.append("PON artifact records are not coordinate-sorted")
                previous_location = location

                try:
                    variant = _variant_from_record(record, source_label="PON artifact")
                except ValueError as exc:
                    errors.append(str(exc))
                    continue

                if variant in seen_variants:
                    errors.append(
                        "PON artifact contains duplicate evidence for "
                        f"{variant.contig}:{variant.ref_pos0 + 1} {variant.ref}>{variant.alt}"
                    )
                seen_variants.add(variant)
                variants.append(variant)

                for sample_name in sample_names:
                    try:
                        _evidence_from_sample(
                            record.samples[sample_name],
                            sample_name=sample_name,
                            variant=variant,
                        )
                    except (KeyError, TypeError, ValueError) as exc:
                        errors.append(str(exc))

    except (OSError, ValueError, pysam.SamtoolsError) as exc:
        errors.append(f"Could not read PON artifact: {exc}")

    if not variants:
        errors.append("PON artifact must contain at least one target record")

    variant_tuple = tuple(variants)
    if reference_path is not None and variant_tuple:
        _validate_reference(reference_path, variant_tuple, errors)

    if target_vcf_path is not None:
        try:
            if variant_tuple != _target_variants(target_vcf_path):
                errors.append("PON artifact targets do not match the supplied target VCF")
        except (OSError, ValueError) as exc:
            errors.append(f"Could not validate target VCF: {exc}")

    return PonValidationResult(inspection=inspection, errors=tuple(errors))


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
        output_path_obj = Path(output_path)
        if str(output_path) != "-":
            output_path_obj.unlink(missing_ok=True)
            Path(f"{output_path_obj}.csi").unlink(missing_ok=True)
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
        if pon_file.format != "BCF":
            raise ValueError("PON artifact must be BCF")
        metadata = _parse_metadata(pon_file.header)
        for record in pon_file:
            variant = _variant_from_record(record, source_label="PON artifact")

            yield variant, tuple(
                _evidence_from_sample(
                    record.samples[sample_name],
                    sample_name=sample_name,
                    variant=variant,
                )
                for sample_name in metadata.sample_names
            )
