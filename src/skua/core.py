"""Core public API for skua."""

import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import sys
from typing import Any
from typing import Callable
from typing import Iterable
from typing import Iterator

import pysam

from .evidence import AggregatedEvidence, collect_evidence_from_alignment
from .stats import aggregate_evidence, compute_stats, DEFAULT_TRUNCATE, truncated_normal_evidences
from .variants import Variant, read_vcf_variant_file


READ_COUNT_FORMAT_FIELD_DEFINITIONS: tuple[tuple[str, str], ...] = (
    ("SKUA_ALT_FWD", "Case ALT-supporting forward reads"),
    ("SKUA_ALT_REV", "Case ALT-supporting reverse reads"),
    ("SKUA_NON_ALT_FWD", "Case non-ALT forward reads"),
    ("SKUA_NON_ALT_REV", "Case non-ALT reverse reads"),
    ("SKUA_USABLE", "Case usable reads at this locus"),
    ("SKUA_UNUSABLE", "Case unusable reads at this locus"),
)

MODEL_SCORE_FORMAT_FIELD_DEFINITIONS: tuple[tuple[str, str, str], ...] = (
    ("SKUA_LOG_BAYES_FACTOR", "Float", "Log Bayes factor artifact-vs-variant"),
    ("SKUA_ARTIFACT_POSTERIOR", "Float", "Posterior probability of the artifact model"),
)

PON_INFO_FIELD_DEFINITIONS: tuple[tuple[str, str, str], ...] = (
    ("SKUA_PON_SAMPLE_COUNT", "Integer", "Number of PON samples included after truncation"),
    ("SKUA_PON_ALT_FWD", "Integer", "PON ALT-supporting forward reads after truncation"),
    ("SKUA_PON_ALT_REV", "Integer", "PON ALT-supporting reverse reads after truncation"),
    ("SKUA_PON_NON_ALT_FWD", "Integer", "PON non-ALT forward reads after truncation"),
    ("SKUA_PON_NON_ALT_REV", "Integer", "PON non-ALT reverse reads after truncation"),
    ("SKUA_PON_USABLE", "Integer", "PON usable reads after truncation"),
    ("SKUA_PON_UNUSABLE", "Integer", "PON unusable reads after truncation"),
    ("SKUA_PON_DISPERSION_FACTOR", "Float", "Estimated dispersion factor"),
)

ANNOTATION_STATUS_INFO_FIELD_DEFINITION = (
    "SKUA_STATUS",
    "String",
    "Skua annotation status for this VCF record",
)


class AnnotationStatus(str, Enum):
    """Outcome of attempting to annotate one VCF record."""

    ANNOTATED = "ANNOTATED"
    UNSUPPORTED_MULTIALLELIC = "UNSUPPORTED_MULTIALLELIC"
    UNSUPPORTED_SYMBOLIC_ALLELE = "UNSUPPORTED_SYMBOLIC_ALLELE"
    UNSUPPORTED_BREAKEND = "UNSUPPORTED_BREAKEND"
    UNSUPPORTED_SPANNING_DELETION = "UNSUPPORTED_SPANNING_DELETION"
    UNSUPPORTED_COMPLEX_ALLELE = "UNSUPPORTED_COMPLEX_ALLELE"
    UNSUPPORTED_NON_STANDARD_ALLELE = "UNSUPPORTED_NON_STANDARD_ALLELE"


@dataclass(frozen=True)
class VcfRecordAnnotation:
    """Supported variant or explicit reason why a VCF record was not annotated."""

    status: AnnotationStatus
    variant: Variant | None


@dataclass
class AnnotationSummary:
    """Counts of VCF records processed during one Skua annotation run."""

    record_count: int = 0
    annotated_record_count: int = 0
    unsupported_record_count_by_status: dict[AnnotationStatus, int] = field(default_factory=dict)

    def record(self, status: AnnotationStatus) -> None:
        """Record one VCF annotation outcome."""
        self.record_count += 1
        if status == AnnotationStatus.ANNOTATED:
            self.annotated_record_count += 1
            return
        self.unsupported_record_count_by_status[status] = (
            self.unsupported_record_count_by_status.get(status, 0) + 1
        )

    @property
    def unsupported_record_count(self) -> int:
        """Return the number of records not annotated by Skua."""
        return self.record_count - self.annotated_record_count

    def format_for_cli(self) -> str:
        """Render a concise, deterministic command-line summary."""
        message = (
            f"skua: records={self.record_count} annotated={self.annotated_record_count} "
            f"unsupported={self.unsupported_record_count}"
        )
        if not self.unsupported_record_count_by_status:
            return message
        details = ", ".join(
            f"{status.value}={count}"
            for status, count in sorted(
                self.unsupported_record_count_by_status.items(),
                key=lambda item: item[0].value,
            )
        )
        return f"{message} ({details})"


@dataclass(frozen=True)
class PonAnnotation:
    """Evidence collected for one case variant and its panel of normals.

    ``normal_evidences`` preserves one :class:`AggregatedEvidence` object per
    normal alignment. ``normal_aggregate_evidence`` is their unfiltered sum;
    callers that need a truncated panel can apply ``truncated_normal_evidences``.
    """

    case_evidence: AggregatedEvidence
    normal_evidences: tuple[AggregatedEvidence, ...]
    normal_aggregate_evidence: AggregatedEvidence


@dataclass(frozen=True)
class CaseSampleSelection:
    """Resolved case sample and any read-group restriction needed to isolate it."""

    sample_name: str
    allowed_read_group_ids: frozenset[str] | None


def _alignment_header_dict(alignment_file: Any) -> dict[str, Any] | None:
    """Return an alignment header as a dictionary, when the object exposes one."""
    header = getattr(alignment_file, "header", None)
    if header is None:
        return None

    if hasattr(header, "to_dict"):
        return header.to_dict()
    if isinstance(header, dict):
        return header
    raise ValueError("Alignment file header does not expose read-group metadata")


def _alignment_sample_names(alignment_file: Any) -> tuple[str, ...]:
    """Return distinct read-group sample names in header order."""
    header_dict = _alignment_header_dict(alignment_file)
    if header_dict is None:
        return ()

    sample_names: list[str] = []
    for read_group in header_dict.get("RG", []):
        if not isinstance(read_group, dict):
            continue
        sample_name = read_group.get("SM")
        if sample_name:
            sample_names.append(str(sample_name))
    return tuple(dict.fromkeys(sample_names))


def _read_group_ids_for_sample(alignment_file: Any, sample_name: str) -> frozenset[str]:
    """Return read-group IDs belonging to one alignment sample."""
    header_dict = _alignment_header_dict(alignment_file)
    if header_dict is None:
        return frozenset()

    return frozenset(
        str(read_group["ID"])
        for read_group in header_dict.get("RG", [])
        if isinstance(read_group, dict)
        and read_group.get("SM") == sample_name
        and read_group.get("ID")
    )


def _alignment_sample_name(alignment_file: Any) -> str:
    """Return the single usable read-group sample name for an alignment file."""
    if _alignment_header_dict(alignment_file) is None:
        raise ValueError("Alignment file does not expose a header with read-group sample names")

    sample_names = _alignment_sample_names(alignment_file)
    if not sample_names:
        raise ValueError("Alignment file must contain exactly one usable read-group SM tag")
    if len(sample_names) > 1:
        raise ValueError(
            "Alignment file contains multiple distinct read-group SM tags: "
            + ", ".join(sample_names)
        )
    return sample_names[0]


def _resolve_case_sample(
    vcf_header: Any,
    alignment_file: Any,
    *,
    requested_sample_name: str | None,
) -> CaseSampleSelection:
    """Resolve the VCF case sample and isolate it in multi-sample alignments."""
    vcf_sample_names = tuple(vcf_header.samples)
    alignment_sample_names = _alignment_sample_names(alignment_file)

    def selection_for(sample_name: str) -> CaseSampleSelection:
        if len(alignment_sample_names) <= 1:
            return CaseSampleSelection(sample_name=sample_name, allowed_read_group_ids=None)

        read_group_ids = _read_group_ids_for_sample(alignment_file, sample_name)
        if not read_group_ids:
            raise ValueError(
                f"Case sample {sample_name!r} has no read-group IDs in a multi-sample alignment"
            )
        return CaseSampleSelection(
            sample_name=sample_name,
            allowed_read_group_ids=read_group_ids,
        )

    if not vcf_sample_names:
        if requested_sample_name is not None:
            if requested_sample_name not in alignment_sample_names:
                raise ValueError(
                    f"Requested sample {requested_sample_name!r} is not present in the case alignment"
                )
            return selection_for(requested_sample_name)
        if len(alignment_sample_names) == 1:
            return selection_for(alignment_sample_names[0])
        if not alignment_sample_names:
            raise ValueError("Site-only VCF input requires a usable read-group SM tag or --sample")
        raise ValueError("Case alignment contains multiple samples; specify --sample")

    if requested_sample_name is not None:
        if requested_sample_name not in vcf_sample_names:
            raise ValueError(f"Requested sample {requested_sample_name!r} is not present in the VCF")
        if alignment_sample_names and requested_sample_name not in alignment_sample_names:
            raise ValueError(
                f"Requested sample {requested_sample_name!r} is not present in the case alignment"
            )
        if not alignment_sample_names and len(vcf_sample_names) > 1:
            raise ValueError(
                "Case alignment has no usable read-group SM tag to select among VCF samples"
            )
        return selection_for(requested_sample_name)

    matching_sample_names = tuple(
        sample_name for sample_name in vcf_sample_names if sample_name in alignment_sample_names
    )
    if len(matching_sample_names) == 1:
        return selection_for(matching_sample_names[0])
    if len(matching_sample_names) > 1:
        raise ValueError("Multiple case alignment samples match the VCF; specify --sample")
    if len(vcf_sample_names) == 1 and not alignment_sample_names:
        return selection_for(vcf_sample_names[0])
    if len(vcf_sample_names) == 1:
        raise ValueError(
            "The sole VCF sample does not match a usable read-group SM tag in the case alignment"
        )
    raise ValueError("No case alignment sample matches the VCF; specify --sample")


def _validate_normal_alignment_samples(normal_alignments: list[Any]) -> None:
    """Require one read-group sample per normal alignment when metadata is available."""
    for index, normal_alignment in enumerate(normal_alignments, start=1):
        if _alignment_header_dict(normal_alignment) is None:
            continue
        try:
            _alignment_sample_name(normal_alignment)
        except ValueError as exc:
            raise ValueError(f"Normal alignment {index}: {exc}") from exc


def _validate_annotation_parameters(
    *,
    min_baseq: int,
    min_mapq: int,
    truncate: float | None = None,
    pseudocount: float | None = None,
    prior_variant_probability: float | None = None,
) -> None:
    """Reject parameter values whose semantics are undefined for annotation."""
    if min_baseq < 0:
        raise ValueError("min_baseq must be >= 0")
    if min_mapq < 0:
        raise ValueError("min_mapq must be >= 0")
    if truncate is not None and not 0.0 < truncate <= 1.0:
        raise ValueError("truncate must be greater than 0 and no greater than 1")
    if pseudocount is not None and pseudocount <= 0:
        raise ValueError("pseudocount must be > 0")
    if prior_variant_probability is not None and not 0.0 < prior_variant_probability < 1.0:
        raise ValueError("prior_variant_probability must be between 0 and 1")


def _validate_alignment_indexes(alignment_files: list[tuple[str, Any]]) -> None:
    """Fail before output when an alignment exposes an unavailable index."""
    for label, alignment_file in alignment_files:
        has_index = getattr(alignment_file, "has_index", None)
        if has_index is not None and not has_index():
            raise ValueError(f"{label} must be indexed")


def _validate_vcf_against_inputs(
    vcf_path: str | Path,
    *,
    alignment_files: list[tuple[str, Any]],
    reference_path: str | Path | None,
    strict: bool = False,
) -> None:
    """Validate supported VCF records against alignment contigs and an optional FASTA."""
    _validate_alignment_indexes(alignment_files)
    alignment_contigs = [
        (
            label,
            frozenset(contigs) if contigs is not None else None,
        )
        for label, alignment_file in alignment_files
        for contigs in (getattr(alignment_file, "references", None),)
    ]

    fasta_file: Any | None = None
    if reference_path is not None:
        fasta_file = pysam.FastaFile(str(reference_path))

    try:
        with pysam.VariantFile(str(vcf_path)) as source_vcf:
            for record in source_vcf:
                assessment = _assess_vcf_record(record)
                if strict and assessment.status != AnnotationStatus.ANNOTATED:
                    raise ValueError(
                        f"Unsupported VCF record at {record.contig}:{record.pos}: "
                        f"{assessment.status.value}"
                    )
                variant = assessment.variant
                if variant is None:
                    continue

                for label, contigs in alignment_contigs:
                    if contigs is not None and variant.contig not in contigs:
                        raise ValueError(f"{label} does not contain contig {variant.contig!r}")

                if fasta_file is None:
                    continue
                if variant.contig not in fasta_file.references:
                    raise ValueError(f"Reference FASTA does not contain contig {variant.contig!r}")

                reference_bases = fasta_file.fetch(
                    variant.contig,
                    variant.ref_pos0,
                    variant.ref_pos0 + len(variant.ref),
                ).upper()
                if reference_bases != variant.ref.upper():
                    raise ValueError(
                        f"VCF REF allele at {variant.contig}:{variant.ref_pos0 + 1} "
                        f"is {variant.ref!r}, but the reference FASTA contains {reference_bases!r}"
                    )
    finally:
        if fasta_file is not None:
            fasta_file.close()


def _ensure_skua_vcf_header_fields(header: Any, *, include_pon_info: bool) -> Any:
    """Ensure SKUA FORMAT/INFO definitions exist on the active VCF header."""
    annotated_header = header

    for field_id, description in READ_COUNT_FORMAT_FIELD_DEFINITIONS:
        if field_id not in annotated_header.formats:
            annotated_header.add_line(
                f'##FORMAT=<ID={field_id},Number=1,Type=Integer,Description="{description}">'
            )

    status_field_id, status_field_type, status_description = ANNOTATION_STATUS_INFO_FIELD_DEFINITION
    if status_field_id not in annotated_header.info:
        annotated_header.add_line(
            f'##INFO=<ID={status_field_id},Number=1,Type={status_field_type},'
            f'Description="{status_description}">'
        )

    if include_pon_info:
        for field_id, field_type, description in MODEL_SCORE_FORMAT_FIELD_DEFINITIONS:
            if field_id not in annotated_header.formats:
                annotated_header.add_line(
                    f'##FORMAT=<ID={field_id},Number=1,Type={field_type},Description="{description}">'
                )

        for field_id, field_type, description in PON_INFO_FIELD_DEFINITIONS:
            if field_id not in annotated_header.info:
                annotated_header.add_line(
                    f'##INFO=<ID={field_id},Number=1,Type={field_type},Description="{description}">'
                )

    return annotated_header


def _assess_vcf_record(record: Any) -> VcfRecordAnnotation:
    """Return a supported variant or an explicit unsupported-record status."""
    alts = record.alts or ()
    if len(alts) != 1:
        return VcfRecordAnnotation(AnnotationStatus.UNSUPPORTED_MULTIALLELIC, None)

    alt = alts[0]
    if alt == "*":
        return VcfRecordAnnotation(AnnotationStatus.UNSUPPORTED_SPANNING_DELETION, None)
    if alt.startswith("<") and alt.endswith(">"):
        return VcfRecordAnnotation(AnnotationStatus.UNSUPPORTED_SYMBOLIC_ALLELE, None)
    if "[" in alt or "]" in alt:
        return VcfRecordAnnotation(AnnotationStatus.UNSUPPORTED_BREAKEND, None)
    if any(base not in {"A", "C", "G", "T"} for base in record.ref.upper() + alt.upper()):
        return VcfRecordAnnotation(AnnotationStatus.UNSUPPORTED_NON_STANDARD_ALLELE, None)

    try:
        variant = Variant.from_vcf_fields(
            contig=record.contig,
            pos1=record.pos,
            ref=record.ref,
            alt=alt,
        )
    except ValueError:
        return VcfRecordAnnotation(AnnotationStatus.UNSUPPORTED_COMPLEX_ALLELE, None)
    return VcfRecordAnnotation(AnnotationStatus.ANNOTATED, variant)


def _variant_from_vcf_record(record: Any) -> Variant | None:
    """Build a Variant from a pysam VCF record when supported, else None."""
    return _assess_vcf_record(record).variant


def _copy_vcf_record_with_sample(record: Any, out_vcf: Any) -> Any:
    """Copy a site-only VCF record into an output header that has one sample."""
    copied_record = out_vcf.new_record(
        contig=record.contig,
        start=record.start,
        stop=record.stop,
        id=record.id,
        alleles=record.alleles,
        qual=record.qual,
    )
    for filter_id in record.filter.keys():
        copied_record.filter.add(filter_id)
    for key, value in record.info.items():
        copied_record.info[key] = value
    return copied_record


def _annotate_read_count_format_fields(
    record: Any,
    evidence: AggregatedEvidence,
    *,
    sample_name: str,
) -> None:
    """Set read-count FORMAT annotations for the selected case sample."""
    sample = record.samples[sample_name]
    sample["SKUA_ALT_FWD"] = evidence.alt_forward
    sample["SKUA_ALT_REV"] = evidence.alt_reverse
    sample["SKUA_NON_ALT_FWD"] = evidence.non_alt_forward
    sample["SKUA_NON_ALT_REV"] = evidence.non_alt_reverse
    sample["SKUA_USABLE"] = evidence.usable
    sample["SKUA_UNUSABLE"] = evidence.unusable


def _annotate_pon_sample_format_fields(
    record: Any,
    *,
    sample_name: str,
    artifact_posterior: float,
    log_bayes_factor: float,
) -> None:
    """Set PON model output FORMAT annotations for the selected case sample."""
    sample = record.samples[sample_name]
    sample["SKUA_LOG_BAYES_FACTOR"] = float(log_bayes_factor)
    sample["SKUA_ARTIFACT_POSTERIOR"] = float(artifact_posterior)


def _vcf_write_mode(output_path: str | Path) -> str:
    """Return the pysam VariantFile write mode for VCF output path."""
    if str(output_path).endswith(".gz"):
        return "wz"
    return "w"


def annotate_vcf(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    output_path: str | Path,
    sample_name: str | None = None,
    reference_path: str | Path | None = None,
    strict: bool = False,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> AnnotationSummary:
    """Annotate an input VCF with read-count FORMAT fields for variants."""
    _validate_annotation_parameters(min_baseq=min_baseq, min_mapq=min_mapq)
    _validate_vcf_against_inputs(
        vcf_path,
        alignment_files=[("Case alignment", alignment_file)],
        reference_path=reference_path,
        strict=strict,
    )

    summary = AnnotationSummary()
    with pysam.VariantFile(str(vcf_path)) as source_vcf:
        header = _ensure_skua_vcf_header_fields(source_vcf.header, include_pon_info=False)
        case_selection = _resolve_case_sample(
            source_vcf.header,
            alignment_file,
            requested_sample_name=sample_name,
        )
        site_only_sample_name: str | None = None
        if len(source_vcf.header.samples) == 0:
            site_only_sample_name = case_selection.sample_name
            header.add_sample(site_only_sample_name)

        with pysam.VariantFile(
            str(output_path),
            _vcf_write_mode(output_path),
            header=header,
        ) as out_vcf:
            for record in source_vcf:
                if site_only_sample_name is not None:
                    record = _copy_vcf_record_with_sample(record, out_vcf)
                assessment = _assess_vcf_record(record)
                summary.record(assessment.status)
                record.info["SKUA_STATUS"] = assessment.status.value
                if assessment.variant is not None:
                    evidence = annotate_variant(
                        alignment_file,
                        assessment.variant,
                        min_baseq=min_baseq,
                        min_mapq=min_mapq,
                        allowed_read_group_ids=case_selection.allowed_read_group_ids,
                    )
                    _annotate_read_count_format_fields(
                        record,
                        evidence,
                        sample_name=case_selection.sample_name,
                    )
                out_vcf.write(record)

    return summary


def annotate_vcf_with_normals(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    normal_alignments: list[Any] | None = None,
    output_path: str | Path,
    sample_name: str | None = None,
    reference_path: str | Path | None = None,
    strict: bool = False,
    min_baseq: int = 20,
    min_mapq: int = 20,
    truncate: float = DEFAULT_TRUNCATE,
    pseudocount: float = sys.float_info.epsilon,
    prior_variant_probability: float = 0.5,
) -> AnnotationSummary:
    """Annotate an input VCF with read-count FORMAT and PON INFO fields."""
    if normal_alignments is None:
        normal_alignments = []

    _validate_annotation_parameters(
        min_baseq=min_baseq,
        min_mapq=min_mapq,
        truncate=truncate,
        pseudocount=pseudocount,
        prior_variant_probability=prior_variant_probability,
    )
    _validate_normal_alignment_samples(normal_alignments)
    _validate_vcf_against_inputs(
        vcf_path,
        alignment_files=[("Case alignment", alignment_file)]
        + [
            (f"Normal alignment {index}", normal_alignment)
            for index, normal_alignment in enumerate(normal_alignments, start=1)
        ],
        reference_path=reference_path,
        strict=strict,
    )

    summary = AnnotationSummary()
    with pysam.VariantFile(str(vcf_path)) as source_vcf:
        header = _ensure_skua_vcf_header_fields(source_vcf.header, include_pon_info=True)
        case_selection = _resolve_case_sample(
            source_vcf.header,
            alignment_file,
            requested_sample_name=sample_name,
        )
        site_only_sample_name: str | None = None
        if len(source_vcf.header.samples) == 0:
            site_only_sample_name = case_selection.sample_name
            header.add_sample(site_only_sample_name)

        with pysam.VariantFile(
            str(output_path),
            _vcf_write_mode(output_path),
            header=header,
        ) as out_vcf:
            for record in source_vcf:
                if site_only_sample_name is not None:
                    record = _copy_vcf_record_with_sample(record, out_vcf)
                assessment = _assess_vcf_record(record)
                summary.record(assessment.status)
                record.info["SKUA_STATUS"] = assessment.status.value
                if assessment.variant is not None:
                    pon_result = annotate_variant_with_normals(
                        alignment_file,
                        assessment.variant,
                        normal_alignments=normal_alignments,
                        min_baseq=min_baseq,
                        min_mapq=min_mapq,
                        allowed_read_group_ids=case_selection.allowed_read_group_ids,
                    )
                    case_evidence = pon_result.case_evidence
                    normal_samples_included = truncated_normal_evidences(
                        list(pon_result.normal_evidences),
                        truncate=truncate,
                    )
                    normal_output_evidence = aggregate_evidence(normal_samples_included)
                    stats = compute_stats(
                        case_evidence,
                        normal_output_evidence,
                        per_sample_evidences=list(pon_result.normal_evidences),
                        truncate=truncate,
                        pseudocount=pseudocount,
                        prior_variant_probability=prior_variant_probability,
                    )

                    _annotate_read_count_format_fields(
                        record,
                        case_evidence,
                        sample_name=case_selection.sample_name,
                    )
                    _annotate_pon_sample_format_fields(
                        record,
                        sample_name=case_selection.sample_name,
                        artifact_posterior=stats.artifact_posterior,
                        log_bayes_factor=stats.log_bayes_factor_artifact_vs_variant,
                    )
                    record.info["SKUA_PON_SAMPLE_COUNT"] = len(normal_samples_included)
                    record.info["SKUA_PON_ALT_FWD"] = normal_output_evidence.alt_forward
                    record.info["SKUA_PON_ALT_REV"] = normal_output_evidence.alt_reverse
                    record.info["SKUA_PON_NON_ALT_FWD"] = normal_output_evidence.non_alt_forward
                    record.info["SKUA_PON_NON_ALT_REV"] = normal_output_evidence.non_alt_reverse
                    record.info["SKUA_PON_USABLE"] = normal_output_evidence.usable
                    record.info["SKUA_PON_UNUSABLE"] = normal_output_evidence.unusable
                    record.info["SKUA_PON_DISPERSION_FACTOR"] = float(stats.dispersion_rho)

                out_vcf.write(record)

    return summary


def annotate_variant(
    alignment_file: Any,
    variant: Variant,
    *,
    min_baseq: int = 20,
    min_mapq: int = 20,
    allowed_read_group_ids: frozenset[str] | None = None,
) -> AggregatedEvidence:
    """Collect strand-aware evidence for one variant from one alignment."""
    return collect_evidence_from_alignment(
        alignment_file,
        contig=variant.contig,
        ref_pos0=variant.ref_pos0,
        ref_base=variant.ref,
        alt_base=variant.alt,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
        allowed_read_group_ids=allowed_read_group_ids,
    )


def annotate_variant_with_normals(
    alignment_file: Any,
    variant: Variant,
    *,
    normal_alignments: list[Any] | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
    allowed_read_group_ids: frozenset[str] | None = None,
) -> PonAnnotation:
    """Collect case and normal evidence for one variant."""
    if normal_alignments is None:
        normal_alignments = []

    case_evidence = annotate_variant(
        alignment_file,
        variant,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
        allowed_read_group_ids=allowed_read_group_ids,
    )

    normal_evidences: list[AggregatedEvidence] = []
    normal_aggregate_evidence = AggregatedEvidence(
        alt_forward=0,
        alt_reverse=0,
        non_alt_forward=0,
        non_alt_reverse=0,
        usable=0,
        unusable=0,
        unusable_by_reason={},
    )

    for normal_alignment in normal_alignments:
        normal_evidence = annotate_variant(
            normal_alignment,
            variant,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        )
        normal_evidences.append(normal_evidence)

        normal_unusable_by_reason = dict(normal_aggregate_evidence.unusable_by_reason)
        for reason, count in normal_evidence.unusable_by_reason.items():
            normal_unusable_by_reason[reason] = normal_unusable_by_reason.get(reason, 0) + count

        normal_aggregate_evidence = AggregatedEvidence(
            alt_forward=normal_aggregate_evidence.alt_forward + normal_evidence.alt_forward,
            alt_reverse=normal_aggregate_evidence.alt_reverse + normal_evidence.alt_reverse,
            non_alt_forward=normal_aggregate_evidence.non_alt_forward + normal_evidence.non_alt_forward,
            non_alt_reverse=normal_aggregate_evidence.non_alt_reverse + normal_evidence.non_alt_reverse,
            usable=normal_aggregate_evidence.usable + normal_evidence.usable,
            unusable=normal_aggregate_evidence.unusable + normal_evidence.unusable,
            unusable_by_reason=normal_unusable_by_reason,
        )

    return PonAnnotation(
        case_evidence=case_evidence,
        normal_evidences=tuple(normal_evidences),
        normal_aggregate_evidence=normal_aggregate_evidence,
    )


def annotate_variants_from_vcf(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> Iterator[tuple[Variant, AggregatedEvidence]]:
    """Yield per-variant evidence for variant records from a VCF file."""
    for variant in read_vcf_variant_file(vcf_path):
        yield (
            variant,
            annotate_variant(
                alignment_file,
                variant,
                min_baseq=min_baseq,
                min_mapq=min_mapq,
            ),
        )


def format_annotation_results(
    results: Iterable[tuple[Variant, AggregatedEvidence]],
) -> list[dict[str, Any]]:
    """Convert annotation results to JSON/tabular-ready row dictionaries."""
    rows: list[dict[str, Any]] = []
    for variant, evidence in results:
        rows.append(
            {
                "contig": variant.contig,
                "pos1": variant.ref_pos0 + 1,
                "ref": variant.ref,
                "alt": variant.alt,
                "counts": {
                    "case": {
                        "alt_forward": evidence.alt_forward,
                        "alt_reverse": evidence.alt_reverse,
                        "non_alt_forward": evidence.non_alt_forward,
                        "non_alt_reverse": evidence.non_alt_reverse,
                        "usable": evidence.usable,
                        "unusable": evidence.unusable,
                        "unusable_by_reason": {
                            reason.value: count
                            for reason, count in evidence.unusable_by_reason.items()
                        },
                    },
                },
            }
        )
    return rows


def render_annotation_results_json(rows: Iterable[dict[str, Any]]) -> str:
    """Render formatted annotation rows as JSON text."""
    return json.dumps(list(rows), indent=2)


def write_annotation_results_json(
    rows: Iterable[dict[str, Any]],
    output_path: str | Path,
) -> None:
    """Write formatted annotation rows to a JSON file."""
    Path(output_path).write_text(
        render_annotation_results_json(rows),
        encoding="utf-8",
    )


def _build_annotation_rows(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    min_baseq: int,
    min_mapq: int,
) -> list[dict[str, Any]]:
    """Build formatted annotation rows from one alignment and one VCF."""
    return format_annotation_results(
        annotate_variants_from_vcf(
            alignment_file,
            vcf_path,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        )
    )


def _render_and_optionally_write(
    rows: Iterable[dict[str, Any]],
    *,
    renderer: Callable[[Iterable[dict[str, Any]]], str],
    output_path: str | Path | None,
) -> str:
    """Render rows and optionally persist the payload to disk."""
    payload = renderer(rows)
    if output_path is not None:
        Path(output_path).write_text(payload, encoding="utf-8")
    return payload


def annotate_vcf_to_json(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    output_path: str | Path | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> str:
    """Run variant annotation from VCF and return JSON output, optionally writing to file."""
    rows = _build_annotation_rows(
        alignment_file,
        vcf_path,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
    )
    return _render_and_optionally_write(
        rows,
        renderer=render_annotation_results_json,
        output_path=output_path,
    )

def annotate_variants_from_vcf_with_normals(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    normal_alignments: list[Any] | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> Iterator[tuple[Variant, PonAnnotation]]:
    """Yield per-variant case+normal evidence for variant records from a VCF file."""
    if normal_alignments is None:
        normal_alignments = []

    for variant in read_vcf_variant_file(vcf_path):
        yield (
            variant,
            annotate_variant_with_normals(
                alignment_file,
                variant,
                normal_alignments=normal_alignments,
                min_baseq=min_baseq,
                min_mapq=min_mapq,
            ),
        )


def format_annotation_results_with_normals(
    results: Iterable[tuple[Variant, PonAnnotation]],
    *,
    truncate: float = DEFAULT_TRUNCATE,
    pseudocount: float = sys.float_info.epsilon,
    prior_variant_probability: float = 0.5,
) -> list[dict[str, Any]]:
    """Convert PON annotation results to JSON/tabular-ready row dictionaries."""
    rows: list[dict[str, Any]] = []
    for variant, pon_result in results:
        evidence = pon_result.case_evidence
        per_sample_evidences = list(pon_result.normal_evidences)

        normal_samples_included = truncated_normal_evidences(
            per_sample_evidences,
            truncate=truncate,
        )
        normal_output_evidence = aggregate_evidence(normal_samples_included)

        stats = compute_stats(
            evidence,
            normal_output_evidence,
            per_sample_evidences=per_sample_evidences,
            truncate=truncate,
            pseudocount=pseudocount,
            prior_variant_probability=prior_variant_probability,
        )
        normal_samples_used = len(normal_samples_included)
        rows.append(
            {
                "contig": variant.contig,
                "pos1": variant.ref_pos0 + 1,
                "ref": variant.ref,
                "alt": variant.alt,
                "stats": {
                    "artifact_posterior": stats.artifact_posterior,
                    "log_bayes_factor_artifact_vs_variant": stats.log_bayes_factor_artifact_vs_variant,
                    "dispersion_factor": stats.dispersion_rho,
                    "pon_sample_count": normal_samples_used,
                },
                "counts": {
                    "case": {
                        "alt_forward": evidence.alt_forward,
                        "alt_reverse": evidence.alt_reverse,
                        "non_alt_forward": evidence.non_alt_forward,
                        "non_alt_reverse": evidence.non_alt_reverse,
                        "usable": evidence.usable,
                        "unusable": evidence.unusable,
                        "unusable_by_reason": {
                            reason.value: count
                            for reason, count in evidence.unusable_by_reason.items()
                        },
                    },
                    "normal": {
                        "alt_forward": normal_output_evidence.alt_forward,
                        "alt_reverse": normal_output_evidence.alt_reverse,
                        "non_alt_forward": normal_output_evidence.non_alt_forward,
                        "non_alt_reverse": normal_output_evidence.non_alt_reverse,
                        "usable": normal_output_evidence.usable,
                        "unusable": normal_output_evidence.unusable,
                        "unusable_by_reason": {
                            reason.value: count
                            for reason, count in normal_output_evidence.unusable_by_reason.items()
                        },
                    },
                },
            }
        )
    return rows


def _build_annotation_rows_with_normals(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    normal_alignments: list[Any] | None,
    min_baseq: int,
    min_mapq: int,
    truncate: float,
    pseudocount: float,
    prior_variant_probability: float,
) -> list[dict[str, Any]]:
    """Build formatted PON annotation rows from case + normal alignments and one VCF."""
    return format_annotation_results_with_normals(
        annotate_variants_from_vcf_with_normals(
            alignment_file,
            vcf_path,
            normal_alignments=normal_alignments,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        ),
        truncate=truncate,
        pseudocount=pseudocount,
        prior_variant_probability=prior_variant_probability,
    )


def annotate_vcf_to_json_with_normals(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    normal_alignments: list[Any] | None = None,
    output_path: str | Path | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
    truncate: float = DEFAULT_TRUNCATE,
    pseudocount: float = sys.float_info.epsilon,
    prior_variant_probability: float = 0.5,
) -> str:
    """Run PON variant annotation from VCF and return JSON output, optionally writing to file."""
    if normal_alignments is None:
        normal_alignments = []

    rows = _build_annotation_rows_with_normals(
        alignment_file,
        vcf_path,
        normal_alignments=normal_alignments,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
        truncate=truncate,
        pseudocount=pseudocount,
        prior_variant_probability=prior_variant_probability,
    )
    return _render_and_optionally_write(
        rows,
        renderer=render_annotation_results_json,
        output_path=output_path,
    )
