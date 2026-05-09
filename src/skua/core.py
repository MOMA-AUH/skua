"""Core public API for skua."""

import gzip
import json
from pathlib import Path
import sys
import tempfile
from typing import Any
from typing import Callable
from typing import Iterable
from typing import Iterator

import pysam

from .evidence import AggregatedEvidence, collect_snv_evidence_from_alignment
from .stats import aggregate_evidence, compute_stats, DEFAULT_TRUNCATE, truncated_normal_evidences
from .variants import Variant, read_vcf_variant_file


CASE_FORMAT_FIELD_DEFINITIONS: tuple[tuple[str, str], ...] = (
    ("SKUA_ALT_FWD", "Case ALT-supporting forward reads"),
    ("SKUA_ALT_REV", "Case ALT-supporting reverse reads"),
    ("SKUA_NON_ALT_FWD", "Case non-ALT forward reads"),
    ("SKUA_NON_ALT_REV", "Case non-ALT reverse reads"),
    ("SKUA_USABLE", "Case usable reads at this locus"),
    ("SKUA_UNUSABLE", "Case unusable reads at this locus"),
)

PON_FORMAT_FIELD_DEFINITIONS: tuple[tuple[str, str, str], ...] = (
    ("SKUA_ARTIFACT_POSTERIOR", "Float", "Posterior probability of the artifact model"),
    ("SKUA_BAYES_FACTOR", "Float", "Bayes factor artifact-vs-variant"),
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


def _ensure_skua_vcf_header_fields(header: Any, *, include_pon_info: bool) -> Any:
    """Ensure SKUA FORMAT/INFO definitions exist on the active VCF header."""
    annotated_header = header

    for field_id, description in CASE_FORMAT_FIELD_DEFINITIONS:
        if field_id not in annotated_header.formats:
            annotated_header.add_line(
                f'##FORMAT=<ID={field_id},Number=1,Type=Integer,Description="{description}">'
            )

    if include_pon_info:
        for field_id, field_type, description in PON_FORMAT_FIELD_DEFINITIONS:
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


def _variant_from_vcf_record(record: Any) -> Variant | None:
    """Build a Variant from a pysam VCF record when supported, else None."""
    if len(record.alts or ()) != 1:
        return None

    alt = record.alts[0]
    try:
        return Variant.from_vcf_fields(
            contig=record.contig,
            pos1=record.pos,
            ref=record.ref,
            alt=alt,
        )
    except ValueError:
        return None


def _annotate_case_format_fields(record: Any, evidence: AggregatedEvidence) -> None:
    """Set case-count FORMAT annotations on all sample columns."""
    for sample in record.samples.values():
        sample["SKUA_ALT_FWD"] = evidence.alt_forward
        sample["SKUA_ALT_REV"] = evidence.alt_reverse
        sample["SKUA_NON_ALT_FWD"] = evidence.non_alt_forward
        sample["SKUA_NON_ALT_REV"] = evidence.non_alt_reverse
        sample["SKUA_USABLE"] = evidence.usable
        sample["SKUA_UNUSABLE"] = evidence.unusable


def _annotate_pon_sample_format_fields(record: Any, *, artifact_posterior: float, bayes_factor: float) -> None:
    """Set PON model output FORMAT annotations on all sample columns."""
    for sample in record.samples.values():
        sample["SKUA_ARTIFACT_POSTERIOR"] = float(artifact_posterior)
        sample["SKUA_BAYES_FACTOR"] = float(bayes_factor)


def _render_annotated_vcf_payload(output_path: Path) -> str:
    """Read and return VCF text written to disk."""
    if output_path.suffix == ".gz":
        with gzip.open(output_path, "rt", encoding="utf-8") as handle:
            return handle.read()
    return output_path.read_text(encoding="utf-8")


def _vcf_write_mode(output_path: Path) -> str:
    """Return the pysam VariantFile write mode for VCF output path."""
    if output_path.suffix == ".gz":
        return "wz"
    return "w"


def verify_snv_vcf_to_annotated_vcf(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    output_path: str | Path | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> str:
    """Annotate an input VCF with case-count FORMAT fields for variants."""
    destination_path: Path
    created_temp = False
    if output_path is None:
        tmp = tempfile.NamedTemporaryFile(prefix="skua_", suffix=".vcf", delete=False)
        tmp.close()
        destination_path = Path(tmp.name)
        created_temp = True
    else:
        destination_path = Path(output_path)

    try:
        with pysam.VariantFile(str(vcf_path)) as source_vcf:
            header = _ensure_skua_vcf_header_fields(source_vcf.header, include_pon_info=False)
            with pysam.VariantFile(
                str(destination_path),
                _vcf_write_mode(destination_path),
                header=header,
            ) as out_vcf:
                for record in source_vcf:
                    variant = _variant_from_vcf_record(record)
                    if variant is not None:
                        evidence = verify_snv_variant(
                            alignment_file,
                            variant,
                            min_baseq=min_baseq,
                            min_mapq=min_mapq,
                        )
                        _annotate_case_format_fields(record, evidence)
                    out_vcf.write(record)

        return _render_annotated_vcf_payload(destination_path)
    finally:
        if created_temp and destination_path.exists():
            destination_path.unlink()


def verify_snv_vcf_to_annotated_vcf_with_normals(
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
    """Annotate an input VCF with case FORMAT and PON INFO fields for variants."""
    if normal_alignments is None:
        normal_alignments = []

    destination_path: Path
    created_temp = False
    if output_path is None:
        tmp = tempfile.NamedTemporaryFile(prefix="skua_", suffix=".vcf", delete=False)
        tmp.close()
        destination_path = Path(tmp.name)
        created_temp = True
    else:
        destination_path = Path(output_path)

    try:
        with pysam.VariantFile(str(vcf_path)) as source_vcf:
            header = _ensure_skua_vcf_header_fields(source_vcf.header, include_pon_info=True)
            with pysam.VariantFile(
                str(destination_path),
                _vcf_write_mode(destination_path),
                header=header,
            ) as out_vcf:
                for record in source_vcf:
                    variant = _variant_from_vcf_record(record)
                    if variant is not None:
                        pon_result = verify_snv_variant_with_normals(
                            alignment_file,
                            variant,
                            normal_alignments=normal_alignments,
                            min_baseq=min_baseq,
                            min_mapq=min_mapq,
                        )
                        case_evidence = pon_result["case_evidence"]
                        normal_samples_included = truncated_normal_evidences(
                            pon_result["normal_evidences"],
                            truncate=truncate,
                        )
                        normal_output_evidence = aggregate_evidence(normal_samples_included)
                        stats = compute_stats(
                            case_evidence,
                            normal_output_evidence,
                            per_sample_evidences=pon_result["normal_evidences"],
                            truncate=truncate,
                            pseudocount=pseudocount,
                            prior_variant_probability=prior_variant_probability,
                        )

                        _annotate_case_format_fields(record, case_evidence)
                        _annotate_pon_sample_format_fields(
                            record,
                            artifact_posterior=stats.artifact_posterior,
                            bayes_factor=stats.bayes_factor,
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

        return _render_annotated_vcf_payload(destination_path)
    finally:
        if created_temp and destination_path.exists():
            destination_path.unlink()


def verify_snv_variant(
    alignment_file: Any,
    variant: Variant,
    *,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> AggregatedEvidence:
    """Collect strand-aware evidence for one variant from one alignment."""
    return collect_snv_evidence_from_alignment(
        alignment_file,
        contig=variant.contig,
        ref_pos0=variant.ref_pos0,
        ref_base=variant.ref,
        alt_base=variant.alt,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
    )


def verify_snv_variant_with_normals(
    alignment_file: Any,
    variant: Variant,
    *,
    normal_alignments: list[Any] | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> dict[str, Any]:
    """Collect case and normal evidence for one variant."""
    if normal_alignments is None:
        normal_alignments = []

    case_evidence = verify_snv_variant(
        alignment_file,
        variant,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
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
        normal_evidence = verify_snv_variant(
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

    return {
        "case_evidence": case_evidence,
        "normal_evidences": normal_evidences,
        "normal_aggregate_evidence": normal_aggregate_evidence,
    }


def verify_snv_variants_from_vcf(
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
            verify_snv_variant(
                alignment_file,
                variant,
                min_baseq=min_baseq,
                min_mapq=min_mapq,
            ),
        )


def format_verification_results(
    results: Iterable[tuple[Variant, AggregatedEvidence]],
) -> list[dict[str, Any]]:
    """Convert verification results to JSON/tabular-ready row dictionaries."""
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


def render_verification_results_json(rows: Iterable[dict[str, Any]]) -> str:
    """Render formatted verification rows as JSON text."""
    return json.dumps(list(rows), indent=2)


def write_verification_results_json(
    rows: Iterable[dict[str, Any]],
    output_path: str | Path,
) -> None:
    """Write formatted verification rows to a JSON file."""
    Path(output_path).write_text(
        render_verification_results_json(rows),
        encoding="utf-8",
    )


def _build_verification_rows(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    min_baseq: int,
    min_mapq: int,
) -> list[dict[str, Any]]:
    """Build formatted verification rows from one alignment and one VCF."""
    return format_verification_results(
        verify_snv_variants_from_vcf(
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


def verify_snv_vcf_to_json(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    output_path: str | Path | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> str:
    """Run SNV verification from VCF and return JSON output, optionally writing to file."""
    rows = _build_verification_rows(
        alignment_file,
        vcf_path,
        min_baseq=min_baseq,
        min_mapq=min_mapq,
    )
    return _render_and_optionally_write(
        rows,
        renderer=render_verification_results_json,
        output_path=output_path,
    )


def verify_snv_variants_from_vcf_with_normals(
    alignment_file: Any,
    vcf_path: str | Path,
    *,
    normal_alignments: list[Any] | None = None,
    min_baseq: int = 20,
    min_mapq: int = 20,
) -> Iterator[tuple[Variant, dict[str, Any]]]:
    """Yield per-variant case+normal evidence for SNV records from a VCF file."""
    if normal_alignments is None:
        normal_alignments = []

    for variant in read_vcf_variant_file(vcf_path):
        yield (
            variant,
            verify_snv_variant_with_normals(
                alignment_file,
                variant,
                normal_alignments=normal_alignments,
                min_baseq=min_baseq,
                min_mapq=min_mapq,
            ),
        )


def format_verification_results_with_normals(
    results: Iterable[tuple[Variant, dict[str, Any]]],
    *,
    truncate: float = DEFAULT_TRUNCATE,
    pseudocount: float = sys.float_info.epsilon,
    prior_variant_probability: float = 0.5,
) -> list[dict[str, Any]]:
    """Convert PON verification results to JSON/tabular-ready row dictionaries."""
    rows: list[dict[str, Any]] = []
    for variant, pon_result in results:
        evidence = pon_result["case_evidence"]
        per_sample_evidences = pon_result["normal_evidences"]

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
                    "bayes_factor": stats.bayes_factor,
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


def _build_verification_rows_with_normals(
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
    """Build formatted PON verification rows from case + normal alignments and one VCF."""
    return format_verification_results_with_normals(
        verify_snv_variants_from_vcf_with_normals(
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


def verify_snv_vcf_to_json_with_normals(
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
    """Run PON SNV verification from VCF and return JSON output, optionally writing to file."""
    if normal_alignments is None:
        normal_alignments = []

    rows = _build_verification_rows_with_normals(
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
        renderer=render_verification_results_json,
        output_path=output_path,
    )


