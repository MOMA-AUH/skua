"""Public Python API for annotating simple variants with read evidence."""

from .core import (
    AnnotationStatus,
    PonAnnotation,
    annotate_variants_from_pon,
    annotate_variant,
    annotate_variant_with_normals,
    annotate_variants,
    annotate_variants_with_normals,
    annotate_variants_from_vcf,
    annotate_variants_from_vcf_with_normals,
    annotate_vcf,
    annotate_vcf_to_json,
    annotate_vcf_to_json_with_normals,
    annotate_vcf_with_pon,
    annotate_vcf_with_normals,
    build_pon,
)
from .evidence import AggregatedEvidence, AlleleSupport, ReadAlleleCall, UnusableReason
from .pon import PonArtifactMetadata, read_pon_evidence, read_pon_metadata
from .stats import Stats, compute_stats
from .variants import Variant, VariantKind
from ._version import __version__

__all__ = [
    "AggregatedEvidence",
    "AnnotationStatus",
    "AlleleSupport",
    "PonAnnotation",
    "PonArtifactMetadata",
    "ReadAlleleCall",
    "Stats",
    "UnusableReason",
    "Variant",
    "VariantKind",
    "annotate_variant",
    "annotate_variant_with_normals",
    "annotate_variants",
    "annotate_variants_with_normals",
    "annotate_variants_from_vcf",
    "annotate_variants_from_vcf_with_normals",
    "annotate_variants_from_pon",
    "annotate_vcf",
    "annotate_vcf_to_json",
    "annotate_vcf_to_json_with_normals",
    "annotate_vcf_with_pon",
    "annotate_vcf_with_normals",
    "build_pon",
    "compute_stats",
    "read_pon_evidence",
    "read_pon_metadata",
    "__version__",
]
