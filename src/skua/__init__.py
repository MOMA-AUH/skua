"""Public Python API for annotating simple variants with read evidence."""

from .core import (
    AnnotationStatus,
    AnnotationSummary,
    PonAnnotation,
    annotate_variant,
    annotate_variant_with_normals,
    annotate_variants_from_vcf,
    annotate_variants_from_vcf_with_normals,
    annotate_vcf,
    annotate_vcf_to_json,
    annotate_vcf_to_json_with_normals,
    annotate_vcf_with_normals,
)
from .evidence import AggregatedEvidence, AlleleSupport, ReadAlleleCall, UnusableReason
from .stats import Stats, compute_stats
from .variants import Variant, VariantKind
from ._version import __version__

__all__ = [
    "AggregatedEvidence",
    "AnnotationStatus",
    "AnnotationSummary",
    "AlleleSupport",
    "PonAnnotation",
    "ReadAlleleCall",
    "Stats",
    "UnusableReason",
    "Variant",
    "VariantKind",
    "annotate_variant",
    "annotate_variant_with_normals",
    "annotate_variants_from_vcf",
    "annotate_variants_from_vcf_with_normals",
    "annotate_vcf",
    "annotate_vcf_to_json",
    "annotate_vcf_to_json_with_normals",
    "annotate_vcf_with_normals",
    "compute_stats",
    "__version__",
]
