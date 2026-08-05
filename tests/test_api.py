"""Tests for the supported top-level Python API."""

import skua


def test_top_level_api_exports_documented_symbols() -> None:
    expected = {
        "AggregatedEvidence",
        "AnnotationStatus",
        "PonAnnotation",
        "PonArtifactMetadata",
        "Stats",
        "Variant",
        "VariantKind",
        "annotate_variant",
        "annotate_variant_with_normals",
        "annotate_variants",
        "annotate_variants_with_normals",
        "annotate_variants_from_vcf",
        "annotate_variants_from_pon",
        "annotate_vcf",
        "annotate_vcf_with_pon",
        "annotate_vcf_with_normals",
        "build_pon",
        "compute_stats",
    }

    assert expected <= set(skua.__all__)
    assert all(hasattr(skua, name) for name in expected)
