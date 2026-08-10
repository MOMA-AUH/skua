"""Variant parsing and normalization helpers."""

import gzip
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Iterator


class VariantKind(str, Enum):
    """Supported simple VCF allele classes."""

    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"


def _normalize_simple_alleles(ref: str, alt: str) -> tuple[str, str]:
    """Return canonical supported alleles or raise ValueError."""
    if not isinstance(ref, str) or not isinstance(alt, str) or not ref or not alt:
        raise ValueError("REF and ALT must be non-empty strings")

    normalized_ref = ref.upper()
    normalized_alt = alt.upper()
    if any(
        base not in {"A", "C", "G", "T"}
        for base in normalized_ref + normalized_alt
    ):
        raise ValueError("REF and ALT must contain only A, C, G, or T")
    if normalized_ref == normalized_alt:
        raise ValueError("REF and ALT must be different")

    is_substitution = len(normalized_ref) == len(normalized_alt)
    is_simple_indel = len(normalized_ref) == 1 or len(normalized_alt) == 1
    if not is_substitution and not is_simple_indel:
        raise ValueError("Only simple substitutions and simple indels are supported")
    if not is_substitution and normalized_ref[0] != normalized_alt[0]:
        raise ValueError("Only left-anchored simple indels are supported")

    return normalized_ref, normalized_alt


@dataclass(frozen=True)
class Variant:
    """Validated simple-variant model using a 0-based reference position."""

    contig: str
    ref_pos0: int
    ref: str
    alt: str

    def __post_init__(self) -> None:
        """Canonicalize and validate the supported variant representation."""
        if not isinstance(self.contig, str) or not self.contig:
            raise ValueError("contig must be a non-empty string")
        if (
            not isinstance(self.ref_pos0, int)
            or isinstance(self.ref_pos0, bool)
            or self.ref_pos0 < 0
        ):
            raise ValueError("ref_pos0 must be a non-negative integer")
        ref, alt = _normalize_simple_alleles(self.ref, self.alt)
        object.__setattr__(self, "ref", ref)
        object.__setattr__(self, "alt", alt)

    @property
    def kind(self) -> VariantKind:
        """Return the simple allele class for this variant."""
        if len(self.ref) == len(self.alt):
            return VariantKind.SUBSTITUTION
        if len(self.ref) == 1 and len(self.alt) > 1:
            return VariantKind.INSERTION
        if len(self.ref) > 1 and len(self.alt) == 1:
            return VariantKind.DELETION
        raise ValueError("Only simple substitutions and simple indels are supported")

    @classmethod
    def from_vcf_fields(cls, *, contig: str, pos1: int, ref: str, alt: str) -> "Variant":
        """Build a Variant from basic VCF fields."""
        if not isinstance(pos1, int) or isinstance(pos1, bool) or pos1 < 1:
            raise ValueError("VCF POS must be >= 1")
        return cls(contig=contig, ref_pos0=pos1 - 1, ref=ref, alt=alt)


def parse_vcf_variant_line(line: str) -> Variant | None:
    """Parse one VCF line and return a Variant when applicable."""
    line = line.strip()
    if not line or line.startswith("#"):
        return None

    fields = line.split("\t")
    if len(fields) < 5:
        return None

    contig, pos_str, _id, ref, alt = fields[:5]
    if "," in alt:
        return None

    try:
        pos1 = int(pos_str)
    except ValueError:
        return None

    try:
        return Variant.from_vcf_fields(contig=contig, pos1=pos1, ref=ref, alt=alt)
    except ValueError:
        return None


def read_vcf_variant_file(path: str | Path) -> Iterator[Variant]:
    """Yield variants from a VCF file, skipping unsupported records."""
    path_obj = Path(path)
    if path_obj.suffix == ".gz":
        handle_cm = gzip.open(path_obj, "rt", encoding="utf-8")
    else:
        handle_cm = path_obj.open("r", encoding="utf-8")

    with handle_cm as handle:
        for line in handle:
            variant = parse_vcf_variant_line(line)
            if variant is not None:
                yield variant
