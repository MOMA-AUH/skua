"""Variant parsing and normalization helpers."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Variant:
    """Minimal SNV variant model using 0-based reference position."""

    contig: str
    ref_pos0: int
    ref: str
    alt: str

    @classmethod
    def from_vcf_fields(cls, *, contig: str, pos1: int, ref: str, alt: str) -> "Variant":
        """Build a Variant from basic VCF fields."""
        if pos1 < 1:
            raise ValueError("VCF POS must be >= 1")
        if len(ref) != 1 or len(alt) != 1:
            raise ValueError("SNV variants require single-base REF and ALT")

        return cls(contig=contig, ref_pos0=pos1 - 1, ref=ref, alt=alt)


def parse_vcf_snv_line(line: str) -> Variant | None:
    """Parse one VCF line and return a SNV Variant when applicable."""
    line = line.strip()
    if not line or line.startswith("#"):
        return None

    fields = line.split("\t")
    if len(fields) < 5:
        return None

    contig, pos_str, _id, ref, alt = fields[:5]
    if "," in alt:
        return None
    if len(ref) != 1 or len(alt) != 1:
        return None

    try:
        pos1 = int(pos_str)
    except ValueError:
        return None
    if pos1 < 1:
        return None

    return Variant(contig=contig, ref_pos0=pos1 - 1, ref=ref, alt=alt)
