"""Variant parsing and normalization helpers."""

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


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

    try:
        pos1 = int(pos_str)
    except ValueError:
        return None

    try:
        return Variant.from_vcf_fields(contig=contig, pos1=pos1, ref=ref, alt=alt)
    except ValueError:
        return None


def read_vcf_snv_file(path: str | Path) -> Iterator[Variant]:
    """Yield SNV variants from a VCF file, skipping unsupported records."""
    path_obj = Path(path)
    if path_obj.suffix == ".gz":
        handle_cm = gzip.open(path_obj, "rt", encoding="utf-8")
    else:
        handle_cm = path_obj.open("r", encoding="utf-8")

    with handle_cm as handle:
        for line in handle:
            variant = parse_vcf_snv_line(line)
            if variant is not None:
                yield variant
