"""Reference identity and compatibility for the contigs used by a target set."""

from dataclasses import asdict, dataclass
import hashlib
import re
from typing import Any, Iterable
from ._headers import remove_header_records


@dataclass(frozen=True)
class ContigReference:
    name: str
    length: int | None
    md5: str | None = None
    verified: bool = False


@dataclass(frozen=True)
class ReferenceIdentity:
    contigs: tuple[ContigReference, ...]

    @property
    def status(self) -> str:
        return "VERIFIED" if self.contigs and all(c.verified for c in self.contigs) else "INSUFFICIENT_METADATA"

    def as_dict(self) -> dict[str, Any]:
        return {"status": self.status, "contigs": [asdict(c) for c in self.contigs]}


def write_reference_header(header: Any, identity: ReferenceIdentity) -> None:
    """Replace only reference metadata owned by Skua."""
    remove_header_records(header, ("SKUA_REFERENCE", "SKUA_REFERENCE_STATUS"))
    header.add_meta("SKUA_REFERENCE_STATUS", value=identity.status)
    for contig in identity.contigs:
        items = [("ID", contig.name), ("Verified", "1" if contig.verified else "0")]
        if contig.length is not None:
            items.append(("Length", str(contig.length)))
        if contig.md5 is not None:
            items.append(("MD5", contig.md5))
        header.add_meta("SKUA_REFERENCE", items=items)


def read_reference_header(header: Any) -> ReferenceIdentity:
    """Read required reference metadata without upgrading incomplete identities."""
    statuses = [record.value for record in header.records if record.key == "SKUA_REFERENCE_STATUS"]
    if len(statuses) != 1:
        raise ValueError("PON artifact must contain exactly one SKUA_REFERENCE_STATUS")
    contigs = []
    seen = set()
    for record in header.records:
        if record.key != "SKUA_REFERENCE":
            continue
        values = {key: str(value).strip('"') for key, value in record.items()}
        name = values.get("ID")
        if not name or name in seen:
            raise ValueError("PON artifact has missing or duplicate reference contig identity")
        seen.add(name)
        length = None
        if "Length" in values:
            try:
                length = int(values["Length"])
            except ValueError as exc:
                raise ValueError(f"Invalid reference length for {name!r}") from exc
            if length <= 0:
                raise ValueError(f"Invalid reference length for {name!r}")
        md5 = values.get("MD5")
        if md5 is not None:
            md5 = md5.lower()
            if re.fullmatch(r"[0-9a-f]{32}", md5) is None:
                raise ValueError(f"Invalid reference checksum for {name!r}")
        verified = values.get("Verified")
        if verified not in {"0", "1"} or (verified == "1" and (length is None or md5 is None)):
            raise ValueError(f"Invalid reference verification metadata for {name!r}")
        contigs.append(ContigReference(name, length, md5, verified == "1"))
    identity = ReferenceIdentity(tuple(contigs))
    if statuses[0] != identity.status:
        raise ValueError("PON reference status is inconsistent with its contig identities")
    return identity


def check_reference_compatibility(
    contigs: Iterable[str], *, alignment_files: list[tuple[str, Any]],
    fasta_file: Any | None = None,
    pon_reference: ReferenceIdentity | None = None,
    vcf_header: Any | None = None,
) -> ReferenceIdentity:
    """Reject conflicts in the available reference dictionaries."""
    names = sorted(set(contigs))
    dictionaries = []
    for label, alignment in alignment_files:
        header = getattr(alignment, "header", None)
        values = header if isinstance(header, dict) else header.to_dict() if header is not None else {}
        dictionaries.append((label, {item["SN"]: item for item in values.get("SQ", [])}, True))
    panel_contigs = {} if pon_reference is None else {c.name: c for c in pon_reference.contigs}
    if pon_reference is not None:
        for name in names:
            if name not in panel_contigs:
                raise ValueError(f"PON artifact is missing reference identity for {name!r}")
        dictionaries.append(("PON artifact", {
            c.name: {"LN": c.length, "M5": c.md5} for c in pon_reference.contigs
        }, True))
    if vcf_header is not None:
        dictionaries.append(("Target VCF", {
            record.get("ID"): {"LN": record.get("length"), "M5": record.get("md5")}
            for record in vcf_header.records if record.key == "contig"
        }, False))
    if fasta_file is not None:
        dictionary = {}
        for name in names:
            if name not in fasta_file.references:
                raise ValueError(f"Reference FASTA does not contain contig {name!r}")
            length = fasta_file.get_reference_length(name)
            digest = hashlib.md5(usedforsecurity=False)
            for start in range(0, length, 1024 * 1024):
                bases = fasta_file.fetch(name, start, min(start + 1024 * 1024, length))
                digest.update(bases.upper().encode("ascii").translate(
                    None, bytes(range(33)) + bytes(range(127, 256)),
                ))
            dictionary[name] = {"LN": length, "M5": digest.hexdigest()}
        dictionaries.append(("Reference FASTA", dictionary, True))
    result = []
    for name in names:
        length = None
        md5 = None
        verified = any(required for _label, _dictionary, required in dictionaries)
        if pon_reference is not None:
            verified = verified and panel_contigs[name].verified
        for label, dictionary, required in dictionaries:
            observed = dictionary.get(name, {}).get("LN")
            if observed is not None:
                if isinstance(observed, bool) or re.fullmatch(r"[0-9]+", str(observed)) is None or int(observed) <= 0:
                    raise ValueError(f"Invalid reference length for {name!r} in {label}")
                observed = int(observed)
                if length is not None and observed != length:
                    raise ValueError(f"Conflicting reference length for {name!r} in {label}")
                length = observed
            checksum = dictionary.get(name, {}).get("M5")
            if required:
                verified = verified and observed is not None and checksum is not None
            if checksum is not None:
                checksum = str(checksum).lower()
                if re.fullmatch(r"[0-9a-f]{32}", checksum) is None:
                    raise ValueError(f"Invalid reference checksum for {name!r} in {label}")
                if md5 is not None and checksum != md5:
                    raise ValueError(f"Conflicting reference checksum for {name!r} in {label}")
                md5 = checksum
        result.append(ContigReference(name, length, md5, verified))
    return ReferenceIdentity(tuple(result))
