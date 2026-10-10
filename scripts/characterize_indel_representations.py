"""Reproduce a small BWA-MEM indel-representation experiment (not an assay benchmark)."""

import argparse
from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path
import random
import subprocess

import pysam

from skua import annotate_variant, annotate_variants
from skua.evidence import is_accepted_sam_flag
from skua.variants import Variant


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bwa", default="bwa")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    reference = args.output_dir / "reference.fa"
    fastq1, fastq2 = (args.output_dir / name for name in ("reads.1.fq", "reads.2.fq"))
    rng = random.Random(19)
    variants = []
    reverse_complement = str.maketrans("ACGT", "TGCA")
    with reference.open("w") as fasta, fastq1.open("w") as mate1, fastq2.open("w") as mate2:
        for motif_name, motif in (("homopolymer", "A"), ("tandem", "AT")):
            for kind in ("insertion", "deletion"):
                contig = f"{motif_name}_{kind}"
                sequence = "".join(rng.choices("ACGT", k=499)) + "C" + motif * 6 + "G"
                sequence += "".join(rng.choices("ACGT", k=1000 - len(sequence)))
                variant = Variant(contig, 499, "C" if kind == "insertion" else "C" + motif,
                                  "C" + motif if kind == "insertion" else "C")
                variants.append(variant)
                haplotype = sequence[:499] + variant.alt + sequence[499 + len(variant.ref):]
                fasta.write(f">{contig}\n{sequence}\n")
                for index in range(40):
                    start = 375 + index % 20
                    fragment = haplotype[start:start + 220]
                    for handle, bases in ((mate1, fragment[:150]),
                                          (mate2, fragment[-150:].translate(reverse_complement)[::-1])):
                        handle.write(f"@{contig}_{index}\n{bases}\n+\n{'I' * len(bases)}\n")
    pysam.faidx(str(reference))
    version = subprocess.run([args.bwa], capture_output=True, text=True, check=False).stderr
    command = [args.bwa, "mem", "-t", "1", "-I", "220,10,250,190",
               "-R", "@RG\\tID:synthetic\\tSM:synthetic", str(reference), str(fastq1), str(fastq2)]
    with (args.output_dir / "bwa.log").open("w") as log:
        subprocess.run([args.bwa, "index", str(reference)], stderr=log, check=True)
        with (args.output_dir / "reads.sam").open("w") as sam:
            subprocess.run(command, stdout=sam, stderr=log, check=True)
    bam = args.output_dir / "reads.bam"
    pysam.sort("-o", str(bam), str(args.output_dir / "reads.sam"))
    pysam.index(str(bam))
    results = []
    with pysam.AlignmentFile(bam, "rb") as alignment:
        for variant in variants:
            offsets: Counter[int] = Counter()
            accepted = 0
            for read in alignment.fetch(variant.contig, variant.ref_pos0, variant.ref_pos0 + 1):
                if not is_accepted_sam_flag(read.flag) or read.mapping_quality < 20:
                    continue
                accepted += 1
                position = read.reference_start
                for op, length in read.cigartuples:
                    if op in (1, 2):
                        offsets[position - 1 - variant.ref_pos0] += 1
                    if op in (0, 2, 3, 7, 8):
                        position += length
            exact = annotate_variant(alignment, variant)
            equivalent = annotate_variant(alignment, variant, reference_path=reference)
            with pysam.FastaFile(reference) as fasta:
                shifted_ref = fasta.fetch(variant.contig, 502, 502 + len(variant.ref))
            shifted = Variant(
                variant.contig, 502, shifted_ref,
                shifted_ref[0] + ("A" if variant.contig.startswith("homopolymer") else "TA")
                if len(variant.alt) > 1 else shifted_ref[0],
            )
            shifted_exact = annotate_variant(alignment, shifted)
            shifted_equivalent = annotate_variant(alignment, shifted, reference_path=reference)
            # Duplicate targets force the batch path on the same observed allele.
            batch = list(annotate_variants(alignment, [variant, variant], reference_path=reference))
            assert batch == [(variant, equivalent), (variant, equivalent)]
            assert equivalent.alt_forward + equivalent.alt_reverse == 40
            assert shifted_equivalent.alt_forward + shifted_equivalent.alt_reverse == 40
            results.append({"variant": asdict(variant), "accepted_alignment_records": accepted,
                            "cigar_anchor_offsets": dict(sorted(offsets.items())),
                            "exact_anchor": asdict(exact), "reference": asdict(equivalent),
                            "shifted_target": asdict(shifted),
                            "shifted_target_exact_anchor": asdict(shifted_exact),
                            "shifted_target_reference": asdict(shifted_equivalent)})
    report = {"bwa_version": next(line for line in version.splitlines() if line.startswith("Version:")),
              "command": command, "seed": 19, "read_length": 150, "fragment_length": 220,
              "fragments_per_allele": 40, "results": results}
    payload = json.dumps(report, indent=2) + "\n"
    (args.output_dir / "results.json").write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
