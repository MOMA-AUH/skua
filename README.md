# skua

[![Conda Version](https://img.shields.io/conda/vn/MOMA-AUH/skua?cacheSeconds=300)](https://anaconda.org/MOMA-AUH/skua) [![Conda Downloads](https://img.shields.io/conda/dn/MOMA-AUH/skua?cacheSeconds=300)](https://anaconda.org/MOMA-AUH/skua)

Implementation of the [shearwater](https://doi.org/10.1093/bioinformatics/btt750) statistical model to assess somatic variant evidence in aligned reads, with support for substitutions, MNVs, and simple insertions and deletions. The **shearwater** authors named their algorithm after seabirds that fly long distances over the ocean, watching the water closely and eventually dive into the water to catch prey. Due to the heavy reuse of the algorithmic core, it is only natural to name this **skua** — a seabird that hunts and steals from other birds.

## Installation

The recommended way to install **skua** is via [conda](https://docs.conda.io/), using the `MOMA-AUH` channel:

```bash
conda install MOMA-AUH::skua
```

## Commands

### `annotate`

Annotate a VCF file with read counts, quality metrics, and artifact posteriors.

```bash
skua annotate \
  --vcf input.vcf.gz \
  --alignment case.bam \
  --normal-list normals.lst \
  --output output.vcf.gz
```

Key input parameters:
- `--vcf`: Input VCF file to annotate
- `--alignment`: Case BAM or CRAM file
- `--normal-list`: Text file with one normal BAM or CRAM path per line
- `--sample`: Case sample to annotate when VCF/BAM sample matching is ambiguous
- `--reference`: Reference FASTA file, required when any input alignment is CRAM
- `--output`: Optional output VCF path; if omitted, output is written to `stdout`

Skua resolves a single case sample from the VCF sample names and alignment read-group `SM` tags. Use `--sample` when that resolution is ambiguous; it must name a VCF sample and an alignment sample. For a site-only VCF, skua adds the selected alignment sample as the sole output sample. The selected sample must have at least one read-group `ID` in the alignment header. In all cases, only reads whose `RG` tag names one of those read groups contribute case evidence. Untagged reads and reads from unknown or unassigned read groups are excluded.

Other optional parameters:
- `--min-baseq` (default `20`): Minimum base quality for read bases
- `--min-mapq` (default `20`): Minimum mapping quality for reads
- `--truncate` (default `0.1`): Truncation percentile for PON sample inclusion
- `--pseudocount` (default `sys.float_info.epsilon`): Pseudocount for beta-binomial rate estimates
- `--prior-variant-probability` (default `0.5`): Prior probability for variant model
- `--strict`: Fail before writing output if any VCF record cannot be annotated

Alignment records must be mapped primary records from a proper pair. Records
whose mate is unmapped, or which are marked secondary, supplementary, failed
quality control, or duplicate, are excluded before evidence classification. If
both mates overlap a variant, only one record is counted for their shared read
group and query name. Agreeing usable mates count once on the first mate's
strand; one usable mate takes precedence over an unusable mate; conflicting
usable mates count once as unusable; and two unusable mates count once as
unusable.

Truncation controls how conservative the panel-of-normals aggregation is at each site. A normal sample is included only if its ALT fraction is strictly less than `--truncate`. With `--truncate 0.1`, normals with ALT fraction `< 0.1` are kept and normals with ALT fraction `>= 0.1` are excluded.

Output FORMAT fields:
- `SKUA_ALT_FWD`: Count of ALT-supporting reads on forward strand
- `SKUA_ALT_REV`: Count of ALT-supporting reads on reverse strand
- `SKUA_NON_ALT_FWD`: Count of non-ALT reads on forward strand
- `SKUA_NON_ALT_REV`: Count of non-ALT reads on reverse strand
- `SKUA_USABLE`: Total usable reads at this locus
- `SKUA_UNUSABLE`: Total unusable reads (low quality, INDELs at locus, etc.)
- `SKUA_ARTIFACT_POSTERIOR`: Posterior probability of artifact model (0–1)
- `SKUA_LOG_BAYES_FACTOR`: Log Bayes factor comparing artifact vs. variant models

Output INFO fields:
- `SKUA_STATUS`: Annotation outcome for every record. `ANNOTATED` records receive Skua evidence; unsupported records are retained with an `UNSUPPORTED_*` status and are not assigned new Skua evidence fields. `UNSUPPORTED_RECORD` includes records with no alternate allele (`ALT=.`).
- `SKUA_PON_SAMPLE_COUNT`: Number of normal samples included after truncation
- `SKUA_PON_ALT_FWD`, `SKUA_PON_ALT_REV`, `SKUA_PON_NON_ALT_FWD`, `SKUA_PON_NON_ALT_REV`: Aggregated read counts across normals
- `SKUA_PON_USABLE`, `SKUA_PON_UNUSABLE`: Aggregated usable/unusable counts
- `SKUA_PON_DISPERSION_FACTOR`: Beta-binomial dispersion parameter estimate

By default, unsupported records do not stop the run. Use `--strict` to reject any input containing one before an output file is created. VCF output is written to `--output` or standard output.

## Python API

The supported library API is available directly from `skua`. It accepts
substitutions, MNVs, and simple insertions and deletions.

```python
import pysam
from skua import Variant, annotate_variant, annotate_variant_with_normals

variant = Variant.from_vcf_fields(contig="chr1", pos1=106, ref="A", alt="T")

with pysam.AlignmentFile("case.bam", "rb") as case_bam:
    evidence = annotate_variant(case_bam, variant)
    print(evidence.alt_forward, evidence.alt_reverse)

    with pysam.AlignmentFile("normal.bam", "rb") as normal_bam:
        annotation = annotate_variant_with_normals(
            case_bam, variant, normal_alignments=[normal_bam]
        )
        print(annotation.case_evidence.usable)
```

For batch work, open each alignment once and use
`annotate_variants_from_vcf()`; this avoids repeatedly opening the same BAM or
CRAM.

## Requirements

- Python ≥ 3.11
- pysam ≥ 0.22

## License

MIT. See [LICENSE](LICENSE) for details.
