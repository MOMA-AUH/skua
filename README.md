# skua

[![Conda Version](https://img.shields.io/conda/vn/MOMA-AUH/skua?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/skua) [![Conda Downloads](https://img.shields.io/conda/dn/MOMA-AUH/skua?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/skua)

Implementation of the [shearwater](https://doi.org/10.1093/bioinformatics/btt750) statistical model to assess somatic variant evidence in aligned reads, with support for substitutions, MNVs, and left-anchored simple insertions and deletions. The **shearwater** authors named their algorithm after seabirds that fly long distances over the ocean, watching the water closely and eventually dive into the water to catch prey. Due to the heavy reuse of the algorithmic core, it is only natural to name this **skua** — a seabird that hunts and steals from other birds.

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
- `--vcf`: Input VCF file to annotate; required with `--normal-list` and optional with `--pon`
- `--alignment`: Case BAM or CRAM file
- `--normal-list`: Text file with one normal BAM or CRAM path per line
- `--pon`: Precomputed PON BCF; mutually exclusive with `--normal-list`
- `--sample`: Case sample to annotate when VCF/BAM sample matching is ambiguous
- `--reference`: Reference FASTA file, required when any input alignment is CRAM
- `--output`: Optional output VCF path; if omitted, output is written to `stdout`
- `--force`: Replace an existing output and recompute existing Skua annotations

Skua resolves a single case sample from the VCF sample names and alignment read-group `SM` tags. Use `--sample` when that resolution is ambiguous; it must name a VCF sample and an alignment sample. For a site-only VCF, skua adds the selected alignment sample as the sole output sample. The selected sample must have at least one read-group `ID` in the alignment header. In all cases, only reads whose `RG` tag names one of those read groups contribute case evidence. Untagged reads and reads from unknown or unassigned read groups are excluded.

Other optional parameters:
- `--min-baseq` (default `20`): Minimum base quality for read bases
- `--min-mapq` (default `20`): Minimum mapping quality for reads
- `--truncate` (default `0.1`): Truncation percentile for PON sample inclusion
- `--pseudocount` (default `sys.float_info.epsilon`): Pseudocount for beta-binomial rate estimates
- `--prior-artifact-probability` (default `0.5`): Fallback artifact prior when the input record has no `SKUA_ARTIFACT_PRIOR`
- `--strict`: Fail before writing output if any VCF record cannot be annotated

Supported biallelic records may provide an allele-specific artifact prior in
INFO:

```vcf
##INFO=<ID=SKUA_ARTIFACT_PRIOR,Number=A,Type=Float,Description="Prior probability that the ALT allele is an artifact before Skua evidence">
```

The value must be finite and strictly between `0` and `1`. An absent field or
`.` uses `--prior-artifact-probability`; an explicit record value takes
precedence. Skua writes the effective value to every annotated output record so
the artifact posterior can be reproduced. The artifact prior affects
`SKUA_ARTIFACT_POSTERIOR` but does not affect `SKUA_LOG_BAYES_FACTOR`.

Artifact priors must be calibrated and independent of the case-read and PON
evidence evaluated by Skua. Deriving a prior from the same evidence would count
that evidence twice. INFO is appropriate for an allele- or site-specific prior
shared by samples; a case- or sample-specific prior belongs in FORMAT instead.

Alignment records must be mapped primary records from a proper pair. Records
whose mate is unmapped, or which are marked secondary, supplementary, failed
quality control, or duplicate, are excluded before evidence classification. If
both mates overlap a variant, only one record is counted for their shared read
group and query name. Agreeing usable mates count once on the first mate's
strand; one usable mate takes precedence over an unusable mate; conflicting
usable mates count once as unusable; and two unusable mates count once as
unusable.

Indel support requires an explicit CIGAR insertion or deletion immediately
after the left anchor and an aligned base on the right. Soft clips, reference
skips, terminal or adjacent complex events, and reads without sequence or base
qualities are unusable rather than ALT or reference evidence.

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
- `SKUA_ARTIFACT_PRIOR`: Effective prior probability that the ALT allele is an artifact before Skua evidence
- `SKUA_STATUS`: Annotation outcome for every record. `ANNOTATED` records receive Skua evidence; unsupported records are retained with an `UNSUPPORTED_*` status and are not assigned new Skua evidence fields. `UNSUPPORTED_RECORD` includes records with no alternate allele (`ALT=.`).
- `SKUA_PON_SAMPLE_COUNT`: Number of normal samples included after truncation
- `SKUA_PON_ALT_FWD`, `SKUA_PON_ALT_REV`, `SKUA_PON_NON_ALT_FWD`, `SKUA_PON_NON_ALT_REV`: Aggregated read counts across normals
- `SKUA_PON_USABLE`, `SKUA_PON_UNUSABLE`: Aggregated usable/unusable counts
- `SKUA_PON_DISPERSION_FACTOR`: Beta-binomial dispersion parameter estimate

By default, unsupported records do not stop the run. Use `--strict` to reject any input containing one before an output file is created. VCF output is written to `--output` or standard output.

File outputs are transactional and no-clobber by default. Skua writes a sibling
temporary VCF and publishes it atomically only after annotation finishes;
`--force` permits replacement of an existing output. Without `--force`, input
VCFs that already define generated `SKUA_*` INFO or FORMAT annotations are
rejected. With `--force`, Skua removes those definitions and all associated
record and sample values before installing the canonical definitions and
computing fresh annotations. Unsupported records therefore retain only their
new `SKUA_STATUS` among generated Skua fields. `SKUA_ARTIFACT_PRIOR` is the
exception: it is a validated input field whose effective value is intentionally
preserved and written to the result. Non-Skua fields are preserved unchanged.

### `pon`

Build reusable panel-of-normals artifacts, inspect their provenance, and
validate them before annotation.

#### `build`

For a fixed target set, normal evidence can be collected once and reused across
case samples:

```bash
skua pon build \
  --vcf hotspots.vcf.gz \
  --normal-list normals.lst \
  --output hotspots.pon.bcf

skua annotate \
  --pon hotspots.pon.bcf \
  --alignment case.bam \
  --output calls.vcf.gz
```

An input VCF can also define a subset of case-specific targets while the PON
supplies their cached normal evidence:

```bash
skua annotate \
  --vcf case-candidates.vcf.gz \
  --pon hotspots.pon.bcf \
  --alignment case.bam \
  --output calls.vcf.gz
```

The indexed PON BCF contains the target records and six strand-aware evidence counts
for every normal sample. During annotation, skua reads those cached counts and
only accesses the case alignment. Per-sample counts are retained so that
`--truncate` and dispersion estimation are still evaluated at annotation time.
Construction also writes a companion `.bcf.csi` index; target records must be
coordinate-sorted. The BCF and CSI are built under sibling temporary names and
published only after both are complete. Existing artifacts are not replaced by
default; pass `--force` to `pon build` to replace an existing pair and remove
pre-existing Skua annotations from its target records before rebuilding the
PON. The target's validated `SKUA_ARTIFACT_PRIOR` is retained.

#### `inspect`

Inspect a PON header without scanning every target record:

```bash
skua pon inspect hotspots.pon.bcf
skua pon inspect hotspots.pon.bcf --json
```

`inspect` reports the artifact format, index presence, schema and evidence
policy versions, quality thresholds, producer version, and normal samples. It
is intentionally permissive: it can describe an incompatible artifact without
claiming that the installed Skua can use it.

#### `validate`

Validate that a PON is structurally sound and compatible with the installed
Skua version:

```bash
skua pon validate hotspots.pon.bcf
skua pon validate hotspots.pon.bcf --reference reference.fa --targets hotspots.vcf.gz
skua pon validate hotspots.pon.bcf --json
```

`validate` checks the BCF and CSI index, provenance metadata, required FORMAT
fields, target alleles, duplicate targets, coordinate order, and per-sample
evidence counts. `--reference` additionally checks every PON REF allele;
`--targets` requires the supplied VCF to define exactly the PON targets. It
returns exit status `0` when valid and `1` when any check fails.

#### Using a PON with `annotate`

When `--vcf` is omitted in cached mode, the PON records define the targets. When
`--vcf` is supplied, its records define the output and their existing IDs,
quality values, filters, INFO annotations, and samples are preserved. Every
supported input allele must have an exact `CHROM`, `POS`, `REF`, and `ALT` match
in the PON; extra PON targets are ignored, while a missing match fails before
the output is created. Alleles should therefore be represented and normalized
consistently when the input VCF and PON are produced.

The record-defining source also owns the artifact prior. In PON-only mode, a
`SKUA_ARTIFACT_PRIOR` retained from the target VCF takes precedence over the CLI
fallback. With `--vcf --pon`, the separate input VCF takes precedence: its
explicit value is used, while an absent or missing value uses the CLI fallback
rather than a prior stored in the reusable PON. Store target-intrinsic priors in
the target VCF before PON construction. A prior prepared specifically for one
case run should live in that run's separate `--vcf`, not in a reusable PON. If
samples within one VCF require different priors, those values belong in FORMAT
rather than this sample-shared INFO field.

Cached annotation always uses the `--min-baseq` and `--min-mapq` values stored
in the PON; those options cannot be supplied together with `--pon`. PON
construction rejects unsupported or multiallelic target records and requires a
unique `CHROM`, `POS`, `REF`, and `ALT` target allele, as well as a unique
read-group `SM` name in each normal alignment. A PON should be rebuilt when the
reference assembly, alignment/evidence policy, or quality thresholds change.

## Python API

The supported library API is available directly from `skua`. It accepts
substitutions, MNVs, and left-anchored simple insertions and deletions.

```python
import pysam
from skua import (
    Variant,
    annotate_variant,
    annotate_variant_with_normals,
    annotate_variants,
)

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
CRAM. Skua automatically shares an alignment fetch across nearby variants while
retaining direct per-site fetches for sparse records. For batches already
represented as `Variant` objects, use `annotate_variants()` directly:

```python
variants = [
    Variant.from_vcf_fields(contig="chr1", pos1=106, ref="A", alt="T"),
    Variant.from_vcf_fields(contig="chr1", pos1=109, ref="A", alt="C"),
]

with pysam.AlignmentFile("case.bam", "rb") as case_bam:
    for variant, evidence in annotate_variants(case_bam, variants):
        print(variant, evidence.usable)
```

Use `annotate_variants_with_normals()` for the corresponding case-plus-PON
workflow.

## Requirements

- Python ≥ 3.11
- pysam ≥ 0.22

## License

MIT. See [LICENSE](LICENSE) for details.
