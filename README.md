# skua

[![Conda Version](https://img.shields.io/conda/vn/MOMA-AUH/skua?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/skua) [![Conda Downloads](https://img.shields.io/conda/dn/MOMA-AUH/skua?style=for-the-badge&cacheSeconds=300)](https://anaconda.org/MOMA-AUH/skua)

Skua annotates candidate somatic variants with read evidence and an artifact
probability using the [shearwater](https://doi.org/10.1093/bioinformatics/btt750)
statistical model. It compares a case sample against a panel of normals (PON),
which can be collected once and reused across cases.

The supported workflow is paired-end targeted DNA mapped with `bwa mem`.
Skua handles biallelic SNVs, MNVs, and simple left-anchored insertions and
deletions. It annotates supplied candidates and preserves their VCF FILTER
values; variant discovery and filtering decisions belong to your workflow.

## Installation

```bash
conda install MOMA-AUH::skua
```

Requires Python ≥ 3.11 and pysam ≥ 0.22, installed automatically by Conda.

## Quick start

Prepare a coordinate-sorted VCF of unique target alleles (`targets.vcf.gz`)
and a text file with one normal BAM or CRAM path per line (`normals.lst`).
Build and validate the panel, then annotate a case at those targets:

```bash
skua pon build \
  --vcf targets.vcf.gz \
  --normal-list normals.lst \
  --reference reference.fa \
  --output targets.pon.bcf

skua pon validate targets.pon.bcf --reference reference.fa

skua annotate \
  --pon targets.pon.bcf \
  --alignment case.bam \
  --reference reference.fa \
  --output calls.vcf.gz
```

Keep the panel's `.bcf.csi` index alongside the BCF. Subsequent cases need only
the panel, reference, and case alignment; the normal alignments are not reread.
Use `skua pon inspect targets.pon.bcf` to view the panel's settings and samples.

To annotate a case-specific subset, add `--vcf case-candidates.vcf.gz` to the
annotation command. Each supported allele must match a panel target exactly by
contig, position, REF, and ALT, so normalize both inputs consistently.

For a run without a saved panel, supply the normal alignments directly:

```bash
skua annotate \
  --vcf candidates.vcf.gz \
  --alignment case.bam \
  --normal-list normals.lst \
  --reference reference.fa \
  --output calls-direct.vcf.gz
```

### Input requirements

- **Alignments:** coordinate-sorted, indexed BAM or CRAM files with read-group
  IDs and sample names (`SM`). Only reads assigned to the selected sample's
  read groups contribute evidence. Use `--sample NAME` if case selection is
  ambiguous; a site-only VCF receives the selected case sample.
- **Normals:** one biological sample per file, with distinct sample names.
  Merge libraries from the same sample upstream. The case must not be in the
  panel.
- **Reference:** use the same assembly and contig names throughout. An indexed
  FASTA is required for CRAM and enables matching of equivalent indel placements.
  A panel built with a FASTA requires it during annotation. Without a FASTA,
  BAM workflows use exact-anchor indel matching.

## Reading the output

Skua preserves input records and adds counts, model scores, and assessment
fields. Case values are in FORMAT; panel summaries are in INFO.

| Field | Location | Meaning |
| --- | --- | --- |
| `SKUA_STATUS` | INFO | `ANNOTATED`, or an `UNSUPPORTED_*` reason for an allele that cannot be annotated. |
| `SKUA_ASSESSMENT` | FORMAT | `ASSESSED` or `INSUFFICIENT_EVIDENCE` for the selected case. |
| `SKUA_REASONS` | FORMAT | Unmet evidence requirements; `.` when assessed. |
| `SKUA_ARTIFACT_POSTERIOR` | FORMAT | Artifact probability from 0 to 1; lower values favor a variant. |
| `SKUA_LBF` | FORMAT | Natural log Bayes factor: positive favors artifact, negative favors variant. |
| `SKUA_ARTIFACT_PRIOR` | INFO | Effective artifact probability before evaluating the evidence. |
| `SKUA_PON_RHO` | INFO | Estimated beta-binomial dispersion parameter. |

Strand counts use `SKUA_ALT_FWD`, `SKUA_ALT_REV`, `SKUA_NON_ALT_FWD`, and
`SKUA_NON_ALT_REV`; totals use `SKUA_USABLE` and `SKUA_UNUSABLE`. Panel counts
use the corresponding `SKUA_PON_*` names, and `SKUA_PON_SAMPLE_COUNT` reports
the number of retained normals. Overlapping mates count once per fragment.

**Interpret scores only when `SKUA_ASSESSMENT=ASSESSED`.** When evidence is
insufficient, both scores are missing (`.`); counts and reasons remain
available. `SKUA_STATUS=ANNOTATED` alone does not mean the evidence is sufficient.
Unsupported alleles remain in the output without scores; `--strict` rejects
them instead.

By default, assessment requires at least one usable case observation and one
usable pooled normal observation. These defaults only exclude absent evidence.
Set coverage requirements and posterior cutoffs appropriate for your assay;
`ASSESSED` means the configured requirements were met.

## Common options

| Option | Default | Purpose |
| --- | --- | --- |
| `--min-baseq`, `--min-mapq` | `20` | Minimum base and mapping qualities for case and normal evidence. |
| `--truncate` | `0.1` | Exclude normals whose ALT fraction is at or above this threshold at a site. |
| `--min-case-depth`, `--min-normal-depth` | `1` | Minimum usable case and pooled normal depths. |
| `--min-normal-samples` | `0` | Minimum number of retained normals; `0` disables this requirement. |
| `--min-case-strand-depth`, `--min-normal-strand-depth` | `0` | Minimum usable depth on each strand; `0` disables this requirement. |
| `--prior-artifact-probability` | `0.5` | Fallback prior when the input has no `SKUA_ARTIFACT_PRIOR`. |
| `--force` | Off | Replace existing output and recompute existing Skua annotations. |

With `--pon`, quality thresholds come from the panel and cannot be overridden.
Normal depth and sample requirements apply after truncation. Record-specific
priors use `INFO/SKUA_ARTIFACT_PRIOR` (`Number=A`, `Type=Float`) and must be
strictly between 0 and 1. Priors should be independent of the case and panel
evidence being scored.

Omit `--output` to write VCF to standard output. Compressed VCF outputs are not
indexed automatically; create an index after annotation if needed. For all
options, run `skua annotate --help` or `skua pon build --help`.

## Python API

Import supported functions directly from `skua`. For example, annotate a case
with an existing panel:

```python
import pysam
from skua import annotate_vcf_with_pon

with pysam.AlignmentFile("case.bam", "rb") as case:
    annotate_vcf_with_pon(
        case,
        "targets.pon.bcf",
        reference_path="reference.fa",
        output_path="calls-python.vcf.gz",
    )
```

The library also provides per-variant and batch evidence collection, model
assessment, PON construction, and JSON output. See the
[Python API reference](docs/v1-contract.md#python-api) for functions and sample
selection rules.

## Further reading

- [Input, output, and API reference](docs/v1-contract.md)
- [Indel matching examples and limits](docs/indel-representations.md)

## License

MIT. See [LICENSE](LICENSE).
