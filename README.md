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

Each normal alignment represents **one biological sample**, identified by its
single distinct read-group `SM` value. Merge multiple libraries or files for the
same sample upstream, retaining their distinct read-group IDs. Skua rejects
repeated handles, repeated files (including local symlinks/hardlinks), and repeated
normal `SM` values instead of counting them as independent normals. Normal
read-group IDs must be unambiguous. As for the case, only reads whose `RG` belongs
to that sample contribute normal evidence; untagged, unknown, and SM-less groups
are excluded from both usable and unusable counts. These rules apply to live
normals and PON construction.

The selected case must not be a panel member: live annotation rejects a shared
input file or matching sample name, and cached annotation checks the PON's sample
names. Supply a panel excluding the case; Skua does not perform automatic
leave-one-out analysis. Identity checks rely on truthful `SM` metadata, not an
inference of biological identity from reads. The low-level Python evidence APIs
also reject known duplicate inputs and case membership, but permit alignment-like
objects without headers; for those objects sample identity and RG assignment
cannot be verified, and the caller must supply distinct, already isolated samples.
Building a reusable PON always requires named normal samples.

Other optional parameters:
- `--min-baseq` (default `20`): Minimum base quality for read bases
- `--min-mapq` (default `20`): Minimum mapping quality for reads
- `--truncate` (default `0.1`): Truncation percentile for PON sample inclusion
- `--pseudocount` (default `sys.float_info.epsilon`): Pseudocount for beta-binomial rate estimates
- `--prior-artifact-probability` (default `0.5`): Fallback artifact prior when the input record has no `SKUA_ARTIFACT_PRIOR`
- `--strict`: Fail before writing output if any VCF record cannot be annotated

Statistical parameters must be finite. `pseudocount` must be positive;
`truncate` must be in `(0, 1]`; and `prior_artifact_probability` must be in
`(0, 1)`. CLI annotation and the Python APIs reject invalid values before
publishing output, including when replacing an existing file.

The exported `compute_stats` API also requires `rho` in `(0, 1)` and
`0 < mu_min <= mu_max < 1`. These checks apply even at zero depth or when
per-sample evidence will replace the supplied `rho`. The statistical helpers
`estimate_rho` and `truncated_normal_evidences` use the same `truncate` contract;
their `pseudo` and `epsilon` parameters must be finite and positive, and
dispersion bounds must satisfy `0 < rho_min <= rho_max < 1`. Existing internal
numerical clipping remains in place for valid parameters.

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

Records with unavailable query names (`*` in SAM/BAM/CRAM, or `None` in the
Python alignment interface) cannot be assigned to fragments. Both singleton and
batch alignment collectors exclude them from usable evidence and count each
record as unusable with reason `missing_query_name`. These diagnostic counts
are per record because the number of fragments is unknown; unnamed records are
never merged together or treated as independent usable fragments. Named
fragments remain scoped to their read group.

MNV support requires consecutive aligned query bases across the entire reference
interval. An internal insertion, deletion, or reference skip makes that read
unusable for the simple MNV, regardless of the inserted bases' quality. Insertions
outside the MNV interval do not affect this continuity check.

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
- `SKUA_ASSESSMENT_STATUS`: `ASSESSED` or `INSUFFICIENT_EVIDENCE` for the selected case sample
- `SKUA_ASSESSMENT_REASONS`: Unmet evidence requirements; `.` when assessed

Output INFO fields:
- `SKUA_ARTIFACT_PRIOR`: Effective prior probability that the ALT allele is an artifact before Skua evidence
- `SKUA_STATUS`: Annotation outcome for every record. `ANNOTATED` records receive Skua evidence; unsupported records are retained with an `UNSUPPORTED_*` status and are not assigned new Skua evidence fields. `UNSUPPORTED_RECORD` includes records with no alternate allele (`ALT=.`).
- `SKUA_PON_SAMPLE_COUNT`: Number of normal samples included after truncation
- `SKUA_PON_ALT_FWD`, `SKUA_PON_ALT_REV`, `SKUA_PON_NON_ALT_FWD`, `SKUA_PON_NON_ALT_REV`: Aggregated read counts across normals
- `SKUA_PON_USABLE`, `SKUA_PON_UNUSABLE`: Aggregated usable/unusable counts
- `SKUA_PON_DISPERSION_FACTOR`: Beta-binomial dispersion parameter estimate

By default, unsupported records do not stop the run. Use `--strict` to reject any input containing one before an output file is created. VCF output is written to `--output` or standard output.

#### Assessment eligibility

`SKUA_STATUS=ANNOTATED` means the allele is supported and its evidence was
collected. Model eligibility is reported separately in the selected sample's
`SKUA_ASSESSMENT_STATUS`. An unsupported allele retains its `UNSUPPORTED_*`
record status and receives no new assessment or score fields. Unselected samples
have missing assessment values. The evidence-only Python APIs do not assess a
model; absence of an assessment must not be interpreted as `ASSESSED`.

**Scores and counts are retained for `INSUFFICIENT_EVIDENCE` records.** In
particular, a zero-case record still has log Bayes factor `0` and posterior equal
to its effective prior. Such a score is not an evidence-supported assessment.
PON counts continue to describe only the normals retained after truncation;
empty panels and all-normals-truncated sites have zero pooled usable depth.
`--strict` concerns allele support and does not reject insufficient evidence.

The following inclusive assessment minima can be configured with either live
normals or a cached PON. Depth means usable ALT plus non-ALT evidence, not ALT
support alone. Normal depths are pooled across retained samples after truncation.

| CLI option | Default | Requirement |
| --- | --- | --- |
| `--min-case-depth` | `1` | Total usable case depth |
| `--min-normal-depth` | `1` | Total usable retained-normal depth |
| `--min-normal-samples` | `0` | Retained normal sample count; `0` disables this extra requirement |
| `--min-case-strand-depth` | `0` | Usable case depth on **each** strand; `0` disables |
| `--min-normal-strand-depth` | `0` | Pooled retained-normal depth on **each** strand; `0` disables |

All limits must be integers. Total-depth minima must be at least `1`, so zero
evidence can never be made eligible; sample-count and strand minima may be `0`.
The defaults only exclude absent evidence and are **not assay-validated coverage
thresholds**. One-strand coverage is eligible by default; set the strand minima
to positive, assay-validated values to require both strands. Agree on production
depth, sample-count, and strand requirements through assay validation.
`ASSESSED` means those configured requirements were met, not that the assay or
variant call has been validated.

All failing requirements are reported in `SKUA_ASSESSMENT_REASONS`, in this
order: `CASE_DEPTH`, `NORMAL_DEPTH`, `NORMAL_SAMPLE_COUNT`,
`CASE_STRAND_DEPTH`, `NORMAL_STRAND_DEPTH`. No scores or counts are changed by
tightening the assessment thresholds.

For example, a downstream Python filter for a selected VCF sample can require
eligibility before applying an illustrative posterior cutoff:

```python
sample = record.samples["CASE"]
posterior = sample.get("SKUA_ARTIFACT_POSTERIOR")
accept = (
    record.info.get("SKUA_STATUS") == "ANNOTATED"
    and sample.get("SKUA_ASSESSMENT_STATUS") == "ASSESSED"
    and posterior is not None
    and posterior < 0.01  # Example only; validate the cutoff for your assay.
)
```

This rule cannot accept a zero-evidence record solely because its prior is low.
Posterior-only filters must be updated to check eligibility.

The exported `AssessmentThresholds` dataclass uses the same option names with
underscores. Pass it as `assessment_thresholds=` to `compute_stats`,
`annotate_vcf_with_normals`, `annotate_vcf_with_pon`, or
`annotate_vcf_to_json_with_normals`. `Stats.assessment_status` is an exported
`AssessmentStatus` string enum; `Stats.assessment_reasons` is a tuple of reason
strings. Python JSON output includes the same status and a reasons list under
`stats` (empty when assessed); it continues to contain supported alleles only.

When `compute_stats` receives `per_sample_evidences`, those samples are
authoritative: the retained pool supplies its normal counts, background summary,
scores, and eligibility, even for an empty list. With aggregate-only input, the
supplied pool is used as-is. If `min_normal_samples > 0`, aggregate-only input is
ineligible with reason `NORMAL_SAMPLE_COUNT_UNAVAILABLE`, because the number of
samples cannot be inferred from pooled counts. Supply per-sample evidence to
evaluate that requirement.

File outputs are transactional and no-clobber by default. Skua writes a sibling
temporary VCF and publishes it atomically only after annotation finishes;
`--force` permits replacement of an existing output. Existing `.tbi` and `.csi`
companions also require `--force`, even if the output VCF itself is absent.
Forced publication removes these indexes before replacing the VCF and restores
the previous files if publication fails. Skua does not rebuild them or index
new VCF outputs automatically; create a fresh index after annotation if you
need regional access to compressed output.

Without `--force`, input VCFs that already define generated `SKUA_*` INFO or
FORMAT annotations are
rejected. With `--force`, Skua removes those definitions and all associated
record and sample values before installing the canonical definitions and
computing fresh annotations. Unsupported records therefore retain only their
new `SKUA_STATUS` among generated Skua fields. `SKUA_ARTIFACT_PRIOR` is the
exception: it is a validated input field whose effective value is intentionally
preserved and written to the result. Non-Skua fields are preserved unchanged,
including per-allele genotype phasing (such as `0|1/2`) and missing genotype
alleles in both selected and unselected samples, even on retained unsupported
records.

### `pon`

Build reusable panel-of-normals artifacts, inspect their metadata, and
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
policy versions, quality thresholds, producer version, normal samples, and
stored reference identity and verification status. It
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

`validate` checks the BCF and CSI index, compatibility metadata, required FORMAT
fields, target alleles, duplicate targets, coordinate order, and per-sample
evidence counts. `--reference` additionally checks every PON REF allele;
`--targets` requires the supplied VCF to define exactly the PON targets. It
returns exit status `0` when valid and `1` when any check fails.

All PON evidence readers, including both cached annotation modes, require the
six count fields to declare `Number=1,Type=Integer`. Counts must be present,
scalar, nonnegative integers, and usable totals must equal the four strand-aware
ALT/non-ALT counts. Invalid counts are rejected without numeric coercion and
before annotation publishes or replaces output. Header compatibility checks are
shared with `validate` and do not require an additional full-panel scan.

Index validation queries every distinct target start through the CSI index and
compares the complete records starting there (including order and multiplicity)
with a sequential BCF scan. Overlapping records that start earlier are excluded
from that comparison. This checks retrieval at the panel's target sites; it does
not prove cryptographic identity of the index or test every possible interval.

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

The current evidence policy is **version 7**. Normal evidence is restricted to
read groups assigned to the normal sample. Records with unavailable query names
are excluded from usable fragment evidence, and internal insertions are excluded
from simple-MNV evidence. PON annotation and validation require the current
evidence-policy version.

The supported production workflow is paired-end targeted DNA mapped with
`bwa mem`. Mapping quality is filtered only by the numeric `--min-mapq` threshold:
values at or above the threshold pass, including MAPQ 255. Although the
[SAM specification, section 1.4](https://samtools.github.io/hts-specs/SAMv1.pdf)
defines 255 as unavailable mapping quality, Skua applies no special exclusion.
This rule applies to case and normal evidence, single-variant and batched
collection, and PON construction. Passing reads still undergo all other evidence
filters and normal fragment collapsing.

The PON schema is **version 2**, with required reference metadata.

#### Reference compatibility

For the contigs used by supported target alleles, Skua compares alignment `@SQ`
lengths and available `M5` checksums, plus available VCF contig `length` and `md5`
metadata. Conflicting values fail before output publication, including forced
replacement. Unused contigs do not need to match; separate-VCF cached annotation
checks the contigs used by that subset. Contig names must match exactly.

With `--reference`, Skua retains the target REF-base checks and also calculates
whole-contig sequence checksums in bounded chunks, following the
[SAM reference MD5 convention](https://samtools.github.io/hts-specs/SAMv1.pdf).
This identifies the sequence rather than the FASTA file's path or formatting.

VCF outputs and PONs contain `SKUA_REFERENCE_STATUS` and per-contig
`SKUA_REFERENCE` header records. `VERIFIED` means that all participating
alignment/PON identities have matching lengths and sequence checksums on the
used contigs, including the supplied FASTA when present. This verifies reference
metadata; it does not independently establish that the reads were aligned to
the declared sequence. Without a FASTA, target REF bases are not checked against
the full sequence.

Missing alignment checksums or reference metadata are allowed and reported as
`INSUFFICIENT_METADATA`, even when lengths agree or a FASTA supplies a checksum.
Names or matching target REF bases alone are never reported as verified
reference identity. PON construction retains this limitation per contig, and a
later case/FASTA cannot upgrade unverifiable normal metadata. `pon inspect`
reports the stored build status; `pon validate` checks available supplied
reference and target metadata and reports definite conflicts as validation
errors. Structural validity does not imply verified reference identity.

#### Run summary and PON metadata

Annotated VCFs contain one human-readable `SKUA_RUN` header record with
`SchemaVersion=1`. It records the Skua version, mode, evidence-policy version,
effective base/mapping-quality thresholds and read-selection
policies. Normal-model runs also record truncation, pseudocount, prior policy
and fallback, and all five assessment thresholds. For example, a case-only run
with default quality thresholds writes:

```text
##SKUA_RUN=<SchemaVersion="1",SkuaVersion="0.7.4",Mode="case_only",EvidencePolicyVersion="7",MinBaseQ="20",MinMapQ="20",CaseReadGroups="assigned_to_sample",NormalReadGroups="not_applicable">
```

`SkuaVersion` reflects the installed version. Cached annotation records the
PON's effective quality thresholds. VCF model annotation keeps each effective
prior in `INFO/SKUA_ARTIFACT_PRIOR`; `PriorPolicy=record_info_then_fallback`
identifies the selection rule and `PriorFallback` records its fallback value.
The existing reference metadata and per-variant annotations remain unchanged.

The summary contains no input paths or filenames, input-file hashes or sizes,
sample names, read-group IDs, or embedded PON build history. It adds no full-file checksum reads. Exact
input traceability belongs in your workflow records; Skua does not create a
separate manifest. Reference sequence checksums used by reference compatibility
checks are unchanged.

PONs contain `SKUA_PON` compatibility metadata, reference metadata, normal
sample columns, and per-normal evidence. These normal identities are
required for membership checks and are not copied into annotated VCF metadata.
PON schema 2 and evidence policy 7 define the supported format. Artifacts with
unsupported versions must be rebuilt from their original targets and alignments.

Forced reannotation replaces Skua annotations and run summaries.
Unrelated input VCF headers and sample columns are preserved, including any
names or paths already written by upstream tools. This change does not anonymize
an input VCF.

## Python API

The supported library API is available directly from `skua`. It accepts
substitutions, MNVs, and left-anchored simple insertions and deletions.

`annotate_vcf_to_json()` and `annotate_vcf_to_json_with_normals()` return a
JSON object with `run_summary` and `records` keys. The summary has
`schema_version: 1` and includes the reference compatibility result. Read the
result rows with `json.loads(payload)["records"]`. The standalone row formatters
and `render_annotation_results_json()` return a list of records.

The JSON wrappers retain their evidence-API semantics: they use all case
alignment reads without VCF sample selection, recorded as
`case_read_groups: "all_alignment_reads"`. The normal-model wrapper records
`model.prior.policy: "constant"` and includes each row's effective `artifact_prior`.
VCF-writing APIs instead select case read groups and use per-record VCF priors.
`evidence.normal_read_groups` records `assigned_to_sample`, or
`per_normal_selection` when headerless Python normals use all reads, without
listing individual normal identities. Case-only summaries use `not_applicable`.
Compare these policy fields before comparing outputs from the two APIs.

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

## Release validation

Pull requests and manual runs of **Publish Conda Package** validate the source
and build and exercise the Conda package without publishing. The source suite
runs on Python 3.11–3.14 with current pysam, and on Python 3.11 with the minimum
supported pysam 0.22.0. Tag publication uses the same tests on the exact tagged
revision and rejects disagreements between the tag, source, and recipe versions.

The installed-package check verifies the distribution and CLI versions, then
generates real indexed BAM and reference-backed CRAM fixtures. It exercises PON
build/validate and direct-normal, PON-only, and separate-VCF cached annotation,
checking known fragment counts and identical evidence and scores. It runs with
isolated Python imports and rejects a package loaded from the source checkout.
Upload requires every validation job to pass and retains the `conda` environment
approval controls. Only a tag run can upload; this workflow does not create tags
or GitHub releases.

## License

MIT. See [LICENSE](LICENSE) for details.
