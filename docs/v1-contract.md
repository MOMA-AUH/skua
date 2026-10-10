# Skua v1 contract

This is the supported interface and compatibility policy for Skua 1.x. The
initial baseline is version **1.0.0**, PON schema **2**, evidence policy **8**,
and run-summary schema **1**. See the [migration notes](releases/v1.0.0.md) for
the transition from 0.x and the [README](../README.md) for command examples.

## Workflow and inputs

Skua assesses supplied alleles using case evidence and a panel of normals. The
intended workflow is paired-end targeted DNA mapped with `bwa mem`: prepare a
fixed, coordinate-sorted target VCF and distinct normal samples, build and
validate a reusable PON, then annotate each case at all or a subset of those
targets. Live-normal annotation is also supported. The CLI requires a normal
source; case-only evidence annotation is available through Python.

Supported alleles are biallelic A/C/G/T substitutions (SNVs and equal-length
MNVs) and simple left-anchored insertions/deletions with one allele of length
one and a shared first base. `Variant` uppercases alleles; it does not left-align
or minimize them. Contig names are exact. PON target matching requires the same
contig, position, REF and ALT representation, even when two indels are
biologically equivalent.

VCF readers accept files readable by pysam (VCF, compressed VCF and BCF). The
CLI writes `.vcf` or `.vcf.gz`, or VCF to standard output when `--output` is
omitted. Alignments must support indexed regional access (BAM or CRAM); CRAM
requires a compatible FASTA through CLI `--reference`. Python callers open and
close their own alignment handles, including configuring CRAM decoding. Keep
handles open until iterators are consumed. Target VCFs need no index for
sequential reading; reusable PONs are BCF with a companion `.bcf.csi` index.

Multiallelic records, symbolic alleles, breakends, spanning deletion `*`,
missing ALT, non-ACGT alleles and complex replacement indels are not scored.
VCF writers retain unsupported records with a status unless `strict=True` /
`--strict` requests rejection before output. Evidence iterators and JSON
wrappers skip them. PON construction always rejects unsupported targets and
requires unique alleles in coordinate order (using the header's contig order).

The initial release does not promise candidate discovery, automatic FILTER or
genotype decisions, multiallelic scoring, UMI consensus generation, or internal
parallel execution. Single-end, long-read and RNA workflows are outside the
supported production workflow. Batched fetches are an optimization, not a
parallel API. Caller-side workflow orchestration remains possible. Stable
interfaces do not establish assay-specific validity or coverage cutoffs.

## Python API

Import supported symbols directly from `skua`. The complete boundary is
`skua.__all__`, listed below. Other submodule functions, classes and constants
are implementation details, even when their names lack an underscore. Their
import paths are not covered by the 1.x promise.

| Export | Contract |
| --- | --- |
| `Variant`, `VariantKind` | Validated allele model: `contig`, zero-based `ref_pos0`, `ref`, `alt`, and derived `kind`. `from_vcf_fields(contig=..., pos1=..., ref=..., alt=...)` accepts one-based VCF positions. Kinds are `substitution`, `insertion`, `deletion`. Invalid alleles/positions raise `ValueError`. |
| `AlleleSupport`, `ReadAlleleCall`, `UnusableReason` | Read-call data: support (`alt`, `non_alt`, `unusable`), `is_reverse`, optional `reason`, `observed_base`, `base_quality`. Reasons are `low_mapq`, `low_baseq`, `missing_baseq`, `no_base_at_site`, `invalid_base`, `conflicting_mates`, `missing_query_name`. |
| `AggregatedEvidence` | Counts `alt_forward`, `alt_reverse`, `non_alt_forward`, `non_alt_reverse`, `usable`, `unusable`, and `unusable_by_reason`. |
| `PonAnnotation` | `case_evidence`, ordered per-normal `normal_evidences`, and their **untruncated** `normal_aggregate_evidence`; this is evidence, not a model assessment. |
| `AnnotationStatus` | Record outcomes listed in the VCF contract below. |
| `AssessmentStatus`, `AssessmentThresholds`, `Stats` | Model eligibility, its five configurable minima, and scores/counts returned by `compute_stats`; see below. |
| `PonArtifactMetadata` | `schema_version`, `evidence_policy_version`, `min_baseq`, `min_mapq`, producer `skua_version`, ordered `sample_names`, `reference_identity`, `indel_matching`. Read the nested reference result through `.status`, `.contigs` and `.as_dict()`; its concrete submodule constructors are not a public API. Each contig exposes `name`, `length`, `md5`, `verified`. |
| `annotate_variant` | Open alignment + `Variant` → `AggregatedEvidence`. |
| `annotate_variants`, `annotate_variants_from_vcf` | Open alignment + iterable of variants or VCF path → iterator of `(Variant, AggregatedEvidence)` in input order. |
| `annotate_variant_with_normals` | Open case alignment + `Variant`, with `normal_alignments=` → `PonAnnotation`. |
| `annotate_variants_with_normals`, `annotate_variants_from_vcf_with_normals` | Batch/VCF equivalents → iterator of `(Variant, PonAnnotation)` in input order. |
| `annotate_variants_from_pon` | Open case alignment + PON path → iterator of `(Variant, PonAnnotation)` in PON order, using cached thresholds and matching policy. |
| `annotate_vcf` | Case-only VCF writer; returns `None`. |
| `annotate_vcf_with_normals`, `annotate_vcf_with_pon` | Live/cached model VCF writers; return `None`. Cached mode accepts optional `vcf_path=` for a target subset. All writers require `output_path=`; `"-"` means stdout. |
| `annotate_vcf_to_json`, `annotate_vcf_to_json_with_normals` | Return JSON text with `run_summary` and `records`; optionally write the same text to `output_path=`. |
| `build_pon` | Target VCF + `normal_alignments=` + `output_path=` → BCF/CSI pair; returns `None`. |
| `read_pon_metadata` | PON path → compatible `PonArtifactMetadata`; validates header, not the full artifact/index. |
| `read_pon_evidence` | PON path → iterator of `(Variant, tuple[AggregatedEvidence, ...])` in artifact order and metadata sample order; validates counts as consumed. |
| `compute_stats` | Case and pooled normal evidence, optionally `per_sample_evidences=` → `Stats`. |
| `__version__` | Installed package version string. |

### Sample selection

VCF writers resolve exactly one case sample from VCF sample names and alignment
read-group `SM` tags. `sample_name=` / CLI `--sample` resolves ambiguity and
must identify a sample in both inputs. For site-only inputs, it selects an
alignment sample and adds that sample to output. A selected sample must have
assigned read-group IDs; only reads whose `RG` names one of those IDs count.
Untagged, unknown and unassigned groups are excluded. Other VCF samples keep
their original values and receive missing generated Skua FORMAT values.
PON-only annotation removes normal sample columns and outputs the selected case.

Evidence APIs (`annotate_variant*` / `annotate_variants*`) instead use **all case
read groups by default**, without VCF sample selection. Pass
`allowed_read_group_ids=frozenset({...})` to isolate case reads. JSON wrappers
use all case read groups and expose no sample/RG selection argument; use an
already isolated case alignment. This difference is recorded in run summaries.

Each normal handle represents one biological sample with one distinct `SM` and
unambiguous read-group IDs. Only its assigned groups count. Duplicate handles,
local files (including symlinks/hardlinks), normal sample names, and case
membership in the panel are rejected when identifiable. Merge libraries for
one sample upstream. There is no automatic leave-one-out analysis. Identity
checks rely on truthful metadata; they do not infer identity from sequence.

Low-level Python calls also allow headerless alignment-like normal objects;
all their reads are used and sample identity cannot be verified. The caller
must supply distinct, isolated samples. Building reusable PONs requires named
normal samples. Empty live-normal lists are allowed by Python, yielding
insufficient normal evidence; CLI normal lists and PON builds require at least
one normal.

### Validation and evidence

VCF-based entry points check available alignment indexes, target contigs,
reference metadata conflicts and, when supplied, FASTA REF bases. Direct
`Variant` evidence collectors are lower-level: they do not perform the complete
VCF preflight or promise the same eager parameter validation. Supply valid
thresholds and compatible inputs. Generator errors can arise on iteration;
consume an iterator before treating a run as successful. Python validation
errors generally use `ValueError`, output conflicts use `FileExistsError`, and
underlying I/O/pysam errors can propagate. Exact exception text is not stable.

Evidence thresholds default to `min_baseq=20`, `min_mapq=20` and must be
nonnegative. MAPQ is tested numerically, including 255. Accepted alignment
records are mapped primary records from proper pairs, with mapped mates, and
without secondary, supplementary, QC-fail or duplicate flags. Records rejected
by flags or sample selection do not contribute usable or unusable counts.

Overlapping mates collapse by `(RG, query name)`: agreeing usable mates count
once on the first mate's strand; a usable mate wins over an unusable mate;
conflicting usable mates count once as unusable. Unnamed records each contribute
an unusable `missing_query_name` diagnostic because fragment identity is unknown.
Counts therefore describe collapsed evidence, with per-record diagnostics for
missing names. Non-ALT means usable support other than the requested ALT, not
necessarily reference-only support.

MNV evidence must span consecutive aligned query bases without internal indels
or reference skips. Simple indel evidence requires explicit CIGAR events and
aligned flanks. With `reference_path=` / `--reference`, equivalent placements
within 100 reference bases on either side of the anchor can count when the read
spans both representations and all comparison bases pass quality checks. The
read must still overlap the VCF anchor. Without a FASTA, matching is exact-anchor.
JSON wrappers expose no reference argument and use exact-anchor evidence.
See [indel limits and examples](indel-representations.md).

### Model and JSON results

`compute_stats` returns case/normal channel counts, `background_rate_by_channel`,
`expected_case_counts`, `log_bayes_factor_artifact_vs_variant`,
`artifact_posterior`, `dispersion_rho`, `pseudocount`, `assessment_status`, and
`assessment_reasons`. The two model score fields have type `float | None`:
they are `None` whenever status is `INSUFFICIENT_EVIDENCE`. Counts and diagnostic
summaries remain available. For assessed evidence, the log Bayes factor is the
natural log of artifact-versus-variant evidence; lower posterior means stronger
variant support.

Statistical inputs must be finite: `0 < truncate <= 1`, `pseudocount > 0`,
`0 < prior_artifact_probability < 1`, `0 < rho < 1`, and
`0 < mu_min <= mu_max < 1`. Defaults are truncation `0.1`, pseudocount
`sys.float_info.epsilon`, prior `0.5`, rho `1e-4`, and mu bounds `1e-6` and
`1 - 1e-6`. Validation applies even at zero depth. When evidence is assessed, the
model bounds the prior to `[1e-12, 1 - 1e-12]` for numerical stability; the recorded
input prior remains unchanged. Supplied per-sample evidence
is authoritative, including an empty list: retained normals determine pooled
counts, dispersion, scores and eligibility, replacing the aggregate argument.
Truncation retains normals whose ALT fraction, adjusted with machine epsilon
in the numerator and denominator, is strictly below `truncate`.

`AssessmentThresholds` has integer `min_case_depth=1`, `min_normal_depth=1`,
`min_normal_samples=0`, `min_case_strand_depth=0`, `min_normal_strand_depth=0`.
The first two must be at least one; the others at least zero (zero disables).
Normal thresholds apply after truncation. These are inclusive usable-depth
minima, not ALT-count thresholds. An aggregate-only panel cannot satisfy a
positive sample minimum because its sample count is unknown.

JSON is an object, not a bare row list: use `json.loads(payload)["records"]`.
Each supported row contains `contig`, one-based `pos1`, `ref`, `alt`, and
`counts.case` with the seven evidence fields. Normal-model rows additionally
contain post-truncation `counts.normal`, `artifact_prior`, and `stats` with
`artifact_posterior`, `log_bayes_factor_artifact_vs_variant`, `dispersion_factor`,
`pon_sample_count`, `assessment_status`, `assessment_reasons`. Both score keys
remain present with JSON `null` values when evidence is insufficient. Reasons
are a list, empty when assessed. JSON model priors are constant API arguments; record INFO
priors are not used. Cached PON counts do not store `unusable_by_reason`; their
Python evidence objects return an empty reason map.

`run_summary` contains schema version 1, producer version, mode, effective
evidence/read-selection policy, model settings (or `null` for case-only), and
reference compatibility. It is not a full input manifest. JSON output files
are ordinary writes that overwrite existing paths; the VCF transaction and
`force` contract below does **not** apply to JSON wrappers.

## CLI and VCF contract

The supported commands are `skua annotate`, `skua pon build`, `skua pon inspect`,
`skua pon validate`, and `skua --version`. Long option names, requiredness,
defaults and documented meaning are part of the compatibility boundary.
`--help` lists each command's complete options. Model options map to the Python
parameters above, using hyphens instead of underscores for the five assessment
thresholds. Cached annotation rejects explicit `--min-baseq` / `--min-mapq`,
even if they equal the stored values: the PON owns those thresholds.

| Exit status | Meaning |
| --- | --- |
| `0` | Successful command, help/version display, or successful PON validation. `inspect` success only means metadata was read, not that a PON is compatible. Unsupported records and insufficient evidence do not fail ordinary annotation. |
| `1` | `pon validate` reports invalid/incompatible input. Unhandled runtime/I/O errors can also terminate nonzero (normally 1); this is not an exclusive validation-error category. |
| `2` | Argument parsing or validation reported through the CLI parser, including handled annotation/build `ValueError` and output conflicts. |

Use exit status and documented structured output, not exact diagnostics, help
layout, human-readable inspection text or stderr from pysam. `pon inspect
--json` exposes header facts; `pon validate --json` exposes `valid`, `errors`,
and `inspection`. Treat every nonzero exit as failure.

VCF writers preserve input record order, alleles, IDs, QUAL, FILTER and unrelated
INFO/FORMAT values, including genotype phasing. Without `force`, existing
generated Skua annotations are rejected. With `force`, they are removed from
all samples/records and recomputed for the selected sample. The validated input
`SKUA_ARTIFACT_PRIOR` is retained. File output uses sibling temporary files and
publication after successful annotation; existing output or `.tbi`/`.csi`
companions require `force`. Forced publication removes stale indexes and restores
previous files if publication fails. New VCFs are not indexed automatically.
Input and output paths must differ. Stdout cannot be rolled back after writes.

The following generated fields are stable by name, location, Number and Type.
Case-only VCF annotation emits the six count fields and `SKUA_STATUS`; model
fields require live normals or a cached PON. Counts are nonnegative integers.

| Location / fields | Number | Type | Meaning |
| --- | --- | --- | --- |
| FORMAT `SKUA_ALT_FWD`, `SKUA_ALT_REV`, `SKUA_NON_ALT_FWD`, `SKUA_NON_ALT_REV` | 1 | Integer | Selected case strand counts. |
| FORMAT `SKUA_USABLE`, `SKUA_UNUSABLE` | 1 | Integer | Case usable/unusable evidence; usable equals the four strand counts' sum. |
| FORMAT `SKUA_LBF` | 1 | Float | Natural log Bayes factor, artifact versus variant; positive favors artifact, negative favors variant. Missing (`.`) when evidence is insufficient. |
| FORMAT `SKUA_ARTIFACT_POSTERIOR` | 1 | Float | Artifact probability; missing (`.`) when evidence is insufficient. |
| FORMAT `SKUA_ASSESSMENT` | 1 | String | `ASSESSED` or `INSUFFICIENT_EVIDENCE`. |
| FORMAT `SKUA_REASONS` | . | String | All unmet requirements; present on every model-annotated record, with value `.` when assessed. |
| INFO `SKUA_STATUS` | 1 | String | Annotation outcome for every record. |
| INFO `SKUA_ARTIFACT_PRIOR` | A | Float | Effective allele prior, strictly between 0 and 1. |
| INFO `SKUA_PON_SAMPLE_COUNT` | 1 | Integer | Number of retained normals. |
| INFO `SKUA_PON_ALT_FWD`, `SKUA_PON_ALT_REV`, `SKUA_PON_NON_ALT_FWD`, `SKUA_PON_NON_ALT_REV`, `SKUA_PON_USABLE`, `SKUA_PON_UNUSABLE` | 1 | Integer | Pooled normal evidence after truncation. |
| INFO `SKUA_PON_RHO` | 1 | Float | Beta-binomial dispersion parameter rho estimated at annotation time. |

Four field names changed at the v1 boundary; see the
[migration table](releases/v1.0.0.md#vcf-field-renames). V1 emits only the new
names, and forced reannotation removes the old generated fields.

`AnnotationStatus` values are `ANNOTATED`, `UNSUPPORTED_RECORD`,
`UNSUPPORTED_MULTIALLELIC`, `UNSUPPORTED_SYMBOLIC_ALLELE`, `UNSUPPORTED_BREAKEND`,
`UNSUPPORTED_SPANNING_DELETION`, `UNSUPPORTED_COMPLEX_ALLELE`, and
`UNSUPPORTED_NON_STANDARD_ALLELE`. Unsupported records receive no new evidence,
model scores or assessment fields. `--strict` rejects unsupported alleles,
not insufficient evidence.

Assessment reasons are `CASE_DEPTH`, `NORMAL_DEPTH`, `NORMAL_SAMPLE_COUNT`,
`NORMAL_SAMPLE_COUNT_UNAVAILABLE`, `CASE_STRAND_DEPTH`, `NORMAL_STRAND_DEPTH`.
The reasons FORMAT field is present on every supported record annotated with
normals: `.` for an assessed case or an unselected sample, otherwise the list
of unmet requirements. Python/JSON reasons remain an empty tuple/list when
assessed. Unsupported records and case-only outputs have no assessment fields.

`ANNOTATED` does not imply `ASSESSED`. For every `INSUFFICIENT_EVIDENCE` result,
both model scores are missing (`.`) in VCF, `None` in Python, and `null` in JSON.
Counts, the input prior, status and reasons remain available. This includes
zero case/normal depth and any unmet configured depth, sample-count or strand
requirement. Tightening thresholds can suppress scores but does not change
counts or the scores of records that remain assessed. Forced reannotation also
replaces old numeric scores with missing values when the new assessment fails.
Missing assessment is not an implicit pass.

An explicit valid record prior overrides the fallback; absent/missing prior
uses the fallback. In separate-VCF cached mode that VCF owns the prior, whereas
PON-only mode uses the prior retained in the PON. Priors must be independent of
the evidence being scored. Skua does not derive case-specific priors or alter
FILTER based on scores.

`SKUA_RUN` header schema 1 records the effective settings and versions, without
input paths, sample identities or full-file hashes. `SKUA_REFERENCE_STATUS`
is `VERIFIED` or `INSUFFICIENT_METADATA`; per-contig `SKUA_REFERENCE` records
contain `ID`, `Verified`, and available `Length`/`MD5`. Verification compares
declared alignment/PON identities on used contigs and a supplied FASTA; it does
not prove how reads were aligned. Missing normal metadata cannot be upgraded by
later supplying case/FASTA metadata. See [reference compatibility](../README.md#reference-compatibility).

## PON compatibility and rebuild decision

The initial supported pair is **schema 2 / evidence policy 8 only**. Older
schema versions, policy 7 and earlier, unknown future versions and missing
required metadata are rejected by annotation/evidence readers and validation.
There is no legacy reader, automatic migration or header-only upgrade. `pon
inspect` remains permissive for diagnosis. Changing version numbers by hand
cannot repair evidence. Rebuild incompatible PONs from original targets and
normal alignments with the installed Skua.

`SKUA_PON` contains `SchemaVersion`, `EvidencePolicyVersion`, `MinBaseQ`,
`MinMapQ`, `SkuaVersion`, `IndelMatching`. Schema 2 also requires reference
metadata. `SkuaVersion` records the producer; compatibility is checked using
schema/policy and input metadata, not equality of producer package versions.
Thus valid v0.8.2 schema-2/policy-8 artifacts remain usable in v1.0.0.

Each normal sample column stores FORMAT `SKUA_PON_AF`, `SKUA_PON_AR`,
`SKUA_PON_NF`, `SKUA_PON_NR`, `SKUA_PON_U`, `SKUA_PON_X`, all
`Number=1,Type=Integer`: ALT forward/reverse, non-ALT forward/reverse, usable,
unusable. Values must be scalar nonnegative integers, with U = AF + AR + NF + NR.
Missing or inconsistent counts fail; they are not coerced. Reason breakdowns
and model scores are not cached. Construction publishes BCF and CSI after both
are complete. `pon validate` scans records and compares indexed retrieval at
every distinct target start with sequential records; metadata/evidence readers
do not replace this full validation. `--targets` requires the same targets in
the same order, whereas cached annotation can select a subset.

`IndelMatching=reference` requires a compatible FASTA at case annotation.
`IndelMatching=exact_anchor` stays exact-anchor even if a FASTA is later supplied
for validation or CRAM decoding. Case and normal matching policies must agree.

| Change | Required action |
| --- | --- |
| Evidence-classification fix or changed read selection, fragment collapsing, indel matching or quality policy | Increment evidence policy when releasing a change that can alter stored counts; rebuild affected panels before reuse under that policy. Current readers require an exact policy match for every panel, even one containing only unaffected allele classes. |
| PON storage/schema version change | Rebuild unless a future release explicitly supplies and documents a supported migration. Initial v1 supplies none. |
| Reference assembly/sequence or contig identity changes | Rebuild from compatible targets and alignments. Renaming/reformatting an identical FASTA alone does not require rebuilding; identity uses sequence checksums, not paths. |
| Normal membership, RG/SM assignments, alignment processing or quality thresholds change | Rebuild; these determine cached counts and identities. |
| Enable reference indel matching for an exact-anchor panel | Rebuild with the FASTA for both normal and case evidence. |
| Add targets or change an allele's representation | Build a panel containing those exact target keys. An existing panel can serve an unchanged subset. |
| Change truncation, pseudocount, prior, assessment minima, or scoring-only implementation | Reannotate; no rebuild solely for these changes because raw per-normal counts are stored. The panel must still satisfy current schema/policy checks. |
| Upgrade package version with unchanged schema, evidence policy and relevant inputs | Validate and reuse the PON; producer version alone does not require rebuilding. |

## Compatibility within 1.x

Existing documented top-level imports, call signatures, return-field meanings,
CLI options/defaults and VCF field names/types retain their meaning within 1.x.
Removal, renaming or incompatible reinterpretation requires a major version.
Minor releases may add optional APIs/options, fields and status/reason values;
consumers should tolerate unknown fields and handle unknown statuses
conservatively instead of treating them as success. JSON object key order,
whitespace, VCF header ordering and diagnostic wording are not interfaces.

Correctness fixes may change evidence or numeric scores in patch/minor releases;
bitwise-identical results across releases are not promised. Release notes must
identify changed results and whether to reannotate or rebuild. A count-changing
fix increments evidence policy; changing PON layout/meaning increments schema.
Changing the run-summary structure incompatibly increments its schema. These
artifact versions are independent of the package major version: the 1.x API
promise does not promise that every historical PON remains readable forever.
Pin the producer version, input/reference identities, settings and PON versions
in workflow records when reproducing results.
