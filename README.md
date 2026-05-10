# skua

Implementation of the [Shearwater](https://doi.org/10.1093/bioinformatics/btt750) statistical model to assess somatic variant evidence in aligned reads with support for SNV, MNV, and INDEL variants. The Shearwater authors named their algorithm after seabirds that fly long distances over the ocean, watching the water closely and eventually dive into the water to pick up prey. Due to the heavy reuse of the algorithmic core, it is only natural to name this **skua** — a seabird that hunts and steals from other birds.

It takes a VCF file of candidate variants, a case alignment file (BAM/CRAM), and a panel of normal alignments, and outputs an annotated VCF with per-sample read counts, quality metrics, and artifact posteriors.

## Installation

Since **skua** is still in active development, installing from source is currently the only option.

```bash
git clone https://github.com/micknudsen/skua.git
cd skua
pip install -e .
```

## Quick start

### Verify with panel of normals
```bash
skua verify \
  --vcf input.vcf \
  --alignment case.bam \
  --normal-list normals.txt \
  --output output.vcf.gz
```

Where `normals.txt` is a text file with one BAM/CRAM path per line:
```
normal1.bam
normal2.bam
normal3.bam
```

### Print to stdout
```bash
skua verify --vcf input.vcf --alignment case.bam --normal-list normals.txt
```

## Input requirements

- **VCF file** (`--vcf`): Single-ALT VCF with candidate variants (SNVs, indels, MNVs supported)
- **Alignment file** (`--alignment`): BAM or CRAM file with aligned reads
- **Reference** (`--reference`): FASTA reference genome; required only if using CRAM input
- **Normal list** (`--normal-list`): Required; text file listing one BAM/CRAM path per line

## Output fields

All output is in VCF format with added FORMAT and INFO fields:

### FORMAT fields

Read-count fields:
- `SKUA_ALT_FWD`: Count of ALT-supporting reads on forward strand
- `SKUA_ALT_REV`: Count of ALT-supporting reads on reverse strand
- `SKUA_NON_ALT_FWD`: Count of non-ALT reads on forward strand
- `SKUA_NON_ALT_REV`: Count of non-ALT reads on reverse strand
- `SKUA_USABLE`: Total usable reads at this locus
- `SKUA_UNUSABLE`: Total unusable reads (low quality, indels at locus, etc.)

Model-score fields:
- `SKUA_ARTIFACT_POSTERIOR`: Posterior probability of artifact model (0–1)
- `SKUA_BAYES_FACTOR`: Bayes factor comparing artifact vs. variant models

### PON INFO fields
- `SKUA_PON_SAMPLE_COUNT`: Number of normal samples included after truncation
- `SKUA_PON_ALT_FWD`, `SKUA_PON_ALT_REV`, `SKUA_PON_NON_ALT_FWD`, `SKUA_PON_NON_ALT_REV`: Aggregated read counts across normals
- `SKUA_PON_USABLE`, `SKUA_PON_UNUSABLE`: Aggregated usable/unusable counts
- `SKUA_PON_DISPERSION_FACTOR`: Beta-binomial dispersion parameter estimate

## Options

```bash
skua verify --help
```

Key parameters:
- `--min-baseq` (default 20): Minimum base quality for read bases
- `--min-mapq` (default 20): Minimum mapping quality for reads
- `--truncate` (default 0.1): Truncation percentile for PON sample inclusion
- `--pseudocount` (default ε): Pseudocount for beta-binomial rate estimates
- `--prior-variant-probability` (default 0.5): Prior probability for variant model

Truncation controls how conservative the panel-of-normals aggregation is at each site. A normal sample is included only if its ALT fraction is strictly less than `--truncate`. With `--truncate 0.1`, normals with ALT fraction `< 0.1` are kept and normals with ALT fraction `>= 0.1` are excluded.

## Requirements

- Python ≥ 3.10
- pysam ≥ 0.22

## License

MIT. See [LICENSE](LICENSE) for details.