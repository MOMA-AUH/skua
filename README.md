# skua

Implementation of the [shearwater](https://doi.org/10.1093/bioinformatics/btt750) statistical model to assess somatic variant evidence in aligned reads with support for SNV, MNV, and INDEL variants. The **shearwater** authors named their algorithm after seabirds that fly long distances over the ocean, watching the water closely and eventually dive into the water to catch prey. Due to the heavy reuse of the algorithmic core, it is only natural to name this **skua** — a seabird that hunts and steals from other birds.

## Installation

The recommended way to install **skua** is via [conda](https://docs.conda.io/), using the `micknudsen` channel:

```bash
conda install micknudsen::skua
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
- `--reference`: Reference FASTA file, required when any input alignment is CRAM
- `--output`: Optional output VCF path; if omitted, output is written to `stdout`

Other optional parameters:
- `--min-baseq` (default `20`): Minimum base quality for read bases
- `--min-mapq` (default `20`): Minimum mapping quality for reads
- `--truncate` (default `0.1`): Truncation percentile for PON sample inclusion
- `--pseudocount` (default `sys.float_info.epsilon`): Pseudocount for beta-binomial rate estimates
- `--prior-variant-probability` (default `0.5`): Prior probability for variant model

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
- `SKUA_PON_SAMPLE_COUNT`: Number of normal samples included after truncation
- `SKUA_PON_ALT_FWD`, `SKUA_PON_ALT_REV`, `SKUA_PON_NON_ALT_FWD`, `SKUA_PON_NON_ALT_REV`: Aggregated read counts across normals
- `SKUA_PON_USABLE`, `SKUA_PON_UNUSABLE`: Aggregated usable/unusable counts
- `SKUA_PON_DISPERSION_FACTOR`: Beta-binomial dispersion parameter estimate

## Requirements

- Python ≥ 3.10
- pysam ≥ 0.22

## License

MIT. See [LICENSE](LICENSE) for details.