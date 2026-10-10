# Indel representation support

Issue #19 is addressed with bounded reference-based equivalence matching.
Supply the alignment reference FASTA through `--reference` or `reference_path=`.
Exact CIGAR matches remain supported. A shifted match requires one explicit
insertion or deletion within 100 bases of the VCF anchor, the same alternate
sequence in the supplied reference, and a clean, quality-passing read span
covering both representations and their flanks. This is not local assembly.

For reference `CATATATG`, inserting `AT` after `C` gives `CATATATATG`.
The same sequence can be represented by CIGAR `2M2I6M`, which inserts `TA`
after the first `A`. Deleting `AT` after `C` gives `CATATG`; CIGAR `2M2D4M`
instead deletes the equivalent `TA`. The integration tests include both
directions of shifting, both strands, and a deletion covering the VCF anchor.

## BWA-MEM experiment

Run from an environment with Skua, pysam, and BWA installed:

```bash
python scripts/characterize_indel_representations.py \
  --bwa bwa --output-dir /tmp/skua-indel-experiment
```

The output directory must be new. The script retains the synthetic reference,
paired FASTQs, SAM, sorted/indexed BAM, BWA log, and JSON results. It uses random
seed 19, four independent 1,000-base contigs, `A` and `AT` repeats of six units,
and one-unit insertions/deletions. Each allele has 40 error-free fragments of
220 bases, sequenced as overlapping 150-base mates at base quality 40. Commands:

```text
bwa index reference.fa
bwa mem -t 1 -I 220,10,250,190 -R '@RG\tID:synthetic\tSM:synthetic' reference.fa reads.1.fq reads.2.fq
```

Observed with **BWA 0.7.19-r1273**, pysam 0.24.1, and Python 3.14.7:

| Allele | Accepted records | CIGAR anchor offset from leftmost target | ALT fragments, exact / reference mode | ALT fragments for equivalent target shifted +3 bases, exact / reference mode |
| --- | ---: | ---: | ---: | ---: |
| Homopolymer insertion | 80 | 0 for all 80 | 40 / 40 | 0 / 40 |
| Homopolymer deletion | 80 | 0 for all 80 | 40 / 40 | 0 / 40 |
| Tandem insertion | 80 | 0 for all 80 | 40 / 40 | 0 / 40 |
| Tandem deletion | 80 | 0 for all 80 | 40 / 40 | 0 / 40 |

BWA placed all four clean examples at the leftmost anchor; this experiment did
not find spontaneous shifting against the leftmost target. Expressing the target
three bases to the right illustrates the representation dependency using those
same alignments. The new mode recognizes all 40 fragments in either target
representation. Overlapping mates collapse to 40 fragments rather than 80
observations. The script also checks singleton/batch agreement.

These observations characterize these synthetic inputs and this BWA build.
They do not measure sensitivity on real targeted data, other BWA options,
sequencing errors, complex haplotypes, long repeats, or reads ending inside the
repeat. Synthetic BAM tests independently exercise shifted CIGARs relative to
leftmost targets, quality failures, complex boundaries, the search bound, and
direct-normal/cached-PON parity.

## Compatibility and limits

PON evidence policy is now 8; policy-7 artifacts must be rebuilt. PONs record
`IndelMatching=reference` or `exact_anchor`, and fresh case evidence uses the
stored mode. Reference-mode PONs require a compatible FASTA at annotation time.
Exact-mode PONs retain exact matching even when a FASTA is supplied later.
Target identifiers in PON joins are still exact, so equivalent VCF records at
different positions are not interchangeable PON keys.

The read must overlap the VCF anchor and span the comparison. Matching does not
rescue shifts beyond 100 bases, soft-clipped events, reference skips, ambiguous
reference sequence, or combinations of multiple events. Existing exact-anchor
classification applies when no shifted equivalence is established. This bounded
behavior is the selected implementation scope for v1; assay-specific sensitivity
remains a separate validation task.
