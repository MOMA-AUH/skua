"""Measure CIGAR gap scaling and ordinary paired-read classification.

Run with the desired Skua checkout on PYTHONPATH, e.g.:
    PYTHONPATH=src python scripts/benchmark_cigar_gaps.py > cigar-gaps.json
Timing is measured without tracemalloc; peak Python allocation is measured
separately. These synthetic cases are not an end-to-end throughput benchmark.
"""

import argparse
from collections.abc import Callable, Iterator
import gc
import json
import platform
from statistics import median
from timeit import repeat
import tracemalloc

import pysam

from skua import __version__
from skua.evidence import (
    AlleleSupport, classify_variant_read, collect_evidence_from_alignment,
)


class InMemoryAlignment:
    """Exercise fragment counting without including BAM I/O in the measurement."""

    def __init__(self, reads: list[pysam.AlignedSegment]) -> None:
        self.reads = reads

    def fetch(self, contig: str, start: int, stop: int) -> Iterator[pysam.AlignedSegment]:
        return iter(self.reads)


def make_read(cigar: str, sequence: str, *, reverse: bool = False) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = "fragment"
    read.query_sequence = sequence
    read.query_qualities = [35] * len(sequence)
    read.flag = 147 if reverse else 99
    read.reference_id = 0
    read.reference_start = 100
    read.mapping_quality = 60
    read.cigarstring = cigar
    return read


def measure(
    run: Callable[[], None], *, iterations: int, repeats: int,
) -> dict[str, float | int]:
    run()  # Warm caches and check the expected classification before timing.
    seconds = median(repeat(run, number=iterations, repeat=repeats)) / iterations
    gc.collect()
    tracemalloc.start()
    try:
        run()
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    return {"seconds_per_call": seconds, "peak_python_bytes": peak}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spans", type=int, nargs="+", default=[10, 1000, 200_000])
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--ordinary-iterations", type=int, default=10_000)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if min(*args.spans, args.iterations, args.ordinary_iterations, args.repeats) < 1:
        parser.error("spans, iterations and repeats must be positive")

    results = []
    for operation in ("N", "D"):
        for span in args.spans:
            cigar = f"1M{span}{operation}1M"
            read = make_read(cigar, "TG")

            def classify() -> None:
                call = classify_variant_read(
                    read, ref_pos0=100, ref_base="A", alt_base="T",
                )
                assert call.support == AlleleSupport.ALT

            results.append({
                "case": cigar, "read_bases": 2,
                **measure(classify, iterations=args.iterations, repeats=args.repeats),
            })

    mates = [
        make_read("150M", "T" + "A" * 149, reverse=reverse)
        for reverse in (False, True)
    ]
    alignment = InMemoryAlignment(mates)

    def classify_pair() -> None:
        evidence = collect_evidence_from_alignment(
            alignment, contig="chr1", ref_pos0=100, ref_base="A", alt_base="T",
        )
        assert evidence.alt_forward == evidence.usable == 1

    results.append({
        "case": "overlapping_150M_mates", "read_bases": 300,
        **measure(classify_pair, iterations=args.ordinary_iterations, repeats=args.repeats),
    })
    print(json.dumps({
        "python": platform.python_version(), "platform": platform.platform(),
        "pysam": pysam.__version__, "skua": __version__, "settings": vars(args),
        "results": results,
    }, indent=2))


if __name__ == "__main__":
    main()
