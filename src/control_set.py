"""
Control set construction, for FDR estimation.

Paper reference (Methods, "Sequence data sets"):
"The Control set was constructed from modern sequencing libraries to have the same
number of reads and length distribution as the Saqqaq library. Reads were randomly
sampled and truncated to match Saqqaq reads."

This builds a new BAM-like read set (as a list of pysam.AlignedSegment-derived
records written to a new BAM) by: (1) sampling reads at random from a modern/source
BAM, (2) truncating each sampled read to a length drawn from the target (Saqqaq)
read-length distribution. Because truncation changes alignment, this script
re-derives the truncated read's genomic placement from its original 5' mapped
position (keeping the same start coordinate and strand -- consistent with the
paper's aim of preserving positional/coverage behavior while matching length stats,
not literally re-aligning truncated reads from scratch).

If you want a stricter reproduction (re-align truncated reads rather than just
clipping coordinates), swap the write step for an actual re-alignment (e.g. via
bwa) on the truncated FASTQ sequences -- flagged here as the more "literal" option
if fidelity matters more to you than speed.
"""

import random
from typing import Iterable, List

import numpy as np
import pysam


def get_read_length_distribution(bam_path: str, chroms=None) -> List[int]:
    """Collect mapped read lengths from the target (e.g. Saqqaq) BAM."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    lengths = []
    iterator = bam.fetch() if chroms is None else (
        r for c in chroms for r in bam.fetch(c)
    )
    for read in iterator:
        if read.is_unmapped:
            continue
        if read.query_alignment_length:
            lengths.append(read.query_alignment_length)
    bam.close()
    return lengths


def build_control_bam(source_bam_path: str, target_length_dist: List[int],
                       out_bam_path: str, n_reads: int = None, seed: int = 0) -> None:
    """
    Sample `n_reads` (default: len(target_length_dist)) reads at random from
    source_bam_path, truncate each to a length drawn (with replacement) from
    target_length_dist, and write to out_bam_path preserving each read's original
    reference start position and strand (see module docstring for rationale/caveat).
    """
    rng = random.Random(seed)
    n_reads = n_reads or len(target_length_dist)

    src = pysam.AlignmentFile(source_bam_path, "rb")
    out = pysam.AlignmentFile(out_bam_path, "wb", template=src)

    # reservoir-style random sample of mapped reads from the source
    all_reads = [r for r in src.fetch() if not r.is_unmapped and r.query_alignment_length]
    if len(all_reads) < n_reads:
        raise ValueError(
            f"Source BAM only has {len(all_reads)} usable mapped reads, need {n_reads}."
        )
    sampled = rng.sample(all_reads, n_reads)

    for read in sampled:
        target_len = rng.choice(target_length_dist)
        orig_len = read.query_alignment_length
        new_len = min(target_len, orig_len)  # can't truncate to *longer* than original

        # truncate query sequence/qualities and shrink the alignment; approximate
        # cigar as a single match block of new_len (fine for depth/GC purposes,
        # not meant to preserve indel structure)
        seq = read.query_sequence
        qual = read.query_qualities
        if seq is None:
            continue

        read.query_sequence = seq[:new_len]
        if qual is not None:
            read.query_qualities = qual[:new_len]
        read.cigartuples = [(0, new_len)]  # 0 = BAM_CMATCH

        out.write(read)

    src.close()
    out.close()

    pysam.sort("-o", out_bam_path + ".sorted.bam", out_bam_path)
    pysam.index(out_bam_path + ".sorted.bam")
    print(f"Control BAM written to {out_bam_path}.sorted.bam ({n_reads} reads).")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Build a length-matched Control BAM from a source BAM.")
    p.add_argument("source_bam", help="Modern WGS BAM to sample reads from")
    p.add_argument("target_bam", help="Target (e.g. Saqqaq) BAM to match read-length distribution to")
    p.add_argument("out_bam")
    p.add_argument("--n-reads", type=int, default=None)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    dist = get_read_length_distribution(args.target_bam)
    build_control_bam(args.source_bam, dist, args.out_bam, n_reads=args.n_reads, seed=args.seed)