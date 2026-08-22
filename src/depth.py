"""
Read depth extraction from a BAM file.

Paper reference (Methods, "Sequence data sets" / Results "Nucleosome occupancy signal"):
read depth is computed per-base from mapped, deduplicated ("collapsed") reads.
The paper's own data had PCR-clone collapsing done upstream (Sesam pipeline, see
Rasmussen 2010 SI); if your BAM is not already deduplicated, dedupe it first
(e.g. `samtools markdup -r` or your aligner's collapsing step) before calling this.
"""

import numpy as np
import pysam


def get_region_depth(bam_path: str, chrom: str, start: int, end: int,
                      min_mapq: int = 0) -> np.ndarray:
    """
    Per-base read depth over [start, end) (0-based, half-open), matching
    samtools coordinate convention.

    Returns an array of length (end - start).
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    depth = np.zeros(end - start, dtype=np.int32)

    # pysam.count_coverage sums per-base per-nucleotide coverage; summing across
    # the 4 bases gives total depth. quality_threshold=0 avoids silently dropping
    # low-quality bases unless you want that (aDNA reads are often short/lower qual).
    a, c, g, t = bam.count_coverage(
        chrom, start, end, quality_threshold=0, read_callback=lambda r: r.mapping_quality >= min_mapq
    )
    depth = np.array(a) + np.array(c) + np.array(g) + np.array(t)
    bam.close()
    return depth


def get_chrom_depth(bam_path: str, chrom: str, min_mapq: int = 0) -> np.ndarray:
    """Convenience wrapper: full-chromosome depth track."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    length = bam.get_reference_length(chrom)
    bam.close()
    return get_region_depth(bam_path, chrom, 0, length, min_mapq=min_mapq)


def genomic_average_depth(bam_path: str, chroms=None) -> float:
    """
    Genomic average (GA) depth, used throughout the paper as the normalization
    baseline (e.g. "0.9x GA in intergenic regions", "6.5x GA at CGIs").
    Computed as total mapped bases / total reference length considered.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    if chroms is None:
        chroms = list(bam.references)
    total_len = sum(bam.get_reference_length(c) for c in chroms)
    bam.close()

    total_depth = 0
    for c in chroms:
        total_depth += get_chrom_depth(bam_path, c).sum()

    return total_depth / total_len


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Extract per-base read depth from a BAM region.")
    p.add_argument("bam")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--out", default=None, help="Optional .npy output path")
    args = p.parse_args()

    d = get_region_depth(args.bam, args.chrom, args.start, args.end)
    if args.out:
        np.save(args.out, d)
        print(f"Saved depth array of length {len(d)} to {args.out}")
    else:
        print(d)