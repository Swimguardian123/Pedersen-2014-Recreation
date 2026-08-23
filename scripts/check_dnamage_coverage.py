"""
Step 3 (the decisive check): for each of the 353 Horvath clock CpGs (now in
hg18, after liftover), check whether our real test BAM has ANY read coverage
within a 2kb window centered on that position -- per Hanghoj et al. 2016's
own stated method (regional Ms score in 2kb windows around each clock CpG).

This is the actual feasibility answer: if most/all windows come back with
zero coverage, building the full Ms-to-beta pipeline further isn't worth it
at this data scale, same reasoning as the full-genome-coverage decision
earlier in this project. If a meaningful number DO have coverage, that
changes the calculus.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from src.depth import get_region_depth


def load_hg18_bed(bed_path: str):
    """Yields (chrom, pos, cpg_id) for each line of a BED file (post-liftOver)."""
    with open(bed_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            chrom, start, end, cpg_id = fields[0], int(fields[1]), int(fields[2]), fields[3]
            yield chrom, start, cpg_id


def check_coverage(bam_path: str, hg18_bed_path: str, window: int = 2000):
    """
    For each CpG, checks depth in a `window`-bp region centered on it.
    Returns (results, n_with_coverage) where results is a list of
    (cpg_id, chrom, pos, has_coverage, max_depth).
    """
    half = window // 2
    results = []
    n_with_coverage = 0

    for chrom, pos, cpg_id in load_hg18_bed(hg18_bed_path):
        region_start = max(pos - half, 0)
        region_end = pos + half
        try:
            depth = get_region_depth(bam_path, chrom, region_start, region_end)
            has_coverage = bool(depth.max() > 0)
            max_depth = int(depth.max())
        except Exception as e:
            has_coverage = False
            max_depth = 0
        if has_coverage:
            n_with_coverage += 1
        results.append((cpg_id, chrom, pos, has_coverage, max_depth))

    return results, n_with_coverage


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Check real coverage at the 353 Horvath clock CpGs (hg18).")
    p.add_argument("bam_path")
    p.add_argument("hg18_bed", help="post-liftOver BED of the 353 clock CpGs, hg18 coordinates")
    p.add_argument("--window", type=int, default=2000, help="window width, bp (Hanghoj 2016: 2kb)")
    args = p.parse_args()

    results, n_covered = check_coverage(args.bam_path, args.hg18_bed, args.window)
    total = len(results)

    print(f"=== Coverage check: {total} clock CpGs, {args.window}bp windows ===\n")
    covered = [r for r in results if r[3]]
    for cpg_id, chrom, pos, has_cov, max_depth in covered:
        print(f"  {cpg_id}  {chrom}:{pos}  max_depth={max_depth}")

    print(f"\n{n_covered}/{total} clock CpGs ({100*n_covered/total:.1f}%) have real coverage "
          f"in a {args.window}bp window, in this test BAM.")

    if n_covered == 0:
        print("\nZero coverage: at this data scale, no basis for a real DNAmAge estimate.")
        print("Recommend treating this as a documented limitation, same as full-genome coverage.")
    elif n_covered < total * 0.05:
        print(f"\n{n_covered} of 353 is a small fraction -- any resulting age estimate would be")
        print("built on well under 5% of the clock's intended CpG panel. Worth discussing before")
        print("investing further in the Ms-to-beta conversion step.")
    else:
        print(f"\n{n_covered} covered CpGs is enough to be worth pursuing further.")