"""
Driver script: smoke-test ctcf_analysis.py's Ms and depth anchor profiles against
real BAM/FASTA data (distinct from your own run_ctcf.py -- this is the automated
smoke-test counterpart to the other run_*_pipeline.py scripts).

Filters the Fu et al. 2008 site table to a single chromosome by default, since a
small test BAM/reference (like our chr12-only setup) has no way to find signal at
sites on other chromosomes.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pysam

from src.ctcf_analysis import parse_fu2008_sites, ctcf_ms_profile, ctcf_depth_profile
from src.gc_correction import build_gc_rate_model, gc_correct_depth


def run(sites_tsv: str, bam_path: str, fasta_path: str, chrom_filter: str = None,
        half_window: int = 1000, bin_size: int = 25) -> None:
    sites = parse_fu2008_sites(sites_tsv, which="occupied")
    print(f"[1/3] Loaded {len(sites)} occupied CTCF sites total")

    if chrom_filter:
        sites = [s for s in sites if s[0] == chrom_filter]
        print(f"  filtered to {len(sites)} site(s) on {chrom_filter}\n")

    if not sites:
        print("No sites to test against after filtering -- nothing further to run.")
        return

    print("[2/3] Ms profile around CTCF sites (methylation.py + anchor_profile)")
    ms_profile = ctcf_ms_profile(bam_path, fasta_path, sites, half_window, bin_size)
    valid = (~np.isnan(ms_profile)).sum()
    print(f"  profile shape={ms_profile.shape}, valid (non-NaN) bins={valid}")
    if valid:
        print(f"  first few values: {ms_profile[:5]}")
    print()

    print("[3/3] GC-corrected depth profile around CTCF sites")
    bam = pysam.AlignmentFile(bam_path, "rb")
    chrom_len = bam.get_reference_length(chrom_filter) if chrom_filter else None
    bam.close()

    if chrom_filter is None:
        print("  skipped: pass --chrom to fit one GC model over that chromosome "
              "(fitting genome-wide is possible but slow/unnecessary for a smoke test)")
        return

    gc_model = build_gc_rate_model(bam_path, fasta_path, chrom_filter, 0, chrom_len)

    def depth_fn(chrom, start, end):
        return gc_correct_depth(bam_path, fasta_path, chrom, max(start, 0),
                                 min(end, chrom_len), gc_model)

    depth_profile = ctcf_depth_profile(depth_fn, sites, half_window, bin_size)
    valid_d = (~np.isnan(depth_profile)).sum()
    print(f"  profile shape={depth_profile.shape}, valid (non-NaN) bins={valid_d}")

    print(f"\n  (at this coverage, most bins will be NaN or near-zero -- this confirms "
          f"the site-parsing, strand-orientation, and anchor-averaging logic runs "
          f"correctly against real data, not that the profile shows real CTCF signal.)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Smoke-test the CTCF anchor pipeline against real data.")
    p.add_argument("sites_tsv")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("--chrom", default=None, help="Filter sites to one chromosome (recommended for small test data)")
    p.add_argument("--half-window", type=int, default=1000)
    p.add_argument("--bin-size", type=int, default=25)
    args = p.parse_args()

    run(args.sites_tsv, args.bam, args.fasta, args.chrom, args.half_window, args.bin_size)