"""
Driver script: run WPS and NucleoMap (call_nucleosomes) on the same real
region and compare -- the actual comparison Hanghoj et al. 2016 ran (chr12
alpha-satellite array), reproduced here at whatever scale your real data
supports. Smoke-tests wps.py's BAM-fetching path (compute_wps), which has
never been run against real data before now, same as every other core module.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from src.wps import compute_wps, call_peaks_from_wps
from src.gc_correction import build_gc_rate_model, gc_correct_depth
from src.nucleosome_calling import call_nucleosomes


def run(bam_path: str, fasta_path: str, chrom: str, start: int, end: int,
        window: int = 120, min_peak_width: int = 50) -> None:
    print(f"=== Region: {chrom}:{start}-{end} ({end - start} bp), WPS window={window}bp, "
          f"min_peak_width={min_peak_width}bp ===\n")

    print("[1/2] WPS (Snyder et al. 2016)")
    wps = compute_wps(bam_path, chrom, start, end, window=window)
    print(f"  shape={wps.shape}, mean={wps.mean():.4f}, min={wps.min():.1f}, max={wps.max():.1f}")
    wps_peaks = call_peaks_from_wps(wps, min_peak_width=min_peak_width)
    print(f"  {len(wps_peaks)} WPS peak region(s)")
    wps_centers = []
    if wps_peaks:
        widths = [e - s for s, e in wps_peaks]
        print(f"  peak widths: min={min(widths)}, max={max(widths)}, mean={np.mean(widths):.1f}")
        wps_centers = [start + (s + e) // 2 for s, e in wps_peaks]
        print(f"  peak center(s), absolute genome position: {wps_centers}")
    print()

    print("[2/2] NucleoMap (nucleosome_calling.py) on the same region, for comparison")
    gc_model = build_gc_rate_model(bam_path, fasta_path, chrom, start, end)
    corrected_depth = gc_correct_depth(bam_path, fasta_path, chrom, start, end, gc_model)
    nucleomap_calls = call_nucleosomes(corrected_depth)
    print(f"  {len(nucleomap_calls)} NucleoMap call(s)")
    nucleomap_centers = [start + c.center_pos for c in nucleomap_calls]
    if nucleomap_calls:
        print(f"  call center(s), absolute genome position: {nucleomap_centers}")
    print()

    print(f"Comparison: WPS found {len(wps_peaks)} peak region(s), NucleoMap found "
          f"{len(nucleomap_calls)} call(s).")
    if wps_centers and nucleomap_centers:
        print("\n  Position agreement check (this, not the raw counts, is the meaningful")
        print("  comparison -- matches Hanghoj et al. 2016's own approach of checking")
        print("  physical proximity between methods' called dyads, not just counting):")
        for wc in wps_centers:
            nearest = min(nucleomap_centers, key=lambda nc: abs(nc - wc))
            print(f"    WPS peak at {wc} -> nearest NucleoMap call at {nearest} "
                  f"(distance: {abs(nearest - wc)}bp)")
    print("(Hanghoj et al. 2016's own finding: NucleoMap outperforms WPS on ancient data --")
    print(" this isn't expected to show WPS winning, the comparison itself is the point.)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Compare WPS against NucleoMap on real data.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--window", type=int, default=120)
    p.add_argument("--min-peak-width", type=int, default=50,
                    help="Lower this (e.g. --min-peak-width 5) to test peak-calling "
                         "mechanics on sparse data, separate from whether there's "
                         "enough coverage for a real nucleosome-scale peak")
    args = p.parse_args()

    run(args.bam, args.fasta, args.chrom, args.start, args.end, args.window, args.min_peak_width)