"""
Driver script: download the real hg18 mappability track and use it to filter
real nucleosome calls -- smoke-tests mappability_filter.py end to end.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.mappability_filter import (
    download_mappability_track, compute_block_mappability, filter_calls_by_mappability,
    AVAILABLE_KMERS,
)
from src.gc_correction import build_gc_rate_model, gc_correct_depth
from src.nucleosome_calling import call_nucleosomes


def run(bam_path: str, fasta_path: str, chrom: str, start: int, end: int,
        mappability_destdir: str = "./data/mappability", threshold: float = 0.9,
        kmer: int = 40) -> None:
    print(f"=== Region: {chrom}:{start}-{end} ({end - start} bp) ===\n")

    print("[1/3] Calling nucleosomes (unfiltered)")
    gc_model = build_gc_rate_model(bam_path, fasta_path, chrom, start, end)
    corrected_depth = gc_correct_depth(bam_path, fasta_path, chrom, start, end, gc_model)
    calls = call_nucleosomes(corrected_depth)
    print(f"  {len(calls)} call(s) before mappability filtering\n")

    print(f"[2/3] Mappability track (real hg18 CRG {kmer}mer, see AVAILABLE_KMERS for sizes)")
    track_path = download_mappability_track(mappability_destdir, kmer=kmer)
    blocks = compute_block_mappability(track_path, chrom, start, end)
    print("  Per-block mean mappability (raw values, not just pass/fail):")
    for b_start, b_end, mean_val in blocks:
        if mean_val is None:
            print(f"    {b_start}-{b_end}: no data in track (assembly gap or track boundary)")
        else:
            status = "PASS" if mean_val >= threshold else "fail"
            print(f"    {b_start}-{b_end}: {mean_val:.4f}  [{status}]")
    regions = [(s, e) for s, e, m in blocks if m is not None and m >= threshold]
    print(f"  {len(regions)}/{len(blocks)} block(s) pass mappability >= {threshold}\n")

    print("[3/3] Filtering calls")
    filtered = filter_calls_by_mappability(calls, start, regions)
    print(f"  {len(filtered)}/{len(calls)} call(s) survive mappability filtering")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Smoke-test mappability filtering against real data.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--mappability-destdir", default="./data/mappability")
    p.add_argument("--kmer", type=int, default=40, choices=AVAILABLE_KMERS,
                    help="40 (default, 1.1G, closest to paper) down to 100 (94M, smaller download)")
    p.add_argument("--threshold", type=float, default=0.9)
    args = p.parse_args()

    run(args.bam, args.fasta, args.chrom, args.start, args.end,
        args.mappability_destdir, args.threshold, args.kmer)