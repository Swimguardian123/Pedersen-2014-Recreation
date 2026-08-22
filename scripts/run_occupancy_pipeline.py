"""
Driver script: chains depth -> GC correction -> nucleosome calling -> periodicity
-> nucleotide patterns against a real BAM/FASTA, for smoke-testing the occupancy
side of the pipeline end to end.

This is deliberately NOT trying to produce a biologically meaningful result on a
sparse test BAM (10k-read subsets won't have real local depth signal -- see the
math on why in-conversation). The goal here is mechanical: does each stage run
without crashing, and do the shapes/types of what comes out look sane. Read the
printed diagnostics, not the biology, until you're running against real coverage.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from src.depth import get_region_depth
from src.gc_correction import build_gc_rate_model, gc_correct_depth
from src.nucleosome_calling import call_nucleosomes
from src.periodicity import welch_periodogram, dominant_period
from src.nucleotide_patterns import mononucleotide_matrix, purine_pyrimidine_dinucleotide_signal


def run(bam_path: str, fasta_path: str, chrom: str, start: int, end: int) -> None:
    print(f"=== Region: {chrom}:{start}-{end} ({end - start} bp) ===\n")

    print("[1/5] Raw depth (depth.py)")
    raw_depth = get_region_depth(bam_path, chrom, start, end)
    print(f"  shape={raw_depth.shape}, mean={raw_depth.mean():.4f}, "
          f"max={raw_depth.max()}, nonzero_bp={(raw_depth > 0).sum()}\n")

    print("[2/5] GC-corrected depth (gc_correction.py)")
    gc_model = build_gc_rate_model(bam_path, fasta_path, chrom, start, end)
    print(f"  fitted rate model for {len(gc_model)} distinct read length(s)")
    corrected_depth = gc_correct_depth(bam_path, fasta_path, chrom, start, end, gc_model)
    print(f"  shape={corrected_depth.shape}, mean={np.nanmean(corrected_depth):.4f}, "
          f"max={np.nanmax(corrected_depth):.4f}\n")

    print("[3/5] Nucleosome calling (nucleosome_calling.py)")
    calls = call_nucleosomes(corrected_depth)
    print(f"  {len(calls)} nucleosome call(s)")
    if calls:
        top = sorted(calls, key=lambda c: c.score, reverse=True)[0]
        print(f"  top call: center_pos={top.center_pos} (abs genome pos "
              f"{start + top.center_pos}), score={top.score:.2f}, "
              f"peak_depth={top.peak_depth:.2f}")
    print()

    print("[4/5] Periodicity (periodicity.py)")
    if (raw_depth > 0).sum() > 10:
        period = dominant_period(raw_depth.astype(np.float64))
        print(f"  dominant period in [50,300]bp window: {period:.1f} bp "
              f"(paper's TSS peak: 193bp -- expect noise, not a match, on sparse data)")
    else:
        print("  skipped: fewer than 10 nonzero positions, periodogram would be meaningless")
    print()

    print("[5/5] Nucleotide patterns (nucleotide_patterns.py)")
    if calls:
        mono = mononucleotide_matrix(fasta_path, chrom, calls, genome_offset=start, halfwidth=50)
        rr = purine_pyrimidine_dinucleotide_signal(fasta_path, chrom, calls, genome_offset=start, halfwidth=50)
        print(f"  mononucleotide matrix shape={mono.shape}, "
              f"purine-purine signal shape={rr.shape}")
        print(f"  (with only {len(calls)} call(s), expect noisy/NaN-heavy output -- "
              f"this checks the code runs, not that the signal is real yet)")
    else:
        print("  skipped: no nucleosome calls to align nucleotide patterns to")
    print()

    print("=== Done. All stages that could run, ran. Check shapes/counts above for sanity. ===")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Smoke-test the occupancy pipeline against a real BAM.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    args = p.parse_args()

    run(args.bam, args.fasta, args.chrom, args.start, args.end)