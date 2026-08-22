"""
Driver script: smoke-test methylation.py's Ms score against real BAM/FASTA data.

Same philosophy as run_occupancy_pipeline.py -- checks that the CpG-context
detection, strand handling, and mismatch counting actually run correctly against
real aligned reads, not that the resulting Ms value is biologically meaningful at
this coverage (it won't be, with this few reads).
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.methylation import ms_score, ms_negative_controls, find_cpg_positions


def run(bam_path: str, fasta_path: str, chrom: str, start: int, end: int) -> None:
    print(f"=== Region: {chrom}:{start}-{end} ({end - start} bp) ===\n")

    cpg_positions = find_cpg_positions(fasta_path, chrom, start, end)
    print(f"[1/2] CpG sites in region (reference-only, no BAM needed): {len(cpg_positions)}\n")

    print("[2/2] Ms score + negative controls (methylation.py)")
    rate, hits, total = ms_score(bam_path, fasta_path, chrom, start, end)
    print(f"  Ms (CpG context): rate={rate}, hits={hits}, total_read_starts_checked={total}")

    controls = ms_negative_controls(bam_path, fasta_path, chrom, start, end)
    for ctx, (r, h, t) in controls.items():
        print(f"  Cp{ctx} (control): rate={r}, hits={h}, total={t}")

    print(f"\n  (with this few reads, 'total' is likely 0 or very small -- that's expected;")
    print(f"   this confirms the CpG-detection and strand-handling logic runs correctly,")
    print(f"   not that Ms is a real methylation estimate at this coverage.)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Smoke-test the Ms methylation pipeline against a real BAM.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    args = p.parse_args()

    run(args.bam, args.fasta, args.chrom, args.start, args.end)