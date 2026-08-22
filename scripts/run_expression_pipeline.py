"""
Driver script: smoke-test expression_inference.py's proxies.

Rs and phasing_strength/first_nucleosome_occupancy are testable against real BAM
depth right now; the GSE3058 fetch is testable too but needs real internet access
that this environment doesn't reliably have when built for you -- run --fetch-gse3058
yourself to confirm that half independently.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from src.expression_inference import (
    compute_rs, first_nucleosome_occupancy, phasing_strength,
    fetch_gse3058_expression, rank_into_quantiles,
)
from src.gc_correction import build_gc_rate_model, gc_correct_depth
from src.methylation import ms_score


def run(bam_path: str, fasta_path: str, chrom: str, tss_pos: int,
        fetch_gse3058: bool = False) -> None:
    print(f"=== TSS test position: {chrom}:{tss_pos} ===\n")

    print("[1/4] Rs (gene body / promoter Ms ratio)")
    promoter_ms, _, promoter_total = ms_score(bam_path, fasta_path, chrom,
                                               tss_pos - 500, tss_pos)
    genebody_ms, _, genebody_total = ms_score(bam_path, fasta_path, chrom,
                                               tss_pos, tss_pos + 2000)
    rs = compute_rs(genebody_ms if genebody_total else 0.0,
                     promoter_ms if promoter_total else 0.0)
    print(f"  promoter Ms={promoter_ms} (n={promoter_total}), "
          f"gene body Ms={genebody_ms} (n={genebody_total})")
    print(f"  Rs={rs:.4f}\n")

    print("[2/4] +1 nucleosome occupancy")
    model = build_gc_rate_model(bam_path, fasta_path, chrom, tss_pos, tss_pos + 300)
    corrected = gc_correct_depth(bam_path, fasta_path, chrom, tss_pos, tss_pos + 300, model)
    occ = first_nucleosome_occupancy(corrected)
    print(f"  +1 nucleosome occupancy={occ:.4f}\n")

    print("[3/4] Phasing strength")
    window_depth = gc_correct_depth(bam_path, fasta_path, chrom, tss_pos - 1000,
                                     tss_pos + 1000, model)
    strength = phasing_strength(window_depth)
    print(f"  phasing strength (power at ~193bp)={strength:.4f}\n")

    print("[4/4] GSE3058 expression data")
    if fetch_gse3058:
        expr = fetch_gse3058_expression()
        print(f"  fetched expression for {len(expr)} probes")
        groups = rank_into_quantiles(expr, n_quantiles=10)
        print(f"  assigned to 10 quantile groups")
    else:
        print("  skipped (pass --fetch-gse3058 to actually pull it -- needs your own internet)")

    print(f"\n  (all three proxies ran against real data -- values are effectively "
          f"noise at this coverage, this confirms mechanics not biology.)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Smoke-test expression inference proxies against real data.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("tss_pos", type=int, help="A test TSS-like position (any real covered position works)")
    p.add_argument("--fetch-gse3058", action="store_true")
    args = p.parse_args()

    run(args.bam, args.fasta, args.chrom, args.tss_pos, args.fetch_gse3058)