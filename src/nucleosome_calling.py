"""
Nucleosome calling from GC-corrected read depth.

Paper reference (Methods, "Nucleosome calls"):
"Nucleosomes were called using a sliding window of 147 bp over regions with a
positive GC-corrected read depth. A given position was called as the center of a
nucleosome when showing the highest GC-corrected read depth within a 147-bp window
centered on that position. A score incorporating both occupancy and positioning was
then calculated for each nucleosome defined as the read depth over the peak
(occupancy) minus the mean read depth of the flanking regions (positioning)."

And (Results, Fig 3A caption):
"score = p - (lf + rf) / 2"
where p = peak depth, lf/rf = mean depth of left/right flanking regions.

Flank width is not explicitly stated as a number in the main text (only "flanking
regions" / "linker regions", elsewhere characterized as ~50 bp given 147 bp
nucleosome core + ~200 bp total periodicity => ~50 bp linker). We default flank
width to 50 bp on each side as the most literature-consistent choice; override via
--flank if you have supplement-derived confirmation.

This is a fully exact reimplementation of the described algorithm modulo that one
unstated parameter (flank width).
"""

from dataclasses import dataclass
from typing import List

import numpy as np


NUCLEOSOME_WIDTH = 147  # bp, stated explicitly in the paper throughout


@dataclass
class NucleosomeCall:
    center_pos: int  # 0-based genomic position, relative to the depth array's start
    peak_depth: float
    score: float


def call_nucleosomes(depth: np.ndarray, flank: int = 50,
                      nuc_width: int = NUCLEOSOME_WIDTH) -> List[NucleosomeCall]:
    """
    depth: 1D array of GC-corrected read depth over a contiguous region.
    flank: bp on each side of the nucleosome window used as "linker" for scoring.

    A position i is called a nucleosome center if depth[i] is the max value within
    the window [i - half, i + half] (147 bp wide, so half = 73 on the smaller side).
    Only positions with positive depth are eligible (per "regions with a positive
    GC-corrected read depth").
    """
    half = nuc_width // 2  # 73
    n = len(depth)
    calls = []

    for i in range(half, n - half):
        if depth[i] <= 0:
            continue
        window = depth[i - half:i + half + 1]
        if depth[i] != window.max():
            continue
        # tie-breaking: if multiple positions in the window share the max, only
        # call the first occurrence to avoid duplicate/adjacent calls at plateaus
        if np.argmax(window) != half:
            continue

        left_start = max(i - half - flank, 0)
        left_flank = depth[left_start:i - half]
        right_end = min(i + half + 1 + flank, n)
        right_flank = depth[i + half + 1:right_end]

        if len(left_flank) == 0 or len(right_flank) == 0:
            continue  # can't score at the very edge of the region

        lf = left_flank.mean()
        rf = right_flank.mean()
        score = depth[i] - (lf + rf) / 2

        calls.append(NucleosomeCall(center_pos=i, peak_depth=float(depth[i]), score=float(score)))

    return calls


def score_threshold_for_fdr(saqqaq_scores: np.ndarray, control_scores: np.ndarray,
                             target_fdr: float = 0.05) -> float:
    """
    Paper reference (Results): "We found that the 25% top-scoring calls (2.66M;
    score >= 22.5) ... have an FDR of 5%" and "Repeating the same procedure over the
    Control set enabled us to calculate the FDR for any given score threshold"
    (Methods).

    FDR(threshold) is estimated by treating every Control call above `threshold` as
    a false positive: FDR = n_control_calls_above_t / n_saqqaq_calls_above_t.
    Returns the lowest score threshold achieving FDR <= target_fdr.
    """
    thresholds = np.unique(np.concatenate([saqqaq_scores, control_scores]))
    thresholds.sort()

    best_t = None
    for t in thresholds:
        n_s = (saqqaq_scores >= t).sum()
        n_c = (control_scores >= t).sum()
        if n_s == 0:
            continue
        fdr = n_c / n_s
        if fdr <= target_fdr:
            best_t = t
            break  # thresholds sorted ascending -> first hit is the lowest qualifying t

    return best_t


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Call nucleosomes from a GC-corrected depth .npy array.")
    p.add_argument("depth_npy")
    p.add_argument("--flank", type=int, default=50)
    p.add_argument("--out", default=None, help="Optional CSV output path")
    args = p.parse_args()

    depth = np.load(args.depth_npy)
    calls = call_nucleosomes(depth, flank=args.flank)
    print(f"Called {len(calls)} nucleosomes.")

    if args.out:
        import csv
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["center_pos", "peak_depth", "score"])
            for c in calls:
                w.writerow([c.center_pos, c.peak_depth, c.score])
        print(f"Wrote calls to {args.out}")