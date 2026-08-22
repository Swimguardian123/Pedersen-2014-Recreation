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

UPDATE (flank width + offset now confirmed, not guessed): Pedersen 2014's own
text doesn't give the flank width as a number. Hanghøj et al. 2016 (Mol Biol
Evol, epiPALEOMIX -- an independent reimplementation of these exact methods,
citing Pedersen 2014 directly) states it precisely in Fig. 1B's caption: "the
nucleosome score is calculated as the coverage at the center minus the mean
read depth of the two 25-bp flanking regions... defined with a 12-bp offset
from 147-bp nucleosome window coordinates." That is: flanks are 25bp wide (not
our earlier guess of 50bp), and are NOT adjacent to the 147bp window -- there's
a 12bp gap between the window edge and where each flank region starts. Layout,
left to right: [25bp flank][12bp gap][147bp window][12bp gap][25bp flank].

This is a real correction, not just added precision: the earlier default (50bp
flank, no gap) used a different, larger, and immediately-adjacent region --
flagged here as history, not hidden.
"""

from dataclasses import dataclass
from typing import List

import numpy as np


NUCLEOSOME_WIDTH = 147  # bp, stated explicitly in the paper throughout
FLANK_WIDTH = 25        # bp, from Hanghøj et al. 2016 Fig. 1B
FLANK_OFFSET = 12       # bp gap between window edge and flank start, same source


@dataclass
class NucleosomeCall:
    center_pos: int  # 0-based genomic position, relative to the depth array's start
    peak_depth: float
    score: float


def call_nucleosomes(depth: np.ndarray, flank: int = FLANK_WIDTH,
                      offset: int = FLANK_OFFSET,
                      nuc_width: int = NUCLEOSOME_WIDTH) -> List[NucleosomeCall]:
    """
    depth: 1D array of GC-corrected read depth over a contiguous region.
    flank: bp width of each flanking region used for scoring (default 25bp,
        Hanghøj et al. 2016 Fig. 1B).
    offset: bp gap between the 147bp window's edge and where each flank region
        starts (default 12bp, same source). Set to 0 to reproduce the earlier
        adjacent-flank behavior if you need to compare against it.

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

        left_flank_end = max(i - half - offset, 0)
        left_flank_start = max(left_flank_end - flank, 0)
        left_flank = depth[left_flank_start:left_flank_end]

        right_flank_start = min(i + half + 1 + offset, n)
        right_flank_end = min(right_flank_start + flank, n)
        right_flank = depth[right_flank_start:right_flank_end]

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
    p.add_argument("--flank", type=int, default=FLANK_WIDTH)
    p.add_argument("--offset", type=int, default=FLANK_OFFSET)
    p.add_argument("--out", default=None, help="Optional CSV output path")
    args = p.parse_args()

    depth = np.load(args.depth_npy)
    calls = call_nucleosomes(depth, flank=args.flank, offset=args.offset)
    print(f"Called {len(calls)} nucleosomes.")

    if args.out:
        import csv
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["center_pos", "peak_depth", "score"])
            for c in calls:
                w.writerow([c.center_pos, c.peak_depth, c.score])
        print(f"Wrote calls to {args.out}")