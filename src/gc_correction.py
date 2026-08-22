"""
GC-content bias correction, following Benjamini & Speed 2012 (NAR 40:e72), the
method Pedersen et al. 2014 cite as "GCcorrect".

Pedersen et al. 2014 Methods: "GCcorrect ... was used to calculate the association
of GC-content with read depth and, in turn, to estimate the expected read depth for
all unique positions over all relevant read lengths in the human genome. The
observed read depth at a given position was normalized for GC content effects by
subtracting the expected read depth, summed across all read lengths."

Benjamini & Speed 2012 Methods (single-position fragment-length model W^s_{a,l}):
  - For each fragment length s, define a GC window of length l = s - a - m, offset
    by margin `a` from the 5' end and margin `m` from the 3' end (small margins are
    trimmed from each end "to reduce the impact of local biases"; the paper doesn't
    give exact a/m values -- default a=m=0 here, override via `margin`).
  - Take a large random sample of *genomic positions* (independent of the sequencing
    data), stratify them by gc = GC(x+a, l), giving N_gc = # sampled positions per
    stratum.
  - Count F_gc = # fragments whose 5' end falls in stratum gc.
  - rate_gc = F_gc / N_gc  <-- this is the key step my first draft got wrong: the
    denominator must be genomic positions, not observed reads, or the model can
    never detect bias (a GC stratum with zero coverage would just be invisible).
  - Sparse strata are pooled and loess-smoothed (smoothness 0.2); we approximate
    this with a simple moving-average smoother over sparse bins (no loess dependency
    required) -- a lighter-weight stand-in, flagged as such.
  - Predicted count per length is summed across lengths (with a scale factor c,
    here fixed at 1 since we're not attempting the unknown-length correction case).
  - Final correction: Benjamini & Speed's own tool divides (observed / predicted).
    Pedersen's paper explicitly says they *subtracted*. Since Pedersen 2014 is the
    reproduction target, default is 'subtract'; pass method='divide' for the
    original GCcorrect behavior.

Caveats vs. the literal method, still open:
  - No mappability track is used to restrict the background sample to mappable
    positions (Benjamini & Speed sample "mappable locations"); we sample uniformly
    across the region instead. Fine for typical unique regions, will understate
    bias correction in low-mappability regions.
  - Margins (a, m) for local-bias trimming default to 0; the paper doesn't specify
    their values.
"""

from collections import defaultdict
from typing import Dict, Tuple

import numpy as np
import pysam


def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    denom = gc + at
    return gc / denom if denom > 0 else np.nan


def _smooth_sparse_bins(rate: np.ndarray, support: np.ndarray, min_support: int = 20,
                         window: int = 5) -> np.ndarray:
    """
    Lightweight stand-in for the paper's loess smoothing of pooled sparse strata:
    bins with fewer than `min_support` background positions get replaced by a
    moving average of their neighbors (weighted by neighbor support), rather than
    left as noisy/undefined single-bin estimates.
    """
    smoothed = rate.copy()
    n = len(rate)
    for i in range(n):
        if support[i] >= min_support:
            continue
        lo, hi = max(0, i - window), min(n, i + window + 1)
        w = support[lo:hi]
        v = rate[lo:hi]
        valid = w > 0
        smoothed[i] = np.average(v[valid], weights=w[valid]) if valid.any() else 0.0
    return smoothed


def build_gc_rate_model(bam_path: str, fasta_path: str, chrom: str, start: int, end: int,
                         n_gc_bins: int = 50, margin: int = 0, n_background: int = 20000,
                         min_mapq: int = 0, seed: int = 0) -> Dict[int, np.ndarray]:
    """
    Fit rate_gc = F_gc / N_gc per fragment length, using a genome-wide (region-wide)
    random background sample of positions as N_gc, and observed 5'-end fragment
    counts as F_gc.

    Returns: {length: rate_array} where rate_array[bin] = fragments per background
    position in that GC bin, for fragments of that length.
    """
    rng = np.random.default_rng(seed)
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)

    # observed fragment 5'-end counts, stratified by (length, gc_bin)
    frag_counts: Dict[int, np.ndarray] = defaultdict(lambda: np.zeros(n_gc_bins))
    lengths_seen = set()

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue
        L = read.query_alignment_length
        if not L or L <= 2 * margin:
            continue
        lengths_seen.add(L)
        w_start = read.reference_start + margin
        w_len = L - 2 * margin
        try:
            seq = fasta.fetch(chrom, w_start, w_start + w_len)
        except (ValueError, IndexError):
            continue
        gc = gc_fraction(seq)
        if np.isnan(gc):
            continue
        b = min(int(gc * n_gc_bins), n_gc_bins - 1)
        frag_counts[L][b] += 1

    # background sample: N random genomic positions in [start, end), independent of
    # the read data, stratified the same way per length
    region_len = end - start
    model: Dict[int, np.ndarray] = {}

    for L in lengths_seen:
        w_len = L - 2 * margin
        if w_len <= 0 or region_len <= w_len:
            continue
        sample_starts = rng.integers(start, end - w_len, size=min(n_background, region_len))
        support = np.zeros(n_gc_bins)
        for s in sample_starts:
            try:
                seq = fasta.fetch(chrom, int(s) + margin, int(s) + margin + w_len)
            except (ValueError, IndexError):
                continue
            gc = gc_fraction(seq)
            if np.isnan(gc):
                continue
            b = min(int(gc * n_gc_bins), n_gc_bins - 1)
            support[b] += 1

        with np.errstate(divide="ignore", invalid="ignore"):
            rate = np.where(support > 0, frag_counts[L] / support, 0.0)
        model[L] = _smooth_sparse_bins(rate, support)

    bam.close()
    fasta.close()
    return model


def gc_correct_depth(bam_path: str, fasta_path: str, chrom: str, start: int, end: int,
                      model: Dict[int, np.ndarray], n_gc_bins: int = 50, margin: int = 0,
                      min_mapq: int = 0, method: str = "subtract",
                      divide_pseudocount: float = 0.1) -> np.ndarray:
    """
    Apply the fitted rate model to compute GC-corrected depth.

    method='subtract' (default, matches Pedersen et al. 2014's stated approach):
        corrected = observed - expected
    method='divide' (matches Benjamini & Speed's own GCcorrect tool):
        corrected = observed / (expected + divide_pseudocount)
    """
    try:
        from .depth import get_region_depth
    except ImportError:
        from depth import get_region_depth

    observed = get_region_depth(bam_path, chrom, start, end, min_mapq=min_mapq).astype(np.float64)
    expected = np.zeros_like(observed)

    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue
        L = read.query_alignment_length
        if not L or L not in model or L <= 2 * margin:
            continue
        w_start = read.reference_start + margin
        w_len = L - 2 * margin
        try:
            seq = fasta.fetch(chrom, w_start, w_start + w_len)
        except (ValueError, IndexError):
            continue
        gc = gc_fraction(seq)
        if np.isnan(gc):
            continue
        b = min(int(gc * n_gc_bins), n_gc_bins - 1)
        expected_val = model[L][b]

        rel_start = max(read.reference_start - start, 0)
        rel_end = min(read.reference_start + L - start, end - start)
        if rel_end > rel_start:
            expected[rel_start:rel_end] += expected_val

    bam.close()
    fasta.close()

    if method == "subtract":
        return observed - expected
    elif method == "divide":
        return observed / (expected + divide_pseudocount)
    else:
        raise ValueError("method must be 'subtract' or 'divide'")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="GC-correct read depth over a region.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--method", choices=["subtract", "divide"], default="subtract")
    p.add_argument("--out", default=None)
    args = p.parse_args()

    m = build_gc_rate_model(args.bam, args.fasta, args.chrom, args.start, args.end)
    corrected = gc_correct_depth(args.bam, args.fasta, args.chrom, args.start, args.end,
                                  m, method=args.method)

    if args.out:
        np.save(args.out, corrected)
        print(f"Saved GC-corrected depth ({len(corrected)} bp, method={args.method}) to {args.out}")
    else:
        print(corrected)