"""
Windowed Protection Score (WPS) -- Snyder et al. 2016 (Cell 164:57-68).

A genuinely different nucleosome-calling algorithm from NucleoMap
(nucleosome_calling.py), included here as a real comparison method, per
Hanghoj et al. 2016's own benchmark (epiPALEOMIX vs. WPS on the chr12
alpha-satellite array). Not a Pedersen 2014 method -- flagged clearly as a
different paper's algorithm, not a refinement of our existing pipeline.

Formula (Snyder et al. 2016, Results, confirmed across multiple independent
citations of the paper, not just one figure):
"a windowed protection score (WPS), which is the number of DNA fragments
completely spanning a [W] bp window centered at a given genomic coordinate
minus the number of fragments with an endpoint within that same window."

Hanghoj et al. 2016 swept window sizes W in {10, 20, 30, 50, 80, 120} bp when
comparing WPS against NucleoMap on the chr12 alpha-satellite array; their own
finding was that NucleoMap outperformed WPS for aDNA specifically (WPS called
fewer peaks than the true nucleosome array, and NucleoMap's calls sat closer
to the true MNase-derived dyad positions) -- this module exists for that
comparison, not because WPS is expected to outperform our existing pipeline
on ancient data.

Snyder et al. 2016 additionally distinguish "L-WPS" (long-fragment WPS, using
120bp windows and typically 120-180bp fragments, targeting nucleosome-scale
protection) from "S-WPS" (short-fragment WPS, much smaller windows, targeting
sub-nucleosomal/transcription-factor footprints in modern cfDNA). Ancient DNA
fragments are already naturally short, so Hanghoj et al. 2016's own comparison
didn't need this L/S split -- this module implements the general single-window
WPS Hanghoj actually used, with fragment-length filtering exposed as an
optional parameter if you want to reproduce the modern-cfDNA L/S distinction.
"""

from typing import List, Optional, Tuple

import numpy as np
import pysam


def wps_from_fragments(fragments: List[Tuple[int, int]], start: int, end: int,
                        window: int = 120) -> np.ndarray:
    """
    Pure WPS counting logic, decoupled from BAM fetching -- takes a plain list
    of (frag_start, frag_end) tuples. Kept separate from compute_wps() so the
    actual counting rule (spanning vs. endpoint-in-window) is independently
    testable without needing pysam or a real BAM file.
    """
    half_w = window // 2
    n = end - start
    wps = np.zeros(n, dtype=np.float64)

    for pos_idx in range(n):
        pos = start + pos_idx
        win_start = pos - half_w
        win_end = pos + half_w

        spanning = 0
        endpoint_in_window = 0
        for frag_start, frag_end in fragments:
            if frag_start <= win_start and frag_end >= win_end:
                spanning += 1
            elif win_start <= frag_start <= win_end or win_start <= frag_end <= win_end:
                endpoint_in_window += 1

        wps[pos_idx] = spanning - endpoint_in_window

    return wps


def compute_wps(bam_path: str, chrom: str, start: int, end: int,
                 window: int = 120, min_fragment_len: Optional[int] = None,
                 max_fragment_len: Optional[int] = None,
                 min_mapq: int = 0) -> np.ndarray:
    """
    Per-base WPS across [start, end).

    window: protection-window width in bp (Snyder et al. 2016 default: 120bp
        for nucleosome-scale L-WPS; Hanghoj et al. 2016 swept 10-120bp).
    min_fragment_len / max_fragment_len: optional fragment-length filter, for
        reproducing Snyder's L-WPS vs S-WPS distinction on modern data. Left
        as None (no filtering) by default, appropriate for naturally-short
        aDNA fragments.

    For single-end aDNA reads (no true fragment/insert-size info the way
    paired-end cfDNA has), the read's own aligned span (reference_start to
    reference_end) is used as the fragment -- reasonable for aDNA where the
    sequenced read typically corresponds to the whole (already short)
    molecule, but flagged explicitly since Snyder's original method assumes
    real paired-end fragment boundaries.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # widen the fetch region so fragments that partially overlap our region's
    # edges are still considered for windows near those edges
    fetch_start = max(start - window, 0)
    fetch_end = end + window

    fragments = []
    for read in bam.fetch(chrom, fetch_start, fetch_end):
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue
        frag_start = read.reference_start
        frag_end = read.reference_end
        if frag_end is None:
            continue
        frag_len = frag_end - frag_start
        if min_fragment_len is not None and frag_len < min_fragment_len:
            continue
        if max_fragment_len is not None and frag_len > max_fragment_len:
            continue
        fragments.append((frag_start, frag_end))

    bam.close()

    return wps_from_fragments(fragments, start, end, window)


def call_peaks_from_wps(wps: np.ndarray, min_peak_width: int = 50) -> List[Tuple[int, int]]:
    """
    Simple peak-region caller on a WPS track: contiguous runs of positive WPS
    (protected regions) at least min_peak_width long. Snyder et al. 2016 use a
    more elaborate heuristic for their genome-wide 10M+ peak set; this is a
    straightforward stand-in sufficient for comparing against NucleoMap calls
    on a single region (the actual comparison Hanghoj et al. 2016 ran), not a
    claim of reproducing Snyder's exact production peak-caller.
    """
    peaks = []
    in_peak = False
    peak_start = None

    for i, v in enumerate(wps):
        if v > 0 and not in_peak:
            in_peak = True
            peak_start = i
        elif v <= 0 and in_peak:
            in_peak = False
            if i - peak_start >= min_peak_width:
                peaks.append((peak_start, i))

    if in_peak and len(wps) - peak_start >= min_peak_width:
        peaks.append((peak_start, len(wps)))

    return peaks


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Compute Windowed Protection Score (Snyder et al. 2016).")
    p.add_argument("bam")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--window", type=int, default=120)
    p.add_argument("--min-fragment-len", type=int, default=None)
    p.add_argument("--max-fragment-len", type=int, default=None)
    args = p.parse_args()

    wps = compute_wps(args.bam, args.chrom, args.start, args.end,
                       window=args.window, min_fragment_len=args.min_fragment_len,
                       max_fragment_len=args.max_fragment_len)
    peaks = call_peaks_from_wps(wps)
    print(f"WPS computed over {len(wps)} positions (window={args.window}bp)")
    print(f"{len(peaks)} peak region(s) found")