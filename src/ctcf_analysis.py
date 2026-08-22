"""
CTCF anchor analysis: Ms methylation and GC-corrected occupancy around CTCF sites.

Paper reference (Methods, "Methylation signal", final paragraph):
"The CTCF nucleosomes patterns are derived from a set of 12,864 published CTCF
binding sites (Fu et al. 2008). We calculated (1) the Ms score, and (2) the
nucleosome occupancy (GC-corrected read depth) of 25-bp sliding windows within 1 kb
of the sites."

Site data: Fu et al. 2008 (PLoS Genetics 4:e1000138) Table S13 lists paired
occupied/unoccupied CTCF motif sites as "chrN:start-end" + strand, hg18 coordinates.
Expected input here is that table exported to tab-delimited text with columns:
    Occupied_site  Occupied_site_strand  Unoccupied_site  Unoccupied_site_strand
(export the table from the supplementary .doc as plain/tab-delimited text -- the
doc's table didn't survive automated .doc->.txt conversion cleanly, so this is a
manual one-time step on your end).

Strand handling: sites are oriented so that profiles align in a consistent 5'->3'
direction relative to the CTCF motif -- for '-' strand sites, the bin order is
reversed before averaging, matching how Fig. 4D's symmetric anchor plot is built.
"""

import re
from typing import Callable, List, Tuple

import numpy as np

SiteTuple = Tuple[str, int, str]  # (chrom, center_pos, strand)

_COORD_RE = re.compile(r"^(chr\w+):(\d+)-(\d+)$")

def _parse_coord(s: str) -> Tuple[str, int, int]:
    m = _COORD_RE.match(s.strip())
    if not m:
        raise ValueError(f"Unrecognized coordinate format: {s!r}")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    return chrom, start, end


def parse_fu2008_sites(path: str, which: str = "occupied",
                        skip_random: bool = True) -> List[SiteTuple]:
    """
    Parse the Fu et al. 2008 occupied/unoccupied CTCF site table (tab-delimited,
    4 columns: Occupied_site, Occupied_site_strand, Unoccupied_site,
    Unoccupied_site_strand -- header line starting with '#' is skipped).

    which: "occupied" (the actual CTCF-bound sites Pedersen et al. used) or
        "unoccupied" (Fu's matched negative-control sites -- useful if you want
        your own FDR-style comparison, analogous to control_set.py's role
        elsewhere in this pipeline).
    skip_random: drop chrN_random / chrN_hap contigs (not part of the standard
        assembly, will fail to fetch from most reference FASTAs/BAMs).

    Returns list of (chrom, center_pos, strand), where center_pos is the midpoint
    of the site interval (0-based).
    """
    col_idx = 0 if which == "occupied" else 2
    strand_idx = col_idx + 1

    sites: List[SiteTuple] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) <= strand_idx:
                continue
            coord_str, strand = fields[col_idx], fields[strand_idx].strip()
            try:
                chrom, start, end = _parse_coord(coord_str)
            except ValueError:
                continue
            if skip_random and ("_random" in chrom or "_hap" in chrom):
                continue
            center = (start + end) // 2
            sites.append((chrom, center, strand))

    return sites

def _binned_profile_for_site(values: np.ndarray, strand: str) -> np.ndarray:
    """Reverse bin order for '-' strand sites so all profiles share orientation."""
    return values[::-1] if strand == "-" else values


def anchor_profile(get_value_fn: Callable[[str, int, int], float],
                    sites: List[SiteTuple], half_window: int = 1000,
                    bin_size: int = 25) -> np.ndarray:
    """
    Generic binned anchor-point averaging, reused for both depth and Ms profiles.

    get_value_fn(chrom, bin_start, bin_end) -> a single scalar for that bin
        (e.g. mean GC-corrected depth, or an Ms score from methylation.py).
    sites: from parse_fu2008_sites().
    half_window: bp on each side of the anchor (paper: 1 kb).
    bin_size: bp per bin (paper: 25 bp sliding windows).

    Returns: array of shape (2*half_window//bin_size,) -- mean value per bin
    across all sites, strand-oriented.
    """
    n_bins = (2 * half_window) // bin_size
    sums = np.zeros(n_bins)
    counts = np.zeros(n_bins)

    for chrom, center, strand in sites:
        bin_vals = np.full(n_bins, np.nan)
        for i in range(n_bins):
            b_start = center - half_window + i * bin_size
            b_end = b_start + bin_size
            try:
                bin_vals[i] = get_value_fn(chrom, b_start, b_end)
            except Exception:
                continue  # skip unfetchable bins (e.g. off-contig edges) rather than aborting
        bin_vals = _binned_profile_for_site(bin_vals, strand)
        valid = ~np.isnan(bin_vals)
        sums[valid] += bin_vals[valid]
        counts[valid] += 1

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(counts > 0, sums / counts, np.nan)


def ctcf_depth_profile(depth_fn: Callable[[str, int, int], np.ndarray],
                        sites: List[SiteTuple], half_window: int = 1000,
                        bin_size: int = 25) -> np.ndarray:
    """
    depth_fn(chrom, start, end) -> array of per-base GC-corrected depth for that
    range (e.g. a wrapper around gc_correction.gc_correct_depth against a
    precomputed model -- build the GC model ONCE for the whole region/chromosome
    and reuse it here; rebuilding it per-bin would be extremely slow).
    """
    def mean_depth(chrom, start, end):
        d = depth_fn(chrom, start, end)
        return float(np.nanmean(d)) if len(d) else np.nan

    return anchor_profile(mean_depth, sites, half_window, bin_size)


def ctcf_ms_profile(bam_path: str, fasta_path: str, sites: List[SiteTuple],
                     half_window: int = 1000, bin_size: int = 25,
                     min_mapq: int = 0) -> np.ndarray:
    """Ms score per 25-bp bin around each CTCF site, averaged across sites."""
    try:
        from .methylation import ms_score
    except ImportError:
        from methylation import ms_score

    def ms_for_bin(chrom, start, end):
        rate, hits, total = ms_score(bam_path, fasta_path, chrom, start, end, min_mapq=min_mapq)
        return rate if total > 0 else np.nan

    return anchor_profile(ms_for_bin, sites, half_window, bin_size)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Compute Ms profile around CTCF sites.")
    p.add_argument("sites_tsv", help="Fu et al. 2008 occupied/unoccupied site table, tab-delimited")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("--which", choices=["occupied", "unoccupied"], default="occupied")
    p.add_argument("--half-window", type=int, default=1000)
    p.add_argument("--bin-size", type=int, default=25)
    args = p.parse_args()

    sites = parse_fu2008_sites(args.sites_tsv, which=args.which)
    print(f"Loaded {len(sites)} {args.which} CTCF sites.")

    profile = ctcf_ms_profile(args.bam, args.fasta, sites, args.half_window, args.bin_size)
    print("Ms profile (per 25bp bin, -1kb to +1kb):")
    print(profile)