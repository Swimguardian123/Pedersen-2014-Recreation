"""
Mono- and dinucleotide positional patterns across nucleosome calls (Fig. 3C).

Paper reference (Methods, "Nucleotide patterns across nucleosomes"):
"Mono- and dinucleotide distributions across nucleosomes were produced by ranking
and stratifying nucleosome calls by score. All nucleosomes within each group were
aligned at the center (dyad) position, and the average usage across each position
was calculated."

Results text gives concrete numbers to validate against once you have real data:
- 10.2-bp periodicity in single-nucleotide frequencies (Brogaard-style helical turn).
- G/C frequencies vary sharply near the dyad: ~35% G at position -2 (C at +2),
  dropping to ~14% G at +1 (C at -1).
- Pattern is reverse-complement symmetric across the dyad (C vs A, G vs T mirrored).
- Dinucleotide analysis: no yeast-style strong/weak (W/S) dinucleotide pattern;
  instead a 10-bp periodic purine/pyrimidine (R/Y) dinucleotide signal, amplitude
  increasing toward the center, reverse-complemented there.

No stratification threshold (e.g. "top 25%") is a fixed universal constant --
Pedersen used it for one specific comparison ("top-scoring nucleosome calls");
`top_fraction` here defaults to None (use all calls) so you choose per-analysis.
"""

from typing import List, Optional

import numpy as np
import pysam

try:
    from .nucleosome_calling import NucleosomeCall
except ImportError:
    from nucleosome_calling import NucleosomeCall

BASES = ["A", "C", "G", "T"]
PURINES = set("AG")


def _is_purine(base: str) -> bool:
    return base.upper() in PURINES


def select_top_calls(calls: List[NucleosomeCall], top_fraction: Optional[float] = None
                      ) -> List[NucleosomeCall]:
    """Stratify by score; top_fraction=0.25 reproduces the paper's 'top-scoring
    nucleosome calls' quartile comparisons. None = use all calls unranked."""
    if top_fraction is None:
        return calls
    ranked = sorted(calls, key=lambda c: c.score, reverse=True)
    n_keep = max(1, int(len(ranked) * top_fraction))
    return ranked[:n_keep]


def mononucleotide_matrix(fasta_path: str, chrom: str, calls: List[NucleosomeCall],
                           genome_offset: int = 0, halfwidth: int = 100) -> np.ndarray:
    """
    Align every call at its dyad (center_pos) and tally per-position base
    frequencies across the aligned stack. `genome_offset` converts a call's
    center_pos (relative to whatever depth array it came from) back to absolute
    genome coordinates -- pass the same `start` you gave call_nucleosomes' source
    region.

    Returns: array of shape (2*halfwidth+1, 4), columns ordered as BASES, giving
    per-position base frequency (rows sum to 1, ignoring any N/ambiguous bases).
    """
    fasta = pysam.FastaFile(fasta_path)
    counts = np.zeros((2 * halfwidth + 1, 4), dtype=np.float64)
    n_used = np.zeros(2 * halfwidth + 1, dtype=np.int64)

    for call in calls:
        center = call.center_pos + genome_offset
        try:
            seq = fasta.fetch(chrom, center - halfwidth, center + halfwidth + 1).upper()
        except (ValueError, IndexError):
            continue
        if len(seq) != 2 * halfwidth + 1:
            continue
        for i, base in enumerate(seq):
            if base in BASES:
                counts[i, BASES.index(base)] += 1
                n_used[i] += 1

    fasta.close()
    with np.errstate(divide="ignore", invalid="ignore"):
        freqs = np.where(n_used[:, None] > 0, counts / n_used[:, None], np.nan)
    return freqs


def purine_pyrimidine_dinucleotide_signal(fasta_path: str, chrom: str,
                                           calls: List[NucleosomeCall],
                                           genome_offset: int = 0,
                                           halfwidth: int = 100) -> np.ndarray:
    """
    Per-position frequency of "purine-purine" dinucleotide steps (RR), the
    strand-specific signal the paper reports as showing 10-bp periodicity with
    amplitude increasing toward the dyad and reverse-complementing there (Fig. 3C).
    A parallel YY (pyrimidine-pyrimidine) track is the natural reverse-complement
    counterpart -- compute both and compare/subtract as needed downstream.

    Returns: array of shape (2*halfwidth,) -- one fewer than mononucleotide since
    dinucleotides span position i, i+1. Value = P(purine at i AND purine at i+1).
    """
    fasta = pysam.FastaFile(fasta_path)
    rr_counts = np.zeros(2 * halfwidth, dtype=np.float64)
    n_used = np.zeros(2 * halfwidth, dtype=np.int64)

    for call in calls:
        center = call.center_pos + genome_offset
        try:
            seq = fasta.fetch(chrom, center - halfwidth, center + halfwidth + 1).upper()
        except (ValueError, IndexError):
            continue
        if len(seq) != 2 * halfwidth + 1:
            continue
        for i in range(len(seq) - 1):
            b1, b2 = seq[i], seq[i + 1]
            if b1 in BASES and b2 in BASES:
                if _is_purine(b1) and _is_purine(b2):
                    rr_counts[i] += 1
                n_used[i] += 1

    fasta.close()
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(n_used > 0, rr_counts / n_used, np.nan)


if __name__ == "__main__":
    import argparse
    import csv

    p = argparse.ArgumentParser(description="Compute nucleotide patterns across nucleosome calls.")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("calls_csv", help="CSV from nucleosome_calling.py (center_pos,peak_depth,score)")
    p.add_argument("--offset", type=int, default=0, help="genome_offset to add to center_pos")
    p.add_argument("--halfwidth", type=int, default=100)
    p.add_argument("--top-fraction", type=float, default=None)
    p.add_argument("--out", default=None)
    args = p.parse_args()

    calls = []
    with open(args.calls_csv) as f:
        for row in csv.DictReader(f):
            calls.append(NucleosomeCall(int(row["center_pos"]), float(row["peak_depth"]), float(row["score"])))

    selected = select_top_calls(calls, args.top_fraction)
    print(f"Using {len(selected)}/{len(calls)} calls.")

    mono = mononucleotide_matrix(args.fasta, args.chrom, selected, args.offset, args.halfwidth)
    rr = purine_pyrimidine_dinucleotide_signal(args.fasta, args.chrom, selected, args.offset, args.halfwidth)

    if args.out:
        np.savez(args.out, mononucleotide=mono, purine_purine_dinuc=rr, bases=BASES)
        print(f"Saved to {args.out}.npz")
    else:
        print("Mononucleotide matrix shape:", mono.shape)
        print("Purine-purine dinucleotide signal shape:", rr.shape)