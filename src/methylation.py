"""
Cytosine methylation proxy (Ms score) from C->T misincorporation at read starts.

Paper reference (Methods, "Methylation signal"):
"We defined a proxy for regional methylation levels, Ms as the fraction of CpG
dinucleotides giving rise to TpG misincorporations at read starts."

Results section, mechanism and validation targets:
- "Focusing on the first position where deamination rates are maximal, we observed
  a 5.04-fold increase in C->T errors for Phusion reads starting at CpG" (vs. non-CpG).
- "This pattern was absent for Phusion reads starting with CpA, CpT, and CpC" --
  i.e. Ms should be computed the same way for those contexts as a negative control.
- "CpG->TpG conversions at read ends and other positions within reads were
  disregarded" -- ONLY the very first sequenced base (biological 5' end) counts,
  not any other read position.
- Run this only on Phusion-amplified-library reads; HiFi reads don't carry the
  methylation signal the same way (see paper's polymerase-bypass logic) -- that's
  a BAM-selection step upstream of this module, not something this code enforces.

Strand handling (this is the part worth getting right, not stated explicitly in the
paper but required for correctness): BAM always stores query_sequence in reference-
forward orientation, regardless of which strand a read came from.
  - Forward-strand read: biological 5' end = leftmost aligned base = reference_start.
    Damage signature there is the classic C->T (ref C, read T).
  - Reverse-strand read: biological 5' end = rightmost aligned base = reference_end-1,
    but since the stored sequence is the reverse-complement of what was actually
    sequenced, the damage signature at that position shows up as G->A (ref G, read A)
    when read off the reference-forward-oriented BAM sequence.
Because CpG is a palindrome (revcomp("CG") == "CG"), the reference-level test for
"is this a CpG site" is identical in both cases -- ref[p]=='C' and ref[p+1]=='G' --
only the strand and which base of the pair to check for a mismatch differs:
  - forward read starting at p       -> check position p for C->T
  - reverse read ending at p+1        -> check position p+1 for G->A
"""

from typing import Optional, Tuple

import numpy as np
import pysam


def find_cpg_positions(fasta_path: str, chrom: str, start: int, end: int) -> np.ndarray:
    """Reference positions p (0-based) where ref[p:p+2] == 'CG', within [start, end)."""
    fasta = pysam.FastaFile(fasta_path)
    # pad by 1 so a CpG straddling the window edge is still found
    seq = fasta.fetch(chrom, start, end + 1).upper()
    fasta.close()
    positions = [start + i for i in range(len(seq) - 1) if seq[i:i + 2] == "CG"]
    return np.array(positions, dtype=np.int64)


def mismatch_rate_at_read_start(bam_path: str, fasta_path: str, chrom: str,
                                 start: int, end: int, ref_dinuc: str = "CG",
                                 min_mapq: int = 0) -> Tuple[float, int, int]:
    """
    General version of the Ms calculation, reproducing Fig. 4B's comparison across
    CpG / CpA / CpT / CpC contexts (ref_dinuc = "CG", "CA", "CT", "CC").

    Only the read's biological first sequenced base is checked (matches the paper's
    explicit exclusion of read-end and internal positions). Only counts read-starts
    whose reference context is ref_dinuc; a C->T (or complement G->A, reverse-strand)
    mismatch there is the "hit".

    Returns: (rate, n_hits, n_total) where rate = n_hits / n_total (nan if n_total==0).
    """
    if ref_dinuc[0] != "C":
        raise ValueError("ref_dinuc must start with 'C' (this measures deamination at C)")
    second_base = ref_dinuc[1]

    fasta = pysam.FastaFile(fasta_path)
    ref_seq = fasta.fetch(chrom, start, end + 1).upper()  # +1 pad for dinucleotide lookups at edge
    fasta.close()

    def ref_base(pos: int) -> str:
        idx = pos - start
        if 0 <= idx < len(ref_seq):
            return ref_seq[idx]
        return "N"

    bam = pysam.AlignmentFile(bam_path, "rb")
    n_hits = 0
    n_total = 0

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue
        seq = read.query_sequence
        if not seq:
            continue

        if not read.is_reverse:
            p = read.reference_start
            if ref_base(p) != "C" or ref_base(p + 1) != second_base:
                continue
            n_total += 1
            if seq[0].upper() == "T":
                n_hits += 1
        else:
            p = read.reference_end - 1  # last aligned reference base = biological 5' end
            if ref_base(p) != second_base or ref_base(p - 1) != "C":
                continue
            n_total += 1
            if seq[-1].upper() == "A":
                n_hits += 1

    bam.close()
    rate = n_hits / n_total if n_total > 0 else float("nan")
    return rate, n_hits, n_total


def ms_score(bam_path: str, fasta_path: str, chrom: str, start: int, end: int,
             min_mapq: int = 0) -> Tuple[float, int, int]:
    """Ms = fraction of CpG-context read-starts showing the deamination mismatch."""
    return mismatch_rate_at_read_start(bam_path, fasta_path, chrom, start, end,
                                        ref_dinuc="CG", min_mapq=min_mapq)


def ms_negative_controls(bam_path: str, fasta_path: str, chrom: str, start: int, end: int,
                          min_mapq: int = 0) -> dict:
    """
    Reproduces the paper's validation check (Fig. 4B): Ms-style rate computed at
    CpA/CpT/CpC contexts, which should show ~no elevated signal versus CpG.
    """
    return {
        ctx: mismatch_rate_at_read_start(bam_path, fasta_path, chrom, start, end,
                                          ref_dinuc="C" + ctx, min_mapq=min_mapq)
        for ctx in ["A", "T", "C"]
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Compute Ms methylation proxy for a region.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--controls", action="store_true", help="Also compute CpA/CpT/CpC negative controls")
    args = p.parse_args()

    rate, hits, total = ms_score(args.bam, args.fasta, args.chrom, args.start, args.end)
    print(f"Ms = {rate:.4f} ({hits}/{total} CpG read-starts showing C->T/G->A)")

    if args.controls:
        ctrl = ms_negative_controls(args.bam, args.fasta, args.chrom, args.start, args.end)
        for ctx, (r, h, t) in ctrl.items():
            print(f"Cp{ctx}: rate = {r:.4f} ({h}/{t})")