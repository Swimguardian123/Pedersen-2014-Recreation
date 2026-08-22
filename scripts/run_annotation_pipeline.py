"""
Driver script: smoke-test genomic_annotations.py's TSS/splice-site/CGI/gene-body
data feeding into the existing periodicity + methylation machinery.

TSS and splice sites are point anchors -> reuse ctcf_analysis.py's generic
anchor_profile() directly (it was built generic for exactly this reuse, not
CTCF-specific despite living in that file). Gene bodies and CpG islands are
region spans -> genomic_annotations.periodogram_across_regions().
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pysam

from src.genomic_annotations import (
    download_ucsc_table, parse_refgene, parse_cpg_islands,
    get_tss_sites, get_splice_sites, get_gene_body_regions,
    periodogram_across_regions,
)
from src.ctcf_analysis import anchor_profile
from src.methylation import ms_score
from src.gc_correction import build_gc_rate_model, gc_correct_depth


def run(bam_path: str, fasta_path: str, chrom: str, annotations_dir: str = "./data/annotations",
        half_window: int = 1000, bin_size: int = 25, max_tss: int = None,
        max_splice: int = None, max_gene_bodies: int = None, max_cgis: int = None) -> None:
    print(f"=== Loading annotations for {chrom} ===")
    refgene_path = download_ucsc_table("refGene", annotations_dir)
    genes = [g for g in parse_refgene(refgene_path) if g.chrom == chrom]
    cgi_path = download_ucsc_table("cpgIslandExt", annotations_dir)
    cgis = [c for c in parse_cpg_islands(cgi_path) if c[0] == chrom]
    print(f"  {len(genes)} genes, {len(cgis)} CpG islands on {chrom}\n")

    tss_sites = get_tss_sites(genes)
    splice_sites = get_splice_sites(genes)
    gene_bodies = get_gene_body_regions(genes)

    # ms_score() opens a fresh BAM/FASTA handle per call -- with 40+ bins per
    # anchor and thousands of anchors, that's millions of file opens (this is
    # what made splice sites hang: 78,068 sites x 40 bins = ~3.1M opens). For a
    # smoke test (confirming mechanics, not exhaustive computation), subsampling
    # is the practical fix rather than a bigger refactor of methylation.py.
    import random
    random.seed(0)
    if max_tss and len(tss_sites) > max_tss:
        tss_sites = random.sample(tss_sites, max_tss)
    if max_splice and len(splice_sites) > max_splice:
        splice_sites = random.sample(splice_sites, max_splice)
    if max_gene_bodies and len(gene_bodies) > max_gene_bodies:
        gene_bodies = random.sample(gene_bodies, max_gene_bodies)
    if max_cgis and len(cgis) > max_cgis:
        cgis = random.sample(cgis, max_cgis)
    print(f"  Using {len(tss_sites)} TSS, {len(splice_sites)} splice sites, "
          f"{len(gene_bodies)} gene bodies, {len(cgis)} CGIs after any subsampling\n")

    print("[1/3] Ms profile anchored at TSS (methylation.py + anchor_profile)")
    def ms_for_bin(c, s, e):
        rate, hits, total = ms_score(bam_path, fasta_path, c, s, e)
        return rate if total > 0 else np.nan
    tss_ms_profile = anchor_profile(ms_for_bin, tss_sites, half_window, bin_size)
    print(f"  profile shape={tss_ms_profile.shape}, "
          f"valid bins={(~np.isnan(tss_ms_profile)).sum()}\n")

    print("[2/3] Ms profile anchored at splice sites")
    splice_ms_profile = anchor_profile(ms_for_bin, splice_sites, half_window=200, bin_size=10)
    print(f"  profile shape={splice_ms_profile.shape}, "
          f"valid bins={(~np.isnan(splice_ms_profile)).sum()}\n")

    print("[3/3] Periodicity across gene bodies and CpG islands (non-anchored)")
    bam = pysam.AlignmentFile(bam_path, "rb")
    chrom_len = bam.get_reference_length(chrom)
    bam.close()
    gc_model = build_gc_rate_model(bam_path, fasta_path, chrom, 0, chrom_len)

    def depth_fn(c, s, e):
        return gc_correct_depth(bam_path, fasta_path, c, max(s, 0), min(e, chrom_len), gc_model)

    gene_body_periods = periodogram_across_regions(depth_fn, gene_bodies)
    valid_gb = sum(1 for p in gene_body_periods if not np.isnan(p))
    print(f"  gene bodies: {len(gene_body_periods)} regions, {valid_gb} with a usable estimate")

    cgi_regions = [(c[0], c[1], c[2]) for c in cgis]
    cgi_periods = periodogram_across_regions(depth_fn, cgi_regions)
    valid_cgi = sum(1 for p in cgi_periods if not np.isnan(p))
    print(f"  CpG islands: {len(cgi_periods)} regions, {valid_cgi} with a usable estimate")

    print(f"\n(At this coverage, expect most anchors/regions to show no usable signal --")
    print(f" this confirms the annotation-loading + generic anchor/periodogram machinery")
    print(f" runs correctly against real data, same caveat as every other smoke test today.)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Smoke-test TSS/splice/CGI/gene-body analyses against real data.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("--annotations-dir", default="./data/annotations")
    p.add_argument("--half-window", type=int, default=1000)
    p.add_argument("--bin-size", type=int, default=25)
    p.add_argument("--max-tss", type=int, default=300)
    p.add_argument("--max-splice", type=int, default=300)
    p.add_argument("--max-gene-bodies", type=int, default=500)
    p.add_argument("--max-cgis", type=int, default=500)
    args = p.parse_args()

    run(args.bam, args.fasta, args.chrom, args.annotations_dir, args.half_window, args.bin_size,
        args.max_tss, args.max_splice, args.max_gene_bodies, args.max_cgis)