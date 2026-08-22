"""
Gene/CpG-island annotations (hg18, UCSC), for the periodicity and methylation
anchor analyses the paper runs beyond CTCF sites.

Paper reference (Methods, "Read-depth periodicity"): "Spectral density plots
(periodograms)... across CpG islands, TSSs +/- 1000 bp, gene bodies, and CTCF
sites +/- 1000 bp." CTCF is covered by ctcf_analysis.py; this module supplies the
other three site/region sources, all public UCSC downloads for hg18 -- no
paywalled or hard-to-find data involved, unlike the CTCF table.

TSS and splice sites are point anchors (used with ctcf_analysis.anchor_profile,
which is already fully generic -- built for exactly this reuse). CpG islands and
gene bodies are variable-length REGIONS, not point anchors: the paper runs the
periodogram directly across each region's own span, not a fixed window around a
point -- a different, simpler analysis pattern implemented here separately.

Source: UCSC refGene table (TSS/exons/gene bodies) and cpgIslandExt table (CGIs),
both at hgdownload.cse.ucsc.edu/goldenPath/hg18/database/ -- public, no accession
hunting needed, unlike the ancient-genome data.
"""

import gzip
from dataclasses import dataclass
from typing import List, Tuple
from urllib.request import urlretrieve

try:
    from .ctcf_analysis import SiteTuple
except ImportError:
    from ctcf_analysis import SiteTuple

UCSC_HG18_DB = "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database"


@dataclass
class GeneRecord:
    name: str
    chrom: str
    strand: str
    tx_start: int
    tx_end: int
    exon_starts: List[int]
    exon_ends: List[int]


def download_ucsc_table(table_name: str, destdir: str) -> str:
    """Download a UCSC hg18 database table (e.g. 'refGene', 'cpgIslandExt')."""
    import os
    os.makedirs(destdir, exist_ok=True)
    dest_path = os.path.join(destdir, f"{table_name}.txt.gz")
    if os.path.exists(dest_path):
        print(f"{dest_path} already exists, skipping download.")
        return dest_path
    url = f"{UCSC_HG18_DB}/{table_name}.txt.gz"
    print(f"Downloading {url} ...")
    urlretrieve(url, dest_path)
    return dest_path


def parse_refgene(path: str) -> List[GeneRecord]:
    """
    Parse a refGene.txt.gz table. Columns (UCSC refGene schema): bin, name,
    chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts,
    exonEnds, ... (exonStarts/Ends are comma-terminated lists).
    """
    records = []
    with gzip.open(path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            name, chrom, strand = fields[1], fields[2], fields[3]
            tx_start, tx_end = int(fields[4]), int(fields[5])
            exon_starts = [int(x) for x in fields[9].rstrip(",").split(",") if x]
            exon_ends = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            records.append(GeneRecord(name, chrom, strand, tx_start, tx_end,
                                       exon_starts, exon_ends))
    return records


def get_tss_sites(genes: List[GeneRecord]) -> List[SiteTuple]:
    """TSS = txStart for '+' strand genes, txEnd for '-' strand (transcription
    start relative to gene orientation, not genome coordinate order)."""
    sites = []
    for g in genes:
        tss = g.tx_start if g.strand == "+" else g.tx_end
        sites.append((g.chrom, tss, g.strand))
    return sites


def get_splice_sites(genes: List[GeneRecord]) -> List[SiteTuple]:
    """
    Every exon boundary (5' and 3' splice site) as its own anchor point, strand-
    tagged. Internal exons only (skips the transcript's outermost tx_start/
    tx_end, which aren't splice sites).
    """
    sites = []
    for g in genes:
        n_exons = len(g.exon_starts)
        for i in range(n_exons):
            if i > 0:  # 3' end of intron / acceptor site before this exon
                sites.append((g.chrom, g.exon_starts[i], g.strand))
            if i < n_exons - 1:  # 5' end of intron / donor site after this exon
                sites.append((g.chrom, g.exon_ends[i], g.strand))
    return sites


def get_gene_body_regions(genes: List[GeneRecord]) -> List[Tuple[str, int, int, str]]:
    """(chrom, start, end, strand) for each gene's full transcript span."""
    return [(g.chrom, g.tx_start, g.tx_end, g.strand) for g in genes]


def parse_cpg_islands(path: str) -> List[Tuple[str, int, int]]:
    """
    Parse cpgIslandExt.txt.gz.

    Schema note (confirmed against the real hg18 file, not assumed): this table's
    flat-file dump does NOT include the leading 'bin' column that its SQL schema
    lists -- columns start directly at chrom, chromStart, chromEnd, name, length,
    cpgNum, gcNum, perCpg, perGc, obsExp. An earlier version of this parser
    assumed a bin column (matching refGene, which DOES have one) and was off by
    one field as a result -- caught by an actual run against the real file.
    """
    islands = []
    with gzip.open(path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            islands.append((chrom, start, end))
    return islands


def periodogram_across_regions(depth_fn, regions: List[Tuple[str, int, int, str]]
                                ) -> List[float]:
    """
    Run periodicity.dominant_period() directly across each region's own span
    (not an anchored window) -- the pattern the paper uses for CpG islands and
    gene bodies, as opposed to the point-anchored TSS/splice-site/CTCF pattern.

    depth_fn(chrom, start, end) -> depth array for that region (e.g. a
    GC-corrected depth wrapper, same pattern as ctcf_analysis.ctcf_depth_profile).
    Returns one dominant-period estimate per region (NaN where a region has no
    usable signal).
    """
    try:
        from .periodicity import dominant_period
    except ImportError:
        from periodicity import dominant_period

    import numpy as np
    results = []
    for region in regions:
        chrom, start, end = region[0], region[1], region[2]
        try:
            depth = depth_fn(chrom, start, end)
            if len(depth) < 20:
                results.append(float("nan"))
                continue
            results.append(dominant_period(depth.astype(np.float64)))
        except Exception:
            results.append(float("nan"))
    return results


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Download and parse UCSC hg18 gene/CGI annotations.")
    p.add_argument("--destdir", default="./data/annotations")
    p.add_argument("--chrom-filter", default=None, help="Only keep records on this chromosome")
    args = p.parse_args()

    refgene_path = download_ucsc_table("refGene", args.destdir)
    genes = parse_refgene(refgene_path)
    if args.chrom_filter:
        genes = [g for g in genes if g.chrom == args.chrom_filter]
    print(f"Parsed {len(genes)} refGene record(s)"
          + (f" on {args.chrom_filter}" if args.chrom_filter else ""))

    tss = get_tss_sites(genes)
    splice = get_splice_sites(genes)
    bodies = get_gene_body_regions(genes)
    print(f"  TSS sites: {len(tss)}, splice sites: {len(splice)}, gene bodies: {len(bodies)}")

    cgi_path = download_ucsc_table("cpgIslandExt", args.destdir)
    cgis = parse_cpg_islands(cgi_path)
    if args.chrom_filter:
        cgis = [c for c in cgis if c[0] == args.chrom_filter]
    print(f"  CpG islands: {len(cgis)}")