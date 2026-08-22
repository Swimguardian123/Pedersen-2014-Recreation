"""
Regression tests for genomic_annotations.py -- specifically pins the real
cpgIslandExt column-shift bug found during this project (the file has no 'bin'
column, unlike refGene which does) so it can't silently regress. Pure file
parsing -- no pysam, no network, just synthetic gzipped files matching the
real UCSC schemas.
"""

import sys
import os
import gzip
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.genomic_annotations import (
    parse_refgene, parse_cpg_islands, get_tss_sites, get_splice_sites,
    get_gene_body_regions, periodogram_across_regions,
)
import numpy as np


def _write_gz(path, lines):
    with gzip.open(path, "wt") as f:
        for line in lines:
            f.write(line + "\n")


def test_parse_cpg_islands_uses_correct_column_offsets(tmp_path):
    """Regression test for the real bug: the actual UCSC file has NO leading
    'bin' column, so columns are chrom(0), start(1), end(2), name(3)... A
    'name'-shaped string like 'CpG: 116' appearing where a coordinate was
    expected is exactly the symptom that caught the original bug."""
    path = tmp_path / "cpgIslandExt.txt.gz"
    _write_gz(path, [
        "chr12\t1000\t1500\tCpG: 42\t500\t45\t300\t60.0\t65.0\t0.85",
        "chr12\t5000\t5300\tCpG: 20\t300\t22\t150\t55.0\t62.0\t0.80",
    ])
    islands = parse_cpg_islands(str(path))
    assert islands == [("chr12", 1000, 1500), ("chr12", 5000, 5300)]


def test_parse_refgene_tss_matches_strand_convention(tmp_path):
    """TSS = txStart for '+' genes, txEnd for '-' genes (transcription start is
    strand-relative, not always the lower genomic coordinate)."""
    path = tmp_path / "refGene.txt.gz"
    _write_gz(path, [
        # bin  name    chrom   strand  txStart txEnd  cdsStart cdsEnd exonCount exonStarts exonEnds
        "0\tNM_001\tchr12\t+\t1000\t5000\t1200\t4800\t2\t1000,3000,\t2000,5000,",
        "0\tNM_002\tchr12\t-\t2000\t6000\t2200\t5800\t2\t2000,4000,\t3000,6000,",
    ])
    genes = parse_refgene(str(path))
    tss = get_tss_sites(genes)

    assert ("chr12", 1000, "+") in tss  # + strand: TSS = txStart
    assert ("chr12", 6000, "-") in tss  # - strand: TSS = txEnd


def test_splice_sites_skip_outermost_boundaries(tmp_path):
    """Splice sites are internal exon boundaries only -- the transcript's own
    outer txStart/txEnd aren't splice sites and shouldn't appear."""
    path = tmp_path / "refGene.txt.gz"
    _write_gz(path, [
        "0\tNM_001\tchr12\t+\t1000\t5000\t1200\t4800\t2\t1000,3000,\t2000,5000,",
    ])
    genes = parse_refgene(str(path))
    splice = get_splice_sites(genes)
    positions = {s[1] for s in splice}

    assert 2000 in positions  # donor site (end of first exon) -- real splice site
    assert 3000 in positions  # acceptor site (start of second exon) -- real splice site
    assert 1000 not in positions  # txStart -- not a splice site
    assert 5000 not in positions  # txEnd -- not a splice site


def test_periodogram_across_regions_handles_a_mix_of_signal_and_empty():
    """With an injectable depth_fn (no real BAM needed), confirm regions with
    real periodic signal get a real estimate and empty regions correctly
    return NaN (exercising the min_nonzero fix from periodicity.py together
    with this module's own region-iteration logic)."""
    def fake_depth_fn(chrom, start, end):
        length = end - start
        if chrom == "chr_signal":
            x = np.arange(length)
            return 50 + 40 * np.sin(2 * np.pi * x / 193) + np.random.normal(0, 3, length)
        return np.zeros(length)  # chr_empty: no real signal

    regions = [("chr_signal", 0, 3000), ("chr_empty", 0, 3000)]
    results = periodogram_across_regions(fake_depth_fn, regions)

    assert not np.isnan(results[0]), "region with real periodic signal should get an estimate"
    assert np.isnan(results[1]), "empty region should return NaN, not a fake number"