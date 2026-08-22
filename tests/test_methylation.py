"""
Regression tests for methylation.py -- specifically for the reverse-strand
C->T/G->A logic that was hand-verified by careful reasoning earlier in this
project (see PROJECT_SUMMARY.md), not previously covered by an automated test.
Requires the synthetic_ref_and_bam fixture (see conftest.py) -- needs pysam,
not run in the environment that wrote this file.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.methylation import ms_score, find_cpg_positions


def test_find_cpg_positions_locates_both_planted_cpgs(synthetic_ref_and_bam):
    fasta_path, bam_path = synthetic_ref_and_bam
    cpgs = find_cpg_positions(fasta_path, "chrTest", 0, 300)
    assert 50 in cpgs
    assert 100 in cpgs


def test_forward_strand_deamination_is_detected(synthetic_ref_and_bam):
    """The forward-strand read has a deliberate C->T mismatch at the CpG it
    starts on -- Ms score for a region containing just that read+CpG should
    show exactly 1 hit out of 1 total."""
    fasta_path, bam_path = synthetic_ref_and_bam
    rate, hits, total = ms_score(fasta_path=fasta_path, bam_path=bam_path,
                                  chrom="chrTest", start=50, end=51)
    assert total == 1, f"expected exactly 1 read-start checked, got {total}"
    assert hits == 1, f"expected the deamination to be detected, got {hits} hits"
    assert rate == 1.0


def test_reverse_strand_deamination_is_detected(synthetic_ref_and_bam):
    """The reverse-strand read has a deliberate G->A mismatch at the CpG its
    alignment ends on -- this is the trickier of the two cases (BAM-forward-
    orientation storage means the deamination shows as G->A, not C->T, at the
    RIGHTMOST position of a reverse-strand read)."""
    fasta_path, bam_path = synthetic_ref_and_bam
    rate, hits, total = ms_score(fasta_path=fasta_path, bam_path=bam_path,
                                  chrom="chrTest", start=100, end=101)
    assert total == 1, f"expected exactly 1 read-start checked, got {total}"
    assert hits == 1, f"expected the reverse-strand deamination to be detected, got {hits} hits"
    assert rate == 1.0


def test_clean_read_region_shows_no_false_hits(synthetic_ref_and_bam):
    """The clean read at position 200 doesn't start on a CpG at all -- should
    show 0 total (nothing to check), not a false hit."""
    fasta_path, bam_path = synthetic_ref_and_bam
    rate, hits, total = ms_score(fasta_path=fasta_path, bam_path=bam_path,
                                  chrom="chrTest", start=200, end=201)
    assert total == 0
    assert hits == 0


def test_excluded_sites_logic():
    """Pure-logic test (no BAM needed) for the site-exclusion refinement
    (Hanghoj et al. 2016): >50% flip rate AND >=5x coverage -> excluded."""
    from src.methylation import excluded_sites
    stats = {
        100: [6, 10],   # 60% flip, cov 10 -> excluded
        200: [2, 10],   # 20% flip, cov 10 -> kept
        300: [3, 4],    # 75% flip, cov 4 (< 5) -> kept (insufficient coverage)
    }
    excluded = excluded_sites(stats, min_coverage=5, max_flip_fraction=0.5)
    assert excluded == {100}


def test_ms_score_filtered_excludes_polymorphism_like_site(polymorphism_site_bam):
    """End-to-end: a site with 6/10 reads showing the flip (60%, real
    polymorphism-like pattern) should be excluded by ms_score_filtered, leaving
    total=0 for the region -- versus plain ms_score(), which pools it in
    regardless and reports total=10."""
    from src.methylation import ms_score, ms_score_filtered

    fasta_path, bam_path = polymorphism_site_bam

    plain_rate, plain_hits, plain_total = ms_score(bam_path, fasta_path, "chrPoly", 100, 101)
    assert plain_total == 10  # unfiltered: pools the polymorphism-like site in

    filt_rate, filt_hits, filt_total, excluded = ms_score_filtered(
        bam_path, fasta_path, "chrPoly", 100, 101)
    assert 100 in excluded
    assert filt_total == 0  # the only site in this tiny region was excluded
    assert np.isnan(filt_rate)