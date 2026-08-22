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