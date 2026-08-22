"""
Regression tests for gc_correction.py -- specifically targets the real bug
found in this project: the original rate model used observed reads as its own
denominator, so a GC stratum with zero real coverage couldn't even appear in
the model. The fix samples background positions independently, so a
zero-coverage stratum shows up as an explicit rate of 0, not an absence.
Requires pysam -- not run in the environment that wrote this file.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.gc_correction import build_gc_rate_model, gc_correct_depth


def test_rate_model_represents_all_gc_bins_not_just_covered_ones(gc_test_bam_and_fasta):
    """All 20 reads were placed only in the low-GC half. The rate model should
    still cover the FULL range of GC bins (since background sampling spans the
    whole region), with near-zero rate in the high-GC bins where no reads
    landed -- not simply missing those bins."""
    fasta_path, bam_path = gc_test_bam_and_fasta
    model = build_gc_rate_model(bam_path, fasta_path, "chrGC", 0, 1000, n_gc_bins=50)

    assert len(model) == 1, f"expected 1 read length in the model, got {len(model)}"
    rate_array = next(iter(model.values()))
    assert len(rate_array) == 50, "rate array should cover all 50 GC bins, not just covered ones"

    low_gc_bins = rate_array[:10]   # near-0% GC end
    high_gc_bins = rate_array[-10:]  # near-100% GC end
    assert low_gc_bins.mean() > high_gc_bins.mean(), (
        "low-GC bins (where all reads actually are) should show higher rate "
        "than high-GC bins (where zero reads are) -- if this fails, the model "
        "may have regressed to the original observed-reads-as-denominator bug"
    )


def test_gc_correct_depth_subtract_method_reduces_to_observed_when_no_model_entry(gc_test_bam_and_fasta):
    """If a read length isn't in the model at all, expected stays 0 for it, so
    corrected == observed exactly for that contribution (documented fallback
    behavior, not a crash)."""
    fasta_path, bam_path = gc_test_bam_and_fasta
    empty_model = {}  # simulates "no model built" / unknown read length
    corrected = gc_correct_depth(bam_path, fasta_path, "chrGC", 0, 1000, empty_model)
    from src.depth import get_region_depth
    observed = get_region_depth(bam_path, "chrGC", 0, 1000).astype(float)
    assert np.allclose(corrected, observed)


def test_divide_method_runs_without_crashing(gc_test_bam_and_fasta):
    fasta_path, bam_path = gc_test_bam_and_fasta
    model = build_gc_rate_model(bam_path, fasta_path, "chrGC", 0, 1000)
    corrected = gc_correct_depth(bam_path, fasta_path, "chrGC", 0, 1000, model, method="divide")
    assert corrected.shape == (1000,)
    assert not np.any(np.isnan(corrected))