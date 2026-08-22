"""
Regression tests for control_set.py -- the one core module that had never
been tested via pytest before, only through one live end-to-end run.
Requires pysam -- not run in the environment that wrote this file.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from src.control_set import get_read_length_distribution, build_control_bam
import pysam


def test_read_length_distribution_matches_target_bam(control_set_bams):
    source_path, target_path = control_set_bams
    dist = get_read_length_distribution(target_path)
    assert len(dist) == 5
    assert all(length == 15 for length in dist), f"expected all target reads to be 15bp, got {dist}"


def test_control_bam_matches_target_read_count_and_length(control_set_bams, tmp_path):
    """Per the paper's stated method: same read COUNT as target, each read
    truncated toward the target's length distribution."""
    source_path, target_path = control_set_bams
    target_dist = get_read_length_distribution(target_path)  # all 15bp

    out_path = tmp_path / "control_out.bam"
    build_control_bam(source_path, target_dist, str(out_path), n_reads=5, seed=0)

    result_bam = str(out_path) + ".sorted.bam"
    bam = pysam.AlignmentFile(result_bam, "rb")
    reads = list(bam.fetch())
    bam.close()

    assert len(reads) == 5, f"expected 5 reads (matching target count), got {len(reads)}"
    for r in reads:
        # source reads are 40bp, target length is 15bp -> truncated length should be
        # min(target_len, orig_len) = 15 for every read
        assert r.query_alignment_length == 15, (
            f"expected truncation to 15bp (target length), got {r.query_alignment_length}"
        )


def test_requesting_more_reads_than_source_has_raises_clear_error(control_set_bams, tmp_path):
    source_path, target_path = control_set_bams  # source only has 20 reads
    out_path = tmp_path / "control_out2.bam"
    with pytest.raises(ValueError):
        build_control_bam(source_path, [15] * 100, str(out_path), n_reads=100, seed=0)