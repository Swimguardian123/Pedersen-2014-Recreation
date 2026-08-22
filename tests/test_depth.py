"""
Regression tests for depth.py against the synthetic fixture. Needs pysam --
not run in the environment that wrote this file.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.depth import get_region_depth


def test_depth_matches_known_read_placement(synthetic_ref_and_bam):
    """Three reads at known, non-overlapping positions (50, 72, 200) should
    each contribute exactly 1x depth over their own 30bp span, and 0 elsewhere."""
    fasta_path, bam_path = synthetic_ref_and_bam
    depth = get_region_depth(bam_path, "chrTest", 0, 300)

    assert depth[50] == 1  # forward read starts here
    assert depth[200] == 1  # clean read starts here
    assert depth[0] == 0    # nothing placed at the very start
    assert depth[290] == 0  # nothing placed at the very end
    assert depth.sum() == 90  # 3 reads x 30bp each, no overlaps