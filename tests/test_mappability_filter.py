"""
Regression tests for mappability_filter.py.

filter_calls_by_mappability() is pure logic (no pyBigWig needed) -- tested
directly. compute_block_mappability()/get_high_mappability_regions() need a
real .bigWig file and pyBigWig -- tested with a small synthetic bigWig built
via pyBigWig's own write support, same self-contained-fixture pattern used
for synthetic BAMs elsewhere in this project.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dataclasses import dataclass
from src.mappability_filter import filter_calls_by_mappability, get_high_mappability_regions


@dataclass
class FakeCall:
    center_pos: int


def test_filter_calls_by_mappability_keeps_only_qualifying_positions():
    calls = [FakeCall(center_pos=100), FakeCall(center_pos=500), FakeCall(center_pos=25000)]
    genome_offset = 1000  # absolute positions: 1100, 1500, 26000
    regions = [(1000, 2000), (30000, 40000)]  # only the first region is relevant here

    kept = filter_calls_by_mappability(calls, genome_offset, regions)
    assert len(kept) == 2
    assert {c.center_pos for c in kept} == {100, 500}


def test_get_high_mappability_regions_against_synthetic_bigwig(tmp_path):
    import pyBigWig

    bw_path = str(tmp_path / "test_mappability.bigWig")
    bw = pyBigWig.open(bw_path, "w")
    bw.addHeader([("chrTest", 100000)])
    # first 20kb block: mappability 1.0 (should pass threshold 0.9)
    # second 20kb block: mappability 0.5 (should NOT pass)
    bw.addEntries(["chrTest", "chrTest"], [0, 20000], ends=[20000, 40000], values=[1.0, 0.5])
    bw.close()

    regions = get_high_mappability_regions(bw_path, "chrTest", 0, 40000,
                                            block_size=20000, threshold=0.9)
    assert regions == [(0, 20000)], f"expected only the first block to pass, got {regions}"


if __name__ == "__main__":
    test_filter_calls_by_mappability_keeps_only_qualifying_positions()
    print("(bigWig test needs pytest fixtures + pyBigWig -- run via `pytest tests/`)")