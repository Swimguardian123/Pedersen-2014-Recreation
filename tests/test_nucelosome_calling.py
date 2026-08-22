"""
Regression tests for nucleosome_calling.py. Pure numpy -- no BAM/FASTA needed.

These specifically encode behaviors validated by hand earlier in this project
(see PROJECT_SUMMARY.md), so a future edit that breaks one of them will be
caught automatically instead of requiring another manual round of "run it and
eyeball the output."
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.nucleosome_calling import call_nucleosomes, score_threshold_for_fdr


def test_periodic_signal_produces_roughly_expected_call_count():
    """A signal built with exact 193bp periodicity over 5000bp should produce
    roughly one call per period (~26), not wildly more or fewer."""
    np.random.seed(0)
    n = 5000
    x = np.arange(n)
    signal = 50 + 40 * np.sin(2 * np.pi * x / 193) + np.random.normal(0, 5, n)
    signal = np.clip(signal, 0, None)

    calls = call_nucleosomes(signal, flank=25, offset=12)
    assert 15 <= len(calls) <= 35, f"expected ~26 calls, got {len(calls)}"


def test_empty_signal_produces_no_calls_and_does_not_crash():
    calls = call_nucleosomes(np.zeros(300))
    assert calls == []


def test_score_equals_peak_minus_flank_average():
    """Directly checks the paper's stated formula: score = peak - avg(flanks)."""
    depth = np.zeros(300)
    depth[150] = 10.0  # isolated peak, zero flanks
    calls = call_nucleosomes(depth, flank=25, offset=12)
    assert len(calls) == 1
    assert calls[0].peak_depth == 10.0
    assert calls[0].score == 10.0  # flanks are 0, so score == peak exactly


def test_flank_geometry_matches_hanghoj_2016_figure_1b():
    """Regression test for the flank-width/offset correction: flanks should be
    25bp wide, starting 12bp beyond the 147bp window's edge -- NOT 50bp and
    adjacent, which was the earlier (unconfirmed) guess. Places a nonzero value
    exactly in the gap (should be excluded from the flank) and exactly in the
    real flank region (should be included), to directly test the geometry.

    Note: this checks the specific call at center_pos=200, rather than
    asserting the total call count -- an isolated nonzero value placed anywhere
    in the array (like the gap-region probe value below) is trivially its own
    local max and legitimately becomes its own separate call. That's correct
    behavior, not a bug; asserting len(calls)==1 was an overly strict
    assumption in an earlier version of this test, not a real requirement."""
    depth = np.zeros(400)
    center = 200
    half = 73  # 147 // 2

    # a value placed in the 12bp GAP (just outside the window, before the flank
    # starts) should NOT be included in the flank average for the center=200 call
    gap_position = center - half - 5  # 5bp into the 12bp gap
    depth[gap_position] = 100.0

    # a value placed inside the real 25bp flank region SHOULD be included
    flank_position = center - half - 12 - 10  # 10bp into the 25bp flank
    depth[flank_position] = 8.0

    depth[center] = 20.0  # the peak itself

    calls = call_nucleosomes(depth, flank=25, offset=12)
    center_calls = [c for c in calls if c.center_pos == center]
    assert len(center_calls) == 1, f"expected exactly one call at center_pos={center}"
    call = center_calls[0]

    # left flank mean should reflect ONLY the 8.0 value (25bp wide, one nonzero
    # entry) -- if the gap-region 100.0 leaked in, this would be wildly higher
    expected_left_flank_mean = 8.0 / 25
    expected_score = 20.0 - (expected_left_flank_mean + 0.0) / 2
    assert abs(call.score - expected_score) < 1e-6, (
        f"expected score ~{expected_score:.4f} (gap value correctly excluded), "
        f"got {call.score:.4f} -- flank geometry may have regressed"
    )


def test_fdr_threshold_is_monotonic_with_target():
    """A stricter (lower) target FDR should never produce a lower score threshold
    than a looser one, for the same score distributions."""
    np.random.seed(1)
    saqqaq_scores = np.random.exponential(5, 200) + 10  # real calls, higher scores
    control_scores = np.random.exponential(5, 200)      # null calls, lower scores

    strict = score_threshold_for_fdr(saqqaq_scores, control_scores, target_fdr=0.01)
    loose = score_threshold_for_fdr(saqqaq_scores, control_scores, target_fdr=0.20)

    # explicit, not conditional: if either came back None, that's a real failure,
    # not something to silently skip past (a previous version of this test
    # only checked monotonicity INSIDE an `if both not None` guard, which meant
    # a completely broken implementation returning None for everything would
    # have reported PASS without checking anything -- fixed here)
    assert strict is not None, "expected a real threshold for target_fdr=0.01, got None"
    assert loose is not None, "expected a real threshold for target_fdr=0.20, got None"
    assert strict >= loose


if __name__ == "__main__":
    test_periodic_signal_produces_roughly_expected_call_count()
    test_empty_signal_produces_no_calls_and_does_not_crash()
    test_score_equals_peak_minus_flank_average()
    test_fdr_threshold_is_monotonic_with_target()
    print("All nucleosome_calling tests passed.")