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

    calls = call_nucleosomes(signal, flank=50)
    assert 15 <= len(calls) <= 35, f"expected ~26 calls, got {len(calls)}"


def test_empty_signal_produces_no_calls_and_does_not_crash():
    calls = call_nucleosomes(np.zeros(300))
    assert calls == []


def test_score_equals_peak_minus_flank_average():
    """Directly checks the paper's stated formula: score = peak - avg(flanks)."""
    depth = np.zeros(300)
    depth[150] = 10.0  # isolated peak, zero flanks
    calls = call_nucleosomes(depth, flank=50)
    assert len(calls) == 1
    assert calls[0].peak_depth == 10.0
    assert calls[0].score == 10.0  # flanks are 0, so score == peak exactly


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