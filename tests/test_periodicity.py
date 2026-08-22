"""
Regression tests for periodicity.py -- specifically encoding the two real bugs
found and fixed during this project (nperseg resolution, min_nonzero NaN gate),
so they can't silently regress.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.periodicity import dominant_period, phasogram, phasogram_dominant_periods


def test_dominant_period_finds_known_periodicity_within_tolerance():
    """Regression test for the nperseg=256 bug: a signal built with exact 193bp
    periodicity should be detected within a reasonable tolerance, not snapped to
    an unrelated value like 256bp (which is what the pre-fix default did)."""
    np.random.seed(0)
    n = 5000
    x = np.arange(n)
    signal = 50 + 40 * np.sin(2 * np.pi * x / 193) + np.random.normal(0, 3, n)
    signal = np.clip(signal, 0, None)

    period = dominant_period(signal)
    assert abs(period - 193) < 25, f"expected ~193bp, got {period}bp"


def test_sparse_signal_returns_nan_not_a_fake_number():
    """Regression test for the min_nonzero gap: a signal with almost no real
    coverage should return NaN, not an arbitrary noise-floor 'answer'."""
    signal = np.zeros(5000)
    signal[100] = 1  # single isolated nonzero point
    signal[2000] = 1
    period = dominant_period(signal)
    assert np.isnan(period), f"expected NaN on near-empty signal, got {period}"


def test_dense_signal_does_not_return_nan():
    """Sanity check the min_nonzero gate doesn't over-trigger on real signal."""
    np.random.seed(0)
    signal = np.random.uniform(5, 15, 500)  # dense, all-nonzero
    period = dominant_period(signal)
    assert not np.isnan(period)


def test_phasogram_detects_known_spacing():
    """Positions spaced exactly 200bp apart, each with sufficient depth (5 reads
    per position, matching the min_depth default), should produce a phasogram
    with a long-range peak near 200bp."""
    unique_positions = np.arange(0, 10000, 200)  # perfectly regular 200bp spacing
    # simulate 5 reads at each position (min_depth default) by repeating each
    # position 5 times, matching the RAW-per-read input the function now expects
    raw_positions = np.repeat(unique_positions, 5)
    hist = phasogram(raw_positions, min_depth=5, max_distance=1000)
    long_range, short_range = phasogram_dominant_periods(hist)
    assert abs(long_range - 200) < 20, f"expected ~200bp long-range peak, got {long_range}"


def test_phasogram_min_depth_actually_filters():
    """Regression test for the dead-parameter bug: positions below min_depth
    must be excluded, not silently included regardless of the parameter."""
    # position 500 has only 2 reads (below min_depth=5) -> should be excluded
    # position 1000 has 6 reads (above min_depth=5) -> should be included
    raw_positions = np.concatenate([
        np.full(2, 500),
        np.full(6, 1000),
    ])
    hist_filtered = phasogram(raw_positions, min_depth=5, max_distance=2000)
    hist_unfiltered = phasogram(raw_positions, min_depth=1, max_distance=2000)

    # with min_depth=5, only position 1000 survives -> no pairwise distances at
    # all exist (only one qualifying position), so the histogram should be empty
    assert hist_filtered.sum() == 0, (
        "expected an empty histogram (only one position passes min_depth=5, "
        "so there are no pairs) -- if this fails, min_depth may not be filtering"
    )
    # with min_depth=1, both positions survive -> at least one real pairwise
    # distance (500) should show up
    assert hist_unfiltered.sum() > 0, "expected real pairwise distances once both positions pass min_depth=1"


if __name__ == "__main__":
    test_dominant_period_finds_known_periodicity_within_tolerance()
    test_sparse_signal_returns_nan_not_a_fake_number()
    test_dense_signal_does_not_return_nan()
    test_phasogram_detects_known_spacing()
    test_phasogram_min_depth_actually_filters()
    print("All periodicity tests passed.")