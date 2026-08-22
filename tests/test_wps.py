"""
Regression tests for wps.py's core counting logic (wps_from_fragments,
call_peaks_from_wps). Pure numpy, no pysam or BAM needed -- compute_wps()
itself (the pysam-fetching wrapper) is exercised via real data separately,
same split as anchor_profile() elsewhere in this codebase.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.wps import wps_from_fragments, call_peaks_from_wps


def test_wps_hand_computed_scenario():
    """Five fragments with deliberately distinct relationships to a window
    centered at position 100 (window=20, so window=[90,110]):
      (50,150)  fully spans        -> spanning
      (95,200)  starts inside      -> endpoint in window
      (0,105)   ends inside        -> endpoint in window
      (200,300) nowhere near       -> neither
      (91,109)  both ends inside   -> counted once as endpoint (elif branch)
    Expected: spanning=1, endpoint_in_window=3 -> WPS = 1 - 3 = -2."""
    fragments = [(50, 150), (95, 200), (0, 105), (200, 300), (91, 109)]
    wps = wps_from_fragments(fragments, start=95, end=106, window=20)
    assert wps[5] == -2  # position 100 is index 5 (100-95)


def test_wps_single_fully_spanning_fragment():
    """One fragment spanning the entire test region should give WPS=1 at
    every position (1 spanning, 0 endpoints, everywhere)."""
    fragments = [(0, 200)]
    wps = wps_from_fragments(fragments, start=50, end=150, window=20)
    assert np.all(wps == 1)


def test_wps_no_fragments_is_all_zero():
    wps = wps_from_fragments([], start=0, end=100, window=20)
    assert np.all(wps == 0)


def test_peak_caller_filters_by_min_width():
    """A run of positive WPS >= min_peak_width should be called; a shorter
    run should not."""
    synthetic_wps = np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, 1, 1, 0])
    peaks = call_peaks_from_wps(synthetic_wps, min_peak_width=3)
    assert peaks == [(2, 8)]  # only the length-6 run qualifies; length-2 run doesn't


if __name__ == "__main__":
    test_wps_hand_computed_scenario()
    test_wps_single_fully_spanning_fragment()
    test_wps_no_fragments_is_all_zero()
    test_peak_caller_filters_by_min_width()
    print("All wps tests passed.")