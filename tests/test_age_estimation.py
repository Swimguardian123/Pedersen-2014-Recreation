"""
Regression tests for age_estimation.py -- specifically pins the real Koch &
Wagner (2011) Figure 3A constants so a future edit can't silently drift them.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.age_estimation import (
    KOCH_WAGNER_MODELS, ModernDonor, fit_ms_to_beta, ms_to_beta, beta_to_age,
    estimate_age_at_death,
)


def test_koch_wagner_constants_match_figure_3a():
    """Pins the real published values (not placeholders) -- percent-scale
    equation from the figure, converted to 0-1 fraction scale."""
    trim58 = KOCH_WAGNER_MODELS["cg07533148"]
    assert abs(trim58["slope"] - 0.39 / 100) < 1e-9
    assert abs(trim58["intercept"] - (-8.02 / 100)) < 1e-9
    assert trim58["used_by_pedersen"] is True

    kcnq1dn = KOCH_WAGNER_MODELS["cg01530101"]
    assert abs(kcnq1dn["slope"] - 0.57 / 100) < 1e-9
    assert abs(kcnq1dn["intercept"] - 0.0117) < 1e-9
    assert kcnq1dn["used_by_pedersen"] is True


def test_ms_to_beta_perfect_linear_fit():
    """A perfectly linear Ms-vs-beta relationship should fit with R^2 == 1."""
    donors = [ModernDonor(ms=x, beta=2 * x + 0.1) for x in [0.1, 0.2, 0.3, 0.4, 0.5]]
    intercept, slope, r2 = fit_ms_to_beta(donors)
    assert abs(r2 - 1.0) < 1e-9
    assert abs(slope - 2.0) < 1e-6
    assert abs(intercept - 0.1) < 1e-6


def test_beta_to_age_inverts_the_model_correctly():
    """If beta = intercept + slope*age exactly, recovering age from beta should
    return the original age (round-trip correctness of the inversion)."""
    model = KOCH_WAGNER_MODELS["cg07533148"]
    true_age = 65.0
    beta = model["intercept"] + model["slope"] * true_age
    recovered_age = beta_to_age(beta, "cg07533148")
    assert abs(recovered_age - true_age) < 1e-6


def test_full_pipeline_runs_and_returns_sane_range():
    """End-to-end pipeline (same shape as scripts/test_age_estimation.py) should
    return an age in a biologically plausible human range, not something wild."""
    donors = {
        "cg07533148": [ModernDonor(ms=m, beta=1.5 * m) for m in [0.02, 0.05, 0.09, 0.14, 0.20]],
        "cg01530101": [ModernDonor(ms=m, beta=1.3 * m + 0.05) for m in [0.10, 0.15, 0.22, 0.28, 0.35]],
    }
    saqqaq_ms = {"cg07533148": 0.11, "cg01530101": 0.19}

    results = estimate_age_at_death(saqqaq_ms, donors)
    mean_age = results["mean_predicted_age"]
    assert mean_age is not None
    assert 0 < mean_age < 120, f"predicted age {mean_age} outside plausible human range"


if __name__ == "__main__":
    test_koch_wagner_constants_match_figure_3a()
    test_ms_to_beta_perfect_linear_fit()
    test_beta_to_age_inverts_the_model_correctly()
    test_full_pipeline_runs_and_returns_sane_range()
    print("All age_estimation tests passed.")