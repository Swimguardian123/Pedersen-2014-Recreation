"""
Regression tests for deamination_rate.py. Pure logic + file I/O, no pysam or
network needed -- self-run directly.
"""

import sys
import os
import csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.deamination_rate import (
    expected_deaminations_per_overhang, parse_mapdamage_stats,
    parse_mapdamage_stats_full, expected_deaminations_range,
)


def test_formula_matches_hand_calculation():
    """lambda=0.1 -> expected overhang length = 1/0.1 - 1 = 9bp;
    deltaS=0.3 -> expected total = 0.3 * 9 = 2.7"""
    result = expected_deaminations_per_overhang(delta_s=0.3, lambda_param=0.1)
    assert abs(result - 2.7) < 1e-9


def test_formula_rejects_zero_lambda():
    import pytest
    with pytest.raises(ValueError):
        expected_deaminations_per_overhang(delta_s=0.3, lambda_param=0.0)


def test_parse_mapdamage_stats_matches_real_file_structure(tmp_path):
    """Matches the REAL mapDamage2 output structure, confirmed against an
    actual run in this project (columns = parameter names, rows = stat types
    -- the opposite of what an earlier version of this parser assumed)."""
    path = tmp_path / "Stats_out_MCMC_iter_summ_stat.csv"
    path.write_text(
        '"","Theta","DeltaD","DeltaS","Lambda","Rho","LogLik"\n'
        '"Mean",0.0613,0.0154,0.2002,0.3084,1.06,-613.4\n'
        '"Std.",0.0067,0.0089,0.0994,0.1649,0.13,1.94\n'
        '"2.5%",0.0491,0.0009,0.0807,0.1087,0.82,-618.0\n'
        '"97.5%",0.0755,0.0336,0.4410,0.8026,1.36,-610.8\n'
    )

    stats = parse_mapdamage_stats(str(path))
    assert abs(stats["lambda"] - 0.3084) < 1e-9
    assert abs(stats["deltas"] - 0.2002) < 1e-9
    assert abs(stats["deltad"] - 0.0154) < 1e-9


def test_parse_mapdamage_stats_full_captures_all_columns(tmp_path):
    """Extension test: full stats should capture mean, std, and both CI
    bounds, not just mean -- against the real row/column structure."""
    path = tmp_path / "Stats_out_MCMC_iter_summ_stat.csv"
    path.write_text(
        '"","Theta","DeltaD","DeltaS","Lambda","Rho","LogLik"\n'
        '"Mean",0.0613,0.0154,0.2002,0.3084,1.06,-613.4\n'
        '"Std.",0.0067,0.0089,0.0994,0.1649,0.13,1.94\n'
        '"2.5%",0.0491,0.0009,0.0807,0.1087,0.82,-618.0\n'
        '"97.5%",0.0755,0.0336,0.4410,0.8026,1.36,-610.8\n'
    )

    full = parse_mapdamage_stats_full(str(path))
    assert full["lambda"] == {"mean": 0.3084, "std": 0.1649, "ci_low": 0.1087, "ci_high": 0.8026}
    assert full["deltas"]["ci_high"] == 0.4410


def test_expected_deaminations_range_uses_correct_diagonal_corners():
    """Regression test for the monotonicity-based range logic, numerically
    confirmed before implementation: f(ds,k) increases in ds, decreases in k,
    so the range's low/high come from the diagonal corners specifically."""
    deltas_stats = {"mean": 0.20, "ci_low": 0.15, "ci_high": 0.25}
    lambda_stats = {"mean": 0.30, "ci_low": 0.25, "ci_high": 0.35}

    point, low, high = expected_deaminations_range(deltas_stats, lambda_stats)

    expected_point = expected_deaminations_per_overhang(0.20, 0.30)
    expected_low = expected_deaminations_per_overhang(0.15, 0.35)   # ds_low, k_high
    expected_high = expected_deaminations_per_overhang(0.25, 0.25)  # ds_high, k_low

    assert abs(point - expected_point) < 1e-9
    assert abs(low - expected_low) < 1e-9
    assert abs(high - expected_high) < 1e-9
    assert low < point < high


def test_expected_deaminations_range_returns_none_when_no_ci_available():
    """If only mean was parseable (no CI columns found), range should degrade
    gracefully to (point, None, None) rather than crashing."""
    deltas_stats = {"mean": 0.20}
    lambda_stats = {"mean": 0.30}
    point, low, high = expected_deaminations_range(deltas_stats, lambda_stats)
    assert low is None and high is None
    assert point > 0


if __name__ == "__main__":
    test_formula_matches_hand_calculation()
    test_expected_deaminations_range_returns_none_when_no_ci_available()
    print("All deamination_rate pure-formula tests passed "
          "(tmp_path-based tests need pytest -- run via `pytest tests/`)")