"""
Regression tests for expression_inference.py's pure-math functions (compute_rs,
first_nucleosome_occupancy, phasing_strength, rank_into_quantiles, evaluate_proxy).
No pysam or network needed -- fetch_gse3058_expression is deliberately not
covered here (needs live GEO access, same reasoning as other network-dependent
functions in this project).
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.expression_inference import (
    compute_rs, first_nucleosome_occupancy, phasing_strength,
    rank_into_quantiles, evaluate_proxy,
)


def test_compute_rs_basic_ratio():
    assert abs(compute_rs(genebody_ms=0.4, promoter_ms=0.2) - 2.0) < 1e-4


def test_compute_rs_pseudocount_avoids_divide_by_zero():
    # promoter_ms=0 would crash on a naive division; pseudocount should prevent that
    result = compute_rs(genebody_ms=0.5, promoter_ms=0.0)
    assert np.isfinite(result)
    assert result > 0


def test_first_nucleosome_occupancy_is_a_windowed_mean():
    depth = np.array([0, 10, 20, 30, 0, 0, 0])
    result = first_nucleosome_occupancy(depth, tss_offset=1, window=3)
    assert abs(result - np.mean([10, 20, 30])) < 1e-9


def test_phasing_strength_higher_for_matching_period_than_random_noise():
    """A signal built with real 193bp periodicity should show more power at
    that target period than pure random noise of the same scale."""
    np.random.seed(0)
    n = 3000
    x = np.arange(n)
    periodic = 50 + 40 * np.sin(2 * np.pi * x / 193) + np.random.normal(0, 3, n)
    noise = np.random.normal(50, 40, n)
    noise = np.clip(noise, 0, None)

    periodic_strength = phasing_strength(np.clip(periodic, 0, None))
    noise_strength = phasing_strength(noise)

    assert periodic_strength > noise_strength


def test_rank_into_quantiles_assigns_extremes_correctly():
    expression = {f"gene{i}": float(i) for i in range(10)}  # gene0=0.0 ... gene9=9.0
    groups = rank_into_quantiles(expression, n_quantiles=5)
    assert groups["gene0"] == 0   # lowest expression -> lowest bin
    assert groups["gene9"] == 4   # highest expression -> highest bin
    assert len(set(groups.values())) == 5  # all 5 bins actually used


def test_evaluate_proxy_perfect_monotonic_relationship():
    proxy = {f"g{i}": float(i) for i in range(20)}
    groups = {f"g{i}": i // 4 for i in range(20)}  # monotonically increasing groups
    rho, p = evaluate_proxy(proxy, groups)
    assert rho > 0.95, f"expected near-perfect positive correlation, got rho={rho}"


def test_evaluate_proxy_too_few_common_genes_returns_nan():
    rho, p = evaluate_proxy({"a": 1.0, "b": 2.0}, {"c": 0, "d": 1})  # no overlap at all
    assert np.isnan(rho) and np.isnan(p)


if __name__ == "__main__":
    test_compute_rs_basic_ratio()
    test_compute_rs_pseudocount_avoids_divide_by_zero()
    test_first_nucleosome_occupancy_is_a_windowed_mean()
    test_phasing_strength_higher_for_matching_period_than_random_noise()
    test_rank_into_quantiles_assigns_extremes_correctly()
    test_evaluate_proxy_perfect_monotonic_relationship()
    test_evaluate_proxy_too_few_common_genes_returns_nan()
    print("All expression_inference tests passed.")