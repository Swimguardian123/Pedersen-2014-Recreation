"""
Expected deamination rate per overhang, from mapDamage2 output.

Paper reference (Hanghoj et al. 2016, Methods, "Regional Methylation Levels"):
"The expected number of deaminations per overhang was estimated within CTCF
regions using ds and k parameters as retrieved from mapDamage2 (Jonsson et al.
2013)... and the following formula (ds*((1/k)-1))."

Parameter mapping (confirmed, not guessed): mapDamage2's Bayesian model
estimates four named parameters -- Lambda (average overhang length /
per-position termination probability), DeltaD (deamination rate in
double-stranded regions), DeltaS (deamination rate in single-stranded
overhangs), Theta (background mismatch rate unrelated to damage). Hanghoj's
"ds" = DeltaS, "k" = Lambda. This mapping isn't just name-matching: the
formula ds*((1/k)-1) is literally "deamination rate x expected overhang
length" -- (1/k)-1 is the expected length of a geometrically-distributed
overhang with per-site termination probability k, which is exactly the
quantity mapDamage2's own model assumes. The formula and the parameter
definitions are mutually consistent, not just plausibly matched.

Source file: mapDamage2 writes `Stats_out_MCMC_iter_summ_stat.csv` --
"summary statistics for the damage parameters estimated posterior
distributions" (confirmed from mapDamage2's own documented output file list).
One honest residual gap: the EXACT column header text/casing in that CSV
(e.g. "Lambda" vs "lambda", "DeltaS" vs "Delta_ss") isn't confirmed
byte-for-byte here -- parse_mapdamage_stats() below matches flexibly
(case-insensitive substring match) rather than assuming one exact spelling,
specifically to avoid silently failing on a real file due to a casing
mismatch. Verify against your actual output file if this needs tightening.
"""

import csv
from typing import Dict, Optional, Tuple


def parse_mapdamage_stats_full(csv_path: str) -> Dict[str, Dict[str, float]]:
    """
    Parse Stats_out_MCMC_iter_summ_stat.csv.

    REAL structure, confirmed against an actual mapDamage2 run in this project
    (a previous version of this function assumed the opposite orientation and
    was wrong -- fixed here against real data, not guessed again):
    columns ARE the parameter names (Theta, DeltaD, DeltaS, Lambda, plus Rho/
    LogLik which aren't damage parameters and are ignored here); each ROW is
    one summary statistic, with the stat name as a row label in column 0:
    "Mean", "Std.", "Acceptance ratio", then a full percentile ladder
    ("0%", "2.5%", "5%", ... "97.5%", "100%").

    Returns {param_name: {"mean": ..., "std": ..., "ci_low": ..., "ci_high": ...}}
    (ci_low/ci_high from the "2.5%"/"97.5%" rows specifically, not the full
    percentile ladder -- those two are what expected_deaminations_range() uses).
    """
    target_params = ["theta", "deltad", "deltas", "lambda"]

    with open(csv_path, newline="") as f:
        rows = list(csv.reader(f))

    if not rows:
        return {}

    header_lower = [h.strip().lower() for h in rows[0]]

    # column index -> parameter name, for columns that match our target params
    param_col: Dict[str, int] = {}
    for i, h in enumerate(header_lower):
        if h in target_params:
            param_col[h] = i

    if not param_col:
        return {}

    results: Dict[str, Dict[str, float]] = {p: {} for p in param_col}

    for row in rows[1:]:
        if not row:
            continue
        row_label = row[0].strip().lower()

        if row_label == "mean":
            stat_role = "mean"
        elif row_label.startswith("std"):
            stat_role = "std"
        elif row_label == "2.5%":
            stat_role = "ci_low"
        elif row_label == "97.5%":
            stat_role = "ci_high"
        else:
            continue  # not a stat row we track (e.g. "Acceptance ratio", other percentiles)

        for param, col_idx in param_col.items():
            try:
                results[param][stat_role] = float(row[col_idx])
            except (ValueError, IndexError):
                pass

    return {p: v for p, v in results.items() if v}


def parse_mapdamage_stats(csv_path: str) -> Dict[str, float]:
    """
    Backward-compatible: {parameter_name: mean_value} only. Kept unchanged in
    behavior/signature so existing callers (and the existing regression test)
    don't break -- built on top of parse_mapdamage_stats_full() now, plus a
    fallback for the column-header layout that full-stats parsing doesn't
    attempt (that layout only ever had one row of values, so "full stats" -- 
    std/CI -- don't meaningfully apply to it anyway).
    """
    full = parse_mapdamage_stats_full(csv_path)
    if full:
        return {p: v["mean"] for p, v in full.items() if "mean" in v}

    # fallback: column-header layout (parameter names ARE the headers)
    target_params = ["lambda", "deltad", "deltas", "theta"]
    results = {}
    with open(csv_path, newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return results
    header = [h.strip().lower() for h in rows[0]]
    for i, h in enumerate(header):
        for param in target_params:
            if param in h and len(rows) > 1:
                try:
                    results[param] = float(rows[1][i])
                except (ValueError, IndexError):
                    pass
    return results


def expected_deaminations_per_overhang(delta_s: float, lambda_param: float) -> float:
    """
    Hanghoj et al. 2016's formula: ds * ((1/k) - 1).
    delta_s: DeltaS (single-strand/overhang deamination rate per site).
    lambda_param: Lambda (overhang-length/termination-probability parameter).
    """
    if lambda_param <= 0:
        raise ValueError("lambda_param must be > 0 (it's a probability)")
    return delta_s * ((1.0 / lambda_param) - 1.0)


def expected_deaminations_range(deltas_stats: Dict[str, float],
                                 lambda_stats: Dict[str, float]
                                 ) -> Tuple[float, Optional[float], Optional[float]]:
    """
    Point estimate (from each parameter's mean) plus a naive sanity range using
    each parameter's own 2.5%/97.5% credible interval bounds.

    Numerically confirmed before implementing: f(ds,k) = ds*(1/k - 1) is
    monotonically INCREASING in ds and DECREASING in k for ds>0, k in (0,1) --
    so the range's low/high bounds come from the diagonal corners, not all four:
        low  = f(ds_ci_low,  k_ci_high)
        high = f(ds_ci_high, k_ci_low)

    Honest limitation, stated plainly rather than left implicit: this treats ds
    and k as independent, but they're jointly estimated by the same MCMC fit and
    are very likely correlated in the real posterior. Plugging in each
    parameter's own marginal CI bound independently like this is a common,
    useful sanity-check technique, but it is NOT a rigorous propagated joint
    confidence interval, and most likely overstates the true uncertainty width.
    Treat this as "here's a plausible spread," not a statistically rigorous CI.

    Returns (point_estimate, naive_low, naive_high) -- the latter two are None
    if the CI fields weren't found in the parsed stats (e.g. only mean was
    available).
    """
    point = expected_deaminations_per_overhang(deltas_stats["mean"], lambda_stats["mean"])

    required = ("ci_low", "ci_high")
    if not all(k in deltas_stats for k in required) or not all(k in lambda_stats for k in required):
        return point, None, None

    low = expected_deaminations_per_overhang(deltas_stats["ci_low"], lambda_stats["ci_high"])
    high = expected_deaminations_per_overhang(deltas_stats["ci_high"], lambda_stats["ci_low"])
    return point, low, high


def expected_deaminations_from_mapdamage_output(csv_path: str) -> Optional[float]:
    """Convenience: parse the file and compute the formula in one call."""
    stats = parse_mapdamage_stats(csv_path)
    if "deltas" not in stats or "lambda" not in stats:
        missing = [p for p in ("deltas", "lambda") if p not in stats]
        raise ValueError(f"Could not find required parameter(s) in {csv_path}: {missing}. "
                          f"Found: {list(stats.keys())}. Check the file's actual column "
                          f"layout against parse_mapdamage_stats()'s assumptions.")
    return expected_deaminations_per_overhang(stats["deltas"], stats["lambda"])


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Compute expected deaminations per overhang from mapDamage2 output.")
    p.add_argument("stats_csv", help="Path to Stats_out_MCMC_iter_summ_stat.csv")
    args = p.parse_args()

    full_stats = parse_mapdamage_stats_full(args.stats_csv)
    print("Parsed parameters (full stats):")
    for param, vals in full_stats.items():
        print(f"  {param}: {vals}")

    if "deltas" in full_stats and "lambda" in full_stats:
        point, low, high = expected_deaminations_range(full_stats["deltas"], full_stats["lambda"])
        print(f"\nExpected deaminations per overhang: {point:.4f}")
        if low is not None:
            print(f"  naive sanity range (see docstring caveat): [{low:.4f}, {high:.4f}]")
    else:
        print("Could not find DeltaS and/or Lambda -- check the file's actual layout "
              "against parse_mapdamage_stats_full()'s assumptions.")