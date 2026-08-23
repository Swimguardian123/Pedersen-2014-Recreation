"""
Driver script: load the real Horvath coefficient CSV (exported from
methylclockData via R) and verify it parses correctly -- surfaces exactly
what was found, same diagnostic-first pattern as run_deamination_rate.py,
so a column-format mismatch is immediately visible rather than silently wrong.

NOTE: this verifies the coefficient table loads and the formula applies
correctly -- it does NOT produce a real predicted age for the Saqqaq sample.
Horvath's clock needs per-CpG-probe beta values (modern bisulfite-array
style measurements at 353 exact named CpG sites); our aDNA methylation.py
produces a region-aggregated deamination-based Ms score, a different kind of
measurement entirely. This script demonstrates the mechanism is correct,
using synthetic beta values for a concrete example -- not a real age result.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.dnamage_clock import load_coefficients_from_csv, predict_age, trafo, anti_trafo


def run(csv_path: str) -> None:
    print(f"=== Loading {csv_path} ===\n")

    if not os.path.exists(csv_path):
        print(f"File not found: {csv_path}")
        return

    coefficients = load_coefficients_from_csv(csv_path)

    if not coefficients:
        print("WARNING: parsed zero coefficients. This means the real file's column")
        print("names don't match what load_coefficients_from_csv() expects")
        print("('CpGmarker'/'CoefficientTraining', or 'cpg'/'coefficient').")
        print("Run this to see the real header before adjusting the parser:")
        print(f"  head -3 {csv_path}")
        return

    has_intercept = "intercept" in coefficients
    n_cpgs = len(coefficients) - (1 if has_intercept else 0)

    print(f"Parsed {n_cpgs} CpG coefficient(s)")
    print(f"Intercept found: {has_intercept}" +
          (f" (value: {coefficients['intercept']})" if has_intercept else ""))

    sample_cpgs = [c for c in coefficients if c != "intercept"][:5]
    print(f"Sample CpG IDs: {sample_cpgs}")
    print(f"Sample coefficients: {[coefficients[c] for c in sample_cpgs]}\n")

    if n_cpgs < 300:
        print(f"NOTE: expected ~353 CpGs, found {n_cpgs}. Real Horvath clock has 353 --")
        print("if this is far off, double check the R export captured the full table.\n")

    if not has_intercept:
        print("WARNING: no intercept found -- predict_age() will need one passed")
        print("explicitly, since it's not in the parsed coefficients.\n")
        return

    # Demonstrate the mechanism with SYNTHETIC beta values (0.5 for every CpG --
    # not real methylation data, just proves the formula pipeline works end to end)
    print("[Mechanism check] Using synthetic beta=0.5 for every CpG (NOT real data,")
    print("just confirms the linear-combination + anti_trafo pipeline runs correctly):")
    real_cpgs = {c: v for c, v in coefficients.items() if c != "intercept"}
    synthetic_betas = {cpg: 0.5 for cpg in real_cpgs}
    result = predict_age(synthetic_betas, real_cpgs, coefficients["intercept"])
    print(f"  Result: {result:.2f} (a real number came out -- pipeline mechanically works)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Verify Horvath clock coefficient loading.")
    p.add_argument("csv_path", help="Path to the CSV exported from R (methylclockData)")
    args = p.parse_args()

    run(args.csv_path)