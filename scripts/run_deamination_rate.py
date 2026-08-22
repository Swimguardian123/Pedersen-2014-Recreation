"""
Driver script: run deamination_rate.py against a REAL mapDamage2 output file.

This is the actual test of parse_mapdamage_stats()'s flexible column-matching
-- it was written defensively (case-insensitive, tries two plausible layouts)
because the exact column format wasn't independently confirmed, only the file
name and the parameter names it should contain. This script surfaces exactly
what was found vs. expected, so a format mismatch is diagnosable rather than
a silent wrong number or a bare crash.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.deamination_rate import parse_mapdamage_stats_full, expected_deaminations_range


def run(stats_csv_path: str) -> None:
    print(f"=== Parsing {stats_csv_path} ===\n")

    if not os.path.exists(stats_csv_path):
        print(f"File not found: {stats_csv_path}")
        print("This should be mapDamage2's Stats_out_MCMC_iter_summ_stat.csv, found in")
        print("mapDamage2's results folder after running it on a real BAM.")
        return

    full_stats = parse_mapdamage_stats_full(stats_csv_path)
    print(f"Parameters found: {list(full_stats.keys())}")
    for name, vals in full_stats.items():
        print(f"  {name}: {vals}")
    print()

    expected = ["lambda", "deltad", "deltas", "theta"]
    missing = [p for p in expected if p not in full_stats]
    if missing:
        print(f"WARNING: expected parameters not found: {missing}")
        print("This means the layout-guessing didn't match this file's actual structure --")
        print("open the CSV directly and compare its real column headers/layout against")
        print("what parse_mapdamage_stats_full() tries (see its docstring). This is exactly")
        print("the kind of format mismatch this script exists to catch, rather than")
        print("silently producing a wrong number.")
        return

    if "deltad" in full_stats and "deltas" in full_stats:
        deltad_mean = full_stats["deltad"].get("mean")
        deltas_mean = full_stats["deltas"].get("mean")
        if deltad_mean and deltas_mean:
            ratio = deltas_mean / deltad_mean
            print(f"[Sanity check] DeltaS/DeltaD ratio: {ratio:.1f}x")
            print("  (real aDNA damage signature: single-strand overhangs deaminate much")
            print("   faster than double-stranded interior DNA -- a large ratio here is a")
            print("   good sign the fit found genuine damage signal, not noise)\n")

    print("[Formula] Expected deaminations per overhang: ds * ((1/k) - 1)")
    point, low, high = expected_deaminations_range(full_stats["deltas"], full_stats["lambda"])
    print(f"  DeltaS (ds) mean = {full_stats['deltas']['mean']}")
    print(f"  Lambda (k)  mean = {full_stats['lambda']['mean']}")
    print(f"  Point estimate   = {point:.4f} expected deaminations per overhang")
    if low is not None:
        print(f"  Naive sanity range = [{low:.4f}, {high:.4f}]")
        print("  (NOT a rigorous joint confidence interval -- ds and k are jointly")
        print("   estimated and likely correlated in the real posterior; this plugs in")
        print("   each parameter's own marginal 2.5%/97.5% bound independently, which")
        print("   is a useful sanity check but probably overstates true uncertainty.)")
    print()
    print("(Per Hanghoj et al. 2016: this value correlates with the minimal read count")
    print(" needed to recover a real out-of-phase CTCF methylation/occupancy signal --")
    print(" higher deamination rate needs fewer reads to show the pattern. A single")
    print(" sample's number alone doesn't test that correlation, just computes the")
    print(" input to it correctly.)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Compute deamination rate from real mapDamage2 output.")
    p.add_argument("stats_csv", help="Path to mapDamage2's Stats_out_MCMC_iter_summ_stat.csv")
    args = p.parse_args()

    run(args.stats_csv)