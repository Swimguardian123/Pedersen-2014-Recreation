"""
Orchestrator: runs every pipeline stage that can run from just a BAM+FASTA+
region, in sequence, with a consolidated summary at the end. Stages needing
an external prerequisite you generate separately (mapDamage2 output, real
Horvath coefficients) are clearly reported as skipped with instructions,
rather than silently omitted or falsely claimed to run.

This does NOT replace the individual scripts -- each remains independently
runnable for focused debugging/iteration, same as throughout this project.
This is a convenience wrapper on top, not a new source of truth.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.gc_correction import build_gc_rate_model, gc_correct_depth
from src.nucleosome_calling import call_nucleosomes
from src.periodicity import welch_periodogram, dominant_period, phasogram, phasogram_dominant_periods
from src.nucleotide_patterns import mononucleotide_matrix, purine_pyrimidine_dinucleotide_signal
from src.methylation import ms_score
from src.wps import compute_wps, call_peaks_from_wps
from src.mappability_filter import download_mappability_track, get_high_mappability_regions, filter_calls_by_mappability
from src.deamination_rate import parse_mapdamage_stats_full, expected_deaminations_range
from src.dnamage_clock import load_coefficients_from_csv, predict_age


def section(title):
    print(f"\n{'=' * 60}\n{title}\n{'=' * 60}")


def run(bam, fasta, chrom, start, end, mapdamage_csv=None, dnamage_csv=None,
        skip_mappability=False, mappability_kmer=40):
    results = {}

    section("[1/7] Nucleosome occupancy (depth -> GC-correct -> call -> periodicity)")
    try:
        gc_model = build_gc_rate_model(bam, fasta, chrom, start, end)
        corrected_depth = gc_correct_depth(bam, fasta, chrom, start, end, gc_model)
        calls = call_nucleosomes(corrected_depth)
        period = dominant_period(corrected_depth)
        print(f"  {len(calls)} nucleosome call(s), dominant period: {period}")
        results["occupancy"] = {"n_calls": len(calls), "dominant_period": period}
    except Exception as e:
        print(f"  FAILED: {e}")
        results["occupancy"] = None

    section("[2/7] Methylation (Ms score)")
    try:
        rate, hits, total = ms_score(bam, fasta, chrom, start, end)
        print(f"  Ms = {rate} ({hits}/{total})")
        results["methylation"] = {"rate": rate, "hits": hits, "total": total}
    except Exception as e:
        print(f"  FAILED: {e}")
        results["methylation"] = None

    section("[3/7] WPS vs NucleoMap comparison")
    try:
        wps = compute_wps(bam, chrom, start, end, window=50)
        wps_peaks = call_peaks_from_wps(wps, min_peak_width=3)
        print(f"  WPS: {len(wps_peaks)} peak region(s); NucleoMap: "
              f"{results['occupancy']['n_calls'] if results['occupancy'] else '?'} call(s)")
        results["wps"] = {"n_peaks": len(wps_peaks)}
    except Exception as e:
        print(f"  FAILED: {e}")
        results["wps"] = None

    section("[4/7] Mappability filtering")
    if skip_mappability:
        print("  SKIPPED (--skip-mappability set -- real ~1-2GB download otherwise)")
        results["mappability"] = None
    else:
        try:
            track_path = download_mappability_track(kmer=mappability_kmer)
            regions = get_high_mappability_regions(track_path, chrom, start, end)
            filtered = filter_calls_by_mappability(calls, start, regions) if results["occupancy"] else []
            print(f"  {len(regions)} high-mappability block(s); "
                  f"{len(filtered)}/{results['occupancy']['n_calls'] if results['occupancy'] else '?'} call(s) survive")
            results["mappability"] = {"n_high_mappability_blocks": len(regions), "n_calls_surviving": len(filtered)}
        except Exception as e:
            print(f"  FAILED: {e}")
            results["mappability"] = None

    section("[5/7] Deamination rate (needs real mapDamage2 output)")
    if mapdamage_csv is None:
        print("  SKIPPED -- no --mapdamage-csv given.")
        print("  To run this stage: run mapDamage on your BAM, then pass")
        print("  --mapdamage-csv results_.../Stats_out_MCMC_iter_summ_stat.csv")
        results["deamination_rate"] = None
    else:
        try:
            stats = parse_mapdamage_stats_full(mapdamage_csv)
            point, low, high = expected_deaminations_range(stats["deltas"], stats["lambda"])
            print(f"  Expected deaminations/overhang: {point:.4f}"
                  + (f"  [{low:.4f}, {high:.4f}]" if low is not None else ""))
            results["deamination_rate"] = {"point_estimate": point}
        except Exception as e:
            print(f"  FAILED: {e}")
            results["deamination_rate"] = None

    section("[6/7] DNAmAge clock (needs real Horvath coefficients CSV)")
    if dnamage_csv is None:
        print("  SKIPPED -- no --dnamage-csv given.")
        print("  To run this stage: export Horvath coefficients from R (methylclockData),")
        print("  then pass --dnamage-csv horvath_coefficients.csv")
        results["dnamage"] = None
    else:
        try:
            coefficients = load_coefficients_from_csv(dnamage_csv)
            print(f"  Loaded {len(coefficients) - 1} CpG coefficient(s) + intercept.")
            print("  (Mechanism-only -- a real age estimate needs per-CpG-probe beta values")
            print("   this aDNA pipeline doesn't produce; see PROJECT_SUMMARY.md.)")
            results["dnamage"] = {"n_cpgs": len(coefficients) - 1}
        except Exception as e:
            print(f"  FAILED: {e}")
            results["dnamage"] = None

    section("[7/7] Summary")
    for stage, result in results.items():
        status = "ran" if result is not None else "skipped/failed"
        print(f"  {stage:20s} {status}")

    return results


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Run every runnable pipeline stage against real data.")
    p.add_argument("bam")
    p.add_argument("fasta")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--mapdamage-csv", default=None,
                    help="Path to a real mapDamage2 Stats_out_MCMC_iter_summ_stat.csv")
    p.add_argument("--dnamage-csv", default=None,
                    help="Path to a real Horvath coefficients CSV (exported via R)")
    p.add_argument("--skip-mappability", action="store_true",
                    help="Skip the mappability stage (real ~1-2GB download otherwise)")
    p.add_argument("--mappability-kmer", type=int, default=40)
    args = p.parse_args()

    run(args.bam, args.fasta, args.chrom, args.start, args.end,
        args.mapdamage_csv, args.dnamage_csv, args.skip_mappability, args.mappability_kmer)