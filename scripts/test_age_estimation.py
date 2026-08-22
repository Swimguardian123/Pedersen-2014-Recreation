"""
Reusable demo/smoke-test for age_estimation.py's full two-layer pipeline.

Pure math -- no BAM/FASTA/network needed, so this is fully runnable standalone at
any point. Uses illustrative synthetic Ms/beta values for the 5 "modern donor"
calibration points and a made-up ancient Ms value; replace REPLACE_WITH_* below
with your real Saqqaq Ms values (from methylation.py) and real modern-donor
(Ms, beta) pairs when you have them, or just run as-is to confirm the pipeline
still produces sane output.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.age_estimation import ModernDonor, estimate_age_at_death, KOCH_WAGNER_MODELS

# Synthetic modern-donor calibration data: 5 donors per CpG, Ms values loosely
# correlated with known beta so the fit isn't degenerate. REPLACE with real
# Ms-vs-beta pairs from your own modern comparison samples for a real estimate.
SYNTHETIC_MODERN_DONORS = {
    "cg07533148": [  # TRIM58
        ModernDonor(ms=0.02, beta=0.03),
        ModernDonor(ms=0.05, beta=0.08),
        ModernDonor(ms=0.09, beta=0.15),
        ModernDonor(ms=0.14, beta=0.22),
        ModernDonor(ms=0.20, beta=0.31),
    ],
    "cg01530101": [  # KCNQ1DN
        ModernDonor(ms=0.10, beta=0.18),
        ModernDonor(ms=0.15, beta=0.25),
        ModernDonor(ms=0.22, beta=0.34),
        ModernDonor(ms=0.28, beta=0.41),
        ModernDonor(ms=0.35, beta=0.50),
    ],
}

# Synthetic ancient Ms values (REPLACE with real ms_score() output from
# methylation.py at each CpG's 2000bp window in your actual data)
SYNTHETIC_SAQQAQ_MS = {
    "cg07533148": 0.11,
    "cg01530101": 0.19,
}


def run(saqqaq_ms=None, modern_donors=None) -> None:
    saqqaq_ms = saqqaq_ms or SYNTHETIC_SAQQAQ_MS
    modern_donors = modern_donors or SYNTHETIC_MODERN_DONORS

    print("Registered Koch & Wagner (2011) models:")
    for cpg, m in KOCH_WAGNER_MODELS.items():
        flag = "used by Pedersen" if m["used_by_pedersen"] else "reference only"
        print(f"  {cpg} ({m['gene']}): intercept={m['intercept']:.4f}, "
              f"slope={m['slope']:.4f}  [{flag}]")
    print()

    results = estimate_age_at_death(saqqaq_ms, modern_donors)

    for cpg, r in results.items():
        if cpg == "mean_predicted_age":
            continue
        print(f"{cpg}:")
        print(f"  Ms->beta calibration: intercept={r['ms_to_beta_intercept']:.4f}, "
              f"slope={r['ms_to_beta_slope']:.4f}, R^2={r['ms_to_beta_r2']:.3f}")
        print(f"  predicted beta={r['predicted_beta']:.4f}")
        print(f"  predicted age={r['predicted_age']:.1f} years\n")

    print(f"Mean predicted age at death: {results['mean_predicted_age']:.1f} years")
    print("\n(Using SYNTHETIC calibration/Ms data unless you passed real values --  "
          "this confirms the pipeline runs correctly, not a real age estimate.)")


if __name__ == "__main__":
    run()