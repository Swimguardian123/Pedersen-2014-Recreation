"""
Age-at-death estimation from ancient methylation levels.

Paper reference (Methods, "Methylation signal", final lines; Results, "Prediction of
the age at death"):
"Estimates of age at death were derived using linear models from the literature
(Koch and Wagner 2011) that relate age and methylation levels at given CpG sites."
"...we focused on four particular CpGs for which age-methylation linear models have
been established... Two CpGs (cg23571857 and cg25148589) showed large differences
between predicted and real age... and were disregarded. However, two other CpGs
(cg07533148 and cg01530101) provided reliable age estimates... Estimating Ms for a
2000-bp-wide region centered on each CpG from the Illumina 450k array, we built a
linear model relating Ms and the methylation levels measured in hairs of five modern
donors... in order to convert Ms into absolute methylation levels for Saqqaq at the
two loci of interest... the absolute levels were used to infer age."

So the pipeline has TWO linear-model layers, not one:
  1. Ms (our aDNA misincorporation-based proxy) -> absolute beta-value methylation
     level, calibrated using Ms + known beta-values from ~5 modern hair donors
     (Slieker et al. 2013) at the SAME two CpGs. This calibration is something YOU
     fit yourself from your own modern comparison data (Pedersen's Supplemental
     Table S3.1 R^2 range was 0.620-0.785 -- expect similar).
  2. Beta-value -> age, using Koch & Wagner (2011)'s published per-CpG linear model
     (beta = A + B*age, inverted to age = (beta - A) / B).

Source confirmed (Koch & Wagner 2011, "Epigenetic-aging-signature to determine age
in different tissues", Aging 3:1018-1027): the two CpGs Pedersen kept are
cg07533148 (TRIM58) and cg01530101 (KCNQ1DN). Figure 3A prints each CpG's regression
line directly on the scatter plot as "y = SLOPE*x + INTERCEPT", where y is CpG
methylation on a 0-100% axis and x is age in years. All five panels' equations:
    TRIM58   (cg07533148): y = 0.39x -  8.02
    KCNQ1DN  (cg01530101): y = 0.57x +  1.17
    NPTX2    (cg12799895): y = 0.39x -  3.36
    GRIA2    (cg25148589): y = 0.32x +  9.27
    BIRC4BP  (cg23571857): y = -0.47x + 76.62
Pedersen only used TRIM58 and KCNQ1DN (the other three -- NPTX2, GRIA2, BIRC4BP --
were part of Koch & Wagner's 5-CpG panel but Pedersen's own modern-donor validation
found only these two gave reliable predictions on hair specifically); all five are
recorded below for completeness/reference.

Unit conversion: the figure's y is percent (0-100 scale, matching the "100%...0%"
axis labels), but Ms and Illumina beta-values conventionally run 0-1. All constants
below are stored pre-converted to the 0-1 fraction scale (slope/100, intercept/100)
so beta_to_age() takes and expects a plain 0-1 beta-value, consistent with the rest
of this pipeline and with ms_to_beta()'s output.
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

# --- Layer 2: Koch & Wagner (2011) per-CpG beta-vs-age linear models -------
# beta (0-1 fraction) = INTERCEPT + SLOPE * age  =>  age = (beta - INTERCEPT) / SLOPE
# Source: Figure 3A regression equations, converted from percent to fraction scale.
KOCH_WAGNER_MODELS: Dict[str, Dict[str, float]] = {
    "cg07533148": {"gene": "TRIM58",  "intercept": -8.02 / 100, "slope": 0.39 / 100,
                   "used_by_pedersen": True},
    "cg01530101": {"gene": "KCNQ1DN", "intercept": 1.17 / 100,  "slope": 0.57 / 100,
                   "used_by_pedersen": True},
    "cg12799895": {"gene": "NPTX2",   "intercept": -3.36 / 100, "slope": 0.39 / 100,
                   "used_by_pedersen": False},
    "cg25148589": {"gene": "GRIA2",   "intercept": 9.27 / 100,  "slope": 0.32 / 100,
                   "used_by_pedersen": False},
    "cg23571857": {"gene": "BIRC4BP", "intercept": 76.62 / 100, "slope": -0.47 / 100,
                   "used_by_pedersen": False},
}


@dataclass
class ModernDonor:
    """One modern comparison sample used to calibrate Ms -> beta (Layer 1)."""
    ms: float    # Ms score computed the same way as for Saqqaq, same 2000bp window
    beta: float  # known methylation beta-value at that CpG (e.g. from Slieker et al. 2013)


def fit_ms_to_beta(donors: List[ModernDonor]) -> Tuple[float, float, float]:
    """
    Layer 1: linear regression of beta (y) on Ms (x) across modern donors, per CpG.
    Returns (intercept, slope, r_squared).
    """
    ms = np.array([d.ms for d in donors])
    beta = np.array([d.beta for d in donors])
    if len(ms) < 2:
        raise ValueError("Need at least 2 modern donors to fit Ms->beta calibration")

    slope, intercept = np.polyfit(ms, beta, 1)
    pred = intercept + slope * ms
    ss_res = np.sum((beta - pred) ** 2)
    ss_tot = np.sum((beta - beta.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return float(intercept), float(slope), float(r2)


def ms_to_beta(ms_value: float, ms_beta_intercept: float, ms_beta_slope: float) -> float:
    """Apply the fitted Layer-1 model to convert an ancient Ms score to a beta-value."""
    return ms_beta_intercept + ms_beta_slope * ms_value


def beta_to_age(beta_value: float, cpg_id: str) -> float:
    """
    Layer 2: invert Koch & Wagner's published beta-vs-age model for the given CpG.
    beta_value is expected on the 0-1 fraction scale (see module docstring).
    """
    model = KOCH_WAGNER_MODELS.get(cpg_id)
    if model is None:
        raise KeyError(f"No Koch & Wagner model registered for {cpg_id}")
    return (beta_value - model["intercept"]) / model["slope"]


def estimate_age_at_death(saqqaq_ms_by_cpg: Dict[str, float],
                           modern_donors_by_cpg: Dict[str, List[ModernDonor]]
                           ) -> Dict[str, dict]:
    """
    Full two-layer pipeline for each CpG: Ms -> beta (fit per-CpG from modern donors)
    -> age (Koch & Wagner model). Mirrors the paper's approach of predicting age
    independently per CpG, then averaging (paper: "Both CpGs considered provided
    strikingly similar age estimates... mean of the predictions").

    Returns per-CpG dict with intermediate values, plus a combined mean estimate.
    """
    results = {}
    ages = []

    for cpg_id, ms_value in saqqaq_ms_by_cpg.items():
        donors = modern_donors_by_cpg.get(cpg_id, [])
        intercept, slope, r2 = fit_ms_to_beta(donors)
        beta = ms_to_beta(ms_value, intercept, slope)
        age = beta_to_age(beta, cpg_id)
        results[cpg_id] = {
            "ms_to_beta_intercept": intercept,
            "ms_to_beta_slope": slope,
            "ms_to_beta_r2": r2,
            "predicted_beta": beta,
            "predicted_age": age,
        }
        ages.append(age)

    results["mean_predicted_age"] = float(np.mean(ages)) if ages else None
    return results


if __name__ == "__main__":
    print("Registered Koch & Wagner (2011) CpG models (beta-fraction scale):")
    for cpg, m in KOCH_WAGNER_MODELS.items():
        flag = "used by Pedersen" if m["used_by_pedersen"] else "not used by Pedersen"
        print(f"  {cpg} ({m['gene']}): intercept={m['intercept']:.4f}, "
              f"slope={m['slope']:.4f}  [{flag}]")