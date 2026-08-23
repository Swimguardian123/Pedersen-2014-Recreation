"""
Horvath 2013 DNAmAge clock (353-CpG multi-tissue epigenetic age predictor).

A DIFFERENT, larger clock from age_estimation.py's Koch & Wagner (2011) model
(2 CpGs) -- kept as its own module rather than merged in, since it's a
genuinely separate external model with its own data source and math, not a
refinement of the same one.

Source, confirmed (not recollected) from real, independently-maintained
source code -- two separate active bioinformatics R packages (wateRmelon:
github.com/schalkwyk/wateRmelon/blob/master/R/horv.R; sesame:
rdrr.io/bioc/sesame/src/R/age.R) give byte-identical formulas for the
non-linear age transform, cross-confirming each other:

    trafo(age, adult_age=20):       # forward: age -> model response
        x = (age + 1) / (adult_age + 1)
        return log(x) if x <= 1 else x - 1

    anti_trafo(response, adult_age=20):   # inverse: response -> age (this is what we need)
        if response < 0: return (1 + adult_age) * exp(response) - 1
        else:            return (1 + adult_age) * response + adult_age

Full prediction pipeline (from wateRmelon's agep() function):
    linear_response = sum(coefficient[cpg] * beta[cpg] for cpg in the 353 CpGs) + intercept
    age = anti_trafo(linear_response, adult_age=20)

Coefficient table (the 353 CpG IDs + their weights + intercept): NOT YET
SOURCED into this module. Confirmed real location: methylclockData
(Bioconductor ExperimentHub, entry EH6071, "Coefficients Horvath's clock"),
an actual RDA data file at github.com/isglobal-brge/methylclock/blob/master/data.
Since this project already has a working R environment (installed for
mapDamage2), the most direct path is exporting that RDA to a plain CSV from R
directly -- see load_coefficients_from_csv() below for the expected format
once you have that file.
"""

import math
from typing import Dict, Optional


def trafo(age: float, adult_age: float = 20) -> float:
    """Forward transform: age -> model response space."""
    x = (age + 1) / (adult_age + 1)
    return math.log(x) if x <= 1 else x - 1


def anti_trafo(response: float, adult_age: float = 20) -> float:
    """Inverse transform: model response -> age. This is what turns a raw
    elastic-net linear prediction into an actual predicted age."""
    if response < 0:
        return (1 + adult_age) * math.exp(response) - 1
    else:
        return (1 + adult_age) * response + adult_age


def predict_age(betas: Dict[str, float], coefficients: Dict[str, float],
                 intercept: float, adult_age: float = 20,
                 require_all_cpgs: bool = False) -> float:
    """
    Full Horvath clock prediction: linear combination of available CpG betas
    (matched against the coefficient table by CpG ID) + intercept, then
    anti_trafo to get the actual age.

    betas: {cpg_id: beta_value} for whatever CpGs you have measured.
    coefficients: {cpg_id: weight} for the clock's 353 CpGs (missing CpGs in
        `betas` are simply skipped, matching wateRmelon's own behavior of
        working with whatever subset of the 353 is actually available --
        unless require_all_cpgs=True, in which case missing CpGs raise).
    """
    if require_all_cpgs:
        missing = set(coefficients) - set(betas)
        if missing:
            raise ValueError(f"Missing {len(missing)} required CpG(s): {sorted(missing)[:5]}...")

    linear_response = intercept
    n_used = 0
    for cpg_id, coef in coefficients.items():
        if cpg_id in betas:
            linear_response += coef * betas[cpg_id]
            n_used += 1

    if n_used == 0:
        raise ValueError("No overlapping CpGs between betas and coefficients -- "
                          "check that CpG IDs match (e.g. 'cg00000029' format)")

    return anti_trafo(linear_response, adult_age=adult_age)


def load_coefficients_from_csv(csv_path: str) -> Dict[str, float]:
    """
    Load a {cpg_id: coefficient} table from a plain CSV, expected format:
        CpGmarker,CoefficientTraining
        (Intercept),0.696
        cg00075967,0.0123
        cg00374717,-0.0056
        ...
    (This is the standard column naming Horvath's own supplementary data and
    most R-package exports use.) The intercept, if present as a row, is
    returned under the key "intercept" rather than as a regular CpG.
    """
    import csv as csv_module

    coefficients = {}
    intercept = None
    with open(csv_path, newline="") as f:
        reader = csv_module.DictReader(f)
        for row in reader:
            cpg_id = row.get("CpGmarker", row.get("cpg", "")).strip()
            coef_str = row.get("CoefficientTraining", row.get("coefficient", ""))
            try:
                coef = float(coef_str)
            except (ValueError, TypeError):
                continue
            if cpg_id.lower() in ("(intercept)", "intercept"):
                intercept = coef
            elif cpg_id:
                coefficients[cpg_id] = coef

    if intercept is not None:
        coefficients["intercept"] = intercept
    return coefficients


if __name__ == "__main__":
    print("dnamage_clock.py -- transform functions ready, coefficient table not yet sourced.")
    print("See module docstring for how to get the real 353-CpG coefficients via R.")
    print(f"\nSanity check: anti_trafo(trafo(45)) should recover 45.0")
    print(f"  Result: {anti_trafo(trafo(45.0)):.6f}")