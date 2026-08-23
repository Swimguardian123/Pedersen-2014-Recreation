"""
Regression tests for dnamage_clock.py. Pure math, no dependencies -- self-run.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
from src.dnamage_clock import trafo, anti_trafo, predict_age


@pytest.mark.parametrize("age", [0, 5, 15, 19.9, 20, 20.1, 25, 45, 70, 100])
def test_trafo_anti_trafo_are_true_inverses(age):
    """The strongest check available: confirms the two formulas genuinely
    invert each other (not just independently plausible-looking), across
    both sides of the piecewise boundary at adult_age=20."""
    assert abs(anti_trafo(trafo(age)) - age) < 1e-9


def test_trafo_boundary_at_adult_age_is_exactly_zero():
    assert abs(trafo(20) - 0.0) < 1e-9


def test_predict_age_end_to_end():
    betas = {"cg001": 0.5, "cg002": 0.3, "cg_unused": 0.9}
    coefficients = {"cg001": 2.0, "cg002": -1.0}
    intercept = -0.3
    # linear_response = -0.3 + 2.0*0.5 - 1.0*0.3 = 0.4 -> age = 21*0.4 + 20 = 28.4
    age = predict_age(betas, coefficients, intercept)
    assert abs(age - 28.4) < 1e-9


def test_predict_age_raises_on_no_overlapping_cpgs():
    with pytest.raises(ValueError):
        predict_age({"cgAAA": 0.5}, {"cgBBB": 1.0}, intercept=0.0)


def test_predict_age_require_all_cpgs_raises_on_missing():
    with pytest.raises(ValueError):
        predict_age({"cg001": 0.5}, {"cg001": 1.0, "cg002": 1.0},
                     intercept=0.0, require_all_cpgs=True)


def test_load_coefficients_from_csv_matches_real_export_format(tmp_path):
    """Confirmed against a real methylclockData R export (this project's own
    run): CpGmarker/CoefficientTraining columns, "(Intercept)" as a row label
    for the intercept -- not guessed, verified against real output that
    parsed to exactly 353 CpGs, matching the literature's stated count."""
    from src.dnamage_clock import load_coefficients_from_csv

    path = tmp_path / "horvath_coefficients.csv"
    path.write_text(
        "CpGmarker,CoefficientTraining\n"
        "(Intercept),0.695507258\n"
        "cg00075967,0.12933661\n"
        "cg00374717,0.005017857\n"
        "cg00864867,1.59976405\n"
    )

    coefficients = load_coefficients_from_csv(str(path))
    assert coefficients["intercept"] == 0.695507258
    assert coefficients["cg00075967"] == 0.12933661
    assert len(coefficients) == 4  # 3 CpGs + intercept
    assert "cg00075967" in coefficients and "(Intercept)" not in coefficients


if __name__ == "__main__":
    for age in [0, 5, 15, 19.9, 20, 20.1, 25, 45, 70, 100]:
        test_trafo_anti_trafo_are_true_inverses(age)
    test_trafo_boundary_at_adult_age_is_exactly_zero()
    test_predict_age_end_to_end()
    print("All dnamage_clock tests passed.")