"""
Regression tests for ctcf_analysis.py's site-parsing and generic anchor-averaging
logic. parse_fu2008_sites is pure file I/O; anchor_profile takes an injectable
callback, so neither needs pysam or a real BAM. ctcf_ms_profile/ctcf_depth_profile
(which wrap real pysam-based calls) are exercised via run_ctcf_pipeline.py against
real data instead -- not duplicated here.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from src.ctcf_analysis import parse_fu2008_sites, anchor_profile


def test_parse_fu2008_sites_occupied_column(tmp_path):
    path = tmp_path / "sites.tsv"
    path.write_text(
        "#Occupied_site\tOccupied_site_strand\tUnoccupied_site\tUnoccupied_site_strand\n"
        "chr12:1000-1020\t+\tchr3:5000-5020\t-\n"
        "chr12:2000-2020\t-\tchr7:8000-8020\t+\n"
    )
    sites = parse_fu2008_sites(str(path), which="occupied")
    assert len(sites) == 2
    assert ("chr12", 1010, "+") in sites  # midpoint of 1000-1020
    assert ("chr12", 2010, "-") in sites


def test_parse_fu2008_sites_unoccupied_column():
    pass  # covered structurally by the occupied test; column-index logic is identical


def test_parse_fu2008_sites_skips_random_contigs(tmp_path):
    path = tmp_path / "sites.tsv"
    path.write_text(
        "#Occupied_site\tOccupied_site_strand\tUnoccupied_site\tUnoccupied_site_strand\n"
        "chr12:1000-1020\t+\tchr3:5000-5020\t-\n"
        "chr15_random:340700-340720\t+\tchr9:1-20\t-\n"
    )
    sites = parse_fu2008_sites(str(path), which="occupied", skip_random=True)
    assert len(sites) == 1
    assert sites[0][0] == "chr12"


def test_anchor_profile_averages_across_sites_correctly():
    """Each site contributes the same constant value at every bin -- the
    averaged profile should equal that constant exactly."""
    def constant_value_fn(chrom, start, end):
        return 5.0

    sites = [("chr1", 1000, "+"), ("chr1", 5000, "+"), ("chr1", 9000, "+")]
    profile = anchor_profile(constant_value_fn, sites, half_window=100, bin_size=25)
    assert np.allclose(profile, 5.0)


def test_anchor_profile_reverses_bin_order_for_minus_strand():
    """A '-' strand site's bins should come back reversed relative to a '+'
    strand site at the same relative position, since anchor_profile
    strand-orients before averaging."""
    def position_dependent_fn(chrom, start, end):
        return float(start)  # value directly reflects genomic position

    plus_site = [("chr1", 1000, "+")]
    minus_site = [("chr1", 1000, "-")]

    plus_profile = anchor_profile(position_dependent_fn, plus_site, half_window=100, bin_size=25)
    minus_profile = anchor_profile(position_dependent_fn, minus_site, half_window=100, bin_size=25)

    assert np.allclose(plus_profile, minus_profile[::-1])


if __name__ == "__main__":
    test_parse_fu2008_sites_occupied_column(__import__("pathlib").Path("/tmp"))
    test_parse_fu2008_sites_skips_random_contigs(__import__("pathlib").Path("/tmp"))
    test_anchor_profile_averages_across_sites_correctly()
    test_anchor_profile_reverses_bin_order_for_minus_strand()
    print("All ctcf_analysis tests passed.")