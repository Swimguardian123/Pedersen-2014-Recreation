"""
Regression tests for functional_enrichment.py's pure gene-list logic
(top_and_bottom_genes). run_enrichr/enrichment_for_proxy need live Enrichr API
access and are deliberately not covered here, same reasoning as other
network-dependent functions in this project.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.functional_enrichment import top_and_bottom_genes


def test_top_and_bottom_split_correctly():
    proxy = {f"gene{i}": float(i) for i in range(10)}  # gene9 highest, gene0 lowest
    result = top_and_bottom_genes(proxy, n=3)

    assert result["top"] == ["gene9", "gene8", "gene7"]
    assert result["bottom"] == ["gene2", "gene1", "gene0"]


def test_n_larger_than_gene_list_returns_everything():
    proxy = {"a": 1.0, "b": 2.0, "c": 3.0}
    result = top_and_bottom_genes(proxy, n=100)
    assert len(result["top"]) == 3
    assert len(result["bottom"]) == 3


if __name__ == "__main__":
    test_top_and_bottom_split_correctly()
    test_n_larger_than_gene_list_returns_everything()
    print("All functional_enrichment tests passed.")