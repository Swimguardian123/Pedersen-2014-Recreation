"""
Functional enrichment of expression-correlated genes.

Paper reference (Results, expression section): genes ranked by expression proxy
(Rs, +1 occupancy, or phasing strength) were checked for functional enrichment
among top/bottom-ranked genes using DAVID.

DAVID itself isn't a library -- it's a specific web tool. Per the earlier
assessment (first message of this project), the honest substitute is Enrichr
(same general purpose: gene-list -> GO/pathway enrichment), via the `gseapy`
package, which wraps Enrichr's API. This is a different backend producing
comparable output, not a reimplementation of DAVID itself -- flagged plainly,
not silently swapped.
"""

from typing import Dict, List

import numpy as np


def top_and_bottom_genes(proxy_values: Dict[str, float], n: int = 1000
                          ) -> Dict[str, List[str]]:
    """
    Split genes into top-n and bottom-n by a given expression proxy (Rs, +1
    occupancy, or phasing strength from expression_inference.py), matching the
    paper's "top/bottom ranked genes" framing for enrichment testing.
    """
    ranked = sorted(proxy_values.items(), key=lambda kv: kv[1], reverse=True)
    genes = [g for g, _ in ranked]
    return {
        "top": genes[:n],
        "bottom": genes[-n:] if n <= len(genes) else genes,
    }


def run_enrichr(gene_list: List[str], gene_sets: List[str] = None,
                 description: str = "saqqaq_repro") -> "pandas.DataFrame":
    """
    Query Enrichr via gseapy for functional enrichment of a gene list -- the
    practical substitute for DAVID. Requires internet access (run on your end).

    gene_sets: Enrichr library names, e.g. ["GO_Biological_Process_2021",
        "KEGG_2021_Human", "GO_Molecular_Function_2021"]. Defaults to a small
        reasonable set if not given.
    """
    import gseapy as gp  # local import: optional dependency, only needed here

    if gene_sets is None:
        gene_sets = ["GO_Biological_Process_2021", "KEGG_2021_Human"]

    enr = gp.enrichr(gene_list=gene_list, gene_sets=gene_sets,
                      description=description, outdir=None)
    return enr.results


def enrichment_for_proxy(proxy_values: Dict[str, float], n: int = 1000,
                          gene_sets: List[str] = None) -> Dict[str, "pandas.DataFrame"]:
    """Run Enrichr separately on the top-n and bottom-n gene sets for one proxy."""
    split = top_and_bottom_genes(proxy_values, n=n)
    return {
        "top": run_enrichr(split["top"], gene_sets=gene_sets, description="top_expression"),
        "bottom": run_enrichr(split["bottom"], gene_sets=gene_sets, description="bottom_expression"),
    }


if __name__ == "__main__":
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="Run Enrichr (DAVID substitute) on top/bottom genes by an expression proxy."
    )
    p.add_argument("proxy_json", help="JSON file: {gene_id: proxy_value, ...}")
    p.add_argument("--n", type=int, default=1000)
    p.add_argument("--gene-sets", nargs="+", default=None)
    args = p.parse_args()

    with open(args.proxy_json) as f:
        proxy_values = json.load(f)

    results = enrichment_for_proxy(proxy_values, n=args.n, gene_sets=args.gene_sets)
    for direction, df in results.items():
        print(f"\n=== {direction} genes ===")
        print(df.head(10) if df is not None and len(df) else "No results / empty.")