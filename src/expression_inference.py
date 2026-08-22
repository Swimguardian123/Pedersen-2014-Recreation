"""
Gene expression inference from nucleosome/methylation proxies.

Paper reference (Methods, "Expression analysis"):
"Three proxies for expression were defined: (1) the level of gene body to promoter
methylation, Rs; (2) occupancy of the +1 nucleosome; and (3) strength of nucleosome
phasing. Their respective performance was evaluated using expression data from 10
samples of modern hair follicles (Kim et al. 2006) (GSE3058)... The modern expression
data was used to define groups of genes with increasing expression levels (using 10,
20, or 50 quantiles), which were ranked by each of the three proxies, and the
Spearman correlation coefficient evaluated."

Results section, defining the three proxies precisely:
- Rs: "we further calculated the ratio Rs of gene body to promoter methylation as a
  proxy for gene expression, with low (high) Rs values indicating low (high)
  expression levels" (following Ball et al. 2009). Rs = Ms(gene body) / Ms(promoter)
  -- active genes have LOW promoter methylation and HIGHER gene-body methylation, so
  this ratio rises with expression.
- First-nucleosome occupancy: "the presence of a strongly positioned +1 nucleosome" ->
  average GC-corrected read depth over the +1 nucleosome region downstream of the TSS.
- Phasing strength: "the level of downstream regularly spaced phasing of nucleosomes"
  via Fourier transform analysis across the TSS region -- reuses periodicity.py's
  Welch periodogram, taking power at the ~193bp nucleosome-periodicity frequency as
  the phasing-strength score.

Data: GSE3058 (Kim et al. 2006), 10 modern hair-follicle expression microarray
samples, public on GEO. Fetched here via GEOparse (`pip install GEOparse`) -- this
needs your own internet access to run, not something done on my end.
"""

from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.stats import spearmanr

try:
    from .periodicity import welch_periodogram
except ImportError:
    from periodicity import welch_periodogram


# ---------------------------------------------------------------------------
# Proxy 1: Rs, gene body / promoter methylation ratio
# ---------------------------------------------------------------------------

def compute_rs(genebody_ms: float, promoter_ms: float, pseudocount: float = 1e-6) -> float:
    """
    Rs = Ms(gene body) / Ms(promoter). Higher Rs -> higher inferred expression.
    pseudocount avoids divide-by-zero when promoter Ms happens to be exactly 0
    (fully unmethylated promoter, common at active CGI promoters).
    """
    return genebody_ms / (promoter_ms + pseudocount)


# ---------------------------------------------------------------------------
# Proxy 2: +1 nucleosome occupancy
# ---------------------------------------------------------------------------

def first_nucleosome_occupancy(gc_corrected_depth: np.ndarray, tss_offset: int = 0,
                                window: int = 147) -> float:
    """
    Mean GC-corrected depth over the +1 nucleosome window downstream of the TSS.

    gc_corrected_depth: depth array covering (at least) [TSS, TSS + tss_offset +
        window). Index 0 of the array should correspond to the TSS itself.
    tss_offset: bp from TSS to the start of the +1 nucleosome window. The paper
        doesn't give an exact number for this main-text proxy (Fig. 2A shows +1
        occupancy peaking shortly downstream of the TSS) -- 0 is a reasonable
        default (start the window right at the TSS), adjust if your own peak-calling
        on nucleosome_calling.py output finds the +1 dyad sitting further downstream.
    window: nucleosome width, 147bp per the rest of this pipeline.
    """
    start = tss_offset
    end = tss_offset + window
    segment = gc_corrected_depth[start:end]
    return float(np.nanmean(segment)) if len(segment) else float("nan")


# ---------------------------------------------------------------------------
# Proxy 3: nucleosome phasing strength
# ---------------------------------------------------------------------------

def phasing_strength(depth_around_tss: np.ndarray, fs: float = 1.0,
                      target_period: float = 193.0, tolerance: float = 15.0) -> float:
    """
    Power of the Welch periodogram at the nucleosome-periodicity frequency band
    (paper's peak TSS periodicity: 193bp, Results/Fig. 2B), reused directly from
    periodicity.py rather than reimplementing the FFT machinery.
    """
    freqs, power = welch_periodogram(depth_around_tss, fs=fs)
    with np.errstate(divide="ignore"):
        periods = np.where(freqs > 0, 1.0 / freqs, np.inf)
    mask = np.abs(periods - target_period) <= tolerance
    return float(power[mask].max()) if mask.any() else float("nan")


# ---------------------------------------------------------------------------
# GSE3058 expression data + validation against the three proxies
# ---------------------------------------------------------------------------

def fetch_gse3058_expression(destdir: str = "./geo_cache") -> Dict[str, float]:
    """
    Download and parse GSE3058 (Kim et al. 2006, 10 modern hair-follicle samples)
    via GEOparse, returning a per-gene mean expression value averaged across all
    10 samples (matching the paper's use of this as a single expression ranking,
    not a per-sample analysis).

    Requires internet access and `pip install GEOparse` -- run this yourself.
    """
    import GEOparse  # local import: optional dependency, only needed for this fetch

    gse = GEOparse.get_GEO(geo="GSE3058", destdir=destdir)

    # Each GSM (sample) has a table with probe IDs and expression values; average
    # across samples per probe, then let the caller map probe -> gene via whatever
    # platform annotation (GPL) file matches this series.
    all_values: Dict[str, List[float]] = {}
    for gsm_name, gsm in gse.gsms.items():
        table = gsm.table
        # typical GEO sample table columns: ID_REF, VALUE (column names can vary by series)
        id_col = "ID_REF" if "ID_REF" in table.columns else table.columns[0]
        val_col = "VALUE" if "VALUE" in table.columns else table.columns[1]
        for probe_id, value in zip(table[id_col], table[val_col]):
            all_values.setdefault(str(probe_id), []).append(float(value))

    return {probe_id: float(np.mean(vals)) for probe_id, vals in all_values.items()}


def rank_into_quantiles(expression: Dict[str, float], n_quantiles: int = 10
                         ) -> Dict[str, int]:
    """
    Bin genes into n_quantiles groups of increasing expression (paper tests 10, 20,
    and 50 quantile granularities). Returns {gene_id: quantile_index} (0 = lowest
    expression group, n_quantiles-1 = highest).
    """
    genes = list(expression.keys())
    values = np.array([expression[g] for g in genes])
    # rankdata-style quantile assignment: percentile rank -> bin index
    ranks = np.argsort(np.argsort(values))  # 0-based rank, ties broken by order
    bin_idx = (ranks * n_quantiles) // len(values)
    bin_idx = np.clip(bin_idx, 0, n_quantiles - 1)
    return dict(zip(genes, bin_idx.tolist()))


def evaluate_proxy(proxy_values: Dict[str, float], expression_groups: Dict[str, int]
                    ) -> Tuple[float, float]:
    """
    Spearman correlation between a proxy (Rs, +1 occupancy, or phasing strength) and
    expression-quantile group, over genes present in both. Returns (rho, p_value).
    """
    common_genes = sorted(set(proxy_values) & set(expression_groups))
    if len(common_genes) < 3:
        return float("nan"), float("nan")
    x = [proxy_values[g] for g in common_genes]
    y = [expression_groups[g] for g in common_genes]
    rho, p = spearmanr(x, y)
    return float(rho), float(p)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Fetch GSE3058 and report per-quantile-granularity setup (no proxy "
                     "values here -- those come from your own Rs/occupancy/phasing runs)."
    )
    p.add_argument("--destdir", default="./geo_cache")
    p.add_argument("--n-quantiles", type=int, default=10)
    args = p.parse_args()

    expr = fetch_gse3058_expression(args.destdir)
    print(f"Fetched expression for {len(expr)} probes from GSE3058.")
    groups = rank_into_quantiles(expr, args.n_quantiles)
    print(f"Assigned to {args.n_quantiles} quantile groups.")