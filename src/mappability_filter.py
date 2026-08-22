"""
Mappability-based filtering of nucleosome calls.

Paper reference (Hanghoj et al. 2016, Methods, "Nucleosome Calls"):
"we restricted nucleosome calls to uniquely mappable genomic regions, which
were defined following the procedure described in Jonsson et al. (2014),
where we computed the sequence mappability within 20-kb genomic blocks using
a 41-mers (Derrien et al. 2012) and restricted nucleosome calls to blocks
showing mappability uniqueness >= 0.9."

Data source: UCSC ENCODE mappability tracks, real filenames confirmed via
UCSC's own files.txt listing at
hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/files.txt
(an earlier version of this module guessed a plain .bigWig filename that
404'd -- the real files are .bw.gz, fixed here against the actual listing).
40mer (not 41mer) used by default: closest available k-mer length to the
paper's stated 41-mers -- a 1bp difference, a much smaller approximation than
most other substitutions in this project. Smaller/larger k-mers available via
the `kmer` parameter if bandwidth is a concern (see AVAILABLE_KMERS).

Caveat: this uses genome-wide sequence mappability, not full replication of
Jonsson et al. 2014's own exact block-mappability procedure -- a reasonable,
clearly-flagged approximation, same spirit as the other paper-refinement
substitutions in this project.
"""

from typing import List, Optional, Tuple
import gzip
import shutil
from urllib.request import urlretrieve

HG18_MAPPABILITY_BASE = (
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability"
)

# Real file sizes, confirmed directly from UCSC's own files.txt listing (not
# inferred/guessed) -- offered so you can trade fidelity (closer k-mer to the
# paper's stated 41-mers) for bandwidth if 40mer's 1.1G is more than you want:
#   36mer: 1.3G (5bp off)   40mer: 1.1G (1bp off, default)   50mer: 679M (9bp off)
#   75mer: 266M (34bp off)  100mer: 94M (59bp off)
AVAILABLE_KMERS = [36, 40, 50, 75, 100]


def download_mappability_track(destdir: str = "./data/mappability", kmer: int = 40) -> str:
    """
    Downloads and decompresses the real UCSC hg18 CRG mappability track.
    Real filenames end in .bw.gz (gzip-compressed bigWig), confirmed from
    UCSC's own files.txt -- an earlier version of this function assumed a
    plain (uncompressed) .bigWig filename, which 404'd. Fixed here against
    the real listing, with an explicit decompression step pyBigWig needs.
    """
    import os
    if kmer not in AVAILABLE_KMERS:
        raise ValueError(f"kmer must be one of {AVAILABLE_KMERS}, got {kmer}")

    os.makedirs(destdir, exist_ok=True)
    gz_filename = f"wgEncodeCrgMapabilityAlign{kmer}mer.bw.gz"
    bw_filename = f"wgEncodeCrgMapabilityAlign{kmer}mer.bw"
    gz_path = os.path.join(destdir, gz_filename)
    bw_path = os.path.join(destdir, bw_filename)

    if os.path.exists(bw_path):
        print(f"{bw_path} already exists, skipping download.")
        return bw_path

    if not os.path.exists(gz_path):
        url = f"{HG18_MAPPABILITY_BASE}/{gz_filename}"
        print(f"Downloading {url} ... (real bandwidth cost, see AVAILABLE_KMERS sizes)")

        def _progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            if total_size > 0:
                pct = min(100, downloaded * 100 / total_size)
                mb_done = downloaded / 1e6
                mb_total = total_size / 1e6
                print(f"\r  {pct:5.1f}%  ({mb_done:,.1f} / {mb_total:,.1f} MB)", end="", flush=True)

        urlretrieve(url, gz_path, reporthook=_progress)
        print()  # newline after the progress line finishes

    print(f"Decompressing {gz_path} -> {bw_path} ...")
    with gzip.open(gz_path, "rb") as f_in, open(bw_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    return bw_path


def compute_block_mappability(bigwig_path: str, chrom: str, start: int, end: int,
                               block_size: int = 20000) -> List[Tuple[int, int, Optional[float]]]:
    """
    Mean mappability per block_size-bp block across [start, end).
    Returns [(block_start, block_end, mean_mappability_or_None), ...] -- None
    (not a silent 0.0) when the bigWig track has no data for that block at all
    (e.g. an assembly gap), distinguishing "genuinely low mappability" from
    "no data available here" -- these were previously conflated, which could
    make a data-sparse region look identical to a real low-mappability region.
    """
    import pyBigWig

    bw = pyBigWig.open(bigwig_path)
    blocks = []
    pos = start
    while pos < end:
        block_end = min(pos + block_size, end)
        try:
            mean_val = bw.stats(chrom, pos, block_end, type="mean")[0]
        except RuntimeError:
            mean_val = None
        blocks.append((pos, block_end, mean_val))  # None stays None, not coerced to 0.0
        pos = block_end
    bw.close()
    return blocks


def get_high_mappability_regions(bigwig_path: str, chrom: str, start: int, end: int,
                                  block_size: int = 20000,
                                  threshold: float = 0.9) -> List[Tuple[int, int]]:
    """Blocks whose mean mappability >= threshold (paper: 0.9). Blocks with no
    data (None) are excluded, same as failing the threshold, but for a
    different, distinguishable reason -- see compute_block_mappability()."""
    blocks = compute_block_mappability(bigwig_path, chrom, start, end, block_size)
    return [(s, e) for s, e, m in blocks if m is not None and m >= threshold]


def filter_calls_by_mappability(calls: list, genome_offset: int,
                                 high_mappability_regions: List[Tuple[int, int]]) -> list:
    """
    Filter a list of NucleosomeCall objects (from nucleosome_calling.py) to only
    those whose absolute genomic position falls in a high-mappability region.
    Kept generic (duck-typed on .center_pos) rather than importing
    NucleosomeCall directly, to avoid a hard dependency loop between modules.
    """
    kept = []
    for call in calls:
        abs_pos = call.center_pos + genome_offset
        if any(s <= abs_pos < e for s, e in high_mappability_regions):
            kept.append(call)
    return kept


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Compute mappability-based nucleosome-call filtering.")
    p.add_argument("chrom")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--destdir", default="./data/mappability")
    p.add_argument("--kmer", type=int, default=40, choices=AVAILABLE_KMERS,
                    help="k-mer size: 40 (default, closest to paper's 41-mer, 1.1G) "
                         "down to 100 (94M, less precise) for smaller downloads")
    p.add_argument("--block-size", type=int, default=20000)
    p.add_argument("--threshold", type=float, default=0.9)
    args = p.parse_args()

    track_path = download_mappability_track(args.destdir, kmer=args.kmer)
    regions = get_high_mappability_regions(track_path, args.chrom, args.start, args.end,
                                            args.block_size, args.threshold)
    total_blocks = (args.end - args.start) // args.block_size + 1
    print(f"{len(regions)}/{total_blocks} block(s) pass mappability >= {args.threshold}")