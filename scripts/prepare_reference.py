"""
Step 2: download and index the hg18 reference (UCSC), matching the reference build
Pedersen et al. 2014 and both source papers (Rasmussen 2010/2011) used.

UCSC serves hg18 as one gzipped FASTA per chromosome. This concatenates the
requested set into a single multi-FASTA, then indexes it for both pysam
(samtools faidx, used throughout src/) and bwa (used by align_reads.py).
"""

import gzip
import shutil
import subprocess
from pathlib import Path
from typing import List
from urllib.request import urlretrieve

UCSC_HG18_BASE = "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/chromosomes"

STANDARD_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def download_chromosome(chrom: str, destdir: Path) -> Path:
    destdir.mkdir(parents=True, exist_ok=True)
    gz_path = destdir / f"{chrom}.fa.gz"
    if not gz_path.exists():
        url = f"{UCSC_HG18_BASE}/{chrom}.fa.gz"
        print(f"Downloading {url} ...")
        urlretrieve(url, gz_path)
    return gz_path


def build_reference(destdir: str = "./data/reference", chroms: List[str] = None,
                     output_name: str = "hg18.fa") -> Path:
    """
    Download the requested chromosomes (default: standard chr1-22,X,Y,M -- skip
    _random/_hap contigs unless explicitly requested) and concatenate into one
    FASTA, then index with samtools faidx and bwa index.
    """
    destdir = Path(destdir)
    destdir.mkdir(parents=True, exist_ok=True)
    chroms = chroms or STANDARD_CHROMS

    output_path = destdir / output_name
    with open(output_path, "wb") as out_f:
        for chrom in chroms:
            gz_path = download_chromosome(chrom, destdir)
            with gzip.open(gz_path, "rb") as in_f:
                shutil.copyfileobj(in_f, out_f)

    print(f"Concatenated {len(chroms)} chromosome(s) into {output_path}")

    subprocess.run(["samtools", "faidx", str(output_path)], check=True)
    print("samtools faidx index built.")

    subprocess.run(["bwa", "index", str(output_path)], check=True)
    print("bwa index built.")

    return output_path


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Download and index hg18 reference.")
    p.add_argument("--destdir", default="./data/reference")
    p.add_argument("--chroms", nargs="+", default=None,
                    help="Specific chromosomes only, e.g. --chroms chr12 chr1 "
                         "(default: standard chr1-22,X,Y,M)")
    args = p.parse_args()

    ref_path = build_reference(args.destdir, args.chroms)
    print(f"Reference ready at: {ref_path}")