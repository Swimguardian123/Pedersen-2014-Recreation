"""
Resolve an old-style SRA study accession (e.g. SRA010102, SRA035301.1) to its
child run accessions (SRRxxxxxx), and download reads via sra-tools.

These two accessions predate the modern SRR/SRX naming convention, so
`prefetch`/`fasterq-dump` can't be pointed at them directly -- they need to be
resolved to run-level accessions first via NCBI's Entrez eutils.

Requires: pysradb (pip install pysradb) OR falls back to raw eutils XML if
pysradb can't resolve the old-format accession; sra-tools (prefetch,
fasterq-dump, fastq-dump) must be installed and on PATH separately (not
pip-installable).

Output naming: --split-files always appends _1 (and _2 for paired data) to the
run accession, even for genuinely single-end libraries like these -- so the
output is e.g. SRR030833_1.fastq, never a bare SRR030833.fastq. Point
align_reads.py at that _1 file.
"""

import subprocess
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List
from urllib.request import urlopen
from urllib.parse import urlencode

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _require_tool(name: str, install_hint: str) -> None:
    """Fail fast with an actionable message instead of an opaque FileNotFoundError
    traceback if a required external binary isn't on PATH."""
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Required tool '{name}' not found on PATH. Install it first: {install_hint}"
        )


def resolve_via_pysradb(study_accession: str) -> List[str]:
    from pysradb.sraweb import SRAweb
    db = SRAweb()
    df = db.sra_metadata(study_accession)
    if df is None or df.empty or "run_accession" not in df.columns:
        return []
    return df["run_accession"].dropna().unique().tolist()


def resolve_via_eutils(study_accession: str) -> List[str]:
    """Fallback: raw esearch (sra db) -> efetch (docsum XML) -> extract SRR runs."""
    search_url = f"{EUTILS_BASE}/esearch.fcgi?" + urlencode(
        {"db": "sra", "term": study_accession, "retmax": 500}
    )
    with urlopen(search_url) as resp:
        search_xml = ET.fromstring(resp.read())
    ids = [el.text for el in search_xml.findall(".//Id")]
    if not ids:
        return []

    fetch_url = f"{EUTILS_BASE}/efetch.fcgi?" + urlencode(
        {"db": "sra", "id": ",".join(ids), "rettype": "runinfo", "retmode": "text"}
    )
    with urlopen(fetch_url) as resp:
        text = resp.read().decode("utf-8", errors="ignore")

    # runinfo CSV: first column is Run (SRRxxxxxx)
    runs = []
    for line in text.splitlines()[1:]:
        if line.strip():
            runs.append(line.split(",")[0])
    return [r for r in runs if r.startswith("SRR") or r.startswith("ERR") or r.startswith("DRR")]


def resolve_run_accessions(study_accession: str) -> List[str]:
    try:
        runs = resolve_via_pysradb(study_accession)
        if runs:
            return runs
    except ImportError:
        print("pysradb not installed, falling back to raw eutils", file=sys.stderr)
    except Exception as e:
        print(f"pysradb resolution failed ({e}), falling back to raw eutils", file=sys.stderr)

    return resolve_via_eutils(study_accession)


def download_run(run_accession: str, destdir: Path, max_reads: int = None) -> None:
    _require_tool("prefetch", "brew install sra-tools  (macOS)  /  "
                               "conda install -c bioconda sra-tools  (conda)")
    destdir.mkdir(parents=True, exist_ok=True)

    subprocess.run(["prefetch", run_accession, "--output-directory", str(destdir)], check=True)

    # prefetch's own convention: <destdir>/<run_accession>/<run_accession>.sra --
    # fasterq-dump needs this explicit path when the file isn't in its default
    # cache location (~/ncbi/public/sra/), otherwise it tries to re-resolve the
    # accession over the network instead of using the file already on disk.
    sra_path = destdir / run_accession / f"{run_accession}.sra"
    if not sra_path.exists():
        candidates = list(destdir.rglob(f"{run_accession}.sra"))
        if not candidates:
            raise RuntimeError(
                f"prefetch reported success but no {run_accession}.sra file was found "
                f"under {destdir} -- check prefetch's own output above for its actual save location."
            )
        sra_path = candidates[0]

    if max_reads:
        # fasterq-dump (3.x) dropped the older fastq-dump's -X/--maxSpotId spot-range
        # limiting in favor of full-extraction speed -- use classic fastq-dump instead
        # when a read cap is requested; it's slower but is the tool that actually
        # supports capping, which fasterq-dump's newer design does not.
        _require_tool("fastq-dump", "same sra-tools package as prefetch/fasterq-dump -- "
                                     "should already be installed alongside them")
        cmd = ["fastq-dump", str(sra_path), "--outdir", str(destdir), "--split-files",
               "-X", str(max_reads)]
    else:
        cmd = ["fasterq-dump", str(sra_path), "--outdir", str(destdir), "--split-files"]

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Resolve an SRA study accession to run accessions and download reads."
    )
    p.add_argument("study_accession", help="e.g. SRA010102 (Saqqaq) or SRA035301.1 (Aboriginal)")
    p.add_argument("--destdir", default="./data/reads")
    p.add_argument("--max-reads", type=int, default=None,
                    help="Limit reads per run for a quick pipeline-validation download "
                         "(omit for a full download)")
    p.add_argument("--max-runs", type=int, default=None,
                    help="Only download the first N resolved run accessions -- use this "
                         "for a sanity check against a study with many sequencing lanes, "
                         "since --max-reads alone still downloads every run's full .sra file")
    p.add_argument("--list-only", action="store_true",
                    help="Just print resolved run accessions, don't download")
    args = p.parse_args()

    runs = resolve_run_accessions(args.study_accession)
    if not runs:
        print(f"No run accessions resolved for {args.study_accession}. "
              "Try searching https://www.ncbi.nlm.nih.gov/sra directly.")
        sys.exit(1)

    print(f"Resolved {len(runs)} run accession(s): {runs}")
    if args.list_only:
        sys.exit(0)

    if args.max_runs:
        runs = runs[:args.max_runs]
        print(f"Limiting to first {len(runs)} run(s) per --max-runs.")

    for run in runs:
        print(f"Downloading {run} ...")
        download_run(run, Path(args.destdir), max_reads=args.max_reads)