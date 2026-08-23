"""
hg19 -> hg18 liftover for the 353 Horvath clock CpG coordinates, using
pyliftover (pure Python, pip-installable) instead of the compiled UCSC
liftOver binary -- avoids the x86_64/arm64 architecture mismatch hit with
the native binary on Apple Silicon. pyliftover is independently verified by
its own maintainers to produce identical output to the real UCSC tool for
single-point coordinates (not full ranges) -- an exact fit here, since each
CpG is a single genomic position, not a region.

Reuses the SAME chain file already downloaded for the native-binary attempt
(hg19ToHg18.over.chain.gz) -- no new download needed.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def liftover_bed(in_bed_path: str, chain_path: str, out_bed_path: str) -> None:
    from pyliftover import LiftOver

    lo = LiftOver(chain_path)

    n_total = 0
    n_lifted = 0
    n_unmapped = 0

    with open(in_bed_path) as fin, open(out_bed_path, "w") as fout:
        for line in fin:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            chrom, start, end, cpg_id = fields[0], int(fields[1]), int(fields[2]), fields[3]
            strand = fields[5] if len(fields) > 5 else "+"
            n_total += 1

            result = lo.convert_coordinate(chrom, start, strand)
            if not result:
                n_unmapped += 1
                continue

            new_chrom, new_pos, new_strand, _ = result[0]
            fout.write(f"{new_chrom}\t{new_pos}\t{new_pos + 1}\t{cpg_id}\t0\t{new_strand}\n")
            n_lifted += 1

    print(f"Lifted {n_lifted}/{n_total} CpGs from hg19 to hg18 ({n_unmapped} unmapped)")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Liftover the 353 clock CpG BED, hg19 -> hg18.")
    p.add_argument("in_bed", help="hg19 BED (from prepare_dnamage_coordinates.py)")
    p.add_argument("chain_file", help="hg19ToHg18.over.chain.gz")
    p.add_argument("--out", default="horvath_cpgs_hg18.bed")
    args = p.parse_args()

    liftover_bed(args.in_bed, args.chain_file, args.out)