"""
Step 3: align FASTQ reads to the hg18 reference, producing a sorted/indexed BAM.

Paper reference: for the ancient horse and polar bear comparison samples, Pedersen
et al. 2014's Methods ("Sequence data sets") state reads were "mapped using BWA
with standard parameters and disabling the seed". Disabling the seed is standard
ancient-DNA practice (also used by mapDamage/PALEOMIX-style pipelines): short,
damaged aDNA fragments don't benefit from BWA's default seed-based heuristic and
can lose sensitivity to true alignments if it's left on. We use `bwa aln` (not
`bwa mem`) for this reason -- `bwa aln` is the algorithm that supports seed
disabling via a very large `-l` value, and is the more historically-appropriate
choice for short single-end aDNA reads like these (Saqqaq/Aboriginal libraries
were sequenced single-end, per the SI's "Sequencing" sections).

This is a real fidelity note, not a hidden detail: modern pipelines increasingly
default to `bwa mem`, which behaves differently on short/damaged reads. If you'd
rather match current best-practice aDNA pipelines instead of the paper's literal
choice, that's a reasonable adjustment -- flagged here so it's a deliberate choice,
not an accidental substitution.
"""

import subprocess
from pathlib import Path


def align_single_end(fastq_path: str, reference_fasta: str, output_bam: str,
                      disable_seed: bool = True, seed_length: int = 1024,
                      max_edit_distance: float = 0.04, threads: int = 4) -> Path:
    """
    bwa aln -> bwa samse -> samtools sort -> samtools index, matching the paper's
    stated "BWA, standard parameters, disabling the seed" approach.

    disable_seed: sets -l to seed_length (bwa treats -l >= read length as
        effectively disabling seeding, since no read will exceed it).
    max_edit_distance: bwa aln's -n parameter (fraction, default 0.04 = bwa's own
        default -- not stated explicitly in the paper's main text for Saqqaq/
        Aboriginal alignment, only "standard parameters" -- left at bwa's default
        rather than guessed at a different value).
    """
    fastq_path = Path(fastq_path)
    reference_fasta = Path(reference_fasta)
    output_bam = Path(output_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    sai_path = output_bam.with_suffix(".sai")
    sam_path = output_bam.with_suffix(".sam")

    aln_cmd = ["bwa", "aln", "-t", str(threads), "-n", str(max_edit_distance)]
    if disable_seed:
        aln_cmd += ["-l", str(seed_length)]
    aln_cmd += [str(reference_fasta), str(fastq_path)]

    print(f"Running: {' '.join(aln_cmd)}")
    with open(sai_path, "wb") as sai_out:
        subprocess.run(aln_cmd, stdout=sai_out, check=True)

    samse_cmd = ["bwa", "samse", str(reference_fasta), str(sai_path), str(fastq_path)]
    print(f"Running: {' '.join(samse_cmd)}")
    with open(sam_path, "wb") as sam_out:
        subprocess.run(samse_cmd, stdout=sam_out, check=True)

    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(output_bam), str(sam_path)]
    print(f"Running: {' '.join(sort_cmd)}")
    subprocess.run(sort_cmd, check=True)

    subprocess.run(["samtools", "index", str(output_bam)], check=True)

    sai_path.unlink(missing_ok=True)
    sam_path.unlink(missing_ok=True)

    return output_bam


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Align single-end FASTQ reads to a reference (aDNA-appropriate BWA settings).")
    p.add_argument("fastq")
    p.add_argument("reference_fasta")
    p.add_argument("output_bam")
    p.add_argument("--no-disable-seed", action="store_true",
                    help="Use BWA's normal seeding behavior instead of the aDNA-style disabled seed")
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()

    bam = align_single_end(args.fastq, args.reference_fasta, args.output_bam,
                            disable_seed=not args.no_disable_seed, threads=args.threads)
    print(f"BAM ready at: {bam}")