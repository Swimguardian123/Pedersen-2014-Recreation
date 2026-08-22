"""
Regression tests for nucleotide_patterns.py. Requires pysam -- not run in the
environment that wrote this file.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pysam
from src.nucleosome_calling import NucleosomeCall
from src.nucleotide_patterns import (
    mononucleotide_matrix, select_top_calls, purine_pyrimidine_dinucleotide_signal, BASES,
)


def _write_fasta(path, seq, name="chrNP"):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(path))


def test_mononucleotide_matrix_exact_frequency_with_single_call(tmp_path):
    """With exactly one call, the frequency at every position is trivially 1.0
    for whichever base is actually there -- a direct, exact check that the
    windowing/indexing logic pulls the right bases from the right positions."""
    ref_seq = "A" * 100 + "G" + "C" * 100  # known base ('G') at position 100
    fasta_path = tmp_path / "np_ref.fa"
    _write_fasta(fasta_path, ref_seq)

    call = NucleosomeCall(center_pos=100, peak_depth=5.0, score=5.0)
    matrix = mononucleotide_matrix(str(fasta_path), "chrNP", [call], genome_offset=0, halfwidth=5)

    # center row (index halfwidth=5) should be 100% G, since only one call and
    # the reference base at the exact center is 'G'
    center_row = matrix[5]
    assert center_row[BASES.index("G")] == 1.0
    assert center_row.sum() == 1.0  # all other bases 0

    # a position clearly in the all-A region should show 100% A
    left_row = matrix[0]  # offset -5 from center = position 95, still in the "A" run
    assert left_row[BASES.index("A")] == 1.0


def test_select_top_calls_stratifies_by_score_correctly():
    calls = [
        NucleosomeCall(center_pos=i, peak_depth=float(i), score=float(i))
        for i in [10, 30, 50, 70, 90]
    ]
    top = select_top_calls(calls, top_fraction=0.4)  # top 40% of 5 = 2 calls
    assert len(top) == 2
    assert {c.score for c in top} == {90.0, 70.0}  # the two highest scores


def test_select_top_calls_none_means_all():
    calls = [NucleosomeCall(center_pos=i, peak_depth=1.0, score=1.0) for i in range(5)]
    assert select_top_calls(calls, top_fraction=None) == calls


def test_purine_purine_signal_exact_on_known_sequence(tmp_path):
    """A reference built as alternating purine-purine / pyrimidine-pyrimidine
    blocks lets us predict the exact purine-purine dinucleotide signal at each
    offset."""
    # positions 95-104 (10bp): all purines (A/G) -> every adjacent pair is RR
    purine_block = "AGAGAGAGAG"
    ref_seq = "C" * 95 + purine_block + "C" * 95  # 200bp total, purine block at 95-104
    fasta_path = tmp_path / "np_ref2.fa"
    _write_fasta(fasta_path, ref_seq)

    call = NucleosomeCall(center_pos=100, peak_depth=5.0, score=5.0)  # center inside purine block
    signal = purine_pyrimidine_dinucleotide_signal(str(fasta_path), "chrNP", [call],
                                                     genome_offset=0, halfwidth=4)
    # positions fully inside the purine block should show RR=1.0 (100%, since
    # only one call, the pair is either purine-purine or not, deterministically)
    center_idx = 4  # corresponds to position 100-101, both purines
    assert signal[center_idx] == 1.0


if __name__ == "__main__":
    print("(This file needs pysam and pytest fixtures -- run via `pytest tests/` "
          "rather than directly.)")