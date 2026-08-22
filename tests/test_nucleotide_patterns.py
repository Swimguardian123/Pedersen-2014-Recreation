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
    """A reference built with the purine block flanked TIGHTLY by pyrimidines on
    BOTH sides (not just one) lets this test catch an off-by-one indexing bug
    in either direction -- an earlier version of this test used a longer purine
    block with margin on only one side, which was concretely demonstrated (by
    directly simulating both shift directions) to only catch a rightward
    off-by-one error, not a leftward one. Fixed here, not just noted."""
    # purine block EXACTLY at positions 96-104 (9bp, matching the 2*halfwidth+1
    # window width below) -- flanked immediately by 'C' on both sides, so ANY
    # 1bp shift in the fetch window (either direction) pulls in a 'C'.
    purine_block = "AGAGAGAGA"
    ref_seq = "C" * 96 + purine_block + "C" * 95  # 200bp total
    fasta_path = tmp_path / "np_ref2.fa"
    _write_fasta(fasta_path, ref_seq)

    call = NucleosomeCall(center_pos=100, peak_depth=5.0, score=5.0)  # center inside purine block
    signal = purine_pyrimidine_dinucleotide_signal(str(fasta_path), "chrNP", [call],
                                                     genome_offset=0, halfwidth=4)
    # the entire fetched window (positions 96-104) is purine-only, so every
    # adjacent pair should show RR=1.0 across the WHOLE signal array, not just
    # the center -- checking all positions (not just index 4) is what actually
    # makes a 1bp shift in either direction detectable
    assert np.allclose(signal, 1.0), (
        f"expected all-1.0 signal (tight purine window), got {signal} -- "
        f"a 1bp off-by-one shift in either direction would show up as a <1.0 "
        f"value somewhere in this array"
    )


if __name__ == "__main__":
    print("(This file needs pysam and pytest fixtures -- run via `pytest tests/` "
          "rather than directly.)")