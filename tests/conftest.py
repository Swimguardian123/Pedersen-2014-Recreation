"""
Shared pytest fixtures: builds small synthetic BAM/FASTA files for tests that
need real pysam objects, not just numpy arrays. Requires pysam (already in
requirements.txt) -- these fixtures were NOT run in the environment that wrote
this file (no pysam available there); run pytest yourself to confirm they work.
"""

import random
import pytest
import pysam


@pytest.fixture
def synthetic_ref_and_bam(tmp_path):
    """
    Builds a 300bp synthetic reference with two deliberately-placed CpG sites
    (at position 50 and position 100), and a BAM with:
      - one forward-strand read starting exactly at the position-50 CpG, with a
        simulated C->T deamination mismatch at its first base
      - one reverse-strand read ending exactly at the position-100 CpG, with a
        simulated G->A deamination mismatch at its last (BAM-stored) base
      - one clean (undamaged) forward-strand read elsewhere, as a negative control

    Returns (fasta_path, bam_path).
    """
    random.seed(42)
    bases = "ACGT"
    ref_len = 300
    ref = list("".join(random.choice(bases) for _ in range(ref_len)))

    # force two known CpG sites, avoiding accidental extra CpGs nearby
    ref[50], ref[51] = "C", "G"
    ref[100], ref[101] = "C", "G"
    ref_seq = "".join(ref)

    fasta_path = tmp_path / "ref.fa"
    with open(fasta_path, "w") as f:
        f.write(">chrTest\n")
        for i in range(0, ref_len, 60):
            f.write(ref_seq[i:i + 60] + "\n")
    pysam.faidx(str(fasta_path))

    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": ref_len, "SN": "chrTest"}]}
    unsorted_bam = tmp_path / "unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_bam), "wb", header=header) as outf:
        read_len = 30

        # 1. Forward-strand read at the position-50 CpG, C->T mismatch at position 0
        seq1 = list(ref_seq[50:50 + read_len])
        seq1[0] = "T"  # deamination: ref C -> observed T
        a1 = pysam.AlignedSegment()
        a1.query_name = "fwd_damaged"
        a1.query_sequence = "".join(seq1)
        a1.flag = 0
        a1.reference_id = 0
        a1.reference_start = 50
        a1.mapping_quality = 40
        a1.cigartuples = [(0, read_len)]
        a1.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        outf.write(a1)

        # 2. Reverse-strand read ending at the position-100 CpG (reference_end-1=101),
        # G->A mismatch at the LAST stored base (BAM-forward-orientation convention)
        rev_start = 101 - read_len + 1  # so reference_end - 1 == 101
        seq2 = list(ref_seq[rev_start:rev_start + read_len])
        seq2[-1] = "A"  # deamination on the reverse strand: ref G -> observed A
        a2 = pysam.AlignedSegment()
        a2.query_name = "rev_damaged"
        a2.query_sequence = "".join(seq2)
        a2.flag = 16  # reverse strand
        a2.reference_id = 0
        a2.reference_start = rev_start
        a2.mapping_quality = 40
        a2.cigartuples = [(0, read_len)]
        a2.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        outf.write(a2)

        # 3. Clean forward-strand read elsewhere (no CpG at its start), negative control
        a3 = pysam.AlignedSegment()
        a3.query_name = "clean_elsewhere"
        a3.query_sequence = ref_seq[200:200 + read_len]
        a3.flag = 0
        a3.reference_id = 0
        a3.reference_start = 200
        a3.mapping_quality = 40
        a3.cigartuples = [(0, read_len)]
        a3.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        outf.write(a3)

    bam_path = tmp_path / "test.bam"
    pysam.sort("-o", str(bam_path), str(unsorted_bam))
    pysam.index(str(bam_path))

    return str(fasta_path), str(bam_path)


@pytest.fixture
def gc_test_bam_and_fasta(tmp_path):
    """
    Reference with a low-GC first half and a high-GC second half. All reads are
    deliberately placed ONLY in the low-GC half -- this directly targets the
    real bug found in this project: the original (buggy) rate model used
    observed reads as its own denominator, meaning a GC stratum with zero
    observed reads simply couldn't appear in the model at all. The fix samples
    background positions independently, so it CAN represent "this GC stratum
    exists in the genome but has zero real coverage" as an explicit rate of 0,
    rather than silently omitting it.
    """
    random.seed(7)
    read_len = 40

    low_gc_half = "".join(random.choice("AATT") for _ in range(500))   # ~0% GC
    high_gc_half = "".join(random.choice("GGCC") for _ in range(500))  # ~100% GC
    ref_seq = low_gc_half + high_gc_half
    ref_len = len(ref_seq)

    fasta_path = tmp_path / "gc_ref.fa"
    with open(fasta_path, "w") as f:
        f.write(">chrGC\n")
        for i in range(0, ref_len, 60):
            f.write(ref_seq[i:i + 60] + "\n")
    pysam.faidx(str(fasta_path))

    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": ref_len, "SN": "chrGC"}]}
    unsorted_bam = tmp_path / "gc_unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_bam), "wb", header=header) as outf:
        for i in range(20):  # 20 reads, all placed in the low-GC half only
            start = random.randint(0, 500 - read_len)
            a = pysam.AlignedSegment()
            a.query_name = f"lowgc_read{i}"
            a.query_sequence = ref_seq[start:start + read_len]
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 40
            a.cigartuples = [(0, read_len)]
            a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            outf.write(a)

    bam_path = tmp_path / "gc_test.bam"
    pysam.sort("-o", str(bam_path), str(unsorted_bam))
    pysam.index(str(bam_path))

    return str(fasta_path), str(bam_path)


@pytest.fixture
def control_set_bams(tmp_path):
    """
    Source BAM: 20 reads, all 40bp (simulating a modern library).
    Target BAM: 5 reads, all 15bp (simulating a short ancient-DNA library).
    """
    random.seed(11)
    ref_seq = "".join(random.choice("ACGT") for _ in range(1000))

    fasta_path = tmp_path / "cs_ref.fa"
    with open(fasta_path, "w") as f:
        f.write(">chrCS\n")
        for i in range(0, len(ref_seq), 60):
            f.write(ref_seq[i:i + 60] + "\n")
    pysam.faidx(str(fasta_path))

    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": len(ref_seq), "SN": "chrCS"}]}

    def _build_bam(path, n_reads, read_len):
        unsorted = str(path) + ".unsorted.bam"
        with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
            for i in range(n_reads):
                start = random.randint(0, len(ref_seq) - read_len)
                a = pysam.AlignedSegment()
                a.query_name = f"read{i}"
                a.query_sequence = ref_seq[start:start + read_len]
                a.flag = 0
                a.reference_id = 0
                a.reference_start = start
                a.mapping_quality = 40
                a.cigartuples = [(0, read_len)]
                a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
                outf.write(a)
        pysam.sort("-o", str(path), unsorted)
        pysam.index(str(path))

    source_path = tmp_path / "source.bam"
    target_path = tmp_path / "target.bam"
    _build_bam(source_path, n_reads=20, read_len=40)
    _build_bam(target_path, n_reads=5, read_len=15)

    return str(source_path), str(target_path)