# Pedersen et al. 2014 — Computational Reconstruction

A from-scratch reimplementation of the computational pipeline from Pedersen
et al. 2014, *Genome Research* 24:454-466 ("Genome-wide nucleosome map and
cytosine methylation levels of an ancient human genome") — built by reading
the paper's Methods closely enough to translate its stated algorithms into
working, tested code, and validated against real (small-scale) public data.

**Start here, depending on what you need:**
- New to this repo / explaining it to someone else → **`PROJECT_SUMMARY.md`**
- Actually running something → **`QUICKSTART.md`** (the real, tested command sequence)
- Which statistical methods are genuinely computed vs. consumed from external
  papers → **`BIOSTATISTICS.md`**
- Module-by-module fidelity to the paper → this file, below

This is a mechanically-validated pipeline, not a biological reproduction of
the paper's findings — real conclusions need the paper's actual ~16-20x
genome coverage (~2.5 billion reads), which is well beyond what's practical
here. See `PROJECT_SUMMARY.md` for the full reasoning.

## Repository layout

```
src/        12 modules, one per paper algorithm (see Status table below)
scripts/    9 driver scripts: fetch data, build a reference, align, run each
            pipeline stage against real data
tests/      Automated regression suite (pytest) — 44 tests across 12 files,
            including tests that specifically pin the real bugs found and
            fixed during this project so they can't silently regress
```

## Status

| Module | File | Fidelity | Notes |
|---|---|---|---|
| Read depth extraction | `src/depth.py` | Exact | Directly from BAM, no ambiguity |
| GC-content correction | `src/gc_correction.py` | Good | F_gc/N_gc rate model per Benjamini & Speed 2012, sparse-bin smoothing, subtract/divide both supported. No mappability-track filtering on the background sample; one aDNA-margin parameter left at an unstated default |
| Nucleosome calling + scoring | `src/nucleosome_calling.py` | Exact | Formula fully specified in Methods. Flank width/offset refined against Hanghøj et al. 2016 Fig. 1B: 25bp flanks with a 12bp gap from the 147bp window (was an unconfirmed 50bp/adjacent guess) |
| Control set (FDR baseline) | `src/control_set.py` | Exact logic | Needs a second modern WGS BAM to sample from (any small one works — see QUICKSTART.md) |
| Periodicity (Welch, STFT, phasograms) | `src/periodicity.py` | Exact | Standard DSP, fully specified. `phasogram()`'s min-depth filter (a previously-dead parameter) is now actually implemented; default 5 (Pedersen 2014's stated value) with 3 available to match Hanghøj et al. 2016's independent reimplementation, a real confirmed discrepancy between the two papers |
| Nucleotide/dinucleotide patterns | `src/nucleotide_patterns.py` | Exact | Mono + purine/pyrimidine dinucleotide positional frequencies across aligned dyads; stratify-by-score built in |
| Methylation (Ms score) | `src/methylation.py` | Exact | Direct reimplementation of Pedersen's own Ms formula (not mapDamage's damage-curve model, deliberately); correct forward/reverse-strand C→T vs G→A handling; CpA/CpT/CpC negative controls for Fig 4B-style validation. Site-exclusion refinement added (`ms_score_filtered`, per Hanghøj et al. 2016): sites with >50% flip rate at ≥5x coverage — likely genuine polymorphisms, not deamination — are excluded; original unfiltered `ms_score()` kept unchanged for literal Pedersen 2014 behavior |
| CTCF anchor analysis | `src/ctcf_analysis.py` | Exact | Real Fu et al. 2008 occupied/unoccupied site table (strand-paired); strand-oriented 25bp-bin averaging matches Methods |
| Age-at-death | `src/age_estimation.py` | Exact | Two-layer pipeline (Ms→beta calibrated from modern donors, then Koch & Wagner 2011 beta→age) with real Figure 3A regression constants for both CpGs Pedersen used (cg07533148/TRIM58, cg01530101/KCNQ1DN) |
| Expression inference | `src/expression_inference.py` | Exact | Rs (gene body/promoter Ms ratio), +1 nucleosome occupancy, phasing strength (reuses `periodicity.py`), plus GSE3058 fetch + Spearman-correlation evaluation |
| Functional enrichment | `src/functional_enrichment.py` | Different backend | Enrichr (via `gseapy`) substituting DAVID — same purpose, not a reimplementation of DAVID itself |
| Genomic annotations | `src/genomic_annotations.py` | Exact | UCSC hg18 refGene/cpgIslandExt parsing: TSS + splice-site anchors (reuses `ctcf_analysis.py`'s generic `anchor_profile`), gene-body/CpG-island regions for the paper's non-anchored periodogram pattern |

## Scripts

```
fetch_sra_data.py            SRA accession -> raw FASTQ reads
prepare_reference.py         hg18 reference FASTA, downloaded + indexed
align_reads.py                FASTQ + reference -> sorted, indexed BAM
run_occupancy_pipeline.py    depth -> GC-correct -> call -> periodicity -> nucleotide patterns
run_methylation_pipeline.py  Ms score + negative controls
run_ctcf_pipeline.py         Ms/depth profiles around real CTCF sites
run_expression_pipeline.py   Rs, +1 occupancy, phasing strength
run_annotation_pipeline.py   TSS/splice-site/CGI/gene-body analyses
test_age_estimation.py       Standalone age-estimation demo (pure math, no data needed)
```

Real accessions confirmed and working: Saqqaq = `SRA010102`, Aboriginal
Australian = `SRA035301.1` (old-format SRA study accessions — resolved to
run-level SRR accessions automatically). Modern comparison data (for
`control_set.py`) — 1000 Genomes low-coverage, e.g. `ERR001268`.

## Install

```bash
brew install sra-tools bwa samtools    # external tools, not pip-installable
pip install -r requirements.txt
```

## Known, deliberate limitations

- **Full-genome coverage was not attempted** — ~2.5 billion raw reads / ~75GB
  minimum download, a hard bandwidth constraint, not a code limitation.
- **Real-data test results are near-empty/noise-level throughout** — an
  expected, explained consequence of validating at small scale (thousands of
  reads) against a dataset whose real scale is billions. See
  `PROJECT_SUMMARY.md` for the actual math on why this is still a meaningful
  validation, not a failure.
- Full per-module caveats are in the Status table above; nothing here is
  hidden or glossed over.

## Bugs found and fixed during this project (see `tests/` for the regression tests)

1. GC-correction rate model had a structurally wrong denominator (used
   observed reads for both numerator and denominator — couldn't detect bias).
2. Periodicity module's default FFT resolution couldn't resolve the paper's
   central 193bp periodicity, silently reporting a wrong value instead.
3. Inconsistent null-result handling between modules (`methylation.py`
   correctly returned NaN on no-data; `periodicity.py` didn't, until fixed).
4. UCSC `cpgIslandExt` table schema mismatch (no leading `bin` column, unlike
   `refGene`).

All four are now locked in by named regression tests in `tests/`.

## Refinements from Hanghøj et al. 2016 (epiPALEOMIX)

This is an independent reimplementation of Pedersen et al. 2014's exact
methods (same GCcorrect, NucleoMap, MethylMap, phasogram approach), published
by an overlapping author group. Rather than duplicate our own pipeline, we
used it to resolve several parameters our own modules had left as documented,
unconfirmed guesses:

- **Nucleosome flank geometry** (`nucleosome_calling.py`) — 25bp flanks with a
  12bp gap from the 147bp window (Fig. 1B), replacing an earlier 50bp/adjacent
  guess.
- **Phasogram min-depth** (`periodicity.py`) — confirmed a genuine discrepancy:
  Pedersen 2014 states 5 reads minimum; Hanghøj 2016's own reimplementation of
  the same method uses 3. Both are supported, defaulting to Pedersen's 5.
- **Ms site-exclusion rule** (`methylation.py`) — sites with >50% C→T flip rate
  at ≥5x coverage are excluded (likely genuine polymorphisms, not deamination).
  Added as `ms_score_filtered()`, alongside the original unfiltered `ms_score()`.

**Scoped but not yet built** (real, feasible, but each needs a new external
data/tool dependency — a scope decision, not started without discussion):
- **Mappability filtering** — the paper restricts nucleosome calls to 20kb
  blocks with ≥0.9 mappability uniqueness (41-mers, Derrien et al. 2012).
  Needs a public mappability track (e.g. UCSC's ENCODE mappability tracks) as
  a practical substitute for the exact 41-mer tool used.
- **Deamination-rate module** — `ds * ((1/k) - 1)` using mapDamage2's own `ds`/`k`
  output. A bounded, specific use of mapDamage2 (just two summary parameters),
  distinct from the full damage-curve model this project deliberately avoided
  earlier.

## Novel modules considered, not built

- **WPS (Windowed Protection Score)** — a genuinely different nucleosome-
  calling algorithm (Snyder et al. 2016) this paper compares against NucleoMap.
  This paper confirms the window sizes tested (10-120bp) but not the exact
  scoring formula (cited from Snyder et al. 2016, not reproduced in full here)
  — would need that source paper before building, same approach as finding
  Koch & Wagner's Figure 3A earlier.
- **Horvath's DNAmAge clock** (353-CpG epigenetic age predictor) — likely
  genuinely infeasible the way Koch & Wagner's 2-CpG model was tractable: this
  needs a full external supplementary coefficient table (353 rows), not a
  couple of numbers readable off one figure. Tentatively a documented
  limitation, same category as full-genome coverage, unless that table turns
  out to be findable.