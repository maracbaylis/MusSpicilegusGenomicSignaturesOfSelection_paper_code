# 07 — Transcription-factor motif substitution analysis

## Purpose

Test whether *M. spicilegus*-specific substitutions in accelerated brain-active
cCREs alter predicted transcription-factor binding architecture and whether
motif-associated changes exceed expectations from matched background cCREs and
local sequence composition.

The analysis compares:

- 285 *M. spicilegus*-accelerated cCRE alignments; and
- 285 matched non-accelerated brain-active cCRE alignments.

The matched background set should be the same set used in the manuscript,
matched for cCRE length, GC content, and *M. spicilegus*-specific substitution
burden.

## Motif scoring

For each cCRE alignment, the non-*M. spicilegus* sequences are used to infer a
consensus ancestral/non-*M. spicilegus* state. Sites at which
*M. spicilegus* differs from this consensus are evaluated in local sequence
windows with JASPAR position-weight matrices.

Motif PWMs are converted to log-odds PSSMs using a uniform nucleotide
background and pseudocounts of 0.5. Motif scores are converted to relative
scores using each motif's minimum and maximum possible score.

The supplied motif-scanning code uses a relative-score threshold of:

```text
0.85
```

A motif-switch site is one at which the best-scoring motif differs between the
ancestral/non-*M. spicilegus* and *M. spicilegus* sequence states and at least
one state passes the relative-score threshold.

## Core analyses

### 1. Motif-family and motif-identity switches

`run_motif_analysis.py` scans the accelerated and matched-background MSA
sets, constructs the site-level substitution/motif table, identifies motif
switches, and performs GC/base-change-matched permutation tests of
motif-to-motif transitions.

```bash
python run_motif_analysis.py
```

The analysis uses 20,000 permutations and applies Benjamini-Hochberg
multiple-testing correction to the matched permutation P-values.

`scan_motif_family_switches_msa_sets.py` contains the shared motif-scoring and
TF-family helper functions used by this analysis.

### 2. Motif-overlap enrichment

`test_gc_corrected_motif_overlap.py` tests whether
*M. spicilegus*-specific substitutions occur within predicted motif bases more
often than expected after matching for GC context.

Local GC is calculated in a window around each substitution and whole-cCRE GC
is calculated from the inferred ancestral sequence. These variables are used
to define GC-matched permutation strata.

```bash
python test_gc_corrected_motif_overlap.py
```

### 3. Destination-motif enrichment

`destination_motif_gc_matched_null.py` tests whether particular
destination motifs are overrepresented among accelerated cCRE motif-switch
sites relative to the GC/base-change-matched background.

```bash
python destination_motif_gc_matched_null.py
```

### 4. De novo motif gains

`de_novo_motif_gains.py` identifies substitutions that generate a
motif match in the *M. spicilegus* state that is absent from the inferred
ancestral state and tests destination motifs with the same matched permutation
framework.

```bash
python de_novo_motif_gains.py
```

### 5. Motif-score effect size

`test_mutation_motif_score_delta.py` quantifies the change in motif score caused
by individual substitutions and compares accelerated cCREs with the matched
background set.

```bash
python test_mutation_motif_score_delta.py
```

## Required inputs

The analysis expects:

```text
selected_MSAs/
background_MSAs/
JASPAR motif file
```

The public repository should document the exact JASPAR release/file used.
The supplied scripts reference a JASPAR combined-matrix file.

The selected and background MSA directories should contain the exact 285 + 285
alignment sets used for the reported analysis, or the repository should provide
a reproducible script/table defining those sets.

## GC/base-change matching

The shared helper `test_family_transition_gc_matched_null.py` defines the
matching strata used by the permutation analyses.

Sites are stratified by:

- substitution base-change class;
- quartile of local GC content in the inferred ancestral/non-*M. spicilegus*
  sequence; and
- quartile of whole-cCRE GC content.

GC quartiles are generated with `pandas.qcut(..., q=4, duplicates="drop")`.
Within each stratum, accelerated/background labels are permuted while
preserving the number of selected sites. The motif analyses use 20,000
permutations and fixed random seeds for reproducibility.

The same helper also provides the Benjamini-Hochberg FDR correction used by
the motif-transition tests.

## Repository contents

```text
07_motif_substitution/
├── README.md
├── run_motif_analysis.py
├── scan_motif_family_switches_msa_sets.py
├── test_gc_corrected_motif_overlap.py
├── test_mutation_motif_score_delta.py
├── de_novo_motif_gains.py
├── destination_motif_gc_matched_null.py
└── test_family_transition_gc_matched_null.py
```

## Requirements

- Python >= 3.8
- numpy
- pandas
- scipy
- biopython
- statsmodels
- openpyxl (optional; only for Excel output)
