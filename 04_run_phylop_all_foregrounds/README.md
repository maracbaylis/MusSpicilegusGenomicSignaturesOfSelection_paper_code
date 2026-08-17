# 04 — phyloP regulatory-acceleration analysis

## Purpose


Run lineage-specific phyloP tests on cCRE multiple-sequence alignments and
apply Benjamini-Hochberg multiple-testing correction within each tested
foreground/tissue/timepoint set.

## Statistical framework

The analysis uses the Siepel-Pollard-Haussler (SPH) method implemented in
PHAST/phyloP with a neutral phylogenetic model fitted separately with
`phyloFit`.

For a focal lineage, phyloP partitions the tree using `--subtree` and evaluates
substitution behavior on the focal subtree relative to the complementary
supertree.

The primary analysis used the SPH subtree test in CONACC mode (--method SPH --mode CONACC --subtree). With --subtree, phyloP partitions the phylogeny into the focal subtree and its complementary supertree and tests for conservation or acceleration in the focal subtree conditional on the complementary supertree under the neutral model.

## Repository contents

```text
04_phylop/
04_phylop/
├── README.md
├── run_phylop_directory.py
├── add_phylop_fdr.py
└── models/
    ├── neutral_model.mod
    └── neutral_model_wo_spic.mod
```

## Inputs

- per-cCRE FASTA multiple-sequence alignments from the alignment-preparation
  workflow
- a PHAST `.mod` neutral model produced by `phyloFit`
- the node name corresponding to the focal lineage in that model/tree

## Neutral models

Two PHAST neutral models are included:

- `models/neutral_model.mod` — neutral model used for the primary
  all-species phyloP analyses.
- `models/neutral_model_wo_spic.mod` — neutral model used for the
  control analysis in which *M. spicilegus* was excluded before testing
  acceleration in the remaining lineages.

The neutral models were generated with `phyloFit` from fourfold-degenerate
sites as described in the manuscript Methods.

## Running phyloP

`run_phylop_directory.py` runs phyloP across all cCRE alignments for a
specified foreground lineage.

The M. spicilegus foreground analysis can be run as:

```bash
python run_phylop_directory.py \
    cCRE_alignments/ \
    models/neutral_model.mod \
    mus_spicilegus \
    outputs/spicilegus/
```

## Multiple-testing correction

`add_phylop_fdr.py` reads phyloP result CSV files, applies
Benjamini-Hochberg correction independently to each `p_*` column, and writes
corresponding `q_*` columns plus an FDR summary.

```bash
python add_phylop_fdr.py \
    phyloP_results_csv/ \
    phyloP_results_fdr/
```

Significant accelerated elements are defined using the acceleration q-value
corresponding to the primary test at q < 0.05.

## Requirements

- PHAST/phyloP/phyloFit 1.5
- Python >= 3.8
