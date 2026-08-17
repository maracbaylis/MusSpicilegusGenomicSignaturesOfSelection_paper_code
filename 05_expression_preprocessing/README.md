# 05 — Expression preprocessing

## Purpose

Prepare cross-species gene-expression values for the downstream
Brownian-motion expression analysis.

The manuscript describes the following preprocessing sequence:

1. align RNA-seq reads from each species to its corresponding genome assembly;
2. project the mouse reference annotation onto each species genome using
   Liftoff v1.6.3 with minimap2 v2.28;
3. generate gene-level counts using HTSeq v2.0.3;
4. retain genes successfully identified across all analyzed taxa;
5. normalize the shared count matrix with the `voom` procedure in limma
   v3.48.3; and
6. average voom-normalized log-expression across biological samples within
   each species.

The output used by the Brownian-motion analysis is a matrix containing one row
per gene and one species-mean expression column per taxon.

## Cross-species annotation

Mouse reference annotations were projected onto the corresponding species
assemblies with Liftoff v1.6.3, using minimap2 v2.28.

## Gene-level counting

Gene-level counts were generated from the species-specific RNA-seq alignments
using HTSeq v2.0.3 and the projected species annotation.

## Shared gene universe

Genes were matched across species using orthologous gene symbols, and only
genes successfully identified in all analyzed taxa were retained for
comparative expression analyses.

The resulting count matrix should contain one row per shared gene and one
column per RNA-seq sample.

## voom normalization and species means

`voom_species_means.R` takes:

```text
counts.tsv
sample_metadata.tsv
```

where `counts.tsv` has gene IDs/names in the first column and sample counts in
the remaining columns, and `sample_metadata.tsv` contains:

```text
sample    species
```

Run:

```bash
Rscript voom_species_means.R \
    shared_gene_counts.tsv \
    sample_metadata.tsv \
    expression
```

The script writes:

```text
expression_voom_expression.tsv
expression_species_mean_voom_expression.csv
```

The species-mean CSV is the direct input to the Brownian-motion analysis in
`06_brownian_expression/`.

## Repository contents

```text
05_expression_preprocessing/
├── README.md
└── voom_species_means.R
```

## Requirements

- R 4.1.1
- limma 3.48.3
- Liftoff 1.6.3
- minimap2 2.28
- HTSeq 2.0.3
