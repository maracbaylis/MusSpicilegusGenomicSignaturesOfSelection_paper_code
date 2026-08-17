# 09 — GO Biological Process enrichment

## Purpose

Test whether genes linked to *M. spicilegus* evolutionary candidate loci are
enriched for Gene Ontology Biological Process (GO:BP) terms.

The primary manuscript analysis uses a custom gene universe consisting of genes
linked to the analyzed cCRE background rather than the entire mouse genome.

## Primary analysis

`run_combined_spicilegus_fg_go_bp_allcCRE_background.R` implements the primary
combined-foreground GO:BP analysis used for the manuscript result set,
including the reported `memory` enrichment.

The script:

1. defines the combined *M. spicilegus* foreground gene set;
2. reads the all-cCRE linked-gene background;
3. maps mouse gene symbols to GO Biological Process annotations using
   `org.Mm.eg.db`;
4. performs one-sided Fisher exact enrichment tests for GO:BP terms; and
5. applies Benjamini-Hochberg correction across tested GO terms.

Run with project/input directories supplied through the environment variables
used throughout the repository:

```bash
MSPIC_PROJECT_ROOT=<repository-root> \
MSPIC_INPUT_DIR=<input-data-directory> \
Rscript run_combined_spicilegus_fg_go_bp_allcCRE_background.R
```

## Reusable custom gene-list analysis

`run_go_bp_enrichment_allccre_neuro_background.R` provides the corresponding
custom-input workflow for an arbitrary gene list while retaining the all-cCRE
background.

It expects:

```text
go_bp_input_genes.txt
go_bg_allcCRE_genes.txt
```

and performs one-sided Fisher exact GO:BP enrichment followed by BH correction.


## Repository contents

```text
09_GO_enrichment/
├── README.md
├── run_combined_spicilegus_fg_go_bp_allcCRE_background.R
└── run_go_bp_enrichment_allccre_neuro_background.R
```

## Requirements

- R 4.1.1
- AnnotationDbi
- org.Mm.eg.db
- GO.db
- dplyr
- readr
- stringr
- ggplot2
