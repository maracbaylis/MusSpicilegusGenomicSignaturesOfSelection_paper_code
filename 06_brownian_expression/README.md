# 06 — Brownian-motion expression divergence

## Purpose

Identify genes whose species-mean expression deviates from the expectation
under a Brownian-motion model of expression evolution.

The input to this folder is produced by
`05_expression_preprocessing/voom_species_means.R`.

## Input

`run_brownian_expression.py` expects one row per gene with these species-mean
voom-expression columns:

```text
gene_name
musculus
castaneus
domesticus
spicilegus
spretus
```

## Brownian-motion model

For each gene and each focal species, the observed expression value is compared
with its conditional expectation given expression in the remaining species and
the fixed phylogenetic covariance matrix used in the analysis.

The script reports:

- observed expression;
- Brownian conditional expected expression;
- observed-minus-expected residual;
- conditional variance and standard deviation;
- standardized residual z-score;
- analytic two-sided normal P-value and BH-adjusted q-value; and
- empirical P-value from 200,000 draws from the absolute standard normal
  distribution and its BH-adjusted q-value.

Multiple-testing correction is performed separately within each focal species.

## Run the Brownian-motion analysis

```bash
python run_brownian_expression.py \
    ../05_expression_preprocessing/expression_species_mean_voom_expression.csv \
    bm_conditional_residuals_all_species_all_genes.csv
```

## Stringent focal-lineage outliers

The manuscript's stringent expression-divergence criteria are:

```text
q < 0.01
|z| >= 3
focal |z| - next-most-extreme species |z| >= 2
```

For *M. spicilegus*:

```bash
python filter_expression_outliers.py \
    bm_conditional_residuals_all_species_all_genes.csv \
    spicilegus_expression_outliers.csv \
    --foreground spicilegus
```

The filter script can use either analytic or empirical BH-adjusted q-values.
The repository and manuscript should explicitly identify which q-value
definition produced the reported final candidate set.

## Repository contents

```text
06_brownian_expression/
├── README.md
├── run_brownian_expression.py
└── filter_expression_outliers.py
```

## Requirements

- Python >= 3.8
- numpy
- pandas
- scipy
