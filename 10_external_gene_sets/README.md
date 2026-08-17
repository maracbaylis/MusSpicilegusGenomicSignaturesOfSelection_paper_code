# 10 — External behavior and plasticity gene-set enrichment

## Purpose

Compare *M. spicilegus* evolutionary candidate genes with published gene sets
related to behavioral evolution, neural plasticity, synaptic maturation, and
other relevant neurobehavioral phenotypes.

This workflow is separate from GO enrichment because the tested sets are
curated from external publications rather than Gene Ontology.

## Workflow

### 1. Build normalized external gene sets

`build_external_gene_sets.py` reads the published supplementary datasets, paths and citations are found in the Methods Section 11 (Comparison of M. spicilegus candidate genes with known conserved developmental programs underlying complex behavior)

and normalizes them into a common long-format table containing dataset, subset,
gene symbol, source, and evidence metadata.

```bash
python build_external_gene_sets.py
```

The script includes explicit source information for the external datasets it
parses. 

### 2. Build candidate-set overlaps

`compare_candidate_overlaps.py` reads the final *M. spicilegus* PAML, phyloP,
and Brownian-expression hit sets and reports the genes shared with each
external dataset/subset.

It also constructs combined candidate sets such as:

```text
Any_candidate_evidence
Multi_evidence_any_two_or_more
PhyloP_plus_BM
PAML_plus_PhyloP
PAML_plus_BM
All_three_evidence_types
```

The final hit workbook is supplied through `MSPIC_FINAL_HITS_WORKBOOK` or the
repository's configured input directory.

```bash
MSPIC_FINAL_HITS_WORKBOOK=<final-hits-workbook.xlsx> \
python compare_candidate_overlaps.py
```

### 3. Test enrichment

`test_external_gene_set_enrichment.py` performs one-sided Fisher/hypergeometric
enrichment tests for the candidate/external-set overlaps using the appropriate
tested-gene backgrounds and applies Benjamini-Hochberg correction.

```bash
MSPIC_FINAL_HITS_WORKBOOK=<final-hits-workbook.xlsx> \
python test_external_gene_set_enrichment.py
```

For the phyloP candidate set, the background should consist of genes linked to
valid cCREs eligible for the phyloP analysis. For Brownian-expression
candidates, the background should consist of genes tested in the corresponding
expression analysis. Combined candidate analyses should use the corresponding
common/testable background implemented by the script.

## Repository contents

```text
10_external_gene_sets/
├── README.md
├── build_external_gene_sets.py
├── compare_candidate_overlaps.py
├── test_external_gene_set_enrichment.py
└── raw/
    └── README.md
```

## Requirements

- Python >= 3.10
- pandas
- openpyxl
