# 03 — PAML branch-model analysis

## Purpose

Run lineage-specific PAML `codeml` branch-model tests and identify genes with evidence of elevated foreground dN/dS.

## Model

For each codon alignment:

- **Null model:** one-ratio branch model (`model = 0`, `NSsites = 0`), in which all branches share a single omega.
- **Alternative model:** two-ratio foreground branch model (`model = 2`, `NSsites = 0`), in which the focal foreground branch has a separate omega.

The likelihood-ratio statistic is:

```text
LRT = 2 * (lnL_alternative - lnL_null)
```

P-values are calculated from a chi-square distribution with 1 degree of freedom and corrected across genes within each foreground using the Benjamini-Hochberg procedure.

Candidate genes satisfy:

```text
q_value_bh <= 0.05
omega_foreground > 1
delta_lnL > 0
```

This is a branch-model analysis, not a branch-site analysis.

## Workflow

### 1. Prepare codon alignments

Codon alignments are supplied to PAML in PHYLIP format (`.phy`).

### 2. Run the null model

The null model is run using `codeml_template_null.ctl` and the all-background tree:

```bash
python run_codeml_directory.py \
    PAML_input_genes \
    codeml_output_null \
    trees/null_all_background.tree \
    templates/codeml_template_null.ctl
```

### 3. Run foreground branch models

The alternative model is run separately for each foreground lineage using `codeml_template_model2.ctl` and the corresponding foreground-labeled tree.

For example, for *Mus spicilegus*:

```bash
python run_codeml_directory.py \
    PAML_input_genes \
    codeml_output_spicilegus \
    trees/mus_spicilegus.tree \
    templates/codeml_template_model2.ctl
```

The same procedure is repeated for each focal lineage using its corresponding tree file.

### 4. Process likelihood-ratio tests

Null and alternative `codeml` results are summarized and processed with:

```bash
python process_paml_lrt_bh.py
```

For each foreground, the script calculates:

- the difference in log likelihood between alternative and null models;
- the likelihood-ratio statistic;
- chi-square P-values with 1 degree of freedom;
- Benjamini-Hochberg-adjusted q-values across genes; and
- the final candidate classification based on q-value and foreground omega.

## Repository contents

```text
03_paml/
├── README.md
├── run_codeml_directory.py
├── process_paml_lrt_bh.py
├── templates/
│   ├── codeml_template_null.ctl
│   └── codeml_template_model2.ctl
└── trees/
    ├── null_all_background.tree
    ├── mus_caroli.tree
    ├── mus_musculus_casteij.tree
    ├── mus_musculus_pwkphj.tree
    ├── mus_musculus_wsbeij.tree
    ├── mus_pahari.tree
    ├── mus_spicilegus.tree
    └── mus_spretus.tree
```

The foreground tree files specify the focal lineage for each two-ratio branch-model analysis. `null_all_background.tree` contains the corresponding tree without a foreground branch designation.

## Requirements

- PAML/codeml 4.10.7
- Python >= 3.8
- pandas

`codeml` must be available on `PATH`.

## Citation

Yang Z. 2007. PAML 4: Phylogenetic Analysis by Maximum Likelihood. *Molecular Biology and Evolution* 24:1586–1591.