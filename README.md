# Mus spicilegus paper code

Code supporting the manuscript:

> Genomic signatures of selection in the mound-building mouse *Mus spicilegus*

This repository contains the analysis and figure-generation scripts used to
identify coding, regulatory, and expression-divergence candidates in
*Mus spicilegus* and related *Mus* lineages.

## Repository organization

The workflow is organized as numbered analysis modules:

| Directory | Purpose |
|---|---|
| `01_prepare_orthologs_and_alignments/` | Extract and stitch locus-level alignments from Ensembl EPO multiple alignments. |
| `02_alignment_qc/` | Apply alignment artifact filters before PAML and phyloP analyses. |
| `03_run_paml_all_foregrounds/` | Run codeml branch-model tests and Benjamini-Hochberg correction. |
| `04_run_phylop_all_foregrounds/` | Run phyloP regulatory acceleration tests and FDR correction. |
| `05_expression_preprocessing/` | Normalize RNA-seq counts and calculate species-mean expression. |
| `06_brownian_expression/` | Test for lineage-specific expression divergence under a Brownian-motion model. |
| `07_motif_substitution/` | Test transcription-factor motif effects of *M. spicilegus*-specific substitutions. |
| `08_enformer_predictions/` | Build Enformer inputs and summarize predicted regulatory effects. |
| `09_GO_enrichment/` | Run GO Biological Process enrichment tests. |
| `10_external_gene_sets/` | Compare candidate genes with external behavior and plasticity gene sets. |
| `11_figures/` | Recreate manuscript and supplementary figure panels from analysis outputs. |

Each subdirectory has its own README with inputs, commands, and expected
outputs for that stage.

## Quick start

Create the software environment:

```bash
conda env create -f environment.yml
conda activate mspic-paper
```

The Conda file installs the Python stack. Exact R/Bioconductor package versions
used for the manuscript are recorded in `SOFTWARE_VERSIONS.md`; older
Bioconductor builds are not consistently available from Conda across platforms.

Some analyses also require external command-line tools that are not installed
by the Conda environment:

- Ensembl Compara `emf2maf.pl`
- PAML/codeml 4.10.7
- PHAST/phyloP/phyloFit 1.5
- MAFFT 7.526
- PAL2NAL 14
- Liftoff 1.6.3
- minimap2 2.28
- HTSeq 2.0.3

See `DATA_MANIFEST.md` for the input data expected by each analysis section.
See `SOFTWARE_VERSIONS.md` for the exact manuscript software versions.

## Data availability

Large public reference files and manuscript input/output tables are not stored
directly in this lightweight code repository. The expected locations and source
descriptions are listed in `DATA_MANIFEST.md`.

By convention, scripts that need shared manuscript inputs look for:

```bash
MSPIC_PROJECT_ROOT=<repository-root>
MSPIC_INPUT_DIR=<input-data-directory>
```

If these variables are not set, many scripts fall back to `data/inputs` under
the current working directory.

## Reproducibility notes

- Coding-sequence analyses use GRCm39-reference alignments.
- cCRE analyses use mm10/GRCm38 coordinates from ENCODE SCREEN.
- PAML results are corrected within foreground lineage using the
  Benjamini-Hochberg procedure.
- phyloP results are corrected within each tested foreground/tissue/timepoint
  set using the Benjamini-Hochberg procedure.
- Brownian-expression tests use analytic normal p-values and empirical
  p-values from 200,000 absolute standard-normal draws.
- Motif analyses use GC/base-change-matched permutation tests with 20,000
  permutations.

## Release notes

- Ensembl Compara EPO alignments are from release-112.
- The phyloP wrapper uses `--method SPH --mode CONACC --subtree`.
- Motif-scanning scripts read a JASPAR-format PWM file from `MSPIC_JASPAR_PATH`
  or from `MSPIC_INPUT_DIR/JASPAR.jaspar`.

## License

Code in this repository is released under the MIT License. See `LICENSE`.
