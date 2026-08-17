# Software versions

This file records the versions used or targeted for the manuscript analyses.
`environment.yml` is kept as the installable Conda environment; some older
Bioconductor package builds are not always available from Conda with exact
version pins on every platform.

## Core Python environment

| Package | Version |
|---|---:|
| Python | 3.10 |
| pandas | 2.1.1 |
| numpy | 1.25.2 |
| scipy | 1.15.3 |
| statsmodels | 0.14.0 |
| openpyxl | 3.1.5 |
| matplotlib | 3.8.0 |
| seaborn | 0.13.0 |
| requests | 2.32.3 |
| biopython | 1.81 |
| pillow | 11.3.0 |
| certifi | 2025.4.26 |

## Core R/Bioconductor environment

| Package | Version |
|---|---:|
| R | 4.1.1 |
| limma | 3.48.3 |
| AnnotationDbi | 1.54.1 |
| org.Mm.eg.db | 3.13.0 |
| GO.db | 3.13.0 |
| dplyr | 1.1.2 |
| readr | 2.1.5 |
| stringr | 1.5.2 |
| ggplot2 | 4.0.0 |

## External command-line tools

| Tool | Version / source |
|---|---|
| Ensembl Compara EPO alignments | release-112, 21-way Murinae EPO |
| Ensembl Compara `emf2maf.pl` | from Ensembl Compara scripts |
| PAML/codeml | 4.10.7 |
| PHAST/phyloP/phyloFit | 1.5 |
| MAFFT | 7.526 |
| PAL2NAL | 14 |
| Liftoff | 1.6.3 |
| minimap2 | 2.28 |
| HTSeq | 2.0.3 |
| TensorFlow | 2.21.0 |
| TensorFlow Hub | 0.16.1 |
