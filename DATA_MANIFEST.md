# Data manifest

This manifest explains which files support the manuscript and how those files
relate to the code in this repository.

The repository contains analysis code, small configuration files, model files,
tree files, and documentation. 

## Data availability

Supplemental Tables S1-S9 are provided with the manuscript. Large public
reference datasets should be obtained from the original sources listed below.
Intermediate working files are available from the corresponding author upon
reasonable request, unless they are already included in the manuscript
supplements or can be regenerated from public sources.

If additional intermediate files are deposited in a public archive, add the DOI
or URL here:

```text
Additional data archive: not applicable
```

## Repository input convention

Several scripts use environment variables to locate data without relying on a
specific local directory structure:

```bash
MSPIC_PROJECT_ROOT=<repository-root>
MSPIC_INPUT_DIR=<input-data-directory>
MSPIC_JASPAR_PATH=<jaspar-pwm-file>
```

If `MSPIC_INPUT_DIR` is not set, many scripts fall back to `data/inputs`.

Recommended local layout for rerunning analyses:

```text
data/
├── inputs/
├── outputs/
└── supplemental_tables/
```

These directories are ignored by git so that large data files do not get
committed accidentally.

## Supplemental result files

These are the manuscript-supporting supplemental tables.

| Table | Filename | Sheets / contents | Description |
|---|---|---|---|
| Supplementary Table 1 | `Table_S1_full_PAML_all_foregrounds.xlsx` | `full_PAML_all_foregrounds` | PAML tests across lineages. Rows report one gene and foreground species, with gene IDs/names, foreground labels, null and alternative log-likelihoods, delta lnL, LRT, nominal p-value, Benjamini-Hochberg q-value, background and foreground omega, q < 0.05 status, and final reference-hit status. |
| Supplementary Table 2 | `Table_S2_Allen_Brain_Atlas_expression_summaries.xlsx` | `gene_coverage_top_region`; `gene_top_structures` | Allen Brain Atlas expression summaries. Includes gene coverage/top-region metadata and expression patterns across brain regions. |
| Supplementary Table 3a | `Table_S3_a_phyloP_whole_brain.xlsx` | `All_accelerated_CREs_genenames`; full phyloP tabs for spicilegus, spretus, domesticus, musculus, and castaneus | Whole-brain cCRE phyloP results for the full species set, including significant accelerated cCREs and complete per-foreground outputs. |
| Supplementary Table 3b | `Table_S3_b_phyloP_whole_brain_without_spicilegus.xlsx` | `phyloP summary` | Whole-brain cCRE phyloP control results after removing the *M. spicilegus* sequence from each alignment. |
| Supplementary Table 4 | `Table_S4_GO_gProfiler_enrichment.xlsx` | `combined_GO_all_terms` | Gene Ontology Biological Process enrichment for genes associated with whole-brain cCREs accelerated in *M. spicilegus*. |
| Supplementary Table 5 | `Table_S5_phyloP_all_foregrounds_across_brain_tissues_and_timepoints.zip` | contains `Combined_PhyloP_Summary.xlsx` | Brain-region and developmental-timecourse cCRE phyloP results across foreground species, tissues, and timepoints. |
| Supplementary Table 6 | `Table_S6_Enformer_allele_replacement_predictions.xlsx` | `top10_targets_by_cCRE`; `phyloP_accelerated_CREs`; `matched_background_CREs` | Enformer allele-replacement predictions for accelerated cCREs and matched-background cCREs. |
| Supplementary Table 7 | `Table_S7_TF_motif_substitution_analysis.xlsx` | `combined_substitution_context` | Transcription-factor motif substitution analysis for cCREs with *M. spicilegus* variants and matched controls. |
| Supplementary Table 8 | `Table_S8_full_BM_expression_all_species.xlsx` | species-by-tissue Brownian-motion tabs | Brownian-motion expression outlier results across species and tissues. |
| Supplementary Table 9 | `Table_S9_external_gene_set_overlap_enrichment.xlsx` | `wang2023_synapse_enrichment`; `spark_enrichment`; `johnson_cichlid_overlap` | External gene-set overlap and enrichment analyses for synaptic maturation, autism-associated CNV genes, and cichlid behavior gene sets. |

## Public reference inputs

| Data type | Species/resource | Version/accession/URL/DOI | Date accessed | Used by |
|---|---|---|---|---|
| Multiple alignment | 21 Murinae EPO alignment | Ensembl Compara release 112 MLSS 9598: `https://useast.ensembl.org/info/genome/compara/mlss.html?mlss=9598`; FTP README: `https://ftp.ensembl.org/pub/release-112/emf/ensembl-compara/multiple_alignments/21_murinae.epo/README.21_murinae.epo`; file: `Compara.112.ncrna_murinae.aln.emf.gz` | 2026-07-20 | phyloP alignment inputs; PAML/ortholog alignment preprocessing where applicable |
| Genome assembly | *Mus musculus* | Ensembl Compara release 112 21 Murinae EPO source assembly; see MLSS 9598 metadata above | 2026-07-20 | orthologs, phyloP, PAML |
| Genome assembly | *Mus spicilegus* | Ensembl Compara release 112 21 Murinae EPO source assembly; see MLSS 9598 metadata above | 2026-07-20 | orthologs, phyloP, PAML |
| Genome assembly | *Mus spretus* | Ensembl Compara release 112 21 Murinae EPO source assembly; see MLSS 9598 metadata above | 2026-07-20 | orthologs, phyloP, PAML |
| Genome assembly | *Mus caroli* | Ensembl Compara release 112 21 Murinae EPO source assembly; see MLSS 9598 metadata above | 2026-07-20 | PAML |
| Genome assembly | *Mus pahari* | Ensembl Compara release 112 21 Murinae EPO source assembly; see MLSS 9598 metadata above | 2026-07-20 | PAML |
| Ortholog table/resource | all species | Ensembl BioMart, Ensembl release 112; one-to-one/ortholog gene mapping tables exported from BioMart for the species in the analysis | 2026-07-20 | PAML/BM gene mapping |
| cCRE/enhancer annotations | mouse developmental tissues | SCREEN Registry/Weng Lab: `https://screen.wenglab.org/` | 2026-07-20 | phyloP, motif, Enformer |
| RNA-seq count/expression tables | wild mouse tissues/species | Edmond dataset DOI: `https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.IMARWT`; ENA projects PRJEB50011 and PRJEB54000: `https://www.ebi.ac.uk/ena/browser/view/PRJEB50011`, `https://www.ebi.ac.uk/ena/browser/view/PRJEB54000` | 2026-07-20 | BM expression |
| GO annotations | mouse GO Biological Process | Gene Ontology enrichment analysis framework: `https://geneontology.org/docs/go-enrichment-analysis/`; implemented locally with R/Bioconductor packages `org.Mm.eg.db` and `GO.db` | 2026-07-20 | GO enrichment |
| TF motif database | vertebrate transcription factor motifs | JASPAR motif matrices from `https://jaspar.elixir.no/`; local reruns should provide this as `MSPIC_JASPAR_PATH` or as `JASPAR.jaspar` under `MSPIC_INPUT_DIR` | 2026-07-20 | motif/substitution analysis |
| Enformer model/assets | model and target metadata | DeepMind Enformer implementation: `https://github.com/google-deepmind/deepmind-research/blob/master/enformer/enformer.py`; Avsec Z. et al., "Effective gene expression prediction from sequence by integrating long-range interactions" | 2026-07-20 | Enformer predictions |
| Allen Mouse Brain Atlas | adult mouse brain expression atlas | Allen Mouse Brain Atlas adult mouse ISH resource: `https://mouse.brain-map.org/static/atlas`; Allen RMA API endpoint used for final matrices: `https://api.brain-map.org/api/v2/data/query.json?criteria=` with `Gene`, `SectionDataSet`, and `StructureUnionize` queries (`products[id$eq1]`, `[failed$eqfalse]`, adult structure graph 1); primary citation: Lein et al. 2007 Nature, `https://doi.org/10.1038/nature05453` | 2026-07-20 | Allen brain expression |
| External gene sets | SPARK/autism, PSD/synaptic maturation, songbird neurogenesis, cichlid behavior | SPARK gene list: `https://sparkforautism.org/discover_article/spark-gene-list/`; Wang, Li, Kaifang Pang, Li Zhou, et al. 2023. "A Cross-Species Proteomic Map Reveals Neoteny of Human Synapse Development." *Nature* 622:112-119; York, Ryan A., Chinar Patil, Kawther Abdilleh, et al. 2018. "Behavior-Dependent Cis Regulation Reveals Genes and Pathways Associated with Bower Building in Cichlid Fishes." *PNAS* 115:E11081-E11090 | 2026-07-20 | external gene-set overlap |

## Regeneration inputs and intermediates

The files below are needed to rerun the full workflow but are not intended to
be committed to git.

| Analysis section | Input or intermediate | How it is produced or supplied |
|---|---|---|
| `01_prepare_orthologs_and_alignments` | Coding-gene interval files | Derived from gene annotations used for the coding-sequence analysis. |
| `01_prepare_orthologs_and_alignments` | `species_map.tsv` | Included in this repository. |
| `02_alignment_qc` | Stitched locus FASTA alignments | Generated by section 01. |
| `03_run_paml_all_foregrounds` | Codon alignments in PHYLIP format | Generated after translation, amino-acid alignment, PAL2NAL conversion, and QC. |
| `03_run_paml_all_foregrounds` | PAML tree and control templates | Included in this repository under `trees/` and `templates/`. |
| `04_run_phylop_all_foregrounds` | cCRE FASTA alignments | Generated by sections 01 and 02. |
| `04_run_phylop_all_foregrounds` | Neutral PHAST models | Included in this repository under `04_run_phylop_all_foregrounds/models/`. |
| `05_expression_preprocessing` | Gene-level count matrix and sample metadata | Generated from aligned RNA-seq reads and projected species annotations. |
| `06_brownian_expression` | Species-mean voom expression matrix | Generated by section 05. |
| `07_motif_substitution` | `selected_MSAs/` and `background_MSAs/` | Accelerated and matched-background MSA sets used for motif analyses. |
| `08_enformer_predictions` | Enformer window manifests and FASTA files | Generated by the scripts in section 08. |
| `09_GO_enrichment` | all-cCRE linked-gene background | Derived from the cCRE-to-gene annotation background used in the manuscript. |
| `11_figures` | Final result tables from sections 03-10 | Use Supplemental Tables S1-S9 or equivalent regenerated outputs. |

## Git tracking policy

The `.gitignore` excludes local data directories, public reference genomes,
raw alignments, sequencing files, binary workbooks, and generated figure files.
This keeps the GitHub repository lightweight and focused on code.
