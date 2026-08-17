#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript voom_species_means.R ",
    "<counts.tsv> <sample_metadata.tsv> <output_prefix>"
  )
}

counts_file <- args[[1]]
metadata_file <- args[[2]]
output_prefix <- args[[3]]

counts <- read.delim(
  counts_file,
  row.names = 1,
  check.names = FALSE
)

metadata <- read.delim(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata <- c("sample", "species")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

missing_samples <- setdiff(metadata$sample, colnames(counts))
if (length(missing_samples) > 0) {
  stop(
    "Samples in metadata but not count matrix: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Reorder counts to match metadata.
counts <- counts[, metadata$sample, drop = FALSE]

# The manuscript describes voom-normalized log-expression values.
# A simple intercept-only design is used here to estimate the voom mean-variance transformation without imposing species contrasts.
design <- matrix(1, nrow = ncol(counts), ncol = 1)
rownames(design) <- colnames(counts)
colnames(design) <- "Intercept"

voom_obj <- voom(counts, design = design, plot = FALSE)
voom_expression <- voom_obj$E

write.table(
  cbind(gene_name = rownames(voom_expression), voom_expression),
  paste0(output_prefix, "_voom_expression.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

species_order <- unique(metadata$species)

species_means <- sapply(
  species_order,
  function(species_name) {
    sample_names <- metadata$sample[metadata$species == species_name]
    rowMeans(voom_expression[, sample_names, drop = FALSE])
  }
)

species_means <- as.data.frame(species_means, check.names = FALSE)
species_means <- cbind(
  gene_name = rownames(voom_expression),
  species_means
)

write.csv(
  species_means,
  paste0(output_prefix, "_species_mean_voom_expression.csv"),
  row.names = FALSE,
  quote = FALSE
)

message("Wrote voom matrix and species-mean expression files.")
