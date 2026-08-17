options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(GO.db)
  library(org.Mm.eg.db)
  library(dplyr)
  library(readr)
  library(stringr)
})

base_dir <- file.path(Sys.getenv("MSPIC_PROJECT_ROOT", "."), "go_bp_custom_gene_list")
gene_file <- file.path(base_dir, "go_bp_input_genes.txt")
background_file <- file.path(Sys.getenv("MSPIC_INPUT_DIR", "data/inputs"), "go_bg_allcCRE_genes.txt")
cluster_file <- file.path(Sys.getenv("MSPIC_INPUT_DIR", "data/inputs"), "GO_terms_with_redundancy_clusters.csv")
out_dir <- file.path(base_dir, "allcCRE_neuro_background")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_gene_list <- function(path) {
  genes <- unique(trimws(readLines(path, warn = FALSE)))
  genes[nzchar(genes)]
}

test_symbols <- read_gene_list(gene_file)
background_original <- read_gene_list(background_file)
fg_not_in_background <- setdiff(test_symbols, background_original)
background <- sort(unique(c(background_original, fg_not_in_background)))

annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = background,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GENENAME", "ENTREZID", "GOALL", "ONTOLOGYALL")
)

go_annot <- annot %>%
  distinct(SYMBOL, GOALL, ONTOLOGYALL) %>%
  filter(!is.na(SYMBOL), !is.na(GOALL), ONTOLOGYALL == "BP")

bg_genes <- sort(unique(go_annot$SYMBOL))
query_genes <- intersect(test_symbols, bg_genes)
unmapped <- setdiff(test_symbols, query_genes)

if (length(query_genes) == 0) {
  stop("No input genes mapped to all-cCRE background BP annotations.")
}

go_bp <- bind_rows(lapply(sort(unique(go_annot$GOALL)), function(go_id) {
  term_genes <- unique(go_annot$SYMBOL[go_annot$GOALL == go_id])
  overlap <- intersect(query_genes, term_genes)
  a <- length(overlap)
  if (a == 0) return(NULL)
  b <- length(query_genes) - a
  c <- length(term_genes) - a
  d <- length(bg_genes) - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  term_name <- tryCatch(as.character(AnnotationDbi::Term(GOTERM[[go_id]])), error = function(e) NA_character_)
  data.frame(
    GOID = go_id,
    Description = term_name,
    ontology = "BP",
    Count = a,
    query_n = length(query_genes),
    bg_term_n = length(term_genes),
    bg_n = length(bg_genes),
    GeneRatio = a / length(query_genes),
    BgRatio = length(term_genes) / length(bg_genes),
    pvalue = ft$p.value,
    odds_ratio = unname(ft$estimate),
    geneID = paste(sort(overlap), collapse = "/"),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
  arrange(p.adjust, pvalue, desc(Count), desc(odds_ratio))

if (file.exists(cluster_file)) {
  clusters <- read.csv(cluster_file, check.names = FALSE) %>%
    select(
      GOID = term_id,
      redundancy_cluster,
      cluster_rep_term_id
    ) %>%
    distinct(GOID, .keep_all = TRUE)
  go_bp <- go_bp %>% left_join(clusters, by = "GOID")
}

behavior_pattern <- regex(
  "behavior|behaviour|locomot|motor|movement|learning|memory|sleep|circadian|rhythm|entrain|sensory|response to stimulus|social|reproductive|mating|maternal|parental|aggression|anxiety|fear|vocal",
  ignore_case = TRUE
)
synaptic_excitability_pattern <- regex(
  "synap|axon|dendrit|action potential|ion channel|calcium channel|neurotrans|neuron migration|neuron projection|neuron remodeling|nervous system",
  ignore_case = TRUE
)
neuro_pattern <- regex(
  paste(
    c(
      "nervous", "neuron", "neuronal", "neurogenesis", "brain", "forebrain", "hindbrain",
      "midbrain", "pallium", "cortex", "cortical", "hippocamp", "axon", "dendrit",
      "synap", "glia", "astrocy", "oligodendro", "myelin", "neurotrans",
      "action potential", "ion channel", "calcium channel", "behavior", "locomot",
      "learning", "memory", "sensory", "motor"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

neuro_only <- go_bp %>%
  filter(str_detect(Description, neuro_pattern)) %>%
  mutate(
    neuro_targeted_q_value = p.adjust(pvalue, method = "BH"),
    significant_neuro_targeted_q05 = neuro_targeted_q_value < 0.05
  ) %>%
  arrange(neuro_targeted_q_value, pvalue, desc(Count))

behavior_neuro <- go_bp %>%
  mutate(
    behavior_related = str_detect(Description, behavior_pattern),
    synaptic_excitability_related = str_detect(Description, synaptic_excitability_pattern),
    behavior_neuro_family = case_when(
      behavior_related & synaptic_excitability_related ~ "Behavior + neural circuit",
      behavior_related ~ "Behavior",
      synaptic_excitability_related ~ "Neural circuit/excitability",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(behavior_neuro_family)) %>%
  mutate(
    behavior_neuro_q_value = p.adjust(pvalue, method = "BH"),
    significant_behavior_neuro_q05 = behavior_neuro_q_value < 0.05
  ) %>%
  arrange(behavior_neuro_q_value, pvalue, desc(Count))

summary_df <- data.frame(
  metric = c(
    "input_genes",
    "original_allcCRE_background_genes",
    "foreground_genes_added_to_background",
    "background_genes_with_bp_annotation",
    "foreground_genes_with_bp_annotation",
    "foreground_genes_unmapped_to_bp",
    "GO_BP_terms_with_foreground_overlap",
    "GO_BP_terms_nominal_p_lt_0.05",
    "GO_BP_terms_q_lt_0.05",
    "neuro_only_terms_tested",
    "neuro_only_terms_targeted_q_lt_0.05",
    "behavior_neuro_terms_tested",
    "behavior_neuro_terms_targeted_q_lt_0.05"
  ),
  value = c(
    length(test_symbols),
    length(background_original),
    length(fg_not_in_background),
    length(bg_genes),
    length(query_genes),
    length(unmapped),
    nrow(go_bp),
    sum(go_bp$pvalue < 0.05),
    sum(go_bp$p.adjust < 0.05),
    nrow(neuro_only),
    sum(neuro_only$significant_neuro_targeted_q05),
    nrow(behavior_neuro),
    sum(behavior_neuro$significant_behavior_neuro_q05)
  )
)

mapping_df <- data.frame(
  gene = test_symbols,
  in_original_allcCRE_background = test_symbols %in% background_original,
  mapped_to_bp = test_symbols %in% query_genes,
  stringsAsFactors = FALSE
)

write.csv(go_bp, file.path(out_dir, "go_bp_allcCRE_background_all_overlap_terms.csv"), row.names = FALSE)
write.csv(go_bp %>% filter(pvalue < 0.05), file.path(out_dir, "go_bp_allcCRE_background_nominal_p05.csv"), row.names = FALSE)
write.csv(go_bp %>% filter(p.adjust < 0.05), file.path(out_dir, "go_bp_allcCRE_background_sig_q05.csv"), row.names = FALSE)
write.csv(head(go_bp, 50), file.path(out_dir, "go_bp_allcCRE_background_top50.csv"), row.names = FALSE)
write.csv(neuro_only, file.path(out_dir, "go_bp_allcCRE_background_neuro_only_targeted_all.csv"), row.names = FALSE)
write.csv(neuro_only %>% filter(significant_neuro_targeted_q05), file.path(out_dir, "go_bp_allcCRE_background_neuro_only_targeted_q05.csv"), row.names = FALSE)
write.csv(behavior_neuro, file.path(out_dir, "go_bp_allcCRE_background_behavior_neuro_targeted_all.csv"), row.names = FALSE)
write.csv(behavior_neuro %>% filter(significant_behavior_neuro_q05), file.path(out_dir, "go_bp_allcCRE_background_behavior_neuro_targeted_q05.csv"), row.names = FALSE)
write.csv(mapping_df, file.path(out_dir, "go_bp_allcCRE_background_input_mapping_summary.csv"), row.names = FALSE)
write.csv(data.frame(gene = fg_not_in_background), file.path(out_dir, "foreground_genes_added_to_allcCRE_background.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "go_bp_allcCRE_background_run_summary.csv"), row.names = FALSE)

cat(capture.output(print(summary_df)), sep = "\n")
cat("\nTop all-cCRE background BP terms:\n")
print(head(go_bp[, c("GOID", "Description", "Count", "pvalue", "p.adjust", "geneID")], 20), row.names = FALSE)
cat("\nTop neuro-only targeted BP terms:\n")
print(head(neuro_only[, c("GOID", "Description", "Count", "pvalue", "neuro_targeted_q_value", "geneID")], 20), row.names = FALSE)
cat("\nTop behavior/neural-circuit targeted BP terms:\n")
print(head(behavior_neuro[, c("GOID", "Description", "behavior_neuro_family", "Count", "pvalue", "behavior_neuro_q_value", "geneID")], 20), row.names = FALSE)
if (length(unmapped) > 0) {
  cat("\nUnmapped from all-cCRE BP annotations:\n")
  cat(paste(unmapped, collapse = ", "), "\n")
}
