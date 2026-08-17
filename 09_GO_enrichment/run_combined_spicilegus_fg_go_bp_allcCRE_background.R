#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(GO.db)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

root <- Sys.getenv("MSPIC_PROJECT_ROOT", ".")
background_path <- file.path(Sys.getenv("MSPIC_INPUT_DIR", "data/inputs"), "go_bg_allcCRE_genes.txt")
out_dir <- file.path(root, "phyloP_MSA_QC", "go_enrichment_spicilegus_combined_fg_allcCRE_background")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes_raw <- c(
  "Sipa1l3", "Slc25a14", "Crebrf", "Ikbkg", "Alkbh6", "Lncpint", "Ing1",
  "Parn", "Plekhh3", "Raver1", "Sdha", "Srsf4", "Zfp574", "A330040F15Rik",
  "Stbd1", "Atp5mc3", "Enc1", "Ctnnb1", "Grin2d", "Npas4", "Nrxn2",
  "2610507B11Rik", "Abcg4", "Abhd11os", "Adam23", "Adgrd2-ps", "Adrb1",
  "Ajuba", "Ank", "Ankrd11", "Avl9", "B4galnt1", "Bbx", "Bcl11a",
  "Bcl9", "Bcl9", "Bmpr2", "Cd1d2", "Clybl", "Cnn3", "Coro2b",
  "Dync1li1", "E130114P18Rik", "Emx2os", "Gm13838", "Gm50348",
  "Gm50394", "Gna15", "Gpr173", "Gse1", "Gtf2ird1", "Hdac7", "Hddc3",
  "Hoxd1", "Hsd17b2", "Igf2os", "Jade1", "Jade2", "Kcnh2", "Kdm7a",
  "Lims2", "Lrrc42", "Lrrc42", "Midn", "Mir19a", "Mir21a", "Mir7045",
  "Mmp15", "Nadsyn1", "Nars", "Nectin3", "Nfam1", "Nfat5", "Nfe2l1",
  "Npas3", "Nsd1", "Ormdl3", "Palmd", "Pcdhb12", "Phtf2", "Pon2",
  "Ptbp3", "R3hdm4", "Rab11fip3", "Rai1", "Rapgef4", "Rassf5", "Rbpj",
  "Rian", "Sema6b", "Sh3bgrl", "Slc10a4", "Slc49a4", "Smad3", "Sox5",
  "Srgap1", "Tacr1", "Tet3", "Tmem238l", "Tsr1", "Ube2o", "Vegfc",
  "Vezf1", "Wnk1", "Wnk2", "Zbtb16", "Zfp512b", "Zfp609",
  "Vxn", "Rims1", "Stk17b", "Dis3l2", "Gm29461", "Gm23782", "Hhat",
  "Gm48132", "Pttg1ip", "Mtmr3", "Ugp2", "Cfap36", "Nkiras2",
  "Gm12589", "Ywhaq", "Twistnb", "Ifrd1", "Foxg1", "Gpr137b",
  "Phactr1", "Gm36101", "Iqgap2", "Fam124a", "Gm47010", "Egr3",
  "C1qtnf3", "Gm48956", "Mroh5", "Cenpm", "Gm49499", "Tprg", "Btg3",
  "Gm28118", "Gm46610", "Rpl36-ps4", "Gm5839", "Dtwd2", "Dntt",
  "Gm14488", "Gm13485", "Pamr1", "Gm43008", "Mab21l1", "Syt11",
  "Strip1", "Elovl6", "Gm43619", "Gm11795", "Ptpn3", "Astn2",
  "Gm12696", "Ece1", "Rell1", "Tesc", "4930413E15Rik", "Hpd", "Psph",
  "Gng12", "Gm25961", "Chl1", "Chl1", "Chl1", "Ybx3", "Snx19",
  "Gm37987", "Fundc1", "Cdk16", "Pgrmc1", "Maged1", "Rlim", "Taf9b",
  "Pof1b", "Prps1", "Tsr2", "Fam120c", "Cdkl5", "Tmsb4x", "Gm25096",
  "Fhl2", "Cd34", "Gm6390", "Sbno2", "Btbd2", "Mir6916", "Myocd",
  "Slc2a4", "Axin2", "Fasn", "Pgf", "Gpr137b-ps", "Gm33489", "Gmpr",
  "Zfp503", "Lrmda", "Zcchc24", "Gja3", "Gm24981", "Bmp1", "Itm2b",
  "Fgf14", "Gm38563", "Eif4b", "Eif2b5", "Bcl6", "Cpox", "D16Ertd472e",
  "Prss30", "Tmem8", "Tbcc", "Kif6", "Wdr43", "Gpr151", "Gstp3",
  "Gcnt1", "Borcs7", "Itga8", "Cacnb2", "Mbd5", "Plcb1", "Bcas1",
  "Polr3k", "Armc1", "Armc1", "Gm11825", "Galt", "Dbf4", "Gm10220",
  "Pgm2", "Fras1", "Wbscr25", "Gm32222", "Aass", "Fbxl14", "Gm44124",
  "Slco1a4", "Slco1a4", "Hpn", "Nudt19", "Vrk3", "Whamm", "Maz",
  "Ppp1r3b", "Adcy7", "Gm24602", "Exoc8", "Cep126", "Robo3", "Lipc",
  "Gm20477", "Dock3", "Trim71", "Porcn", "Bcor", "Bcor", "Maob",
  "Mcts1", "Arhgef6", "F9", "Pabpc5", "Prps1", "Gpr173", "Sh3kbp1",
  "Gpm6b", "Jazf1", "Lif", "Gm20644", "Plk2"
)

foreground_rows <- tibble(gene = genes_raw) %>%
  mutate(input_order = row_number(), duplicated_in_input = duplicated(gene) | duplicated(gene, fromLast = TRUE))
foreground <- foreground_rows %>%
  distinct(gene) %>%
  arrange(gene) %>%
  pull(gene)

background_original <- read_lines(background_path) %>%
  str_trim() %>%
  .[. != ""] %>%
  unique()
fg_not_in_background <- setdiff(foreground, background_original)
background <- sort(unique(c(background_original, fg_not_in_background)))

write_csv(foreground_rows, file.path(out_dir, "foreground_input_rows_with_duplicates.csv"))
write_lines(foreground, file.path(out_dir, "foreground_unique_genes.txt"))
write_csv(tibble(gene = fg_not_in_background), file.path(out_dir, "foreground_genes_added_to_background.csv"))
write_lines(background, file.path(out_dir, "background_allcCRE_genes_plus_foreground.txt"))

annot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = background,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GENENAME", "ENTREZID", "GOALL", "ONTOLOGYALL")
)

gene_desc <- annot %>%
  distinct(SYMBOL, GENENAME, ENTREZID) %>%
  filter(!is.na(SYMBOL))

foreground_annot <- tibble(gene = foreground) %>%
  left_join(gene_desc, by = c("gene" = "SYMBOL")) %>%
  mutate(in_original_background = gene %in% background_original)
write_csv(foreground_annot, file.path(out_dir, "foreground_unique_gene_annotations.csv"))

go_annot <- annot %>%
  distinct(SYMBOL, GOALL, ONTOLOGYALL) %>%
  filter(!is.na(SYMBOL), !is.na(GOALL), ONTOLOGYALL == "BP")

bg_genes <- unique(go_annot$SYMBOL)
query_genes <- intersect(foreground, bg_genes)

go_bp <- bind_rows(lapply(sort(unique(go_annot$GOALL)), function(go_id) {
  term_genes <- unique(go_annot$SYMBOL[go_annot$GOALL == go_id])
  overlap <- intersect(query_genes, term_genes)
  a <- length(overlap)
  if (a == 0) return(NULL)
  b <- length(query_genes) - a
  c <- length(term_genes) - a
  d <- length(bg_genes) - a - b - c
  p <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
  term_name <- tryCatch(as.character(Term(GOTERM[[go_id]])), error = function(e) NA_character_)
  tibble(
    go_id = go_id,
    term_name = term_name,
    ontology = "BP",
    overlap_n = a,
    query_n = length(query_genes),
    bg_term_n = length(term_genes),
    bg_n = length(bg_genes),
    overlap_genes = paste(sort(overlap), collapse = ";"),
    p_value = p
  )
})) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(q_value, p_value, desc(overlap_n))

go_sig <- go_bp %>% filter(q_value < 0.05)
go_bp_top <- go_bp %>% slice_head(n = 25)

write_csv(go_bp, file.path(out_dir, "combined_fg_GO_BP_fisher_all_terms.csv"))
write_csv(go_sig, file.path(out_dir, "combined_fg_GO_BP_fisher_q05.csv"))
write_csv(go_bp_top, file.path(out_dir, "combined_fg_GO_BP_top25.csv"))

behavior_pattern <- regex("behavior|behaviour|locomot|motor|movement|learning|memory|sleep|circadian|rhythm|entrain|sensory|response to stimulus|social|reproductive|mating|maternal|parental|aggression|anxiety|fear|vocal", ignore_case = TRUE)
synaptic_excitability_pattern <- regex("synap|axon|dendrit|action potential|ion channel|calcium channel|neurotrans|neuron migration|neuron projection|neuron remodeling|nervous system", ignore_case = TRUE)

behavior_neuro <- go_bp %>%
  mutate(
    behavior_related = str_detect(term_name, behavior_pattern),
    synaptic_excitability_related = str_detect(term_name, synaptic_excitability_pattern),
    behavior_neuro_family = case_when(
      behavior_related & synaptic_excitability_related ~ "Behavior + neural circuit",
      behavior_related ~ "Behavior",
      synaptic_excitability_related ~ "Neural circuit/excitability",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(behavior_neuro_family)) %>%
  mutate(
    behavior_neuro_q_value = p.adjust(p_value, method = "BH"),
    significant_behavior_neuro_q05 = behavior_neuro_q_value < 0.05
  ) %>%
  arrange(behavior_neuro_q_value, p_value, desc(overlap_n))

behavior_neuro_sig <- behavior_neuro %>% filter(significant_behavior_neuro_q05)
behavior_neuro_plot <- bind_rows(
    behavior_neuro_sig,
    behavior_neuro %>% filter(go_id == "GO:0007399")
  ) %>%
  distinct(go_id, .keep_all = TRUE) %>%
  arrange(behavior_neuro_q_value, p_value) %>%
  mutate(
    term_label = recode(
      term_name,
      "positive regulation of high voltage-gated calcium channel activity" =
        "Positive regulation of\nhigh voltage-gated calcium channel activity",
      .default = str_wrap(term_name, width = 34)
    ),
    term_label = factor(term_label, levels = rev(term_label)),
    neg_log10_q = -log10(behavior_neuro_q_value),
    gene_count_label = paste0("n=", overlap_n),
    star_label = if_else(significant_behavior_neuro_q05, "*", ""),
    star_x = neg_log10_q + 0.035,
    count_x = neg_log10_q + 0.18
  )

write_csv(behavior_neuro, file.path(out_dir, "combined_fg_behavior_relevant_neuro_targeted_all.csv"))
write_csv(behavior_neuro_sig, file.path(out_dir, "combined_fg_behavior_relevant_neuro_targeted_q05.csv"))
write_csv(behavior_neuro_plot, file.path(out_dir, "combined_fg_behavior_relevant_neuro_plot_data.csv"))

summary_df <- tibble(
  metric = c(
    "input_rows",
    "unique_foreground_genes",
    "duplicated_input_rows",
    "background_original_genes",
    "foreground_genes_added_to_background",
    "background_genes_used",
    "foreground_genes_with_bp_go_annotation",
    "GO_BP_terms_tested",
    "GO_BP_terms_q_lt_0.05",
    "behavior_neuro_terms_tested",
    "behavior_neuro_terms_q_lt_0.05"
  ),
  value = c(
    length(genes_raw),
    length(foreground),
    sum(foreground_rows$duplicated_in_input),
    length(background_original),
    length(fg_not_in_background),
    length(background),
    length(query_genes),
    nrow(go_bp),
    nrow(go_sig),
    nrow(behavior_neuro),
    nrow(behavior_neuro_sig)
  )
)
write_csv(summary_df, file.path(out_dir, "combined_fg_GO_BP_summary.csv"))

if (nrow(behavior_neuro_plot) > 0) {
  p <- ggplot(behavior_neuro_plot, aes(x = term_label, y = neg_log10_q)) +
    geom_col(fill = "#002F6C", width = 0.68) +
    geom_text(aes(y = star_x, label = star_label), hjust = 0.5, vjust = 0.5, size = 6.6, color = "black") +
    geom_text(aes(y = count_x, label = gene_count_label), hjust = -0.05, vjust = 0.5, size = 4.0, color = "black") +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(
      x = "Neurological GO:BP",
      y = expression(-log[10]("BH q-value"))
    ) +
    theme_bw(base_size = 14) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 12, lineheight = 0.95),
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.title.y = element_text(margin = margin(r = 12)),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.45),
      plot.margin = margin(8, 56, 8, 8)
    )

  ggsave(file.path(out_dir, "combined_fg_behavior_neuro_GO_allcCRE.png"), p, width = 8.2, height = 4.8, dpi = 600)
  ggsave(file.path(out_dir, "combined_fg_behavior_neuro_GO_allcCRE.pdf"), p, width = 8.2, height = 4.8)
  ggsave(file.path(out_dir, "combined_fg_behavior_neuro_GO_allcCRE.tiff"), p, width = 8.2, height = 4.8, dpi = 600, compression = "lzw")
}

cat(capture.output(print(summary_df)), sep = "\n")
if (nrow(behavior_neuro_sig) > 0) {
  cat("\nTop behavior-neuro targeted terms:\n")
  print(head(behavior_neuro_sig, 20))
}
