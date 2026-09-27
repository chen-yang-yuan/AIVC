library(dplyr)
library(msigdbr)
library(stringr)

here::i_am("data/_utils/hallmark_pathways_test.R")

# helper function
overlap_filter <- function(df, min_overlap, gene_panel) {
  df %>%
    group_by(gs_name) %>%
    summarise(
      overlap_genes = list(sort(intersect(unique(gene_symbol), gene_panel))),
      overlap_n = length(overlap_genes[[1]]),
      .groups = "drop"
    ) %>%
    filter(overlap_n >= min_overlap)
}

# filtering criteria
min_overlap_hallmark <- 5
min_overlap_reactome <- 5
min_overlap_gobp <- 10

# main operation
for (data in c("Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC", "Xenium_5K_LC", "Xenium_5K_Prostate", "Xenium_5K_Skin")){
  
  message("==================== Processing ", data, " ====================")
  
  # gene panel
  gene_panel <- read.csv(here::here(paste0("code/test/SG_vs_nuclei/", data, "/sg_gene_detection_fraction.csv")), header = TRUE)
  gene_panel <- gene_panel$gene[1:100]
  
  # ==================== hallmark pathways ==================== #
  hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
    select(gs_name, gene_symbol)
  
  hallmark_filtered <- overlap_filter(hallmark_sets, min_overlap = min_overlap_hallmark, gene_panel = gene_panel) %>%
    mutate(source = "H:hallmark")
  
  message("Kept ", nrow(hallmark_filtered), " out of ", length(unique(hallmark_sets$gs_name)), " hallmark pathways with overlap ≥ ", min_overlap_hallmark, " genes.")
  
  hallmark_filtered_chr <- hallmark_filtered %>%
    mutate(overlap_genes = purrr::map_chr(overlap_genes, ~ paste(unlist(.x), collapse = ", ")),
           gs_name = gs_name %>%
             str_remove("^HALLMARK_") %>%
             str_to_lower() %>%
             str_replace_all("_", " ") %>%
             str_to_title()) %>%
    arrange(desc(overlap_n))
  
  write.csv(hallmark_filtered_chr, here::here(paste0("data/_utils/test_filtered_pathways_hallmark_", data, ".csv")), row.names = FALSE)
  
  # ==================== reactome pathways ==================== #
  reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
    select(gs_name, gene_symbol)
  
  reactome_filtered <- overlap_filter(reactome_sets, min_overlap = min_overlap_reactome, gene_panel = gene_panel) %>%
    mutate(source = "C2:REACTOME")
  
  message("Kept ", nrow(reactome_filtered), " out of ", length(unique(reactome_sets$gs_name)), " reactome pathways with overlap ≥ ", min_overlap_reactome, " genes.")
  
  reactome_filtered_chr <- reactome_filtered %>%
    mutate(overlap_genes = purrr::map_chr(overlap_genes, ~ paste(unlist(.x), collapse = ", ")),
           gs_name = gs_name %>%
             str_remove("^REACTOME_") %>%
             str_to_lower() %>%
             str_replace_all("_", " ") %>%
             str_to_title()) %>%
    arrange(desc(overlap_n))
  
  write.csv(reactome_filtered_chr, here::here(paste0("data/_utils/test_filtered_pathways_reactome_", data, ".csv")), row.names = FALSE)
  
  # ==================== GO pathways ==================== #
  gobp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
    select(gs_name, gene_symbol)
  
  gobp_filtered <- overlap_filter(gobp_sets, min_overlap = min_overlap_gobp, gene_panel = gene_panel) %>%
    mutate(source = "C5:GO_BP")
  
  message("Kept ", nrow(gobp_filtered), " out of ", length(unique(gobp_sets$gs_name)), " GO pathways with overlap ≥ ", min_overlap_gobp, " genes.")
  
  gobp_filtered_chr <- gobp_filtered %>%
    mutate(overlap_genes = purrr::map_chr(overlap_genes, ~ paste(unlist(.x), collapse = ", ")),
           gs_name = gs_name %>%
             str_remove("^GOBP_") %>%
             str_to_lower() %>%
             str_replace_all("_", " ") %>%
             str_to_title()) %>%
    arrange(desc(overlap_n))
  
  write.csv(gobp_filtered_chr, here::here(paste0("data/_utils/test_filtered_pathways_gobp_", data, ".csv")), row.names = FALSE)
  
}










