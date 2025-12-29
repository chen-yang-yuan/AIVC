library(dplyr)
library(msigdbr)
library(stringr)

here::i_am("data/_utils/hallmark_pathways.R")

# gene panel
gene_panel <- read.csv(here::here("data/_utils/shared_genes.csv"), header = FALSE)
gene_panel <- gene_panel$V1

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
min_overlap_hallmark <- 50
min_overlap_others <- 50

# ==================== hallmark pathways ==================== #
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

hallmark_filtered <- overlap_filter(hallmark_sets, min_overlap = min_overlap_hallmark, gene_panel = gene_panel) %>%
  mutate(source = "H:hallmark")

message("Kept ", nrow(hallmark_filtered), " out of ", nrow(hallmark_sets), " hallmark pathways with overlap ≥ ", min_overlap_hallmark, " genes.")

# ==================== reactome pathways ==================== #
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, gene_symbol)

kw_reactome <- c("TRANSLATION", "RIBOSOME", "INITIATION", "ELONGATION", "TERMINATION",
                 "SRP", "CO-TRANSLATIONAL", "ER", "ENDOPLASMIC_RETICULUM",
                 "RNA", "MRNA", "SPLIC", "SPLICE", "SPLICING", "PROCESSING",
                 "EXPORT", "TRANSPORT", "NUCLEAR", "CYTOPLAS", "RNP",
                 "NMD", "NONSENSE_MEDIATED", "DEADENYL", "DECAPPING", "EXOSOME",
                 "UNFOLDED_PROTEIN", "UPR", "EIF2", "ISR", "ATF4", "XBP1",
                 "HEAT_SHOCK", "CHAPERONE", "PROTEASOME", "UBIQUITIN",
                 "AUTOPHAG", "LYSOSOME", "OXIDATIVE_STRESS", "REACTIVE_OXYGEN")

reactome_focus <- reactome_sets %>%
  filter(str_detect(gs_name, str_c(kw_reactome, collapse = "|"))) %>%
  distinct()

reactome_filtered <- overlap_filter(reactome_focus, min_overlap = min_overlap_others, gene_panel = gene_panel) %>%
  mutate(source = "C2:REACTOME")

message("Kept ", nrow(reactome_filtered), " out of ", nrow(reactome_sets), " reactome pathways with overlap ≥ ", min_overlap_others, " genes.")

# ==================== GO pathways ==================== #
gobp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  select(gs_name, gene_symbol)

kw_gobp <- c("RNA", "MRNA", "RIBONUCLEO", "SPLIC", "SPLICE", "SPLICING",
             "TRANSLATION", "RIBOSOME", "INITIATION", "ELONGATION",
             "RNA LOCALIZATION", "RNA TRANSPORT", "NUCLEAR EXPORT",
             "STRESS", "RESPONSE TO STRESS", "OXIDATIVE", "HYPOXIA",
             "UNFOLDED PROTEIN", "PROTEIN FOLDING", "CHAPERONE",
             "PROTEASOME", "UBIQUITIN", "AUTOPHAG", "LYSOSOME")

gobp_focus <- gobp_sets %>%
  filter(str_detect(gs_name, str_c(kw_gobp, collapse = "|"))) %>%
  distinct()

gobp_filtered <- overlap_filter(gobp_focus, min_overlap = min_overlap_others, gene_panel = gene_panel) %>%
  mutate(source = "C5:GO_BP")

message("Kept ", nrow(gobp_filtered), " out of ", nrow(gobp_sets), " GO pathways with overlap ≥ ", min_overlap_others, " genes.")

# ==================== Combine ==================== #

all_sets <- bind_rows(hallmark_filtered, reactome_filtered, gobp_filtered) %>%
  filter(overlap_n > 0)

message("Kept:")
message("  Hallmark:  ", nrow(hallmark_filtered))
message("  Reactome:  ", nrow(reactome_filtered))
message("  GO BP:     ", nrow(gobp_filtered))
message("  Total:     ", nrow(all_sets))

jaccard <- function(a, b){
  inter <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(0)
  inter / uni
}

dedup_by_jaccard <- function(df, threshold = 0.7){
  df <- df %>% arrange(desc(overlap_n))
  keep <- rep(TRUE, nrow(df))
  
  for (i in seq_len(nrow(df))) {
    if (!keep[i]) next
    ai <- df$overlap_genes[[i]]
    if (length(ai) == 0) { keep[i] <- FALSE; next }
    for (j in (i+1):nrow(df)) {
      if (j > nrow(df)) break
      if (!keep[j]) next
      aj <- df$overlap_genes[[j]]
      if (jaccard(ai, aj) >= threshold) {
        keep[j] <- FALSE
      }
    }
  }
  df[keep, , drop = FALSE]
}

all_sets_dedup <- dedup_by_jaccard(all_sets, threshold = 0.7)
message("After dedup (Jaccard>=0.7 removed): ", nrow(all_sets_dedup))

# write GMT
# gmt_path <- here::here("data/_utils/hallmark_pathways_filtered.gmt")
# con <- file(gmt_path, open = "wt")
# for (i in seq_len(nrow(hallmark_filtered))) {
#   gs <- hallmark_filtered$gs_name[i]
#   genes <- hallmark_filtered$overlap_genes[[i]]
#   line <- paste(c(gs, "msigdbr_filtered", genes), collapse = "\t")
#   writeLines(line, con)
# }
# close(con)

gmt_path <- here::here("data/_utils/all_pathways_filtered.gmt")
con <- file(gmt_path, open = "wt")
for (i in seq_len(nrow(hallmark_filtered))) {
  gs <- hallmark_filtered$gs_name[i]
  genes <- hallmark_filtered$overlap_genes[[i]]
  line <- paste(c(gs, "msigdbr_filtered", genes), collapse = "\t")
  writeLines(line, con)
}
close(con)

