library(dplyr)
library(msigdbr)

here::i_am("data/_utils/hallmark_pathways.R")

# gene panel
gene_panel <- read.csv(here::here("data/_utils/shared_genes.csv"), header = FALSE)
gene_panel <- gene_panel$V1

# hallmark pathways
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

# filtering
min_overlap <- 50
hallmark_filtered <- hallmark_sets %>%
  group_by(gs_name) %>%
  summarise(overlap_genes = list(sort(intersect(unique(gene_symbol), gene_panel))),
            overlap_n = length(overlap_genes[[1]]),
            .groups = "drop") %>%
  filter(overlap_n >= min_overlap)

message("Kept ", nrow(hallmark_filtered), " hallmark pathways with overlap ≥ ", min_overlap)

# write GMT
gmt_path <- here::here("data/_utils/hallmark_pathways_filtered.gmt")
con <- file(gmt_path, open = "wt")
for (i in seq_len(nrow(hallmark_filtered))) {
  gs <- hallmark_filtered$gs_name[i]
  genes <- hallmark_filtered$overlap_genes[[i]]
  line <- paste(c(gs, "msigdbr_filtered", genes), collapse = "\t")
  writeLines(line, con)
}
close(con)

