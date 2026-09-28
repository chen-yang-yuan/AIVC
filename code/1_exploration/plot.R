library(tidyverse)
library(ggalluvial)
library(ggnewscale)

here::i_am("code/1_exploration/plot.R")
nuc_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
nuc_colors_light <- c("#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
cyto_colors <- c("#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3", "#6D1A9C", "#15821E", "#3A84E6", "#997273", "#D1D1D1")
dpi <- 500

sample <- "Xenium_5K_BC"

df <- read.csv(here::here(paste0("output/1_exploration/", sample, "/clustering_results.csv")))

df_sankey <- df %>%
  count(subtype_louvain_nuc, subtype_louvain_cyto) %>%
  filter(!is.na(subtype_louvain_nuc), !is.na(subtype_louvain_cyto)) %>%
  mutate(subtype_louvain_nuc = factor(subtype_louvain_nuc, levels = sort(unique(subtype_louvain_nuc))),
         subtype_louvain_cyto = factor(subtype_louvain_cyto, levels = sort(unique(subtype_louvain_cyto))),
         nuc_stratum  = paste0("nuc_", subtype_louvain_nuc),
         cyto_stratum = paste0("cyto_", subtype_louvain_cyto))

n_level_nuc <- min(length(levels(df_sankey$subtype_louvain_nuc)), length(nuc_colors))
n_level_cyto <- min(length(levels(df_sankey$subtype_louvain_cyto)), length(cyto_colors))
subdomain_levels <- levels(df_sankey$subtype_louvain_nuc)[1:n_level_nuc]
layer_levels  <- levels(df_sankey$subtype_louvain_cyto)[1:n_level_cyto]
stratum_colors <- c(setNames(nuc_colors_light[1:n_level_nuc], subdomain_levels), setNames(cyto_colors[1:n_level_cyto], layer_levels))

p <- ggplot(df_sankey, aes(axis1 = nuc_stratum, axis2 = cyto_stratum, y = n)) +
  geom_alluvium(aes(fill = subtype_louvain_nuc), width = 0.08, alpha = 0.9) +
  scale_fill_manual(values = setNames(nuc_colors_light, subdomain_levels)) +
  new_scale_fill() +  # 👈 reset fill scale
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.08, color = "black") +
  scale_fill_manual(values = c(setNames(nuc_colors[1:n_level_nuc], paste0("nuc_", subdomain_levels)), setNames(cyto_colors[1:n_level_cyto], paste0("cyto_", layer_levels)))) +
  scale_x_discrete(limits = c("Nuclear expression", "Cytoplasmic expression"), expand = c(.1, .1)) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none")

ggsave(here::here(paste0("output/1_exploration/", sample, "/clustering_results.jpeg")), width = 8, height = 6, dpi = dpi)
