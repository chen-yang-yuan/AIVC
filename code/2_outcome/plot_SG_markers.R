library(ggplot2)
library(readxl)
library(tidyverse)

here::i_am("code/2_outcome/plot_SG_markers.R")
dpi <- 500

# Settings
all_data <- c("Xenium_5K_BC", "Xenium_5K_OC", "Xenium_5K_CC", "Xenium_5K_LC", "Xenium_5K_Prostate", "Xenium_5K_Skin")

# All genes
genes <- read.csv(here::here("data/_utils/shared_genes.csv"), header = FALSE)
genes <- genes$V1

# SG markers
sg_markers_df <- read_excel(here::here("data/_utils/SG_markers.xlsx"))

thr <- 0.25
sg_markers_df <- sg_markers_df %>%
  arrange(desc(`Fraction of RNA molecules in SGs`)) %>%
  filter(`Fraction of RNA molecules in SGs` > thr)
sg_marker_genes <- sg_markers_df$gene

overlap_genes <- sg_marker_genes[sg_marker_genes %in% genes]

# In-cytoplasm ratio
for (data in all_data){
  
  df <- read.csv(here::here(paste0("output/", data, "/in_cytoplasm_ratio.csv")))
  df$is_SG_marker <- ifelse(df$gene %in% overlap_genes, "SG markers", "Non SG markers")
  df$is_SG_marker <- factor(df$is_SG_marker, levels = c("SG markers", "Non SG markers"))
  
  test.res <- t.test(df[df$is_SG_marker == "SG markers", ]$in_cytoplasm_ratio, df[df$is_SG_marker == "Non SG markers", ]$in_cytoplasm_ratio)
  print(paste("Test data:", data, "test statistic:", test.res$statistic, "p-value:", test.res$p.value))
  
  p <- ggplot(df, aes(x = is_SG_marker, y = in_cytoplasm_ratio, fill = is_SG_marker)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.5, outlier.size = 0.5) +
    scale_fill_manual(values = c("SG markers" = "#a0ccec", "Non SG markers" = "#f48488")) +
    labs(x = " ", y = "In-cytoplasm ratio", fill = "Gene type") +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))
  ggsave(here::here(paste0("output/", data, "/in_cytoplasm_ratio.jpeg")), p, width = 6, height = 6, dpi = dpi)

}

# In-SG ratio
df <- read.csv(here::here(paste0("output/", data, "/in_SG_ratio.csv")))
df$is_SG_marker <- factor(df$is_SG_marker, levels = c(1, 0), labels = c("SG markers", "Non SG markers"))
df$is_SG_marker <- factor(df$is_SG_marker, levels = c("SG markers", "Non SG markers"))

t.test(df[df$is_SG_marker == "SG markers", ]$in_SG_ratio, df[df$is_SG_marker == "Non SG markers", ]$in_SG_ratio)

p <- ggplot(df, aes(x = is_SG_marker, y = in_SG_ratio, fill = is_SG_marker)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5, outlier.size = 0.5) +
  scale_fill_manual(values = c("SG markers" = "#a0ccec", "Non SG markers" = "#f48488")) +
  labs(x = " ", y = "In-SG ratio", fill = "Gene type") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12))
ggsave(here::here(paste0("output/", data, "/in_SG_ratio.jpeg")), p, width = 6, height = 6, dpi = dpi)

