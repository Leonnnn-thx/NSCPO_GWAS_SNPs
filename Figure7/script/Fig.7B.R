library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scales)
setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/24.TFbindingSNP/Results/HEPM/2.TOBIAS_results/2.Footprintscore/1.P70_TF_region")
data <- read.table("CPenh3381_Top30_vs_Bot70_comparison.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Combined_CPGene <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/CP_combined_with_mouse_orthologs.csv", sep=",", header = T)
Combined_CPGene_HEPM <- Combined_CPGene$CP_Gene %>% unique()
# ----- bound sites comparison ----------------------------

data_clean <- data %>%
  filter(!is.na(bound_Top30) & 
           !is.na(bound_Bot70) & 
           !is.na(bound_diff))
data_clean <- data_clean %>%
  mutate(is_cp_gene = name %in% Combined_CPGene_HEPM)
n_cp_genes <- sum(data_clean$is_cp_gene)
cat(sprintf("\nFound %d TFs in Combined_CPGene_HEPM (out of %d total TFs)\n", 
            n_cp_genes, nrow(data_clean)))  # Found 46 TFs in Combined_CPGene_HEPM (out of 602 total TFs)
cp_gene_data <- data_clean %>%
  filter(is_cp_gene) %>%
  arrange(desc(abs(bound_diff)))
# Bound sites range
min_bound <- min(c(data_clean$bound_Top30, data_clean$bound_Bot70))
max_bound <- max(c(data_clean$bound_Top30, data_clean$bound_Bot70))
bound_range <- max_bound - min_bound
min_bound <- max(0, min_bound - bound_range * 0.05) 
max_bound <- max_bound + bound_range * 0.05


# Bound difference range
min_diff <- min(data_clean$bound_diff)
max_diff <- max(data_clean$bound_diff)
abs_max <- max(abs(min_diff), abs(max_diff))

# threshold = 0
diff_threshold <- 0 

data_clean <- data_clean %>%
  mutate(
    significance = case_when(
      bound_diff > diff_threshold ~ "Top30-enriched",
      bound_diff < -diff_threshold ~ "Bot70-enriched",
      TRUE ~ "Not significant"
    )
  )
genes_to_label <- c("CTCF", "OSR2", "TWIST1", "IKZF2", "TFAP2A")
labeled_genes_data <- data_clean %>%
  filter(name %in% genes_to_label)
p2_symmetric <- ggplot(data_clean, aes(x = bound_Bot70, y = bound_Top30, color = bound_diff)) +
  geom_point(alpha = 0.4, size = 5, shape = 16) +
  geom_point(data = cp_gene_data,
             aes(x = bound_Bot70, y = bound_Top30),
             size = 5, alpha = 0.8, shape = 21, 
             stroke = 1, fill = NA, color = "black") +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "gray30", size = 0.8) +
  geom_text_repel(
    data = labeled_genes_data,
    aes(label = name),
    size = 4,
    fontface = "bold",
    box.padding = 1,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 1,
    arrow = arrow(length = unit(0.02, "npc")),
    max.overlaps = 20,
    color = "black",
    bg.color = "white",
    bg.r = 0.15
  ) +
  scale_color_gradientn(
    colors = rev(brewer.pal(11, "RdYlBu")),
    limits = c(-abs_max, abs_max),
    breaks = c(-abs_max, 0, abs_max),
    labels = c(sprintf("%.0f", -abs_max), "0", sprintf("%.0f", abs_max)),
    name = "Bound Diff\n(Top30 - Bot70)",
    guide = guide_colorbar(
      barwidth = 1.5,
      barheight = 15,
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  coord_fixed(ratio = 1,
              xlim = c(min_bound, max_bound),
              ylim = c(min_bound, max_bound)) +

  labs(
    title = "TF Binding Sites: Top30 vs Bot70",
    subtitle = sprintf("Color indicates binding difference (n = %d TFs)", nrow(data_clean)),
    x = "Bound Sites in Bot70",
    y = "Bound Sites in Top30"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 1
  )

print(p2_symmetric)

ggsave("Fig8B.TF_bound_sites_scatter_Top30_vs_Bot70.pdf", 
       plot = p2_symmetric, width = 8, height = 8, dpi = 300)


# ------------------------------------------------------------------------------------------

# -------- Plot by score_log2fc---------------------------------
data_clean <- data %>%
  filter(!is.na(score_Top30) & 
           !is.na(score_Bot70) & 
           !is.na(score_log2fc) &
           !is.infinite(score_log2fc))

min_score <- min(c(data_clean$score_Top30, data_clean$score_Bot70))
max_score <- max(c(data_clean$score_Top30, data_clean$score_Bot70))
score_range <- max_score - min_score
min_score <- min_score - score_range * 0.05
max_score <- max_score + score_range * 0.05

# Log2FC range
min_log2fc <- min(data_clean$score_log2fc)
max_log2fc <- max(data_clean$score_log2fc)
max_abs_log2fc <- max(abs(min_log2fc), abs(max_log2fc))

# label significant tf-motif
log2fc_threshold <- 2 
data_clean <- data_clean %>%
  mutate(
    significance = case_when(
      score_log2fc > log2fc_threshold ~ "Top30-enriched",
      score_log2fc < -log2fc_threshold ~ "Bot70-enriched",
      TRUE ~ "Not significant"
    )
  )

# plot 
p1 <- ggplot(data_clean, aes(x = score_Bot70, y = score_Top30, color = score_log2fc)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "gray30", size = 0.8) +
  scale_color_gradientn(
    colors = rev(brewer.pal(11, "RdYlBu")),
    limits = c(-max_abs_log2fc, max_abs_log2fc),
    breaks = c(-max_abs_log2fc, 0, max_abs_log2fc),
    labels = c(sprintf("%.2f", -max_abs_log2fc), "0", sprintf("%.2f", max_abs_log2fc)),
    name = "Score Log2FC\n(Top30/Bot70)",
    guide = guide_colorbar(
      barwidth = 1.5,
      barheight = 15,
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  coord_fixed(ratio = 1, 
              xlim = c(min_score, max_score), 
              ylim = c(min_score, max_score)) +
  labs(
    title = "TF Binding Scores: Top30 vs Bot70",
    subtitle = sprintf("Color indicates log2 fold change (n = %d TFs)", nrow(data_clean)),
    x = "Score in Bot70",
    y = "Score in Top30"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 1
  )

ggsave("Fig.8B.TF_scores_scatter_Top30_vs_Bot70.pdf", 
       plot = p1, width = 8, height = 7, dpi = 300)

save(list = ls(), file = "Fig.8B.Comparison_top30_bot70_TFmotif.Rdata")
