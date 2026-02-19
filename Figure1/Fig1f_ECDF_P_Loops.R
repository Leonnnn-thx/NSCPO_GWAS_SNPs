library(ggplot2)
library(dplyr)
library(tidyr)

# ========== 1. read data ==========
Data1 <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results/1.Loop_annotation/Table1.HEPM_MicroC_annotatedLoops.csv")
Data2 <- read.csv("path/to/MEPM_data.csv")  
# ========== 2. gene list ==========
extract_genes <- function(data, gene_col) {
  data %>%
    filter(.data[[gene_col]] != "") %>%
    pull(.data[[gene_col]]) %>%
    strsplit(split = ",") %>%
    unlist() %>%
    trimws()
}

# ========== 3. gene frequency  ==========
count_gene_frequency <- function(loops_data, loop_type_name) {
  P_L_genes <- extract_genes(loops_data, "P_L_genes")
  P_R_genes <- extract_genes(loops_data, "P_R_genes")
  all_genes <- c(P_L_genes, P_R_genes)
  gene_freq <- as.data.frame(table(all_genes))
  colnames(gene_freq) <- c("Gene", "Frequency")
  gene_freq$Loop_Type <- loop_type_name
  
  return(gene_freq)
}

# ========== 4. HEPM data ==========

## 4.1 only E-P loops
HEPM_EP_loops <- Data1 %>%
  filter(interaction_type == "Enhancer-Promoter")

HEPM_EP_freq <- count_gene_frequency(HEPM_EP_loops, "E-P only")
HEPM_EP_freq$Group <- "HEPM"

## 4.2 all promoter relavent loops (E-P + P-P + P-U)
HEPM_all_promoter_loops <- Data1 %>%
  filter(interaction_type %in% c("Enhancer-Promoter", 
                                 "Promoter-Promoter", 
                                 "Promoter-Unannotated"))

HEPM_all_freq <- count_gene_frequency(HEPM_all_promoter_loops, "All loops")
HEPM_all_freq$Group <- "HEPM"

# ========== 5. MEPM data ==========

## 5.1 only 仅E-P loops
MEPM_EP_loops <- Data2 %>%
  filter(interaction_type == "Enhancer-Promoter")

MEPM_EP_freq <- count_gene_frequency(MEPM_EP_loops, "E-P only")
MEPM_EP_freq$Group <- "MEPM"

## 5.2 all promoter relavent loops (E-P + P-P + P-U)

MEPM_all_promoter_loops <- Data2 %>%
  filter(interaction_type %in% c("Enhancer-Promoter", 
                                 "Promoter-Promoter", 
                                 "Promoter-Unannotated"))

MEPM_all_freq <- count_gene_frequency(MEPM_all_promoter_loops, "All loops")
MEPM_all_freq$Group <- "MEPM"

# ========== 6. merge data ==========
combined_data <- rbind(
  HEPM_EP_freq,
  HEPM_all_freq,
  MEPM_EP_freq,
  MEPM_all_freq
)

# 
combined_data$Group_Type <- paste(combined_data$Group, combined_data$Loop_Type, sep = " - ")

# ========== 7. statistical summary ==========
summary_stats <- combined_data %>%
  group_by(Group, Loop_Type) %>%
  summarise(
    Mean = mean(Frequency),
    Median = median(Frequency),
    SD = sd(Frequency),
    Min = min(Frequency),
    Max = max(Frequency),
    N_genes = n(),
    .groups = "drop"
  )

print(summary_stats)

# ========== 8. plot ==========

colors <- c(
  "HEPM - E-P only" = "#E8C9C9",
  "HEPM - All loops" = "#D4A5A5",
  "MEPM - E-P only" = "#B8D4E4",
  "MEPM - All loops" = "#8FB9D4"
)

linetypes <- c(
  "HEPM - E-P only" = "solid",
  "HEPM - All loops" = "dashed",
  "MEPM - E-P only" = "solid",
  "MEPM - All loops" = "dashed"
)

# plot
p_combined <- ggplot(combined_data, aes(x = Frequency, 
                                        color = Group_Type,
                                        linetype = Group_Type)) +
  stat_ecdf(size = 1.2) +

  geom_hline(yintercept = c(0.25, 0.50, 0.75, 0.90), 
             linetype = "dotted", 
             color = "gray50", 
             alpha = 0.4,
             linewidth = 0.3) +
  
  scale_color_manual(values = colors,
                     name = "",
                     labels = c("HEPM E-P loops",
                                "HEPM All promoter loops",
                                "MEPM E-P loops",
                                "MEPM All promoter loops")) +
  
  scale_linetype_manual(values = linetypes,
                        name = "",
                        labels = c("HEPM E-P loops",
                                   "HEPM All promoter loops",
                                   "MEPM E-P loops",
                                   "MEPM All promoter loops")) +
  
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1.00)) +
  
  labs(title = "Cumulative Distribution of Loops per Gene",
       subtitle = "Comparison between E-P loops only vs All promoter-involved loops",
       x = "Number of Loops per Gene",
       y = "Cumulative Probability") +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm")
  )

print(p_combined)


ggsave("cumulative_distribution_EP_vs_AllLoops.pdf", 
       p_combined, 
       width = 12, 
       height = 7)



# color by tissues
colors <- c("HEPM" = "#87CEEB", "MEPM" = "#EC8C12")

p_facet <- ggplot(combined_data, aes(x = Frequency, color = Group)) +
  stat_ecdf(size = 1.2) +

  facet_wrap(~Loop_Type, ncol = 2, 
             labeller = labeller(Loop_Type = c(
               "E-P only" = "E-P Loops Only",
               "All loops" = "All Promoter-involved Loops"
             ))) +

  geom_hline(yintercept = c(0.25, 0.50, 0.75, 0.90), 
             linetype = "dotted", 
             color = "gray50", 
             alpha = 0.4,
             linewidth = 0.3) +

  scale_color_manual(values = colors,
                     name = "Tissue") +

  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1.00)) +

  labs(title = "Cumulative Distribution of Loops per Gene",
       subtitle = "Comparison between E-P loops and All promoter-involved loops",
       x = "Number of Loops per Gene",
       y = "Cumulative Probability") +

  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "gray90", color = "gray50"),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(1.5, "lines")  
  )

print(p_facet)

ggsave("cumulative_distribution_faceted_by_LoopType.pdf", 
       p_facet, 
       width = 12, 
       height = 8)

# ========== 10. save results ==========

write.csv(HEPM_EP_freq, "HEPM_EP_loops_gene_frequency.csv", row.names = FALSE)
write.csv(HEPM_all_freq, "HEPM_all_promoter_loops_gene_frequency.csv", row.names = FALSE)
write.csv(MEPM_EP_freq, "MEPM_EP_loops_gene_frequency.csv", row.names = FALSE)
write.csv(MEPM_all_freq, "MEPM_all_promoter_loops_gene_frequency.csv", row.names = FALSE)

write.csv(summary_stats, "loops_frequency_summary_stats.csv", row.names = FALSE)
