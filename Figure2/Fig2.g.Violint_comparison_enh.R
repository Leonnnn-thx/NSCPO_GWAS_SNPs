#------------ Plot 3. Enhancer number of Genes ---------------------------------
setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/3.Results/2.Enhancer_Hierary_Network/Fig3G_enhancer_comparison")
load("Fig3g_violint_comparisonplot.Rdata")
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(patchwork)
library(dplyr) 

# load data
Gene_regulatory_layers_HEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/HEPM/HEPM_gene_regulatory_layers.csv")
Gene_regulatory_layers_MEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/MEPM/MEPM_gene_regulatory_layers.csv")

Combined_CPGene <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/CP_combined_with_mouse_orthologs.csv", sep=",", header = T)
Combined_CPGene_HEPM <- Combined_CPGene$CP_Gene %>% unique()
Combined_CPGene_MEPM <- Combined_CPGene$mm_gene_symbol %>% 
  .[!is.na(.)] %>% 
  unique()

# ========== 1. Calculate the number and proportion of Tier 2. ==========

# HEPM - All Genes
HEPM_analysis <- Gene_regulatory_layers_HEPM %>%
  mutate(
    Additional_Enhancers = Total_Enhancers - Tier1_Count,
    Group = "HEPM_All"
  )

# MEPM - All Genes
MEPM_analysis <- Gene_regulatory_layers_MEPM %>%
  mutate(
    Additional_Enhancers = Total_Enhancers - Tier1_Count,
    Group = "MEPM_All"
  )

# HEPM - CP Genes
HEPM_CP_analysis <- HEPM_analysis %>%
  filter(Gene %in% Combined_CPGene_HEPM) %>%
  mutate(Group = "HEPM_CP")

# MEPM - CP Genes
MEPM_CP_analysis <- MEPM_analysis %>%
  filter(Gene %in% Combined_CPGene_MEPM) %>%
  mutate(Group = "MEPM_CP")

# stacked bar plot
stacked_data <- bind_rows(
  HEPM_analysis %>% dplyr::select(Gene, Tier1_Count, Additional_Enhancers, Group),
  MEPM_analysis %>% dplyr::select(Gene, Tier1_Count, Additional_Enhancers, Group),
  HEPM_CP_analysis %>% dplyr::select(Gene, Tier1_Count, Additional_Enhancers, Group),
  MEPM_CP_analysis %>% dplyr::select(Gene, Tier1_Count, Additional_Enhancers, Group)
) %>%
  group_by(Group) %>%
  summarise(
    Mean_Tier1 = mean(Tier1_Count),
    Mean_Additional = mean(Additional_Enhancers),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(Mean_Tier1, Mean_Additional),
               names_to = "Tier",
               values_to = "Count") %>%
  mutate(
    Category = ifelse(grepl("CP", Group), "CP Genes", "All Genes"),
    Species = ifelse(grepl("HEPM", Group), "Human", "Mouse"),
    Tier = case_when(
      Tier == "Mean_Tier1" ~ "Tier1",
      Tier == "Mean_Additional" ~ "Additional (Network)",
      TRUE ~ Tier
    )
  )

# check
print("Data structure:")
print(head(stacked_data))
print(unique(stacked_data$Tier))

# plot

p1 <- ggplot(stacked_data, aes(x = Species, y = Count, fill = Tier)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = round(Count, 1)), 
            position = position_stack(vjust = 0.5), 
            size = 3.5, color = "white", fontface = "bold") +
  facet_wrap(~Category) +
  scale_fill_manual(
    values = c("Tier1" = "#999999", "Additional (Network)" = "#E64B35"),
    breaks = c("Tier1", "Additional (Network)")
  ) +
  labs(title = "Network Found More Enhancers: Composition of Total",
       subtitle = "Tier1 alone vs Additional enhancers found by network",
       x = "", y = "Mean Number of Enhancers per Gene",
       fill = "Enhancer Type") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 11),
        plot.title = element_text(face = "bold", size = 14))

print(p1)
ggsave("Network_Enhancer_Composition.pdf", p1, width = 10, height = 6)

# ========== other version ==========

comparison_data <- bind_rows(
  HEPM_analysis %>% dplyr::select(Gene, Tier1_Count, Total_Enhancers, Group),
  MEPM_analysis %>% dplyr::select(Gene, Tier1_Count, Total_Enhancers, Group),
  HEPM_CP_analysis %>% dplyr::select(Gene, Tier1_Count, Total_Enhancers, Group),
  MEPM_CP_analysis %>% dplyr::select(Gene, Tier1_Count, Total_Enhancers, Group)
) %>%
  group_by(Group) %>%
  summarise(
    Mean_Tier1 = mean(Tier1_Count),
    Mean_Total = mean(Total_Enhancers),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(Mean_Tier1, Mean_Total),
               names_to = "Type",
               values_to = "Count") %>%
  mutate(
    Category = ifelse(grepl("CP", Group), "CP Genes", "All Genes"),
    Species = ifelse(grepl("HEPM", Group), "Human", "Mouse"),
    Type = case_when(
      Type == "Mean_Tier1" ~ "Tier1 Only",
      Type == "Mean_Total" ~ "Total (Network)",
      TRUE ~ Type
    )
  )

p2 <- ggplot(comparison_data, aes(x = Species, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = round(Count, 1)), 
            position = position_dodge(0.8), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  facet_wrap(~Category) +
  scale_fill_manual(
    values = c("Tier1 Only" = "#999999", "Total (Network)" = "#E64B35")
  ) +
  labs(title = "Network Identified More Regulated Enhancers",
       subtitle = "Comparison: Tier1 alone vs Total enhancers found by network",
       x = "", y = "Mean Number of Enhancers per Gene",
       fill = "") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 14))

print(p2)
ggsave("Network_vs_Tier1_Comparison.pdf", p2, width = 10, height = 6)

# ------ violint plot

# HEPM - All Genes
HEPM_analysis <- Gene_regulatory_layers_HEPM %>%
  mutate(
    Additional_Enhancers = Total_Enhancers - Tier1_Count,
    Fold_Increase = Total_Enhancers / Tier1_Count,
    Percent_Increase = (Total_Enhancers - Tier1_Count) / Tier1_Count * 100,
    Category = "All Genes",
    Species = "Human"
  )

# MEPM - All Genes
MEPM_analysis <- Gene_regulatory_layers_MEPM %>%
  mutate(
    Additional_Enhancers = Total_Enhancers - Tier1_Count,
    Fold_Increase = Total_Enhancers / Tier1_Count,
    Percent_Increase = (Total_Enhancers - Tier1_Count) / Tier1_Count * 100,
    Category = "All Genes",
    Species = "Mouse"
  )

# HEPM - CP Genes
HEPM_CP_analysis <- Gene_regulatory_layers_HEPM %>%
  filter(Gene %in% Combined_CPGene_HEPM) %>%
  mutate(
    Additional_Enhancers = Total_Enhancers - Tier1_Count,
    Fold_Increase = Total_Enhancers / Tier1_Count,
    Percent_Increase = (Total_Enhancers - Tier1_Count) / Tier1_Count * 100,
    Category = "CP Genes",
    Species = "Human"
  )

# MEPM - CP Genes
MEPM_CP_analysis <- Gene_regulatory_layers_MEPM %>%
  filter(Gene %in% Combined_CPGene_MEPM) %>%
  mutate(
    Additional_Enhancers = Total_Enhancers - Tier1_Count,
    Fold_Increase = Total_Enhancers / Tier1_Count,
    Percent_Increase = (Total_Enhancers - Tier1_Count) / Tier1_Count * 100,
    Category = "CP Genes",
    Species = "Mouse"
  )

# combined data
all_data <- bind_rows(
  HEPM_analysis,
  MEPM_analysis,
  HEPM_CP_analysis,
  MEPM_CP_analysis
)

# ========== violin plot：Tier1 vs Total ==========

violin_data <- all_data %>%
  dplyr::select(Gene, Category, Species, Tier1_Count, Total_Enhancers) %>%
  pivot_longer(cols = c(Tier1_Count, Total_Enhancers),
               names_to = "Type",
               values_to = "Count") %>%
  mutate(Type = case_when(
    Type == "Tier1_Count" ~ "Tier1 Only",
    Type == "Total_Enhancers" ~ "Total enhancers (Network)",
    TRUE ~ Type
  ))
p5 <- ggplot(violin_data, aes(x = Species, y = Count, fill = Type)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.10, position = position_dodge(0.9), 
               outlier.size = 0.5, alpha = 0.8) +
  stat_compare_means(aes(group = Type), method = "wilcox.test", 
                     label = "p.signif", label.y.npc = 0.95) +
  facet_wrap(~Category, scales = "free_y") +
  scale_fill_manual(values = c("Tier1 Only" = "#999999", "Total enhancers (Network)" = "#C9A063")) +
  labs(title = "Network Found Significantly More Enhancers",
       subtitle = "Statistical comparison between Tier1 and Total",
       x = "", y = "Number of Enhancers per Gene",
       fill = "") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 14))

print(p5)
ggsave("Violin_With_Statistics.pdf", p5, width = 12, height = 9)


# ====================
descriptive_stats <- violin_data %>%
  group_by(Category, Species, Type) %>%
  summarise(
    N = n(),
    Mean = mean(Count, na.rm = TRUE),
    Median = median(Count, na.rm = TRUE),
    SD = sd(Count, na.rm = TRUE),
    SE = SD / sqrt(N),
    Min = min(Count, na.rm = TRUE),
    Max = max(Count, na.rm = TRUE),
    Q1 = quantile(Count, 0.25, na.rm = TRUE),
    Q3 = quantile(Count, 0.75, na.rm = TRUE),
    IQR = IQR(Count, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    CI_lower = Mean - 1.96 * SE,
    CI_upper = Mean + 1.96 * SE
  )

print(descriptive_stats)
write.csv(descriptive_stats, "Violin_Descriptive_Statistics.csv", row.names = FALSE)


# ========== summary enhancer linking to CP Genes  ==========

# Human CP Genes
cat("\n========== Human CP Genes ==========\n")
cat("Total CP Genes:", nrow(HEPM_CP_analysis), "\n")
cat("Total Tier1 Enhancers:", sum(HEPM_CP_analysis$Tier1_Count, na.rm = TRUE), "\n") # 1765
cat("Total Network Enhancers:", sum(HEPM_CP_analysis$Total_Enhancers, na.rm = TRUE), "\n") # 3381

# Mouse CP Genes
cat("\n========== Mouse CP Genes ==========\n")
cat("Total CP Genes:", nrow(MEPM_CP_analysis), "\n")
cat("Total Tier1 Enhancers:", sum(MEPM_CP_analysis$Tier1_Count, na.rm = TRUE), "\n") # 2076
cat("Total Network Enhancers:", sum(MEPM_CP_analysis$Total_Enhancers, na.rm = TRUE), "\n") # 3717


save(list = ls(), file = "Fig3g_violint_comparisonplot.Rdata")
