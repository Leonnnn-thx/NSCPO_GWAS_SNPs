setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/3.Results/2.Enhancer_Hierary_Network/Fig3f_h")
load("Fig3_SEs_network.Rdata")
library(dplyr)

# -- load data ---
Gene_regulatory_layers_HEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/HEPM/Fig3G_enhancer_comparison_violint/HEPM_gene_regulatory_layers.csv")
Gene_regulatory_layers_MEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/MEPM/MEPM_gene_regulatory_layers.csv")

Deep_gene_list_HEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/HEPM/deep_regulation_genes.txt", sep = " ", header = F)
Deep_gene_list_MEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/MEPM/deep_regulation_genes.txt", sep = " ", header = F)

Deep_gene_list_HEPM <- Deep_gene_list_HEPM %>% dplyr::rename(Gene = V1)
Deep_gene_list_MEPM <- Deep_gene_list_MEPM %>% dplyr::rename(Gene = V1)

Shadow_Gene_list_HEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/HEPM/shallow_regulation_genes.txt", sep = " ", header = F)
Shadow_Gene_list_MEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/MEPM/shallow_regulation_genes.txt", sep = " ", header = F)

Shallow_gene_list_HEPM <- Shadow_Gene_list_HEPM %>% dplyr::rename(Gene = V1)
Shallow_gene_list_MEPM <- Shadow_Gene_list_MEPM %>% dplyr::rename(Gene = V1)


Combined_CPGene <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/CP_combined_with_mouse_orthologs.csv", sep=",", header = T)

# ---plot1 : CP gene distribution of shallow and deep group ------
cp_genes <- Combined_CPGene$CP_Gene %>% unique()
deep_genes_HEPM <- Deep_gene_list_HEPM$Gene
deep_genes_MEPM <- Deep_gene_list_MEPM$Gene

shallow_genes_HEPM <- Shallow_gene_list_HEPM$Gene
shallow_genes_MEPM <- Shallow_gene_list_MEPM$Gene

# - stack bar plot 
all_genes_data <- data.frame(
  Group = rep(c("HEPM", "MEPM"), each = 2),
  Type = rep(c("Shallow", "Deep"), 2),
  Count = c(
    length(shallow_genes_HEPM),
    length(deep_genes_HEPM),
    length(shallow_genes_MEPM),
    length(deep_genes_MEPM)
  ),
  Category = "All genes"
)

# distribution of CP genes in HEPM 
cp_deep_HEPM <- sum(cp_genes %in% deep_genes_HEPM)
cp_shallow_HEPM <- sum(cp_genes %in% shallow_genes_HEPM)


# Distribution of CP genes in MEPM
cp_genes_mm <- Combined_CPGene %>%
  filter(!is.na(mm_gene_symbol)) %>%
  pull(mm_gene_symbol)  %>% unique()

cp_deep_MEPM <- sum(cp_genes_mm %in% deep_genes_MEPM)
cp_shallow_MEPM <- sum(cp_genes_mm %in% shallow_genes_MEPM)

cp_genes_data <- data.frame(
  Group = rep(c("HEPM", "MEPM"), each = 2),
  Type = rep(c("Shallow", "Deep"), 2),
  Count = c(
    cp_shallow_HEPM,
    cp_deep_HEPM,
    cp_shallow_MEPM,
    cp_deep_MEPM
  ),
  Category = "CP genes"
)

# 3. merge data
plot_data <- rbind(all_genes_data, cp_genes_data)


plot_data$Type <- factor(plot_data$Type, levels = c("Shallow", "Deep"))
plot_data$Group <- factor(plot_data$Group, levels = c("HEPM", "MEPM"))
plot_data$Category <- factor(plot_data$Category, levels = c("All genes", "CP genes"))

# 4. plot
p <- ggplot(plot_data, aes(x = Group, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~Category, scales = "free_y") +
  scale_fill_manual(values = c("Shallow" = "#87CEEB", "Deep" = "#FF6B6B")) +
  labs(
    title = "Distribution of Shallow and Deep Genes",
    x = "",
    y = "Number of Genes",
    fill = "Gene Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  )

print(p)
print(plot_data)

# calculate
plot_data_pct <- plot_data %>%
  group_by(Category, Group) %>%
  mutate(
    Total = sum(Count),
    Percentage = Count / Total * 100
  )

# plot
p_pct <- ggplot(plot_data_pct, aes(x = Group, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~Category) +
  scale_fill_manual(values = c("Shallow" = "#87CEEB", "Deep" = "#FF6B6B")) +
  labs(
    title = "Percentage Distribution of Shallow and Deep Genes",
    x = "",
    y = "Percentage (%)",
    fill = "Gene Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  )

print(p_pct)
print(plot_data_pct)
ggsave("Gene_Distribution_Stacked_Barplot_number.pdf", p, width = 10, height = 6, dpi = 300)
ggsave("Gene_Distribution_Stacked_Barplot_percent.pdf", p_pct, width = 10, height = 6, dpi = 300)
write.csv(plot_data_pct, "statis_summary_Fig3F.csv")
# --------- Fig3h: plot Tier1 and Tier2 enhancer number of top10 CP genes in HEPM -----

# load HEPM CP genes
hepm_cp_genes <- Combined_CPGene$CP_Gene %>% unique()

cat("HEPM CP genes:", length(hepm_cp_genes), "\n")

# combined deep and shallow data in HEPM
hepm_combined <- Gene_regulatory_layers_HEPM %>%
  mutate(Type = Regulation_Type) %>%
  dplyr::select(Gene, Tier1_Count, Tier2_Count, Total_Enhancers, Type)

# label CP genes in HEPM
hepm_combined <- hepm_combined %>%
  mutate(
    Is_CP = Gene %in% hepm_cp_genes,
    Dataset = "HEPM"
  )

hepm_cp_data <- hepm_combined %>%
  filter(Is_CP)

cat("CP gene numbers of HEPM:", nrow(hepm_cp_data), "\n") # 424

# Top10 CP genes in HPEM cp data
hepm_top10 <- hepm_cp_data %>%
  arrange(desc(Total_Enhancers)) %>%
  head(10)

cat("HEPM Top 10 CP gene:\n")
print(hepm_top10 %>% dplyr::select(Gene, Total_Enhancers, Type))

# label Top10 CP genes in HEPM combined data
hepm_combined <- hepm_combined %>%
  mutate(
    Is_Top10_CP = Gene %in% hepm_top10$Gene,
    Gene_Label = ifelse(Is_Top10_CP, Gene, ""),
    Tier1_Plot = Tier1_Count
  )
# load gene in MEPM

mepm_cp_genes <- Combined_CPGene$mm_gene_symbol %>% unique()
mepm_cp_genes <- mepm_cp_genes[!is.na(mepm_cp_genes)]

cat("\nMEPM CP gene:", length(mepm_cp_genes), "\n")

# combined Deep/Shallow data in MEPM
mepm_combined <- Gene_regulatory_layers_MEPM %>%
  mutate(Type = Regulation_Type) %>%
  dplyr::select(Gene, Tier1_Count, Tier2_Count, Total_Enhancers, Type)

# label CP genes in MEPM
mepm_combined <- mepm_combined %>%
  mutate(
    Is_CP = Gene %in% mepm_cp_genes,
    Dataset = "MEPM"
  )

# filter CP gene in MEPM
mepm_cp_data <- mepm_combined %>%
  filter(Is_CP)

cat("CP genes number of MEPM:", nrow(mepm_cp_data), "\n") #462

# find Top10 CP genes in MEPM
mepm_top10 <- mepm_cp_data %>%
  arrange(desc(Total_Enhancers)) %>%
  head(10)

cat("MEPM Top 10 CP gene:\n")
print(mepm_top10 %>% dplyr::select(Gene, Total_Enhancers, Type))

# label MEPM Top10 gene
mepm_combined <- mepm_combined %>%
  mutate(
    Is_Top10_CP = Gene %in% mepm_top10$Gene,
    Gene_Label = ifelse(Is_Top10_CP, Gene, ""),
    Tier1_Plot = -Tier1_Count 
  )

# combined MEPM and HEPM data
combined_plot_data <- rbind(hepm_combined, mepm_combined)


# X axis = max|abs|
x_max <- max(abs(combined_plot_data$Tier1_Plot), na.rm = TRUE)
x_limits <- c(-x_max * 1.1, x_max * 1.1)

# Create X-axis breaks and labels (all displayed as positive values)
x_breaks <- seq(-ceiling(x_max), ceiling(x_max), by = 5)
x_labels <- abs(x_breaks)

cat("\nX axis:", x_limits, "\n")

# mirror dotplot
# custom gene
custom_genes <- c("SNAI2", "Snai2", "HHAT", "Hhat")
combined_plot_data <- combined_plot_data %>% mutate(Is_Custom = Gene %in% custom_genes)
library(ggrepel)
p_mirror <- ggplot(combined_plot_data %>% filter(Is_CP), aes(x = Tier1_Plot, y = Tier2_Count)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8, linetype = "solid") +
  # background color
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#FDE6CC", alpha = 0.3) +  # MEPM
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#E6F4FA", alpha = 0.3) +  # HEPM 
  geom_point(data = combined_plot_data %>% filter(Is_CP & !Is_Top10_CP),
             aes(color = Type), 
             alpha = 0.5, size = 3) +
  geom_point(data = combined_plot_data %>% filter(Is_Top10_CP),
             aes(color = Type), 
             size = 3.5, alpha = 1) +
  geom_point(data = combined_plot_data %>% filter(Is_Top10_CP),
             color = "black", size = 3.5, shape = 1, stroke = 0.6) +
  geom_point(data = combined_plot_data %>% filter(Is_Custom),
             aes(color = Type), 
             size = 3.5, alpha = 1) +
  geom_point(data = combined_plot_data %>% filter(Is_Custom),
             color = "red", size = 3.5, shape = 1, stroke = 0.8) +
  
  scale_color_manual(
    values = c(
      "Deep" = "#C9A063",     
      "Shallow" = "black"      
    ),
    labels = c(
      "Deep" = "Deep",
      "Shallow" = "Shallow"
    ),
    name = "Regulation Type"
  ) +
  # HEPM gene (right)
  geom_text_repel(
    data = combined_plot_data %>% filter(Is_Top10_CP & Dataset == "HEPM"),
    aes(label = Gene_Label),
    size = 6,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    arrow = arrow(length = unit(0.01, "npc")),
    max.overlaps = 30,
    force = 2,
    min.segment.length = 0,
    nudge_x = 2,
    direction = "y"
  ) +
  # MEPM gene (left)
  geom_text_repel(
    data = combined_plot_data %>% filter(Is_Top10_CP & Dataset == "MEPM"),
    aes(label = Gene_Label),
    size = 6,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    arrow = arrow(length = unit(0.01, "npc")),
    max.overlaps = 30,
    force = 2,
    min.segment.length = 0,
    nudge_x = -2,
    direction = "y"
  ) +
  scale_x_continuous(
    breaks = seq(-35, 35, by = 5),
    labels = c(35, 30, 25, 20, 15, 10, 5, 0, 5, 10, 15, 20, 25, 30, 35),
    limits = c(-37, 37)
  ) +
  labs(
    title = "Mirror Comparison: HEPM vs MEPM (CP Genes Only)",
    subtitle = "Tier 1 vs Tier 2 Enhancers for Top 10 CP Genes",
    x = "Tier 1 Enhancers Count",
    y = "Tier 2 Enhancers Count"
  ) +
  annotate("text", x = -28, y = max(combined_plot_data$Tier2_Count) * 0.8,
           label = "MEPM", size = 6, fontface = "bold", color = "#E78712") +
  annotate("text", x = 28, y = max(combined_plot_data$Tier2_Count) * 0.8,
           label = "HEPM", size = 6, fontface = "bold", color = "#82C9E6") +

  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8)
  )

print(p_mirror)

ggsave("HEPM_MEPM_Mirror_CPGenes_Only.pdf", p_mirror, width = 14, height = 8)
# check top 10 CP gene
cat("\n========================================\n")
cat("HEPM Top 10 CP Genes\n")
cat("========================================\n")
hepm_top10_detail <- hepm_combined %>% 
  filter(Is_Top10_CP) %>% 
  dplyr::select(Gene, Tier1_Count, Tier2_Count, Total_Enhancers, Type) %>%
  arrange(desc(Total_Enhancers))
print(hepm_top10_detail)

cat("\n========================================\n")
cat("MEPM Top 10 CP Genes\n")
cat("========================================\n")
mepm_top10_detail <- mepm_combined %>% 
  filter(Is_Top10_CP) %>% 
  dplyr::select(Gene, Tier1_Count, Tier2_Count, Total_Enhancers, Type) %>%
  arrange(desc(Total_Enhancers))
print(mepm_top10_detail)

# ----- save Rdata --
save(list = ls(), file = "Fig3_SEs_network.Rdata")
load("Fig3_SEs_network.Rdata")
