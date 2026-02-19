setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results")

library(ggplot2)
library(ggpubr)
library(dplyr)

# load data 
load("./OMIM_cleftpalateGene_microC/Fig2J_CombinedCPGenes_chromatinLoops.Rdata")

# caculate frequency for all gene in HEPM
all_genes <- unique(unlist(strsplit(
  apply(Data1, 1, function(row) {
    genes <- c(row["P_L_genes"], row["P_R_genes"])
    genes <- genes[!is.na(genes) & genes != ""]
    paste(genes, collapse = ",")
  }), 
  ","
)))
all_genes <- all_genes[all_genes != "" & !is.na(all_genes)]
All_Gene <- data.frame(CP_Gene = all_genes)

HEPM_frequency <- detailed_CPgene_analysis(Data1, All_Gene, gene_col = "CP_Gene")

# select conserved genes
conserved_loops <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/processed_results/3.compare_loop/Total_conservedLoop/HEPM_conserved_chromatinLoops_CPGenes.csv")
conserved_Genes <- unique(unlist(strsplit(
  apply(conserved_loops, 1, function(row) {
    genes <- c(row["P_L_genes"], row["P_R_genes"])
    genes <- genes[!is.na(genes) & genes != ""]
    paste(genes, collapse = ",")
  }), 
  ","
)))
conserved_Genes <- conserved_Genes[conserved_Genes != "" & !is.na(conserved_Genes)]

# select conserved CP genes
conserved_CPGenes <- Combined_CPGene[Combined_CPGene$CP_Gene %in% conserved_Genes, ]

# prepare data for plot
hepm_genes <- HEPM_frequency$gene_frequency$gene
detailed_genes <- detailed_result$gene_frequency$gene

# group 1：Conserved CP Genes
conserved_cp_data <- detailed_result$gene_frequency[
  detailed_result$gene_frequency$gene %in% conserved_CPGenes$CP_Gene, 
]

# group 2：CP Genes
all_cp_data <- detailed_result$gene_frequency

# group 3：Other Genes
hepm_filtered <- HEPM_frequency$gene_frequency[!hepm_genes %in% detailed_genes, ]

# frequency for three groups
cat("group1 - Conserved CP Genes frequency:\n")
print(summary(conserved_cp_data$count))
cat("\ngroup2 - CP Genes frequency:\n")
print(summary(all_cp_data$count))
cat("\ngroup3 - Other Genes frequency:\n")
print(summary(hepm_filtered$count))
cat("\n")

write.csv(conserved_cp_data, "Group1_Conserved_CP_Genes.csv", row.names = FALSE)
write.csv(all_cp_data, "Group2_All_CP_Genes.csv", row.names = FALSE)
write.csv(hepm_filtered, "Group3_Other_Genes.csv", row.names = FALSE)

# combine data
data_combined <- data.frame(
  Frequency = c(conserved_cp_data$count,
                all_cp_data$count,
                hepm_filtered$count),
  Group = c(rep("Conserved CP Genes", nrow(conserved_cp_data)),
            rep("CP Genes", nrow(all_cp_data)),
            rep("Other Genes", nrow(hepm_filtered)))
)
data_combined$Group <- factor(data_combined$Group, 
                              levels = c("Conserved CP Genes", "CP Genes", "Other Genes"))


# statistic summary
detailed_stats <- data_combined %>%
  group_by(Group) %>%
  summarise(
    样本数 = n(),
    最小值 = min(Frequency),
    第一四分位数 = quantile(Frequency, 0.25),
    中位数 = median(Frequency),
    平均数 = mean(Frequency),
    第三四分位数 = quantile(Frequency, 0.75),
    最大值 = max(Frequency),
    标准差 = sd(Frequency),
    方差 = var(Frequency),
    变异系数 = sd(Frequency) / mean(Frequency) * 100,
    标准误 = sd(Frequency) / sqrt(n()),
    .groups = 'drop'
  )

detailed_stats[, 2:12] <- round(detailed_stats[, 2:12], 2)
print(detailed_stats)

# violin plot
p <- ggplot(data_combined, aes(x = Group, y = Frequency, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  stat_compare_means(comparisons = list(
    c("Conserved CP Genes", "CP Genes"),
    c("Conserved CP Genes", "Other Genes"),
    c("CP Genes", "Other Genes")
  ), method = "wilcox.test", label = "p.format", size = 4) +
  theme_bw() +
  labs(title = "Gene Frequency Distribution",
       x = "",
       y = "Gene Frequency (Count)") +
  scale_fill_manual(values = c(
    "Conserved CP Genes" = "#D55E00",
    "CP Genes" = "#E69F00",
    "Other Genes" = "#56B4E9"
  )) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 13))

print(p)
ggsave("./OMIM_cleftpalateGene_microC/Three_groups_violin_stats_HEPM_loopsfrequency.pdf", width = 6, height = 6)

cat(paste("Total", nrow(conserved_CPGenes), "Conserved CP Gene:\n\n"))

write.csv(conserved_CPGenes, "Conserved_CP_Genes_list.csv", row.names = FALSE)

gene_classification <- data.frame(
  Gene = c(conserved_cp_data$gene, all_cp_data$gene, hepm_filtered$gene),
  Frequency = c(conserved_cp_data$count, all_cp_data$count, hepm_filtered$count),
  Category = c(rep("Conserved CP Genes", nrow(conserved_cp_data)),
               rep("CP Genes", nrow(all_cp_data)),
               rep("Other Genes", nrow(hepm_filtered)))
)

write.csv(gene_classification, "Fig3L.Gene_Classification_Summary.csv", row.names = FALSE)

save(list = ls(), file = "Fig2L_chromatinLoops_frequency_3group.Rdata")
