library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(tidyr)
library(ggpubr) 
save(list = ls(), file = "Fig.8.Identify_alleleSpecific_TFbSNP.Rdata")

setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/24.TFbindingSNP/Results/HEPM/2.TOBIAS_results/3.TFbSNP")
source("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/24.TFbindingSNP/script/4.TFbSNP_deltaScore_analysis.R")
result <- analyze_TFbSNP(
   input_csv = "TFbSNP_CPae_footprint.csv",
   genome_version = "hg19",
   flank_size = 9,
   svm_model_prefix = "/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/13_gkmSVM/Final_input_atac/full_sequence/atac_full",
   output_prefix = "./HEPM_CPae_TFbSNP_score"
 )


#----- Plot Fig 8C -----------------------
TfbSNP_gkmSVM <- read.csv("./HEPM_CPae_TFbSNP_score_results.csv")

# -------------------------------------------------- Identify AEs and targeted genes for allele-specific TFbSNP ---------------------------------------
library(GenomicRanges)
CP_enhancers <- read.delim("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/24.TFbindingSNP/Results/HEPM/1.CP_Enhancer/HEPM_CP_Enhancers_Combined.bed", header = T)
HEPM_genes_Enhancers <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/HEPM/HEPM_gene_regulatory_layers.csv")
Combined_CPGene <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/CP_combined_with_mouse_orthologs.csv", sep=",", header = T)
Combined_CPGene_HEPM <- Combined_CPGene$CP_Gene %>% unique()
# mapping enhancer to CP genes
HEPM_CPGene_Enhancers <- HEPM_genes_Enhancers %>%
  filter(Gene %in% Combined_CPGene_HEPM)

CP_enhancers_CPgenes <- CP_enhancers %>%
  left_join(
    HEPM_CPGene_Enhancers %>%
      pivot_longer(cols = c(Tier1_Enhancers, Tier2_Enhancers),
                   values_to = "Original_ID") %>%
      separate_longer_delim(Original_ID, delim = ";") %>%
      filter(!is.na(Original_ID) & Original_ID != "") %>%
      mutate(Original_ID = trimws(Original_ID)) %>%
      group_by(Original_ID) %>%
      summarise(CP_Genes = paste(sort(unique(Gene)), collapse = ";"), .groups = "drop"),
    by = "Original_ID"
  ) %>%
  mutate(CP_Genes = replace_na(CP_Genes, ""))



# SNP GRanges
snp_gr <- GRanges(
  seqnames = TfbSNP_gkmSVM$CHROM,
  ranges = IRanges(start = TfbSNP_gkmSVM$POS, end = TfbSNP_gkmSVM$POS),
  strand = "*"
)
mcols(snp_gr) <- TfbSNP_gkmSVM

# Enhancer GRanges
enhancer_gr <- GRanges(
  seqnames = CP_enhancers_CPgenes$chr,
  ranges = IRanges(start = CP_enhancers_CPgenes$start, end = CP_enhancers_CPgenes$end),
  strand = "*"
)
mcols(enhancer_gr) <- CP_enhancers_CPgenes

# 2. Finding overlaps
cat("Finding overlaps...\n")
overlaps <- findOverlaps(snp_gr, enhancer_gr)

# 3. mapping SNP to genes
snp_to_genes <- data.frame(
  snp_idx = queryHits(overlaps),
  enhancer_idx = subjectHits(overlaps)
) %>%
  mutate(
    SNP_ID = TfbSNP_gkmSVM$ID[snp_idx],
    CP_Genes = CP_enhancers_CPgenes$CP_Genes[enhancer_idx],
    Enhancer_ID = CP_enhancers_CPgenes$Enhancer_ID[enhancer_idx],
    Tier = CP_enhancers_CPgenes$Tier[enhancer_idx]
  ) %>%
  group_by(SNP_ID) %>%
  summarise(
    CP_Genes = paste(unique(CP_Genes[CP_Genes != ""]), collapse = ";"),
    Enhancer_IDs = paste(unique(Enhancer_ID), collapse = ";"),
    Tiers = paste(unique(Tier), collapse = ";"),
    n_enhancers = n(),
    .groups = "drop"
  )

# 4. return to raw data
TfbSNP_gkmSVM <- TfbSNP_gkmSVM %>%
  left_join(snp_to_genes, by = c("ID" = "SNP_ID")) %>%
  mutate(
    CP_Genes = replace_na(CP_Genes, ""),
    Enhancer_IDs = replace_na(Enhancer_IDs, ""),
    Tiers = replace_na(Tiers, ""),
    n_enhancers = replace_na(n_enhancers, 0)
  )

snps_no_cpgenes <- TfbSNP_gkmSVM %>%
  filter(CP_Genes == "" | is.na(CP_Genes))

snps_with_cpgenes <- TfbSNP_gkmSVM %>%
  filter(CP_Genes != "" & !is.na(CP_Genes))

# ====== plot 8C =============
freq_stats <- snps_with_cpgenes %>%
  dplyr::count(FREQUENCY_CLASS, name = "Count") %>%
  dplyr::mutate(
    Percentage = Count / sum(Count) * 100,
    PercentLabel = paste0(round(Percentage, 1), "%"),
    Ypos = cumsum(Count) - 0.5 * Count,
    CumPercentage = cumsum(Percentage)
  ) %>%
  dplyr::arrange(desc(FREQUENCY_CLASS))

print(freq_stats)

# pie 
pie_plot <- ggplot(freq_stats, aes(x = "", y = Count, fill = FREQUENCY_CLASS)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.1) +
  coord_polar("y", start = 0) +
  geom_text(
    aes(y = Ypos, 
        label = paste0(FREQUENCY_CLASS, "\n", 
                       format(Count, big.mark = ","))),
    size = 5, 
    fontface = "bold", 
    color = "white"
  ) +
  geom_text(
    aes(y = Ypos, label = PercentLabel),
    position = position_nudge(x = 0.7),
    size = 6, 
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "RARE" = "#e74c3c",    
      "COMMON" = "#3498db"   
    ),
    name = "Frequency Class"
  ) +
  theme_void() +
  theme(plot.title = element_text(
    hjust = 0.5, 
    size = 20, 
    face = "bold", 
    margin = margin(b = 10)
  ),
  plot.subtitle = element_text(
    hjust = 0.5, 
    size = 14, 
    color = "gray40",
    margin = margin(b = 20)
  ),
  legend.position = "bottom",
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 13),
  legend.key.size = unit(1, "cm"),
  plot.margin = margin(20, 20, 20, 20),
  plot.background = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "SNP Frequency Class Distribution",
    subtitle = paste0("Total SNPs: ", 
                      format(nrow(snps_with_cpgenes), big.mark = ","))
  )

print(pie_plot)

ggsave("Fig.8C.TFbSNP_Frequency_Distribution_Pie.pdf", 
       pie_plot, 
       width = 10, height = 10, dpi = 300)

# ============================================================ plot figure 8D ===========================================
snps_with_cpgenes <- snps_with_cpgenes %>%
  mutate(RELATIVE_FOOTPRINT_SCORE_norm = (RELATIVE_FOOTPRINT_SCORE - min(RELATIVE_FOOTPRINT_SCORE, na.rm = TRUE)) / 
           (max(RELATIVE_FOOTPRINT_SCORE, na.rm = TRUE) - min(RELATIVE_FOOTPRINT_SCORE, na.rm = TRUE)))

plot_data <- snps_with_cpgenes %>%
  filter(!is.na(deltaScore), 
         !is.na(RELATIVE_FOOTPRINT_SCORE),
         FREQUENCY_CLASS %in% c("RARE", "COMMON"))

# 
q90_deltaScore <- quantile(plot_data$deltaScore, 0.90, na.rm = TRUE)
q10_deltaScore <- quantile(plot_data$deltaScore, 0.10, na.rm = TRUE)

# label SNP list
custom_snp_list <- c("rs592164", "rs56268268", "rs193044641")
custom_labeled <- plot_data %>%
  filter(ID %in% custom_snp_list) %>%
  mutate(label_type = "Custom")
#top10 deltaSVM score
top10_snps <- plot_data %>%
  arrange(desc(deltaScore)) %>%
  head(10) %>%
  mutate(label_type = "Top10")
# Bottom 10 deltaSVM score 
bottom10_snps <- plot_data %>%
  arrange(deltaScore) %>%
  head(10) %>%
  mutate(label_type = "Bottom10")
all_labeled_snps <- bind_rows(
  custom_labeled
) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(
    #  Custom > Top10 > Bottom10
    label_type = case_when(
      ID %in% custom_snp_list ~ "Custom",
      ID %in% top10_snps$ID ~ "Top10",
      ID %in% bottom10_snps$ID ~ "Bottom10",
      TRUE ~ label_type
    )
  )


# plot
p_custom <- ggplot(plot_data, aes(x = deltaScore, y = RELATIVE_FOOTPRINT_SCORE)) +
  geom_point(aes(color = FREQUENCY_CLASS), alpha = 0.6, size = 4, shape = 16) +
  geom_point(data = all_labeled_snps,
             aes(fill = FREQUENCY_CLASS),    
             size = 6,                      
             shape = 21,                     
             stroke = 2,                     
             color = "black") +             
  geom_text_repel(
    data = all_labeled_snps,
    aes(label = ID),                       
    size = 3.5, 
    fontface = "bold", 
    color = "black",                        
    box.padding = 0.8, 
    point.padding = 0.5,
    segment.color = "black",                
    segment.size = 0.8,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed"),
    max.overlaps = 100, 
    min.segment.length = 0,
    force = 15,
    force_pull = 2
  ) +
  scale_color_manual(
    values = c("RARE" = "#e74c3c", "COMMON" = "#3498db"),
    name = "Frequency Class"
  ) +
  scale_fill_manual(
    values = c("RARE" = "#e74c3c", "COMMON" = "#3498db"),
    name = "Frequency Class"
  ) +
  geom_vline(xintercept = q90_deltaScore, 
             linetype = "dashed", color = "gray30", linewidth = 0.8) +
  geom_vline(xintercept = q10_deltaScore, 
             linetype = "dashed", color = "gray30", linewidth = 0.8) +
  annotate("text", 
           x = q90_deltaScore, 
           y = max(plot_data$RELATIVE_FOOTPRINT_SCORE_norm) * 0.95,
           label = paste0("Q90 = ", round(q90_deltaScore, 2)),
           hjust = -0.1, vjust = 1, 
           color = "gray30", size = 4, fontface = "bold") +
  annotate("text", 
           x = q10_deltaScore, 
           y = max(plot_data$RELATIVE_FOOTPRINT_SCORE_norm) * 0.95,
           label = paste0("Q10 = ", round(q10_deltaScore, 2)),
           hjust = 1.1, vjust = 1, 
           color = "gray30", size = 4, fontface = "bold") +
  labs(
    title = "Relationship between deltaSVM Score and Footprint Score",
    subtitle = paste0("Total SNPs: ", format(nrow(plot_data), big.mark = ","),
                      " | RARE: ", sum(plot_data$FREQUENCY_CLASS == "RARE"),
                      " | COMMON: ", sum(plot_data$FREQUENCY_CLASS == "COMMON"),
                      " | Labeled: ", nrow(all_labeled_snps)),
    x = "deltaSVM Score",
    y = "Relative Footprint Score"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11),
    plot.margin = margin(20, 20, 20, 20)
  )


print(p_custom)

ggsave("Fig.8D.TfbSNP_scatter_custom_labeled.pdf", p_custom, 
       width = 13, height = 10, dpi = 300)

# Filter SNP >= Q10/Q90
Top_candidate_TFbSNP <- plot_data %>%
  filter(deltaScore > q90_deltaScore | deltaScore < q10_deltaScore) %>%
  mutate(
    candidate_type = case_when(
      deltaScore > q90_deltaScore ~ "High_deltaScore (>Q90)",
      deltaScore < q10_deltaScore ~ "Low_deltaScore (<Q10)",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(desc(abs(deltaScore)))  

write.csv(Top_candidate_TFbSNP, 
          "Top_candidate_TFbSNP.csv", 
          row.names = FALSE)

write.csv(snps_with_cpgenes, "Table17.TFbSNP_allele_gkmSVM_Enhancer_CPGene.csv", row.names = FALSE)


############################ Figure 8 E. ##################################################
neg_SNP <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/24.TFbindingSNP/Results/HEPM/2.TOBIAS_results/3.TFbSNP/Negcontrol/control_SNPs_full.csv")

library(GenomicRanges)

# 将SNP转换为GRanges对象
snp_gr <- GRanges(
  seqnames = neg_SNP$CHROM,
  ranges = IRanges(start = neg_SNP$POS, width = 1),
  strand = "*"
)

# 将enhancer转换为GRanges对象
enhancer_gr <- GRanges(
  seqnames = CP_enhancers$chr,
  ranges = IRanges(start = CP_enhancers$start, end = CP_enhancers$end),
  strand = "*"
)

# 找到重叠的SNP
overlaps <- findOverlaps(snp_gr, enhancer_gr)

# 获取在enhancer内的SNP索引
snp_in_enhancer_idx <- unique(queryHits(overlaps))

# 过滤掉在enhancer内的SNP
neg_SNP_filtered <- neg_SNP[-snp_in_enhancer_idx, ]

# 查看被过滤掉的SNP（可选）
snp_removed <- neg_SNP[snp_in_enhancer_idx, ]


#-- 
source("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/24.TFbindingSNP/Results/HEPM/2.TOBIAS_results/3.TFbSNP/Negcontrol/Control_SNP_deltaSCORE.R")

result <- analyze_control_SNPs_gkmSVM(
  control_snps = neg_SNP,
  genome_version = "hg19",
  flank_size = 9,
  ref_fasta = "control_reference_alleles.fa",
  alt_fasta = "control_alt_alleles.fa",
  outfn_snp = "control_deltaSVM_score.out",
  svm_model_prefix = "/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/13_gkmSVM/Final_input_atac/full_sequence/atac_full",
  output_prefix = "control_SNP_analysis",
  verbose = TRUE
)

# 
## input file
testfn_ref= 'control_reference_alleles_unique.fa'   
testfn_alt= 'control_alt_alleles_unique.fa'    
svmfnprfx = '/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/13_gkmSVM/Final_input_atac/full_sequence/atac_full'
## output file 
outfn_snp = 'control_score.out'

## deltaSVM get score for SNPs
gkmsvm_delta(testfn_ref,testfn_alt, svmfnprfx, outfn_snp)


#---  Plot Fig8 TFbSNP/control -----

select_controls <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/24.TFbindingSNP/Results/HEPM/2.TOBIAS_results/3.TFbSNP/Negcontrol/Enhancer_Activity_Data.csv")
select_controls <- select_controls %>% filter(Type == "Control")

candidate_data <- Top_candidate_TFbSNP %>%
  select(ID, deltaScore) %>%
  mutate(
    Type = "Candidate",
    Color = "Candidate SNP"
  )

# conbined data
combined_data <- bind_rows(candidate_data, select_controls) %>%
  arrange(deltaScore) %>%
  mutate(x_position = row_number())

#  SNP group ：Enhancing、Disrupting、Neutral
combined_data <- combined_data %>%
  mutate(
    abs_deltaScore = abs(deltaScore),
    Effect_Direction = case_when(
      deltaScore > 0 ~ "Enhancing",
      deltaScore < 0 ~ "Disrupting",
      TRUE ~ "Neutral"
    )
  )


stats_summary <- combined_data %>%
  group_by(Type) %>%
  summarise(
    n = n(),
    mean_delta = mean(deltaScore),
    median_delta = median(deltaScore),
    sd_delta = sd(deltaScore),
    mean_abs_delta = mean(abs_deltaScore),
    median_abs_delta = median(abs_deltaScore),
    n_enhancing = sum(deltaScore > 0),
    n_disrupting = sum(deltaScore < 0),
    pct_enhancing = sum(deltaScore > 0) / n() * 100,
    pct_disrupting = sum(deltaScore < 0) / n() * 100,
    .groups = "drop"
  )


cat(rep("-", 80), "\n", sep = "")
print(stats_summary)
cat("\n")


stats_by_direction <- combined_data %>%
  group_by(Type, Effect_Direction) %>%
  summarise(
    n = n(),
    mean_delta = mean(deltaScore),
    median_delta = median(deltaScore),
    sd_delta = sd(deltaScore),
    min_delta = min(deltaScore),
    max_delta = max(deltaScore),
    .groups = "drop"
  ) %>%
  arrange(Type, desc(Effect_Direction))


cat(rep("-", 80), "\n", sep = "")
print(stats_by_direction)
cat("\n")


# plot 8E
direction_summary <- combined_data %>%
  group_by(Type, Effect_Direction) %>%
  summarise(
    n = n(),
    mean_delta = mean(deltaScore),
    se_delta = sd(deltaScore) / sqrt(n()),
    .groups = "drop"
  ) %>%
  filter(Effect_Direction != "Neutral")
my_comparisons <- list(c("Candidate", "Control"))

p_facet <- ggplot(
  combined_data %>% filter(Effect_Direction != "Neutral"),
  aes(x = Type, y = deltaScore, fill = Type)
) +
  geom_violin(alpha = 0.4, width = 0.8) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  
  # wilcox.test
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",     
    label = "p.signif",        
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  
  stat_summary(
    fun = mean, geom = "point", shape = 23, size = 4,
    fill = "white", color = "black", stroke = 1.5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  
  facet_wrap(
    ~ Effect_Direction,
    scales = "free_y",
    labeller = labeller(
      Effect_Direction = c(
        "Enhancing" = "Enhancing Activity (deltaSVM > 0)",
        "Disrupting" = "Disrupting Activity (deltaSVM < 0)"
      )
    )
  ) +
  
  scale_fill_manual(values = c("Candidate" = "#E41A1C", "Control" = "#808080")) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Enhancer Activity Modulation: Candidate vs Control SNPs",
    x = "SNP Type",
    y = "deltaSVM Score"
  )


print(p_facet)

ggsave(
  "Fig.8E.Enhancer_Activity_Faceted.pdf",
  plot = p_facet,
  width = 10,
  height = 10,
  dpi = 300
)

write.csv(
  combined_data,
  "Enhancer_Activity_Data.csv",
  row.names = FALSE
)
