
setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results")
getwd()
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotly)
library(scales)
library(stringr)
library(GenomicRanges)

# --- HEPM chromatin loops with  VISTA positive human facial messenchymal tissue enhancer  ---
Data1 <- read.csv("./1.Loop_annotation/HEPM_microC_6400bp_annotate_loops_20251204.csv", header = T, sep = ",")
fm_enh <- read.delim("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/13_gkmSVM/VistaEnh/VISTA_ENH_250222/VISTA_positive_fm_hg19.bed", header = F)

fm_enh_gr <- GRanges(
  seqnames = fm_enh$V1,
  ranges = IRanges(start = fm_enh$V2, end = fm_enh$V3)
)


micro_gr1 <- GRanges(
  seqnames = Data1$chr1,
  ranges = IRanges(start = Data1$start1, end = Data1$end1)
)
mcols(micro_gr1)$original_row <- 1:nrow(Data1)
mcols(micro_gr1)$region <- "region1"

# convert anchor 2 in HEPM Micro-C to GRanges object
micro_gr2 <- GRanges(
  seqnames = Data1$chr2,
  ranges = IRanges(start = Data1$start2, end = Data1$end2)
)
mcols(micro_gr2)$original_row <- 1:nrow(Data1)
mcols(micro_gr2)$region <- "region2"

# find fm_enh overlapped with micro_gr1/micro_gr2
overlaps1 <- findOverlaps(fm_enh_gr, micro_gr1)
overlaps2 <- findOverlaps(fm_enh_gr, micro_gr2)

# dataframe for overlapped data
overlap_info1 <- data.frame(
  fm_idx = queryHits(overlaps1),
  micro_row = mcols(micro_gr1)$original_row[subjectHits(overlaps1)],
  region = "region1"
)

overlap_info2 <- data.frame(
  fm_idx = queryHits(overlaps2),
  micro_row = mcols(micro_gr2)$original_row[subjectHits(overlaps2)],
  region = "region2"
)

all_overlaps <- rbind(overlap_info1, overlap_info2)
all_original_rows <- unique(c(overlap_info1$micro_row, overlap_info2$micro_row))

# select HEPM Micro-C data
Data1_overlapping <- Data1[all_original_rows, ]

# add fm_enh data
Data1_overlapping$fm_chr <- ""
Data1_overlapping$fm_start <- ""
Data1_overlapping$fm_end <- ""
Data1_overlapping$fm_ID <- "" 
Data1_overlapping$Intersected_Anchor <- ""  
Data1_overlapping$Targeted_Gene <- ""       

for (i in 1:nrow(Data1_overlapping)) {
  row_idx <- all_original_rows[i]
  fm_indices_region1 <- overlap_info1$fm_idx[overlap_info1$micro_row == row_idx]
  fm_indices_region2 <- overlap_info2$fm_idx[overlap_info2$micro_row == row_idx]
  all_fm_indices <- unique(c(fm_indices_region1, fm_indices_region2))
  if (length(all_fm_indices) > 0) {
    Data1_overlapping$fm_chr[i] <- paste(fm_enh$V1[all_fm_indices], collapse = ",")
    Data1_overlapping$fm_start[i] <- paste(fm_enh$V2[all_fm_indices], collapse = ",")
    Data1_overlapping$fm_end[i] <- paste(fm_enh$V3[all_fm_indices], collapse = ",")
    Data1_overlapping$fm_ID[i] <- paste(fm_enh$V4[all_fm_indices], collapse = ",")
  }
  
  # identify Intersected_Anchor and Targeted_Gene
  has_region1_overlap <- length(fm_indices_region1) > 0
  has_region2_overlap <- length(fm_indices_region2) > 0
  
  if (has_region1_overlap && has_region2_overlap) {
    # case 1: overlap with two anchors
    Data1_overlapping$Intersected_Anchor[i] <- "Both"
    # combined target genes
    p_r_gene <- ifelse(is.na(Data1_overlapping$P_R_genes[i]) || Data1_overlapping$P_R_genes[i] == "", 
                       "", as.character(Data1_overlapping$P_R_genes[i]))
    p_l_gene <- ifelse(is.na(Data1_overlapping$P_L_genes[i]) || Data1_overlapping$P_L_genes[i] == "", 
                       "", as.character(Data1_overlapping$P_L_genes[i]))
    combined_genes <- paste(c(p_r_gene, p_l_gene)[c(p_r_gene, p_l_gene) != ""], collapse = ",")
    Data1_overlapping$Targeted_Gene[i] <- combined_genes
    
  } else if (has_region1_overlap) {
    # case 2: overlap with left anchor  -> Left anchor, target gene from P_R_genes
    Data1_overlapping$Intersected_Anchor[i] <- "Left"
    Data1_overlapping$Targeted_Gene[i] <- ifelse(is.na(Data1_overlapping$P_R_genes[i]) || Data1_overlapping$P_R_genes[i] == "", 
                                                 "", as.character(Data1_overlapping$P_R_genes[i]))
    
  } else if (has_region2_overlap) {
    # case 3: overlap with right anchor  -> Right anchor, target gene from P_L_genes
    Data1_overlapping$Intersected_Anchor[i] <- "Right"
    Data1_overlapping$Targeted_Gene[i] <- ifelse(is.na(Data1_overlapping$P_L_genes[i]) || Data1_overlapping$P_L_genes[i] == "", 
                                                 "", as.character(Data1_overlapping$P_L_genes[i]))
  }
}

# check 
print(paste("found", nrow(Data1_overlapping), "loops overlapped with fm_enh"))

cat("\n=== Intersected_Anchor Distribution ===\n")
print(table(Data1_overlapping$Intersected_Anchor, useNA = "ifany"))

cat("\n=== Targeted_Gene Summary ===\n")
cat("Rows with targeted genes:", sum(Data1_overlapping$Targeted_Gene != "", na.rm = TRUE), "\n")
cat("Rows without targeted genes:", sum(Data1_overlapping$Targeted_Gene == "", na.rm = TRUE), "\n")

# check the subtype distribution: Vista-type.

Data1_overlapping <- Data1_overlapping %>%
  mutate(
    Vista_type = case_when(
      # case 1: if Targeted_Gene != " "，then is "E-P"
      Targeted_Gene != "" & !is.na(Targeted_Gene) ~ "E-P",
      
      # case 2 : if Targeted_Gene = " "，and Intersected_Anchor is "Left"
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Left" & 
        E_R == "TRUE" ~ "E-E",
      
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Left" & 
        E_R == "FALSE" ~ "E-U",
      
      # case 3 : if Targeted_Gene = " "，and Intersected_Anchor is "Right"
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Right" & 
        E_L == "TRUE" ~ "E-E",
      
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Right" & 
        E_L == "FALSE" ~ "E-U",
      
      # case 4 : if Targeted_Gene = " "，and Intersected_Anchor is "Both"
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Both" & 
        (E_L == "TRUE" | E_R == "TRUE") ~ "E-E",
      
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Both" & 
        E_L == "FALSE" & E_R == "FALSE" ~ "E-E",
      
      # other
      TRUE ~ "E-E"
    )
  )

cat("=== Check Vista_type Distribution ===\n")
print(table(Data1_overlapping$Vista_type, useNA = "ifany"))

# plot the distribution of interaction type
interaction_summary <- Data1_overlapping %>%
  group_by(Vista_type) %>%
  summarise(
    loop_count = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    percentage = round(loop_count / sum(loop_count) * 100, 2),
    proportion = round(loop_count / sum(loop_count), 4)
  ) %>%
  arrange(desc(loop_count))

print(interaction_summary)

total_loops <- sum(interaction_summary$loop_count)
cat("\n number of total loops:", total_loops, "\n")

pie_chart <- ggplot(interaction_summary, aes(x = "", y = loop_count, fill = Vista_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("E-P" = "#CC5D4B", 
                               "E-U" = "#B4CFD8", 
                               "E-E" = "#428AB7")) +
  labs(title = "VISTA FM enhancers related chromatin Loop distribution(n=334)",
       fill = "Vista_type") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  geom_text(aes(label = paste0(loop_count, "\n(", percentage, "%)")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black", fontface = "bold")


print(pie_chart)

ggsave("./3.VISTA_Enh_Gene/77FMenh_334chromatin_loop_interaction_types_pie_HEPM.pdf", pie_chart, width = 8, height = 6)
write.csv(Data1_overlapping, "./3.VISTA_Enh_Gene/77FMenh_3834chromatin_loop_interaction_types_HEPM.csv")


# --- Enhancer ID ---
FacialMesenchym_VistaID <- Data1_overlapping %>%
  mutate(temp = fm_ID) %>%
  separate_rows(temp, sep = ",") %>%
  dplyr::select(temp) %>%
  distinct() %>%
  dplyr::select(fm_Enh_ID = temp) %>% arrange((fm_Enh_ID))
write.csv(FacialMesenchym_VistaID, "77FacialMesenchymeVistaID.csv")

# --- find targeted genes regulated by FM enhancers in VISTA database ---
E_P_vista_Enh <- Data1_overlapping %>% dplyr::filter(Vista_type == "E-P")
Enh_ID <- E_P_vista_Enh %>% dplyr::select(fm_ID) %>% distinct()
Targeted_Gene_vista <- E_P_vista_Enh %>% 
  dplyr::select(Targeted_Gene) %>%
  filter(Targeted_Gene != "" & !is.na(Targeted_Gene)) %>%
  # separate by ","
  separate_rows(Targeted_Gene, sep = ",") %>%
  mutate(Targeted_Gene = str_trim(Targeted_Gene)) %>%
  filter(Targeted_Gene != "") %>%
  distinct(Targeted_Gene) %>%
  arrange(Targeted_Gene) 
write.csv(E_P_vista_Enh, "74EPLoops_35FMenh_41Genes.csv")
write.csv(Targeted_Gene_vista, "41Targeted_genes_35VistaFMenhancer.csv")


#-------------------------- MEPM chromatin loops associated with VISTA FM enhancers ------------
Data1 <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/MEPM_MicroC_6400bp_annotatedLoops_20251204.csv", header = T, sep = ",")
fm_enh <- read.delim("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/13_gkmSVM/VistaEnh/VISTA_ENH_250222/VISTA_positive_fm_mm10.bed", header = F)

fm_enh_gr <- GRanges(
  seqnames = fm_enh$V1,
  ranges = IRanges(start = fm_enh$V2, end = fm_enh$V3)
)

#
micro_gr1 <- GRanges(
  seqnames = Data1$chr1,
  ranges = IRanges(start = Data1$start1, end = Data1$end1)
)
mcols(micro_gr1)$original_row <- 1:nrow(Data1)
mcols(micro_gr1)$region <- "region1"
#
micro_gr2 <- GRanges(
  seqnames = Data1$chr2,
  ranges = IRanges(start = Data1$start2, end = Data1$end2)
)
mcols(micro_gr2)$original_row <- 1:nrow(Data1)
mcols(micro_gr2)$region <- "region2"

#  fm_enh overlapped with micro_gr1 and micro_gr2, respectively
overlaps1 <- findOverlaps(fm_enh_gr, micro_gr1)
overlaps2 <- findOverlaps(fm_enh_gr, micro_gr2)


overlap_info1 <- data.frame(
  fm_idx = queryHits(overlaps1),
  micro_row = mcols(micro_gr1)$original_row[subjectHits(overlaps1)],
  region = "region1"
)

overlap_info2 <- data.frame(
  fm_idx = queryHits(overlaps2),
  micro_row = mcols(micro_gr2)$original_row[subjectHits(overlaps2)],
  region = "region2"
)

# merge
all_overlaps <- rbind(overlap_info1, overlap_info2)
all_original_rows <- unique(c(overlap_info1$micro_row, overlap_info2$micro_row))

# 
Data1_overlapping <- Data1[all_original_rows, ]

# add fm_enh for corresponding micro_hg19_overlapping

Data1_overlapping$fm_chr <- ""
Data1_overlapping$fm_start <- ""
Data1_overlapping$fm_end <- ""
Data1_overlapping$fm_ID <- "" 
Data1_overlapping$Intersected_Anchor <- "" 
Data1_overlapping$Targeted_Gene <- ""      

for (i in 1:nrow(Data1_overlapping)) {
  row_idx <- all_original_rows[i]

  fm_indices_region1 <- overlap_info1$fm_idx[overlap_info1$micro_row == row_idx]
  fm_indices_region2 <- overlap_info2$fm_idx[overlap_info2$micro_row == row_idx]
  all_fm_indices <- unique(c(fm_indices_region1, fm_indices_region2))
  
  if (length(all_fm_indices) > 0) {
    Data1_overlapping$fm_chr[i] <- paste(fm_enh$V1[all_fm_indices], collapse = ",")
    Data1_overlapping$fm_start[i] <- paste(fm_enh$V2[all_fm_indices], collapse = ",")
    Data1_overlapping$fm_end[i] <- paste(fm_enh$V3[all_fm_indices], collapse = ",")
    Data1_overlapping$fm_ID[i] <- paste(fm_enh$V4[all_fm_indices], collapse = ",")
  }
  
  has_region1_overlap <- length(fm_indices_region1) > 0
  has_region2_overlap <- length(fm_indices_region2) > 0
  
  if (has_region1_overlap && has_region2_overlap) {
    Data1_overlapping$Intersected_Anchor[i] <- "Both"
    p_r_gene <- ifelse(is.na(Data1_overlapping$P_R_genes[i]) || Data1_overlapping$P_R_genes[i] == "", 
                       "", as.character(Data1_overlapping$P_R_genes[i]))
    p_l_gene <- ifelse(is.na(Data1_overlapping$P_L_genes[i]) || Data1_overlapping$P_L_genes[i] == "", 
                       "", as.character(Data1_overlapping$P_L_genes[i]))
    combined_genes <- paste(c(p_r_gene, p_l_gene)[c(p_r_gene, p_l_gene) != ""], collapse = ",")
    Data1_overlapping$Targeted_Gene[i] <- combined_genes
    
  } else if (has_region1_overlap) {
    Data1_overlapping$Intersected_Anchor[i] <- "Left"
    Data1_overlapping$Targeted_Gene[i] <- ifelse(is.na(Data1_overlapping$P_R_genes[i]) || Data1_overlapping$P_R_genes[i] == "", 
                                                 "", as.character(Data1_overlapping$P_R_genes[i]))
    
  } else if (has_region2_overlap) {
    Data1_overlapping$Intersected_Anchor[i] <- "Right"
    Data1_overlapping$Targeted_Gene[i] <- ifelse(is.na(Data1_overlapping$P_L_genes[i]) || Data1_overlapping$P_L_genes[i] == "", 
                                                 "", as.character(Data1_overlapping$P_L_genes[i]))
  }
}

print(paste("find", nrow(Data1_overlapping), " overlapped with fm_enh"))
print(paste("其中", length(unique(overlap_info1$micro_row)), "region1 overlapped with fm_enh"))
print(paste("其中", length(unique(overlap_info2$micro_row)), "region2 overlapped with fm_enh"))

cat("\n=== Intersected_Anchor Distribution ===\n")
print(table(Data1_overlapping$Intersected_Anchor, useNA = "ifany"))

cat("\n=== Targeted_Gene Summary ===\n")
cat("Rows with targeted genes:", sum(Data1_overlapping$Targeted_Gene != "", na.rm = TRUE), "\n")
cat("Rows without targeted genes:", sum(Data1_overlapping$Targeted_Gene == "", na.rm = TRUE), "\n")

cat("\n=== Sample Results ===\n")
print(head(Data1_overlapping[, c("chr1", "start1", "end1", "chr2", "start2", "end2", 
                                 "P_L_genes", "P_R_genes", "fm_ID", 
                                 "Intersected_Anchor", "Targeted_Gene")], 10))

# check the subtype distribution
Data1_overlapping <- Data1_overlapping %>%
  mutate(
    Vista_type = case_when(
      Targeted_Gene != "" & !is.na(Targeted_Gene) ~ "E-P",
        (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Left" & 
        E_R == "TRUE" ~ "E-E",
      
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Left" & 
        E_R == "FALSE" ~ "E-U",
      
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Right" & 
        E_L == "TRUE" ~ "E-E",
      
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Right" & 
        E_L == "FALSE" ~ "E-U",
      
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Both" & 
        (E_L == "TRUE" | E_R == "TRUE") ~ "E-E",
      
      (Targeted_Gene == "" | is.na(Targeted_Gene)) & 
        Intersected_Anchor == "Both" & 
        E_L == "FALSE" & E_R == "FALSE" ~ "E-E",
      
      # other
      TRUE ~ "E-E"
    )
  )

cat("=== Vista_type Distribution ===\n")
print(table(Data1_overlapping$Vista_type, useNA = "ifany"))

# plot the distribution of interaction type
interaction_summary <- Data1_overlapping %>%
  group_by(Vista_type) %>%
  summarise(
    loop_count = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    percentage = round(loop_count / sum(loop_count) * 100, 2),
    proportion = round(loop_count / sum(loop_count), 4)
  ) %>%
  arrange(desc(loop_count))

print(interaction_summary)

total_loops <- sum(interaction_summary$loop_count)
cat("\n total number of loop :", total_loops, "\n")

pie_chart <- ggplot(interaction_summary, aes(x = "", y = loop_count, fill = Vista_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("E-P" = "#CC5D4B", 
                               "E-U" = "#B4CFD8", 
                               "E-E" = "#428AB7")) +
  labs(title = "VISTA FM enhancers related chromatin Loop distribution(n=334)",
       fill = "Vista_type") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  geom_text(aes(label = paste0(loop_count, "\n(", percentage, "%)")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black", fontface = "bold")


print(pie_chart)

ggsave("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/VISTA_FM_Enh/MEPM_408chromatin_loop_interaction_types_pie.pdf", pie_chart, width = 8, height = 6)
write.csv(Data1_overlapping, "/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/VISTA_FM_Enh/100FMenh_408chromatin_loop_interaction_types_MEPM.csv")


# --- Enhancer ID ---
FacialMesenchym_VistaID <- Data1_overlapping %>%
  mutate(temp = fm_ID) %>%
  separate_rows(temp, sep = ",") %>%
  dplyr::select(temp) %>%
  distinct() %>%
  dplyr::select(fm_Enh_ID = temp) %>% arrange((fm_Enh_ID))
write.csv(FacialMesenchym_VistaID, "/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/VISTA_FM_Enh/100FacialMesenchymeVistaID.csv")

# --- find targeted genes regulated by FM enhancers in VISTA database ---
E_P_vista_Enh <- Data1_overlapping %>% dplyr::filter(Vista_type == "E-P")
Enh_ID <- E_P_vista_Enh %>% dplyr::select(fm_ID) %>% distinct()
Targeted_Gene_vista <- E_P_vista_Enh %>% 
  dplyr::select(Targeted_Gene) %>%
  filter(Targeted_Gene != "" & !is.na(Targeted_Gene)) %>%
  separate_rows(Targeted_Gene, sep = ",") %>%
  mutate(Targeted_Gene = str_trim(Targeted_Gene)) %>%
  filter(Targeted_Gene != "") %>%
  distinct(Targeted_Gene) %>%
  arrange(Targeted_Gene) 
write.csv(E_P_vista_Enh, "/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/VISTA_FM_Enh/92EPLoops_41FMenh_52Genes.csv")
write.csv(Targeted_Gene_vista, "/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/VISTA_FM_Enh/52Targeted_genes_41VistaFMenhancer_MEPM.csv")
# combined VISTA Genes in HEPM and MEPM

Targeted_Gene_vista_HEPM <- read.csv("41Targeted_genes_35VistaFMenhancer.csv", sep = ",", header = T)

library(gprofiler2)
orthologs <- gorth(
  query = Targeted_Gene_vista$Targeted_Gene,
  source_organism = "mmusculus", 
  target_organism = "hsapiens"
)
final_result <- Targeted_Gene_vista %>%
  left_join(orthologs, by = c("Targeted_Gene" = "input")) %>%
  dplyr::select(
    Original_Gene = Targeted_Gene, 
    HSA_Gene_Symbol = ortholog_name, 
    HSA_Ensembl_ID = ortholog_ensg
  ) %>%
  distinct() 

# combined targeted genes in HEPM and MEPM
df1 <- final_result %>% 
  dplyr::select(Gene = HSA_Gene_Symbol) %>% 
  filter(!is.na(Gene))

df2 <- Targeted_Gene_vista_HEPM %>% 
  dplyr::select(Gene = Targeted_Gene)

all_genes_unique <- bind_rows(df1, df2) %>%
  distinct(Gene) %>%  
  arrange(Gene)  

# chech CP Genes
Combined_CPGene <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/CP_combined_with_mouse_orthologs.csv", sep=",", header = T)
intersection <- intersect(all_genes_unique$Gene, Combined_CPGene$CP_Gene)


setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results/3.VISTA_Enh_Gene")
# save
save(list = ls(), file = "./3.VISTA_Enh_Gene/Fig2M.FMVEs_related_Loops_Genes.RData")

load("VISTA_EMenh_HEPM_microC_20251119.RData")
