HEPM_gene_Enhancers <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/HEPM/HEPM_gene_regulatory_layers.csv", header = T, sep = ",")

MEPM_gene_Enhancers <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/23.Enhancer_Hierary_Network/2.Process_data/1.BFS/MEPM/MEPM_gene_regulatory_layers.csv", header = T, sep = ",")


Combined_CPGene <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/CP_combined_with_mouse_orthologs.csv", sep=",", header = T)
Combined_CPGene_HEPM <- Combined_CPGene$CP_Gene %>% unique()
Combined_CPGene_MEPM <- Combined_CPGene$mm_gene_symbol %>% unique() %>% .[!is.na(.)]

HEPM_AEs <- read.delim2("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results/0.HEPM_MicroC/Enhancer/HEPM_NFR_enhancers_sorted_250928.bed", header = F)
MEPM_AEs <- read.delim2("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/1.E145_Mes_AEs_Pros/Final_results/MEPM_AEs_NFR_final.bed", header = F)

library(dplyr)

HEPM_CP_Enhancers <- HEPM_gene_Enhancers %>% filter(Gene %in% Combined_CPGene_HEPM)
MEPM_CP_Enhancers <- MEPM_gene_Enhancers %>% filter(Gene %in% Combined_CPGene_MEPM)
# ---------------------------------------- Select HEPM enhancers -------------------------------------

library(dplyr)
library(tidyr)

# select Tier1 Enhancers
tier1_ids <- HEPM_CP_Enhancers$Tier1_Enhancers %>%
  paste(collapse = ";") %>%  
  strsplit(";") %>%           
  unlist() %>%                
  unique()                    

cat("Tier1 Enhancer numbers:", length(tier1_ids), "\n") # 1766

# select Tier2 Enhancers
tier2_ids <- HEPM_CP_Enhancers$Tier2_Enhancers %>%
  paste(collapse = ";") %>%
  strsplit(";") %>%
  unlist() %>%
  unique()

cat("Tier2 Enhancer numbers:", length(tier2_ids), "\n") # 1534

# Enhancer ID both in Tier1 and Tier2
overlap_ids <- intersect(tier1_ids, tier2_ids)
cat("overlapped Enhancer number:", length(overlap_ids), "\n")

# Tier2 without overlapped ID
tier2_ids_unique <- setdiff(tier2_ids, tier1_ids)

cat("  Tier1 Enhancer number :", length(tier1_ids), "\n") # 1766
cat("  Tier2 Enhancer number (without overlapped):", length(tier2_ids_unique), "\n") # 1499
# obtained coordinate from HEPM_AEs
# add colnames
colnames(HEPM_AEs) <- c("chr", "start", "end", "Enhancer_ID")

# Tier1 gene coordinate
tier1_coords <- HEPM_AEs %>%
  filter(Enhancer_ID %in% tier1_ids) %>%
  arrange(chr, start)

cat("\nTier1 Enhancer:", nrow(tier1_coords), "\n")

# Tier2 gene coordinate
tier2_coords <- HEPM_AEs %>%
  filter(Enhancer_ID %in% tier2_ids_unique) %>%
  arrange(chr, start)

cat("Tier2 Enhancers:", nrow(tier2_coords), "\n")

# Tier1: rename as  T1_Enhancer_1, T1_Enhancer_2, ...
tier1_coords <- tier1_coords %>%
  mutate(
    Original_ID = Enhancer_ID,
    Enhancer_ID = paste0("T1_", Original_ID)
  )

# Tier2: rename as T2_Enhancer_1, T2_Enhancer_2, ...
tier2_coords <- tier2_coords %>%
  mutate(
    Original_ID = Enhancer_ID,
    Enhancer_ID = paste0("T2_", Original_ID)
  )

# combined Tier1 and Tier2
combined_enhancers <- rbind(
  tier1_coords %>% mutate(Tier = "Tier1"),
  tier2_coords %>% mutate(Tier = "Tier2")
)

# BED format: chr, start, end, name
bed_data <- combined_enhancers %>%
  dplyr::select(chr, start, end, Enhancer_ID, Tier, Original_ID)

print(head(bed_data, 5))

# integrated data
write.table(
  bed_data,
  file = "HEPM_CP_Enhancers_Combined.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

# standard format
bed_standard <- bed_data %>%
  dplyr::select(chr, start, end, Enhancer_ID)

write.table(
  bed_standard,
  file = "HEPM_CP_Enhancers_Standard.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE  
)
