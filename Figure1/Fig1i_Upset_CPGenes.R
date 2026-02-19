# -------- 185 OMIM Cleft Palate Genes --------------
setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene")
OMIM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/OMIM-Gene-Map-Search.tsv", header = T, sep = "\t")
head(OMIM, 5)
library("dplyr")
CP_Gene <- OMIM %>% dplyr::select(Approved.Symbol, Phenotype, Mouse.Gene..from.MGI.)

CP_causative_Gene <- CP_Gene %>% 
  dplyr::select(Approved.Symbol) %>% 
  filter(!is.na(Approved.Symbol) & Approved.Symbol != "") %>%
  unique()
OMIM_CP_genes <- CP_causative_Gene$Approved.Symbol
write.csv(CP_Gene, "NSCPO_Gene_list_OMIM.csv")
write.csv(CP_causative_Gene, "NSCPO_Gene_list_OMIM_symbol.csv")

# ----------296 Mouse CP Genes --------------
Mice_CP_Genes <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/Cleftpalate_Genes_MGI.csv", sep = ",", header = T)
Mice_CP_Genes_clean <- Mice_CP_Genes %>%
  pull(gene_symbol.human.) %>%
  .[!is.na(.) & . != ""] %>%
  unique()

# ----------- HEPM RNA-seq Expressed Genes --------
HEPM <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/1.initial_data/4_HEPM_wildtype_RNAseq/HEPM_RNAseq_sleuth_TPM/HEPMresults_TPM/HEPM_exp_transcript.csv", header = T, sep = ",")
# total 11533 genes (TPM >1)
HEPM_exp_genes <- HEPM %>%
  group_by(gene_symbol) %>%
  summarise(max_aveTPM = max(aveTPM, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(max_aveTPM > 1) %>%
  pull(gene_symbol) %>%
  .[!is.na(.) & . != ""]


# ------ E14 MEPM RNA-seq Expressed Genes ---------
setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/1.E145_Mes_AEs_Pros/RNA-seq")
library(biomaRt)
Mm <- read.csv("E14_TPM.tsv", sep = "\t", header = T)
Mm <- Mm %>%
  mutate(aveTPM = rowMeans(dplyr::select(., E14_AP_rep1_tpm, E14_AP_rep_2_tpm, E14_AP_rep3_tpm, 
                                         E14_PP_rep1_tpm, E14_PP_rep2_tpm, E14_PP_rep3_tpm)))
Mm$transcript_id_clean <- gsub("\\..*$", "", Mm$target_id)
options(timeout = 3000)
mouse_mart <- useMart("ensembl", 
                      dataset = "mmusculus_gene_ensembl",
                      host = "https://nov2020.archive.ensembl.org")
transcript_to_gene <- getBM(
  attributes = c("ensembl_transcript_id", "mgi_symbol", "external_gene_name"),
  filters = "ensembl_transcript_id",
  values = Mm$transcript_id_clean,
  mart = mouse_mart
)


#  transfer to symbol

Mm_with_symbol <- Mm %>%
  left_join(transcript_to_gene, 
            by = c("transcript_id_clean" = "ensembl_transcript_id"))


# aveTPM > 1 ; 13978 genes
Mm_exp_genes <- Mm_with_symbol %>%
  group_by(external_gene_name) %>%
  summarise(max_aveTPM = max(aveTPM, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(max_aveTPM > 1) %>%
  pull(external_gene_name)


# human orthologs gene
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mouse_to_human <- getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
                        filters = "external_gene_name",
                        values = Mm_exp_genes,
                        mart = mouse_mart)

# 12736 MEPM-human orthologs genes
MEPM_human_orthologs_genes_clean <- mouse_to_human %>%
  pull(hsapiens_homolog_associated_gene_name) %>%
  .[!is.na(.) & . != ""] %>%
  unique()

# --------- Upset plot -------

library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(grid)
set_list <- list(
  "Human CP" = OMIM_CP_genes,
  "HEPM Expressed" = HEPM_exp_genes,
  "MEPM Expressed" = MEPM_human_orthologs_genes_clean,
  "Mouse CP" = Mice_CP_Genes_clean
)
m <- make_comb_mat(set_list)
# color
set_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00")
names(set_colors) <- names(set_list)
comb_name(m)
existing_combs <- comb_name(m)
comb_info <- data.frame(
  index = 1:length(comb_size(m)),
  size = comb_size(m),
  degree = comb_degree(m),
  name = comb_name(m)
)
print(comb_info)
custom_comb_order <- match(valid_custom_order, existing_combs)
valid_custom_order <- c("1111", "1110", "0111", "1011","1010", "0011", "1101", "1100", "1001", "1000", "0001", "0110", "0010", "0100")
comb_indices <- match(valid_custom_order, comb_name(m))
comb_indices <- comb_indices[!is.na(comb_indices)]  

p <- UpSet(m, 
           top_annotation = upset_top_annotation(
             m, 
             add_numbers = TRUE,
             numbers_rot = 0,
             numbers_gp = gpar(fontsize = 12, fontface = "bold"),
             bar_width = 0.8,
             gp = gpar(fill = "steelblue", alpha = 0.8),
             height = unit(4, "cm")
           ),
           right_annotation = upset_right_annotation(
             m, 
             add_numbers = TRUE,
             numbers_rot = 0,
             numbers_gp = gpar(fontsize = 12, fontface = "bold"),
             bar_width = 0.8,
             gp = gpar(fill = set_colors, alpha = 0.8),
             width = unit(2, "cm")
           ),
           pt_size = unit(2, "mm"),
           lwd = 2,
           bg_col = "white",
           bg_pt_col = "lightgray",
           height = unit(3, "cm"), 
           width = unit(8, "cm"),
           # order
           set_order = order(set_size(m), decreasing = TRUE),
           comb_order = comb_indices,
           column_title = "Gene Set Intersection Analysis",
           column_title_gp = gpar(fontsize = 16, fontface = "bold"),
           row_title = "Gene Sets",
           row_title_gp = gpar(fontsize = 14, fontface = "bold"))

pdf("HEPM_MEPM_CP_Genes_upset.pdf", width = 10, height = 8)
p
dev.off()

# ------------ combined human and mouse CP genes --------
OMIM_CP_df <- data.frame(
  CP_Gene = OMIM_CP_genes,
  Source = "Human_CP_Genes",
  stringsAsFactors = FALSE
)

Mice_CP_df <- data.frame(
  CP_Gene = Mice_CP_Genes_clean,
  Source = "Mouse_CP_Genes",
  stringsAsFactors = FALSE
)

CP_combined <- full_join(OMIM_CP_df, Mice_CP_df, by = "CP_Gene") #CP Genes: 748
# common CP genes in human and mouse : 118
CP_both <- inner_join(OMIM_CP_df, Mice_CP_df, by = "CP_Gene")
nrow(CP_both) 
# convert to mouse gene
library(gprofiler2)
orthologs <- gorth(
  query = CP_combined$CP_Gene,
  source_organism = "hsapiens",
  target_organism = "mmusculus",
  mthreshold = Inf,
  filter_na = TRUE
)

head(orthologs)
str(orthologs)

orthologs_clean <- orthologs %>%
  dplyr::select(
    CP_Gene = input,
    mm_gene_symbol = ortholog_name,
    ortholog_ensg = ortholog_ensg
  ) %>%
  distinct()

CP_combined_with_mouse <- CP_combined %>%
  left_join(orthologs_clean, by = "CP_Gene")


# unmapped gene
unmapped <- CP_combined_with_mouse %>%
  filter(is.na(mm_gene_symbol)) %>%
  dplyr::select(CP_Gene, Source.x, Source.y)

print(unmapped)

write.csv(CP_combined_with_mouse, "CP_combined_with_mouse_orthologs.csv", row.names = FALSE)


# ------ 
library(openxlsx)
wb <- createWorkbook()
summary_df <- data.frame(
  Intersection_Code = comb_name(m),
  Size = comb_size(m),
  Degree = comb_degree(m)
)


set_names <- names(set_list)
for (i in 1:nrow(summary_df)) {
  pattern_bits <- strsplit(summary_df$Intersection_Code[i], "")[[1]]
  included_sets <- set_names[pattern_bits == "1"]
  summary_df$Sets[i] <- paste(included_sets, collapse = " & ")
}

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)


for (i in 1:length(comb_name(m))) {
  comb_pattern <- comb_name(m)[i]
  genes <- extract_comb(m, comb_pattern)
  
  if (length(genes) > 0) {

    sheet_name <- substr(comb_pattern, 1, 31)
    
    df <- data.frame(Gene = genes)
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)
  }
}


saveWorkbook(wb, "upset_all_intersections.xlsx", overwrite = TRUE)

save(list = ls(), file = "Fig2i_Upset_CPGenes.Rdata")
