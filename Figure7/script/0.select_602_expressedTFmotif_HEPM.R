
# ----- homoer motif analysis -----
library(TFBSTools)
library(tidyverse)
pfm_list <- readJASPARMatrix("./TFs_jaspar/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar/All_TFs1019.jaspar", matrixClass = "PFM")
HEPM_rnaseq <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/1.initial_data/4_HEPM_wildtype_RNAseq/HEPM_RNAseq_sleuth_TPM/HEPM_TPM_FINAL.CSV", header = T)
HEPM_rnaseq <- HEPM_rnaseq %>% mutate(Mean_TPM = (HEPM1_TPM + HEPM2_TPM + HEPM3_TPM) / 3)
# TPM >1
expressed_genes <- HEPM_rnaseq %>%
  filter(Mean_TPM > 1)
# -- select motif for TFs ---
tf_info <- data.frame(
  Motif_ID = sapply(pfm_list, ID),
  TF_Name = sapply(pfm_list, name),
  stringsAsFactors = FALSE
)
# upper 
tf_info$TF_upper <- toupper(tf_info$TF_Name)
tf_info_expanded <- tf_info %>%
  rowwise() %>%
  mutate(
    TF_variants = list(toupper(trimws(unlist(strsplit(TF_Name, "::")))))
  ) %>%
  unnest(TF_variants) %>%
  ungroup()
# matching
expressed_genes$Gene_upper <- toupper(expressed_genes$gene_symbol)

matched_tfs <- tf_info_expanded %>%
  filter(TF_variants %in% expressed_genes$Gene_upper) %>%
  left_join(expressed_genes, by = c("TF_variants" = "Gene_upper"))

# summary
expressed_tfs_summary <- matched_tfs %>%
  group_by(Motif_ID, TF_Name) %>%
  summarise(
    Matched_genes = paste(unique(TF_variants), collapse = ", "),
    Mean_TPM = mean(Mean_TPM),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_TPM))


write.csv(expressed_tfs_summary, "HEPM_expressed_TFs.csv", row.names = FALSE)

expressed_motif_ids <- expressed_tfs_summary$Motif_ID
expressed_pfm_list <- pfm_list[sapply(pfm_list, ID) %in% expressed_motif_ids]
# -- output jaspar
write_jaspar_format <- function(pfm_list, output_file) {
  con <- file(output_file, "w")
  for (i in 1:length(pfm_list)) {
    pfm <- pfm_list[[i]]
    # obtained motif data
    motif_id <- ID(pfm)
    motif_name <- name(pfm)
    matrix_data <- Matrix(pfm)
    writeLines(paste0(">", motif_id, "\t", motif_name), con)
    
    # JASPAR format
    for (base in c("A", "C", "G", "T")) {
      values <- matrix_data[base, ]
      # format: "A  [  123  57  0  ... ]"
      line <- paste0(base, "  [ ", 
                     paste(sprintf("%6d", values), collapse = " "), 
                     " ]")
      writeLines(line, con)
    }
  }

  close(con)
  
  cat("write", length(pfm_list), "motifs to", output_file, "\n")
}

write_jaspar_format(expressed_pfm_list, "HEPM_expressed_TFs.jaspar") # 602个motifs




