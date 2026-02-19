
# ---  gene expression level in HEPM cells ---
setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results")
getwd()
library(dplyr)
library(tidyverse)
gene_list <- read.csv("./2.SNP_candidateGenes/SNP_targetedGenes.csv", header = F, sep = ",")
gene_list <- gene_list %>%
  rename_with(~ "gene_list", 1)
tpm_file_path <- "/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/1.initial_data/4_HEPM_wildtype_RNAseq/HEPM_RNAseq_sleuth_TPM/HEPM_TPM_FINAL.CSV"
output_dir <- "./2.SNP_candidateGenes"


result <- plot_gene_expression(
  gene_list = gene_list$gene_list,
  tpm_file_path = tpm_file_path,
  output_dir = output_dir,
  csv_filename = "Exp_CandidateGenes_HEPM.csv",
  pdf_filename = "Exp_CandidateGenes_HEPM.pdf",
  plot_title = "Gene Expression of candidate Genes"
)

plot_gene_expression <- function(gene_list, 
                                 tpm_file_path, 
                                 output_dir = NULL,
                                 csv_filename = "Exp_Candidates_HEPM.csv",
                                 pdf_filename = "exp_Candidates_HEPM.pdf",
                                 pdf_width = 6,
                                 pdf_height = 4,
                                 plot_title = "Gene Expression of candidate Genes") {
  
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  csv_output_path <- file.path(output_dir, csv_filename)
  pdf_output_path <- file.path(output_dir, pdf_filename)
  
  require(dplyr)
  require(tidyverse)
  
  HEPM_TPM_modified_2 <- read.csv(tpm_file_path)
  
  HEPM_TFs_TPM <- HEPM_TPM_modified_2 %>% filter(gene_symbol %in% gene_list)
  
  if (nrow(HEPM_TFs_TPM) == 0) {
    stop("not found genes in target list")
  }
  
  HEPM_TFs_TPM_long <- HEPM_TFs_TPM %>%
    pivot_longer(cols = c(HEPM1_TPM, HEPM2_TPM, HEPM3_TPM),
                 names_to = "Sample",
                 values_to = "TPM") %>%
    mutate(Sample = factor(Sample, levels = c("HEPM1_TPM", "HEPM2_TPM", "HEPM3_TPM")))
  
  HEPM_TFs_TPM_avg <- HEPM_TFs_TPM %>%
    rowwise() %>%
    mutate(Mean_TPM = mean(c(HEPM1_TPM, HEPM2_TPM, HEPM3_TPM))) %>%
    ungroup() %>%
    arrange(desc(Mean_TPM))
  
  write.csv(HEPM_TFs_TPM, csv_output_path, row.names = FALSE)
  
  gene_order <- HEPM_TFs_TPM_avg$gene_symbol
  
  HEPM_TFs_TPM_avg$gene_symbol <- factor(HEPM_TFs_TPM_avg$gene_symbol, levels = gene_order)
  
  HEPM_TFs_TPM_long$gene_symbol <- factor(HEPM_TFs_TPM_long$gene_symbol, levels = gene_order)
  
  plot.format <- theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 7),
    axis.title = element_text(color = "black", size = 7),
    plot.title = element_text(color = "black", size = 7),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(color = "black", size = 7),
    legend.title = element_text(color = "black", size = 7)
  )
  
  p <- ggplot() +
    geom_bar(data = HEPM_TFs_TPM_avg, 
             aes(x = gene_symbol, y = Mean_TPM), 
             stat = "identity", 
             fill = "steelblue", 
             alpha = 0.7,
             width = 0.7) +
    geom_point(data = HEPM_TFs_TPM_long, 
               aes(x = gene_symbol, y = TPM, color = Sample),
               position = position_dodge(width = 0.5),
               size = 1) +
    geom_text(data = HEPM_TFs_TPM_avg,
              aes(x = gene_symbol, y = Mean_TPM, label = sprintf("%.1f", Mean_TPM)),
              vjust = -0.5,
              size = 2.5) +
    labs(x = "Gene Symbol", 
         y = "Gene Expression Level (TPM)",
         title = plot_title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom") +
    scale_color_manual(values = c("HEPM1_TPM" = "black", 
                                  "HEPM2_TPM" = "black", 
                                  "HEPM3_TPM" = "black")) +
    plot.format
  
  ggsave(pdf_output_path, plot = p, width = pdf_width, height = pdf_height)
  
  print(p)
  
  return(list(data = HEPM_TFs_TPM, plot = p))
}

# --- Gene expression level in secondary palate of mice ---

lapply(c("sctransform", "harmony", "data.table", "scales", "plotly", "ggplot2", "ggpointdensity", "monocle3", "magrittr", 
         "dplyr", "tidyverse", "cowplot", "clustree", "gridExtra", "RColorBrewer", "Seurat"), library, character.only = TRUE); 
set.seed(123); theme_set(theme_bw())


sp_combined <- readRDS("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/5_Expression_pattern_of_candidategenes/3_mice/GSE205449_RAW_23NC/scATAC_scRNA/combined/SP_combined_celltype_20241022.rds")
library(gprofiler2)

orthology_result <- gorth(
  query = gene_list$gene_list,
  source_organism = "hsapiens",      
  target_organism = "mmusculus",     
  mthreshold = 1,                    
  filter_na = TRUE                   
)
gene_mm <- orthology_result %>%
  select(
    Human_Gene = input,
    Mouse_Gene = ortholog_name,
    Human_Ensembl = input_ensg,
    Mouse_Ensembl = ortholog_ensg
  )


genes_to_check <- c(
  "Acvr2a",    # ACVR2A
  "Ak1",       # AK1
  "Angptl2",   # ANGPTL2
  "A130010J15Rik",     # C1orf74
  "Cdk9",      # CDK9
  "G0s2",      # G0S2
  "Hhat",      # HHAT
  "Lamb3",     # LAMB3
  "Niban2",    # NIBAN2
  "Sec22c",    # SEC22C
  "Sertad4",   # SERTAD4
  "Slc2a8",    # SLC2A8
  "Snrk",      # SNRK
  "Tmem254",   # TMEM254
  "Tor2a",     # TOR2A
  "Traf3ip3",  # TRAF3IP3
  "Utp25",     # UTP25
  "Zbtb34",    # ZBTB34
  "Zbtb43",    # ZBTB43
  "Zbtb47",     # ZBTB47
  
)



Idents(sp_combined) <- sp_combined$celltype
DotPlot(sp_combined, features=unique(gene_mm$Mouse_Gene))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45,hjust = 0.5, vjust = 0.5))+
  labs(x=NULL, y=NULL)+ guides(size = guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066', '#336699','#66CC66','#FFCC33'))
ggsave("./2.SNP_candidateGenes/mm10_SP_dotplot_NSCPO_CandidateGenes.pdf", width = 30, height =15,units = "cm")

p1 <- FeaturePlot(sp_combined, features = "Hhat",label = T, pt.size  = 0.1, min.cutoff = 0.1, 
                  cols = c("lightgrey", "red"))+  theme(panel.border = element_rect(fill=NA,color="black", 
                                                                                    linewidth = 1, linetype="solid"), 
                                                        legend.position = "right")+  labs(title = "Hhat")

p1

dir.create("./2.SNP_candidateGenes/cleft_palate_gene", showWarnings = FALSE)

for(gene in gene_mm$Mouse_Gene) {
  tryCatch({
    p <- FeaturePlot(sp_combined, features = gene, label = TRUE, pt.size = 0.1, 
                     min.cutoff = 0.1, cols = c("lightgrey", "red")) +
      theme(panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
            legend.position = "right") +
      labs(title = gene)
    
    ggsave(filename = paste0("./2.SNP_candidateGenes/cleft_palate_gene/", gene, ".pdf"), 
           plot = p, width = 8, height = 7, device = "pdf")
    
    cat("save:", gene, "\n")
  }, error = function(e) {
    cat("skipped genes:", gene, "wrong message:", conditionMessage(e), "\n")
  })
}





