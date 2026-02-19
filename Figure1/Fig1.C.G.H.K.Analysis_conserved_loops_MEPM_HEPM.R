setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/processed_results/3.compare_loop/Total_conservedLoop")
getwd()
list.files()
library(readr)
library(dplyr)

file_list <- list.files(pattern = "^HEPM_.*_loop_conserved_matched_pairs_details\\.csv$")
merged_data <- bind_rows(lapply(file_list, read_csv))
write_csv(merged_data, "HEPM_merged_all.csv")
print(dim(merged_data))
HEPM_uniq <- merged_data %>% 
  dplyr::select(HEPM_loop_id) %>%
  unique()
# HEPM_uniq : 100095
mm10_uniq <- merged_data %>%
  dplyr::select(mm10_loop_id) %>%
  unique()
# mm10_uniq: 89465
HEPM_only: 151979
common: 100095
mm10_only: 154819
library(VennDiagram)

if (dev.cur() != 1) {
  dev.off()
}

pdf("Fig.2C.Conserved_loops_HEPM_MEPM.pdf", width=4, height=4)
venn.plot <- draw.pairwise.venn(
  area1 = 252074,
  area2 = 254914,
  cross.area = 100095,
  category = c("HEPM", "mm10"),
  fill = c("skyblue", "lightgreen"),
  lty = "blank",                   
  cex = 1.5,                       
  cat.cex = 1.8,                   
  cat.pos = c(-20, 20),            
  cat.dist = 0.05,                 
  margin = 0.05                    
)
dev.off()


# Figure 2g/h
# --- load HEPM/MEPM loops ---
HEPM_loop <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results/1.Loop_annotation/HEPM_microC_6400bp_annotate_loops_20251204.csv", header = T, sep = ",")
E145_loop <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/MEPM_MicroC_6400bp_annotatedLoops_20251204.csv", header = T, sep = ",")
# start1 > start2 
sum(merged_data$mm10_start1 > merged_data$mm10_start2)
merged_data <- merged_data %>%
  mutate(
    temp_chr1 = mm10_chr1,
    temp_start1 = mm10_start1,
    temp_end1 = mm10_end1,
    temp_chr2 = mm10_chr2,
    temp_start2 = mm10_start2,
    temp_end2 = mm10_end2,
    # make sure start1 < start2
    mm10_chr1 = ifelse(temp_start1 > temp_start2, temp_chr2, temp_chr1),
    mm10_start1 = ifelse(temp_start1 > temp_start2, temp_start2, temp_start1),
    mm10_end1 = ifelse(temp_start1 > temp_start2, temp_end2, temp_end1),
    mm10_chr2 = ifelse(temp_start1 > temp_start2, temp_chr1, temp_chr2),
    mm10_start2 = ifelse(temp_start1 > temp_start2, temp_start1, temp_start2),
    mm10_end2 = ifelse(temp_start1 > temp_start2, temp_end1, temp_end2)
  ) %>%
  dplyr::select(-starts_with("temp_"))  
# return 0
sum(merged_data$mm10_start1 > merged_data$mm10_start2)



HEPM_loop_renamed <- HEPM_loop %>%
  dplyr::select(-X) %>%  
  rename_with(~paste0("HEPM_", .), 
              .cols = -c(chr1, start1, end1, chr2, start2, end2)) 

merged_result <- merged_data %>%
  left_join(HEPM_loop_renamed, 
            by = c("HEPM_chr1" = "chr1", 
                   "HEPM_start1" = "start1", 
                   "HEPM_end1" = "end1",
                   "HEPM_chr2" = "chr2", 
                   "HEPM_start2" = "start2", 
                   "HEPM_end2" = "end2"))


MEPM_loop_renamed <- E145_loop %>%
  dplyr::select(-X) %>%  
  rename_with(~paste0("MEPM_", .), 
              .cols = -c(chr1, start1, end1, chr2, start2, end2)) 

merged_result_V2 <- merged_result %>%
  left_join(MEPM_loop_renamed, 
            by = c("mm10_chr1" = "chr1", 
                   "mm10_start1" = "start1", 
                   "mm10_end1" = "end1",
                   "mm10_chr2" = "chr2", 
                   "mm10_start2" = "start2", 
                   "mm10_end2" = "end2"))
# --- filter conserved Enhancer-Promtoer --------
conserved_EP <- merged_result_V2 %>% 
  filter(HEPM_interaction_type == "Enhancer-Promoter", 
         MEPM_interaction_type == "Enhancer-Promoter")
conserve_HEPM_EP <- conserved_EP %>%
    dplyr::select(HEPM_loop_id) %>%
     distinct()
conserve_MEPM_EP <- conserved_EP %>%
    dplyr::select(mm10_loop_id) %>% distinct()

library(VennDiagram)
pdf("Fig.2g.Conserved_EPloops_HEPM_MEPM.pdf", width=4, height=4)

venn.plot <- draw.pairwise.venn(
  area1 = 28626,
  area2 = 23994,
  cross.area = 5533,
  category = c("HEPM", "MEPM"),
  fill = c("skyblue", "#EC8C12"), 
  lty = "blank",                   
  cex = 1.5,                       
  cat.cex = 1.8,                   
  cat.pos = c(-20, 20),             
  cat.dist = 0.05,                 
  margin = 0.05                    
)
dev.off()

# --- distance length distribution of HEPM conserved E-P loops ---
conserved_EP_distance <- conserved_EP %>%
  dplyr::select(HEPM_loop_id, HEPM_Distance_group) %>%
  distinct()

#
group_counts <- conserved_EP_distance %>%
  group_by(HEPM_Distance_group) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(Distance_group = factor(HEPM_Distance_group, 
                                 levels = c(">1000kb", "500~1000kb", "100~500kb", "<100kb")))


print(group_counts)

colors <- c("<100kb" = "#663399",     
            "100~500kb" = "#339999",  
            "500~1000kb" = "#66CC66",   
            ">1000kb" = "#CCCC33")     

# Fig2H distance for conserved_EPloops in HEPM and MEPM
p <- ggplot(group_counts, aes(x = "", y = percentage, fill = Distance_group)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.01) +
  scale_fill_manual(values = colors,
                    breaks = c("<100kb", "100~500kb", "500~1000kb", ">1000kb")) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12, color = "black", face = "bold"),  
    axis.text.x = element_text(size = 10, color = "black"),   
    axis.ticks.x = element_line(color = "black"),            
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line.x = element_line(color = "black", linewidth = 0.5)  
  ) +
  labs(
    x = "Percentage (%)", 
    fill = "Interaction Distance"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20), 
    expand = c(0, 0)
  ) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5),
            color = "white", 
            size = 5, 
            fontface = "bold")

print(p)
ggsave("Fig.2H.Conserved_HEPM_EP_loops_Distribution_of_LoopDistance.pdf", p, width = 25, height = 5, dpi = 300)


# ------------------- select genes -------------------
genes_from_rule1 <- conserved_EP %>%
  filter(HEPM_E_L == "TRUE") %>%
  pull(HEPM_P_R_genes)

genes_from_rule2 <- conserved_EP %>%
  filter(HEPM_E_R == "TRUE") %>%
  pull(HEPM_P_L_genes)


all_genes_raw <- c(genes_from_rule1, genes_from_rule2)


conserved_genes <- all_genes_raw %>%
  str_split(",") %>%
  unlist() %>%
  .[. != ""] %>%
  na.omit() %>%
  unique() %>%
  sort()

print(conserved_genes)
conserved_genes <- as.data.frame(conserved_genes)
write.csv(conserved_EP, "Conserved_EPloops_HEPM_MEPM.csv")
write.csv(conserved_genes, "conserved_genes_HEPM_MEPM.csv")

# --- mice genes ---
genes_from_MEPM1 <- conserved_EP %>%
  filter(MEPM_E_L == "TRUE") %>%
  pull(MEPM_P_R_genes)

genes_from_MEPM2 <- conserved_EP %>%
  filter(MEPM_E_R == "TRUE") %>%
  pull(MEPM_P_L_genes)


all_genes_MEPM <- c(genes_from_MEPM1, genes_from_MEPM2)

# 
conserved_genes_MEPM <- all_genes_MEPM %>%
  str_split(",") %>%
  unlist() %>%
  .[. != ""] %>%
  na.omit() %>%
  unique() %>%
  sort()
conserved_genes_MEPM <- as.data.frame(conserved_genes_MEPM)
write.csv(conserved_genes_MEPM, "conserved_genes_MEPM_Genes.csv")

##### ------ Fig. 2k conserved chromatin loops associated with CP genes --------
# load combined cp genes
Combined_CPGene <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/CP_combined_with_mouse_orthologs.csv", sep=",", header = T)
# load detailed_CPgene_analysis
detailed_CPgene_analysis <- function(data1, CP_gene, gene_col = "hsa_symbol") {
  cp_genes <- CP_gene[[gene_col]]
  cp_genes <- cp_genes[!is.na(cp_genes) & cp_genes != ""]
  cp_genes <- unique(cp_genes)
  pattern <- paste0("\\b(", paste(cp_genes, collapse = "|"), ")\\b")
  filtered_data <- data1 %>%
    filter(
      (str_detect(P_L_genes, pattern, negate = FALSE) & !is.na(P_L_genes)) | 
        (str_detect(P_R_genes, pattern, negate = FALSE) & !is.na(P_R_genes))
    )
  detailed_data <- filtered_data %>%
    mutate(
      matched_P_L = str_extract_all(P_L_genes, pattern),
      matched_P_R = str_extract_all(P_R_genes, pattern)
    )
  all_genes <- c(
    unlist(detailed_data$matched_P_L),
    unlist(detailed_data$matched_P_R)
  )
  all_genes <- all_genes[all_genes != ""]
  gene_counts <- as.data.frame(table(all_genes), stringsAsFactors = FALSE)
  colnames(gene_counts) <- c("gene", "count")
  gene_counts <- gene_counts[order(gene_counts$count, decreasing = TRUE), ]
  rownames(gene_counts) <- NULL
  return(list(
    filtered_data = filtered_data,
    gene_frequency = gene_counts,
    summary = list(
      total_loops = nrow(filtered_data),
      unique_genes = nrow(gene_counts),
      total_CP_genes = length(cp_genes)  
    )
  ))
}
# load extract_enhancer_ids
extract_enhancer_ids <- function(filtered_df, CP_gene, gene_col = "gene_symbol") {
  # Step 1: Subset Enhancer-Promoter interactions
  ep_interactions <- filtered_df %>%
    filter(interaction_type == "Enhancer-Promoter")
  
  cat("Total Enhancer-Promoter interactions:", nrow(ep_interactions), "\n")
  
  if(nrow(ep_interactions) == 0) {
    cat("No Enhancer-Promoter interactions found.\n")
    return(list(
      summary = list(
        total_ep_interactions = 0,
        ep_interactions_with_omim = 0,
        interactions_with_enhancers = 0,
        total_enhancer_instances = 0,
        unique_enhancer_count = 0,
        unique_enhancer_ids = character(0)
      ),
      detailed_results = data.frame(),
      unique_enhancers = character(0),
      E_P_Loops_OMIMGenes = data.frame()
    ))
  }
  
  # Get CP gene symbols 
  cp_genes <- CP_gene[[gene_col]]
  cp_genes <- cp_genes[!is.na(cp_genes) & cp_genes != ""]
  cp_genes <- unique(cp_genes)
  if(length(cp_genes) == 0) {
    stop(paste0("No valid gene symbols found in column '", gene_col, "'"))
  }
  
  cat("Valid CP genes loaded:", length(cp_genes), "\n")
  
  # construct pattern
  omim_pattern <- paste0("\\b(", paste(cp_genes, collapse = "|"), ")\\b")
  
  # Initialize result vectors
  all_enhancer_ids <- character(0)
  detailed_results <- data.frame()
  
  # Step 2: Process each row based on the criteria
  for(i in 1:nrow(ep_interactions)) {
    row <- ep_interactions[i, ]
    enhancer_ids_from_row <- character(0)
    extraction_reason <- character(0)
    # Convert logical values to character for comparison
    p_l_status <- as.character(row$P_L)
    p_r_status <- as.character(row$P_R)
    # Extract CP genes from P_L_genes and P_R_genes
    omim_cp_genes <- character(0)
    
    # Check P_L_genes for OMIM genes
    if(!is.na(row$P_L_genes) && row$P_L_genes != "") {
      p_l_matches <- str_extract_all(as.character(row$P_L_genes), omim_pattern)[[1]]
      omim_cp_genes <- c(omim_cp_genes, p_l_matches)
    }
    
    # Check P_R_genes for CP genes
    if(!is.na(row$P_R_genes) && row$P_R_genes != "") {
      p_r_matches <- str_extract_all(as.character(row$P_R_genes), omim_pattern)[[1]]
      omim_cp_genes <- c(omim_cp_genes, p_r_matches)
    }
    
    # Remove duplicates and empty strings
    omim_cp_genes <- unique(omim_cp_genes[omim_cp_genes != ""])
    omim_cp_gene_string <- paste(omim_cp_genes, collapse = ",")
    
    # Case 1: P_L is True, P_R is False -> extract from E_R_IDs
    if(p_l_status == "TRUE" && p_r_status == "FALSE") {
      if(!is.na(row$E_R_IDs) && row$E_R_IDs != "") {
        enhancer_ids_from_row <- unlist(strsplit(as.character(row$E_R_IDs), ","))
        extraction_reason <- "P_L=TRUE, P_R=FALSE -> E_R_IDs"
      }
    }
    
    # Case 2: P_R is True, P_L is False -> extract from E_L_IDs
    else if(p_r_status == "TRUE" && p_l_status == "FALSE") {
      if(!is.na(row$E_L_IDs) && row$E_L_IDs != "") {
        enhancer_ids_from_row <- unlist(strsplit(as.character(row$E_L_IDs), ","))
        extraction_reason <- "P_R=TRUE, P_L=FALSE -> E_L_IDs"
      }
    }
    
    # Case 3: Both P_L and P_R are True -> check CP genes
    else if(p_l_status == "TRUE" && p_r_status == "TRUE") {
      p_l_has_omim <- !is.na(row$P_L_genes) && row$P_L_genes != "" && 
        str_detect(as.character(row$P_L_genes), omim_pattern)
      p_r_has_omim <- !is.na(row$P_R_genes) && row$P_R_genes != "" && 
        str_detect(as.character(row$P_R_genes), omim_pattern)
      
      # Both contain CP genes -> extract from both E_L_IDs and E_R_IDs
      if(p_l_has_omim && p_r_has_omim) {
        e_l_ids <- character(0)
        e_r_ids <- character(0)
        if(!is.na(row$E_L_IDs) && row$E_L_IDs != "") {
          e_l_ids <- unlist(strsplit(as.character(row$E_L_IDs), ","))
        }
        if(!is.na(row$E_R_IDs) && row$E_R_IDs != "") {
          e_r_ids <- unlist(strsplit(as.character(row$E_R_IDs), ","))
        }
        
        enhancer_ids_from_row <- c(e_l_ids, e_r_ids)
        extraction_reason <- "Both P_L and P_R contain OMIM -> E_L_IDs + E_R_IDs"
      }
      
      # Only P_L contains CP genes -> extract from E_R_IDs
      else if(p_l_has_omim && !p_r_has_omim) {
        if(!is.na(row$E_R_IDs) && row$E_R_IDs != "") {
          enhancer_ids_from_row <- unlist(strsplit(as.character(row$E_R_IDs), ","))
          extraction_reason <- "P_L contains OMIM -> E_R_IDs"
        }
      }
      
      # Only P_R contains CP genes -> extract from E_L_IDs
      else if(!p_l_has_omim && p_r_has_omim) {
        if(!is.na(row$E_L_IDs) && row$E_L_IDs != "") {
          enhancer_ids_from_row <- unlist(strsplit(as.character(row$E_L_IDs), ","))
          extraction_reason <- "P_R contains OMIM -> E_L_IDs"
        }
      }
    }
    
    # Clean and add enhancer IDs
    if(length(enhancer_ids_from_row) > 0) {
      enhancer_ids_from_row <- str_trim(enhancer_ids_from_row)
      enhancer_ids_from_row <- enhancer_ids_from_row[enhancer_ids_from_row != ""]
      all_enhancer_ids <- c(all_enhancer_ids, enhancer_ids_from_row)
      
      # Store detailed results
      for(enhancer_id in enhancer_ids_from_row) {
        detailed_results <- rbind(detailed_results, data.frame(
          row_index = i,
          chr1 = row$chr1,
          start1 = row$start1,
          end1 = row$end1,
          chr2 = row$chr2,
          start2 = row$start2,
          end2 = row$end2,
          P_L_genes = row$P_L_genes,
          P_R_genes = row$P_R_genes,
          P_L = row$P_L,
          P_R = row$P_R,
          OMIM_CP_Gene = omim_cp_gene_string,
          enhancer_id = enhancer_id,
          extraction_reason = extraction_reason,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Step 3: Sort, distinct, and summarize
  unique_enhancer_ids <- sort(unique(all_enhancer_ids))
  # Create E_P_Loops_OMIMGenes dataframe - all EP interactions containing OMIM genes
  E_P_Loops_OMIMGenes <- ep_interactions %>%
    rowwise() %>%
    mutate(
      has_omim_genes = {
        p_l_match <- FALSE
        p_r_match <- FALSE
        # Check P_L_genes
        if(!is.na(P_L_genes) && P_L_genes != "") {
          p_l_match <- str_detect(as.character(P_L_genes), omim_pattern)
        }
        
        # Check P_R_genes
        if(!is.na(P_R_genes) && P_R_genes != "") {
          p_r_match <- str_detect(as.character(P_R_genes), omim_pattern)
        }
        
        p_l_match || p_r_match
      }
    ) %>%
    ungroup() %>%
    filter(has_omim_genes) %>%
    dplyr::select(-has_omim_genes)
  
  # Create summary
  summary_stats <- list(
    total_ep_interactions = nrow(ep_interactions),
    ep_interactions_with_omim = nrow(E_P_Loops_OMIMGenes),
    interactions_with_enhancers = length(unique(detailed_results$row_index)),
    total_enhancer_instances = length(all_enhancer_ids),
    unique_enhancer_count = length(unique_enhancer_ids),
    unique_enhancer_ids = unique_enhancer_ids
  )
  
  cat("\n=== Summary ===\n")
  cat("Total E-P interactions:", summary_stats$total_ep_interactions, "\n")
  cat("E-P interactions with OMIM genes:", summary_stats$ep_interactions_with_omim, "\n")
  cat("Interactions with extracted enhancers:", summary_stats$interactions_with_enhancers, "\n")
  cat("Total enhancer instances:", summary_stats$total_enhancer_instances, "\n")
  cat("Unique enhancers:", summary_stats$unique_enhancer_count, "\n")
  
  return(list(
    summary = summary_stats,
    detailed_results = detailed_results,
    unique_enhancers = unique_enhancer_ids,
    E_P_Loops_OMIMGenes = E_P_Loops_OMIMGenes
  ))
}

# start analysis: conserved chromatin loops in HEPM associated with CP genes
Conserved_HEPM_loop <- merged_result_V2 %>%
  dplyr::select(contains("HEPM")) %>% distinct()
# 
colnames(Conserved_HEPM_loop) <- gsub("HEPM_", "", colnames(Conserved_HEPM_loop))

conserved_CP_loops <- detailed_CPgene_analysis(Conserved_HEPM_loop, Combined_CPGene, gene_col = "CP_Gene")

if(nrow(conserved_CP_loops$gene_frequency) > 0) {
  p <- ggplot(conserved_CP_loops$gene_frequency, aes(x = reorder(gene, count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "CP genes frequency in conserved loops",
         x = "Conbimed CP Genes", y = "Number") +
    theme_minimal()
  
  print(p)
  
  # 保存图片
  ggsave("CP_genes_in_conservedchromatin_loops.pdf", plot = p, width = 10, height = 20)
}
# plot pie for intersection type
filtered_df <- conserved_CP_loops$filtered_data
write.csv(filtered_df, "HEPM_conserved_chromatinLoops_CPGenes.csv")
interaction_summary <- filtered_df %>%
  group_by(interaction_type) %>%
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
# plot
pie_chart <- ggplot(interaction_summary, aes(x = "", y = loop_count, fill = interaction_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Enhancer-Promoter" = "#CC5D4B", 
                               "Promoter-Promoter" = "#FAD8C4", 
                               "Promoter-Undefined" = "#F5F1BB")) +
  labs(title = "Combined-CP genes related chromatin Loop distribution(n=2688)",
       fill = "Interaction type") +
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
ggsave("AllCPgenes_436CPGenes_2688conservedloop_interaction_types_pie.pdf", pie_chart, width = 8, height = 6)

enhancer_analysis <- extract_enhancer_ids(filtered_df, Combined_CPGene, gene_col = "CP_Gene")
# Total Enhancer-Promoter interactions: 1224
#All CP genes loaded: 748
#Total E-P interactions: 1224
#all total loops: 2688
# unique CP genes: 436
#Total enhancer instances: 1519
#unique enhancers: 1050

# Now save the all chromatin loops and E-P loops :
E_P_Loops_CpGenes <- enhancer_analysis$E_P_Loops_OMIMGenes
write.csv(filtered_df, "2688conservedLoops_436CPGenes_MEPM.csv", row.names = FALSE)

# save 1050 enhancers associated with CP genes
E_ID <- data.frame(
  enhancer_id = enhancer_analysis$summary$unique_enhancer_ids,
  stringsAsFactors = FALSE
)
write.csv(E_ID, "1050_ConservedAEs_CPGenes.csv")

save(list = ls(), file = "Fig2.c.g.h.k.HEPM_MEPM_conserved_loops_genes.Rdata")

