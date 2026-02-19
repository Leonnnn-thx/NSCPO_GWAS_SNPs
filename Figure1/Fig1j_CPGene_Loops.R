save(list = ls(), file = "Fig2J_CombinedCPGenes_chromatinLoops.Rdata")
setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results")
getwd()
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotly)
library(scales)
library(stringr)
library(GenomicRanges)
Combined_CPGene <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/21_NSCPO_causativeGene/CP_combined_with_mouse_orthologs.csv", sep=",", header = T)

#--------- HEPM combined CP genes related chromatin loops -------------------------------------------
Data1 <- read.csv("./1.Loop_annotation/HEPM_microC_6400bp_annotate_loops_20251204.csv", header = T, sep = ",")

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

detailed_result <- detailed_CPgene_analysis(Data1, Combined_CPGene, gene_col = "CP_Gene")


cat("number of combined CP genes associated chromatin loop:", detailed_result$summary$total_loops, "\n") # 4928
cat("the number of CP genes in our chromatin loops:", detailed_result$summary$unique_genes, "\n") # 513
cat("Total CP genes:", detailed_result$summary$total_CP_genes, "\n\n") # 748

print(detailed_result$gene_frequency)

# ------------------Plot for CP gene frequency ----------------
library(ggplot2)
if(nrow(detailed_result$gene_frequency) > 0) {
  p <- ggplot(detailed_result$gene_frequency, aes(x = reorder(gene, count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Combined CP genes frequency in loops",
         x = "HEPM CP Genes", y = "Number") +
    theme_minimal()
  print(p)
  # save picture
  ggsave("Combined_CPgenes_in_chromatin_loops.pdf", plot = p, width = 10, height = 20)
}

#------------------- pie plot for loop intersection types associated wit CP genes --------------------------------
filtered_df <- detailed_result$filtered_data

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
pie_chart <- ggplot(interaction_summary, aes(x = "", y = loop_count, fill = interaction_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Enhancer-Promoter" = "#CC5D4B", 
                               "Promoter-Promoter" = "#FAD8C4", 
                               "Promoter-Undefined" = "#F5F1BB")) +
  labs(title = "All CP genes related chromatin Loop distribution(n=4928)",
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
ggsave("All748CPgenes_513CPGenes_4928chromatin_loop_interaction_types_pie.pdf", pie_chart, width = 8, height = 6)

# ———————————————— Extract E-P intersection and identify AEs linked to combined CP genes -----
# Function to extract enhancer IDs from Enhancer-Promoter interactions
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


# Execute the analysis
enhancer_analysis <- extract_enhancer_ids(filtered_df, Combined_CPGene, gene_col = "CP_Gene")

# Now save the all chromatin loops and E-P loops :
E_P_Loops_CpGenes <- enhancer_analysis$E_P_Loops_OMIMGenes
write.csv(E_P_Loops_CpGenes, "1765AEs_2218EPLoops_CPGenes.csv", row.names = FALSE)
write.csv(filtered_df, "4928Loops_513CPGenes_HEPM.csv", row.names = FALSE)

# save 1765 enhancers associated with CP genes
E_ID <- data.frame(
  enhancer_id = enhancer_analysis$summary$unique_enhancer_ids,
  stringsAsFactors = FALSE
)
write.csv(E_ID, "1765AEs_CPGenes.csv")
######################

#--------- HEPM combined CP genes related chromatin loops -------------------------------------------
Data2 <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/MEPM_MicroC_6400bp_annotatedLoops_20251204.csv", header = T, sep = ",")
MEPM_result <- detailed_CPgene_analysis(Data2, Combined_CPGene, gene_col = "mm_gene_symbol")

cat("number of combined CP genes associated chromatin loop:", MEPM_result$summary$total_loops, "\n") # 4988
cat("the number of CP genes in our chromatin loops:", MEPM_result$summary$unique_genes, "\n") # 576
cat("Total CP genes:", MEPM_result$summary$total_CP_genes, "\n\n") # 759

print(detailed_result$gene_frequency)

#------------------- pie plot for loop intersection types associated wit CP genes --------------------------------
filtered_df <- MEPM_result$filtered_data

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
pie_chart <- ggplot(interaction_summary, aes(x = "", y = loop_count, fill = interaction_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Enhancer-Promoter" = "#CC5D4B", 
                               "Promoter-Promoter" = "#FAD8C4", 
                               "Promoter-Undefined" = "#F5F1BB")) +
  labs(title = "Combined-CP genes related chromatin Loop distribution(n=4988)",
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
# ———————————————— Extract E-P intersection and identify AEs linked to combined CP genes -----

# Execute the analysis
enhancer_analysis <- extract_enhancer_ids(filtered_df, Combined_CPGene, gene_col = "mm_gene_symbol")

#Total Enhancer-Promoter interactions: 2009
#Valid OMIM genes loaded: 759
#=== Summary ===
#Total E-P interactions: 2009
#E-P interactions with OMIM genes: 2009 
#Interactions with extracted enhancers: 1813
#Total enhancer instances: 2889
#Unique enhancers: 2076
# unique genes: 576
# total loops: 4988
# Now save the all chromatin loops and E-P loops :
E_P_Loops_CpGenes <- enhancer_analysis$E_P_Loops_OMIMGenes
write.csv(E_P_Loops_CpGenes, "2076AEs_2009EPLoops_CPGenes_MEPM.csv", row.names = FALSE)
write.csv(filtered_df, "4988Loops_576CPGenes_MEPM.csv", row.names = FALSE)
ggsave("All759mouseCPgenes_576CPGenes_4988chromatin_loop_interaction_types_pie.pdf", pie_chart, width = 8, height = 6)

# save 2076 enhancers associated with CP genes
E_ID <- data.frame(
  enhancer_id = enhancer_analysis$summary$unique_enhancer_ids,
  stringsAsFactors = FALSE
)
write.csv(E_ID, "2076AEs_mouseCPGenes.csv")
######################

