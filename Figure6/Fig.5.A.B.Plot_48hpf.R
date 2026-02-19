#' Description: subset timepoint your interseted from reference scRNA-data and check genes expression in featureplot format.

setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/5_Expression_pattern_of_candidategenes/5_danRe_scRNA/zebra_scRNA_2024nature/tutorial/48hpf")
getwd()
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)
  library(qs)
  library(mascarade)
  library(ggrepel)
  library(showtext)
})
source("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/5_Expression_pattern_of_candidategenes/5_danRe_scRNA/zebra_scRNA_2024nature/tutorial/saunders_srivatsan_2023_utils.R")
ref_cds <- readRDS("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/5_Expression_pattern_of_candidategenes/5_danRe_scRNA/zebra_scRNA_2024nature/tutorial/zscape_reference_cds.RDS")
# select timepoint == 48hpf scRNA data from reference cds data
Time_48hpf_logical <- !is.na(colData(ref_cds)$timepoint) & 
  colData(ref_cds)$timepoint == "48hpf"
Time_48hpf_cds <- ref_cds[, Time_48hpf_cell]

# check the data
print(dim(Time_48hpf_cds))  # input as [Gene count, cell count]
head(colData(Time_48hpf_cds))  
# save 
qsave(Time_48hpf_cds, "zebra_timepoint_48_filtered_by_tissue.qs")
ref_cds <- qread("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/5_Expression_pattern_of_candidategenes/5_danRe_scRNA/zebra_scRNA_2024nature/tutorial/48hpf/zebra_timepoint_48_filtered_by_tissue.qs")

ref_coldata <-  
  ref_cds %>%
  colData() %>%
  as.data.frame()

# Extract UMAP coordinates from ref_cds to ref_coldata.
ref_coldata$umap3d_1 <- reducedDim(ref_cds, type = "UMAP")[, 1]  # UMAP1 
ref_coldata$umap3d_2 <- reducedDim(ref_cds, type = "UMAP")[, 2]  # UMAP2

# Define colors for time points -------------------------------------------

# Figure 1d ---------------------------------------------------------------

vibrant.colors <-  
  c('#EE7733', '#0077BB', '#228833', '#33BBEE', '#EE3377', '#CC3311',
    '#AA3377', '#009988', '#004488', '#DDAA33', '#99CC66','#D590DD')

bright.colors <-  
  c('#4477AA', 
    '#EE6677', 
    '#228833', 
    '#CCBB44', 
    '#66CCEE',
    '#AA3377', 
    '#BBBBBB')


num.colors.tissue <-  
  ref_coldata %>%
  filter(!is.na(tissue)) %>%
  pull(tissue) %>%
  unique() %>%
  sort() %>%
  length()

tissue.colors <- colorRampPalette(c(vibrant.colors,bright.colors))(num.colors.tissue) 

# plot cells by tissue

maskTable <- generateMask(
  dims = ref_coldata[, c("umap3d_1", "umap3d_2")],  
  cluster = ref_coldata$tissue,                     
  minDensity = 10,                                  
  smoothSigma = 0.005                              
)

p1 <- ggplot(
  ref_coldata %>%
    sample_n(size = dim(ref_coldata)[1]) %>%
    filter(!is.na(timepoint)) %>%
    group_by(tissue) %>%
    add_tally() %>%
    arrange(-n)
) +
  geom_point(
    aes(x = umap3d_1, y = umap3d_2),
    color = "black",
    stroke = 0,
    size = 0.2
  ) +
  geom_point(
    aes(
      x = umap3d_1,
      y = umap3d_2,
      color = as.character(tissue)
    ),
    stroke = 0,
    size = 0.25
  ) +
  geom_path(
    data = maskTable[cluster == "Pharyngeal Arch", ],
    aes(
      x = umap3d_1,  
      y = umap3d_2,  
      group = cluster
    ),
    linewidth = 0.5,
    linetype = 2
  ) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = tissue.colors)
p1
##########
z <- ggplot(
  ref_coldata %>%
    sample_n(size = dim(ref_coldata)[1]) %>%
    filter(!is.na(timepoint)) %>%
    group_by(tissue) %>%
    add_tally() %>%
    arrange(-n)
) +
  geom_point(
    aes(x = umap3d_1, y = umap3d_2),
    color = "black",
    stroke = 0,
    size = 0.2
  ) +
  geom_point(
    aes(
      x = umap3d_1,
      y = umap3d_2,
      color = as.character(tissue) 
    ),
    stroke = 0,
    size = 0.25
  ) +
  geom_path(
    data = maskTable[cluster == "Pharyngeal Arch", ],
    aes(
      x = umap3d_1,  
      y = umap3d_2,  
      group = cluster
    ),
    linewidth = 0.5,
    linetype = 2
  ) +
  theme_void() +
  scale_color_manual(
      values = tissue.colors,
      name = "Tissue"  
  ) +
  theme(
    legend.position = "right" 
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3)
    )
  )

z
####

tissue_centers <- ref_coldata %>%
  group_by(tissue) %>%
  summarize(
    umap3d_1 = median(umap3d_1), 
    umap3d_2 = median(umap3d_2)  
  )

G <- p1 +
  geom_text_repel(
    data = tissue_centers,
    aes(
      x = umap3d_1,
      y = umap3d_2,
      label = tissue 
    ),
    size = 2, 
    color = "black"
  )
G

## save
ggsave("zebra_48hpf_umap_v1.tiff", plot = p1,
       dpi = 750,
       height = 5,
       width = 5,
       bg = "transparent")
ggsave("zebra_48hpf_umap_v2.pdf", plot = z,
       dpi = 100,
       height = 5,
       width = 8,
       bg = "transparent")
ggsave("zebra_48hpf_umap_v3.pdf", plot = G,
       dpi = 750,
       height = 5,
       width = 5,
       bg = "transparent")

###### Plot : featureplot of gene (hhat) at 48 hpf##########
genes_list <- c("hhat","acvr2aa","acvr2ab", "prrx1b", "hand2", "zbtb47b")
genes_list <- c("hhat","tlnrd1","tle3a", "tle3b", "spsb1", "Traf3ip2")
genes_list <- c("ENSDARG00000091734")

font_add(family = "Arial", regular = "arial.ttf")
showtext_auto()
for (g in genes_list){
  p <- plot_cells(ref_cds, 
                  x = 1, 
                  y = 2,  
                  genes = g, 
                  label_cell_groups = FALSE, 
                  show_trajectory_graph = FALSE,
                  cell_size = 0.2, 
                  cell_stroke = 0, 
                  alpha = 0.8) +
    scale_color_viridis_c() +
    theme_void()
  p <- p + geom_path(
    data = maskTable[cluster == "Pharyngeal Arch", ],
    aes(
      x = umap3d_1,  
      y = umap3d_2,  
      group = cluster
    ),
    color = "black",
    linewidth = 0.5,
    linetype = 2
  )
  p <- p + theme(
    text = element_text(family = "Arial"), 
    
    legend.position = "right",
    legend.text = element_text(size = 70),   
    legend.title = element_text(size = 70),  
    legend.key.width = unit(0.01, "in"),     
    legend.key.height = unit(0.05, "in")    
  ) +
    guides(
      color = guide_colorbar(
        title = "Expression",  
        barwidth = 0.5,        
        barheight = 3        
      )
    )
  ggsave(
    paste0("Pharyngeal_arch_", g, "_48hpf_umap.tiff"),
    plot = p,
    device = "tiff",
    dpi = 1000,          
    height = 2.5,
    width = 3,
    units = "in"         
  )
}


 showtext_auto(FALSE)

