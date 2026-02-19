setwd("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/5_Expression_pattern_of_candidategenes/5_danRe_scRNA/zebra_scRNA_2024nature/tutorial")
getwd()
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(monocle3)
  library(data.table)
  library(qs)
})
source("saunders_srivatsan_2023_utils.R")
ref_cds <- readRDS("zscape_reference_cds.RDS")
# subset Pharyngeal Arch from ref_cds
pharyngeal_arch_cells <- colData(ref_cds)$cell[colData(ref_cds)$tissue == "Pharyngeal Arch"]
pharyngeal_arch_cds <- ref_cds[, pharyngeal_arch_cells]


print(dim(pharyngeal_arch_cds))  #  [Gene, Cell]
head(colData(pharyngeal_arch_cds)) 
qsave(pharyngeal_arch_cds, "pharyngeal_arch_subset.qs")

PA_cds <- qread("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/5_Expression_pattern_of_candidategenes/5_danRe_scRNA/zebra_scRNA_2024nature/tutorial/pharyngeal_arch_subset.qs")


colData(PA_cds)$sub_umap1 <- reducedDim(x = PA_cds,
                                        type = "UMAP")[,1]
colData(PA_cds)$sub_umap2 <- reducedDim(x = PA_cds,
                                        type = "UMAP")[,2]
colData(PA_cds)$sub_umap3 <- reducedDim(x = PA_cds,
                                        type = "UMAP")[,3]


PA_coldata <- colData(PA_cds) %>% 
  as.data.frame()

######## Plot 1. PharyngealArch-highlight_umap  ###########

# load full reference metadata
global_coldata <- fread("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/5_Expression_pattern_of_candidategenes/5_danRe_scRNA/zebra_scRNA_2024nature/tutorial/reference_cell_metadata.csv", 
                        sep = ",", stringsAsFactors = F, data.table = F) # Download from ZSCAPE
head(global_coldata)

PA_cells <- PA_coldata %>% 
  pull(cell)

PA_cells <- intersect(PA_cells, global_coldata$cell)

p1 <- ggplot() +
  geom_point(data = global_coldata, 
             aes(x = umap3d_1,
                 y = umap3d_2),
             color = "grey90",
             stroke = 0,
             size = 0.1) +
  geom_point(data = global_coldata %>% 
               filter(cell %in% PA_cells),
             aes(x = umap3d_1,
                 y = umap3d_2),
             color = "black",
             stroke = 0,
             size = 0.1) +
  theme_void() +
  theme(legend.position = "none")
ggsave("1.ref_cells_PharyngealArch-highlight_umap.pdf",
       plot = p1, 
       dpi = 750,
       height = 1.5,
       width = 1.5,
       bg = "transparent")

######### Plot2.    Paryngeal Arch UMAP by timepoint  #######
PA_coldata_numeric <- PA_coldata %>%
  filter(!is.na(timepoint)) %>%
  mutate(time_numeric = as.numeric(gsub("h", "", timepoint)))
# get timepoint colors
rainbow_timepoint_colors <- 
  c("18h" = "#62b5d9", 
    "24h" = "#e5a18f",
    "36h" = "#c840a9", 
    "48h" = "#ffe6d8", 
    "72h" = "#7e62d5",
    "96h" = "#00a696"
    )

p2_with_colorbar <- ggplot(
  PA_coldata_numeric %>% sample_n(size = dim(.)[1]),
  aes(x = sub_umap1, y = sub_umap3, color = time_numeric)
) +
  geom_point(color = "black", stroke = 0, size = 0.5) +
  geom_point(stroke = 0, size = 0.4) +
  scale_color_gradientn(
    colors = unname(rainbow_timepoint_colors), 
    name = "Timepoint (hpf)",
    breaks = c(min(PA_coldata_numeric$time_numeric), max(PA_coldata_numeric$time_numeric)) 
  ) +
  theme_void() +
  theme(
    legend.position = c(0.7, 0.9),     
    legend.direction = "horizontal",    
    legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 10, face = "bold")
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",           
      barwidth = unit(8, "lines"),     
      barheight = unit(0.8, "lines")    
    )
  )
p2_with_colorbar





ggsave("2.all_Paryngeal-Arch_timepoint_umap_v2.pdf",
       plot = p2_with_colorbar,
       dpi = 750,
       height = 4,
       width = 6,
       bg = "transparent")

######### plot 3. pharyngeal arch cell type in pesudotime ######
cell_type_sub <- PA_coldata %>% filter(!is.na(cell_type_sub)) %>% pull(cell_type_sub) %>% unique()
cell_type_sub
cell_tpye_broad <- PA_coldata %>% filter(!is.na(cell_type_broad)) %>% pull(cell_type_broad) %>% unique()
cell_tpye_broad

PA_colors = c(
  "head mesenchyme (maybe ventral, hand2+)" = "#AEC6CF",  
  "head mesenchyme/PA cartilage" = "#AEC6CF",            
  "pharyngeal arch (early)" = "#45B7D1",                  
  "chondrocranium" = "#96CEB4",                           
  "pharyngeal arch (contains muscle, early cartilage)" = "#45B7D1",
  "cranial muscle (early)" = "#C23B22",                  
  "head mesenchyme" = "#AEC6CF",                          
  "vascular smooth muscle" = "#B19CD9",                  
  "pharyngeal arch (NC-derived)" = "#45B7D1",            
  "cranial muscle (progenitor)" = "#C23B22",              
  "jaw chondrocyte" = "#F49AC2",                         
  "cranial muscle (mid)" = "#C23B22"                     
)

cell_type_centers <- PA_coldata %>%
  filter(!is.na(cell_type_broad)) %>%
  group_by(cell_type_broad) %>%
  summarize(
    umap1_center = median(sub_umap1),
    umap3_center = median(sub_umap3)
  )


p3 <- ggplot(PA_coldata %>%
              sample_n(size = nrow(PA_coldata))) +
  geom_point(aes(x = sub_umap1, y = sub_umap3),
             color = "black", stroke = 0, size = 0.2) +
  geom_point(aes(x = sub_umap1, y = sub_umap3,
                 color = cell_type_broad),
             stroke = 0, size = 0.2) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = PA_colors)
p3
p4 <- p3 +
  geom_text_repel(
    data = cell_type_centers,  
    aes(x = umap1_center, y = umap3_center, label = cell_type_broad),
    size = 5,                  
    color = "black",           
    fontface = "bold",         
    box.padding = unit(0.5, "lines"), 
    point.padding = unit(0.3, "lines"),
    max.overlaps = Inf         
  )
p4
ggsave("3.Paryngeal-Arch_subcelltype_umap.pdf", plot = p3,
       dpi = 750,
       height = 2,
       width = 2.2,
       bg = "transparent")


##########


######### plot 4. featureplot of gene (hhat)  ##########
genes_list <- c("hhat")
for (g in genes_list){
  p <- plot_cells(PA_cds, 
                  x = 1, 
                  y = 3,  
                  genes = g, 
                  label_cell_groups = FALSE, 
                  show_trajectory_graph = FALSE,
                  cell_size = 0.2, 
                  cell_stroke = 0, 
                  alpha = 0.8) +
    scale_color_viridis_c(
      name = expression("expression"),  
      limits = c(0, 100),                   
      breaks = seq(0, 100, 25),            
      labels = c("0", "25", "50", "75", "100")           
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(
        size = 8, 
        margin = margin(b = 2),      
        hjust = 0.5                  
      ),
      legend.text = element_text(
        size = 6,
        margin = margin(t = 1)       
      ),
      legend.key.height = unit(0.8, "cm"),  
      legend.key.width = unit(0.2, "cm")   
    ) +
    guides(
      color = guide_colorbar(
        ticks.linewidth = 0.1,        
        ticks.colour = "white",       
        barwidth = 0.5,              
        barheight = 1.5,              
        title.position = "top",       
        title.hjust = 0.5            
      )
    )
  
  ggsave(
    paste0("Pharyngeal_arch_", g, "_expr_umap.tiff"),
    plot = p,
    dpi = 750,
    height = 2.5,
    width = 3,
    bg = "transparent"
  )
}

