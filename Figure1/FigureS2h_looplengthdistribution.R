#- plot Figure 2.f loops length distribution of HEPM and MEPM --

library(ggplot2)
library(dplyr)

# - HEPM loop -
Data1 <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/19_microC_data/Results/1.Loop_annotation/HEPM_microC_6400bp_annotate_loops_20251204.csv")
Data_HEPM <- Data1 %>%
  mutate(
    loop_length = abs(((start2 + end2) / 2) - ((start1 + end1) / 2))
  )

plot_data <- Data_HEPM %>%
  filter(interaction_type %in% c("Enhancer-Promoter", "Promoter-Promoter", "Enhancer-Enhancer"))

plot_data$interaction_type <- factor(plot_data$interaction_type, levels = c("Enhancer-Enhancer", "Promoter-Promoter", "Enhancer-Promoter"))

# - MEPM loop -
Data2 <- read.csv("/Users/thx/Desktop/PH.D_candidate/1_Project/B_NSCPO_FunctionalVariants_CRE/2_analysis_process/22_MEPM_MicroC/2.145_Mes_microc/Results/2.annotation_loop_mm/MEPM_MicroC_6400bp_annotatedLoops_20251204.csv")
Data_mm10 <- Data2 %>%
  mutate(
    loop_length = abs(((start2 + end2) / 2) - ((start1 + end1) / 2))
  )

plot_data_mm10 <- Data_mm10 %>%
  filter(interaction_type %in% c("Enhancer-Promoter", "Promoter-Promoter", "Enhancer-Enhancer"))

plot_data_mm10$interaction_type <- factor(plot_data_mm10$interaction_type, levels = c("Enhancer-Enhancer", "Promoter-Promoter", "Enhancer-Promoter"))


plot_data_combined <- plot_data %>%
  mutate(group = "HEPM") %>%
  bind_rows(
    plot_data_mm10 %>%
      mutate(group = "Mice-Mes-Palate")
  )

plot_data_combined$interaction_type <- factor(plot_data_combined$interaction_type, 
                                              levels = c("Enhancer-Promoter", "Promoter-Promoter", "Enhancer-Enhancer"))

plot_data_combined$group <- factor(plot_data_combined$group, 
                                   levels = c("HEPM", "Mice-Mes-Palate"))

median_labels_combined <- plot_data_combined %>%
  group_by(interaction_type, group) %>%
  summarise(
    median_len = median(loop_length, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    median_text = paste0(round(median_len / 1000, 1), "kb"),
    x_pos = case_when(
      group == "HEPM" ~ as.numeric(interaction_type) - 0.2,
      group == "Mice-Mes-Palate" ~ as.numeric(interaction_type) + 0.2
    )
  )

p_combined <- ggplot(plot_data_combined, 
                     aes(x = interaction_type, y = loop_length, 
                         fill = interaction_type, color = interaction_type, alpha = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               width = 0.5,
               size = 0.8) +
  geom_text(
    data = median_labels_combined,
    aes(x = x_pos, y = 2e6, label = median_text), 
    hjust = 0, 
    size = 3.5,
    vjust = -0.5,
    color = "black",
    alpha = 1,
    inherit.aes = FALSE
  ) +
  coord_flip(
    ylim = c(5000, 2e6), 
    clip = 'off'
  ) +
  scale_y_log10(
    name = "Loop length", 
    breaks = c(10000, 50000, 100000, 1000000), 
    labels = c("10kb", "50kb", "100kb", "1Mb")  
  ) +
  scale_fill_manual(
    name = "Interaction Type",
    values = c(
      "Enhancer-Promoter" = "#D6604D",   
      "Promoter-Promoter" = "#FDDBC7",   
      "Enhancer-Enhancer" = "#4393C3"    
    )
  ) +
  scale_color_manual(
    name = "Interaction Type",
    values = c(
      "Enhancer-Promoter" = "#67001F",   
      "Promoter-Promoter" = "#F4A582",   
      "Enhancer-Enhancer" = "#2166AC"    
    )
  ) +
  scale_alpha_manual(
    name = "Dataset",
    values = c(
      "HEPM" = 0.5,                      
      "Mice-Mes-Palate" = 0.0           
    ),
    guide = guide_legend(override.aes = list(
      fill = c("#808080", "white"),      
      color = c("#808080", "#808080"),   
      alpha = c(1, 1)                    
    ))
  ) +
  labs(
    x = "Subtypes",
    y = "Loop length",
    title = "Loop Length Distribution by Interaction Type and Dataset"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), 
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(5.5, 80, 5.5, 5.5, "pt")
  ) +
  guides(
    fill = guide_legend(order = 1),
    alpha = guide_legend(order = 2)
  )
print(p_combined)
ggsave("combined_loop_length_boxplot.pdf", 
       plot = p_combined, 
       width = 12, 
       height = 6,
       dpi = 300)

