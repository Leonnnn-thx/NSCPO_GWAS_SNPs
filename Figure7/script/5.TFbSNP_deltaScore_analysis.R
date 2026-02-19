# ============================================================================
# TFbSNP Analysis Pipeline: 完整功能函数
# ============================================================================
# 作者: THX
# 日期: 2024
# 描述: 分析 TFbSNP 数据，包括频率分布、deltaSVM 评分和可视化
# ============================================================================

#' TFbSNP 完整分析流程
#'
#' @param input_csv 输入的 CSV 文件路径
#' @param genome_version 基因组版本，默认 "hg19"，可选 "hg38"
#' @param flank_size 侧翼序列长度，默认 9
#' @param ref_fasta 参考等位基因 FASTA 文件名，默认 "reference_alleles.fa"
#' @param alt_fasta 替代等位基因 FASTA 文件名，默认 "alt_alleles.fa"
#' @param outfn_snp deltaSVM 输出文件名，默认 "deltaSVM_score.out"
#' @param svm_model_prefix gkmSVM 模型文件前缀路径
#' @param output_prefix 输出文件前缀，默认 "TFbSNP_analysis"
#' @param plot_top_n 在散点图中标记 Top N 个 SNP，默认 10
#' @param run_gkmsvm 是否运行 gkmSVM 分析，默认 TRUE
#' @param verbose 是否显示详细信息，默认 TRUE
#'
#' @return 返回包含分析结果的列表
#' @export
#'
#' @examples
#' result <- analyze_TFbSNP(
#'   input_csv = "AEs_805_TFbSNP_20251127.csv",
#'   genome_version = "hg19",
#'   flank_size = 9,
#'   svm_model_prefix = "/path/to/svm/model",
#'   output_prefix = "my_analysis"
#' )

analyze_TFbSNP <- function(
    input_csv,
    genome_version = "hg19",
    flank_size = 9,
    ref_fasta = "reference_alleles.fa",
    alt_fasta = "alt_alleles.fa",
    outfn_snp = "deltaSVM_score.out",
    svm_model_prefix = NULL,
    output_prefix = "TFbSNP_analysis",
    plot_top_n = 10,
    run_gkmsvm = TRUE,
    verbose = TRUE
) {
  
  # ============================================================================
  # 步骤 0: 加载必需的包
  # ============================================================================
  
  if (verbose) cat("\n", rep("=", 80), "\n", sep = "")
  if (verbose) cat("TFbSNP Analysis Pipeline\n")
  if (verbose) cat(rep("=", 80), "\n\n", sep = "")
  
  required_packages <- c("ggplot2", "dplyr", "scales", "ggrepel")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (verbose) cat("安装缺失的包:", pkg, "\n")
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
  
  # 加载基因组包
  genome_pkg <- paste0("BSgenome.Hsapiens.UCSC.", genome_version)
  if (!require(genome_pkg, character.only = TRUE, quietly = TRUE)) {
    if (verbose) cat("安装基因组包:", genome_pkg, "\n")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(genome_pkg)
    library(genome_pkg, character.only = TRUE)
  }
  
  # 获取基因组对象
  genome <- get(genome_pkg)
  
  if (run_gkmsvm) {
    if (!require("gkmSVM", quietly = TRUE)) {
      if (verbose) cat("安装 gkmSVM 包\n")
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("gkmSVM")
      library(gkmSVM)
    }
  }
  
  # ============================================================================
  # 步骤 1: 读取数据
  # ============================================================================
  
  if (verbose) cat("步骤 1: 读取输入数据\n")
  if (verbose) cat(rep("-", 80), "\n", sep = "")
  
  if (!file.exists(input_csv)) {
    stop("错误: 输入文件不存在: ", input_csv)
  }
  
  TFbSNP <- read.csv(input_csv, header = TRUE, stringsAsFactors = FALSE)
  
  if (verbose) {
    cat("输入文件:", input_csv, "\n")
    cat("总 SNP 数:", nrow(TFbSNP), "\n")
    cat("列名:", paste(names(TFbSNP), collapse = ", "), "\n\n")
  }
  
  # ============================================================================
  # 步骤 2: 频率分布分析
  # ============================================================================
  
  if (verbose) cat("步骤 2: 分析频率分布\n")
  if (verbose) cat(rep("-", 80), "\n", sep = "")
  
  freq_stats <- TFbSNP %>%
    group_by(FREQUENCY_CLASS) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    arrange(desc(Count)) %>%
    mutate(
      Percentage = Count / sum(Count) * 100,
      Fraction = Count / sum(Count),
      Ypos = cumsum(Fraction) - 0.5 * Fraction,
      PercentLabel = paste0(round(Percentage, 1), "%"),
      Label = paste0(FREQUENCY_CLASS, "\n", 
                     format(Count, big.mark = ","), "\n",
                     "(", round(Percentage, 1), "%)")
    )
  
  if (verbose) {
    cat("\n频率分布统计:\n")
    print(freq_stats %>% select(FREQUENCY_CLASS, Count, Percentage))
    cat("\n")
  }
  
  # 绘制饼图
  pie_plot <- ggplot(freq_stats, aes(x = "", y = Count, fill = FREQUENCY_CLASS)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 0) +
    coord_polar("y", start = 0) +
    geom_text(aes(y = Ypos * sum(Count), 
                  label = paste0(FREQUENCY_CLASS, "\n", Count)),
              size = 3, fontface = "bold", color = "white") +
    geom_text(aes(y = Ypos * sum(Count), label = PercentLabel),
              position = position_nudge(x = 0.6),
              size = 5, fontface = "bold") +
    scale_fill_manual(
      values = c("RARE" = "#e74c3c", 
                 "COMMON" = "#3498db", 
                 "UNKNOWN" = "#95a5a6"),
      name = "Frequency Class"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 20)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
      legend.position = "bottom",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = "SNP Frequency Class Distribution",
      subtitle = paste0("Total SNPs: ", format(nrow(TFbSNP), big.mark = ","))
    )
  
  # 保存饼图
  pie_filename_pdf <- paste0(output_prefix, "_frequency_pie.pdf")
  pie_filename_png <- paste0(output_prefix, "_frequency_pie.png")
  
  ggsave(pie_filename_pdf, pie_plot, width = 10, height = 7, dpi = 300)
  ggsave(pie_filename_png, pie_plot, width = 10, height = 7, dpi = 300)
  
  if (verbose) cat("✅ 饼图已保存:", pie_filename_pdf, "和", pie_filename_png, "\n\n")
  
  # ============================================================================
  # 步骤 3: 准备 FASTA 序列
  # ============================================================================
  
  if (verbose) cat("步骤 3: 生成参考和替代等位基因序列\n")
  if (verbose) cat(rep("-", 80), "\n", sep = "")
  
  # 过滤 UNKNOWN 频率的 SNP
  TFbSNP_gkmSVM_list <- TFbSNP %>% 
    dplyr::filter(FREQUENCY_CLASS != "UNKNOWN") %>% 
    dplyr::select(CHROM, POS, ID, REF, ALT)
  
  if (verbose) {
    cat("过滤后的 SNP 数:", nrow(TFbSNP_gkmSVM_list), "\n")
    cat("侧翼长度:", flank_size, "bp\n")
    cat("基因组版本:", genome_version, "\n\n")
  }
  
  # 删除已存在的 FASTA 文件
  if (file.exists(ref_fasta)) file.remove(ref_fasta)
  if (file.exists(alt_fasta)) file.remove(alt_fasta)
  
  # 生成 FASTA 序列
  success_count <- 0
  failed_count <- 0
  failed_snps <- c()
  
  for (i in 1:nrow(TFbSNP_gkmSVM_list)) {
    snp_info <- TFbSNP_gkmSVM_list[i, ]
    chrom <- snp_info$CHROM
    snp_pos <- snp_info$POS
    ref_allele <- snp_info$REF
    alt_allele <- snp_info$ALT
    snp_id <- snp_info$ID
    
    if (verbose && (i %% 100 == 0 || i == 1)) {
      cat(sprintf("处理进度: %d / %d (%.1f%%) - %s\n", 
                  i, nrow(TFbSNP_gkmSVM_list), 
                  i/nrow(TFbSNP_gkmSVM_list)*100, snp_id))
    }
    
    start_pos <- snp_pos - flank_size
    end_pos <- snp_pos + flank_size
    
    reference_sequence <- tryCatch({
      getSeq(genome, chrom, start_pos, end_pos)
    }, error = function(e) {
      if (verbose) cat(sprintf("⚠️  无法提取 SNP %s 的序列\n", snp_id))
      return(NULL)
    })
    
    if (is.null(reference_sequence)) {
      failed_count <- failed_count + 1
      failed_snps <- c(failed_snps, snp_id)
      next 
    }
    
    reference_sequence_char <- as.character(reference_sequence)
    
    # 检查并替换中心位置的碱基
    center_base_pos <- flank_size + 1
    center_base <- substr(reference_sequence_char, center_base_pos, center_base_pos)
    if (center_base != ref_allele) {
      substr(reference_sequence_char, center_base_pos, center_base_pos) <- ref_allele
    }
    
    # 生成替代序列
    alt_sequence_char <- reference_sequence_char
    substr(alt_sequence_char, center_base_pos, center_base_pos) <- alt_allele
    
    # 创建 FASTA header
    ref_fasta_header <- paste0(">", snp_id, "_", chrom, ":", snp_pos, "_", ref_allele, "_", alt_allele)  
    alt_fasta_header <- paste0(">", snp_id, "_", chrom, ":", snp_pos, "_", alt_allele)  
    
    ref_fasta_content <- paste(ref_fasta_header, reference_sequence_char, sep = "\n")
    alt_fasta_content <- paste(alt_fasta_header, alt_sequence_char, sep = "\n")
    
    cat(ref_fasta_content, file = ref_fasta, append = TRUE, sep = "\n")
    cat(alt_fasta_content, file = alt_fasta, append = TRUE, sep = "\n")
    
    success_count <- success_count + 1
  }
  
  if (verbose) {
    cat("\n")
    cat("FASTA 生成完成:\n")
    cat("  成功:", success_count, "\n")
    cat("  失败:", failed_count, "\n")
    cat("  参考序列文件:", ref_fasta, "\n")
    cat("  替代序列文件:", alt_fasta, "\n\n")
  }
  
  # ============================================================================
  # 步骤 4: 运行 gkmSVM (可选)
  # ============================================================================
  
  if (run_gkmsvm) {
    if (verbose) cat("步骤 4: 运行 gkmSVM 分析\n")
    if (verbose) cat(rep("-", 80), "\n", sep = "")
    
    if (is.null(svm_model_prefix)) {
      warning("未提供 svm_model_prefix，跳过 gkmSVM 分析")
      run_gkmsvm <- FALSE
    } else {
      if (verbose) cat("SVM 模型前缀:", svm_model_prefix, "\n")
      
      tryCatch({
        gkmsvm_delta(ref_fasta, alt_fasta, svm_model_prefix, outfn_snp)
        if (verbose) cat("✅ gkmSVM 分析完成，结果保存至:", outfn_snp, "\n\n")
      }, error = function(e) {
        warning("gkmSVM 运行失败: ", e$message)
        run_gkmsvm <<- FALSE
      })
    }
  }
  
  # ============================================================================
  # 步骤 5: 整合 deltaSVM 分数
  # ============================================================================
  
  if (run_gkmsvm && file.exists(outfn_snp)) {
    if (verbose) cat("步骤 5: 整合 deltaSVM 分数\n")
    if (verbose) cat(rep("-", 80), "\n", sep = "")
    
    score <- read.table(outfn_snp, header = FALSE, stringsAsFactors = FALSE)
    names(score) <- c("V1", "V2")
    
    TFbSNP_score <- TFbSNP %>%
      left_join(
        score %>%
          mutate(ID = sub("_.*", "", V1)) %>%
          select(ID, V2),
        by = "ID"
      ) %>%
      dplyr::rename(deltaScore = V2) %>% 
      dplyr::filter(FREQUENCY_CLASS != "UNKNOWN")
    
    if (verbose) {
      cat("deltaSVM 分数统计:\n")
      cat("  有分数的 SNP:", sum(!is.na(TFbSNP_score$deltaScore)), "\n")
      cat("  缺失分数的 SNP:", sum(is.na(TFbSNP_score$deltaScore)), "\n")
      print(summary(TFbSNP_score$deltaScore))
      cat("\n")
    }
    
  } else {
    if (verbose) cat("步骤 5: 跳过 deltaSVM 分数整合（未运行 gkmSVM）\n\n")
    TFbSNP_score <- TFbSNP %>%
      dplyr::filter(FREQUENCY_CLASS != "UNKNOWN")
  }
  
  # ============================================================================
  # 步骤 6: 绘制散点图
  # ============================================================================
  
  if (run_gkmsvm && "deltaScore" %in% names(TFbSNP_score)) {
    if (verbose) cat("步骤 6: 绘制散点图\n")
    if (verbose) cat(rep("-", 80), "\n", sep = "")
    
    # 准备绘图数据
    plot_data <- TFbSNP_score %>%
      filter(!is.na(deltaScore), 
             !is.na(RELATIVE_FOOTPRINT_SCORE),
             FREQUENCY_CLASS %in% c("RARE", "COMMON"))
    
    if (nrow(plot_data) == 0) {
      warning("没有可用于绘图的数据")
    } else {
      # 计算分位数
      q90_deltaScore <- quantile(plot_data$deltaScore, 0.90, na.rm = TRUE)
      q10_deltaScore <- quantile(plot_data$deltaScore, 0.10, na.rm = TRUE)
      q75_footprint <- quantile(plot_data$RELATIVE_FOOTPRINT_SCORE, 0.75, na.rm = TRUE)
      
      if (verbose) {
        cat("分位数统计:\n")
        cat("  deltaScore Q10:", round(q10_deltaScore, 4), "\n")
        cat("  deltaScore Q90:", round(q90_deltaScore, 4), "\n")
        cat("  Footprint Q75:", round(q75_footprint, 4), "\n\n")
      }
      
      # 基础散点图
      p_basic <- ggplot(plot_data, aes(x = deltaScore, y = RELATIVE_FOOTPRINT_SCORE)) +
        geom_point(aes(color = FREQUENCY_CLASS), alpha = 0.6, size = 2) +
        scale_color_manual(
          values = c("RARE" = "#e74c3c", "COMMON" = "#3498db"),
          name = "Frequency Class"
        ) +
        geom_vline(xintercept = q90_deltaScore, linetype = "dashed", color = "gray30", size = 0.8) +
        geom_vline(xintercept = q10_deltaScore, linetype = "dashed", color = "gray30", size = 0.8) +
        geom_hline(yintercept = q75_footprint, linetype = "dashed", color = "gray30", size = 0.8) +
        annotate("text", x = q90_deltaScore, y = max(plot_data$RELATIVE_FOOTPRINT_SCORE) * 0.95,
                 label = paste0("Q90 = ", round(q90_deltaScore, 2)),
                 hjust = -0.1, vjust = 1, color = "gray30", size = 4, fontface = "bold") +
        annotate("text", x = q10_deltaScore, y = max(plot_data$RELATIVE_FOOTPRINT_SCORE) * 0.95,
                 label = paste0("Q10 = ", round(q10_deltaScore, 2)),
                 hjust = 1.1, vjust = 1, color = "gray30", size = 4, fontface = "bold") +
        annotate("text", x = min(plot_data$deltaScore) * 0.95, y = q75_footprint,
                 label = paste0("Q75 = ", round(q75_footprint, 3)),
                 hjust = 0, vjust = -0.5, color = "gray30", size = 4, fontface = "bold") +
        labs(
          title = "Relationship between deltaSVM Score and Footprint Score",
          subtitle = paste0("Total SNPs: ", format(nrow(plot_data), big.mark = ","),
                            " | RARE: ", sum(plot_data$FREQUENCY_CLASS == "RARE"),
                            " | COMMON: ", sum(plot_data$FREQUENCY_CLASS == "COMMON")),
          x = "deltaSVM Score",
          y = "Relative Footprint Score"
        ) +
        theme_bw(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
          legend.position = "top",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(size = 11)
        )
      
      # 保存基础散点图
      basic_plot_pdf <- paste0(output_prefix, "_scatter_basic.pdf")
      basic_plot_png <- paste0(output_prefix, "_scatter_basic.png")
      
      ggsave(basic_plot_pdf, p_basic, width = 10, height = 8, dpi = 300)
      ggsave(basic_plot_png, p_basic, width = 10, height = 8, dpi = 300)
      
      if (verbose) cat("✅ 基础散点图已保存:", basic_plot_pdf, "\n")
      
      # Top N SNPs 散点图
      if (plot_top_n > 0) {
        top_snps <- plot_data %>%
          arrange(desc(RELATIVE_FOOTPRINT_SCORE)) %>%
          head(plot_top_n) %>%
          mutate(rank = row_number())
        
        if (verbose) {
          cat("\nTop", plot_top_n, "SNPs:\n")
          print(top_snps %>% select(rank, ID, RELATIVE_FOOTPRINT_SCORE, deltaScore))
        }
        
        p_top <- ggplot(plot_data, aes(x = deltaScore, y = RELATIVE_FOOTPRINT_SCORE)) +
          geom_point(aes(color = FREQUENCY_CLASS), alpha = 0.6, size = 2) +
          geom_point(data = top_snps, color = "black", size = 4, shape = 21, 
                    stroke = 1.5) +
          geom_text_repel(data = top_snps,
                          aes(label = ID),
                          size = 3.5, fontface = "bold", color = "black",
                          bg.color = "white", bg.r = 0.1,
                          box.padding = 0.5, point.padding = 0.3,
                          segment.color = "gray40", segment.size = 0.5,
                          arrow = arrow(length = unit(0.01, "npc")),
                          max.overlaps = 20, min.segment.length = 0) +
          scale_color_manual(
            values = c("RARE" = "#e74c3c", "COMMON" = "#3498db"),
            name = "Frequency Class"
          ) +
          geom_vline(xintercept = q90_deltaScore, linetype = "dashed", color = "gray30", size = 0.8) +
          geom_vline(xintercept = q10_deltaScore, linetype = "dashed", color = "gray30", size = 0.8) +
          geom_hline(yintercept = q75_footprint, linetype = "dashed", color = "gray30", size = 0.8) +
          annotate("text", x = q90_deltaScore, y = max(plot_data$RELATIVE_FOOTPRINT_SCORE) * 0.95,
                   label = paste0("Q90 = ", round(q90_deltaScore, 2)),
                   hjust = -0.1, vjust = 1, color = "gray30", size = 4, fontface = "bold") +
          annotate("text", x = q10_deltaScore, y = max(plot_data$RELATIVE_FOOTPRINT_SCORE) * 0.95,
                   label = paste0("Q10 = ", round(q10_deltaScore, 2)),
                   hjust = 1.1, vjust = 1, color = "gray30", size = 4, fontface = "bold") +
          annotate("text", x = min(plot_data$deltaScore) * 0.95, y = q75_footprint,
                   label = paste0("Q75 = ", round(q75_footprint, 3)),
                   hjust = 0, vjust = -0.5, color = "gray30", size = 4, fontface = "bold") +
          labs(
            title = "Relationship between deltaSVM Score and Footprint Score",
            subtitle = paste0("Total SNPs: ", format(nrow(plot_data), big.mark = ","),
                              " | Top ", plot_top_n, " SNPs highlighted"),
            x = "deltaSVM Score",
            y = "Relative Footprint Score"
          ) +
          theme_bw(base_size = 14) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
            legend.position = "top",
            legend.title = element_text(face = "bold"),
            legend.text = element_text(size = 12),
            panel.grid.minor = element_blank(),
            axis.title = element_text(face = "bold"),
            axis.text = element_text(size = 11)
          )
        
        # 保存 Top N 散点图
        top_plot_pdf <- paste0(output_prefix, "_scatter_top", plot_top_n, ".pdf")
        top_plot_png <- paste0(output_prefix, "_scatter_top", plot_top_n, ".png")
        
        ggsave(top_plot_pdf, p_top, width = 12, height = 10, dpi = 300)
        ggsave(top_plot_png, p_top, width = 12, height = 10, dpi = 300)
        
        if (verbose) cat("✅ Top", plot_top_n, "散点图已保存:", top_plot_pdf, "\n")
        
        # 保存 Top N SNPs 信息
        top_snps_file <- paste0(output_prefix, "_top", plot_top_n, "_SNPs.csv")
        write.csv(top_snps, top_snps_file, row.names = FALSE)
        if (verbose) cat("✅ Top", plot_top_n, "SNPs 信息已保存:", top_snps_file, "\n")
      }
    }
  }
  
  # ============================================================================
  # 步骤 7: 保存最终结果
  # ============================================================================
  
  if (verbose) cat("\n步骤 7: 保存最终结果\n")
  if (verbose) cat(rep("-", 80), "\n", sep = "")
  
  output_csv <- paste0(output_prefix, "_results.csv")
  write.csv(TFbSNP_score, output_csv, row.names = FALSE)
  
  if (verbose) cat("✅ 最终结果已保存:", output_csv, "\n")
  
  # ============================================================================
  # 返回结果
  # ============================================================================
  
  if (verbose) {
    cat("\n")
    cat(rep("=", 80), "\n", sep = "")
    cat("分析完成！\n")
    cat(rep("=", 80), "\n\n", sep = "")
  }
  
  result_list <- list(
    data = TFbSNP_score,
    freq_stats = freq_stats,
    input_csv = input_csv,
    genome_version = genome_version,
    flank_size = flank_size,
    ref_fasta = ref_fasta,
    alt_fasta = alt_fasta,
    output_prefix = output_prefix,
    success_count = success_count,
    failed_count = failed_count,
    failed_snps = failed_snps
  )
  
  if (run_gkmsvm && "deltaScore" %in% names(TFbSNP_score)) {
    result_list$plot_data <- plot_data
    result_list$quantiles <- list(
      q10_deltaScore = q10_deltaScore,
      q90_deltaScore = q90_deltaScore,
      q75_footprint = q75_footprint
    )
    if (plot_top_n > 0) {
      result_list$top_snps <- top_snps
    }
  }
  
  invisible(result_list)
}


# ============================================================================
# 辅助函数：快速分析（使用默认参数）
# ============================================================================

#' 快速 TFbSNP 分析
#'
#' @param input_csv 输入 CSV 文件
#' @param svm_model_prefix SVM 模型路径
#' @param output_prefix 输出前缀
#'
#' @export
quick_analyze_TFbSNP <- function(input_csv, svm_model_prefix, output_prefix = "TFbSNP") {
  analyze_TFbSNP(
    input_csv = input_csv,
    svm_model_prefix = svm_model_prefix,
    output_prefix = output_prefix,
    verbose = TRUE
  )
}


# ============================================================================
# 辅助函数：只生成 FASTA（不运行 gkmSVM）
# ============================================================================

#' 只生成 FASTA 序列
#'
#' @param input_csv 输入 CSV 文件
#' @param genome_version 基因组版本
#' @param flank_size 侧翼长度
#' @param ref_fasta 参考序列文件名
#' @param alt_fasta 替代序列文件名
#'
#' @export
generate_fasta_only <- function(
    input_csv,
    genome_version = "hg19",
    flank_size = 9,
    ref_fasta = "reference_alleles.fa",
    alt_fasta = "alt_alleles.fa"
) {
  analyze_TFbSNP(
    input_csv = input_csv,
    genome_version = genome_version,
    flank_size = flank_size,
    ref_fasta = ref_fasta,
    alt_fasta = alt_fasta,
    run_gkmsvm = FALSE,
    verbose = TRUE
  )
}


