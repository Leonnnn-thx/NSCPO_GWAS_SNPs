# ============================================================================
# 检查并去重 FASTA 文件中的重复序列 ID
# ============================================================================

check_and_remove_fasta_duplicates <- function(
    fasta_file, 
    output_file = NULL,
    keep_first = TRUE,
    verbose = TRUE
) {
  
  if (verbose) {
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("检查并去重 FASTA 文件中的重复序列 ID\n")
    cat(rep("=", 80), "\n\n", sep = "")
    cat("输入文件:", fasta_file, "\n")
  }
  
  # 检查文件是否存在
  if (!file.exists(fasta_file)) {
    stop("❌ 文件不存在: ", fasta_file)
  }
  
  # 设置输出文件名
  if (is.null(output_file)) {
    output_file <- paste0(tools::file_path_sans_ext(fasta_file), "_unique.fa")
  }
  
  if (verbose) cat("输出文件:", output_file, "\n\n")
  
  # 读取文件
  lines <- readLines(fasta_file)
  
  # 提取所有的序列 ID（以 > 开头的行）
  header_indices <- grep("^>", lines)
  header_lines <- lines[header_indices]
  
  if (verbose) {
    cat("文件统计:\n")
    cat("  总行数:", length(lines), "\n")
    cat("  序列数:", length(header_lines), "\n\n")
  }
  
  if (length(header_lines) == 0) {
    cat("⚠️  文件中没有找到序列 ID（以 > 开头的行）\n")
    return(invisible(NULL))
  }
  
  # 显示前几个序列 ID
  if (verbose) {
    cat("前10个序列 ID:\n")
    cat(paste(head(header_lines, 10), collapse = "\n"), "\n\n")
  }
  
  # 提取序列 ID（去掉 > 符号）
  seq_ids <- sub("^>", "", header_lines)
  
  # 检查重复
  duplicated_ids <- seq_ids[duplicated(seq_ids)]
  
  if (verbose) {
    cat(rep("-", 80), "\n", sep = "")
    cat("重复检查结果:\n")
    cat(rep("-", 80), "\n", sep = "")
  }
  
  if (length(duplicated_ids) == 0) {
    if (verbose) {
      cat("✅ 没有发现重复的序列 ID\n")
      cat("   文件已经是唯一的，无需去重\n\n")
    }
    
    # 复制原文件到输出文件
    file.copy(fasta_file, output_file, overwrite = TRUE)
    
    if (verbose) cat("✅ 文件已复制到:", output_file, "\n\n")
    
    return(invisible(list(
      file = fasta_file,
      output_file = output_file,
      total_sequences = length(seq_ids),
      unique_sequences = length(unique(seq_ids)),
      duplicates = 0,
      duplicated_ids = character(0),
      removed_count = 0
    )))
    
  } else {
    if (verbose) {
      cat("❌ 发现", length(duplicated_ids), "个重复的序列 ID\n\n")
    }
    
    # 统计每个 ID 出现的次数
    id_counts <- table(seq_ids)
    duplicated_id_counts <- id_counts[id_counts > 1]
    duplicated_id_counts <- sort(duplicated_id_counts, decreasing = TRUE)
    
    if (verbose) {
      cat("重复次数统计:\n")
      cat("  唯一序列 ID 数:", length(unique(seq_ids)), "\n")
      cat("  重复的唯一 ID 数:", length(duplicated_id_counts), "\n")
      cat("  总重复次数:", length(duplicated_ids), "\n\n")
      
      cat("出现次数最多的序列 ID (前20个):\n")
      print(head(duplicated_id_counts, 20))
      cat("\n")
    }
    
    # 显示重复 ID 的详细信息（前5个）
    if (verbose) {
      cat("重复序列 ID 详情 (前5个):\n")
      cat(rep("-", 80), "\n", sep = "")
      
      top_dup_ids <- names(head(duplicated_id_counts, 5))
      
      for (dup_id in top_dup_ids) {
        matching_indices <- which(seq_ids == dup_id)
        
        cat("\n序列 ID:", dup_id, "\n")
        cat("出现次数:", length(matching_indices), "\n")
        cat("序列索引:", paste(matching_indices, collapse = ", "), "\n")
        
        cat("完整 header:\n")
        for (idx in matching_indices) {
          cat("  索引", idx, ":", header_lines[idx], "\n")
        }
      }
      cat("\n")
    }
    
    # ========================================================================
    # 去重：保留第一次出现的序列
    # ========================================================================
    
    if (verbose) {
      cat(rep("=", 80), "\n", sep = "")
      cat("开始去重...\n")
      cat(rep("=", 80), "\n\n", sep = "")
    }
    
    # 找到要保留的序列索引
    if (keep_first) {
      # 保留第一次出现
      keep_indices <- !duplicated(seq_ids)
    } else {
      # 保留最后一次出现
      keep_indices <- !duplicated(seq_ids, fromLast = TRUE)
    }
    
    # 获取要保留的序列的行号
    keep_seq_indices <- which(keep_indices)
    
    if (verbose) {
      cat("去重策略:", if(keep_first) "保留第一次出现" else "保留最后一次出现", "\n")
      cat("保留的序列数:", length(keep_seq_indices), "\n")
      cat("移除的序列数:", length(seq_ids) - length(keep_seq_indices), "\n\n")
    }
    
    # 提取要保留的序列及其序列内容
    unique_lines <- character()
    
    for (i in seq_along(keep_seq_indices)) {
      seq_idx <- keep_seq_indices[i]
      header_line_num <- header_indices[seq_idx]
      
      # 确定序列的结束位置
      if (seq_idx < length(header_indices)) {
        next_header_line_num <- header_indices[seq_idx + 1]
        seq_end_line_num <- next_header_line_num - 1
      } else {
        seq_end_line_num <- length(lines)
      }
      
      # 提取 header 和序列
      seq_lines <- lines[header_line_num:seq_end_line_num]
      unique_lines <- c(unique_lines, seq_lines)
      
      if (verbose && (i %% 1000 == 0 || i == 1)) {
        cat(sprintf("处理进度: %d / %d (%.1f%%)\r", 
                    i, length(keep_seq_indices), 
                    i/length(keep_seq_indices)*100))
      }
    }
    
    if (verbose) cat("\n\n")
    
    # 写入去重后的文件
    writeLines(unique_lines, output_file)
    
    if (verbose) {
      cat("✅ 去重后的文件已保存:", output_file, "\n")
      cat("   原始序列数:", length(seq_ids), "\n")
      cat("   去重后序列数:", length(keep_seq_indices), "\n")
      cat("   移除序列数:", length(seq_ids) - length(keep_seq_indices), "\n\n")
    }
    
    # 保存重复 ID 列表
    dup_file <- paste0(tools::file_path_sans_ext(fasta_file), "_duplicates.txt")
    
    writeLines(
      c(
        paste("重复序列 ID 报告 -", fasta_file),
        paste("生成时间:", Sys.time()),
        "",
        paste("总序列数:", length(seq_ids)),
        paste("唯一序列数:", length(unique(seq_ids))),
        paste("重复的唯一 ID 数:", length(duplicated_id_counts)),
        paste("总重复次数:", length(duplicated_ids)),
        paste("去重后保留:", length(keep_seq_indices)),
        paste("去重后移除:", length(seq_ids) - length(keep_seq_indices)),
        "",
        "重复的序列 ID 及其出现次数:",
        "",
        paste(names(duplicated_id_counts), duplicated_id_counts, sep = "\t")
      ),
      dup_file
    )
    
    if (verbose) cat("✅ 重复 ID 列表已保存到:", dup_file, "\n\n")
    
    # 验证去重结果
    if (verbose) {
      cat(rep("-", 80), "\n", sep = "")
      cat("验证去重结果...\n")
      cat(rep("-", 80), "\n", sep = "")
    }
    
    # 读取去重后的文件
    unique_lines_check <- readLines(output_file)
    unique_headers <- grep("^>", unique_lines_check, value = TRUE)
    unique_ids <- sub("^>", "", unique_headers)
    
    if (verbose) {
      cat("去重后文件统计:\n")
      cat("  总行数:", length(unique_lines_check), "\n")
      cat("  序列数:", length(unique_headers), "\n")
      cat("  唯一 ID 数:", length(unique(unique_ids)), "\n")
      
      if (length(unique_ids) == length(unique(unique_ids))) {
        cat("\n✅ 验证通过：去重后的文件没有重复 ID\n\n")
      } else {
        cat("\n⚠️  警告：去重后的文件仍有重复 ID\n\n")
      }
    }
    
    return(invisible(list(
      file = fasta_file,
      output_file = output_file,
      total_sequences = length(seq_ids),
      unique_sequences = length(unique(seq_ids)),
      duplicates = length(duplicated_ids),
      duplicated_ids = names(duplicated_id_counts),
      duplicated_counts = duplicated_id_counts,
      duplicate_file = dup_file,
      removed_count = length(seq_ids) - length(keep_seq_indices),
      kept_count = length(keep_seq_indices)
    )))
  }
}


# ============================================================================
# 批量处理多个 FASTA 文件
# ============================================================================

batch_remove_fasta_duplicates <- function(
    fasta_files,
    output_suffix = "_unique",
    keep_first = TRUE,
    verbose = TRUE
) {
  
  results <- list()
  
  for (fasta_file in fasta_files) {
    
    if (verbose) {
      cat("\n", rep("=", 80), "\n", sep = "")
      cat("处理文件:", fasta_file, "\n")
      cat(rep("=", 80), "\n", sep = "")
    }
    
    # 生成输出文件名
    base_name <- tools::file_path_sans_ext(fasta_file)
    ext <- tools::file_ext(fasta_file)
    output_file <- paste0(base_name, output_suffix, ".", ext)
    
    # 处理文件
    result <- check_and_remove_fasta_duplicates(
      fasta_file = fasta_file,
      output_file = output_file,
      keep_first = keep_first,
      verbose = verbose
    )
    
    results[[fasta_file]] <- result
  }
  
  # 打印总结
  if (verbose) {
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("批量处理完成\n")
    cat(rep("=", 80), "\n\n", sep = "")
    
    cat("处理总结:\n")
    for (file in names(results)) {
      result <- results[[file]]
      cat(sprintf("  %s:\n", basename(file)))
      cat(sprintf("    原始序列: %d\n", result$total_sequences))
      cat(sprintf("    唯一序列: %d\n", result$unique_sequences))
      cat(sprintf("    移除重复: %d\n", result$removed_count))
      cat(sprintf("    输出文件: %s\n\n", basename(result$output_file)))
    }
  }
  
  invisible(results)
}



result_ref <- check_and_remove_fasta_duplicates(
  fasta_file = "control_reference_alleles.fa",
  output_file = "control_reference_alleles_unique.fa",
  keep_first = TRUE,
  verbose = TRUE
)


result_alt <- check_and_remove_fasta_duplicates(
  fasta_file = "control_alt_alleles.fa",
  output_file = "control_alt_alleles_unique.fa",
  keep_first = TRUE,
  verbose = TRUE
)
