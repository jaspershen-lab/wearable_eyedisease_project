# 时间窗口特异性相关分析
# Time Window Specific Correlation Analysis
# 探索不同时间窗口的可穿戴设备数据与OCTA改善的关系

library(tidyverse)
library(corrplot)
library(ggplot2)
library(gridExtra)

# 设置工作目录
setwd(get_project_wd())

# ================== 1. 加载和准备数据 ==================
cat("===== 时间窗口特异性分析 =====\n")

# 加载可穿戴设备原始数据（包含所有时间窗口信息）
ppv_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)

# 加载OCTA改善数据（从之前的分析中获取）
# 需要重新处理OCTA数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# 创建输出目录
dir.create("3_data_analysis/6_clustering_modeling/time_window_correlation", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/time_window_correlation")

# ================== 2. 重新处理OCTA数据 ==================
# 使用之前定义的函数处理OCTA数据
process_octa_improvements <- function(baseline_data, octa_data, id_column = "id") {
  ppv_patients <- baseline_data %>%
    filter(surgery_1..0.PI.1.other. == 1) %>%
    distinct(ID) %>%
    pull(ID)
  
  octa_features <- baseline_data %>%
    filter(ID %in% ppv_patients & !is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  process_patient_octa <- function(patient_data, time_points = c("T0", "T2")) {
    current_eye <- patient_data$surgery_eye_1[1]
    pattern <- if(current_eye == 1) "_OS_" else "_OD_"
    
    result <- patient_data %>% dplyr::select(ID)
    
    for(suffix in time_points) {
      cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
      cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
      
      if(length(cols_to_keep) > 0) {
        time_data <- patient_data %>%
          dplyr::select("ID", all_of(cols_to_keep)) %>%
          rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
        
        result <- result %>% left_join(time_data, by = "ID")
      }
    }
    
    return(result)
  }
  
  patient_list <- split(octa_features, octa_features$ID)
  processed_data <- map_dfr(patient_list, process_patient_octa)
  
  return(processed_data)
}

# 筛选关键OCTA参数
filter_key_octa_params <- function(data, param_type = "bloodflow") {
  if(param_type == "bloodflow") {
    layers <- c("SVP", "ICP", "DCP", "Choroid")
  } else {
    layers <- c("GCL.IPL", "INL", "Retina")
  }
  
  regions <- c("0_21", "0_6")  # 黄斑区和广角区
  
  pattern <- paste0("(", paste(layers, collapse = "|"), ").*(",
                    paste(regions, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  return(list(
    base_params = valid_base_params,
    params_T0 = paste0(valid_base_params, "_T0"),
    params_T2 = paste0(valid_base_params, "_T2")
  ))
}

# 处理OCTA数据
ppv_bloodflow <- process_octa_improvements(baseline_info, octa_bloodflow)
ppv_thickness <- process_octa_improvements(baseline_info, octa_thickness)

# 筛选参数
bloodflow_filtered <- filter_key_octa_params(ppv_bloodflow, "bloodflow")
thickness_filtered <- filter_key_octa_params(ppv_thickness, "thickness")

# 计算改善值
calculate_improvement <- function(data, params_T0, params_T2) {
  result <- data %>% dplyr::select(ID)
  
  for(i in 1:length(params_T0)) {
    t0_param <- params_T0[i]
    t2_param <- params_T2[i]
    base_param <- gsub("_T0$", "", t0_param)
    
    result[[paste0(base_param, "_improvement")]] <- data[[t2_param]] - data[[t0_param]]
  }
  
  return(result)
}

bloodflow_improvements <- calculate_improvement(
  ppv_bloodflow %>% dplyr::select(ID, all_of(c(bloodflow_filtered$params_T0, bloodflow_filtered$params_T2))),
  bloodflow_filtered$params_T0, bloodflow_filtered$params_T2
)

thickness_improvements <- calculate_improvement(
  ppv_thickness %>% dplyr::select(ID, all_of(c(thickness_filtered$params_T0, thickness_filtered$params_T2))),
  thickness_filtered$params_T0, thickness_filtered$params_T2
)

# 合并OCTA改善数据
octa_improvements <- bloodflow_improvements %>%
  full_join(thickness_improvements, by = "ID")

# ================== 3. 提取时间窗口数据 ==================
# 定义时间窗口（与你的聚类代码一致）
time_windows <- list(
  baseline = list(days = -4:-1, name = "baseline"),
  acute_recovery = list(days = 0:3, name = "acute_recovery"),
  early_recovery = list(days = 4:7, name = "early_recovery"),
  mid_recovery = list(days = 8:15, name = "mid_recovery"),
  late_recovery = list(days = 16:30, name = "late_recovery")
)

# 提取时间窗口数据的函数
extract_time_window_data <- function(data, metrics, time_windows) {
  result_data <- data %>% dplyr::select(subject_id)
  
  for(metric in metrics) {
    cat(sprintf("处理指标: %s\n", metric))
    
    for(window_name in names(time_windows)) {
      window <- time_windows[[window_name]]
      
      # 收集该时间窗口内的所有列
      window_cols <- c()
      for(day in window$days) {
        day_str <- paste0("day_", day, "_", metric)
        if(day_str %in% colnames(data)) {
          window_cols <- c(window_cols, day_str)
        }
      }
      
      if(length(window_cols) > 0) {
        # 计算时间窗口均值
        window_data <- data %>%
          dplyr::select(subject_id, all_of(window_cols))
        
        min_valid_points <- max(1, floor(length(window_cols) / 2))
        
        window_mean <- window_data %>%
          mutate(
            valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
            window_mean = ifelse(
              valid_count >= min_valid_points,
              rowMeans(dplyr::select(., -subject_id), na.rm = TRUE),
              NA
            )
          ) %>%
          dplyr::select(subject_id, window_mean)
        
        colnames(window_mean)[2] <- paste0(window$name, "_", metric)
        result_data <- result_data %>%
          left_join(window_mean, by = "subject_id")
      }
    }
  }
  
  return(result_data)
}

# 关键可穿戴设备指标
key_metrics <- c("cv_rhr_1", "steps_max")

# 提取时间窗口数据
wearable_windows <- extract_time_window_data(ppv_data, key_metrics, time_windows)

# ================== 4. 时间窗口特异性相关分析 ==================
# 合并数据
analysis_data <- wearable_windows %>%
  left_join(octa_improvements, by = c("subject_id" = "ID"))

cat("时间窗口分析数据:", nrow(analysis_data), "patients\n")

# 时间窗口特异性相关分析函数
perform_window_specific_correlation <- function(data, window_names, metrics, octa_params) {
  results <- data.frame()
  
  for(window in window_names) {
    cat(sprintf("\n分析时间窗口: %s\n", window))
    
    for(metric in metrics) {
      window_metric_col <- paste0(window, "_", metric)
      
      if(!window_metric_col %in% names(data)) {
        cat(sprintf("跳过: %s (列不存在)\n", window_metric_col))
        next
      }
      
      for(octa_param in octa_params) {
        if(!octa_param %in% names(data)) next
        
        # 创建完整案例数据
        complete_data <- data[!is.na(data[[window_metric_col]]) & !is.na(data[[octa_param]]), ]
        
        if(nrow(complete_data) >= 3) {
          # Pearson相关
          cor_test <- try(cor.test(complete_data[[window_metric_col]], complete_data[[octa_param]], 
                                   method = "pearson"), silent = TRUE)
          
          if(class(cor_test) != "try-error") {
            # 参数分类
            param_type <- case_when(
              grepl("SVP|ICP|DCP|Choroid", octa_param) ~ "BloodFlow",
              grepl("GCL|INL|Retina", octa_param) ~ "Thickness",
              TRUE ~ "Other"
            )
            
            region <- case_when(
              grepl("0_21", octa_param) ~ "Macular",
              grepl("0_6", octa_param) ~ "Widefield",
              TRUE ~ "Other"
            )
            
            results <- rbind(results, data.frame(
              Time_Window = window,
              Wearable_Metric = metric,
              OCTA_Parameter = octa_param,
              Parameter_Type = param_type,
              Region = region,
              N = nrow(complete_data),
              Correlation = as.numeric(cor_test$estimate),
              P_Value = cor_test$p.value,
              CI_Lower = cor_test$conf.int[1],
              CI_Upper = cor_test$conf.int[2],
              Significant = cor_test$p.value < 0.05,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  # FDR校正
  if(nrow(results) > 0) {
    results$P_FDR <- p.adjust(results$P_Value, method = "fdr")
    results$Significant_FDR <- results$P_FDR < 0.05
    results <- results %>% arrange(P_Value)
  }
  
  return(results)
}

# 获取所有OCTA改善参数
octa_improvement_params <- names(octa_improvements)[grep("_improvement$", names(octa_improvements))]

# 执行时间窗口特异性分析
window_correlations <- perform_window_specific_correlation(
  analysis_data, 
  names(time_windows), 
  key_metrics, 
  octa_improvement_params
)

# 显示结果
cat("\n===== 时间窗口特异性相关分析结果 =====\n")
if(nrow(window_correlations) > 0) {
  # 显著结果
  significant_results <- window_correlations %>%
    filter(Significant == TRUE) %>%
    arrange(desc(abs(Correlation)))
  
  if(nrow(significant_results) > 0) {
    cat("显著相关结果 (p < 0.05):\n")
    print(significant_results %>%
            dplyr::select(Time_Window, Wearable_Metric, OCTA_Parameter, 
                          Parameter_Type, Region, Correlation, P_Value, N))
  } else {
    cat("未发现传统意义上的显著相关\n")
  }
  
  # 趋势性显著结果
  trend_results <- window_correlations %>%
    filter(P_Value < 0.10 & abs(Correlation) >= 0.4) %>%
    arrange(P_Value)
  
  if(nrow(trend_results) > 0) {
    cat("\n\n趋势性显著结果 (p < 0.10, |r| ≥ 0.4):\n")
    print(trend_results %>%
            dplyr::select(Time_Window, Wearable_Metric, OCTA_Parameter, 
                          Parameter_Type, Region, Correlation, P_Value, N))
  }
  
  # 按时间窗口总结
  window_summary <- window_correlations %>%
    group_by(Time_Window) %>%
    summarise(
      Total_Tests = n(),
      Significant_p05 = sum(Significant),
      Significant_p10 = sum(P_Value < 0.10),
      Strong_Correlation = sum(abs(Correlation) >= 0.5),
      Mean_Abs_Correlation = round(mean(abs(Correlation)), 3),
      .groups = 'drop'
    ) %>%
    arrange(desc(Significant_p05), desc(Significant_p10))
  
  cat("\n\n各时间窗口预测能力总结:\n")
  print(window_summary)
} else {
  cat("未找到有效的相关性结果\n")
}

# 保存详细结果
write.csv(window_correlations, "time_window_specific_correlations.csv", row.names = FALSE)

# ================== 5. 可视化时间窗口特异性结果 ==================
# 创建热图：时间窗口 vs OCTA参数的相关性
create_window_correlation_heatmap <- function(corr_results) {
  if(nrow(corr_results) == 0) return(NULL)
  
  # 选择最强的相关性进行可视化
  top_results <- corr_results %>%
    group_by(Time_Window, OCTA_Parameter) %>%
    slice_max(abs(Correlation), n = 1) %>%
    ungroup() %>%
    mutate(
      OCTA_Clean = gsub("_improvement", "", OCTA_Parameter),
      OCTA_Clean = gsub("_", " ", OCTA_Clean),
      OCTA_Display = case_when(
        grepl("0_21", OCTA_Parameter) ~ paste0(OCTA_Clean, " (Mac)"),
        grepl("0_6", OCTA_Parameter) ~ paste0(OCTA_Clean, " (WF)"),
        TRUE ~ OCTA_Clean
      ),
      Significance_Label = case_when(
        P_Value < 0.001 ~ "***",
        P_Value < 0.01 ~ "**",
        P_Value < 0.05 ~ "*",
        P_Value < 0.10 ~ ".",
        TRUE ~ ""
      )
    )
  
  p <- ggplot(top_results, aes(x = Time_Window, y = OCTA_Display, fill = Correlation)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = paste0(sprintf("%.2f", Correlation), "\n", Significance_Label)),
              color = "black", size = 3) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = "Correlation",
      limits = c(-1, 1)
    ) +
    labs(
      title = "Time Window Specific Correlations\nWearable Metrics vs OCTA Improvements",
      x = "Time Window",
      y = "OCTA Parameters",
      caption = "*** p<0.001, ** p<0.01, * p<0.05, . p<0.10"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

# 创建时间窗口特异性热图
window_heatmap <- create_window_correlation_heatmap(window_correlations)
if(!is.null(window_heatmap)) {
  ggsave("time_window_correlation_heatmap.pdf", window_heatmap, width = 12, height = 10)
  print(window_heatmap)
}

# ================== 6. 临界时间窗口分析 ==================
# 识别最具预测价值的时间窗口
identify_critical_windows <- function(corr_results) {
  if(nrow(corr_results) == 0) return(NULL)
  
  # 计算每个时间窗口的预测价值得分
  window_scores <- corr_results %>%
    group_by(Time_Window) %>%
    summarise(
      Significant_Count = sum(Significant),
      Trend_Count = sum(P_Value < 0.10),
      Strong_Corr_Count = sum(abs(Correlation) >= 0.5),
      Mean_Abs_Corr = mean(abs(Correlation)),
      Best_Correlation = max(abs(Correlation)),
      # 综合评分
      Prediction_Score = Significant_Count * 3 + Trend_Count * 2 + Strong_Corr_Count * 1 + Mean_Abs_Corr,
      .groups = 'drop'
    ) %>%
    arrange(desc(Prediction_Score))
  
  cat("\n===== 临界时间窗口识别 =====\n")
  cat("各时间窗口预测价值排序:\n")
  print(window_scores)
  
  # 识别最佳时间窗口
  best_window <- window_scores$Time_Window[1]
  cat(sprintf("\n最具预测价值的时间窗口: %s\n", best_window))
  
  # 分析该时间窗口的最佳预测组合
  best_window_results <- corr_results %>%
    filter(Time_Window == best_window) %>%
    arrange(desc(abs(Correlation))) %>%
    head(10)
  
  cat(sprintf("\n%s 时间窗口的最佳预测组合:\n", best_window))
  print(best_window_results %>%
          dplyr::select(Wearable_Metric, OCTA_Parameter, Correlation, P_Value, N))
  
  return(list(
    window_scores = window_scores,
    best_window = best_window,
    best_combinations = best_window_results
  ))
}

critical_windows <- identify_critical_windows(window_correlations)

# ================== 7. 生成改进建议报告 ==================
generate_improvement_report <- function(window_results, critical_results) {
  total_tests <- nrow(window_results)
  significant_count <- sum(window_results$Significant)
  trend_count <- sum(window_results$P_Value < 0.10)
  
  best_window <- if(!is.null(critical_results)) critical_results$best_window else "未确定"
  
  report <- paste0(
    "========================================\n",
    "时间窗口特异性分析改进建议报告\n",
    "========================================\n\n",
    
    "分析概况:\n",
    "- 总相关性测试: ", total_tests, " 个\n",
    "- 显著相关 (p<0.05): ", significant_count, " 个 (", 
    round(significant_count/total_tests*100, 1), "%)\n",
    "- 趋势显著 (p<0.10): ", trend_count, " 个 (", 
    round(trend_count/total_tests*100, 1), "%)\n",
    "- 最佳预测时间窗口: ", best_window, "\n\n",
    
    "改进策略建议:\n\n",
    
    "1. 聚焦最佳时间窗口分析:\n",
    "   - 重点分析 ", best_window, " 时间窗口\n",
    "   - 该窗口显示最强的预测能力\n",
    "   - 可以单独针对此窗口建立预测模型\n\n",
    
    "2. 样本量优化策略:\n",
    "   - 当前样本量较小，影响统计功效\n",
    "   - 建议增加样本至25-30例\n",
    "   - 或考虑meta分析方法整合多中心数据\n\n",
    
    "3. 多变量分析方法:\n",
    "   - 使用主成分分析降维\n",
    "   - 建立多元回归模型\n",
    "   - 考虑机器学习方法（随机森林、支持向量机）\n\n",
    
    "4. 临床意义解释:\n",
    "   - 重视效应量 (|r|≥0.5) 而非仅看p值\n",
    "   - 趋势性结果 (p<0.10) 在探索性研究中有价值\n",
    "   - 结合临床专业知识解释时间窗口特异性\n\n",
    
    "5. 分析方法优化:\n",
    "   - 考虑单侧检验（如有方向性假设）\n",
    "   - 使用Bootstrap方法增强统计稳健性\n",
    "   - 交叉验证评估模型泛化能力\n\n",
    
    "6. 数据预处理改进:\n",
    "   - 异常值检测和处理\n",
    "   - 考虑非线性关系（Spearman相关）\n",
    "   - 时间序列平滑处理\n\n"
  )
  
  if(trend_count > significant_count) {
    report <- paste0(report,
                     "🎯 当前发现:\n",
                     "虽然传统显著性结果有限，但发现了多个趋势性关联，\n",
                     "这在探索性研究中具有重要价值。建议:\n",
                     "- 扩大样本量验证这些趋势\n",
                     "- 重点关注效应量大的关联\n",
                     "- 考虑临床相关性而非仅统计显著性\n\n"
    )
  }
  
  writeLines(report, "Time_Window_Analysis_Improvement_Report.txt")
  cat(report)
}

generate_improvement_report(window_correlations, critical_windows)

# 保存工作空间
save(analysis_data, window_correlations, critical_windows,
     file = "time_window_correlation_analysis.RData")

cat("\n========================================\n")
cat("时间窗口特异性分析完成！\n")
cat("========================================\n")
cat("关键发现:\n")
if(sum(window_correlations$Significant) > 0) {
  cat("✓ 发现显著的时间窗口特异性关联\n")
} else {
  cat("• 未发现传统显著关联，但可能存在趋势性关系\n")
}
cat("✓ 识别了最具预测价值的时间窗口\n")
cat("✓ 提供了详细的改进策略建议\n")
cat("\n请查看生成的报告和可视化结果！\n")
