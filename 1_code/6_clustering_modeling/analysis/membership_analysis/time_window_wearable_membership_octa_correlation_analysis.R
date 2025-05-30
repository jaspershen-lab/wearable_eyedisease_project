# 时间窗口特异性Membership分析
# Time Window Specific Membership vs OCTA Analysis
# 分析不同时间窗口数据产生的membership与OCTA改善的关系

library(tidyverse)
library(Biobase)
library(Mfuzz)
library(ggplot2)
library(gridExtra)

# ================== 1. 为每个时间窗口计算独立的membership ==================
cat("===== 时间窗口特异性Membership分析 =====\n")

# 定义时间窗口（与之前保持一致）
time_windows <- list(
  baseline = list(days = -4:-1, name = "baseline"),
  acute_recovery = list(days = 0:3, name = "acute_recovery"),
  early_recovery = list(days = 4:7, name = "early_recovery"),
  mid_recovery = list(days = 8:15, name = "mid_recovery"),
  late_recovery = list(days = 16:30, name = "late_recovery")
)

# 加载必要数据
ppv_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)

# 关键指标
key_metrics <- c("cv_rhr_1", "steps_max")

# 函数：为单个时间窗口计算membership
calculate_window_membership <- function(data, window_info, metrics) {
  window_name <- window_info$name
  window_days <- window_info$days
  
  cat(sprintf("计算 %s 时间窗口的membership...\n", window_name))
  
  # 提取该时间窗口的数据
  window_cols <- c()
  for(metric in metrics) {
    for(day in window_days) {
      day_str <- paste0("day_", day, "_", metric)
      if(day_str %in% colnames(data)) {
        window_cols <- c(window_cols, day_str)
      }
    }
  }
  
  if(length(window_cols) == 0) {
    cat(sprintf("警告: %s 时间窗口没有可用数据\n", window_name))
    return(NULL)
  }
  
  # 选择数据并计算每个患者在该时间窗口的均值
  window_data <- data %>%
    dplyr::select(subject_id, all_of(window_cols))
  
  # 为每个指标计算时间窗口均值
  processed_data <- data %>% dplyr::select(subject_id)
  
  for(metric in metrics) {
    metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
    if(length(metric_cols) > 0) {
      metric_means <- window_data %>%
        dplyr::select(subject_id, all_of(metric_cols)) %>%
        mutate(
          valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
          metric_mean = ifelse(
            valid_count >= max(1, floor(length(metric_cols)/2)),
            rowMeans(dplyr::select(., -subject_id), na.rm = TRUE),
            NA
          )
        ) %>%
        dplyr::select(subject_id, metric_mean)
      
      names(metric_means)[2] <- paste0(window_name, "_", metric)
      processed_data <- processed_data %>%
        left_join(metric_means, by = "subject_id")
    }
  }
  
  # 移除有太多NA的患者
  complete_patients <- processed_data %>%
    filter(rowSums(is.na(dplyr::select(., -subject_id))) < ncol(dplyr::select(., -subject_id)))
  
  if(nrow(complete_patients) < 5) {
    cat(sprintf("警告: %s 时间窗口有效患者数不足 (%d)\n", window_name, nrow(complete_patients)))
    return(NULL)
  }
  
  # 用均值填充剩余的NA值
  numeric_cols <- names(complete_patients)[-1]
  for(col in numeric_cols) {
    if(sum(!is.na(complete_patients[[col]])) > 0) {
      complete_patients[is.na(complete_patients[[col]]), col] <- 
        mean(complete_patients[[col]], na.rm = TRUE)
    }
  }
  
  # 标准化数据
  scaled_data <- complete_patients
  for(col in numeric_cols) {
    scaled_data[[col]] <- scale(complete_patients[[col]])[,1]
  }
  
  # 准备Mfuzz数据
  data_matrix <- scaled_data %>%
    dplyr::select(-subject_id) %>%
    as.matrix()
  
  rownames(data_matrix) <- scaled_data$subject_id
  
  # 创建ExpressionSet
  eset <- ExpressionSet(assayData = data_matrix)
  eset_std <- standardise(eset)
  
  # 估计最佳参数
  m_value <- mestimate(eset_std)
  
  # 执行聚类 (使用2-3个聚类)
  optimal_c <- min(3, max(2, floor(nrow(complete_patients)/4)))
  
  set.seed(123)
  clustering_result <- mfuzz(eset_std, c = optimal_c, m = m_value)
  
  # 提取membership信息
  main_clusters <- apply(clustering_result$membership, 1, which.max)
  max_memberships <- apply(clustering_result$membership, 1, max)
  
  # 创建结果数据框
  membership_result <- data.frame(
    subject_id = rownames(clustering_result$membership),
    window = window_name,
    max_cluster = main_clusters,
    max_membership = max_memberships,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("%s 聚类完成: %d 患者, %d 聚类, 平均membership = %.3f\n", 
              window_name, nrow(membership_result), optimal_c, mean(max_memberships)))
  
  return(list(
    membership_data = membership_result,
    clustering_result = clustering_result,
    window_name = window_name,
    n_patients = nrow(complete_patients),
    n_clusters = optimal_c
  ))
}

# ================== 2. 为每个时间窗口计算membership ==================
window_memberships <- list()
all_membership_data <- data.frame()

for(window_name in names(time_windows)) {
  window_result <- calculate_window_membership(ppv_data, time_windows[[window_name]], key_metrics)
  
  if(!is.null(window_result)) {
    window_memberships[[window_name]] <- window_result
    all_membership_data <- rbind(all_membership_data, window_result$membership_data)
  }
}

cat(sprintf("\n成功计算了 %d 个时间窗口的membership\n", length(window_memberships)))

# ================== 3. 重塑数据用于分析 ==================
# 将长格式转换为宽格式，每个时间窗口一列membership
membership_wide <- all_membership_data %>%
  dplyr::select(subject_id, window, max_membership) %>%
  pivot_wider(
    names_from = window,
    values_from = max_membership,
    names_prefix = "membership_"
  )

cat("\n时间窗口membership数据结构:\n")
print(names(membership_wide))
print(head(membership_wide))

# ================== 4. 合并OCTA数据 ==================
# 重新加载和处理OCTA数据（使用之前的方法）
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# 处理OCTA数据（重用之前的函数）
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

# 筛选和计算OCTA改善参数（重用之前的函数）
filter_key_octa_params <- function(data, param_type = "bloodflow") {
  if(param_type == "bloodflow") {
    layers <- c("SVP", "ICP", "DCP", "Choroid")
  } else {
    layers <- c("GCL.IPL", "INL", "Retina")
  }
  
  regions <- c("0_21", "0_6")
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

# 处理OCTA数据
ppv_bloodflow <- process_octa_improvements(baseline_info, octa_bloodflow)
ppv_thickness <- process_octa_improvements(baseline_info, octa_thickness)

bloodflow_filtered <- filter_key_octa_params(ppv_bloodflow, "bloodflow")
thickness_filtered <- filter_key_octa_params(ppv_thickness, "thickness")

bloodflow_improvements <- calculate_improvement(
  ppv_bloodflow %>% dplyr::select(ID, all_of(c(bloodflow_filtered$params_T0, bloodflow_filtered$params_T2))),
  bloodflow_filtered$params_T0, bloodflow_filtered$params_T2
)

thickness_improvements <- calculate_improvement(
  ppv_thickness %>% dplyr::select(ID, all_of(c(thickness_filtered$params_T0, thickness_filtered$params_T2))),
  thickness_filtered$params_T0, thickness_filtered$params_T2
)

octa_improvements <- bloodflow_improvements %>%
  full_join(thickness_improvements, by = "ID")

# ================== 5. 时间窗口特异性membership相关分析 ==================
# 合并membership和OCTA数据
window_membership_analysis <- membership_wide %>%
  left_join(octa_improvements, by = c("subject_id" = "ID"))

cat("\n时间窗口membership分析数据:\n")
cat("总患者:", nrow(window_membership_analysis), "\n")

# 检查每个时间窗口的有效membership数量
membership_cols <- grep("^membership_", names(window_membership_analysis), value = TRUE)
for(col in membership_cols) {
  valid_count <- sum(!is.na(window_membership_analysis[[col]]))
  cat(paste0(col, ": ", valid_count, " 有效值\n"))
}

# 获取OCTA参数
octa_improvement_params <- names(octa_improvements)[grep("_improvement$", names(octa_improvements))]

# 执行时间窗口特异性membership相关分析
perform_window_membership_correlation <- function(data, membership_cols, octa_params) {
  results <- data.frame()
  
  for(membership_col in membership_cols) {
    window_name <- gsub("^membership_", "", membership_col)
    cat(sprintf("\n分析 %s 时间窗口的membership...\n", window_name))
    
    for(octa_param in octa_params) {
      if(!octa_param %in% names(data)) next
      
      # 创建完整案例数据
      complete_data <- data[!is.na(data[[membership_col]]) & !is.na(data[[octa_param]]), ]
      
      if(nrow(complete_data) >= 3) {
        # Pearson相关
        cor_test <- try(cor.test(complete_data[[membership_col]], complete_data[[octa_param]], 
                                 method = "pearson"), silent = TRUE)
        
        # Spearman相关
        spearman_test <- try(cor.test(complete_data[[membership_col]], complete_data[[octa_param]], 
                                      method = "spearman"), silent = TRUE)
        
        if(class(cor_test) != "try-error" && class(spearman_test) != "try-error") {
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
          
          effect_size <- case_when(
            abs(cor_test$estimate) >= 0.5 ~ "Large",
            abs(cor_test$estimate) >= 0.3 ~ "Medium",
            abs(cor_test$estimate) >= 0.1 ~ "Small",
            TRUE ~ "Negligible"
          )
          
          results <- rbind(results, data.frame(
            Time_Window = window_name,
            OCTA_Parameter = octa_param,
            Parameter_Type = param_type,
            Region = region,
            N = nrow(complete_data),
            Pearson_r = as.numeric(cor_test$estimate),
            Pearson_p = cor_test$p.value,
            CI_Lower = cor_test$conf.int[1],
            CI_Upper = cor_test$conf.int[2],
            Spearman_rho = as.numeric(spearman_test$estimate),
            Spearman_p = spearman_test$p.value,
            Effect_Size = effect_size,
            Significant_p05 = cor_test$p.value < 0.05,
            Trend_p10 = cor_test$p.value < 0.10,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # FDR校正
  if(nrow(results) > 0) {
    results$Pearson_p_FDR <- p.adjust(results$Pearson_p, method = "fdr")
    results$Spearman_p_FDR <- p.adjust(results$Spearman_p, method = "fdr")
    results$Significant_FDR <- results$Pearson_p_FDR < 0.05
    
    # 按相关性强度排序
    results <- results %>% arrange(desc(abs(Pearson_r)))
  }
  
  return(results)
}

# 执行时间窗口membership相关分析
window_membership_correlations <- perform_window_membership_correlation(
  window_membership_analysis, 
  membership_cols, 
  octa_improvement_params
)

# ================== 6. 显示结果 ==================
cat("\n===== 时间窗口特异性Membership与OCTA相关分析结果 =====\n")

if(nrow(window_membership_correlations) > 0) {
  # 显著结果
  significant_results <- window_membership_correlations %>%
    filter(Significant_p05 == TRUE) %>%
    arrange(desc(abs(Pearson_r)))
  
  if(nrow(significant_results) > 0) {
    cat("🎯 显著相关结果 (p < 0.05):\n")
    print(significant_results %>%
            dplyr::select(Time_Window, OCTA_Parameter, Parameter_Type, Region, 
                          Pearson_r, Pearson_p, Effect_Size, N))
  } else {
    cat("❌ 未发现显著相关 (p < 0.05)\n")
  }
  
  # 趋势性结果
  trend_results <- window_membership_correlations %>%
    filter(Trend_p10 == TRUE & abs(Pearson_r) >= 0.4) %>%
    arrange(Pearson_p)
  
  if(nrow(trend_results) > 0) {
    cat("\n📈 趋势性显著结果 (p < 0.10, |r| ≥ 0.4):\n")
    print(trend_results %>%
            dplyr::select(Time_Window, OCTA_Parameter, Parameter_Type, Region, 
                          Pearson_r, Pearson_p, Effect_Size, N))
  }
  
  # 按时间窗口汇总
  window_summary <- window_membership_correlations %>%
    group_by(Time_Window) %>%
    summarise(
      Total_Tests = n(),
      Significant_p05 = sum(Significant_p05),
      Trend_p10 = sum(Trend_p10),
      Large_Effect = sum(Effect_Size == "Large"),
      Mean_Abs_Correlation = round(mean(abs(Pearson_r)), 3),
      Best_Correlation = round(max(abs(Pearson_r)), 3),
      .groups = 'drop'
    ) %>%
    arrange(desc(Significant_p05), desc(Trend_p10))
  
  cat("\n📊 各时间窗口membership预测能力:\n")
  print(window_summary)
  
  # 最强相关性
  if(nrow(window_membership_correlations) > 0) {
    top_result <- window_membership_correlations[1, ]
    cat(sprintf("\n🏆 最强相关性:\n"))
    cat(sprintf("时间窗口: %s\n", top_result$Time_Window))
    cat(sprintf("OCTA参数: %s\n", top_result$OCTA_Parameter))
    cat(sprintf("相关系数: r = %.3f (p = %.4f)\n", top_result$Pearson_r, top_result$Pearson_p))
    cat(sprintf("效应量: %s\n", top_result$Effect_Size))
  }
  
} else {
  cat("❌ 未找到任何有效的相关性结果\n")
}

# ================== 7. 保存结果 ==================
write.csv(window_membership_correlations, "time_window_membership_correlations.csv", row.names = FALSE)
write.csv(membership_wide, "time_window_membership_data.csv", row.names = FALSE)

# 保存聚类信息
clustering_summary <- data.frame(
  Window = names(window_memberships),
  N_Patients = sapply(window_memberships, function(x) x$n_patients),
  N_Clusters = sapply(window_memberships, function(x) x$n_clusters),
  Mean_Membership = sapply(window_memberships, function(x) mean(x$membership_data$max_membership)),
  stringsAsFactors = FALSE
)

write.csv(clustering_summary, "time_window_clustering_summary.csv", row.names = FALSE)

cat("\n===== 分析完成 =====\n")
cat("这个分析计算了每个时间窗口独立的membership值\n")
cat("然后分析这些时间特异性membership与OCTA改善的关系\n")
cat("结果文件:\n")
cat("- time_window_membership_correlations.csv: 详细相关性结果\n")
cat("- time_window_membership_data.csv: 时间窗口membership数据\n")
cat("- time_window_clustering_summary.csv: 聚类信息汇总\n")

