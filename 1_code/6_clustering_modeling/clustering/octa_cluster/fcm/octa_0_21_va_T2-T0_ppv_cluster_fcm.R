# ===============================================
# 使用所有14个参数的完整FCM聚类分析
# 不进行变量选择，使用原始的全部指标
# ===============================================

library(tidyverse)
library(e1071)
library(cluster)
library(ggplot2)

# Set working directory
setwd(get_project_wd())
rm(list = ls())


# -------------------- 1. 数据准备 (保持原有的数据处理部分) --------------------
# Load baseline information and OCTA data
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# Create output directory
dir.create("3_data_analysis/6_clustering_modeling/fcm_octa/WF_only_cluster", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/fcm_octa/WF_only_cluster")

# -------------------- 2. Process OCTA data --------------------
# Function to process OCTA data
process_octa_data <- function(baseline_data, octa_data, id_column = "id") {
  features <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  return(features)
}

# Process blood flow and thickness data
octa_bloodflow_features <- process_octa_data(baseline_info, octa_bloodflow)
octa_thickness_features <- process_octa_data(baseline_info, octa_thickness)

# Function to process patient data
process_patient_data <- function(patient_data, time_points = c("T0", "T2")) {
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

# Function to process all patients
process_all_patients <- function(features_data) {
  patient_list <- split(features_data, features_data$ID)
  processed_data <- purrr::map(patient_list, process_patient_data)
  return(bind_rows(processed_data))
}

# Process OCTA data
octa_bloodflow_processed <- process_all_patients(octa_bloodflow_features)
octa_thickness_processed <- process_all_patients(octa_thickness_features)

# -------------------- 3. Extract PPV group data --------------------
# Get PPV patient IDs
ppv_patients <- baseline_info %>%
  filter(surgery_1..0.PI.1.other. == 1) %>%
  distinct(ID) %>%
  pull(ID)

# Filter for PPV patients
ppv_bloodflow <- octa_bloodflow_processed %>%
  filter(ID %in% ppv_patients)

ppv_thickness <- octa_thickness_processed %>%
  filter(ID %in% ppv_patients)

# -------------------- 4. Process vision data --------------------
# Create vision dataset
vision_data <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,
      surgery_eye_1 == 1 ~ os_corrected_bas,
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,
      TRUE ~ NA_real_
    ),
    post_vision_1m = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1m,
      surgery_eye_1 == 1 ~ os_corrected_1m,
      surgery_eye_1 == 2 ~ (od_corrected_1m + os_corrected_1m)/2,
      TRUE ~ NA_real_
    ),
    vision_improvement_1m = post_vision_1m - pre_vision
  ) %>%
  dplyr::select(ID, vision_improvement_1m, pre_vision, post_vision_1m, age, gender)

# Filter vision data for PPV patients
ppv_vision <- vision_data %>%
  filter(ID %in% ppv_patients)

# -------------------- 5. Filter OCTA parameters --------------------
# Function to filter blood flow layers (enhanced for 0_21 only)
filter_bloodflow_layers <- function(data) {
  layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
  # MODIFIED: Only include 0_21 region (removing 0_6)
  regions_of_interest <- c("0_21")
  
  # Create pattern for 0_21 region only
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("No target blood flow layer T0 parameters found!")
    return(list(data = data, params = character(0)))
  } else {
    cat("Found", length(params_T0), "target blood flow layer T0 parameters:\n")
    
    # Separate macular parameters for better reporting
    macular_params <- params_T0[grep("0_21_T0$", params_T0)]
    
    cat("- Macular region (0_21):", length(macular_params), "parameters\n")
    cat("- Total:", length(params_T0), "parameters\n\n")
    
    if(length(macular_params) > 0) {
      cat("Macular parameters:\n")
      cat(paste(gsub("_T0$", "", macular_params), collapse = "\n"), "\n\n")
    }
  }
  
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  if(length(params_T2) < length(params_T0)) {
    missing_params <- gsub("_T0$", "_T2", params_T0[!(gsub("_T0$", "_T2", params_T0) %in% params_T2)])
    warning("Missing T2 parameters: ", paste(missing_params, collapse = ", "))
  }
  
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  final_params_T0 <- paste0(valid_base_params, "_T0")
  final_params_T2 <- paste0(valid_base_params, "_T2")
  
  filtered_data <- data %>%
    dplyr::select(ID, all_of(c(final_params_T0, final_params_T2)))
  
  return(list(
    data = filtered_data,
    params_T0 = final_params_T0,
    params_T2 = final_params_T2,
    base_params = valid_base_params,
    macular_params = gsub("_T0$", "", macular_params),
    widefield_params = character(0)  # Empty since we removed 0_6
  ))
}

# Function to filter thickness layers (enhanced for 0_21 only)
filter_thickness_layers <- function(data) {
  layers_of_interest <- c("GCL.IPL", "INL", "Retina")
  # MODIFIED: Only include 0_21 region (removing 0_6)
  regions_of_interest <- c("0_21")
  
  # Create pattern for 0_21 region only
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("No target thickness layer T0 parameters found!")
    return(list(data = data, params = character(0)))
  } else {
    cat("Found", length(params_T0), "target thickness layer T0 parameters:\n")
    
    # Separate macular parameters for better reporting
    macular_params <- params_T0[grep("0_21_T0$", params_T0)]
    
    cat("- Macular region (0_21):", length(macular_params), "parameters\n")
    cat("- Total:", length(params_T0), "parameters\n\n")
    
    if(length(macular_params) > 0) {
      cat("Macular thickness parameters:\n")
      cat(paste(gsub("_T0$", "", macular_params), collapse = "\n"), "\n\n")
    }
  }
  
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  if(length(params_T2) < length(params_T0)) {
    missing_params <- gsub("_T0$", "_T2", params_T0[!(gsub("_T0$", "_T2", params_T0) %in% params_T2)])
    warning("Missing T2 parameters: ", paste(missing_params, collapse = ", "))
  }
  
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  final_params_T0 <- paste0(valid_base_params, "_T0")
  final_params_T2 <- paste0(valid_base_params, "_T2")
  
  filtered_data <- data %>%
    dplyr::select(ID, all_of(c(final_params_T0, final_params_T2)))
  
  return(list(
    data = filtered_data,
    params_T0 = final_params_T0,
    params_T2 = final_params_T2,
    base_params = valid_base_params,
    macular_params = gsub("_T0$", "", macular_params),
    widefield_params = character(0)  # Empty since we removed 0_6
  ))
}

# Apply filters
ppv_bloodflow_filtered <- filter_bloodflow_layers(ppv_bloodflow)
ppv_thickness_filtered <- filter_thickness_layers(ppv_thickness)

# -------------------- 6. Calculate OCTA improvements --------------------
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

# Calculate improvements
ppv_bloodflow_improvement <- calculate_improvement(
  ppv_bloodflow_filtered$data, 
  ppv_bloodflow_filtered$params_T0, 
  ppv_bloodflow_filtered$params_T2
)

ppv_thickness_improvement <- calculate_improvement(
  ppv_thickness_filtered$data, 
  ppv_thickness_filtered$params_T0, 
  ppv_thickness_filtered$params_T2
)

# -------------------- 7. Combine OCTA and Vision data --------------------
# Merge blood flow and thickness improvements
ppv_octa_combined <- ppv_bloodflow_improvement %>%
  full_join(ppv_thickness_improvement, by = "ID", suffix = c("_bloodflow", "_thickness"))

# Merge with vision data
comprehensive_data <- ppv_vision %>%
  dplyr::select(ID, vision_improvement_1m, pre_vision, age) %>%
  full_join(ppv_octa_combined, by = "ID")

# Print data combination summary (enhanced for regions)
cat("\n===== Comprehensive Data Summary =====\n")
cat("Total patients in comprehensive dataset:", nrow(comprehensive_data), "\n")
cat("Vision parameters: vision_improvement_1m, pre_vision, age\n")

if(length(ppv_bloodflow_filtered$macular_params) > 0) {
  cat("OCTA blood flow - Macular region (0_21):", length(ppv_bloodflow_filtered$macular_params), "parameters\n")
}
cat("Total OCTA blood flow parameters:", ncol(ppv_bloodflow_improvement) - 1, "\n")

if(length(ppv_thickness_filtered$macular_params) > 0) {
  cat("OCTA thickness - Macular region (0_21):", length(ppv_thickness_filtered$macular_params), "parameters\n")
}
cat("Total OCTA thickness parameters:", ncol(ppv_thickness_improvement) - 1, "\n")
cat("Total parameters:", ncol(comprehensive_data) - 1, "\n")

# -------------------- 8. Prepare data for clustering --------------------
# Analyze missing values
original_ids <- comprehensive_data$ID
clustering_data <- comprehensive_data %>% dplyr::select(-ID)

# Detailed missing value analysis
na_count_by_row <- rowSums(is.na(clustering_data))
na_count_table <- table(na_count_by_row)
na_count_by_col <- colSums(is.na(clustering_data))

cat("\n===== Missing Values Analysis =====\n")
cat("Total rows (patients):", nrow(clustering_data), "\n")
cat("Total columns (parameters):", ncol(clustering_data), "\n")

cat("\nPatient missing value distribution:\n")
for(na_count in names(na_count_table)) {
  cat("Patients with", na_count, "missing values:", na_count_table[na_count], "\n")
}

# Create summary of missing values by parameter type
na_summary <- data.frame(
  Parameter = names(clustering_data),
  Missing_Count = na_count_by_col,
  Missing_Percent = round(na_count_by_col/nrow(clustering_data)*100, 1),
  Data_Type = case_when(
    grepl("vision_improvement|pre_vision|age", names(clustering_data)) ~ "Vision",
    grepl("SVP|ICP|DCP|Choroid", names(clustering_data)) ~ "BloodFlow",
    grepl("GCL|INL|Retina", names(clustering_data)) ~ "Thickness",
    TRUE ~ "Other"
  )
)

cat("\nTop 15 parameters with most missing values:\n")
na_sorted <- na_summary %>% arrange(desc(Missing_Percent))
print(head(na_sorted, 15))

# Filter for complete cases
complete_rows <- complete.cases(clustering_data)
complete_data <- clustering_data[complete_rows, ]
complete_ids <- original_ids[complete_rows]

cat("\n===== Complete Data Statistics =====\n")
cat("Original patient count:", nrow(clustering_data), "\n")
cat("Complete data patient count:", nrow(complete_data), 
    "(", round(nrow(complete_data)/nrow(clustering_data)*100, 1), "%)\n")
cat("Removed", sum(!complete_rows), "patients with missing values\n")

# Count parameters by type
vision_params <- names(complete_data)[grepl("vision_improvement|pre_vision|age", names(complete_data))]
bloodflow_params <- names(complete_data)[grepl("SVP|ICP|DCP|Choroid", names(complete_data))]
thickness_params <- names(complete_data)[grepl("GCL|INL|Retina", names(complete_data))]

# Count parameters by type and region
bloodflow_macular <- bloodflow_params[grepl("0_21", bloodflow_params)]
bloodflow_widefield <- character(0)  # Empty since we removed 0_6
thickness_macular <- thickness_params[grepl("0_21", thickness_params)]
thickness_widefield <- character(0)  # Empty since we removed 0_6

cat("Vision parameters:", length(vision_params), "\n")
cat("Blood flow parameters - Total:", length(bloodflow_params), "\n")
cat("  - Macular (0_21):", length(bloodflow_macular), "\n")
cat("Thickness parameters - Total:", length(thickness_params), "\n")
cat("  - Macular (0_21):", length(thickness_macular), "\n")

# Check data adequacy
if(nrow(complete_data) < 5) {
  stop("Not enough patients with complete data. Consider imputation strategies.")
} else if(nrow(complete_data) < 10) {
  cat("⚠️ Warning: Small sample size (", nrow(complete_data), " patients)\n")
}

# Standardize data
comprehensive_data_std <- as.data.frame(scale(complete_data))

# -------------------- 1. 使用所有参数进行FCM聚类 --------------------

perform_fcm_all_parameters <- function() {
  cat("===== 使用所有14个参数进行FCM聚类 =====\n")
  
  # 使用完整的标准化数据（所有14个参数）
  all_data <- comprehensive_data_std
  
  cat("完整数据信息:\n")
  cat("患者数量:", nrow(all_data), "\n")
  cat("参数数量:", ncol(all_data), "\n")
  cat("参数列表:\n")
  for(i in 1:ncol(all_data)) {
    param_type <- case_when(
      names(all_data)[i] %in% vision_params ~ "视力",
      names(all_data)[i] %in% bloodflow_params ~ "血流",
      names(all_data)[i] %in% thickness_params ~ "厚度",
      TRUE ~ "其他"
    )
    cat(sprintf("%2d. %s (%s)\n", i, names(all_data)[i], param_type))
  }
  
  # -------------------- 2. FCM参数优化（全参数版本） --------------------
  
  cat("\n===== FCM参数优化（全部14个参数） =====\n")
  
  results_summary <- data.frame(
    k = integer(),
    fuzzy_coeff = numeric(),
    cluster_balance = numeric(),
    avg_membership = numeric(),
    partition_coefficient = numeric(),
    silhouette_score = numeric(),
    composite_score = numeric()
  )
  
  # 测试参数组合
  k_range <- 2:3  # 限制聚类数，避免过度拟合
  fuzzy_coeffs <- c(1.5, 2.0, 2.5)
  
  best_config <- list(score = -Inf)
  
  for(k in k_range) {
    for(m in fuzzy_coeffs) {
      cat(sprintf("测试 k=%d, m=%.1f (全部14个参数)...\n", k, m))
      
      # 执行FCM聚类
      set.seed(123)
      fcm_result <- cmeans(all_data, centers = k, m = m, iter.max = 100)
      
      # 计算评估指标
      hard_clusters <- fcm_result$cluster
      membership_matrix <- fcm_result$membership
      max_memberships <- apply(membership_matrix, 1, max)
      
      # 1. 聚类平衡度
      cluster_sizes <- table(hard_clusters)
      balance_score <- min(cluster_sizes) / max(cluster_sizes)
      
      # 2. 平均隶属度
      avg_membership <- mean(max_memberships)
      
      # 3. 分割系数
      pc <- sum(membership_matrix^2) / nrow(membership_matrix)
      
      # 4. 轮廓系数
      sil_score <- NA
      if(length(unique(hard_clusters)) > 1) {
        tryCatch({
          sil <- silhouette(hard_clusters, dist(all_data))
          sil_score <- mean(sil[, 3])
        }, error = function(e) {
          sil_score <<- balance_score * avg_membership
        })
      }
      
      if(is.na(sil_score)) {
        sil_score <- balance_score * avg_membership
      }
      
      # 5. 综合评分
      composite_score <- (balance_score * 0.4 + 
                            avg_membership * 0.3 + 
                            pc * 0.2 + 
                            sil_score * 0.1)
      
      # 记录结果
      results_summary <- rbind(results_summary, data.frame(
        k = k,
        fuzzy_coeff = m,
        cluster_balance = balance_score,
        avg_membership = avg_membership,
        partition_coefficient = pc,
        silhouette_score = sil_score,
        composite_score = composite_score
      ))
      
      # 更新最佳配置
      if(composite_score > best_config$score) {
        best_config <- list(
          k = k,
          m = m,
          score = composite_score,
          fcm_result = fcm_result,
          metrics = list(
            balance = balance_score,
            membership = avg_membership,
            pc = pc,
            silhouette = sil_score
          )
        )
      }
      
      cat(sprintf("  平衡度=%.3f, 隶属度=%.3f, PC=%.3f, 综合=%.3f\n",
                  balance_score, avg_membership, pc, composite_score))
      cat("  聚类分布:", paste(cluster_sizes, collapse=" vs "), "\n\n")
    }
  }
  
  cat("=== 全参数FCM最优结果 ===\n")
  cat(sprintf("最优聚类数: k=%d\n", best_config$k))
  cat(sprintf("最优模糊系数: m=%.1f\n", best_config$m))
  cat(sprintf("聚类平衡度: %.3f\n", best_config$metrics$balance))
  cat(sprintf("平均隶属度: %.3f\n", best_config$metrics$membership))
  cat(sprintf("分割系数: %.3f\n", best_config$metrics$pc))
  cat(sprintf("综合得分: %.3f\n", best_config$score))
  
  return(list(
    best_config = best_config,
    all_results = results_summary,
    all_data = all_data
  ))
}

# 执行全参数FCM分析
fcm_all_results <- perform_fcm_all_parameters()

# -------------------- 3. 分析全参数FCM结果 --------------------

analyze_fcm_all_results <- function(fcm_result, all_data) {
  cat("\n===== 全参数FCM聚类结果分析 =====\n")
  
  membership_matrix <- fcm_result$membership
  hard_clusters <- fcm_result$cluster
  
  # 创建结果数据框
  results_df <- data.frame(
    Patient_ID = complete_ids,
    Hard_Cluster = hard_clusters,
    Max_Membership = apply(membership_matrix, 1, max)
  )
  
  # 添加每个聚类的隶属度
  for(i in 1:ncol(membership_matrix)) {
    results_df[[paste0("Membership_C", i)]] <- membership_matrix[, i]
  }
  
  cat("患者归属详细分析（全参数）:\n")
  for(i in 1:nrow(results_df)) {
    patient_id <- results_df$Patient_ID[i]
    max_mem <- results_df$Max_Membership[i]
    main_cluster <- results_df$Hard_Cluster[i]
    
    # 获取所有隶属度
    memberships <- membership_matrix[i, ]
    
    cat(sprintf("患者 %s: 主要归属C%d (%.3f)", patient_id, main_cluster, max_mem))
    
    # 显示所有隶属度
    cat(" [")
    for(j in 1:length(memberships)) {
      cat(sprintf("C%d:%.3f", j, memberships[j]))
      if(j < length(memberships)) cat(", ")
    }
    cat("]")
    
    # 分类患者类型
    if(max_mem >= 0.8) {
      patient_type <- "强归属"
    } else if(max_mem >= 0.6) {
      patient_type <- "中等归属"
    } else {
      patient_type <- "模糊归属"
    }
    
    cat(sprintf(" (%s)\n", patient_type))
  }
  
  # 统计归属类型
  strong_patients <- sum(results_df$Max_Membership >= 0.8)
  medium_patients <- sum(results_df$Max_Membership >= 0.6 & results_df$Max_Membership < 0.8)
  fuzzy_patients <- sum(results_df$Max_Membership < 0.6)
  
  cat(sprintf("\n归属类型统计（全参数）:\n"))
  cat(sprintf("- 强归属患者 (≥0.8): %d名\n", strong_patients))
  cat(sprintf("- 中等归属患者 (0.6-0.8): %d名\n", medium_patients))
  cat(sprintf("- 模糊归属患者 (<0.6): %d名\n", fuzzy_patients))
  
  return(results_df)
}

fcm_all_results_df <- analyze_fcm_all_results(fcm_all_results$best_config$fcm_result, 
                                              fcm_all_results$all_data)

# -------------------- 4. 创建显示所有14个参数的完整热图 --------------------

# ===============================================
# 简单修复：只移除热图底部的NA行
# 保留原来的完整参数热图，只去掉NA分面
# ===============================================

create_fcm_heatmap_remove_na_row <- function() {
  cat("===== 移除热图底部NA行（保留所有其他参数） =====\n")
  
  # -------------------- 1. 使用原来的热图数据 --------------------
  
  # 重新读取之前创建的完整热图数据
  if(exists("complete_heatmap_data")) {
    heatmap_data <- complete_heatmap_data
    cat("使用已存在的complete_heatmap_data\n")
  } else {
    # 如果不存在，重新创建
    cat("重新创建热图数据...\n")
    
    fcm_result <- fcm_all_results$best_config$fcm_result
    hard_clusters <- fcm_result$cluster
    
    # 使用原始数据
    original_all_data <- comprehensive_data[1:length(hard_clusters), ]
    
    # 合并数据与聚类结果
    data_with_clusters <- original_all_data %>%
      mutate(FCM_Cluster = hard_clusters)
    
    # 获取所有参数名（排除ID列）
    all_param_names <- names(original_all_data)[names(original_all_data) != "ID"]
    
    # 计算聚类均值
    cluster_means <- data_with_clusters %>%
      group_by(FCM_Cluster) %>%
      summarise(across(all_of(all_param_names), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
    
    # 参数分类
    vision_params_all <- all_param_names[all_param_names %in% vision_params]
    bloodflow_params_all <- all_param_names[all_param_names %in% bloodflow_params]
    thickness_params_all <- all_param_names[all_param_names %in% thickness_params]
    
    # 转换为长格式
    heatmap_data <- cluster_means %>%
      pivot_longer(cols = all_of(all_param_names), 
                   names_to = "Parameter", 
                   values_to = "Mean_Value") %>%
      mutate(
        Parameter_Type = case_when(
          Parameter %in% vision_params_all ~ "Vision",
          Parameter %in% bloodflow_params_all ~ "Blood Flow - WF",
          Parameter %in% thickness_params_all ~ "Thickness - WF",
          TRUE ~ "Other"
        ),
        Parameter_Clean = gsub("_improvement", "", Parameter),
        Parameter_Clean = gsub("_0_21", "", Parameter_Clean),
        Parameter_Clean = gsub("_1m", "", Parameter_Clean),
        Parameter_Clean = gsub("Thickness_", "", Parameter_Clean),
        Parameter_Clean = gsub("_", " ", Parameter_Clean),
        Cluster_Label = paste("Cluster", FCM_Cluster)
      )
  }
  
  cat("原始热图数据统计:\n")
  cat("总数据点:", nrow(heatmap_data), "\n")
  cat("参数类型分布:\n")
  print(table(heatmap_data$Parameter_Type))
  
  # -------------------- 2. 移除NA和Other类型 --------------------
  
  cat("\n移除NA值和Other类型...\n")
  
  # 移除NA值和"Other"类型、"NA"类型的数据
  clean_heatmap_data <- heatmap_data %>%
    filter(!is.na(Mean_Value)) %>%                    # 移除NA值
    filter(Parameter_Type != "Other") %>%             # 移除Other类型
    filter(Parameter_Type != "NA") %>%                # 移除NA类型（如果存在）
    filter(!is.na(Parameter_Type)) %>%                # 移除参数类型为NA的
    filter(Parameter_Type %in% c("Vision", "Blood Flow - WF", "Thickness - WF"))  # 只保留这三种类型
  
  cat("清理后数据统计:\n")
  cat("清理后数据点:", nrow(clean_heatmap_data), "\n")
  cat("保留的参数类型:\n")
  print(table(clean_heatmap_data$Parameter_Type))
  
  # -------------------- 3. 重新设置因子顺序 --------------------
  
  # 设置参数类型顺序
  clean_heatmap_data$Parameter_Type <- factor(clean_heatmap_data$Parameter_Type, 
                                              levels = c("Blood Flow - WF", "Thickness - WF", "Vision"))
  
  # 按类型内部排序参数
  clean_heatmap_data <- clean_heatmap_data %>%
    arrange(Parameter_Type, Parameter_Clean)
  
  # 设置参数顺序
  clean_heatmap_data$Parameter_Clean <- factor(clean_heatmap_data$Parameter_Clean, 
                                               levels = unique(clean_heatmap_data$Parameter_Clean))
  
  # -------------------- 4. 创建清洁的热图（无NA行） --------------------
  
  cat("\n创建无NA行的完整热图...\n")
  
  p_final_heatmap <- ggplot(clean_heatmap_data, aes(x = Parameter_Clean, y = Cluster_Label, fill = Mean_Value)) +
    geom_tile(color = "white", size = 0.5) +
    # 添加数值标签
    geom_text(aes(label = sprintf("%.2f", Mean_Value)), 
              color = "black", size = 3.2, fontface = "bold") +
    # 颜色配置
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 0,
      name = "Mean\nValue",
      guide = guide_colorbar(barwidth = 1, barheight = 10)
    ) +
    # 分面设置
    facet_wrap(~ Parameter_Type, scales = "free_x", ncol = 1, strip.position = "top") +
    # 主题设置
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      strip.text = element_text(size = 12, face = "bold", color = "darkblue"),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.4, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 1, size = 10, color = "gray50"),
      legend.position = "right"
    ) +
    # 标签
    labs(
      title = "Complete FCM Cluster Comparison Heatmap",
      subtitle = "All Valid Parameters: Vision + OCTA WF Parameters | n = 15",
      x = "Parameters",
      y = "FCM Cluster",
      caption = "Red = Higher values, Blue = Lower values"
    )
  
  # 保存最终热图
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  ggsave("plots/fcm_final_heatmap_no_na_row.pdf", p_final_heatmap, 
         width = 18, height = 10, dpi = 300)
  ggsave("plots/fcm_final_heatmap_no_na_row.png", p_final_heatmap, 
         width = 18, height = 10, dpi = 300)
  
  cat("最终FCM热图已保存到:\n")
  cat("- plots/fcm_final_heatmap_no_na_row.pdf\n")
  cat("- plots/fcm_final_heatmap_no_na_row.png\n")
  
  # -------------------- 5. 显示最终统计 --------------------
  
  cat("\n=== 最终热图统计 ===\n")
  
  param_counts <- clean_heatmap_data %>%
    group_by(Parameter_Type) %>%
    summarise(Count = n_distinct(Parameter_Clean), .groups = 'drop')
  
  print(param_counts)
  
  total_params <- sum(param_counts$Count)
  cat("显示的总参数数:", total_params, "\n")
  cat("聚类分布: 7 vs 8\n")
  cat("成功移除了底部的NA行！\n")
  
  # 保存最终数据
  write.csv(clean_heatmap_data, "fcm_final_heatmap_data.csv", row.names = FALSE)
  cat("最终数据已保存到: fcm_final_heatmap_data.csv\n")
  
  return(clean_heatmap_data)
}

# 执行最终修复
final_heatmap_data <- create_fcm_heatmap_remove_na_row()



# 保存全参数FCM结果
write.csv(fcm_all_results_df, "fcm_all_parameters_results.csv", row.names = FALSE)
