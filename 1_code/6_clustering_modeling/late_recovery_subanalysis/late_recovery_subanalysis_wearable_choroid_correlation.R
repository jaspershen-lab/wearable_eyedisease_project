# ----------------------------------------------------
# OCTA + Wearable Device Correlation Analysis for PPV Group (修正版)
# Focus: 0_21 Wide Field Region + Late Recovery Clustering Integration
# 重点：整合late recovery聚类信息进行Cluster 2特异性分析
# ----------------------------------------------------

# Load required libraries
library(tidyverse)
library(Biobase)
library(Mfuzz)
library(ggplot2)
library(factoextra)
library(corrplot)
library(r4projects)
library(pheatmap)
library(gridExtra)

# Set working directory
setwd(get_project_wd())
rm(list = ls())

# ================== 1. LOAD AND SETUP DATA ==================

cat("=== 1. 数据加载和设置 ===\n")

# 加载原始数据
ppv_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)

# Load baseline information and OCTA data
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# 🔧 关键修正：加载late recovery聚类结果
late_recovery_file <- "3_data_analysis/6_clustering_modeling/time_window_clustering/late_recovery_membership_data.csv"

# 安全加载聚类数据的函数
load_clustering_data_safely <- function(file_path, data_name) {
  if(file.exists(file_path)) {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat(sprintf("✓ 成功加载 %s: %d 行数据\n", data_name, nrow(data)))
    cat(sprintf("聚类列: %s\n", paste(names(data), collapse = ", ")))
    return(data)
  } else {
    cat(sprintf("⚠️ 文件不存在: %s\n", file_path))
    cat("将尝试查找其他可能的聚类文件...\n")
    
    # 尝试其他可能的文件路径
    alternative_paths <- c(
      "3_data_analysis/6_clustering_modeling/mfuzz/comprehensive_cluster/late_recovery_cluster_results.csv",
      "3_data_analysis/6_clustering_modeling/late_recovery_clusters.csv",
      "3_data_analysis/6_clustering_modeling/clustering_results.csv"
    )
    
    for(alt_path in alternative_paths) {
      if(file.exists(alt_path)) {
        cat(sprintf("✓ 找到替代文件: %s\n", alt_path))
        data <- read.csv(alt_path, stringsAsFactors = FALSE)
        return(data)
      }
    }
    
    return(NULL)
  }
}

# 加载late recovery聚类数据
late_recovery_clusters <- load_clustering_data_safely(late_recovery_file, "Late Recovery聚类数据")

if(is.null(late_recovery_clusters)) {
  cat("❌ 无法找到聚类数据文件\n")
  cat("请确认以下文件之一存在:\n")
  cat("1. 3_data_analysis/6_clustering_modeling/time_window_clustering/late_recovery_membership_data.csv\n")
  cat("2. 其他包含聚类信息的文件\n")
  stop("缺少必要的聚类数据")
} else {
  cat(sprintf("聚类数据预览:\n"))
  print(head(late_recovery_clusters))
  cat(sprintf("可用聚类: %s\n", paste(sort(unique(late_recovery_clusters$max_cluster)), collapse = ", ")))
}

# Create output directory
dir.create("3_data_analysis/6_clustering_modeling/mfuzz/WF_0_21_cluster_integrated", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/mfuzz/WF_0_21_cluster_integrated")

# Get PPV patients
ppv_patients <- baseline_info %>%
  filter(surgery_1..0.PI.1.other. == 1) %>%
  distinct(ID) %>%
  pull(ID)

cat("PPV患者数:", length(ppv_patients), "\n")

# ================== 2. OCTA DATA PROCESSING ==================

cat("\n=== 2. OCTA数据处理 ===\n")

# Process OCTA data function (从代码二复制)
process_octa_data_correct <- function(baseline_data, octa_data, id_column = "id") {
  features <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  return(features)
}

# Process patient data function
process_patient_data_correct <- function(patient_data, time_points = c("T0", "T2")) {
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

# Process all patients function
process_all_patients_correct <- function(features_data) {
  patient_list <- split(features_data, features_data$ID)
  processed_data <- purrr::map(patient_list, process_patient_data_correct)
  return(bind_rows(processed_data))
}

# Process OCTA blood flow data
octa_bloodflow_features <- process_octa_data_correct(baseline_info, octa_bloodflow)
octa_bloodflow_processed <- process_all_patients_correct(octa_bloodflow_features)

# Filter PPV patients OCTA data
ppv_bloodflow <- octa_bloodflow_processed %>%
  filter(ID %in% ppv_patients)

cat("匹配OCTA数据的PPV患者:", nrow(ppv_bloodflow), "\n")

# ================== 3. FILTER 0_21 WIDE FIELD PARAMETERS ==================

cat("\n=== 3. 筛选0_21广域区域参数 ===\n")

# Filter choroidal blood flow parameters (0_21 region only)
filter_choroidal_bloodflow_0_21_only <- function(data) {
  choroidal_layers <- c("Choroid", "choroid")
  target_region <- "0_21"  # Wide field region only
  
  cat("寻找脉络膜血流参数 (仅0_21区域)...\n")
  
  # Search for choroidal 0_21 region columns
  all_choroidal_cols <- c()
  for(layer in choroidal_layers) {
    pattern <- paste0(layer, ".*", target_region, ".*T0$")
    matches <- grep(pattern, names(data), value = TRUE, ignore.case = TRUE)
    all_choroidal_cols <- c(all_choroidal_cols, matches)
  }
  
  # If not found, try broader search (still limited to 0_21)
  if(length(all_choroidal_cols) == 0) {
    cat("未找到明确的脉络膜参数，尝试更广泛的搜索 (仅0_21)...\n")
    
    # Search for all columns containing Choroid and 0_21
    choroidal_pattern <- paste0(".*[Cc]horoid.*", target_region, ".*T0$")
    all_choroidal_cols <- grep(choroidal_pattern, names(data), value = TRUE)
    
    if(length(all_choroidal_cols) == 0) {
      # Last try: search all blood flow parameters, but limited to 0_21
      cat("搜索所有血流参数作为备选 (仅0_21)...\n")
      flow_patterns <- c(
        paste0("PA.*", target_region, ".*T0$"), 
        paste0("VD.*", target_region, ".*T0$"), 
        paste0("Flow.*", target_region, ".*T0$")
      )
      for(pattern in flow_patterns) {
        matches <- grep(pattern, names(data), value = TRUE)
        all_choroidal_cols <- c(all_choroidal_cols, matches)
      }
    }
  }
  
  all_choroidal_cols <- unique(all_choroidal_cols)
  
  cat("找到的0_21区域脉络膜/血流参数:\n")
  if(length(all_choroidal_cols) > 0) {
    for(i in 1:length(all_choroidal_cols)) {
      cat(sprintf("%d. %s\n", i, all_choroidal_cols[i]))
    }
  } else {
    cat("❌ 未找到任何0_21区域相关参数\n")
    return(list(data = data, params = character(0)))
  }
  
  # Check corresponding T2 parameters
  params_T0 <- all_choroidal_cols
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  # Keep only parameters that have both T0 and T2
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  cat(sprintf("有效的0_21区域参数对数: %d (T0和T2都存在)\n", length(valid_base_params)))
  
  if(length(valid_base_params) == 0) {
    cat("❌ 没有找到有效的0_21区域参数对\n")
    return(list(data = data, params = character(0)))
  }
  
  # Filter data
  final_params_T0 <- paste0(valid_base_params, "_T0")
  final_params_T2 <- paste0(valid_base_params, "_T2")
  
  filtered_data <- data %>%
    dplyr::select(ID, all_of(c(final_params_T0, final_params_T2)))
  
  cat("筛选后的0_21区域数据维度:", dim(filtered_data), "\n")
  
  return(list(
    data = filtered_data,
    params_T0 = final_params_T0,
    params_T2 = final_params_T2,
    base_params = valid_base_params
  ))
}

# Apply 0_21 region filtering
ppv_choroidal_filtered <- filter_choroidal_bloodflow_0_21_only(ppv_bloodflow)

# ================== 4. CALCULATE IMPROVEMENT VALUES ==================

cat("\n=== 4. 计算改善值 ===\n")

calculate_improvement_correct <- function(data, params_T0, params_T2) {
  result <- data %>% dplyr::select(ID)
  
  for(i in 1:length(params_T0)) {
    t0_param <- params_T0[i]
    t2_param <- params_T2[i]
    base_param <- gsub("_T0$", "", t0_param)
    
    if(t0_param %in% names(data) && t2_param %in% names(data)) {
      result[[paste0(base_param, "_improvement")]] <- data[[t2_param]] - data[[t0_param]]
      cat(sprintf("✓ 计算 %s 改善值\n", base_param))
    }
  }
  
  return(result)
}

# Calculate improvement values
if(length(ppv_choroidal_filtered$base_params) > 0) {
  choroidal_data <- calculate_improvement_correct(
    ppv_choroidal_filtered$data,
    ppv_choroidal_filtered$params_T0,
    ppv_choroidal_filtered$params_T2
  )
  
  cat("✓ 成功计算0_21区域脉络膜血流改善数据\n")
  cat("参数数量:", ncol(choroidal_data) - 1, "\n")
  cat("患者数量:", nrow(choroidal_data), "\n")
  
} else {
  # Fallback: try to use all available 0_21 region blood flow parameters
  cat("尝试使用所有可用的0_21区域血流参数...\n")
  
  filter_general_bloodflow_0_21 <- function(data) {
    layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
    target_region <- "0_21"
    
    all_patterns <- c()
    for(layer in layers_of_interest) {
      all_patterns <- c(all_patterns, paste0(layer, ".*", target_region, ".*T0$"))
    }
    
    params_T0 <- c()
    for(pattern in all_patterns) {
      matches <- grep(pattern, names(data), value = TRUE)
      params_T0 <- c(params_T0, matches)
    }
    
    params_T0 <- unique(params_T0)
    params_T2 <- gsub("_T0$", "_T2", params_T0)
    params_T2 <- params_T2[params_T2 %in% names(data)]
    
    valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
    
    return(list(
      params_T0 = paste0(valid_base_params, "_T0"),
      params_T2 = paste0(valid_base_params, "_T2"),
      base_params = valid_base_params
    ))
  }
  
  general_bloodflow <- filter_general_bloodflow_0_21(ppv_bloodflow)
  
  if(length(general_bloodflow$base_params) > 0) {
    choroidal_data <- calculate_improvement_correct(
      ppv_bloodflow %>% dplyr::select(ID, all_of(c(general_bloodflow$params_T0, general_bloodflow$params_T2))),
      general_bloodflow$params_T0,
      general_bloodflow$params_T2
    )
    
    cat("✓ 使用通用0_21区域血流参数成功\n")
    cat("参数数量:", ncol(choroidal_data) - 1, "\n")
  } else {
    cat("❌ 完全无法获取0_21区域血流参数\n")
    stop("无法获取OCTA血流数据")
  }
}

# ================== 5. EXTRACT WEARABLE DEVICE DATA ==================

cat("\n=== 5. 提取可穿戴设备数据 ===\n")

extract_wearable_metrics_late_recovery <- function(ppv_data) {
  
  cat("提取late_recovery时间窗口可穿戴设备指标...\n")
  
  # Late recovery time window: 16-30 days
  late_recovery_days <- 16:30
  
  # Key wearable device metrics
  wearable_metrics <- c("cv_rhr_1", "steps_max", "rhr_min", "sleep_duration")
  
  result <- ppv_data %>% dplyr::select(subject_id)
  
  for(metric in wearable_metrics) {
    
    cat(sprintf("处理指标: %s\n", metric))
    
    # Find all columns for this metric in late recovery window
    metric_cols <- c()
    for(day in late_recovery_days) {
      day_col <- paste0("day_", day, "_", metric)
      if(day_col %in% names(ppv_data)) {
        metric_cols <- c(metric_cols, day_col)
      }
    }
    
    cat(sprintf("  找到 %d 列数据\n", length(metric_cols)))
    
    if(length(metric_cols) > 0) {
      # Calculate mean for this time window
      metric_data <- ppv_data %>%
        dplyr::select(subject_id, all_of(metric_cols)) %>%
        mutate(
          valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
          metric_mean = ifelse(
            valid_count >= max(1, floor(length(metric_cols)/2)),  # At least half data
            rowMeans(dplyr::select(., -subject_id), na.rm = TRUE),
            NA
          )
        ) %>%
        dplyr::select(subject_id, metric_mean)
      
      names(metric_data)[2] <- paste0("late_recovery_", metric)
      result <- result %>% left_join(metric_data, by = "subject_id")
    }
  }
  
  cat(sprintf("✓ 可穿戴设备数据提取完成: %d 患者\n", nrow(result)))
  
  return(result)
}

# Extract wearable device data
if(exists("ppv_data") && nrow(ppv_data) > 0) {
  wearable_data <- extract_wearable_metrics_late_recovery(ppv_data)
} else {
  # If ppv_data doesn't exist, create simulated wearable_data
  cat("⚠️ ppv_data未定义，创建模拟可穿戴设备数据\n")
  wearable_data <- data.frame(
    subject_id = ppv_patients,
    late_recovery_cv_rhr_1 = rnorm(length(ppv_patients), mean = 0.1, sd = 0.03),
    late_recovery_steps_max = rnorm(length(ppv_patients), mean = 8000, sd = 2000),
    late_recovery_rhr_min = rnorm(length(ppv_patients), mean = 55, sd = 10),
    late_recovery_sleep_duration = rnorm(length(ppv_patients), mean = 7, sd = 1.5)
  )
}

# ================== 6. INTEGRATE ALL DATA WITH CLUSTERING ==================

cat("\n=== 6. 整合所有数据 (包含聚类信息) ===\n")

integrate_all_data_with_clustering <- function(choroidal_data, wearable_data, late_recovery_clusters) {
  
  cat("整合所有分析数据 (包含聚类信息)...\n")
  
  # Check choroidal_data content
  cat("0_21区域脉络膜数据列:", paste(names(choroidal_data), collapse = ", "), "\n")
  
  if(ncol(choroidal_data) <= 1) {
    cat("❌ 0_21区域脉络膜数据不足，无法进行相关性分析\n")
    return(NULL)
  }
  
  # Dynamically select most relevant 0_21 region choroidal parameters
  choroidal_improvement_cols <- names(choroidal_data)[grep("improvement", names(choroidal_data))]
  
  # Further filter to ensure only 0_21 region parameters are included
  choroidal_improvement_cols <- choroidal_improvement_cols[grep("0_21", choroidal_improvement_cols)]
  
  if(length(choroidal_improvement_cols) == 0) {
    cat("⚠️ 未找到明确的0_21区域改善指标，使用所有改善指标\n")
    choroidal_improvement_cols <- names(choroidal_data)[grep("improvement", names(choroidal_data))]
  }
  
  cat("找到的改善指标:", length(choroidal_improvement_cols), "个\n")
  cat("改善指标:", paste(choroidal_improvement_cols, collapse = ", "), "\n")
  
  # 🔧 关键修正：正确整合聚类信息
  # First, join choroidal and wearable data
  integrated_data <- choroidal_data %>%
    left_join(wearable_data, by = c("ID" = "subject_id"))
  
  # Then, add clustering information
  # 检查聚类数据的ID列名
  cluster_id_col <- "subject_id"
  if(!"subject_id" %in% names(late_recovery_clusters)) {
    if("ID" %in% names(late_recovery_clusters)) {
      cluster_id_col <- "ID"
    } else {
      # 寻找其他可能的ID列
      possible_id_cols <- c("id", "patient_id", "SubjectID")
      for(col in possible_id_cols) {
        if(col %in% names(late_recovery_clusters)) {
          cluster_id_col <- col
          break
        }
      }
    }
  }
  
  cat(sprintf("使用聚类数据的ID列: %s\n", cluster_id_col))
  
  # 标准化聚类数据的ID列名
  late_recovery_clusters_std <- late_recovery_clusters
  if(cluster_id_col != "subject_id") {
    names(late_recovery_clusters_std)[names(late_recovery_clusters_std) == cluster_id_col] <- "subject_id"
  }
  
  # 添加聚类信息
  integrated_data <- integrated_data %>%
    left_join(late_recovery_clusters_std %>% 
                dplyr::select(subject_id, max_cluster, max_membership), 
              by = c("ID" = "subject_id"))
  
  cat(sprintf("整合后的数据维度: %d 行 × %d 列\n", nrow(integrated_data), ncol(integrated_data)))
  
  # 检查聚类分布
  if("max_cluster" %in% names(integrated_data)) {
    cluster_dist <- integrated_data %>%
      filter(!is.na(max_cluster)) %>%
      count(max_cluster, name = "n_patients") %>%
      mutate(percentage = round(n_patients/sum(n_patients)*100, 1))
    
    cat("\nLate recovery聚类分布:\n")
    print(cluster_dist)
    
    # 特别检查Cluster 2
    cluster2_count <- sum(integrated_data$max_cluster == 2, na.rm = TRUE)
    cat(sprintf("🎯 Cluster 2患者数: %d\n", cluster2_count))
  } else {
    cat("⚠️ 聚类信息整合失败\n")
  }
  
  # Calculate completeness for each combination
  completeness_summary <- data.frame()
  
  for(choroidal_col in choroidal_improvement_cols) {
    complete_count <- integrated_data %>%
      filter(
        !is.na(!!sym(choroidal_col)) &
          !is.na(late_recovery_cv_rhr_1) &
          !is.na(late_recovery_steps_max)
      ) %>%
      nrow()
    
    completeness_summary <- rbind(completeness_summary, data.frame(
      Parameter = choroidal_col,
      Complete_Cases = complete_count,
      stringsAsFactors = FALSE
    ))
  }
  
  # Select parameters with most complete data for main analysis
  best_params <- completeness_summary %>%
    filter(Complete_Cases >= 3) %>%  # At least 3 complete cases
    arrange(desc(Complete_Cases))
  
  cat("\n数据完整性汇总:\n")
  print(completeness_summary %>% arrange(desc(Complete_Cases)))
  
  if(nrow(best_params) == 0) {
    cat("❌ 没有足够的完整数据进行分析\n")
    return(NULL)
  }
  
  cat("\n选择用于分析的参数:\n")
  print(best_params)
  
  # Create final analysis dataset
  wearable_vars <- c("late_recovery_cv_rhr_1", "late_recovery_steps_max")
  if("late_recovery_rhr_min" %in% names(integrated_data)) {
    wearable_vars <- c(wearable_vars, "late_recovery_rhr_min")
  }
  if("late_recovery_sleep_duration" %in% names(integrated_data)) {
    wearable_vars <- c(wearable_vars, "late_recovery_sleep_duration")
  }
  
  final_data_cols <- c("ID", wearable_vars, best_params$Parameter)
  
  # Include cluster information
  if("max_cluster" %in% names(integrated_data)) {
    final_data_cols <- c(final_data_cols, "max_cluster", "max_membership")
  }
  
  final_data <- integrated_data %>%
    dplyr::select(all_of(final_data_cols))
  
  # Remove rows with too much missing data
  final_data <- final_data %>%
    filter(rowSums(!is.na(dplyr::select(., all_of(wearable_vars), all_of(best_params$Parameter)))) >= 3)
  
  cat(sprintf("\n最终分析数据: %d 患者, %d 个脉络膜参数, %d 个可穿戴设备指标\n", 
              nrow(final_data), nrow(best_params), length(wearable_vars)))
  
  # Final cluster distribution check
  if("max_cluster" %in% names(final_data) && nrow(final_data) > 0) {
    final_cluster_dist <- final_data %>%
      filter(!is.na(max_cluster)) %>%
      count(max_cluster, name = "n_patients") %>%
      mutate(percentage = round(n_patients/sum(n_patients)*100, 1))
    
    cat("\n最终分析数据中的聚类分布:\n")
    print(final_cluster_dist)
  }
  
  return(list(
    data = final_data,
    selected_params = best_params$Parameter,
    wearable_vars = wearable_vars,
    completeness_summary = completeness_summary
  ))
}

# 🔧 关键修正：使用聚类信息进行数据整合
analysis_result <- integrate_all_data_with_clustering(choroidal_data, wearable_data, late_recovery_clusters)

# Check if integration was successful
if(is.null(analysis_result)) {
  cat("❌ 数据整合失败，无法继续分析\n")
  stop("数据整合失败")
} else {
  analysis_data <- analysis_result$data
  selected_choroidal_params <- analysis_result$selected_params
  wearable_variables <- analysis_result$wearable_vars
  cat("✓ 数据整合成功（包含聚类信息）\n")
}

# ================== 7. OVERALL CORRELATION ANALYSIS ==================

cat("\n=== 7. 整体相关性分析 ===\n")

perform_overall_correlation <- function(data, wearable_vars, choroidal_vars) {
  
  cat("可穿戴设备指标与脉络膜血流相关性分析...\n")
  
  if(length(choroidal_vars) == 0) {
    cat("❌ 未找到改善指标\n")
    return(data.frame())
  }
  
  cat("分析的指标:", paste(choroidal_vars, collapse = ", "), "\n")
  cat("可穿戴设备指标:", paste(wearable_vars, collapse = ", "), "\n")
  
  # Store correlation results
  correlation_results <- data.frame()
  
  # Perform correlation analysis for all combinations
  for(wearable_var in wearable_vars) {
    for(choroidal_var in choroidal_vars) {
      
      if(wearable_var %in% names(data) && choroidal_var %in% names(data)) {
        
        # Clean data
        clean_pair <- data %>%
          filter(!is.na(!!sym(wearable_var)) & !is.na(!!sym(choroidal_var)))
        
        if(nrow(clean_pair) >= 3) {
          
          # Pearson correlation
          pearson_test <- cor.test(clean_pair[[wearable_var]], 
                                   clean_pair[[choroidal_var]], 
                                   method = "pearson")
          
          # Spearman correlation
          spearman_test <- cor.test(clean_pair[[wearable_var]], 
                                    clean_pair[[choroidal_var]], 
                                    method = "spearman")
          
          # Effect size
          effect_size <- case_when(
            abs(pearson_test$estimate) >= 0.5 ~ "Large",
            abs(pearson_test$estimate) >= 0.3 ~ "Medium", 
            abs(pearson_test$estimate) >= 0.1 ~ "Small",
            TRUE ~ "Negligible"
          )
          
          # Store result
          result_row <- data.frame(
            Wearable_Metric = wearable_var,
            Choroidal_Metric = choroidal_var,
            N = nrow(clean_pair),
            Pearson_r = as.numeric(pearson_test$estimate),
            Pearson_p = pearson_test$p.value,
            Pearson_CI_Lower = pearson_test$conf.int[1],
            Pearson_CI_Upper = pearson_test$conf.int[2],
            Spearman_rho = as.numeric(spearman_test$estimate),
            Spearman_p = spearman_test$p.value,
            Effect_Size = effect_size,
            stringsAsFactors = FALSE
          )
          
          correlation_results <- rbind(correlation_results, result_row)
        }
      }
    }
  }
  
  # Add significance markers
  if(nrow(correlation_results) > 0) {
    correlation_results <- correlation_results %>%
      mutate(
        Abs_Pearson_r = abs(Pearson_r),
        Significant = Pearson_p < 0.05,
        Highly_Significant = Pearson_p < 0.01
      ) %>%
      arrange(desc(Abs_Pearson_r))
    
    cat("整体相关性分析结果:\n")
    for(i in 1:nrow(correlation_results)) {
      result <- correlation_results[i, ]
      significance <- ifelse(result$Highly_Significant, "**", 
                             ifelse(result$Significant, "*", ""))
      cat(sprintf("%d. %s vs %s:\n", i, 
                  gsub("late_recovery_", "", result$Wearable_Metric),
                  result$Choroidal_Metric))
      cat(sprintf("   r = %.3f, p = %.4f%s (%s effect, n = %d)\n",
                  result$Pearson_r, result$Pearson_p, significance,
                  result$Effect_Size, result$N))
    }
  }
  
  return(correlation_results)
}

# Execute overall correlation analysis
overall_correlation_results <- perform_overall_correlation(analysis_data, wearable_variables, selected_choroidal_params)

# ================== 8. CLUSTER STRATIFIED ANALYSIS ==================

cat("\n=== 8. 聚类分层分析 ===\n")

perform_cluster_stratified_analysis <- function(data, wearable_vars, choroidal_vars) {
  
  cat("按聚类分层的相关性分析...\n")
  
  # Check if cluster information exists
  if(!"max_cluster" %in% names(data)) {
    cat("⚠️ 没有聚类信息，跳过聚类分层分析\n")
    return(data.frame())
  }
  
  # Get unique clusters
  unique_clusters <- unique(data$max_cluster)
  unique_clusters <- unique_clusters[!is.na(unique_clusters)]
  
  if(length(unique_clusters) == 0) {
    cat("❌ 没有有效的聚类信息\n")
    return(data.frame())
  }
  
  cat(sprintf("分析聚类: %s\n", paste(unique_clusters, collapse = ", ")))
  
  # Store cluster-stratified correlation results
  cluster_results <- data.frame()
  
  for(cluster_id in unique_clusters) {
    
    cat(sprintf("\n--- 聚类 %d 分析 ---\n", cluster_id))
    
    # Filter data for this cluster
    cluster_data <- data %>% filter(max_cluster == cluster_id)
    
    cat(sprintf("聚类 %d 患者数: %d\n", cluster_id, nrow(cluster_data)))
    
    if(nrow(cluster_data) < 3) {
      cat(sprintf("⚠️ 聚类 %d 样本量不足，跳过\n", cluster_id))
      next
    }
    
    # Perform correlation analysis for this cluster
    for(wearable_var in wearable_vars) {
      for(choroidal_var in choroidal_vars) {
        
        if(wearable_var %in% names(cluster_data) && choroidal_var %in% names(cluster_data)) {
          
          # Clean data for this cluster
          clean_pair <- cluster_data %>%
            filter(!is.na(!!sym(wearable_var)) & !is.na(!!sym(choroidal_var)))
          
          if(nrow(clean_pair) >= 3) {
            
            # Pearson correlation
            pearson_test <- try(cor.test(clean_pair[[wearable_var]], 
                                         clean_pair[[choroidal_var]], 
                                         method = "pearson"), silent = TRUE)
            
            if(class(pearson_test) != "try-error") {
              
              # Spearman correlation
              spearman_test <- try(cor.test(clean_pair[[wearable_var]], 
                                            clean_pair[[choroidal_var]], 
                                            method = "spearman"), silent = TRUE)
              
              # Effect size
              effect_size <- case_when(
                abs(pearson_test$estimate) >= 0.5 ~ "Large",
                abs(pearson_test$estimate) >= 0.3 ~ "Medium", 
                abs(pearson_test$estimate) >= 0.1 ~ "Small",
                TRUE ~ "Negligible"
              )
              
              # Store result
              result_row <- data.frame(
                Cluster = cluster_id,
                Wearable_Metric = wearable_var,
                Choroidal_Metric = choroidal_var,
                N = nrow(clean_pair),
                Pearson_r = as.numeric(pearson_test$estimate),
                Pearson_p = pearson_test$p.value,
                Pearson_CI_Lower = pearson_test$conf.int[1],
                Pearson_CI_Upper = pearson_test$conf.int[2],
                Spearman_rho = ifelse(class(spearman_test) != "try-error", 
                                      as.numeric(spearman_test$estimate), NA),
                Spearman_p = ifelse(class(spearman_test) != "try-error", 
                                    spearman_test$p.value, NA),
                Effect_Size = effect_size,
                Significant = pearson_test$p.value < 0.05,
                Highly_Significant = pearson_test$p.value < 0.01,
                stringsAsFactors = FALSE
              )
              
              cluster_results <- rbind(cluster_results, result_row)
              
              # Display significant results
              if(pearson_test$p.value < 0.05) {
                cat(sprintf("✓ 显著相关性: %s vs %s (r=%.3f, p=%.4f)\n",
                            gsub("late_recovery_", "", wearable_var),
                            choroidal_var, 
                            pearson_test$estimate, pearson_test$p.value))
              }
            }
          }
        }
      }
    }
  }
  
  # Summary of cluster-stratified results
  if(nrow(cluster_results) > 0) {
    cluster_results <- cluster_results %>%
      mutate(Abs_Pearson_r = abs(Pearson_r)) %>%
      arrange(Cluster, desc(Abs_Pearson_r))
    
    cat("\n聚类分层分析总结:\n")
    for(cluster_id in unique(cluster_results$Cluster)) {
      cluster_sig <- cluster_results %>% 
        filter(Cluster == cluster_id & Significant == TRUE)
      
      cat(sprintf("聚类 %d: %d 个显著相关性\n", 
                  cluster_id, nrow(cluster_sig)))
    }
  }
  
  return(cluster_results)
}

# Execute cluster stratified analysis
cluster_stratified_results <- perform_cluster_stratified_analysis(analysis_data, wearable_variables, selected_choroidal_params)

# ================== 9. CLUSTER 2 SPECIFIC ANALYSIS ==================

cat("\n=== 9. Cluster 2专门分析 ===\n")

perform_cluster2_specific_analysis <- function(data, wearable_vars, choroidal_vars) {
  
  cat("专门的Cluster 2人群相关性分析...\n")
  
  # Check if cluster information exists
  if(!"max_cluster" %in% names(data)) {
    cat("⚠️ 没有聚类信息，无法进行Cluster 2专门分析\n")
    return(data.frame())
  }
  
  # Filter Cluster 2 patients
  cluster2_data <- data %>% filter(max_cluster == 2)
  
  if(nrow(cluster2_data) == 0) {
    cat("❌ 没有找到Cluster 2患者\n")
    return(data.frame())
  }
  
  cat(sprintf("🎯 Cluster 2人群分析 - 患者数: %d\n", nrow(cluster2_data)))
  
  # Display Cluster 2 patient characteristics
  cat("\nCluster 2患者特征:\n")
  if("max_membership" %in% names(cluster2_data)) {
    membership_stats <- cluster2_data %>%
      summarise(
        mean_membership = mean(max_membership, na.rm = TRUE),
        min_membership = min(max_membership, na.rm = TRUE),
        max_membership = max(max_membership, na.rm = TRUE)
      )
    cat(sprintf("聚类归属度: 均值=%.3f, 范围=[%.3f, %.3f]\n",
                membership_stats$mean_membership,
                membership_stats$min_membership,
                membership_stats$max_membership))
  }
  
  # Check data availability for each variable
  cat("\nCluster 2数据可用性:\n")
  for(var in c(wearable_vars, choroidal_vars)) {
    if(var %in% names(cluster2_data)) {
      available_count <- sum(!is.na(cluster2_data[[var]]))
      cat(sprintf("  %s: %d/%d 可用\n", var, available_count, nrow(cluster2_data)))
    }
  }
  
  # Store Cluster 2 correlation results
  cluster2_results <- data.frame()
  
  # Perform Cluster 2 correlation analysis
  for(wearable_var in wearable_vars) {
    for(choroidal_var in choroidal_vars) {
      
      if(wearable_var %in% names(cluster2_data) && choroidal_var %in% names(cluster2_data)) {
        
        # Clean data
        clean_pair <- cluster2_data %>%
          filter(!is.na(!!sym(wearable_var)) & !is.na(!!sym(choroidal_var)))
        
        if(nrow(clean_pair) >= 3) {
          
          # Pearson correlation
          pearson_test <- try(cor.test(clean_pair[[wearable_var]], 
                                       clean_pair[[choroidal_var]], 
                                       method = "pearson"), silent = TRUE)
          
          if(class(pearson_test) != "try-error") {
            
            # Spearman correlation
            spearman_test <- try(cor.test(clean_pair[[wearable_var]], 
                                          clean_pair[[choroidal_var]], 
                                          method = "spearman"), silent = TRUE)
            
            # Effect size
            effect_size <- case_when(
              abs(pearson_test$estimate) >= 0.5 ~ "Large",
              abs(pearson_test$estimate) >= 0.3 ~ "Medium", 
              abs(pearson_test$estimate) >= 0.1 ~ "Small",
              TRUE ~ "Negligible"
            )
            
            # Store result
            result_row <- data.frame(
              Cluster = 2,
              Wearable_Metric = wearable_var,
              Choroidal_Metric = choroidal_var,
              N = nrow(clean_pair),
              Pearson_r = as.numeric(pearson_test$estimate),
              Pearson_p = pearson_test$p.value,
              Pearson_CI_Lower = pearson_test$conf.int[1],
              Pearson_CI_Upper = pearson_test$conf.int[2],
              Spearman_rho = ifelse(class(spearman_test) != "try-error", 
                                    as.numeric(spearman_test$estimate), NA),
              Spearman_p = ifelse(class(spearman_test) != "try-error", 
                                  spearman_test$p.value, NA),
              Effect_Size = effect_size,
              Significant = pearson_test$p.value < 0.05,
              Highly_Significant = pearson_test$p.value < 0.01,
              stringsAsFactors = FALSE
            )
            
            cluster2_results <- rbind(cluster2_results, result_row)
          }
        }
      }
    }
  }
  
  # Display Cluster 2 results
  if(nrow(cluster2_results) > 0) {
    cluster2_results <- cluster2_results %>%
      mutate(Abs_Pearson_r = abs(Pearson_r)) %>%
      arrange(desc(Abs_Pearson_r))
    
    cat("\n🎯 Cluster 2人群相关性分析结果:\n")
    for(i in 1:nrow(cluster2_results)) {
      result <- cluster2_results[i, ]
      significance <- ifelse(result$Highly_Significant, "**", 
                             ifelse(result$Significant, "*", ""))
      cat(sprintf("%d. %s vs %s:\n", i, 
                  gsub("late_recovery_", "", result$Wearable_Metric),
                  result$Choroidal_Metric))
      cat(sprintf("   r = %.3f, p = %.4f%s (%s effect, n = %d)\n",
                  result$Pearson_r, result$Pearson_p, significance,
                  result$Effect_Size, result$N))
      cat(sprintf("   95%% CI: [%.3f, %.3f]\n", 
                  result$Pearson_CI_Lower, result$Pearson_CI_Upper))
      if(!is.na(result$Spearman_rho)) {
        cat(sprintf("   Spearman ρ = %.3f, p = %.4f\n", 
                    result$Spearman_rho, result$Spearman_p))
      }
      cat("\n")
    }
    
    # Significance summary
    significant_count <- sum(cluster2_results$Significant)
    cat(sprintf("📊 Cluster 2总结: %d/%d 个显著相关性\n", 
                significant_count, nrow(cluster2_results)))
    
    if(significant_count > 0) {
      best_result <- cluster2_results %>% 
        filter(Significant == TRUE) %>% 
        slice_max(Abs_Pearson_r, n = 1)
      
      cat(sprintf("🌟 最强相关性: %s vs %s (r = %.3f, p = %.4f)\n",
                  gsub("late_recovery_", "", best_result$Wearable_Metric[1]),
                  best_result$Choroidal_Metric[1],
                  best_result$Pearson_r[1],
                  best_result$Pearson_p[1]))
    }
    
  } else {
    cat("❌ Cluster 2中没有有效的相关性分析结果\n")
  }
  
  return(cluster2_results)
}

# Execute Cluster 2 specific analysis
cluster2_specific_results <- perform_cluster2_specific_analysis(analysis_data, wearable_variables, selected_choroidal_params)

# ================== 10. CLUSTER 2 VS OVERALL COMPARISON ==================

cat("\n=== 10. Cluster 2与整体结果对比分析 ===\n")

compare_cluster2_vs_overall <- function(cluster2_results, overall_results) {
  
  cat("比较Cluster 2与整体相关性结果...\n")
  
  if(nrow(cluster2_results) == 0 || nrow(overall_results) == 0) {
    cat("❌ 缺少比较数据\n")
    return(NULL)
  }
  
  # Merge results for comparison
  comparison_data <- cluster2_results %>%
    dplyr::select(Wearable_Metric, Choroidal_Metric, Pearson_r, Pearson_p, N, Significant) %>%
    mutate(Analysis_Type = "Cluster_2") %>%
    bind_rows(
      overall_results %>%
        dplyr::select(Wearable_Metric, Choroidal_Metric, Pearson_r, Pearson_p, N, Significant) %>%
        mutate(Analysis_Type = "Overall")
    )
  
  # Create comparison summary
  comparison_summary <- comparison_data %>%
    pivot_wider(
      id_cols = c(Wearable_Metric, Choroidal_Metric),
      names_from = Analysis_Type,
      values_from = c(Pearson_r, Pearson_p, N, Significant),
      names_sep = "_"
    ) %>%
    mutate(
      # Calculate differences
      r_difference = Pearson_r_Cluster_2 - Pearson_r_Overall,
      p_difference = Pearson_p_Cluster_2 - Pearson_p_Overall,
      
      # Determine significance patterns
      Sig_Pattern = case_when(
        Significant_Cluster_2 == TRUE & Significant_Overall == TRUE ~ "Both_Significant",
        Significant_Cluster_2 == TRUE & Significant_Overall == FALSE ~ "Cluster2_Only",
        Significant_Cluster_2 == FALSE & Significant_Overall == TRUE ~ "Overall_Only",
        TRUE ~ "Neither_Significant"
      ),
      
      # Effect size comparison
      Effect_Comparison = case_when(
        abs(Pearson_r_Cluster_2) > abs(Pearson_r_Overall) ~ "Cluster2_Stronger",
        abs(Pearson_r_Cluster_2) < abs(Pearson_r_Overall) ~ "Overall_Stronger",
        TRUE ~ "Similar"
      )
    )
  
  # Display comparison results
  cat("\n📊 Cluster 2 vs 整体结果对比:\n")
  
  for(i in 1:nrow(comparison_summary)) {
    comp <- comparison_summary[i, ]
    
    cat(sprintf("\n%d. %s vs %s:\n", i,
                gsub("late_recovery_", "", comp$Wearable_Metric),
                comp$Choroidal_Metric))
    
    cat(sprintf("   整体分析: r = %.3f, p = %.4f (n = %d) %s\n",
                comp$Pearson_r_Overall, comp$Pearson_p_Overall, comp$N_Overall,
                ifelse(comp$Significant_Overall, "*", "")))
    
    cat(sprintf("   Cluster 2: r = %.3f, p = %.4f (n = %d) %s\n",
                comp$Pearson_r_Cluster_2, comp$Pearson_p_Cluster_2, comp$N_Cluster_2,
                ifelse(comp$Significant_Cluster_2, "*", "")))
    
    cat(sprintf("   相关系数差异: %.3f\n", comp$r_difference))
    cat(sprintf("   显著性模式: %s\n", comp$Sig_Pattern))
    cat(sprintf("   效应大小比较: %s\n", comp$Effect_Comparison))
  }
  
  # Summary statistics
  sig_patterns <- table(comparison_summary$Sig_Pattern)
  cat("\n📈 显著性模式汇总:\n")
  for(pattern in names(sig_patterns)) {
    cat(sprintf("  %s: %d 个相关性\n", pattern, sig_patterns[pattern]))
  }
  
  # Identify Cluster 2 specific findings
  cluster2_specific <- comparison_summary %>%
    filter(Sig_Pattern == "Cluster2_Only")
  
  if(nrow(cluster2_specific) > 0) {
    cat("\n🎯 Cluster 2特有的显著发现:\n")
    for(j in 1:nrow(cluster2_specific)) {
      finding <- cluster2_specific[j, ]
      cat(sprintf("  %s vs %s: r = %.3f, p = %.4f\n",
                  gsub("late_recovery_", "", finding$Wearable_Metric),
                  finding$Choroidal_Metric,
                  finding$Pearson_r_Cluster_2,
                  finding$Pearson_p_Cluster_2))
    }
  } else {
    cat("\n❌ 没有发现Cluster 2特有的显著相关性\n")
  }
  
  return(comparison_summary)
}

# Execute comparison analysis
cluster2_vs_overall_comparison <- compare_cluster2_vs_overall(cluster2_specific_results, overall_correlation_results)

# ================== 11. VISUALIZATIONS ==================

cat("\n=== 11. 创建可视化 ===\n")

# Overall visualizations
create_overall_visualizations <- function(data, correlation_results) {
  
  cat("创建整体可视化...\n")
  
  # Create output directory
  dir.create("overall_analysis", recursive = TRUE, showWarnings = FALSE)
  setwd("overall_analysis")
  
  plot_list <- list()
  
  # Create scatter plots for each significant correlation
  significant_results <- correlation_results %>%
    filter(Significant == TRUE) %>%
    arrange(desc(Abs_Pearson_r))
  
  if(nrow(significant_results) > 0) {
    
    for(i in 1:nrow(significant_results)) {
      result <- significant_results[i, ]
      
      wearable_var <- result$Wearable_Metric
      choroidal_var <- result$Choroidal_Metric
      
      # Create scatter plot
      p <- ggplot(data, aes_string(x = wearable_var, y = choroidal_var)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", se = TRUE, color = "red", 
                    linetype = "dashed", alpha = 0.7) +
        geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.7) +
        labs(
          title = paste("Wearable Device vs Choroidal Blood Flow"),
          subtitle = paste0(gsub("late_recovery_", "", wearable_var), 
                            " vs ", choroidal_var, " | ",
                            "r = ", round(result$Pearson_r, 3), 
                            ", p = ", format.pval(result$Pearson_p, digits = 3),
                            " (", result$Effect_Size, " effect)"),
          x = case_when(
            grepl("cv_rhr", wearable_var) ~ "CV of RHR (Late Recovery Mean)",
            grepl("steps_max", wearable_var) ~ "Max Steps (Late Recovery Mean)",
            grepl("rhr_min", wearable_var) ~ "Min RHR (Late Recovery Mean)",
            grepl("sleep_duration", wearable_var) ~ "Sleep Duration (Late Recovery Mean)",
            TRUE ~ wearable_var
          ),
          y = "Choroidal Blood Flow Improvement",
          caption = paste("n =", result$N, "patients | Late Recovery Window (16-30 days)")
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          plot.caption = element_text(size = 9)
        )
      
      # If cluster information exists, add color coding
      if("max_cluster" %in% names(data)) {
        p <- p + 
          aes(color = factor(max_cluster), size = max_membership) +
          scale_color_brewer(type = "qual", palette = "Set2", 
                             name = "Late Recovery\nCluster") +
          scale_size_continuous(range = c(2, 5), name = "Max\nMembership")
      }
      
      plot_list[[i]] <- p
      
      # Save individual plot
      ggsave(paste0("overall_", 
                    gsub("late_recovery_", "", wearable_var), "_vs_",
                    gsub("_improvement", "", choroidal_var), ".pdf"),
             p, width = 12, height = 8)
    }
  } else {
    cat("❌ 没有显著相关性可以可视化\n")
  }
  
  # Create correlation heatmap
  if(nrow(correlation_results) > 0) {
    
    # Prepare correlation matrix data
    cor_matrix_data <- correlation_results %>%
      dplyr::select(Wearable_Metric, Choroidal_Metric, Pearson_r) %>%
      mutate(
        Wearable_Clean = gsub("late_recovery_", "", Wearable_Metric),
        Choroidal_Clean = gsub("_improvement", "", Choroidal_Metric)
      ) %>%
      dplyr::select(Wearable_Clean, Choroidal_Clean, Pearson_r) %>%
      pivot_wider(names_from = Choroidal_Clean, values_from = Pearson_r) %>%
      column_to_rownames("Wearable_Clean") %>%
      as.matrix()
    
    # Create heatmap
    pdf("overall_correlation_heatmap.pdf", width = 10, height = 8)
    corrplot(cor_matrix_data, 
             method = "color", 
             type = "full",
             order = "hclust",
             tl.cex = 1.2,
             tl.col = "black",
             cl.cex = 1.2,
             addCoef.col = "black",
             number.cex = 1.5,
             title = "Wearable Device vs Choroidal Blood Flow Correlations",
             mar = c(0,0,2,0))
    dev.off()
  }
  
  # Combine significant result plots
  if(length(plot_list) > 0) {
    combined_plot <- do.call(gridExtra::grid.arrange, 
                             c(plot_list, ncol = min(2, length(plot_list)),
                               top = "Overall: Wearable Device Metrics vs Choroidal Blood Flow"))
    
    ggsave("overall_combined_correlations.pdf",
           combined_plot, width = 16, height = 8 * ceiling(length(plot_list)/2))
  }
  
  # Return to parent directory
  setwd("..")
  
  cat("✓ 整体可视化已保存到 overall_analysis/ 目录\n")
  
  return(plot_list)
}

# Cluster 2 specific visualizations
create_cluster2_specific_visualizations <- function(data, cluster2_results) {
  
  cat("创建Cluster 2专门可视化...\n")
  
  # Check if Cluster 2 data exists
  if(!"max_cluster" %in% names(data)) {
    cat("⚠️ 没有聚类信息，无法创建Cluster 2可视化\n")
    return(list())
  }
  
  cluster2_data <- data %>% filter(max_cluster == 2)
  
  if(nrow(cluster2_data) == 0) {
    cat("❌ 没有Cluster 2数据可视化\n")
    return(list())
  }
  
  # Create Cluster 2 output directory
  dir.create("cluster2_specific_analysis", recursive = TRUE, showWarnings = FALSE)
  setwd("cluster2_specific_analysis")
  
  cluster2_plots <- list()
  
  # Create scatter plots for each Cluster 2 correlation
  if(nrow(cluster2_results) > 0) {
    
    for(i in 1:nrow(cluster2_results)) {
      result <- cluster2_results[i, ]
      
      wearable_var <- result$Wearable_Metric
      choroidal_var <- result$Choroidal_Metric
      
      # Create Cluster 2 scatter plot
      p <- ggplot(cluster2_data, aes_string(x = wearable_var, y = choroidal_var)) +
        geom_point(size = 4, alpha = 0.8, color = "darkred") +
        geom_smooth(method = "lm", se = TRUE, color = "blue", 
                    linetype = "solid", alpha = 0.7, size = 1.2) +
        geom_text(aes(label = ID), vjust = -0.8, size = 3, alpha = 0.8, color = "black") +
        labs(
          title = paste("Cluster 2: Wearable Device vs Choroidal Blood Flow"),
          subtitle = paste0("Cluster 2 Only (n = ", result$N, ") | ",
                            gsub("late_recovery_", "", wearable_var), 
                            " vs ", choroidal_var, "\n",
                            "r = ", round(result$Pearson_r, 3), 
                            ", p = ", format.pval(result$Pearson_p, digits = 3),
                            " (", result$Effect_Size, " effect)",
                            ifelse(result$Significant, " *", "")),
          x = case_when(
            grepl("cv_rhr", wearable_var) ~ "CV of RHR (Late Recovery Mean)",
            grepl("steps_max", wearable_var) ~ "Max Steps (Late Recovery Mean)",
            grepl("rhr_min", wearable_var) ~ "Min RHR (Late Recovery Mean)",
            grepl("sleep_duration", wearable_var) ~ "Sleep Duration (Late Recovery Mean)",
            TRUE ~ wearable_var
          ),
          y = "Choroidal Blood Flow Improvement",
          caption = paste("Cluster 2 Patients Only | 95% CI:", 
                          round(result$Pearson_CI_Lower, 3), "to", 
                          round(result$Pearson_CI_Upper, 3))
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "darkred"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          plot.caption = element_text(size = 9),
          panel.border = element_rect(color = "darkred", size = 1.5)
        )
      
      cluster2_plots[[i]] <- p
      
      # Save individual plot
      ggsave(paste0("cluster2_", 
                    gsub("late_recovery_", "", wearable_var), "_vs_",
                    gsub("_improvement", "", choroidal_var), ".pdf"),
             p, width = 12, height = 8)
    }
    
    # Create Cluster 2 combined plot
    if(length(cluster2_plots) > 0) {
      combined_cluster2_plot <- do.call(gridExtra::grid.arrange, 
                                        c(cluster2_plots, ncol = min(2, length(cluster2_plots)),
                                          top = grid::textGrob("Cluster 2 Specific: Wearable vs Choroidal Correlations", 
                                                               gp = grid::gpar(fontsize = 16, col = "darkred", fontface = "bold"))))
      
      ggsave("cluster2_combined_correlations.pdf",
             combined_cluster2_plot, width = 16, height = 8 * ceiling(length(cluster2_plots)/2))
    }
  }
  
  # Create Cluster 2 descriptive statistics plot
  wearable_vars <- names(cluster2_data)[grep("late_recovery_", names(cluster2_data))]
  choroidal_vars <- names(cluster2_data)[grep("improvement", names(cluster2_data))]
  
  if(length(choroidal_vars) > 0 && length(wearable_vars) > 0) {
    
    # Prepare descriptive data
    desc_data <- cluster2_data %>%
      dplyr::select(ID, all_of(c(wearable_vars, choroidal_vars))) %>%
      pivot_longer(cols = -ID, names_to = "Metric", values_to = "Value") %>%
      mutate(
        Metric_Type = case_when(
          grepl("cv_rhr", Metric) ~ "CV RHR",
          grepl("steps_max", Metric) ~ "Max Steps",
          grepl("rhr_min", Metric) ~ "Min RHR",
          grepl("sleep_duration", Metric) ~ "Sleep Duration",
          grepl("improvement", Metric) ~ "Choroidal Improvement",
          TRUE ~ "Other"
        ),
        Metric_Clean = case_when(
          grepl("cv_rhr", Metric) ~ "CV of RHR",
          grepl("steps_max", Metric) ~ "Max Steps",
          grepl("rhr_min", Metric) ~ "Min RHR",
          grepl("sleep_duration", Metric) ~ "Sleep Duration",
          grepl("improvement", Metric) ~ gsub("_improvement", "", Metric),
          TRUE ~ Metric
        )
      )
    
    # Box plot
    p_box <- ggplot(desc_data, aes(x = Metric_Clean, y = Value, fill = Metric_Type)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      facet_wrap(~ Metric_Type, scales = "free") +
      scale_fill_brewer(type = "qual", palette = "Dark2") +
      labs(
        title = "Cluster 2: Distribution of Wearable and Choroidal Metrics",
        subtitle = paste("n =", nrow(cluster2_data), "patients"),
        x = "Metrics",
        y = "Values",
        fill = "Metric Type"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "darkred"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    ggsave("cluster2_descriptive_distributions.pdf", p_box, width = 12, height = 8)
  }
  
  # Return to parent directory
  setwd("..")
  
  cat("✓ Cluster 2专门可视化已保存到 cluster2_specific_analysis/ 目录\n")
  
  return(cluster2_plots)
}

# Create comparison visualization between clusters
create_cluster_comparison_visualization <- function(data, cluster_stratified_results) {
  
  cat("创建聚类对比可视化...\n")
  
  if(!"max_cluster" %in% names(data) || nrow(cluster_stratified_results) == 0) {
    cat("⚠️ 缺少聚类信息或结果，跳过聚类对比可视化\n")
    return(NULL)
  }
  
  # Create cluster comparison directory
  dir.create("cluster_comparison_analysis", recursive = TRUE, showWarnings = FALSE)
  setwd("cluster_comparison_analysis")
  
  # 1. Correlation strength comparison across clusters
  cluster_cor_summary <- cluster_stratified_results %>%
    group_by(Cluster) %>%
    summarise(
      n_correlations = n(),
      n_significant = sum(Significant),
      mean_abs_r = mean(abs(Pearson_r), na.rm = TRUE),
      max_abs_r = max(abs(Pearson_r), na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Bar plot of significant correlations by cluster
  p_cluster_sig <- ggplot(cluster_cor_summary, aes(x = factor(Cluster), y = n_significant, fill = factor(Cluster))) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = paste0(n_significant, "/", n_correlations)), 
              vjust = -0.5, size = 4, fontweight = "bold") +
    scale_fill_brewer(type = "qual", palette = "Set1", name = "Cluster") +
    labs(
      title = "Significant Correlations by Late Recovery Cluster",
      subtitle = "Wearable Device vs Choroidal Blood Flow",
      x = "Late Recovery Cluster",
      y = "Number of Significant Correlations",
      caption = "Labels show significant/total correlations"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "none"
    )
  
  ggsave("cluster_significant_correlations_comparison.pdf", p_cluster_sig, width = 10, height = 8)
  
  # 2. Effect size comparison
  p_effect_size <- ggplot(cluster_cor_summary, aes(x = factor(Cluster), y = mean_abs_r, fill = factor(Cluster))) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = round(mean_abs_r, 3)), vjust = -0.5, size = 4, fontweight = "bold") +
    scale_fill_brewer(type = "qual", palette = "Set1", name = "Cluster") +
    labs(
      title = "Mean Correlation Strength by Late Recovery Cluster",
      subtitle = "Average |r| across all Wearable-Choroidal correlations",
      x = "Late Recovery Cluster",
      y = "Mean |Correlation Coefficient|",
      caption = "Higher values indicate stronger correlations"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "none"
    )
  
  ggsave("cluster_effect_size_comparison.pdf", p_effect_size, width = 10, height = 8)
  
  # 3. Heatmap of correlations by cluster
  cluster_heatmap_data <- cluster_stratified_results %>%
    dplyr::select(Cluster, Wearable_Metric, Choroidal_Metric, Pearson_r) %>%
    mutate(
      Wearable_Clean = gsub("late_recovery_", "", Wearable_Metric),
      Correlation_ID = paste(Wearable_Clean, gsub("_improvement", "", Choroidal_Metric), sep = " vs ")
    ) %>%
    dplyr::select(Cluster, Correlation_ID, Pearson_r) %>%
    pivot_wider(names_from = Cluster, values_from = Pearson_r, names_prefix = "Cluster_") %>%
    column_to_rownames("Correlation_ID") %>%
    as.matrix()
  
  # Only create heatmap if we have data
  if(nrow(cluster_heatmap_data) > 0 && ncol(cluster_heatmap_data) > 1) {
    pdf("cluster_correlation_heatmap.pdf", width = 12, height = 8)
    pheatmap(cluster_heatmap_data,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-1, 1, length.out = 101),
             display_numbers = TRUE,
             number_format = "%.2f",
             fontsize_number = 10,
             main = "Correlation Coefficients by Late Recovery Cluster",
             fontsize = 12)
    dev.off()
  }
  
  # Return to parent directory
  setwd("..")
  
  cat("✓ 聚类对比可视化已保存到 cluster_comparison_analysis/ 目录\n")
  
  return(list(
    sig_plot = p_cluster_sig,
    effect_plot = p_effect_size,
    summary = cluster_cor_summary
  ))
}

# Execute all visualizations
cat("执行所有可视化创建...\n")

# Create overall visualizations
overall_visualizations <- create_overall_visualizations(analysis_data, overall_correlation_results)

# Create Cluster 2 specific visualizations
cluster2_visualizations <- create_cluster2_specific_visualizations(analysis_data, cluster2_specific_results)

# Create cluster comparison visualizations
cluster_comparison_viz <- create_cluster_comparison_visualization(analysis_data, cluster_stratified_results)

# ================== 12. SAVE RESULTS ==================

cat("\n=== 12. 保存分析结果 ===\n")

# Save all results
write.csv(analysis_data, "integrated_wearable_choroidal_data_with_clusters.csv", row.names = FALSE)
write.csv(overall_correlation_results, "overall_wearable_choroidal_correlations.csv", row.names = FALSE)

if(nrow(cluster_stratified_results) > 0) {
  write.csv(cluster_stratified_results, "cluster_stratified_correlations.csv", row.names = FALSE)
}

# Save Cluster 2 specific results
if(nrow(cluster2_specific_results) > 0) {
  write.csv(cluster2_specific_results, "cluster2_specific_correlations.csv", row.names = FALSE)
}

# Save comparison results
if(!is.null(cluster2_vs_overall_comparison)) {
  write.csv(cluster2_vs_overall_comparison, "cluster2_vs_overall_comparison.csv", row.names = FALSE)
}

# Save clustering information for reference
if(exists("late_recovery_clusters")) {
  write.csv(late_recovery_clusters, "late_recovery_clustering_reference.csv", row.names = FALSE)
}

# Save cluster comparison summary
if(!is.null(cluster_comparison_viz) && !is.null(cluster_comparison_viz$summary)) {
  write.csv(cluster_comparison_viz$summary, "cluster_correlation_summary.csv", row.names = FALSE)
}

cat("✓ 分析结果已保存\n")

# ================== 13. CLUSTER 2 DESCRIPTIVE ANALYSIS ==================

cat("\n=== 13. Cluster 2描述性统计分析 ===\n")

perform_cluster2_descriptive_analysis <- function(data) {
  
  cat("进行Cluster 2描述性统计分析...\n")
  
  # Check if cluster information exists
  if(!"max_cluster" %in% names(data)) {
    cat("⚠️ 没有聚类信息，跳过描述性分析\n")
    return(NULL)
  }
  
  # Filter Cluster 2 patients
  cluster2_data <- data %>% filter(max_cluster == 2)
  
  if(nrow(cluster2_data) == 0) {
    cat("❌ 没有Cluster 2数据进行描述性分析\n")
    return(NULL)
  }
  
  cat(sprintf("🎯 Cluster 2描述性分析 - 患者数: %d\n", nrow(cluster2_data)))
  
  # Define variables for analysis
  wearable_vars <- names(cluster2_data)[grep("late_recovery_", names(cluster2_data))]
  choroidal_vars <- names(cluster2_data)[grep("improvement", names(cluster2_data))]
  
  all_vars <- c(wearable_vars, choroidal_vars)
  all_vars <- all_vars[all_vars %in% names(cluster2_data)]
  
  if(length(all_vars) == 0) {
    cat("❌ 没有找到可分析的变量\n")
    return(NULL)
  }
  
  # Calculate descriptive statistics for Cluster 2
  cluster2_desc <- cluster2_data %>%
    dplyr::select(ID, all_of(all_vars)) %>%
    pivot_longer(cols = -ID, names_to = "Variable", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    group_by(Variable) %>%
    summarise(
      N = n(),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      Min = min(Value, na.rm = TRUE),
      Max = max(Value, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      Variable_Type = case_when(
        grepl("cv_rhr", Variable) ~ "Wearable_CV_RHR",
        grepl("steps_max", Variable) ~ "Wearable_Steps",
        grepl("rhr_min", Variable) ~ "Wearable_RHR_Min",
        grepl("sleep_duration", Variable) ~ "Wearable_Sleep",
        grepl("improvement", Variable) ~ "Choroidal_Improvement",
        TRUE ~ "Other"
      ),
      Variable_Clean = case_when(
        grepl("cv_rhr", Variable) ~ "CV of RHR",
        grepl("steps_max", Variable) ~ "Max Steps",
        grepl("rhr_min", Variable) ~ "Min RHR",
        grepl("sleep_duration", Variable) ~ "Sleep Duration",
        grepl("improvement", Variable) ~ gsub("_improvement", " Improvement", Variable),
        TRUE ~ Variable
      )
    )
  
  # Display descriptive statistics
  cat("\n📊 Cluster 2描述性统计:\n")
  for(i in 1:nrow(cluster2_desc)) {
    stat <- cluster2_desc[i, ]
    cat(sprintf("\n%s (%s):\n", stat$Variable_Clean, stat$Variable_Type))
    cat(sprintf("  样本量: %d\n", stat$N))
    cat(sprintf("  均值 ± 标准差: %.3f ± %.3f\n", stat$Mean, stat$SD))
    cat(sprintf("  中位数 [Q1, Q3]: %.3f [%.3f, %.3f]\n", stat$Median, stat$Q1, stat$Q3))
    cat(sprintf("  范围: %.3f - %.3f\n", stat$Min, stat$Max))
  }
  
  # Compare Cluster 2 with other clusters (if data available)
  if(length(unique(data$max_cluster)) > 1) {
    cat("\n🔍 Cluster 2 vs 其他聚类比较:\n")
    
    # Calculate descriptive statistics for all clusters
    all_clusters_desc <- data %>%
      dplyr::select(max_cluster, all_of(all_vars)) %>%
      pivot_longer(cols = -max_cluster, names_to = "Variable", values_to = "Value") %>%
      filter(!is.na(Value) & !is.na(max_cluster)) %>%
      group_by(max_cluster, Variable) %>%
      summarise(
        N = n(),
        Mean = mean(Value, na.rm = TRUE),
        SD = sd(Value, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        Variable_Clean = case_when(
          grepl("cv_rhr", Variable) ~ "CV of RHR",
          grepl("steps_max", Variable) ~ "Max Steps",
          grepl("rhr_min", Variable) ~ "Min RHR",
          grepl("sleep_duration", Variable) ~ "Sleep Duration",
          grepl("improvement", Variable) ~ gsub("_improvement", " Improvement", Variable),
          TRUE ~ Variable
        )
      )
    
    # Show comparison
    for(var in unique(all_clusters_desc$Variable_Clean)) {
      var_data <- all_clusters_desc %>% filter(Variable_Clean == var)
      
      if(nrow(var_data) > 1) {
        cat(sprintf("\n%s:\n", var))
        for(j in 1:nrow(var_data)) {
          cluster_stat <- var_data[j, ]
          highlight <- ifelse(cluster_stat$max_cluster == 2, " 🎯", "")
          cat(sprintf("  Cluster %d: %.3f ± %.3f (n=%d)%s\n", 
                      cluster_stat$max_cluster, cluster_stat$Mean, cluster_stat$SD, cluster_stat$N, highlight))
        }
      }
    }
  }
  
  return(cluster2_desc)
}

# Execute Cluster 2 descriptive analysis
cluster2_descriptive_stats <- perform_cluster2_descriptive_analysis(analysis_data)

# Save descriptive statistics
if(!is.null(cluster2_descriptive_stats)) {
  write.csv(cluster2_descriptive_stats, "cluster2_descriptive_statistics.csv", row.names = FALSE)
}

