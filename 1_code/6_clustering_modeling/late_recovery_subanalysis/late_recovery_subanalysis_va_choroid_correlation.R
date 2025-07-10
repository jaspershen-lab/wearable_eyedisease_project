# ----------------------------------------------------
# Late Recovery Cluster PA Choroid vs Visual Acuity Analysis
# 分析late recovery聚类中PA Choroid与视力VA的相关性
# ----------------------------------------------------

# Load required libraries
library(tidyverse)
library(ggplot2)
library(corrplot)
library(gridExtra)
library(r4projects)
library(pheatmap)

# Set working directory
setwd(get_project_wd())
rm(list = ls())

# ================== 1. DATA LOADING AND SETUP ==================

cat("=== 1. 数据加载和设置 ===\n")

# Load baseline information and OCTA data
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")



# Load late recovery clustering results
late_recovery_file <- "3_data_analysis/6_clustering_modeling/time_window_clustering/late_recovery_membership_data.csv"

# Safe loading function
load_data_safely <- function(file_path, data_name) {
  if(file.exists(file_path)) {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat(sprintf("✓ 成功加载 %s: %d 行数据\n", data_name, nrow(data)))
    return(data)
  } else {
    cat(sprintf("⚠️ 文件不存在: %s\n", file_path))
    cat("将尝试查找其他可能的聚类文件...\n")
    
    # Try alternative paths
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

# Load late recovery clustering data
late_recovery_clusters <- load_data_safely(late_recovery_file, "Late Recovery聚类数据")

if(is.null(late_recovery_clusters)) {
  cat("❌ 无法找到聚类数据文件\n")
  stop("缺少必要的聚类数据")
} else {
  cat(sprintf("聚类数据预览:\n"))
  print(head(late_recovery_clusters))
  cat(sprintf("可用聚类: %s\n", paste(sort(unique(late_recovery_clusters$max_cluster)), collapse = ", ")))
}

# Create output directory
dir.create("3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/late_recovery_pa_choroid_va_analysis", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/late_recovery_pa_choroid_va_analysis")

# Get PPV patients
ppv_patients <- baseline_info %>%
  filter(surgery_1..0.PI.1.other. == 1) %>%
  distinct(ID) %>%
  pull(ID)

cat("PPV患者数:", length(ppv_patients), "\n")

# ================== 2. PROCESS OCTA AND VA DATA ==================

cat("\n=== 2. 处理OCTA和VA数据 ===\n")

# Process OCTA data function (from reference code)
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

# Process VA data using your definition
process_va_data_custom <- function(baseline_data) {
  
  cat("处理视力数据 (使用自定义vision_improvement_1m)...\n")
  
  # Create vision data using your exact definition
  vision_data <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    mutate(
      pre_vision = case_when(
        surgery_eye_1 == 0 ~ od_corrected_bas,  # Right eye surgery
        surgery_eye_1 == 1 ~ os_corrected_bas,  # Left eye surgery
        surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,  # Both eyes (average)
        TRUE ~ NA_real_
      ),
      post_vision_1w = case_when(
        surgery_eye_1 == 0 ~ od_corrected_1w,   # Right eye post-surgery 1 week
        surgery_eye_1 == 1 ~ os_corrected_1w,   # Left eye post-surgery 1 week
        surgery_eye_1 == 2 ~ (od_corrected_1w + os_corrected_1w)/2,   # Both eyes post-surgery 1 week
        TRUE ~ NA_real_
      ),
      post_vision_1m = case_when(
        surgery_eye_1 == 0 ~ od_corrected_1m,   # Right eye post-surgery 1 month
        surgery_eye_1 == 1 ~ os_corrected_1m,   # Left eye post-surgery 1 month
        surgery_eye_1 == 2 ~ (od_corrected_1m + os_corrected_1m)/2,   # Both eyes post-surgery 1 month
        TRUE ~ NA_real_
      ),
      vision_improvement_1w = post_vision_1w - pre_vision,
      vision_improvement_1m = post_vision_1m - pre_vision,
      vision_improved_1m = if_else(vision_improvement_1m >= 0, 1, 0),  # Binary improvement indicator
      vision_improved_factor_1m = factor(vision_improved_1m, 
                                         levels = c(0, 1), 
                                         labels = c("NoImprovement", "Improved"))  # Factor version
    ) %>%
    dplyr::select(ID, surgery_eye_1, pre_vision, post_vision_1w, post_vision_1m,
                  vision_improvement_1w, vision_improvement_1m, 
                  vision_improved_1m, vision_improved_factor_1m,
                  age, gender)
  
  return(vision_data)
}

# Process VA data for PPV patients using custom definition
ppv_va <- process_va_data_custom(baseline_info) %>%
  filter(ID %in% ppv_patients)

cat("匹配VA数据的PPV患者:", nrow(ppv_va), "\n")
cat("可用的vision_improvement_1m数据:", sum(!is.na(ppv_va$vision_improvement_1m)), "患者\n")

# ================== 3. EXTRACT PA CHOROID PARAMETERS ==================

cat("\n=== 3. 提取PA Choroid参数 ===\n")

# Filter PA Choroid parameters
filter_pa_choroid_parameters <- function(data) {
  
  cat("寻找PA Choroid参数...\n")
  
  # Look for PA Choroid parameters
  pa_choroid_patterns <- c(
    "PA.*Choroid.*0_21.*T0$",  # Wide field
    "PA.*Choroid.*0_6.*T0$",   # Macular
    "PA_Choroid.*T0$"          # General PA Choroid
  )
  
  all_pa_choroid_cols <- c()
  for(pattern in pa_choroid_patterns) {
    matches <- grep(pattern, names(data), value = TRUE, ignore.case = TRUE)
    all_pa_choroid_cols <- c(all_pa_choroid_cols, matches)
  }
  
  all_pa_choroid_cols <- unique(all_pa_choroid_cols)
  
  cat("找到的PA Choroid参数:\n")
  if(length(all_pa_choroid_cols) > 0) {
    for(i in 1:length(all_pa_choroid_cols)) {
      cat(sprintf("%d. %s\n", i, all_pa_choroid_cols[i]))
    }
  } else {
    cat("❌ 未找到PA Choroid参数\n")
    return(list(data = data, params = character(0)))
  }
  
  # Check corresponding T2 parameters
  params_T0 <- all_pa_choroid_cols
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  # Keep only parameters that have both T0 and T2
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  cat(sprintf("有效的PA Choroid参数对数: %d (T0和T2都存在)\n", length(valid_base_params)))
  
  if(length(valid_base_params) == 0) {
    cat("❌ 没有找到有效的PA Choroid参数对\n")
    return(list(data = data, params = character(0)))
  }
  
  # Filter data
  final_params_T0 <- paste0(valid_base_params, "_T0")
  final_params_T2 <- paste0(valid_base_params, "_T2")
  
  filtered_data <- data %>%
    dplyr::select(ID, all_of(c(final_params_T0, final_params_T2)))
  
  cat("筛选后的PA Choroid数据维度:", dim(filtered_data), "\n")
  
  return(list(
    data = filtered_data,
    params_T0 = final_params_T0,
    params_T2 = final_params_T2,
    base_params = valid_base_params
  ))
}

# Apply PA Choroid filtering
ppv_pa_choroid_filtered <- filter_pa_choroid_parameters(ppv_bloodflow)

# ================== 4. CALCULATE PA CHOROID IMPROVEMENTS ==================

cat("\n=== 4. 计算PA Choroid改善值 ===\n")

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

# Calculate PA Choroid improvements
if(length(ppv_pa_choroid_filtered$base_params) > 0) {
  pa_choroid_data <- calculate_improvement_correct(
    ppv_pa_choroid_filtered$data,
    ppv_pa_choroid_filtered$params_T0,
    ppv_pa_choroid_filtered$params_T2
  )
  
  cat("✓ 成功计算PA Choroid改善数据\n")
  cat("参数数量:", ncol(pa_choroid_data) - 1, "\n")
  cat("患者数量:", nrow(pa_choroid_data), "\n")
  
} else {
  cat("❌ 无法获取PA Choroid参数\n")
  stop("无法获取PA Choroid数据")
}

# ================== 5. EXTRACT VA PARAMETERS (SIMPLIFIED) ==================

cat("\n=== 5. 提取视力VA参数 (使用vision_improvement_1m) ===\n")

# Since we're using the custom vision improvement definition, we don't need complex filtering
# We already have vision_improvement_1m as our target VA parameter

va_params <- c("vision_improvement_1m")
cat("使用的VA参数: vision_improvement_1m\n")
cat("有效vision_improvement_1m数据:", sum(!is.na(ppv_va$vision_improvement_1m)), "患者\n")

# Create a simplified VA data structure for consistency with the rest of the code
va_data_processed <- ppv_va %>%
  dplyr::select(ID, vision_improvement_1m) %>%
  filter(!is.na(vision_improvement_1m))

cat("✓ 成功处理VA数据\n")
cat("参数数量: 1 (vision_improvement_1m)\n")
cat("患者数量:", nrow(va_data_processed), "\n")

# ================== 6. INTEGRATE ALL DATA ==================

cat("\n=== 6. 整合所有数据 ===\n")

integrate_cluster_pa_choroid_va_data <- function(pa_choroid_data, va_data, late_recovery_clusters) {
  
  cat("整合聚类、PA Choroid和VA数据...\n")
  
  # Check data availability
  cat("PA Choroid数据列:", paste(names(pa_choroid_data), collapse = ", "), "\n")
  cat("VA数据列:", paste(names(va_data), collapse = ", "), "\n")
  
  # Standardize cluster data ID column
  cluster_id_col <- "subject_id"
  if(!"subject_id" %in% names(late_recovery_clusters)) {
    if("ID" %in% names(late_recovery_clusters)) {
      cluster_id_col <- "ID"
    } else {
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
  
  # Standardize cluster data
  late_recovery_clusters_std <- late_recovery_clusters
  if(cluster_id_col != "subject_id") {
    names(late_recovery_clusters_std)[names(late_recovery_clusters_std) == cluster_id_col] <- "subject_id"
  }
  
  # Integrate data step by step
  # First: PA Choroid + VA
  integrated_data <- pa_choroid_data %>%
    left_join(va_data, by = "ID")
  
  # Then: Add clustering information
  integrated_data <- integrated_data %>%
    left_join(late_recovery_clusters_std %>% 
                dplyr::select(subject_id, max_cluster, max_membership), 
              by = c("ID" = "subject_id"))
  
  cat(sprintf("整合后的数据维度: %d 行 × %d 列\n", nrow(integrated_data), ncol(integrated_data)))
  
  # Check cluster distribution
  if("max_cluster" %in% names(integrated_data)) {
    cluster_dist <- integrated_data %>%
      filter(!is.na(max_cluster)) %>%
      count(max_cluster, name = "n_patients") %>%
      mutate(percentage = round(n_patients/sum(n_patients)*100, 1))
    
    cat("\nLate recovery聚类分布:\n")
    print(cluster_dist)
  } else {
    cat("⚠️ 聚类信息整合失败\n")
  }
  
  # Get improvement parameters
  pa_choroid_improvement_cols <- names(integrated_data)[grep("PA.*Choroid.*improvement", names(integrated_data))]
  va_improvement_cols <- c("vision_improvement_1m")  # Using our specific VA parameter
  
  cat("\nPA Choroid改善参数:", length(pa_choroid_improvement_cols), "个\n")
  if(length(pa_choroid_improvement_cols) > 0) {
    cat("参数列表:", paste(pa_choroid_improvement_cols, collapse = ", "), "\n")
  }
  
  cat("\nVA改善参数: 1个 (vision_improvement_1m)\n")
  
  # Calculate data completeness
  if(length(pa_choroid_improvement_cols) > 0 && length(va_improvement_cols) > 0) {
    
    completeness_summary <- data.frame()
    
    for(pa_col in pa_choroid_improvement_cols) {
      for(va_col in va_improvement_cols) {
        complete_count <- integrated_data %>%
          filter(!is.na(!!sym(pa_col)) & !is.na(!!sym(va_col)) & !is.na(max_cluster)) %>%
          nrow()
        
        completeness_summary <- rbind(completeness_summary, data.frame(
          PA_Choroid_Param = pa_col,
          VA_Param = va_col,
          Complete_Cases = complete_count,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Select best parameter combinations
    best_combinations <- completeness_summary %>%
      filter(Complete_Cases >= 3) %>%
      arrange(desc(Complete_Cases))
    
    cat("\n数据完整性汇总:\n")
    print(best_combinations)
    
    if(nrow(best_combinations) == 0) {
      cat("❌ 没有足够的完整数据进行分析\n")
      return(NULL)
    }
    
    # Create final analysis dataset
    analysis_vars <- c("ID", "max_cluster", "max_membership",
                       unique(c(best_combinations$PA_Choroid_Param, best_combinations$VA_Param)))
    
    final_data <- integrated_data %>%
      dplyr::select(all_of(analysis_vars)) %>%
      filter(!is.na(max_cluster))
    
    cat(sprintf("\n最终分析数据: %d 患者\n", nrow(final_data)))
    
    return(list(
      data = final_data,
      pa_choroid_params = pa_choroid_improvement_cols,
      va_params = va_improvement_cols,
      best_combinations = best_combinations
    ))
    
  } else {
    cat("❌ 缺少PA Choroid或VA改善参数\n")
    return(NULL)
  }
}

# Integrate all data
analysis_result <- integrate_cluster_pa_choroid_va_data(pa_choroid_data, va_data_processed, late_recovery_clusters)

if(is.null(analysis_result)) {
  cat("❌ 数据整合失败\n")
  stop("无法整合分析数据")
} else {
  analysis_data <- analysis_result$data
  pa_choroid_params <- analysis_result$pa_choroid_params
  va_params <- analysis_result$va_params
  best_combinations <- analysis_result$best_combinations
  cat("✓ 数据整合成功\n")
}

# ================== 7. CORRELATION ANALYSIS ==================

cat("\n=== 7. PA Choroid与VA相关性分析 ===\n")

perform_pa_choroid_va_correlation <- function(data, pa_choroid_params, va_params) {
  
  cat("执行PA Choroid与VA相关性分析...\n")
  
  # Store correlation results
  correlation_results <- data.frame()
  
  # Perform correlation analysis for all combinations
  for(pa_col in pa_choroid_params) {
    for(va_col in va_params) {
      
      if(pa_col %in% names(data) && va_col %in% names(data)) {
        
        # Clean data
        clean_pair <- data %>%
          filter(!is.na(!!sym(pa_col)) & !is.na(!!sym(va_col)))
        
        if(nrow(clean_pair) >= 3) {
          
          # Pearson correlation
          pearson_test <- cor.test(clean_pair[[pa_col]], 
                                   clean_pair[[va_col]], 
                                   method = "pearson")
          
          # Spearman correlation
          spearman_test <- cor.test(clean_pair[[pa_col]], 
                                    clean_pair[[va_col]], 
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
            PA_Choroid_Metric = pa_col,
            VA_Metric = va_col,
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
    
    cat("PA Choroid与VA相关性分析结果:\n")
    for(i in 1:nrow(correlation_results)) {
      result <- correlation_results[i, ]
      significance <- ifelse(result$Highly_Significant, "**", 
                             ifelse(result$Significant, "*", ""))
      cat(sprintf("%d. %s vs %s:\n", i, 
                  result$PA_Choroid_Metric, result$VA_Metric))
      cat(sprintf("   r = %.3f, p = %.4f%s (%s effect, n = %d)\n",
                  result$Pearson_r, result$Pearson_p, significance,
                  result$Effect_Size, result$N))
    }
  }
  
  return(correlation_results)
}

# Execute overall correlation analysis
overall_correlation_results <- perform_pa_choroid_va_correlation(analysis_data, pa_choroid_params, va_params)

# ================== 8. CLUSTER-SPECIFIC ANALYSIS ==================

cat("\n=== 8. 聚类特异性分析 ===\n")

perform_cluster_specific_analysis <- function(data, pa_choroid_params, va_params) {
  
  cat("执行基于聚类的分层分析...\n")
  
  # Get unique clusters
  unique_clusters <- unique(data$max_cluster)
  unique_clusters <- unique_clusters[!is.na(unique_clusters)]
  
  if(length(unique_clusters) == 0) {
    cat("❌ 没有有效的聚类信息\n")
    return(data.frame())
  }
  
  cat(sprintf("分析聚类: %s\n", paste(unique_clusters, collapse = ", ")))
  
  # Store cluster-specific correlation results
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
    for(pa_col in pa_choroid_params) {
      for(va_col in va_params) {
        
        if(pa_col %in% names(cluster_data) && va_col %in% names(cluster_data)) {
          
          # Clean data for this cluster
          clean_pair <- cluster_data %>%
            filter(!is.na(!!sym(pa_col)) & !is.na(!!sym(va_col)))
          
          if(nrow(clean_pair) >= 3) {
            
            # Pearson correlation
            pearson_test <- try(cor.test(clean_pair[[pa_col]], 
                                         clean_pair[[va_col]], 
                                         method = "pearson"), silent = TRUE)
            
            if(class(pearson_test) != "try-error") {
              
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
                PA_Choroid_Metric = pa_col,
                VA_Metric = va_col,
                N = nrow(clean_pair),
                Pearson_r = as.numeric(pearson_test$estimate),
                Pearson_p = pearson_test$p.value,
                Pearson_CI_Lower = pearson_test$conf.int[1],
                Pearson_CI_Upper = pearson_test$conf.int[2],
                Effect_Size = effect_size,
                Significant = pearson_test$p.value < 0.05,
                Highly_Significant = pearson_test$p.value < 0.01,
                stringsAsFactors = FALSE
              )
              
              cluster_results <- rbind(cluster_results, result_row)
              
              # Display significant results
              if(pearson_test$p.value < 0.05) {
                cat(sprintf("✓ 显著相关性: %s vs %s (r=%.3f, p=%.4f)\n",
                            pa_col, va_col, 
                            pearson_test$estimate, pearson_test$p.value))
              }
            }
          }
        }
      }
    }
  }
  
  # Summary of cluster-specific results
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

# Execute cluster-specific analysis
cluster_specific_results <- perform_cluster_specific_analysis(analysis_data, pa_choroid_params, va_params)

# ================== 9. CLUSTER COMPARISON ANALYSIS ==================

cat("\n=== 9. 聚类间比较分析 ===\n")

perform_cluster_comparison_analysis <- function(data, pa_choroid_params, va_params) {
  
  cat("执行聚类间PA Choroid和VA差异分析...\n")
  
  # Calculate descriptive statistics by cluster
  cluster_descriptives <- data.frame()
  
  unique_clusters <- sort(unique(data$max_cluster[!is.na(data$max_cluster)]))
  
  for(cluster_id in unique_clusters) {
    cluster_data <- data %>% filter(max_cluster == cluster_id)
    
    for(param in c(pa_choroid_params, va_params)) {
      if(param %in% names(cluster_data)) {
        
        param_stats <- cluster_data %>%
          filter(!is.na(!!sym(param))) %>%
          summarise(
            N = n(),
            Mean = mean(!!sym(param), na.rm = TRUE),
            SD = sd(!!sym(param), na.rm = TRUE),
            Median = median(!!sym(param), na.rm = TRUE),
            Q1 = quantile(!!sym(param), 0.25, na.rm = TRUE),
            Q3 = quantile(!!sym(param), 0.75, na.rm = TRUE),
            .groups = 'drop'
          ) %>%
          mutate(
            Cluster = cluster_id,
            Parameter = param,
            Parameter_Type = ifelse(grepl("PA.*Choroid", param), "PA_Choroid", "VA")
          )
        
        cluster_descriptives <- rbind(cluster_descriptives, param_stats)
      }
    }
  }
  
  # Display cluster differences
  cat("\n聚类间参数差异:\n")
  
  for(param_type in c("PA_Choroid", "VA")) {
    type_data <- cluster_descriptives %>% filter(Parameter_Type == param_type)
    
    if(nrow(type_data) > 0) {
      cat(sprintf("\n=== %s Parameters ===\n", param_type))
      
      for(param in unique(type_data$Parameter)) {
        param_data <- type_data %>% filter(Parameter == param)
        
        if(nrow(param_data) > 1) {
          cat(sprintf("\n%s:\n", param))
          for(i in 1:nrow(param_data)) {
            stat <- param_data[i, ]
            cat(sprintf("  Cluster %d: %.3f ± %.3f (n=%d)\n", 
                        stat$Cluster, stat$Mean, stat$SD, stat$N))
          }
          
          # Perform Kruskal-Wallis test if more than 2 clusters
          if(nrow(param_data) >= 2) {
            param_values <- data %>%
              filter(!is.na(max_cluster) & !is.na(!!sym(param))) %>%
              dplyr::select(max_cluster, !!sym(param))
            
            if(nrow(param_values) >= 5) {
              kw_test <- try(kruskal.test(param_values[[param]], param_values$max_cluster), silent = TRUE)
              if(class(kw_test) != "try-error") {
                cat(sprintf("  Kruskal-Wallis p = %.4f\n", kw_test$p.value))
              }
            }
          }
        }
      }
    }
  }
  
  return(cluster_descriptives)
}

# Execute cluster comparison analysis
cluster_descriptives <- perform_cluster_comparison_analysis(analysis_data, pa_choroid_params, va_params)

# ================== 10. VISUALIZATIONS ==================

cat("\n=== 10. 创建可视化 ===\n")

# Create scatter plots for significant correlations
create_pa_choroid_va_visualizations <- function(data, correlation_results) {
  
  cat("创建PA Choroid与VA相关性可视化...\n")
  
  # Create output directory
  dir.create("pa_choroid_va_analysis", recursive = TRUE, showWarnings = FALSE)
  setwd("pa_choroid_va_analysis")
  
  plot_list <- list()
  
  # Create scatter plots for significant correlations
  significant_results <- correlation_results %>%
    filter(Significant == TRUE) %>%
    arrange(desc(Abs_Pearson_r))
  
  if(nrow(significant_results) > 0) {
    
    for(i in 1:nrow(significant_results)) {
      result <- significant_results[i, ]
      
      pa_var <- result$PA_Choroid_Metric
      va_var <- result$VA_Metric
      
      # Create scatter plot
      p <- ggplot(data, aes_string(x = pa_var, y = va_var)) +
        geom_point(aes(color = factor(max_cluster), size = max_membership),
                   alpha = 0.8) +
        geom_smooth(method = "lm", se = TRUE, color = "black", 
                    linetype = "dashed", alpha = 0.7) +
        geom_text(aes(label = ID), vjust = -0.5, size = 2.5, alpha = 0.7) +
        scale_color_brewer(type = "qual", palette = "Set2", 
                           name = "Late Recovery\nCluster") +
        scale_size_continuous(range = c(2, 5), name = "Max\nMembership") +
        labs(
          title = paste("PA Choroid vs Visual Acuity Correlation"),
          subtitle = paste0(pa_var, " vs ", va_var, " | ",
                            "r = ", round(result$Pearson_r, 3), 
                            ", p = ", format.pval(result$Pearson_p, digits = 3),
                            " (", result$Effect_Size, " effect)"),
          x = "PA Choroid Improvement",
          y = "Visual Acuity Improvement",
          caption = paste("n =", result$N, "patients | Colors by late recovery cluster")
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.position = "right",
          plot.caption = element_text(size = 9)
        )
      
      plot_list[[i]] <- p
      
      # Save individual plot
      ggsave(paste0("pa_choroid_va_correlation_", i, ".pdf"),
             p, width = 12, height = 8)
    }
  } else {
    cat("❌ 没有显著相关性可以可视化\n")
  }
  
  # Create correlation heatmap if we have results
  if(nrow(correlation_results) > 0) {
    
    # Prepare correlation matrix data
    cor_matrix_data <- correlation_results %>%
      dplyr::select(PA_Choroid_Metric, VA_Metric, Pearson_r) %>%
      mutate(
        PA_Clean = gsub("_improvement", "", PA_Choroid_Metric),
        VA_Clean = gsub("_improvement", "", VA_Metric)
      ) %>%
      dplyr::select(PA_Clean, VA_Clean, Pearson_r) %>%
      pivot_wider(names_from = VA_Clean, values_from = Pearson_r) %>%
      column_to_rownames("PA_Clean") %>%
      as.matrix()
    
    # Create heatmap
    pdf("pa_choroid_va_correlation_heatmap.pdf", width = 10, height = 8)
    corrplot(cor_matrix_data, 
             method = "color", 
             type = "full",
             order = "hclust",
             tl.cex = 1.2,
             tl.col = "black",
             cl.cex = 1.2,
             addCoef.col = "black",
             number.cex = 1.5,
             title = "PA Choroid vs Visual Acuity Correlations",
             mar = c(0,0,2,0))
    dev.off()
  }
  
  # Combine significant result plots
  if(length(plot_list) > 0) {
    combined_plot <- do.call(gridExtra::grid.arrange, 
                             c(plot_list, ncol = min(2, length(plot_list)),
                               top = "PA Choroid vs Visual Acuity Correlations"))
    
    ggsave("pa_choroid_va_combined_correlations.pdf",
           combined_plot, width = 16, height = 8 * ceiling(length(plot_list)/2))
  }
  
  # Return to parent directory
  setwd("..")
  
  cat("✓ PA Choroid vs VA可视化已保存到 pa_choroid_va_analysis/ 目录\n")
  
  return(plot_list)
}

# Create cluster-specific visualizations
create_cluster_specific_visualizations <- function(data, cluster_results) {
  
  cat("创建聚类特异性可视化...\n")
  
  # Create cluster-specific output directory
  dir.create("cluster_specific_pa_choroid_va", recursive = TRUE, showWarnings = FALSE)
  setwd("cluster_specific_pa_choroid_va")
  
  cluster_plots <- list()
  
  # Get significant cluster-specific results
  significant_cluster_results <- cluster_results %>%
    filter(Significant == TRUE) %>%
    arrange(Cluster, desc(Abs_Pearson_r))
  
  if(nrow(significant_cluster_results) > 0) {
    
    for(i in 1:nrow(significant_cluster_results)) {
      result <- significant_cluster_results[i, ]
      
      cluster_id <- result$Cluster
      pa_var <- result$PA_Choroid_Metric
      va_var <- result$VA_Metric
      
      # Filter data for this cluster
      cluster_data <- data %>% filter(max_cluster == cluster_id)
      
      # Create cluster-specific scatter plot
      p <- ggplot(cluster_data, aes_string(x = pa_var, y = va_var)) +
        geom_point(size = 4, alpha = 0.8, color = "darkblue") +
        geom_smooth(method = "lm", se = TRUE, color = "red", 
                    linetype = "solid", alpha = 0.7, size = 1.2) +
        geom_text(aes(label = ID), vjust = -0.8, size = 3, alpha = 0.8, color = "black") +
        labs(
          title = paste0("Cluster ", cluster_id, ": PA Choroid vs Visual Acuity"),
          subtitle = paste0("Cluster ", cluster_id, " Only (n = ", result$N, ") | ",
                            pa_var, " vs ", va_var, "\n",
                            "r = ", round(result$Pearson_r, 3), 
                            ", p = ", format.pval(result$Pearson_p, digits = 3),
                            " (", result$Effect_Size, " effect)",
                            ifelse(result$Significant, " *", "")),
          x = "PA Choroid Improvement",
          y = "Visual Acuity Improvement",
          caption = paste("Cluster", cluster_id, "Patients Only | 95% CI:", 
                          round(result$Pearson_CI_Lower, 3), "to", 
                          round(result$Pearson_CI_Upper, 3))
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "darkblue"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          plot.caption = element_text(size = 9),
          panel.border = element_rect(color = "darkblue", size = 1.5)
        )
      
      cluster_plots[[i]] <- p
      
      # Save individual plot
      ggsave(paste0("cluster_", cluster_id, "_pa_choroid_va_", i, ".pdf"),
             p, width = 12, height = 8)
    }
  }
  
  # Create cluster comparison boxplot
  if(nrow(cluster_descriptives) > 0) {
    
    # Prepare data for boxplot
    cluster_long_data <- data %>%
      dplyr::select(ID, max_cluster, all_of(c(pa_choroid_params, va_params))) %>%
      filter(!is.na(max_cluster)) %>%
      pivot_longer(cols = -c(ID, max_cluster), names_to = "Parameter", values_to = "Value") %>%
      filter(!is.na(Value)) %>%
      mutate(
        Parameter_Type = ifelse(grepl("PA.*Choroid", Parameter), "PA Choroid", "Visual Acuity"),
        Parameter_Clean = gsub("_improvement", " Improvement", Parameter)
      )
    
    # Box plot by cluster
    p_cluster_box <- ggplot(cluster_long_data, aes(x = factor(max_cluster), y = Value, fill = factor(max_cluster))) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
      facet_wrap(~ Parameter_Clean, scales = "free_y") +
      scale_fill_brewer(type = "qual", palette = "Set2", name = "Late Recovery\nCluster") +
      labs(
        title = "PA Choroid and Visual Acuity by Late Recovery Cluster",
        subtitle = "Distribution comparison across clusters",
        x = "Late Recovery Cluster",
        y = "Improvement Value",
        caption = "Points show individual patients | Box shows median and quartiles"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10),
        legend.position = "bottom"
      )
    
    ggsave("cluster_comparison_boxplot.pdf", p_cluster_box, width = 14, height = 10)
  }
  
  # Return to parent directory
  setwd("..")
  
  cat("✓ 聚类特异性可视化已保存到 cluster_specific_pa_choroid_va/ 目录\n")
  
  return(cluster_plots)
}

# Execute visualizations
overall_visualizations <- create_pa_choroid_va_visualizations(analysis_data, overall_correlation_results)
cluster_visualizations <- create_cluster_specific_visualizations(analysis_data, cluster_specific_results)

# ================== 11. SAVE RESULTS ==================

cat("\n=== 11. 保存分析结果 ===\n")

# Save all results
write.csv(analysis_data, "integrated_late_recovery_pa_choroid_va_data.csv", row.names = FALSE)
write.csv(overall_correlation_results, "pa_choroid_va_overall_correlations.csv", row.names = FALSE)

if(nrow(cluster_specific_results) > 0) {
  write.csv(cluster_specific_results, "pa_choroid_va_cluster_specific_correlations.csv", row.names = FALSE)
}

if(nrow(cluster_descriptives) > 0) {
  write.csv(cluster_descriptives, "cluster_descriptive_statistics.csv", row.names = FALSE)
}

# Save clustering information for reference
write.csv(late_recovery_clusters, "late_recovery_clustering_reference.csv", row.names = FALSE)

cat("✓ 分析结果已保存\n")

# ================== 12. GENERATE COMPREHENSIVE REPORT ==================

cat("\n=== 12. 生成综合报告 ===\n")

generate_comprehensive_report <- function(analysis_data, overall_results, cluster_results, cluster_descriptives) {
  
  report <- paste0(
    "========================================\n",
    "Late Recovery聚类中PA Choroid与视力VA相关性分析报告\n", 
    "========================================\n\n",
    
    "🔬 研究背景:\n",
    "基于late recovery时间窗口聚类分析，探索不同恢复模式患者中\n",
    "脉络膜血流灌注(PA Choroid)与视力改善(VA)的关联性。\n\n",
    
    "📊 数据概览:\n",
    "- 分析患者数: ", nrow(analysis_data), "\n",
    "- PA Choroid参数: ", length(pa_choroid_params), " 个\n",
    "- VA参数: ", length(va_params), " 个\n"
  )
  
  # Add cluster information
  if("max_cluster" %in% names(analysis_data)) {
    cluster_dist <- analysis_data %>% 
      filter(!is.na(max_cluster)) %>%
      count(max_cluster) %>%
      arrange(max_cluster)
    
    report <- paste0(report, "- 聚类分布: ")
    for(i in 1:nrow(cluster_dist)) {
      report <- paste0(report, "Cluster ", cluster_dist$max_cluster[i], " (n=", cluster_dist$n[i], ")")
      if(i < nrow(cluster_dist)) report <- paste0(report, ", ")
    }
    report <- paste0(report, "\n")
  }
  
  report <- paste0(report, "\n🎯 主要发现:\n")
  
  # Overall correlation results
  if(nrow(overall_results) > 0) {
    
    significant_results <- overall_results %>% filter(Significant == TRUE)
    
    if(nrow(significant_results) > 0) {
      report <- paste0(report, "✅ 整体发现 ", nrow(significant_results), " 个显著相关性:\n\n")
      
      for(i in 1:nrow(significant_results)) {
        result <- significant_results[i, ]
        report <- paste0(report,
                         sprintf("%d. %s vs %s:\n",
                                 i, result$PA_Choroid_Metric, result$VA_Metric),
                         sprintf("   相关系数: r = %.3f\n", result$Pearson_r),
                         sprintf("   显著性: p = %.4f\n", result$Pearson_p),
                         sprintf("   效应大小: %s\n", result$Effect_Size),
                         sprintf("   样本量: n = %d\n\n", result$N))
      }
    } else {
      report <- paste0(report, "❌ 未发现显著的整体相关性\n\n")
    }
  }
  
  # Cluster-specific results
  if(nrow(cluster_results) > 0) {
    cluster_significant <- cluster_results %>% filter(Significant == TRUE)
    
    if(nrow(cluster_significant) > 0) {
      report <- paste0(report, "🎯 聚类特异性相关性发现:\n")
      
      for(cluster in unique(cluster_significant$Cluster)) {
        cluster_results_subset <- cluster_significant %>% filter(Cluster == cluster)
        
        if(nrow(cluster_results_subset) > 0) {
          report <- paste0(report, sprintf("\n聚类 %d (%d 个显著相关性):\n", 
                                           cluster, nrow(cluster_results_subset)))
          
          for(j in 1:nrow(cluster_results_subset)) {
            result <- cluster_results_subset[j, ]
            report <- paste0(report,
                             sprintf("  %d. %s vs %s:\n",
                                     j, result$PA_Choroid_Metric, result$VA_Metric),
                             sprintf("     r = %.3f, p = %.4f (%s effect)\n",
                                     result$Pearson_r, result$Pearson_p, result$Effect_Size))
          }
        }
      }
      report <- paste0(report, "\n")
    } else {
      report <- paste0(report, "❌ 未发现聚类特异性相关性\n\n")
    }
  }
  
  # Cluster differences summary
  if(nrow(cluster_descriptives) > 0) {
    report <- paste0(report, "📊 聚类间差异概览:\n")
    
    # PA Choroid differences
    pa_choroid_stats <- cluster_descriptives %>% 
      filter(Parameter_Type == "PA_Choroid") %>%
      group_by(Parameter) %>%
      summarise(
        mean_diff = max(Mean, na.rm = TRUE) - min(Mean, na.rm = TRUE),
        .groups = 'drop'
      )
    
    if(nrow(pa_choroid_stats) > 0) {
      report <- paste0(report, "- PA Choroid参数聚类间差异较大的指标:\n")
      top_pa_diff <- pa_choroid_stats %>% 
        arrange(desc(mean_diff)) %>% 
        slice_head(n = 3)
      
      for(k in 1:nrow(top_pa_diff)) {
        report <- paste0(report, sprintf("  %d. %s (差异: %.3f)\n", 
                                         k, top_pa_diff$Parameter[k], top_pa_diff$mean_diff[k]))
      }
    }
    
    # VA differences
    va_stats <- cluster_descriptives %>% 
      filter(Parameter_Type == "VA") %>%
      group_by(Parameter) %>%
      summarise(
        mean_diff = max(Mean, na.rm = TRUE) - min(Mean, na.rm = TRUE),
        .groups = 'drop'
      )
    
    if(nrow(va_stats) > 0) {
      report <- paste0(report, "- VA参数聚类间差异较大的指标:\n")
      top_va_diff <- va_stats %>% 
        arrange(desc(mean_diff)) %>% 
        slice_head(n = 3)
      
      for(k in 1:nrow(top_va_diff)) {
        report <- paste0(report, sprintf("  %d. %s (差异: %.3f)\n", 
                                         k, top_va_diff$Parameter[k], top_va_diff$mean_diff[k]))
      }
    }
    report <- paste0(report, "\n")
  }
  
  report <- paste0(report,
                   "💡 临床意义:\n",
                   "1. PA Choroid血流灌注改善与视力恢复可能存在关联\n",
                   "2. 不同late recovery聚类显示不同的血流-视力关系模式\n",
                   "3. 脉络膜血流可能是视力恢复的重要影响因素\n",
                   "4. 聚类特异性相关性提示个体化治疗策略的重要性\n\n",
                   
                   "📈 生成文件:\n",
                   "- integrated_late_recovery_pa_choroid_va_data.csv: 整合数据\n",
                   "- pa_choroid_va_overall_correlations.csv: 整体相关性结果\n",
                   "- pa_choroid_va_cluster_specific_correlations.csv: 聚类特异性结果\n",
                   "- cluster_descriptive_statistics.csv: 聚类描述性统计\n",
                   "- pa_choroid_va_analysis/ 目录: 整体分析可视化\n",
                   "- cluster_specific_pa_choroid_va/ 目录: 聚类特异性分析\n\n",
                   
                   "🔬 下一步建议:\n",
                   "1. 验证PA Choroid与VA关联的因果关系\n",
                   "2. 探索聚类特异性机制的生物学基础\n",
                   "3. 开发基于血流灌注的视力预后预测模型\n",
                   "4. 考虑将脉络膜血流作为治疗监测指标\n\n",
                   
                   "分析完成时间: ", Sys.time(), "\n",
                   "========================================\n")
  
  # Save report
  writeLines(report, "pa_choroid_va_comprehensive_report.txt")
  cat(report)
  
  return(report)
}

# Generate comprehensive report
final_report <- generate_comprehensive_report(analysis_data, overall_correlation_results, 
                                              cluster_specific_results, cluster_descriptives)

# ================== 13. FINAL SUMMARY ==================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("🎯 Late Recovery聚类PA Choroid与视力VA相关性分析 - 完整总结\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Display main findings
cat("📊 主要发现汇总:\n")
cat("1. 数据质量: 成功整合", nrow(analysis_data), "患者的完整数据\n")

if("max_cluster" %in% names(analysis_data)) {
  cluster_count <- length(unique(analysis_data$max_cluster[!is.na(analysis_data$max_cluster)]))
  cat("2. 聚类数: ", cluster_count, "个late recovery聚类\n")
}

cat("3. PA Choroid参数: ", length(pa_choroid_params), "个\n")
cat("4. VA参数: ", length(va_params), "个\n")

if(nrow(overall_correlation_results) > 0) {
  significant_count <- sum(overall_correlation_results$Significant)
  cat("5. 整体相关性: ", significant_count, "/", nrow(overall_correlation_results), "个显著相关性\n")
  
  if(significant_count > 0) {
    best_correlation <- overall_correlation_results %>% 
      filter(Significant == TRUE) %>% 
      slice_max(abs(Pearson_r), n = 1)
    
    cat("6. 最强相关性: ", best_correlation$PA_Choroid_Metric[1], 
        " vs ", best_correlation$VA_Metric[1], 
        " (r = ", round(best_correlation$Pearson_r[1], 3), ")\n")
  }
}

if(nrow(cluster_specific_results) > 0) {
  cluster_significant_count <- sum(cluster_specific_results$Significant)
  cat("7. 聚类特异性相关性: ", cluster_significant_count, "个显著发现\n")
}

cat("\n📁 生成的主要分析文件:\n")
output_files_main <- c(
  "integrated_late_recovery_pa_choroid_va_data.csv",
  "pa_choroid_va_overall_correlations.csv", 
  "pa_choroid_va_cluster_specific_correlations.csv",
  "cluster_descriptive_statistics.csv",
  "late_recovery_clustering_reference.csv",
  "pa_choroid_va_comprehensive_report.txt"
)

for(file in output_files_main) {
  if(file.exists(file)) {
    cat("✅", file, "\n")
  } else {
    cat("❌", file, "(未生成)\n")
  }
}

cat("\n📁 生成的分析目录:\n")
analysis_dirs <- c("pa_choroid_va_analysis", "cluster_specific_pa_choroid_va")

for(dir in analysis_dirs) {
  if(dir.exists(dir)) {
    cat("✅", dir, "/\n")
  } else {
    cat("❌", dir, "/ (未生成)\n")
  }
}

cat("\n🌟 关键发现:\n")
if(nrow(overall_correlation_results) > 0) {
  if(sum(overall_correlation_results$Significant) > 0) {
    cat("- ✅ PA Choroid与VA存在显著相关性\n")
  } else {
    cat("- ❌ 整体上PA Choroid与VA无显著相关性\n")
  }
}

if(nrow(cluster_specific_results) > 0) {
  if(sum(cluster_specific_results$Significant) > 0) {
    cat("- ✅ 特定聚类中发现PA Choroid与VA显著相关性\n")
    
    # Show which clusters have significant correlations
    sig_clusters <- cluster_specific_results %>% 
      filter(Significant == TRUE) %>% 
      distinct(Cluster) %>% 
      pull(Cluster)
    
    if(length(sig_clusters) > 0) {
      cat("- 🎯 发现显著相关性的聚类:", paste(sig_clusters, collapse = ", "), "\n")
    }
  } else {
    cat("- ❌ 各聚类中均无PA Choroid与VA显著相关性\n")
  }
}

cat("\n🔍 研究意义:\n")
cat("- PA Choroid血流灌注与视力改善的关系\n")
cat("- late recovery聚类模式的临床预测价值\n")
cat("- 个体化治疗策略的数据支持\n")
cat("- 脉络膜血流作为预后指标的潜力\n")

cat("\n✅ Late Recovery聚类PA Choroid与视力VA相关性分析完成！\n")
cat("📊 使用vision_improvement_1m作为核心VA指标\n")
cat("📁 所有结果已保存至相应文件和目录\n")
cat(paste(rep("=", 80), collapse = ""), "\n")