# 基于时间窗口硬聚类的组间差异分析
# 针对所有membership=1的硬聚类情况，使用cluster比较分析代替相关性分析

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(corrplot)
library(r4projects)

# ================== 1. 设置工作目录并读取数据 ==================

setwd(get_project_wd())
rm(list = ls())

# 读取代码三生成的max membership数据
max_membership_data <- read.csv("3_data_analysis/6_clustering_modeling/time_window_clustering/time_window_2cluster_membership_data.csv")

cat("✓ 读取时间窗Max Membership数据\n")
cat("数据维度:", dim(max_membership_data), "\n")
cat("列名:", paste(names(max_membership_data), collapse = ", "), "\n\n")

# 检查数据结构
membership_cols <- grep("^membership_", names(max_membership_data), value = TRUE)
cluster_cols <- grep("^cluster_", names(max_membership_data), value = TRUE)

cat("🎯 发现的数据列:\n")
cat("Membership列:", paste(membership_cols, collapse = ", "), "\n")
cat("Cluster列:", paste(cluster_cols, collapse = ", "), "\n\n")

# 检查是否为硬聚类
sample_memberships <- max_membership_data %>%
  dplyr::select(all_of(membership_cols)) %>%
  slice_head(n = 5)

cat("🔍 样本membership值检查:\n")
print(sample_memberships)

all_ones <- all(sapply(max_membership_data[membership_cols], function(x) all(x == 1, na.rm = TRUE)))
cat(sprintf("📊 所有membership值都是1: %s\n", ifelse(all_ones, "是", "否")))

if(all_ones) {
  cat("🎯 检测到硬聚类，将使用cluster组间比较分析\n\n")
} else {
  cat("⚠️ 非硬聚类，可考虑使用传统相关性分析\n\n")
}

# ================== 2. 读取和处理OCTA数据 ==================

if(!exists("octa_improvements") || !exists("va_improvements")) {
  
  cat("===== 处理OCTA和VA数据 =====\n")
  
  # 读取基础数据
  baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
  octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
  octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")
  
  # OCTA数据处理函数
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
  
  cat("✓ OCTA data processed\n")
  
  # 处理VA数据
  process_vision_improvements <- function(baseline_data) {
    ppv_patients <- baseline_data %>%
      filter(surgery_1..0.PI.1.other. == 1) %>%
      distinct(ID) %>%
      pull(ID)
    
    vision_data <- baseline_data %>%
      filter(ID %in% ppv_patients & !is.na(surgery_eye_1)) %>%
      distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
      mutate(
        pre_vision = case_when(
          surgery_eye_1 == 0 ~ od_corrected_bas,
          surgery_eye_1 == 1 ~ os_corrected_bas,
          surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,
          TRUE ~ NA_real_
        ),
        post_vision_1w = case_when(
          surgery_eye_1 == 0 ~ od_corrected_1w,
          surgery_eye_1 == 1 ~ os_corrected_1w,
          surgery_eye_1 == 2 ~ (od_corrected_1w + os_corrected_1w)/2,
          TRUE ~ NA_real_
        ),
        post_vision_1m = case_when(
          surgery_eye_1 == 0 ~ od_corrected_1m,
          surgery_eye_1 == 1 ~ os_corrected_1m,
          surgery_eye_1 == 2 ~ (od_corrected_1m + os_corrected_1m)/2,
          TRUE ~ NA_real_
        ),
        vision_improvement_1w = post_vision_1w - pre_vision,
        vision_improvement_1m = post_vision_1m - pre_vision
      ) %>%
      dplyr::select(ID, vision_improvement_1w, vision_improvement_1m)
    
    return(vision_data)
  }
  
  va_improvements <- process_vision_improvements(baseline_info)
  cat("✓ VA data processed\n")
  
} else {
  cat("✓ 使用现有的OCTA和VA数据\n")
}

# ================== 3. 整合数据 ==================

# 整合所有数据
enhanced_cluster_analysis <- max_membership_data %>%
  left_join(octa_improvements, by = c("subject_id" = "ID")) %>%
  left_join(va_improvements, by = c("subject_id" = "ID"))

cat("✓ 数据整合完成\n")
cat("最终分析数据:", nrow(enhanced_cluster_analysis), "行,", ncol(enhanced_cluster_analysis), "列\n\n")

# ================== 4. 硬聚类组间差异分析函数 ==================

analyze_hard_cluster_differences <- function(data, membership_cols, cluster_cols, outcome_params, outcome_type) {
  
  cat(sprintf("=== 分析硬聚类组间差异 (%s) ===\n", outcome_type))
  
  results <- data.frame()
  
  for(i in 1:length(membership_cols)) {
    membership_col <- membership_cols[i]
    cluster_col <- cluster_cols[i]
    window_name <- gsub("^membership_", "", membership_col)
    
    cat(sprintf("分析 %s 窗口...\n", window_name))
    
    for(param in outcome_params) {
      if(param %in% names(data)) {
        # 清理数据
        clean_data <- data %>%
          filter(!is.na(!!sym(membership_col)) & !is.na(!!sym(param)) & !is.na(!!sym(cluster_col)))
        
        if(nrow(clean_data) >= 3) {
          
          # 聚类描述性统计
          cluster_summary <- clean_data %>%
            group_by(!!sym(cluster_col)) %>%
            summarise(
              n = n(),
              mean_outcome = mean(!!sym(param)),
              sd_outcome = sd(!!sym(param)),
              median_outcome = median(!!sym(param)),
              min_outcome = min(!!sym(param)),
              max_outcome = max(!!sym(param)),
              .groups = 'drop'
            )
          
          # 如果有2个cluster，进行t检验
          if(nrow(cluster_summary) == 2) {
            cluster1_data <- clean_data %>% filter(!!sym(cluster_col) == 1) %>% pull(!!sym(param))
            cluster2_data <- clean_data %>% filter(!!sym(cluster_col) == 2) %>% pull(!!sym(param))
            
            if(length(cluster1_data) >= 2 && length(cluster2_data) >= 2) {
              
              # t检验
              t_test <- try(t.test(cluster1_data, cluster2_data), silent = TRUE)
              
              # Welch t检验（方差不等）
              welch_test <- try(t.test(cluster1_data, cluster2_data, var.equal = FALSE), silent = TRUE)
              
              # Mann-Whitney U检验（非参数）
              mann_whitney <- try(wilcox.test(cluster1_data, cluster2_data), silent = TRUE)
              
              if(class(t_test) != "try-error" && class(welch_test) != "try-error" && class(mann_whitney) != "try-error") {
                
                # 计算效应大小 (Cohen's d)
                pooled_sd <- sqrt(((length(cluster1_data)-1)*var(cluster1_data) + 
                                     (length(cluster2_data)-1)*var(cluster2_data)) / 
                                    (length(cluster1_data) + length(cluster2_data) - 2))
                cohens_d <- abs(mean(cluster1_data) - mean(cluster2_data)) / pooled_sd
                
                # Glass's delta (使用对照组的标准差)
                glass_delta1 <- abs(mean(cluster1_data) - mean(cluster2_data)) / sd(cluster1_data)
                glass_delta2 <- abs(mean(cluster1_data) - mean(cluster2_data)) / sd(cluster2_data)
                
                # Hedges' g (小样本校正)
                hedges_g <- cohens_d * (1 - (3 / (4 * (length(cluster1_data) + length(cluster2_data)) - 9)))
                
                # 确定较好的cluster（outcome更高的）
                better_cluster <- ifelse(mean(cluster1_data) > mean(cluster2_data), 1, 2)
                better_cluster_mean <- max(mean(cluster1_data), mean(cluster2_data))
                worse_cluster_mean <- min(mean(cluster1_data), mean(cluster2_data))
                better_cluster_n <- ifelse(better_cluster == 1, length(cluster1_data), length(cluster2_data))
                worse_cluster_n <- ifelse(better_cluster == 1, length(cluster2_data), length(cluster1_data))
                
                # 效应大小分类
                effect_size_category <- case_when(
                  cohens_d >= 0.8 ~ "Large",
                  cohens_d >= 0.5 ~ "Medium", 
                  cohens_d >= 0.2 ~ "Small",
                  TRUE ~ "Negligible"
                )
                
                # 计算百分比改善
                percent_improvement <- ifelse(worse_cluster_mean != 0, 
                                              (better_cluster_mean - worse_cluster_mean) / abs(worse_cluster_mean) * 100,
                                              NA)
                
                result_row <- data.frame(
                  Time_Window = window_name,
                  Outcome_Type = outcome_type,
                  Outcome_Parameter = param,
                  N_Total = nrow(clean_data),
                  Cluster1_N = length(cluster1_data),
                  Cluster2_N = length(cluster2_data),
                  Cluster1_Mean = mean(cluster1_data),
                  Cluster1_SD = sd(cluster1_data),
                  Cluster1_Median = median(cluster1_data),
                  Cluster2_Mean = mean(cluster2_data),
                  Cluster2_SD = sd(cluster2_data),
                  Cluster2_Median = median(cluster2_data),
                  Mean_Difference = abs(mean(cluster1_data) - mean(cluster2_data)),
                  # t检验结果
                  T_Statistic = t_test$statistic,
                  T_P_Value = t_test$p.value,
                  T_CI_Lower = t_test$conf.int[1],
                  T_CI_Upper = t_test$conf.int[2],
                  # Welch t检验
                  Welch_T_Statistic = welch_test$statistic,
                  Welch_P_Value = welch_test$p.value,
                  # Mann-Whitney U检验
                  MW_U_Statistic = mann_whitney$statistic,
                  MW_P_Value = mann_whitney$p.value,
                  # 效应大小
                  Cohens_D = cohens_d,
                  Hedges_G = hedges_g,
                  Glass_Delta1 = glass_delta1,
                  Glass_Delta2 = glass_delta2,
                  Effect_Size_Category = effect_size_category,
                  # 临床解释
                  Better_Cluster = better_cluster,
                  Better_Cluster_N = better_cluster_n,
                  Better_Cluster_Mean = better_cluster_mean,
                  Worse_Cluster_Mean = worse_cluster_mean,
                  Percent_Improvement = percent_improvement,
                  # 显著性
                  T_Significant = t_test$p.value < 0.05,
                  Welch_Significant = welch_test$p.value < 0.05,
                  MW_Significant = mann_whitney$p.value < 0.05,
                  stringsAsFactors = FALSE
                )
                
                results <- rbind(results, result_row)
              }
            }
          } else if(nrow(cluster_summary) > 2) {
            # 如果有超过2个cluster，进行ANOVA
            cat(sprintf("  警告：%s窗口的%s有%d个clusters，跳过分析\n", 
                        window_name, param, nrow(cluster_summary)))
          }
        }
      }
    }
  }
  
  return(results)
}

# ================== 5. 执行硬聚类分析 ==================

# 获取outcome参数
octa_improvement_params <- names(enhanced_cluster_analysis)[grep("_improvement$", names(enhanced_cluster_analysis))]
octa_improvement_params <- octa_improvement_params[!grepl("vision_improvement", octa_improvement_params)]
va_improvement_params <- c("vision_improvement_1w", "vision_improvement_1m")

cat("🔬 开始硬聚类组间差异分析...\n")
cat("OCTA改善参数:", length(octa_improvement_params), "个\n")
cat("VA改善参数:", length(va_improvement_params), "个\n\n")

# 分析OCTA和VA的cluster差异
octa_cluster_results <- analyze_hard_cluster_differences(
  enhanced_cluster_analysis, membership_cols, cluster_cols, 
  octa_improvement_params, "OCTA"
)

va_cluster_results <- analyze_hard_cluster_differences(
  enhanced_cluster_analysis, membership_cols, cluster_cols,
  va_improvement_params, "VA"
)

# 合并结果
all_cluster_results <- rbind(octa_cluster_results, va_cluster_results)

# ================== 6. 精细化FDR校正 ==================

perform_refined_fdr_correction <- function(all_cluster_results) {
  
  cat("=== 执行精细化的FDR校正 ===\n")
  
  if(nrow(all_cluster_results) == 0) {
    cat("❌ 没有结果进行校正\n")
    return(data.frame())
  }
  
  # 更精细的参数分类
  refined_results <- all_cluster_results %>%
    mutate(
      Parameter_Type_Detailed = case_when(
        # 血流参数 - 按层次分类
        grepl("VD_SVP|PA_SVP", Outcome_Parameter) ~ "BloodFlow_SVP",
        grepl("VD_ICP|PA_ICP", Outcome_Parameter) ~ "BloodFlow_ICP", 
        grepl("VD_DCP|PA_DCP", Outcome_Parameter) ~ "BloodFlow_DCP",
        grepl("VD_Choroid|PA_Choroid", Outcome_Parameter) ~ "BloodFlow_Choroid",
        # 厚度参数
        grepl("Thickness_GCL", Outcome_Parameter) ~ "Thickness_GCL_IPL",
        grepl("Thickness_INL", Outcome_Parameter) ~ "Thickness_INL",
        grepl("Thickness_Retina", Outcome_Parameter) ~ "Thickness_Retina",
        # 视力参数
        grepl("vision_improvement", Outcome_Parameter) ~ "Vision",
        TRUE ~ "Other"
      ),
      # 区域分类
      Region_Type = case_when(
        grepl("_0_6_", Outcome_Parameter) ~ "Macular",
        grepl("_0_21_", Outcome_Parameter) ~ "Widefield", 
        grepl("vision_improvement", Outcome_Parameter) ~ "Vision",
        TRUE ~ "Other"
      )
    )
  
  # 显示分类结果
  cat("\n📋 参数分类统计:\n")
  classification_summary <- refined_results %>%
    count(Time_Window, Parameter_Type_Detailed, Region_Type, name = "Tests") %>%
    arrange(Time_Window, Parameter_Type_Detailed, Region_Type)
  
  print(classification_summary)
  
  # 策略1: 按时间窗口 + 详细参数类型分别校正
  cat("\n🔬 策略1: 按时间窗口 + 详细参数类型分别校正\n")
  
  strategy1_results <- refined_results %>%
    group_by(Time_Window, Parameter_Type_Detailed) %>%
    mutate(
      T_P_Value_FDR_Strategy1 = p.adjust(T_P_Value, method = "fdr"),
      Welch_P_Value_FDR_Strategy1 = p.adjust(Welch_P_Value, method = "fdr"),
      MW_P_Value_FDR_Strategy1 = p.adjust(MW_P_Value, method = "fdr"),
      T_Significant_FDR_Strategy1 = T_P_Value_FDR_Strategy1 < 0.05,
      Welch_Significant_FDR_Strategy1 = Welch_P_Value_FDR_Strategy1 < 0.05,
      MW_Significant_FDR_Strategy1 = MW_P_Value_FDR_Strategy1 < 0.05,
      Tests_in_Group_Strategy1 = n()
    ) %>%
    ungroup()
  
  # 策略2: 按时间窗口 + 参数类型 + 区域分别校正
  cat("🔬 策略2: 按时间窗口 + 参数类型 + 区域分别校正\n")
  
  strategy2_results <- strategy1_results %>%
    group_by(Time_Window, Parameter_Type_Detailed, Region_Type) %>%
    mutate(
      T_P_Value_FDR_Strategy2 = p.adjust(T_P_Value, method = "fdr"),
      Welch_P_Value_FDR_Strategy2 = p.adjust(Welch_P_Value, method = "fdr"),
      MW_P_Value_FDR_Strategy2 = p.adjust(MW_P_Value, method = "fdr"),
      T_Significant_FDR_Strategy2 = T_P_Value_FDR_Strategy2 < 0.05,
      Welch_Significant_FDR_Strategy2 = Welch_P_Value_FDR_Strategy2 < 0.05,
      MW_Significant_FDR_Strategy2 = MW_P_Value_FDR_Strategy2 < 0.05,
      Tests_in_Group_Strategy2 = n()
    ) %>%
    ungroup()
  
  # 策略3: 仅按参数类型校正（不考虑时间窗口）
  cat("🔬 策略3: 仅按详细参数类型校正\n")
  
  final_results <- strategy2_results %>%
    group_by(Parameter_Type_Detailed) %>%
    mutate(
      T_P_Value_FDR_Strategy3 = p.adjust(T_P_Value, method = "fdr"),
      Welch_P_Value_FDR_Strategy3 = p.adjust(Welch_P_Value, method = "fdr"),
      MW_P_Value_FDR_Strategy3 = p.adjust(MW_P_Value, method = "fdr"),
      T_Significant_FDR_Strategy3 = T_P_Value_FDR_Strategy3 < 0.05,
      Welch_Significant_FDR_Strategy3 = Welch_P_Value_FDR_Strategy3 < 0.05,
      MW_Significant_FDR_Strategy3 = MW_P_Value_FDR_Strategy3 < 0.05,
      Tests_in_Group_Strategy3 = n()
    ) %>%
    ungroup()
  
  # 比较不同策略的结果
  cat("\n📊 不同FDR校正策略的结果比较:\n")
  
  strategy_summary <- data.frame(
    Strategy = c(
      "原始显著(t-test)", "策略1(t-test)", "策略2(t-test)", "策略3(t-test)",
      "原始显著(Welch)", "策略1(Welch)", "策略2(Welch)", "策略3(Welch)",
      "原始显著(Mann-Whitney)", "策略1(Mann-Whitney)", "策略2(Mann-Whitney)", "策略3(Mann-Whitney)"
    ),
    Significant_Count = c(
      sum(final_results$T_Significant),
      sum(final_results$T_Significant_FDR_Strategy1),
      sum(final_results$T_Significant_FDR_Strategy2),
      sum(final_results$T_Significant_FDR_Strategy3),
      sum(final_results$Welch_Significant),
      sum(final_results$Welch_Significant_FDR_Strategy1),
      sum(final_results$Welch_Significant_FDR_Strategy2),
      sum(final_results$Welch_Significant_FDR_Strategy3),
      sum(final_results$MW_Significant),
      sum(final_results$MW_Significant_FDR_Strategy1),
      sum(final_results$MW_Significant_FDR_Strategy2),
      sum(final_results$MW_Significant_FDR_Strategy3)
    ),
    Total_Tests = rep(nrow(final_results), 12)
  ) %>%
    mutate(
      Success_Rate = round(Significant_Count / Total_Tests * 100, 1)
    )
  
  print(strategy_summary)
  
  # 推荐策略
  cat("\n💡 推荐使用策略2（按时间窗口 + 参数类型 + 区域分别校正）:\n")
  
  # 选择最严格的Welch t检验作为主要结果
  recommended_significant <- final_results %>%
    filter(Welch_Significant_FDR_Strategy2) %>%
    arrange(Welch_P_Value_FDR_Strategy2)
  
  if(nrow(recommended_significant) > 0) {
    cat(sprintf("\n🎯 推荐策略下的显著结果 (%d个, Welch t-test):\n", nrow(recommended_significant)))
    
    for(i in 1:nrow(recommended_significant)) {
      result <- recommended_significant[i, ]
      cat(sprintf("\n%d. %s - %s Window:\n", i, result$Outcome_Parameter, result$Time_Window))
      cat(sprintf("   📊 Cluster差异: %.3f vs %.3f (差值=%.3f)\n", 
                  result$Cluster1_Mean, result$Cluster2_Mean, result$Mean_Difference))
      cat(sprintf("   📈 Welch t检验: t=%.3f, p=%.4f → p_FDR=%.4f\n", 
                  result$Welch_T_Statistic, result$Welch_P_Value, result$Welch_P_Value_FDR_Strategy2))
      cat(sprintf("   🔬 效应大小: Cohen's d=%.3f (%s), Hedges' g=%.3f\n", 
                  result$Cohens_D, result$Effect_Size_Category, result$Hedges_G))
      cat(sprintf("   🎯 更好cluster: Cluster %s (n=%d, mean=%.3f)\n", 
                  result$Better_Cluster, result$Better_Cluster_N, result$Better_Cluster_Mean))
      cat(sprintf("   📋 校正组: %s_%s_%s (%d个测试)\n", 
                  result$Time_Window, result$Parameter_Type_Detailed, result$Region_Type,
                  result$Tests_in_Group_Strategy2))
      if(!is.na(result$Percent_Improvement)) {
        cat(sprintf("   📈 相对改善: %.1f%%\n", result$Percent_Improvement))
      }
    }
  } else {
    cat("   在推荐策略下仍无显著结果\n")
    
    # 显示最接近显著的结果
    near_significant <- final_results %>%
      filter(T_Significant) %>%
      arrange(Welch_P_Value_FDR_Strategy2) %>%
      slice_head(n = 3)
    
    cat("\n🔍 最接近显著的结果 (Top 3, Welch t-test):\n")
    for(i in 1:nrow(near_significant)) {
      result <- near_significant[i, ]
      cat(sprintf("%d. %s - %s: p_FDR=%.4f (组内%d个测试, Cohen's d=%.3f)\n", 
                  i, result$Outcome_Parameter, result$Time_Window,
                  result$Welch_P_Value_FDR_Strategy2, result$Tests_in_Group_Strategy2, result$Cohens_D))
    }
  }
  
  return(final_results)
}

# 执行精细化校正
if(nrow(all_cluster_results) > 0) {
  refined_results <- perform_refined_fdr_correction(all_cluster_results)
} else {
  cat("❌ 没有cluster差异结果进行FDR校正\n")
  refined_results <- data.frame()
}

# ================== 7. 创建输出目录并保存结果 ==================

# 创建输出目录
output_dir <- "3_data_analysis/6_clustering_modeling/time_window_cluster_compare_analysis"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
setwd(output_dir)

# 保存结果
if(nrow(refined_results) > 0) {
  write.csv(refined_results, "time_window_hard_cluster_analysis_complete.csv", row.names = FALSE)
  
  # 提取推荐策略的显著结果
  recommended_significant_results <- refined_results %>%
    filter(Welch_Significant_FDR_Strategy2) %>%
    arrange(Welch_P_Value_FDR_Strategy2)
  
  write.csv(recommended_significant_results, "time_window_hard_cluster_analysis_refined_significant.csv", row.names = FALSE)
  
  # 提取原始显著结果
  original_significant_results <- refined_results %>%
    filter(T_Significant) %>%
    arrange(T_P_Value)
  
  write.csv(original_significant_results, "time_window_hard_cluster_analysis_original_significant.csv", row.names = FALSE)
  
  cat("✅ 精细化硬聚类分析完成！\n")
  cat("📁 结果已保存：\n")
  cat("   - time_window_hard_cluster_analysis_complete.csv (完整结果)\n")
  cat("   - time_window_hard_cluster_analysis_refined_significant.csv (FDR校正后显著)\n")
  cat("   - time_window_hard_cluster_analysis_original_significant.csv (原始显著)\n")
} else {
  cat("❌ 没有硬聚类分析结果保存\n")
}

# ================== 8. 结果输出和分析 ==================

if(nrow(all_cluster_results) > 0) {
  cat("\n📊 硬聚类分析结果摘要:\n")
  cat(sprintf("- 总测试数: %d\n", nrow(all_cluster_results)))
  cat(sprintf("- t检验显著 (p < 0.05): %d (%.1f%%)\n", 
              sum(all_cluster_results$T_Significant), 
              sum(all_cluster_results$T_Significant)/nrow(all_cluster_results)*100))
  cat(sprintf("- Welch t检验显著: %d (%.1f%%)\n", 
              sum(all_cluster_results$Welch_Significant),
              sum(all_cluster_results$Welch_Significant)/nrow(all_cluster_results)*100))
  cat(sprintf("- Mann-Whitney显著: %d (%.1f%%)\n", 
              sum(all_cluster_results$MW_Significant),
              sum(all_cluster_results$MW_Significant)/nrow(all_cluster_results)*100))
  
  if(nrow(refined_results) > 0) {
    cat(sprintf("- FDR校正后显著(策略2-Welch): %d\n", sum(refined_results$Welch_Significant_FDR_Strategy2)))
  }
  
  # 显示原始显著结果
  if(nrow(original_significant_results) > 0) {
    cat("\n🎯 Top 5 显著的cluster差异 (t-test):\n")
    for(i in 1:min(5, nrow(original_significant_results))) {
      result <- original_significant_results[i, ]
      cat(sprintf("\n%d. %s (%s) - %s Window:\n", i, result$Outcome_Parameter, result$Outcome_Type, result$Time_Window))
      cat(sprintf("   📊 Cluster差异: %.3f vs %.3f (差值=%.3f)\n", 
                  result$Cluster1_Mean, result$Cluster2_Mean, result$Mean_Difference))
      cat(sprintf("   📈 统计: t=%.3f, p=%.4f | Welch t=%.3f, p=%.4f\n", 
                  result$T_Statistic, result$T_P_Value, result$Welch_T_Statistic, result$Welch_P_Value))
      cat(sprintf("   📈 非参数: MW U=%.0f, p=%.4f\n", result$MW_U_Statistic, result$MW_P_Value))
      cat(sprintf("   🔬 效应: Cohen's d=%.3f (%s), Hedges' g=%.3f\n", 
                  result$Cohens_D, result$Effect_Size_Category, result$Hedges_G))
      cat(sprintf("   🎯 较好cluster: Cluster %s (n=%d, mean=%.3f)\n", 
                  result$Better_Cluster, result$Better_Cluster_N, result$Better_Cluster_Mean))
      if(!is.na(result$Percent_Improvement)) {
        cat(sprintf("   📈 相对改善: %.1f%%\n", result$Percent_Improvement))
      }
    }
    
    # 按时间窗口汇总
    cat("\n📈 按时间窗口汇总 (t-test显著):\n")
    window_summary <- original_significant_results %>%
      group_by(Time_Window) %>%
      summarise(
        Count = n(),
        Mean_Effect_Size = round(mean(Cohens_D), 3),
        Max_Effect_Size = round(max(Cohens_D), 3),
        Best_P = min(T_P_Value),
        Large_Effects = sum(Effect_Size_Category == "Large"),
        .groups = 'drop'
      ) %>%
      arrange(desc(Count))
    
    print(window_summary)
    
    # 按参数类型汇总
    cat("\n🔬 按参数类型汇总:\n")
    param_summary <- original_significant_results %>%
      mutate(
        Param_Category = case_when(
          grepl("BloodFlow", Parameter_Type_Detailed) ~ "Blood Flow",
          grepl("Thickness", Parameter_Type_Detailed) ~ "Thickness",
          grepl("Vision", Parameter_Type_Detailed) ~ "Vision",
          TRUE ~ "Other"
        )
      ) %>%
      group_by(Param_Category) %>%
      summarise(
        Count = n(),
        Mean_Effect_Size = round(mean(Cohens_D), 3),
        Large_Effects = sum(Effect_Size_Category == "Large"),
        .groups = 'drop'
      ) %>%
      arrange(desc(Count))
    
    print(param_summary)
  }
}

# 补充缺少的FDR校正可视化函数

# ================== 在你的代码第9节可视化部分替换为以下内容 ==================

# ================== 9. 创建FDR校正后的可视化函数 ==================

create_fdr_corrected_hard_cluster_visualizations <- function(data, analysis_results, fdr_threshold = 0.05) {
  
  cat("\n===== 创建FDR校正后的硬聚类分析可视化 =====\n")
  
  if(nrow(analysis_results) == 0) {
    cat("❌ 没有分析结果可用于可视化\n")
    return(NULL)
  }
  
  # 🎯 关键修改：优先使用FDR校正后显著的结果
  if("Welch_Significant_FDR_Strategy2" %in% names(analysis_results)) {
    # 使用FDR校正后的结果
    significant_for_viz <- analysis_results %>%
      filter(Welch_Significant_FDR_Strategy2) %>%  # FDR校正后显著
      arrange(desc(Cohens_D)) %>%
      slice_head(n = 12)  # 增加显示数量，因为FDR后可能变少
    
    cat(sprintf("📊 使用FDR校正后的显著结果创建可视化 (%d个变量)\n", nrow(significant_for_viz)))
    
    # 如果FDR校正后没有显著结果，回退到原始显著结果但给出警告
    if(nrow(significant_for_viz) == 0) {
      cat("⚠️ FDR校正后无显著结果，回退到原始显著结果\n")
      significant_for_viz <- analysis_results %>%
        filter(T_Significant) %>%
        arrange(desc(Cohens_D)) %>%
        slice_head(n = 6)
      use_fdr_correction <- FALSE
    } else {
      use_fdr_correction <- TRUE
    }
    
  } else {
    # 如果没有FDR校正结果，使用原始显著结果
    cat("⚠️ 没有FDR校正结果，使用原始显著结果\n")
    significant_for_viz <- analysis_results %>%
      filter(T_Significant) %>%
      arrange(desc(Cohens_D)) %>%
      slice_head(n = 6)
    use_fdr_correction <- FALSE
  }
  
  if(nrow(significant_for_viz) == 0) {
    cat("❌ 没有显著结果可用于可视化\n")
    return(NULL)
  }
  
  # 创建箱线图列表
  plot_list <- list()
  
  for(i in 1:nrow(significant_for_viz)) {
    result <- significant_for_viz[i, ]
    param <- result$Outcome_Parameter
    window <- result$Time_Window
    cluster_col <- paste0("cluster_", window)
    
    if(cluster_col %in% names(data) && param %in% names(data)) {
      # 准备绘图数据
      plot_data <- data %>%
        filter(!is.na(!!sym(cluster_col)) & !is.na(!!sym(param))) %>%
        mutate(
          cluster_label = paste0("Cluster ", !!sym(cluster_col)),
          cluster_factor = factor(!!sym(cluster_col))
        )
      
      if(nrow(plot_data) >= 3) {
        
        # 创建参数的清晰名称
        param_clean <- param %>%
          gsub("_improvement", " Improvement", .) %>%
          gsub("_", " ", .) %>%
          gsub("vision improvement", "Vision Improvement", .) %>%
          gsub("1w", "(1 Week)", .) %>%
          gsub("1m", "(1 Month)", .) %>%
          gsub("SVP", "Superficial Vascular Plexus", .) %>%
          gsub("ICP", "Intermediate Capillary Plexus", .) %>%
          gsub("DCP", "Deep Capillary Plexus", .) %>%
          gsub("GCL IPL", "GCL-IPL", .) %>%
          gsub("0 6", "Macular (0_6)", .) %>%
          gsub("0 21", "Widefield (0_21)", .)
        
        window_clean <- gsub("_", " ", window) %>% str_to_title()
        
        # 🎯 关键修改：标题和统计信息显示FDR校正状态
        if(use_fdr_correction) {
          # 使用FDR校正后的p值
          p_value_to_show <- result$Welch_P_Value_FDR_Strategy2
          p_label <- "p_FDR"
          subtitle_text <- paste0("FDR-Corrected Cluster Comparison | ",
                                  "Welch t = ", round(result$Welch_T_Statistic, 3), 
                                  ", ", p_label, " = ", format.pval(p_value_to_show, digits = 3),
                                  " | Cohen's d = ", round(result$Cohens_D, 3),
                                  " (", result$Effect_Size_Category, " effect)",
                                  " | n = ", result$N_Total)
        } else {
          # 使用原始p值
          p_value_to_show <- result$T_P_Value
          p_label <- "p"
          subtitle_text <- paste0("Original Cluster Comparison | ",
                                  "t = ", round(result$T_Statistic, 3), 
                                  ", ", p_label, " = ", format.pval(p_value_to_show, digits = 3),
                                  " | Cohen's d = ", round(result$Cohens_D, 3),
                                  " (", result$Effect_Size_Category, " effect)",
                                  " | n = ", result$N_Total)
        }
        
        # 创建箱线图 + 散点图
        p <- ggplot(plot_data, aes(x = cluster_factor, y = !!sym(param), fill = cluster_factor)) +
          geom_boxplot(alpha = 0.7, outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
          stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                       fill = "red", color = "darkred") +
          scale_fill_manual(
            values = c("1" = "#a488bf", "2" = "#bd992e"),
            labels = c("1" = "Cluster 1", "2" = "Cluster 2"),
            name = "Cluster"
          ) +
          labs(
            title = paste(param_clean, "-", window_clean, "Window"),
            subtitle = subtitle_text,
            x = "Cluster",
            y = "Improvement Value",
            caption = paste("Better cluster:", result$Better_Cluster, 
                            "| Red diamonds = means | Improvement:", 
                            ifelse(!is.na(result$Percent_Improvement), 
                                   paste0(round(result$Percent_Improvement, 1), "%"), "N/A"),
                            ifelse(use_fdr_correction, "| FDR corrected", "| Original p-values"))
          ) +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            legend.position = "bottom",
            plot.caption = element_text(size = 8)
          )
        
        plot_list[[i]] <- p
        
        # 保存单独的图形
        filename_suffix <- ifelse(use_fdr_correction, "_FDR_corrected", "_original")
        ggsave(paste0("hard_cluster_", gsub("[^A-Za-z0-9]", "_", param), "_", window, filename_suffix, "_boxplot.pdf"),
               p, width = 8, height = 6)
      }
    }
  }
  
  # 创建组合图形
  if(length(plot_list) > 0) {
    # 计算布局
    ncol_layout <- ifelse(length(plot_list) >= 4, 2, 1)
    nrow_layout <- ceiling(length(plot_list) / ncol_layout)
    
    # 组合图形标题
    combined_title <- ifelse(use_fdr_correction, 
                             "Time Window Hard Cluster Analysis - FDR Corrected Significant Differences",
                             "Time Window Hard Cluster Analysis - Original Significant Differences")
    
    # 组合图形
    combined_plot <- do.call(gridExtra::grid.arrange, 
                             c(plot_list, 
                               ncol = ncol_layout,
                               top = combined_title))
    
    # 保存组合图形
    filename_suffix <- ifelse(use_fdr_correction, "_FDR_corrected", "_original")
    ggsave(paste0("time_window_hard_cluster_analysis_combined", filename_suffix, ".pdf"),
           combined_plot, width = 14, height = 6 * nrow_layout)
    
    cat("✓ 可视化已保存\n")
    
    # 🎯 创建FDR校正状态报告
    if(use_fdr_correction) {
      cat("\n🎯 FDR校正状态报告:\n")
      cat("✅ 使用FDR校正后的p值 (Welch t-test, Strategy 2)\n")
      cat(sprintf("✅ 显示了 %d 个FDR校正后显著的变量\n", length(plot_list)))
      cat("✅ 所有p值都经过了多重比较校正\n")
    } else {
      cat("\n⚠️ 使用原始p值 (未经FDR校正)\n")
      cat(sprintf("⚠️ 显示了 %d 个原始显著的变量\n", length(plot_list)))
      cat("⚠️ 建议检查FDR校正结果\n")
    }
    
    return(combined_plot)
  }
  
  return(NULL)
}

# ================== 创建FDR校正对比图函数 ==================

create_fdr_correction_comparison <- function(original_results, fdr_results) {
  
  cat("\n📊 创建FDR校正前后对比图...\n")
  
  if(nrow(original_results) == 0 || nrow(fdr_results) == 0) {
    cat("❌ 缺少结果数据，无法创建对比图\n")
    return(NULL)
  }
  
  # 合并原始和FDR结果
  comparison_data <- original_results %>%
    dplyr::select(Time_Window, Outcome_Parameter, T_P_Value, Welch_P_Value, Cohens_D, Effect_Size_Category) %>%
    left_join(
      fdr_results %>% 
        dplyr::select(Time_Window, Outcome_Parameter, Welch_P_Value_FDR_Strategy2, Welch_Significant_FDR_Strategy2),
      by = c("Time_Window", "Outcome_Parameter")
    ) %>%
    filter(!is.na(Welch_P_Value_FDR_Strategy2)) %>%  # 过滤掉没有FDR结果的行
    mutate(
      Original_Significant = T_P_Value < 0.05,
      FDR_Significant = Welch_Significant_FDR_Strategy2,
      Significance_Status = case_when(
        Original_Significant & FDR_Significant ~ "Both Significant",
        Original_Significant & !FDR_Significant ~ "Only Original",
        !Original_Significant & FDR_Significant ~ "Only FDR",
        TRUE ~ "Not Significant"
      )
    )
  
  if(nrow(comparison_data) == 0) {
    cat("❌ 没有可比较的数据\n")
    return(NULL)
  }
  
  # 创建对比图
  p_comparison <- ggplot(comparison_data, aes(x = T_P_Value, y = Welch_P_Value_FDR_Strategy2)) +
    geom_point(aes(color = Significance_Status, size = Cohens_D), alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
    scale_color_manual(
      values = c(
        "Both Significant" = "#2E8B57",
        "Only Original" = "#FF6347", 
        "Only FDR" = "#4169E1",
        "Not Significant" = "#D3D3D3"
      ),
      name = "Significance Status"
    ) +
    scale_size_continuous(name = "Effect Size", range = c(1, 4)) +
    labs(
      title = "FDR Correction Impact on Significance",
      subtitle = "Comparison of original vs FDR-corrected p-values",
      x = "Original p-value",
      y = "FDR-corrected p-value",
      caption = "Dashed lines: significance threshold (p = 0.05)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom"
    )
  
  ggsave("fdr_correction_comparison.pdf", p_comparison, width = 10, height = 8)
  
  # 打印统计摘要
  cat("\n📋 FDR校正影响统计:\n")
  significance_summary <- comparison_data %>%
    count(Significance_Status) %>%
    mutate(Percentage = round(n / sum(n) * 100, 1))
  
  print(significance_summary)
  
  return(p_comparison)
}

# ================== 修正后的可视化调用部分 ==================

# 替换你原来的可视化创建代码
if(nrow(all_cluster_results) > 0) {
  
  cat("\n📊 创建FDR校正后的硬聚类分析可视化...\n")
  
  # 如果有FDR校正结果，使用FDR校正后的结果
  if(exists("refined_results") && nrow(refined_results) > 0) {
    fdr_corrected_plots <- create_fdr_corrected_hard_cluster_visualizations(
      data = enhanced_cluster_analysis,
      analysis_results = refined_results,  # 使用FDR校正后的结果
      fdr_threshold = 0.05
    )
    
    # 创建FDR校正前后对比图
    if(nrow(refined_results) > 0) {
      fdr_comparison_plot <- create_fdr_correction_comparison(all_cluster_results, refined_results)
    }
  } else {
    # 如果没有FDR校正结果，使用原始结果
    cat("⚠️ 没有FDR校正结果，使用原始结果创建可视化\n")
    fdr_corrected_plots <- create_fdr_corrected_hard_cluster_visualizations(
      data = enhanced_cluster_analysis,
      analysis_results = all_cluster_results,  # 使用原始结果
      fdr_threshold = 0.05
    )
  }
  
} else {
  cat("❌ 没有硬聚类分析结果用于可视化\n")
}

# ================== 使用说明 ==================

cat("\n💡 修正说明:\n")
cat("1. ✅ 补充了缺少的 create_fdr_corrected_hard_cluster_visualizations 函数\n")
cat("2. ✅ 补充了缺少的 create_fdr_correction_comparison 函数\n")
cat("3. ✅ 修正了函数调用部分，添加了错误处理\n")
cat("4. ✅ 图表会自动检测是否有FDR校正结果并相应标注\n")
cat("5. ✅ 生成的图表文件名会包含 '_FDR_corrected' 或 '_original' 后缀\n\n")

cat("🎨 现在你的代码应该能正常运行，生成:\n")
cat("- FDR校正后的箱线图（如果有显著结果）\n")
cat("- FDR校正前后的对比图\n")
cat("- 明确标注是否使用了FDR校正的图表\n")

# ================== 10. 效应大小分析 ==================

analyze_effect_sizes <- function(analysis_results) {
  
  if(nrow(analysis_results) == 0) {
    cat("❌ 没有结果进行效应大小分析\n")
    return(NULL)
  }
  
  cat("\n===== 效应大小深度分析 =====\n")
  
  # 效应大小分布
  effect_size_dist <- analysis_results %>%
    count(Effect_Size_Category) %>%
    mutate(Percentage = round(n / sum(n) * 100, 1))
  
  cat("📊 效应大小分布:\n")
  print(effect_size_dist)
  
  # 大效应结果
  large_effects <- analysis_results %>%
    filter(Effect_Size_Category == "Large") %>%
    arrange(desc(Cohens_D)) %>%
    dplyr::select(Time_Window, Outcome_Parameter, Cohens_D, T_P_Value, Better_Cluster, Percent_Improvement)
  
  if(nrow(large_effects) > 0) {
    cat(sprintf("\n🎯 大效应结果 (%d个):\n", nrow(large_effects)))
    print(large_effects)
  }
  
  # 创建效应大小可视化
  p_effect <- ggplot(analysis_results, aes(x = Cohens_D, y = -log10(T_P_Value))) +
    geom_point(aes(color = Effect_Size_Category, size = abs(Percent_Improvement)), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(0.2, 0.5, 0.8), linetype = "dashed", color = "gray") +
    scale_color_manual(
      values = c("Large" = "#D73027", "Medium" = "#FC8D59", "Small" = "#91BFDB", "Negligible" = "#EEEEEE"),
      name = "Effect Size"
    ) +
    scale_size_continuous(name = "% Improvement", range = c(1, 5)) +
    labs(
      title = "Effect Size vs Statistical Significance",
      subtitle = "Volcano plot of cluster differences",
      x = "Cohen's d (Effect Size)",
      y = "-log10(p-value)",
      caption = "Dashed red line: p = 0.05 | Dashed gray lines: Cohen's d thresholds (0.2, 0.5, 0.8)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  
  ggsave("effect_size_volcano_plot.pdf", p_effect, width = 10, height = 8)
  
  return(list(
    distribution = effect_size_dist,
    large_effects = large_effects,
    volcano_plot = p_effect
  ))
}

# 执行效应大小分析
if(nrow(all_cluster_results) > 0) {
  effect_analysis <- analyze_effect_sizes(all_cluster_results)
}

# ================== 11. 生成最终报告 ==================

generate_hard_cluster_analysis_report <- function(analysis_results, refined_results = NULL) {
  
  report <- paste0(
    "========================================\n",
    "TIME WINDOW HARD CLUSTER ANALYSIS REPORT\n",
    "========================================\n\n",
    
    "🎯 ANALYSIS OVERVIEW:\n",
    "✅ 针对硬聚类（所有membership=1）进行组间差异分析\n",
    "✅ 使用t检验、Welch t检验和Mann-Whitney U检验\n",
    "✅ 计算多种效应大小指标（Cohen's d, Hedges' g等）\n",
    "✅ 实施精细化FDR校正策略\n\n",
    
    "📊 ANALYSIS DETAILS:\n",
    "- Analysis Date: ", Sys.Date(), "\n",
    "- Total Comparisons: ", nrow(analysis_results), "\n"
  )
  
  if(nrow(analysis_results) > 0) {
    report <- paste0(report,
                     "- t-test Significant (p < 0.05): ", sum(analysis_results$T_Significant), "\n",
                     "- Welch t-test Significant: ", sum(analysis_results$Welch_Significant), "\n",
                     "- Mann-Whitney Significant: ", sum(analysis_results$MW_Significant), "\n")
    
    if(!is.null(refined_results) && nrow(refined_results) > 0) {
      report <- paste0(report,
                       "- FDR Corrected Significant (Strategy 2): ", sum(refined_results$Welch_Significant_FDR_Strategy2), "\n")
    }
    
    # 效应大小分布
    effect_dist <- table(analysis_results$Effect_Size_Category)
    report <- paste0(report, "\n🔬 EFFECT SIZE DISTRIBUTION:\n")
    for(effect in names(effect_dist)) {
      report <- paste0(report, sprintf("- %s: %d (%.1f%%)\n", 
                                       effect, effect_dist[effect], 
                                       effect_dist[effect]/nrow(analysis_results)*100))
    }
    
    # Top发现
    significant_results <- analysis_results %>%
      filter(T_Significant) %>%
      arrange(desc(Cohens_D)) %>%
      slice_head(n = 3)
    
    if(nrow(significant_results) > 0) {
      report <- paste0(report, "\n🏆 TOP FINDINGS:\n")
      
      for(i in 1:nrow(significant_results)) {
        result <- significant_results[i, ]
        report <- paste0(report,
                         sprintf("\n%d. %s - %s Window:\n", i, result$Outcome_Parameter, result$Time_Window),
                         sprintf("   📊 Difference: %.3f vs %.3f\n", result$Cluster1_Mean, result$Cluster2_Mean),
                         sprintf("   📈 Statistics: t=%.3f, p=%.4f\n", result$T_Statistic, result$T_P_Value),
                         sprintf("   🔬 Effect: Cohen's d=%.3f (%s)\n", result$Cohens_D, result$Effect_Size_Category),
                         sprintf("   🎯 Better Cluster: %s (n=%d)\n", result$Better_Cluster, result$Better_Cluster_N))
        
        if(!is.na(result$Percent_Improvement)) {
          report <- paste0(report, sprintf("   📈 Improvement: %.1f%%\n", result$Percent_Improvement))
        }
      }
    }
    
    # 时间窗口表现
    window_performance <- analysis_results %>%
      filter(T_Significant) %>%
      group_by(Time_Window) %>%
      summarise(
        Count = n(),
        Mean_Effect = round(mean(Cohens_D), 3),
        Large_Effects = sum(Effect_Size_Category == "Large"),
        .groups = 'drop'
      ) %>%
      arrange(desc(Count))
    
    if(nrow(window_performance) > 0) {
      report <- paste0(report, "\n⏰ TIME WINDOW PERFORMANCE:\n")
      for(i in 1:nrow(window_performance)) {
        window_info <- window_performance[i, ]
        report <- paste0(report,
                         sprintf("%d. %s: %d significant differences (Mean d=%.3f, %d large effects)\n",
                                 i, window_info$Time_Window, window_info$Count,
                                 window_info$Mean_Effect, window_info$Large_Effects))
      }
    }
  }
  
  report <- paste0(report,
                   "\n📁 GENERATED FILES:\n",
                   "- time_window_hard_cluster_analysis_complete.csv: Complete results\n",
                   "- time_window_hard_cluster_analysis_refined_significant.csv: FDR significant results\n",
                   "- time_window_hard_cluster_analysis_original_significant.csv: Original significant results\n",
                   "- hard_cluster_*.pdf: Individual cluster comparison plots\n",
                   "- time_window_hard_cluster_analysis_combined.pdf: Combined visualization\n",
                   "- effect_size_volcano_plot.pdf: Effect size visualization\n\n",
                   
                   "🎯 INTERPRETATION GUIDELINES:\n",
                   "- Cohen's d ≥ 0.8: Large effect (clinically very meaningful)\n",
                   "- Cohen's d ≥ 0.5: Medium effect (clinically meaningful)\n",
                   "- Cohen's d ≥ 0.2: Small effect (may be clinically relevant)\n",
                   "- p < 0.05: Statistically significant (before correction)\n",
                   "- p_FDR < 0.05: Significant after multiple comparison correction\n\n",
                   
                   "🔬 STATISTICAL METHODS:\n",
                   "1. Student's t-test: Assumes equal variances\n",
                   "2. Welch's t-test: Does not assume equal variances (recommended)\n",
                   "3. Mann-Whitney U: Non-parametric alternative\n",
                   "4. Cohen's d: Standardized effect size\n",
                   "5. Hedges' g: Small sample corrected effect size\n",
                   "6. FDR correction: Controls false discovery rate\n\n",
                   
                   "💡 CLINICAL IMPLICATIONS:\n",
                   "1. Focus on large effect size differences (d ≥ 0.8)\n",
                   "2. Consider both statistical and clinical significance\n",
                   "3. Better clusters indicate more favorable recovery patterns\n",
                   "4. Time window differences suggest optimal intervention periods\n",
                   "5. Validate findings in independent cohorts\n\n",
                   
                   "🚀 NEXT STEPS:\n",
                   "1. Investigate baseline characteristics of better-performing clusters\n",
                   "2. Develop cluster-specific treatment protocols\n",
                   "3. Validate cluster assignments in prospective studies\n",
                   "4. Consider time-window-specific interventions\n",
                   "5. Analyze cluster stability across different time windows\n\n",
                   
                   "========================================\n",
                   "Hard Cluster Analysis completed successfully! 🎉\n",
                   "All cluster comparison results are now available.\n",
                   "========================================\n"
  )
  
  # 将报告写入文件
  writeLines(report, "time_window_hard_cluster_analysis_report.txt")
  
  # 显示报告
  cat(report)
  
  return(report)
}

# 生成最终报告
if(nrow(all_cluster_results) > 0) {
  final_report <- generate_hard_cluster_analysis_report(all_cluster_results, refined_results)
} else {
  cat("❌ 没有分析结果可生成报告\n")
}



# 显示生成的文件
cat("\n📁 生成的主要文件:\n")
main_files <- c(
  "time_window_hard_cluster_analysis_complete.csv",
  "time_window_hard_cluster_analysis_refined_significant.csv",
  "time_window_hard_cluster_analysis_original_significant.csv",
  "time_window_hard_cluster_analysis_report.txt"
)

for (file in main_files) {
  if (file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  } else {
    cat(sprintf("❌ %s (未找到)\n", file))
  }
}

cat("\n📊 生成的可视化文件:\n")
viz_files <- list.files(pattern = "\\.(pdf|png)$")
if(length(viz_files) > 0) {
  for(file in viz_files) {
    cat(sprintf("✓ %s\n", file))
  }
} else {
  cat("❌ 未找到可视化文件\n")
}

cat("\n🎯 关键发现总结:\n")
if(nrow(all_cluster_results) > 0) {
  cat(sprintf("- 总共测试了 %d 个cluster比较\n", nrow(all_cluster_results)))
  cat(sprintf("- 发现 %d 个显著差异 (p < 0.05)\n", sum(all_cluster_results$T_Significant)))
  if(exists("refined_results") && nrow(refined_results) > 0) {
    cat(sprintf("- FDR校正后剩余 %d 个显著差异\n", sum(refined_results$Welch_Significant_FDR_Strategy2)))
  }
  
  large_effect_count <- sum(all_cluster_results$Effect_Size_Category == "Large")
  if(large_effect_count > 0) {
    cat(sprintf("- 发现 %d 个大效应差异 (Cohen's d ≥ 0.8)\n", large_effect_count))
  }
} else {
  cat("- 未检测到cluster差异\n")
}


