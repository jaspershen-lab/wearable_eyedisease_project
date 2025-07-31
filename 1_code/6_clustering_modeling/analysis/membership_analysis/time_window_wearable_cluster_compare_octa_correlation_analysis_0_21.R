# 基于时间窗口硬聚类的组间差异分析 - 修改版（仅0_21数据）
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

# ================== 2. 读取和处理OCTA数据（修改：仅保留0_21数据）==================

if(!exists("octa_improvements") || !exists("va_improvements")) {
  
  cat("===== 处理OCTA和VA数据（仅0_21区域）=====\n")
  
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
  
  # 修改：仅筛选0_21区域的参数
  filter_key_octa_params <- function(data, param_type = "bloodflow") {
    if(param_type == "bloodflow") {
      layers <- c("SVP", "ICP", "DCP", "Choroid")
    } else {
      layers <- c("GCL.IPL", "INL", "Retina")
    }
    
    # 修改：仅保留0_21区域，去掉0_6
    regions <- c("0_21")  # 只保留0_21
    pattern <- paste0("(", paste(layers, collapse = "|"), ").*(",
                      paste(regions, collapse = "|"), ")_T0$")
    
    params_T0 <- names(data)[grep(pattern, names(data))]
    params_T2 <- gsub("_T0$", "_T2", params_T0)
    params_T2 <- params_T2[params_T2 %in% names(data)]
    
    valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
    
    cat(sprintf("🔍 筛选到 %s 参数: %d个 (仅0_21区域)\n", param_type, length(valid_base_params)))
    
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
  
  cat("✓ OCTA data processed (仅0_21区域)\n")
  
  # 处理VA数据（保持不变）
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

# ================== 4. 修改后的硬聚类组间差异分析函数（去掉区域校正）==================

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

# 获取outcome参数（确保只包含0_21区域的参数）
octa_improvement_params <- names(enhanced_cluster_analysis)[grep("_improvement$", names(enhanced_cluster_analysis))]
octa_improvement_params <- octa_improvement_params[!grepl("vision_improvement", octa_improvement_params)]

# 验证参数确实都是0_21区域的
cat("🔍 验证OCTA参数（应该都包含0_21）:\n")
octa_21_params <- octa_improvement_params[grepl("0_21", octa_improvement_params)]
octa_6_params <- octa_improvement_params[grepl("0_6", octa_improvement_params)]

cat(sprintf("✓ 0_21区域参数: %d个\n", length(octa_21_params)))
cat(sprintf("⚠️ 0_6区域参数: %d个 (应该为0)\n", length(octa_6_params)))

if(length(octa_6_params) > 0) {
  cat("发现0_6参数，将被排除:\n")
  print(octa_6_params)
  octa_improvement_params <- octa_21_params  # 只保留0_21参数
}

va_improvement_params <- c("vision_improvement_1w", "vision_improvement_1m")

cat("🔬 开始硬聚类组间差异分析（仅0_21区域）...\n")
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

# ================== 6. 简化的FDR校正（去掉区域分组）==================

perform_simplified_fdr_correction <- function(all_cluster_results) {
  
  cat("=== 执行简化的FDR校正（不按区域分组）===\n")
  
  if(nrow(all_cluster_results) == 0) {
    cat("❌ 没有结果进行校正\n")
    return(data.frame())
  }
  
  # 简化的参数分类（去掉区域分组）
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
      # 主要参数类型
      Parameter_Major_Type = case_when(
        grepl("BloodFlow", Parameter_Type_Detailed) ~ "BloodFlow",
        grepl("Thickness", Parameter_Type_Detailed) ~ "Thickness",
        grepl("Vision", Parameter_Type_Detailed) ~ "Vision",
        TRUE ~ "Other"
      )
    )
  
  # 显示分类结果
  cat("\n📋 参数分类统计:\n")
  classification_summary <- refined_results %>%
    count(Time_Window, Parameter_Type_Detailed, name = "Tests") %>%
    arrange(Time_Window, Parameter_Type_Detailed)
  
  print(classification_summary)
  
  # 策略1: 按时间窗口 + 详细参数类型分别校正
  cat("\n🔬 策略1: 按时间窗口 + 详细参数类型分别校正\n")
  
  strategy1_results <- refined_results %>%
    group_by(Time_Window, Parameter_Type_Detailed) %>%
    mutate(
      Welch_P_Value_FDR_Strategy1 = p.adjust(Welch_P_Value, method = "fdr"),
      Welch_Significant_FDR_Strategy1 = Welch_P_Value_FDR_Strategy1 < 0.05,
      Tests_in_Group_Strategy1 = n()
    ) %>%
    ungroup()
  
  # 策略2: 按时间窗口 + 主要参数类型分别校正
  cat("🔬 策略2: 按时间窗口 + 主要参数类型分别校正\n")
  
  strategy2_results <- strategy1_results %>%
    group_by(Time_Window, Parameter_Major_Type) %>%
    mutate(
      Welch_P_Value_FDR_Strategy2 = p.adjust(Welch_P_Value, method = "fdr"),
      Welch_Significant_FDR_Strategy2 = Welch_P_Value_FDR_Strategy2 < 0.05,
      Tests_in_Group_Strategy2 = n()
    ) %>%
    ungroup()
  
  # 策略3: 仅按主要参数类型校正（不考虑时间窗口）
  cat("🔬 策略3: 仅按主要参数类型校正\n")
  
  final_results <- strategy2_results %>%
    group_by(Parameter_Major_Type) %>%
    mutate(
      Welch_P_Value_FDR_Strategy3 = p.adjust(Welch_P_Value, method = "fdr"),
      Welch_Significant_FDR_Strategy3 = Welch_P_Value_FDR_Strategy3 < 0.05,
      Tests_in_Group_Strategy3 = n()
    ) %>%
    ungroup() %>%
    # 策略4: 全局校正
    mutate(
      Welch_P_Value_FDR_Global = p.adjust(Welch_P_Value, method = "fdr"),
      Welch_Significant_FDR_Global = Welch_P_Value_FDR_Global < 0.05,
      Tests_in_Global = n()
    )
  
  # 比较不同策略的结果
  cat("\n📊 不同FDR校正策略的结果比较:\n")
  
  strategy_summary <- data.frame(
    Strategy = c(
      "原始显著(Welch)", "策略1(时间窗口+详细类型)", "策略2(时间窗口+主要类型)", 
      "策略3(仅主要类型)", "全局校正"
    ),
    Significant_Count = c(
      sum(final_results$Welch_Significant),
      sum(final_results$Welch_Significant_FDR_Strategy1),
      sum(final_results$Welch_Significant_FDR_Strategy2),
      sum(final_results$Welch_Significant_FDR_Strategy3),
      sum(final_results$Welch_Significant_FDR_Global)
    ),
    Total_Tests = rep(nrow(final_results), 5)
  ) %>%
    mutate(
      Success_Rate = round(Significant_Count / Total_Tests * 100, 1)
    )
  
  print(strategy_summary)
  
  # 推荐策略（考虑小样本量）
  cat("\n💡 推荐使用策略2（按时间窗口 + 主要参数类型校正）:\n")
  
  # 选择Welch t检验作为主要结果
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
      cat(sprintf("   📋 校正组: %s_%s (%d个测试)\n", 
                  result$Time_Window, result$Parameter_Major_Type, result$Tests_in_Group_Strategy2))
      if(!is.na(result$Percent_Improvement)) {
        cat(sprintf("   📈 相对改善: %.1f%%\n", result$Percent_Improvement))
      }
    }
  } else {
    cat("   在推荐策略下仍无显著结果\n")
    
    # 显示最接近显著的结果
    near_significant <- final_results %>%
      filter(Welch_Significant) %>%
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

# 执行简化FDR校正
if(nrow(all_cluster_results) > 0) {
  refined_results <- perform_simplified_fdr_correction(all_cluster_results)
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
  write.csv(refined_results, "time_window_hard_cluster_analysis_0_21_only_complete.csv", row.names = FALSE)
  
  # 提取推荐策略的显著结果
  recommended_significant_results <- refined_results %>%
    filter(Welch_Significant_FDR_Strategy2) %>%
    arrange(Welch_P_Value_FDR_Strategy2)
  
  write.csv(recommended_significant_results, "time_window_hard_cluster_analysis_0_21_only_refined_significant.csv", row.names = FALSE)
  
  # 提取原始显著结果
  original_significant_results <- refined_results %>%
    filter(Welch_Significant) %>%
    arrange(Welch_P_Value)
  
  write.csv(original_significant_results, "time_window_hard_cluster_analysis_0_21_only_original_significant.csv", row.names = FALSE)
  
  cat("✅ 简化硬聚类分析完成（仅0_21区域）！\n")
  cat("📁 结果已保存：\n")
  cat("   - time_window_hard_cluster_analysis_0_21_only_complete.csv (完整结果)\n")
  cat("   - time_window_hard_cluster_analysis_0_21_only_refined_significant.csv (FDR校正后显著)\n")
  cat("   - time_window_hard_cluster_analysis_0_21_only_original_significant.csv (原始显著)\n")
} else {
  cat("❌ 没有硬聚类分析结果保存\n")
}

# ================== 8. 可视化函数（修改后去掉区域标注）==================

create_simplified_hard_cluster_visualizations <- function(data, analysis_results, fdr_threshold = 0.05) {
  
  cat("\n===== 创建简化的硬聚类分析可视化（仅0_21区域）=====\n")
  
  if(nrow(analysis_results) == 0) {
    cat("❌ 没有分析结果可用于可视化\n")
    return(NULL)
  }
  
  # 优先使用FDR校正后显著的结果
  if("Welch_Significant_FDR_Strategy2" %in% names(analysis_results)) {
    # 使用FDR校正后的结果
    significant_for_viz <- analysis_results %>%
      filter(Welch_Significant_FDR_Strategy2) %>%
      arrange(desc(Cohens_D)) %>%
      slice_head(n = 12)
    
    cat(sprintf("📊 使用FDR校正后的显著结果创建可视化 (%d个变量)\n", nrow(significant_for_viz)))
    
    # 如果FDR校正后没有显著结果，回退到原始显著结果
    if(nrow(significant_for_viz) == 0) {
      cat("⚠️ FDR校正后无显著结果，回退到原始显著结果\n")
      significant_for_viz <- analysis_results %>%
        filter(Welch_Significant) %>%
        arrange(desc(Cohens_D)) %>%
        slice_head(n = 8)
      use_fdr_correction <- FALSE
    } else {
      use_fdr_correction <- TRUE
    }
    
  } else {
    # 如果没有FDR校正结果，使用原始显著结果
    cat("⚠️ 没有FDR校正结果，使用原始显著结果\n")
    significant_for_viz <- analysis_results %>%
      filter(Welch_Significant) %>%
      arrange(desc(Cohens_D)) %>%
      slice_head(n = 8)
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
        
        # 创建参数的清晰名称（简化，去掉区域信息）
        param_clean <- param %>%
          gsub("_improvement", " Improvement", .) %>%
          gsub("_0_21_", " ", .) %>%  # 去掉0_21标识
          gsub("_", " ", .) %>%
          gsub("vision improvement", "Vision Improvement", .) %>%
          gsub("1w", "(1 Week)", .) %>%
          gsub("1m", "(1 Month)", .) %>%
          gsub("SVP", "Superficial Vascular Plexus", .) %>%
          gsub("ICP", "Intermediate Capillary Plexus", .) %>%
          gsub("DCP", "Deep Capillary Plexus", .) %>%
          gsub("GCL IPL", "GCL-IPL", .)
        
        window_clean <- gsub("_", " ", window) %>% str_to_title()
        
        # 标题和统计信息显示FDR校正状态
        if(use_fdr_correction) {
          # 使用FDR校正后的p值
          p_value_to_show <- result$Welch_P_Value_FDR_Strategy2
          p_label <- "p_FDR"
          subtitle_text <- paste0("FDR-Corrected Cluster Comparison (0_21 region only) | ",
                                  "Welch t = ", round(result$Welch_T_Statistic, 3), 
                                  ", ", p_label, " = ", format.pval(p_value_to_show, digits = 3),
                                  " | Cohen's d = ", round(result$Cohens_D, 3),
                                  " (", result$Effect_Size_Category, " effect)",
                                  " | n = ", result$N_Total)
        } else {
          # 使用原始p值
          p_value_to_show <- result$Welch_P_Value
          p_label <- "p"
          subtitle_text <- paste0("Original Cluster Comparison (0_21 region only) | ",
                                  "Welch t = ", round(result$Welch_T_Statistic, 3), 
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
                            "| 0_21 region only",
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
        ggsave(paste0("hard_cluster_", gsub("[^A-Za-z0-9]", "_", param), "_", window, "_0_21_only", filename_suffix, "_boxplot.pdf"),
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
                             "Time Window Hard Cluster Analysis - FDR Corrected (0_21 Region Only)",
                             "Time Window Hard Cluster Analysis - Original Significant (0_21 Region Only)")
    
    # 组合图形
    combined_plot <- do.call(gridExtra::grid.arrange, 
                             c(plot_list, 
                               ncol = ncol_layout,
                               top = combined_title))
    
    # 保存组合图形
    filename_suffix <- ifelse(use_fdr_correction, "_FDR_corrected", "_original")
    ggsave(paste0("time_window_hard_cluster_analysis_0_21_only_combined", filename_suffix, ".pdf"),
           combined_plot, width = 14, height = 6 * nrow_layout)
    
    cat("✓ 可视化已保存\n")
    
    # FDR校正状态报告
    if(use_fdr_correction) {
      cat("\n🎯 FDR校正状态报告:\n")
      cat("✅ 使用FDR校正后的p值 (Welch t-test, Strategy 2)\n")
      cat(sprintf("✅ 显示了 %d 个FDR校正后显著的变量\n", length(plot_list)))
      cat("✅ 所有p值都经过了多重比较校正\n")
      cat("✅ 仅包含0_21区域的参数\n")
    } else {
      cat("\n⚠️ 使用原始p值 (未经FDR校正)\n")
      cat(sprintf("⚠️ 显示了 %d 个原始显著的变量\n", length(plot_list)))
      cat("⚠️ 建议检查FDR校正结果\n")
      cat("✅ 仅包含0_21区域的参数\n")
    }
    
    return(combined_plot)
  }
  
  return(NULL)
}

# ================== 9. 创建可视化 ==================

if(nrow(all_cluster_results) > 0) {
  
  cat("\n📊 创建简化的硬聚类分析可视化（仅0_21区域）...\n")
  
  # 如果有FDR校正结果，使用FDR校正后的结果
  if(exists("refined_results") && nrow(refined_results) > 0) {
    simplified_plots <- create_simplified_hard_cluster_visualizations(
      data = enhanced_cluster_analysis,
      analysis_results = refined_results,
      fdr_threshold = 0.05
    )
  } else {
    # 如果没有FDR校正结果，使用原始结果
    cat("⚠️ 没有FDR校正结果，使用原始结果创建可视化\n")
    simplified_plots <- create_simplified_hard_cluster_visualizations(
      data = enhanced_cluster_analysis,
      analysis_results = all_cluster_results,
      fdr_threshold = 0.05
    )
  }
  
} else {
  cat("❌ 没有硬聚类分析结果用于可视化\n")
}

# ================== 10. 效应大小分析 ==================

analyze_effect_sizes_simplified <- function(analysis_results) {
  
  if(nrow(analysis_results) == 0) {
    cat("❌ 没有结果进行效应大小分析\n")
    return(NULL)
  }
  
  cat("\n===== 效应大小深度分析（仅0_21区域）=====\n")
  
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
    dplyr::select(Time_Window, Outcome_Parameter, Cohens_D, Welch_P_Value, Better_Cluster, Percent_Improvement)
  
  if(nrow(large_effects) > 0) {
    cat(sprintf("\n🎯 大效应结果 (%d个):\n", nrow(large_effects)))
    print(large_effects)
  }
  
  # 创建效应大小可视化
  p_effect <- ggplot(analysis_results, aes(x = Cohens_D, y = -log10(Welch_P_Value))) +
    geom_point(aes(color = Effect_Size_Category, size = abs(Percent_Improvement)), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(0.2, 0.5, 0.8), linetype = "dashed", color = "gray") +
    scale_color_manual(
      values = c("Large" = "#D73027", "Medium" = "#FC8D59", "Small" = "#91BFDB", "Negligible" = "#EEEEEE"),
      name = "Effect Size"
    ) +
    scale_size_continuous(name = "% Improvement", range = c(1, 5)) +
    labs(
      title = "Effect Size vs Statistical Significance (0_21 Region Only)",
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
  
  ggsave("effect_size_volcano_plot_0_21_only.pdf", p_effect, width = 10, height = 8)
  
  return(list(
    distribution = effect_size_dist,
    large_effects = large_effects,
    volcano_plot = p_effect
  ))
}

# 执行效应大小分析
if(nrow(all_cluster_results) > 0) {
  effect_analysis <- analyze_effect_sizes_simplified(all_cluster_results)
}

# ================== 11. 生成最终报告 ==================

generate_simplified_cluster_analysis_report <- function(analysis_results, refined_results = NULL) {
  
  report <- paste0(
    "========================================\n",
    "TIME WINDOW HARD CLUSTER ANALYSIS REPORT\n",
    "        (0_21 REGION ONLY)\n",
    "========================================\n\n",
    
    "🎯 ANALYSIS OVERVIEW:\n",
    "✅ 针对硬聚类（所有membership=1）进行组间差异分析\n",
    "✅ 仅分析0_21区域（21mm广域视网膜区）的OCTA参数\n",
    "✅ 去掉了按区域的FDR校正分组\n",
    "✅ 使用Welch t检验作为主要统计方法\n",
    "✅ 实施简化的FDR校正策略\n\n",
    
    "📊 ANALYSIS DETAILS:\n",
    "- Analysis Date: ", Sys.Date(), "\n",
    "- Total Comparisons: ", nrow(analysis_results), "\n",
    "- Region: 0_21 (21mm diameter) only\n"
  )
  
  if(nrow(analysis_results) > 0) {
    report <- paste0(report,
                     "- Welch t-test Significant (p < 0.05): ", sum(analysis_results$Welch_Significant), "\n")
    
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
      filter(Welch_Significant) %>%
      arrange(desc(Cohens_D)) %>%
      slice_head(n = 5)
    
    if(nrow(significant_results) > 0) {
      report <- paste0(report, "\n🏆 TOP FINDINGS (0_21 Region):\n")
      
      for(i in 1:nrow(significant_results)) {
        result <- significant_results[i, ]
        report <- paste0(report,
                         sprintf("\n%d. %s - %s Window:\n", i, result$Outcome_Parameter, result$Time_Window),
                         sprintf("   📊 Difference: %.3f vs %.3f\n", result$Cluster1_Mean, result$Cluster2_Mean),
                         sprintf("   📈 Statistics: Welch t=%.3f, p=%.4f\n", result$Welch_T_Statistic, result$Welch_P_Value),
                         sprintf("   🔬 Effect: Cohen's d=%.3f (%s)\n", result$Cohens_D, result$Effect_Size_Category),
                         sprintf("   🎯 Better Cluster: %s (n=%d)\n", result$Better_Cluster, result$Better_Cluster_N))
        
        if(!is.na(result$Percent_Improvement)) {
          report <- paste0(report, sprintf("   📈 Improvement: %.1f%%\n", result$Percent_Improvement))
        }
      }
    }
    
    # 时间窗口表现
    window_performance <- analysis_results %>%
      filter(Welch_Significant) %>%
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
                   "- time_window_hard_cluster_analysis_0_21_only_complete.csv: Complete results\n",
                   "- time_window_hard_cluster_analysis_0_21_only_refined_significant.csv: FDR significant\n",
                   "- time_window_hard_cluster_analysis_0_21_only_original_significant.csv: Original significant\n",
                   "- hard_cluster_*_0_21_only*.pdf: Individual cluster comparison plots\n",
                   "- time_window_hard_cluster_analysis_0_21_only_combined*.pdf: Combined visualization\n",
                   "- effect_size_volcano_plot_0_21_only.pdf: Effect size visualization\n\n",
                   
                   "🎯 KEY MODIFICATIONS:\n",
                   "1. ✅ Removed 0_6 (macular) region parameters\n",
                   "2. ✅ Focus on 0_21 (widefield) region only\n",
                   "3. ✅ Simplified FDR correction (no region-based grouping)\n",
                   "4. ✅ Reduced parameter complexity while maintaining statistical rigor\n\n",
                   
                   "🔬 FDR CORRECTION STRATEGIES:\n",
                   "- Strategy 1: By Time Window + Detailed Parameter Type\n",
                   "- Strategy 2: By Time Window + Major Parameter Type (Recommended)\n", 
                   "- Strategy 3: By Major Parameter Type only\n",
                   "- Global: All parameters together (Most conservative)\n\n",
                   
                   "💡 CLINICAL IMPLICATIONS:\n",
                   "1. Focus on widefield (0_21) retinal health patterns\n",
                   "2. Better clusters indicate favorable recovery in broader retinal areas\n",
                   "3. Time window differences suggest optimal intervention periods\n",
                   "4. Simplified analysis reduces multiple comparison burden\n",
                   "5. More focused results for clinical interpretation\n\n",
                   
                   "🚀 NEXT STEPS:\n",
                   "1. Validate 0_21 region findings in independent cohorts\n",
                   "2. Investigate baseline predictors of better-performing clusters\n",
                   "3. Develop widefield-focused treatment protocols\n",
                   "4. Consider peripheral retinal intervention strategies\n",
                   "5. Compare 0_21 vs 0_6 region recovery patterns in future studies\n\n",
                   
                   "========================================\n",
                   "Simplified Hard Cluster Analysis completed! 🎉\n",
                   "All 0_21 region cluster comparison results available.\n",
                   "========================================\n"
  )
  
  # 将报告写入文件
  writeLines(report, "time_window_hard_cluster_analysis_0_21_only_report.txt")
  
  # 显示报告
  cat(report)
  
  return(report)
}

# 生成最终报告
if(nrow(all_cluster_results) > 0) {
  final_report <- generate_simplified_cluster_analysis_report(all_cluster_results, refined_results)
} else {
  cat("❌ 没有分析结果可生成报告\n")
}

# 显示生成的文件
cat("\n📁 生成的主要文件（仅0_21区域）:\n")
main_files <- c(
  "time_window_hard_cluster_analysis_0_21_only_complete.csv",
  "time_window_hard_cluster_analysis_0_21_only_refined_significant.csv", 
  "time_window_hard_cluster_analysis_0_21_only_original_significant.csv",
  "time_window_hard_cluster_analysis_0_21_only_report.txt"
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

cat("\n🎯 关键发现总结（仅0_21区域）:\n")
if(nrow(all_cluster_results) > 0) {
  cat(sprintf("- 总共测试了 %d 个cluster比较（仅0_21区域）\n", nrow(all_cluster_results)))
  cat(sprintf("- 发现 %d 个显著差异 (p < 0.05)\n", sum(all_cluster_results$Welch_Significant)))
  if(exists("refined_results") && nrow(refined_results) > 0) {
    cat(sprintf("- FDR校正后剩余 %d 个显著差异\n", sum(refined_results$Welch_Significant_FDR_Strategy2)))
  }
  
  large_effect_count <- sum(all_cluster_results$Effect_Size_Category == "Large")
  if(large_effect_count > 0) {
    cat(sprintf("- 发现 %d 个大效应差异 (Cohen's d ≥ 0.8)\n", large_effect_count))
  }
  
  cat("- ✅ 分析简化：去除了0_6区域，专注于0_21广域区域\n")
  cat("- ✅ 校正优化：移除了区域分组，减少多重比较负担\n")
} else {
  cat("- 未检测到cluster差异\n")
}

cat("\n💡 修改总结:\n")
cat("✅ 1. 去除了所有0_6区域的参数\n")
cat("✅ 2. 只保留0_21区域的OCTA参数\n") 
cat("✅ 3. 简化了FDR校正，去掉了按区域分组\n")
cat("✅ 4. 减少了参数数量，降低了多重比较负担\n")
cat("✅ 5. 保持了统计严格性，专注于广域视网膜分析\n")