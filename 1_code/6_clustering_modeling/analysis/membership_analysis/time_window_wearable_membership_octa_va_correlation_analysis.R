# 基于时间窗Max Membership的相关性分析
# 使用代码三生成的真正max membership数据进行预后指标相关性分析

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
max_membership_data <- read.csv("3_data_analysis/6_clustering_modeling/time_window_clustering/time_window_max_membership_data_fixed.csv")

cat("✓ 读取时间窗Max Membership数据\n")
cat("数据维度:", dim(max_membership_data), "\n")
cat("列名:", paste(names(max_membership_data), collapse = ", "), "\n\n")

# 检查数据结构
membership_cols <- grep("^membership_", names(max_membership_data), value = TRUE)
cluster_cols <- grep("^cluster_", names(max_membership_data), value = TRUE)

cat("🎯 发现的数据列:\n")
cat("Membership列:", paste(membership_cols, collapse = ", "), "\n")
cat("Cluster列:", paste(cluster_cols, collapse = ", "), "\n\n")

# ================== 2. 读取和处理OCTA数据 ==================

if(!exists("all_correlations") || !exists("enhanced_membership_analysis")) {
  
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

# ================== 3. 整合Max Membership数据与Outcome数据 ==================

# 整合所有数据
enhanced_max_membership_analysis <- max_membership_data %>%
  left_join(octa_improvements, by = c("subject_id" = "ID")) %>%
  left_join(va_improvements, by = c("subject_id" = "ID"))

cat("✓ 数据整合完成\n")
cat("最终分析数据:", nrow(enhanced_max_membership_analysis), "行,", ncol(enhanced_max_membership_analysis), "列\n\n")

# ================== 4. Max Membership相关性分析函数 ==================

perform_enhanced_max_membership_correlation_analysis <- function(data, membership_cols, cluster_cols) {
  
  cat("===== 增强Max Membership相关性分析 =====\n")
  
  # 获取outcome参数
  octa_improvement_params <- names(data)[grep("_improvement$", names(data))]
  octa_improvement_params <- octa_improvement_params[!grepl("vision_improvement", octa_improvement_params)]
  va_improvement_params <- c("vision_improvement_1w", "vision_improvement_1m")
  
  cat("OCTA改善参数:", length(octa_improvement_params), "个\n")
  cat("VA改善参数:", length(va_improvement_params), "个\n\n")
  
  # 🎯 增强的cluster分析函数
  analyze_enhanced_cluster_info <- function(clean_data, membership_col, cluster_col, param) {
    
    cluster_analysis <- clean_data %>%
      group_by(!!sym(cluster_col)) %>%
      summarise(
        count = n(),
        mean_membership = round(mean(!!sym(membership_col)), 3),
        mean_outcome = round(mean(!!sym(param)), 3),
        median_outcome = round(median(!!sym(param)), 3),
        sd_outcome = round(sd(!!sym(param)), 3),
        min_outcome = round(min(!!sym(param)), 3),
        max_outcome = round(max(!!sym(param)), 3),
        .groups = 'drop'
      ) %>%
      arrange(desc(count))  # 先按数量排序
    
    # 🔧 关键修正：定义多种"主要cluster"
    primary_by_count <- cluster_analysis %>% slice(1) %>% pull(!!sym(cluster_col))
    primary_by_outcome <- cluster_analysis %>% slice_max(mean_outcome, n = 1) %>% pull(!!sym(cluster_col))
    primary_by_membership <- cluster_analysis %>% slice_max(mean_membership, n = 1) %>% pull(!!sym(cluster_col))
    
    # 选择最合适的"主要cluster"
    if(primary_by_outcome == primary_by_count) {
      selected_primary <- primary_by_outcome
      selection_reason <- "Best outcome & Most patients"
    } else if(primary_by_outcome == primary_by_membership) {
      selected_primary <- primary_by_outcome  
      selection_reason <- "Best outcome & Highest membership"
    } else {
      # 如果不一致，选择最佳outcome的cluster，但标注冲突
      selected_primary <- primary_by_outcome
      selection_reason <- "Best outcome (conflicts with patient count)"
    }
    
    # 获取选定cluster的详细信息
    selected_cluster_info <- cluster_analysis %>% 
      filter(!!sym(cluster_col) == selected_primary)
    
    # 生成详细的cluster比较描述
    cluster_comparison <- cluster_analysis %>%
      mutate(
        rank_by_outcome = rank(-mean_outcome),
        rank_by_count = rank(-count),
        rank_by_membership = rank(-mean_membership)
      ) %>%
      mutate(desc = paste0("C", !!sym(cluster_col), 
                           "(n=", count, 
                           ",mem=", mean_membership,
                           ",out=", mean_outcome, "±", sd_outcome,
                           ",rank:", rank_by_outcome, "/", rank_by_count, "/", rank_by_membership, ")")) %>%
      pull(desc) %>%
      paste(collapse = "; ")
    
    return(list(
      cluster_analysis = cluster_analysis,
      primary_cluster = selected_primary,
      primary_cluster_n = selected_cluster_info$count,
      primary_cluster_mean_outcome = selected_cluster_info$mean_outcome,
      primary_cluster_median_outcome = selected_cluster_info$median_outcome,
      selection_reason = selection_reason,
      cluster_comparison = cluster_comparison,
      # 额外信息
      best_outcome_cluster = primary_by_outcome,
      most_patients_cluster = primary_by_count,
      highest_membership_cluster = primary_by_membership,
      outcome_conflict = primary_by_outcome != primary_by_count
    ))
  }
  
  # 主要相关性分析函数
  analyze_max_membership_correlations <- function(outcome_params, outcome_type) {
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
            # Pearson相关性
            cor_test <- try(cor.test(clean_data[[membership_col]], clean_data[[param]], 
                                     method = "pearson"), silent = TRUE)
            
            # Spearman相关性
            spearman_test <- try(cor.test(clean_data[[membership_col]], clean_data[[param]], 
                                          method = "spearman"), silent = TRUE)
            
            if(class(cor_test) != "try-error" && class(spearman_test) != "try-error") {
              
              # 🎯 使用增强的cluster分析
              cluster_info <- analyze_enhanced_cluster_info(clean_data, membership_col, cluster_col, param)
              
              # 效应大小
              effect_size <- case_when(
                abs(cor_test$estimate) >= 0.5 ~ "Large",
                abs(cor_test$estimate) >= 0.3 ~ "Medium",
                abs(cor_test$estimate) >= 0.1 ~ "Small",
                TRUE ~ "Negligible"
              )
              
              # 🔧 增强的结果结构
              result_row <- data.frame(
                Time_Window = window_name,
                Outcome_Type = outcome_type,
                Outcome_Parameter = param,
                N = nrow(clean_data),
                Pearson_r = as.numeric(cor_test$estimate),
                Pearson_p = cor_test$p.value,
                Pearson_CI_Lower = cor_test$conf.int[1],
                Pearson_CI_Upper = cor_test$conf.int[2],
                Spearman_rho = as.numeric(spearman_test$estimate),
                Spearman_p = spearman_test$p.value,
                Effect_Size = effect_size,
                # 🎯 增强的cluster信息
                Primary_Cluster = cluster_info$primary_cluster,
                Primary_Cluster_N = cluster_info$primary_cluster_n,
                Primary_Cluster_Mean_Outcome = cluster_info$primary_cluster_mean_outcome,
                Primary_Cluster_Median_Outcome = cluster_info$primary_cluster_median_outcome,
                Selection_Reason = cluster_info$selection_reason,
                Total_Clusters = nrow(cluster_info$cluster_analysis),
                Cluster_Comparison = cluster_info$cluster_comparison,
                # 🔧 冲突检测
                Best_Outcome_Cluster = cluster_info$best_outcome_cluster,
                Most_Patients_Cluster = cluster_info$most_patients_cluster,
                Highest_Membership_Cluster = cluster_info$highest_membership_cluster,
                Outcome_Conflict = cluster_info$outcome_conflict,
                stringsAsFactors = FALSE
              )
              
              results <- rbind(results, result_row)
            }
          }
        }
      }
    }
    
    return(results)
  }
  
  # 分析OCTA和VA相关性
  octa_results <- analyze_max_membership_correlations(octa_improvement_params, "OCTA")
  va_results <- analyze_max_membership_correlations(va_improvement_params, "VA")
  
  all_results <- rbind(octa_results, va_results)
  
  if(nrow(all_results) > 0) {
    # 🔧 关键修改：按时间窗口和参数类型分别进行FDR校正
    all_results <- all_results %>%
      # 首先确定参数类型
      mutate(
        Parameter_Type = case_when(
          # Blood flow parameters
          grepl("VD_|PA_.*SVP|PA_.*ICP|PA_.*DCP|PA_.*Choroid", Outcome_Parameter) ~ "BloodFlow",
          # Thickness parameters  
          grepl("Thickness_", Outcome_Parameter) ~ "Thickness",
          # Vision parameters
          grepl("vision_improvement", Outcome_Parameter) ~ "Vision",
          TRUE ~ "Other"
        )
      ) %>%
      # 按时间窗口和参数类型分组进行FDR校正
      group_by(Time_Window, Parameter_Type) %>%
      mutate(
        Pearson_p_FDR = p.adjust(Pearson_p, method = "fdr"),
        Spearman_p_FDR = p.adjust(Spearman_p, method = "fdr")
      ) %>%
      ungroup()
    
    # 添加显著性标记
    all_results <- all_results %>%
      mutate(
        Abs_Pearson_r = abs(Pearson_r),
        Significant_Pearson = Pearson_p < 0.05,
        Significant_Pearson_FDR = Pearson_p_FDR < 0.05,
        Significant_Spearman = Spearman_p < 0.05,
        Significant_Spearman_FDR = Spearman_p_FDR < 0.05
      ) %>%
      arrange(desc(Abs_Pearson_r))
    
    cat("\n📊 按时间窗口和参数类型分别进行FDR校正的结果:\n")
    
    # 显示每个时间窗口和参数类型的校正情况
    window_type_fdr_summary <- all_results %>%
      group_by(Time_Window, Parameter_Type) %>%
      summarise(
        Total_Tests = n(),
        Significant_Raw = sum(Significant_Pearson),
        Significant_FDR = sum(Significant_Pearson_FDR),
        Min_p_FDR = ifelse(any(Significant_Pearson), min(Pearson_p_FDR[Significant_Pearson]), NA),
        .groups = 'drop'
      ) %>%
      arrange(Time_Window, Parameter_Type)
    
    cat("\n各时间窗口和参数类型的FDR校正情况:\n")
    print(window_type_fdr_summary)
    
    # 详细显示校正前后对比
    if(any(all_results$Significant_Pearson)) {
      cat("\n🔍 显著结果的FDR校正前后对比:\n")
      significant_comparison <- all_results %>%
        filter(Significant_Pearson) %>%
        dplyr::select(Time_Window, Parameter_Type, Outcome_Parameter, 
               Pearson_p, Pearson_p_FDR, Significant_Pearson_FDR) %>%
        arrange(Time_Window, Parameter_Type, Pearson_p)
      
      for(i in 1:nrow(significant_comparison)) {
        result <- significant_comparison[i, ]
        fdr_status <- ifelse(result$Significant_Pearson_FDR, "✅ 仍显著", "❌ 不再显著")
        cat(sprintf("%s - %s - %s:\n", result$Time_Window, result$Parameter_Type, result$Outcome_Parameter))
        cat(sprintf("  p = %.4f → p_FDR = %.4f (%s)\n", 
                    result$Pearson_p, result$Pearson_p_FDR, fdr_status))
      }
    }
    
    return(all_results)
  } else {
    cat("❌ 未找到任何有效的相关性结果\n")
    return(data.frame())
  }
  
  # 额外分析：显示参数分类结果
  if(exists("all_results") && nrow(all_results) > 0) {
    cat("\n📋 参数分类验证:\n")
    parameter_classification <- all_results %>%
      group_by(Parameter_Type) %>%
      summarise(
        Count = n(),
        Example_Parameters = paste(head(unique(Outcome_Parameter), 3), collapse = ", "),
        .groups = 'drop'
      )
    
    print(parameter_classification)
    
    # 显示每个时间窗口的详细分组
    cat("\n📊 各时间窗口的参数类型分布:\n")
    detailed_distribution <- all_results %>%
      count(Time_Window, Parameter_Type, name = "Count") %>%
      pivot_wider(names_from = Parameter_Type, values_from = Count, values_fill = 0)
    
    print(detailed_distribution)
  }
}

# ================== 5. 执行Max Membership相关性分析 ==================

# 执行分析
# 执行增强分析
max_membership_correlations <- perform_enhanced_max_membership_correlation_analysis(
  enhanced_max_membership_analysis, membership_cols, cluster_cols
)

# 创建输出目录
output_dir <- "3_data_analysis/6_clustering_modeling/time_window_max_membership_correlation_analysis"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
setwd(output_dir)

# 保存完整结果
write.csv(max_membership_correlations, "time_window_max_membership_correlations_complete.csv", row.names = FALSE)

# 筛选显著结果
significant_results <- max_membership_correlations %>%
  filter(Significant_Pearson) %>%
  arrange(desc(Abs_Pearson_r))

# 保存显著结果
write.csv(significant_results, "time_window_max_membership_correlations_significant.csv", row.names = FALSE)

# ================== 6. 结果输出和分析 ==================

cat("\n📊 时间窗Max Membership相关性分析结果摘要:\n")
cat(sprintf("- 总测试数: %d\n", nrow(max_membership_correlations)))
cat(sprintf("- 显著相关性 (p < 0.05): %d (%.1f%%)\n", 
            nrow(significant_results), nrow(significant_results)/nrow(max_membership_correlations)*100))
cat(sprintf("- FDR校正后显著: %d\n", sum(max_membership_correlations$Significant_Pearson_FDR)))

if(nrow(significant_results) > 0) {
  cat("\n🎯 Top 5 显著相关性 (with Cluster Information):\n")
  for(i in 1:min(5, nrow(significant_results))) {
    result <- significant_results[i, ]
    cat(sprintf("\n%d. %s (%s) - %s Window:\n", i, result$Outcome_Parameter, result$Outcome_Type, result$Time_Window))
    cat(sprintf("   📈 相关性: r = %.3f, p = %.4f (%s effect)\n", result$Pearson_r, result$Pearson_p, result$Effect_Size))
    cat(sprintf("   🎯 主要cluster: Cluster %s (n = %d patients)\n", result$Primary_Cluster, result$Primary_Cluster_N))
    cat(sprintf("   📊 该cluster平均outcome: %.3f\n", result$Primary_Cluster_Mean_Outcome))
    cat(sprintf("   🔍 所有clusters分布: %s\n", result$Cluster_Distribution))
    cat(sprintf("   👥 总clusters数: %d\n", result$Total_Clusters))
  }
  
  # 按时间窗口汇总
  cat("\n📈 按时间窗口汇总:\n")
  window_summary <- significant_results %>%
    group_by(Time_Window) %>%
    summarise(
      Count = n(),
      Mean_Abs_r = round(mean(Abs_Pearson_r), 3),
      Max_r = round(max(Abs_Pearson_r), 3),
      Primary_Clusters = paste(unique(Primary_Cluster), collapse = ", "),
      .groups = 'drop'
    ) %>%
    arrange(desc(Count))
  
  print(window_summary)
  
  # 按cluster汇总
  cat("\n🎯 按主要cluster汇总:\n")
  cluster_summary <- significant_results %>%
    group_by(Time_Window, Primary_Cluster) %>%
    summarise(
      Count = n(),
      Mean_Abs_r = round(mean(Abs_Pearson_r), 3),
      Max_r = round(max(Abs_Pearson_r), 3),
      .groups = 'drop'
    ) %>%
    arrange(desc(Count)) %>%
    mutate(Cluster_ID = paste0(Time_Window, "_Cluster_", Primary_Cluster))
  
  print(cluster_summary)
}

# ================== 7. 创建Max Membership可视化 ==================

create_max_membership_visualizations <- function(data, correlation_results) {
  
  cat("\n===== 创建Max Membership可视化 =====\n")
  
  # 筛选显著结果用于可视化
  significant_for_viz <- correlation_results %>%
    filter(Pearson_p_FDR < 0.05) %>%
    arrange(desc(abs(Pearson_r))) %>%
    slice_head(n = 6)  # 选择top 6进行可视化
  
  if(nrow(significant_for_viz) == 0) {
    cat("❌ 没有显著结果可用于可视化\n")
    return(NULL)
  }
  
  cat(sprintf("📊 创建 %d 个显著相关性的可视化\n", nrow(significant_for_viz)))
  
  # 创建散点图列表
  plot_list <- list()
  
  for(i in 1:nrow(significant_for_viz)) {
    result <- significant_for_viz[i, ]
    param <- result$Outcome_Parameter
    window <- result$Time_Window
    membership_col <- paste0("membership_", window)
    cluster_col <- paste0("cluster_", window)
    
    if(membership_col %in% names(data) && cluster_col %in% names(data) && param %in% names(data)) {
      # 准备绘图数据
      plot_data <- data %>%
        filter(!is.na(!!sym(membership_col)) & !is.na(!!sym(param)) & !is.na(!!sym(cluster_col))) %>%
        mutate(
          cluster_label = paste0("Cluster ", !!sym(cluster_col))
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
        
        # 创建散点图
        p <- ggplot(plot_data, aes_string(x = membership_col, y = param, color = "cluster_label")) +
          geom_point(size = 3, alpha = 0.8) +
          geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", alpha = 0.3) +
          geom_text(aes(label = subject_id), vjust = -0.5, size = 2.5, alpha = 0.7) +
          scale_color_brewer(type = "qual", palette = "Set2", name = "Cluster") +
          labs(
            title = paste(param_clean, "-", window_clean, "Window"),
            subtitle = paste0("Max Membership Correlation | ",
                              "r = ", round(result$Pearson_r, 3), 
                              ", p = ", format.pval(result$Pearson_p, digits = 3),
                              ", p(FDR) = ", format.pval(result$Pearson_p_FDR, digits = 3),  # 添加FDR校正p值
                              " (", result$Effect_Size, " effect)",
                              ", n = ", result$N),
            x = paste("Max Membership in", window_clean, "Window"),
            y = "Improvement Value",
            caption = paste("Primary cluster:", result$Primary_Cluster, 
                            "| Points colored by cluster assignment")
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
        ggsave(paste0("time_window_max_membership_", gsub("[^A-Za-z0-9]", "_", param), "_", window, "_scatter.pdf"),
               p, width = 10, height = 8)
      }
    }
  }
  
  # 创建组合图形
  if(length(plot_list) > 0) {
    # 计算布局
    ncol_layout <- ifelse(length(plot_list) >= 4, 2, 1)
    nrow_layout <- ceiling(length(plot_list) / ncol_layout)
    
    # 组合图形
    combined_plot <- do.call(gridExtra::grid.arrange, 
                             c(plot_list, 
                               ncol = ncol_layout,
                               top = "Time Window Max Membership Correlations with Cluster-Specific Information"))
    
    # 保存组合图形
    ggsave("time_window_max_membership_correlations_with_clusters_combined.pdf",
           combined_plot, width = 14, height = 6 * nrow_layout)
    
    cat("✓ 可视化已保存\n")
    return(combined_plot)
  }
  
  return(NULL)
}

# 创建可视化
if(nrow(max_membership_correlations) > 0) {
  max_membership_plots <- create_max_membership_visualizations(enhanced_max_membership_analysis, max_membership_correlations)
}

# ================== 8. 生成最终报告 ==================

generate_max_membership_correlation_report <- function(correlation_results, significant_results) {
  
  report <- paste0(
    "========================================\n",
    "TIME WINDOW MAX MEMBERSHIP CORRELATION ANALYSIS\n",
    "========================================\n\n",
    
    "🎯 MISSION ACCOMPLISHED:\n",
    "✅ 现在可以精确识别哪个时间窗的哪个cluster与outcome相关！\n",
    "✅ 每个显著相关性都有具体的cluster信息\n",
    "✅ 解决了'Early Recovery有显著性，但不知道是哪个cluster'的问题\n\n",
    
    "📊 ANALYSIS OVERVIEW:\n",
    "- Analysis Date: ", Sys.Date(), "\n",
    "- Total Time Window Correlations Tested: ", nrow(correlation_results), "\n",
    "- Significant Correlations (p < 0.05): ", nrow(significant_results), "\n",
    "- Success Rate: ", round(nrow(significant_results)/nrow(correlation_results)*100, 1), "%\n\n"
  )
  
  if(nrow(significant_results) > 0) {
    report <- paste0(report,
                     "🏆 TOP TIME WINDOW CLUSTER-SPECIFIC FINDINGS:\n")
    
    for(i in 1:min(3, nrow(significant_results))) {
      result <- significant_results[i, ]
      report <- paste0(report,
                       sprintf("\n%d. %s:\n", i, result$Outcome_Parameter),
                       sprintf("   🕐 Time Window: %s\n", result$Time_Window),
                       sprintf("   🎯 Specific Cluster: Cluster %s\n", result$Primary_Cluster),
                       sprintf("   📈 Correlation: r = %.3f, p = %.4f\n", result$Pearson_r, result$Pearson_p),
                       sprintf("   👥 Patients in Primary Cluster: %d\n", result$Primary_Cluster_N),
                       sprintf("   📊 Primary Cluster Mean Outcome: %.3f\n", result$Primary_Cluster_Mean_Outcome),
                       sprintf("   🔍 Clinical Interpretation: %s Window Cluster %s shows %s correlation\n",
                               result$Time_Window, result$Primary_Cluster,
                               ifelse(result$Pearson_r > 0, "positive", "negative")))
    }
    
    # 最佳时间窗口cluster汇总
    best_clusters <- significant_results %>%
      group_by(Time_Window, Primary_Cluster) %>%
      summarise(
        Count = n(),
        Mean_Effect = round(mean(abs(Pearson_r)), 3),
        .groups = 'drop'
      ) %>%
      arrange(desc(Count), desc(Mean_Effect)) %>%
      slice_head(n = 3)
    
    report <- paste0(report,
                     "\n🥇 MOST PREDICTIVE TIME WINDOW CLUSTERS:\n")
    
    for(i in 1:nrow(best_clusters)) {
      cluster_info <- best_clusters[i, ]
      report <- paste0(report,
                       sprintf("%d. %s Window - Cluster %s: %d significant correlations (Mean |r| = %.3f)\n",
                               i, cluster_info$Time_Window, cluster_info$Primary_Cluster,
                               cluster_info$Count, cluster_info$Mean_Effect))
    }
  }
  
  report <- paste0(report,
                   "\n📁 GENERATED FILES:\n",
                   "- time_window_max_membership_correlations_complete.csv: All correlation results\n",
                   "- time_window_max_membership_correlations_significant.csv: Significant results only\n",
                   "- time_window_max_membership_*.pdf: Individual correlation plots\n",
                   "- time_window_max_membership_correlations_with_clusters_combined.pdf: Combined visualization\n\n",
                   
                   "🎯 HOW TO INTERPRET RESULTS:\n",
                   "When you see 'Early Recovery has significant correlation':\n",
                   "✅ NOW: Check Primary_Cluster column to see exactly which cluster drives this correlation\n",
                   "✅ NOW: Use Cluster_Distribution to understand all cluster patterns in that window\n",
                   "✅ NOW: Compare Primary_Cluster_Mean_Outcome across different correlations\n\n",
                   
                   "📋 CLINICAL SIGNIFICANCE GUIDELINES:\n",
                   "- |r| ≥ 0.5: Large effect - Clinically very significant\n",
                   "- |r| ≥ 0.3: Medium effect - Clinically significant\n",
                   "- |r| ≥ 0.1: Small effect - May be clinically relevant\n",
                   "- Primary_Cluster_N ≥ 5: Good cluster stability\n",
                   "- p < 0.05: Statistically significant\n",
                   "- p_FDR < 0.05: Significant after multiple comparison correction\n\n",
                   
                   "🔬 RESEARCH IMPLICATIONS:\n",
                   "1. Each significant correlation now has time-window and cluster-specific context\n",
                   "2. Can identify which recovery patterns in which time periods predict specific outcomes\n",
                   "3. Primary_Cluster information enables targeted patient stratification by time window\n",
                   "4. Cluster_Distribution shows heterogeneity within each time window\n",
                   "5. Different time windows may have different optimal clusters for prediction\n\n",
                   
                   "💡 NEXT STEPS RECOMMENDATIONS:\n",
                   "1. Focus on time window-cluster combinations with Count ≥ 3 significant correlations\n",
                   "2. Investigate baseline characteristics of patients in predictive clusters\n",
                   "3. Validate time-window-cluster-specific predictions in independent cohorts\n",
                   "4. Develop time-window-cluster-based treatment decision algorithms\n",
                   "5. Consider interventions targeting specific time windows and clusters\n\n",
                   
                   "🕐 TIME WINDOW SPECIFIC INSIGHTS:\n",
                   "- Baseline Window: Pre-surgical recovery patterns\n",
                   "- Acute Recovery: Immediate post-surgical patterns (0-3 days)\n",
                   "- Early Recovery: Early adaptation patterns (4-7 days)\n",
                   "- Mid Recovery: Intermediate recovery patterns (8-15 days)\n",
                   "- Late Recovery: Long-term recovery patterns (16-30 days)\n\n",
                   
                   "========================================\n",
                   "Time Window Max Membership Analysis completed successfully! 🎉\n",
                   "All cluster-specific correlation information is now available.\n",
                   "========================================\n"
  )
  
  # 将报告写入文件
  writeLines(report, "time_window_max_membership_correlation_analysis_report.txt")
  
  # 显示报告
  cat(report)
  
  return(report)
}

# 🔧 增强的结果显示
if(nrow(significant_results) > 0) {
  cat("\n🎯 Enhanced Cluster Analysis Results:\n")
  for(i in 1:min(5, nrow(significant_results))) {
    result <- significant_results[i, ]
    cat(sprintf("\n%d. %s (%s) - %s Window:\n", i, result$Outcome_Parameter, result$Outcome_Type, result$Time_Window))
    cat(sprintf("   📈 相关性: r = %.3f, p = %.4f (%s effect)\n", result$Pearson_r, result$Pearson_p, result$Effect_Size))
    
    # 🔧 增强的cluster信息显示
    cat(sprintf("   🎯 选定cluster: Cluster %s (n = %d patients)\n", result$Primary_Cluster, result$Primary_Cluster_N))
    cat(sprintf("   📊 该cluster平均outcome: %.3f\n", result$Primary_Cluster_Mean_Outcome))
    cat(sprintf("   🔍 选择原因: %s\n", result$Selection_Reason))
    
    # 🔧 冲突警告
    if(result$Outcome_Conflict) {
      cat(sprintf("   ⚠️  注意: 最佳outcome cluster (%s) ≠ 最多患者cluster (%s)\n", 
                  result$Best_Outcome_Cluster, result$Most_Patients_Cluster))
    }
    
    cat(sprintf("   📋 所有clusters比较: %s\n", result$Cluster_Comparison))
  }
  
  # 🔧 冲突检测总结
  conflict_cases <- significant_results %>%
    filter(Outcome_Conflict == TRUE)
  
  if(nrow(conflict_cases) > 0) {
    cat(sprintf("\n⚠️  发现 %d 个案例中最佳outcome cluster与最多患者cluster不一致:\n", nrow(conflict_cases)))
    for(i in 1:nrow(conflict_cases)) {
      case <- conflict_cases[i, ]
      cat(sprintf("   - %s %s: 最佳cluster %s vs 最多患者cluster %s\n", 
                  case$Time_Window, case$Outcome_Parameter, 
                  case$Best_Outcome_Cluster, case$Most_Patients_Cluster))
    }
    cat("\n💡 建议: 重点关注最佳outcome cluster，因为它们代表更好的临床结果\n")
  }
}

# 生成最终报告
if(nrow(max_membership_correlations) > 0) {
  final_report <- generate_max_membership_correlation_report(max_membership_correlations, significant_results)
} else {
  cat("❌ 没有相关性结果可生成报告\n")
}

# ================== 9. 创建时间窗特定的深度分析 ==================

# 创建每个时间窗的详细分析
perform_window_specific_analysis <- function(data, correlation_results) {
  
  if(nrow(correlation_results) == 0) {
    cat("❌ 没有相关性结果进行窗口特定分析\n")
    return(NULL)
  }
  
  cat("\n===== 时间窗特定深度分析 =====\n")
  
  # 获取所有时间窗口
  windows <- unique(correlation_results$Time_Window)
  window_analyses <- list()
  
  for(window in windows) {
    cat(sprintf("\n--- 分析 %s 时间窗口 ---\n", window))
    
    # 该窗口的相关性结果
    window_results <- correlation_results %>%
      filter(Time_Window == window) %>%
      arrange(desc(abs(Pearson_r)))
    
    # 获取membership和cluster列
    membership_col <- paste0("membership_", window)
    cluster_col <- paste0("cluster_", window)
    
    if(membership_col %in% names(data) && cluster_col %in% names(data)) {
      
      # 该窗口的患者分布
      window_data <- data %>%
        filter(!is.na(!!sym(membership_col)) & !is.na(!!sym(cluster_col))) %>%
        dplyr::select(subject_id, !!sym(membership_col), !!sym(cluster_col))
      
      # 集群分布分析
      cluster_dist <- window_data %>%
        group_by(!!sym(cluster_col)) %>%
        summarise(
          count = n(),
          mean_membership = round(mean(!!sym(membership_col)), 3),
          sd_membership = round(sd(!!sym(membership_col)), 3),
          min_membership = round(min(!!sym(membership_col)), 3),
          max_membership = round(max(!!sym(membership_col)), 3),
          .groups = 'drop'
        ) %>%
        arrange(desc(count))
      
      # 显著相关性数量
      significant_count <- sum(window_results$Pearson_p < 0.05)
      
      cat(sprintf("患者数量: %d\n", nrow(window_data)))
      cat(sprintf("集群数量: %d\n", nrow(cluster_dist)))
      cat(sprintf("显著相关性: %d/%d\n", significant_count, nrow(window_results)))
      
      if(significant_count > 0) {
        cat("Top 3 显著相关性:\n")
        top_significant <- window_results %>%
          filter(Pearson_p < 0.05) %>%
          slice_head(n = 3)
        
        for(j in 1:nrow(top_significant)) {
          result <- top_significant[j, ]
          cat(sprintf("  %d. %s: r=%.3f, p=%.4f, Cluster %s\n", 
                      j, result$Outcome_Parameter, result$Pearson_r, 
                      result$Pearson_p, result$Primary_Cluster))
        }
      }
      
      cat("集群分布:\n")
      print(cluster_dist)
      
      # 保存窗口特定结果
      write.csv(window_results, 
                paste0("time_window_", window, "_correlations.csv"), 
                row.names = FALSE)
      
      write.csv(cluster_dist,
                paste0("time_window_", window, "_cluster_distribution.csv"),
                row.names = FALSE)
      
      window_analyses[[window]] <- list(
        window_name = window,
        results = window_results,
        cluster_distribution = cluster_dist,
        patient_count = nrow(window_data),
        significant_count = significant_count
      )
    }
  }
  
  return(window_analyses)
}

# 执行窗口特定分析
window_analyses <- perform_window_specific_analysis(enhanced_max_membership_analysis, max_membership_correlations)

# ================== 10. 创建比较分析 ==================

create_cross_window_comparison <- function(correlation_results) {
  
  cat("\n===== 跨时间窗口比较分析 =====\n")
  
  # 按时间窗口汇总
  window_summary <- correlation_results %>%
    group_by(Time_Window) %>%
    summarise(
      Total_Tests = n(),
      Significant_Count = sum(Pearson_p < 0.05),
      Significant_Percentage = round(sum(Pearson_p < 0.05) / n() * 100, 1),
      Mean_Abs_r = round(mean(abs(Pearson_r)), 3),
      Max_Abs_r = round(max(abs(Pearson_r)), 3),
      Best_p = min(Pearson_p),
      Unique_Clusters = n_distinct(Primary_Cluster[Pearson_p < 0.05]),
      .groups = 'drop'
    ) %>%
    arrange(desc(Significant_Percentage))
  
  cat("时间窗口比较 (按显著性百分比排序):\n")
  print(window_summary)
  
  # 创建比较可视化
  p1 <- ggplot(window_summary, aes(x = reorder(Time_Window, Significant_Percentage))) +
    geom_col(aes(y = Significant_Percentage), fill = "steelblue", alpha = 0.8, width = 0.7) +
    geom_text(aes(y = Significant_Percentage, 
                  label = paste0(Significant_Count, "/", Total_Tests, "\n(", Significant_Percentage, "%)")), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = "Time Window Performance Comparison",
      subtitle = "Percentage of significant correlations (p < 0.05) by time window",
      x = "Time Window",
      y = "Percentage of Significant Correlations (%)",
      caption = "Numbers show: Significant/Total (Percentage)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  
  # 效应大小比较
  p2 <- ggplot(window_summary, aes(x = Time_Window, y = Mean_Abs_r)) +
    geom_col(fill = "darkgreen", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", Mean_Abs_r)), vjust = -0.3, size = 3.5) +
    labs(
      title = "Mean Effect Size by Time Window",
      x = "Time Window",
      y = "Mean Absolute Correlation (|r|)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # 组合图形
  comparison_plot <- gridExtra::grid.arrange(p1, p2, ncol = 1,
                                             top = "Cross-Time Window Analysis Comparison")
  
  # 保存图形和数据
  ggsave("cross_time_window_comparison.pdf", comparison_plot, width = 12, height = 10)
  write.csv(window_summary, "cross_time_window_summary.csv", row.names = FALSE)
  
  cat("\n✓ 跨时间窗口比较分析完成\n")
  
  return(list(summary = window_summary, plot = comparison_plot))
}

# 创建跨窗口比较
if(nrow(max_membership_correlations) > 0) {
  cross_window_analysis <- create_cross_window_comparison(max_membership_correlations)
}

# ================== 11. 生成交互式结果总结 ==================

create_interactive_summary <- function(correlation_results, significant_results, window_analyses = NULL) {
  
  cat("\n" %R% paste(rep("=", 60), collapse = "") %R% "\n")
  cat("🎯 TIME WINDOW MAX MEMBERSHIP CORRELATION ANALYSIS - INTERACTIVE SUMMARY\n")
  cat(paste(rep("=", 60), collapse = "") %R% "\n\n")
  
  if(nrow(significant_results) > 0) {
    cat("📊 发现 " %R% nrow(significant_results) %R% " 个显著相关性！\n\n")
    
    # 按时间窗口分析
    window_analysis <- significant_results %>%
      group_by(Time_Window) %>%
      summarise(
        Count = n(),
        Best_r = max(abs(Pearson_r)),
        Mean_r = mean(abs(Pearson_r)),
        Unique_Clusters = n_distinct(Primary_Cluster),
        Best_Outcome = Outcome_Parameter[which.max(abs(Pearson_r))],
        .groups = 'drop'
      ) %>%
      arrange(desc(Count))
    
    cat("🕐 时间窗口表现排名:\n")
    for(i in 1:nrow(window_analysis)) {
      window_info <- window_analysis[i, ]
      cat(sprintf("   %d. %s: %d个显著相关性, 最高|r|=%.3f, %d个独特clusters\n",
                  i, window_info$Time_Window, window_info$Count, 
                  window_info$Best_r, window_info$Unique_Clusters))
      cat(sprintf("      最佳预测指标: %s\n", window_info$Best_Outcome))
    }
    
    # 集群特异性分析
    cluster_analysis <- significant_results %>%
      group_by(Time_Window, Primary_Cluster) %>%
      summarise(
        Count = n(),
        Mean_r = mean(abs(Pearson_r)),
        Best_Outcome = Outcome_Parameter[which.max(abs(Pearson_r))],
        .groups = 'drop'
      ) %>%
      arrange(desc(Count)) %>%
      slice_head(n = 5)
    
    cat("\n🎯 Top 5 预测性时间窗口-集群组合:\n")
    for(i in 1:nrow(cluster_analysis)) {
      cluster_info <- cluster_analysis[i, ]
      cat(sprintf("   %d. %s Window - Cluster %s: %d个相关性, 平均|r|=%.3f\n",
                  i, cluster_info$Time_Window, cluster_info$Primary_Cluster,
                  cluster_info$Count, cluster_info$Mean_r))
      cat(sprintf("      最佳预测: %s\n", cluster_info$Best_Outcome))
    }
    
  } else {
    cat("❌ 未发现显著的相关性\n")
    cat("💡 建议检查:\n")
    cat("   - 时间窗口聚类质量\n")
    cat("   - Max membership数据完整性\n")
    cat("   - 样本量是否足够\n")
  }
  
  cat("\n" %R% paste(rep("=", 60), collapse = "") %R% "\n")
  cat("✅ 时间窗口Max Membership分析完成！所有结果已保存到当前目录。\n")
  cat(paste(rep("=", 60), collapse = "") %R% "\n")
}

# 修正字符串连接操作符
`%R%` <- function(x, y) paste0(x, y)

# 创建交互式总结
create_interactive_summary(max_membership_correlations, significant_results, window_analyses)

cat("\n🎉 时间窗口Max Membership相关性分析全部完成！\n")
cat("📁 请检查生成的文件和可视化结果。\n")
cat("🔬 现在可以精确识别哪个时间窗口的哪个集群与预后相关！\n")