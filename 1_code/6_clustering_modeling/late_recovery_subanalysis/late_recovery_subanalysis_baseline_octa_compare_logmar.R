# Late Recovery时间窗口聚类的术前基线特征差异分析 - logMAR版本
# 分析目标：基于late recovery聚类结果，分析患者群体在术前的基线差异（使用logMAR视力）

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(r4projects)

setwd(get_project_wd())
rm(list = ls())

# ================== 1. 加载数据和聚类结果 ==================

cat("===== Late Recovery时间窗口聚类的术前基线特征差异分析 (logMAR版本) =====\n")
cat("分析目标：确定late recovery阶段的患者群体差异是否在术前就存在（使用logMAR视力）\n\n")

# 加载baseline信息和OCTA数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# 🎯 关键修改：加载late recovery聚类结果
late_recovery_results_file <- "3_data_analysis/6_clustering_modeling/time_window_clustering/late_recovery_detailed_2cluster_membership.csv"

if(file.exists(late_recovery_results_file)) {
  late_recovery_clusters <- read.csv(late_recovery_results_file, stringsAsFactors = FALSE)
  cat("✓ 成功加载late recovery聚类结果:", nrow(late_recovery_clusters), "患者\n")
  
  # 检查数据结构
  cat("Late recovery聚类数据列名:", paste(names(late_recovery_clusters), collapse = ", "), "\n")
  cat("前3行数据:\n")
  print(head(late_recovery_clusters, 3))
  
} else {
  stop("请先运行late recovery时间窗口聚类分析！文件路径:", late_recovery_results_file)
}

# 创建输出目录
dir.create("3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/baseline_analysis", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/baseline_analysis")

# ================== 2. 删除不需要的复杂转换函数 ==================

# 简化：直接使用小数视力转logMAR的简单公式
# logMAR = -log10(decimal_acuity)
# 例如：小数视力1.0 → logMAR = -log10(1.0) = 0
#      小数视力0.5 → logMAR = -log10(0.5) = 0.301
#      小数视力0.1 → logMAR = -log10(0.1) = 1.0

# ================== 3. 处理术前OCTA数据（保持不变）==================

# 重用原代码的数据处理函数
process_octa_data <- function(baseline_data, octa_data, id_column = "id") {
  features <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  return(features)
}

process_patient_data <- function(patient_data, time_points = c("T0")) {  # 只要T0
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

process_all_patients <- function(features_data) {
  patient_list <- split(features_data, features_data$ID)
  processed_data <- purrr::map(patient_list, process_patient_data)
  return(bind_rows(processed_data))
}

# 处理OCTA数据
octa_bloodflow_features <- process_octa_data(baseline_info, octa_bloodflow)
octa_thickness_features <- process_octa_data(baseline_info, octa_thickness)

octa_bloodflow_t0 <- process_all_patients(octa_bloodflow_features)
octa_thickness_t0 <- process_all_patients(octa_thickness_features)

# 获取PPV患者
ppv_patients <- baseline_info %>%
  filter(surgery_1..0.PI.1.other. == 1) %>%
  distinct(ID) %>%
  pull(ID)

ppv_bloodflow_t0 <- octa_bloodflow_t0 %>% filter(ID %in% ppv_patients)
ppv_thickness_t0 <- octa_thickness_t0 %>% filter(ID %in% ppv_patients)

# ================== 4. 提取术前参数（保持0_21区域专注）==================

# 筛选T0参数 - 专注于0_21区域
filter_baseline_bloodflow <- function(data) {
  layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
  regions_of_interest <- c("0_21")  # 只关注0_21区域
  
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("未找到0_21区域T0血流参数！")
    return(list(data = data %>% dplyr::select(ID), params = character(0)))
  }
  
  cat("找到", length(params_T0), "个0_21区域T0血流参数\n")
  
  filtered_data <- data %>%
    dplyr::select(ID, all_of(params_T0))
  
  return(list(
    data = filtered_data,
    params_T0 = params_T0,
    target_region = "0_21"
  ))
}

filter_baseline_thickness <- function(data) {
  layers_of_interest <- c("GCL.IPL", "INL", "Retina")
  regions_of_interest <- c("0_21")  # 只关注0_21区域
  
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("未找到0_21区域T0厚度参数！")
    return(list(data = data %>% dplyr::select(ID), params = character(0)))
  }
  
  cat("找到", length(params_T0), "个0_21区域T0厚度参数\n")
  
  filtered_data <- data %>%
    dplyr::select(ID, all_of(params_T0))
  
  return(list(
    data = filtered_data,
    params_T0 = params_T0,
    target_region = "0_21"
  ))
}

# 提取T0参数
baseline_bloodflow_filtered <- filter_baseline_bloodflow(ppv_bloodflow_t0)
baseline_thickness_filtered <- filter_baseline_thickness(ppv_thickness_t0)

# ================== 5. 创建术前基线数据集（修改为logMAR视力）==================

# 🎯 关键修改：处理术前视力数据并转换为logMAR
baseline_vision <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    # 根据手术眼确定对应的视力值（原始小数视力）
    pre_vision_decimal = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,
      surgery_eye_1 == 1 ~ os_corrected_bas,
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,
      TRUE ~ NA_real_
    ),
    # 🎯 核心修改：将小数视力转换为logMAR视力
    pre_vision = -log10(pre_vision_decimal)  # logMAR = -log10(decimal_acuity)
  ) %>%
  dplyr::select(ID, pre_vision_decimal, pre_vision, age, gender) %>%
  filter(ID %in% ppv_patients)

# 检查logMAR转换结果
cat("\n===== 小数视力转logMAR转换结果检查 =====\n")
cat("原始小数视力范围:", sprintf("%.3f - %.3f", 
                         min(baseline_vision$pre_vision_decimal, na.rm = TRUE),
                         max(baseline_vision$pre_vision_decimal, na.rm = TRUE)), "\n")
cat("转换后logMAR范围:", sprintf("%.3f - %.3f", 
                            min(baseline_vision$pre_vision, na.rm = TRUE),
                            max(baseline_vision$pre_vision, na.rm = TRUE)), "\n")
cat("成功转换患者数:", sum(!is.na(baseline_vision$pre_vision)), "/", nrow(baseline_vision), "\n")

# 显示一些转换示例
conversion_examples <- baseline_vision %>%
  filter(!is.na(pre_vision_decimal) & !is.na(pre_vision)) %>%
  mutate(
    snellen_equiv = paste0("20/", round(pre_vision_decimal * 20)),
    logmar_rounded = round(pre_vision, 3)
  ) %>%
  head(10) %>%
  dplyr::select(ID, pre_vision_decimal, snellen_equiv, logmar_rounded)
cat("转换示例:\n")
print(conversion_examples)

# 统一ID列名
if("subject_id" %in% names(late_recovery_clusters)) {
  # 已经是subject_id，不需要修改
} else if("ID" %in% names(late_recovery_clusters)) {
  names(late_recovery_clusters)[names(late_recovery_clusters) == "ID"] <- "subject_id"
}

# 🎯 合并所有术前数据（现在pre_vision已经是logMAR格式）
baseline_comprehensive <- baseline_vision %>%
  # 只保留需要的列（pre_vision现在是logMAR）
  dplyr::select(ID, pre_vision, age, gender) %>%  # pre_vision_decimal可选择保留或移除
  full_join(baseline_bloodflow_filtered$data, by = "ID") %>%
  full_join(baseline_thickness_filtered$data, by = "ID") %>%
  # 添加late recovery聚类信息
  inner_join(late_recovery_clusters %>% 
               dplyr::select(subject_id, max_cluster, max_membership), 
             by = c("ID" = "subject_id"))

cat("\n===== 术前基线数据集摘要（Late Recovery聚类 - logMAR版本）=====\n")
cat("专注分析：0_21区域（广角区域）+ logMAR视力参数\n")
cat("聚类来源：Late Recovery时间窗口\n")
cat("视力转换：小数视力 → logMAR视力（-log10转换）\n")
cat("患者数量:", nrow(baseline_comprehensive), "\n")
cat("总参数数:", ncol(baseline_comprehensive) - 3, "\n")  # 排除ID, cluster, membership
cat("- logMAR视力/基本信息:", 3, "\n")  # pre_vision(logMAR), age, gender
cat("- OCTA血流 (0_21区域):", length(baseline_bloodflow_filtered$params_T0), "\n")
cat("- OCTA厚度 (0_21区域):", length(baseline_thickness_filtered$params_T0), "\n")

# 检查聚类分布
cat("\nLate Recovery聚类分布:\n")
cluster_distribution <- table(baseline_comprehensive$max_cluster)
print(cluster_distribution)

# 检查数据完整性
baseline_complete_cases <- baseline_comprehensive[complete.cases(baseline_comprehensive), ]
cat("完整数据患者:", nrow(baseline_complete_cases), 
    "(", round(nrow(baseline_complete_cases)/nrow(baseline_comprehensive)*100, 1), "%)\n")

# 检查logMAR视力分布
cat("\nlogMAR视力描述性统计:\n")
logmar_stats <- baseline_comprehensive %>%
  group_by(max_cluster) %>%
  summarise(
    n = n(),
    mean_logmar = mean(pre_vision, na.rm = TRUE),
    sd_logmar = sd(pre_vision, na.rm = TRUE),
    median_logmar = median(pre_vision, na.rm = TRUE),
    min_logmar = min(pre_vision, na.rm = TRUE),
    max_logmar = max(pre_vision, na.rm = TRUE),
    # 添加Snellen等效值
    mean_snellen_equiv = paste0("20/", round(10^mean_logmar * 20)),
    median_snellen_equiv = paste0("20/", round(10^median_logmar * 20)),
    .groups = 'drop'
  )
print(logmar_stats)

cat("\nlogMAR转换说明:\n")
cat("- logMAR = 0: 相当于20/20视力（正常）\n")
cat("- logMAR = 0.3: 相当于20/40视力（轻度下降）\n") 
cat("- logMAR = 0.7: 相当于20/100视力（中度下降）\n")
cat("- logMAR = 1.0: 相当于20/200视力（重度下降）\n")
cat("- logMAR值越小，视力越好\n")

# ================== 6. 术前基线差异统计分析（保持原函数，但现在使用logMAR）==================

analyze_baseline_differences_late_recovery <- function(data) {
  
  cat("\n===== 术前基线差异统计分析（Late Recovery聚类 - logMAR视力）=====\n")
  
  # 确定要分析的参数
  analysis_params <- names(data)[!names(data) %in% c("ID", "max_cluster", "max_membership")]
  
  results <- data.frame(
    Parameter = character(),
    Data_Type = character(),
    Region = character(),
    Variable_Type = character(),
    Cluster1_Mean = numeric(),
    Cluster2_Mean = numeric(),
    Mean_Difference = numeric(),
    P_Value = numeric(),
    Effect_Size = numeric(),
    Test_Method = character(),
    Significant = character(),
    stringsAsFactors = FALSE
  )
  
  for(param in analysis_params) {
    if(param %in% names(data) && !all(is.na(data[[param]]))) {
      
      # 确定参数类型和区域
      data_type <- case_when(
        param %in% c("pre_vision", "age", "gender") ~ "Baseline",
        grepl("SVP|ICP|DCP|Choroid", param) ~ "Blood Flow",
        grepl("GCL|INL|Retina", param) ~ "Thickness",
        TRUE ~ "Other"
      )
      
      region <- case_when(
        param == "pre_vision" ~ "logMAR Vision",
        param %in% c("age", "gender") ~ "Demographics",
        grepl("0_21", param) ~ "Wide-field (0_21)",
        TRUE ~ "Unknown"
      )
      
      # 判断变量类型
      param_values <- data[[param]][!is.na(data[[param]])]
      unique_values <- unique(param_values)
      
      # 判断是否为二分类变量
      is_binary <- length(unique_values) == 2 && all(unique_values %in% c(0, 1))
      
      variable_type <- if(is_binary) "Binary" else "Continuous"
      
      # 使用max_cluster进行分组分析
      param_data <- data[, c("max_cluster", param)]
      param_data <- param_data[!is.na(param_data[[param]]), ]
      
      if(nrow(param_data) > 0 && length(unique(param_data$max_cluster)) >= 2) {
        
        # 计算均值（或比例）
        means <- tapply(param_data[[param]], param_data$max_cluster, mean, na.rm = TRUE)
        
        # 根据变量类型选择检验方法
        if(is_binary) {
          # 二分类变量：使用Fisher精确检验
          cat(sprintf("对二分类变量 %s 使用Fisher精确检验\n", param))
          
          # 创建列联表
          contingency_table <- table(param_data$max_cluster, param_data[[param]])
          
          # Fisher精确检验
          test_result <- try(fisher.test(contingency_table), silent = TRUE)
          test_method <- "Fisher's Exact Test"
          
          # 计算效应量（Cramér's V或φ系数）
          if(class(test_result) != "try-error") {
            # 对于2x2表，φ系数 = √(χ²/n)
            chi_sq <- chisq.test(contingency_table, correct = FALSE)$statistic
            effect_size <- sqrt(as.numeric(chi_sq) / nrow(param_data))
          } else {
            effect_size <- NA
          }
          
        } else {
          # 连续变量：使用t检验（包括logMAR视力）
          cat(sprintf("对连续变量 %s 使用t检验\n", param))
          if(param == "pre_vision") {
            cat("  -> 分析logMAR视力差异\n")
          }
          
          test_result <- try(t.test(reformulate("max_cluster", param), data = param_data), silent = TRUE)
          test_method <- "Independent t-test"
          
          # 计算Cohen's d
          if(class(test_result) != "try-error") {
            cluster1_data <- param_data[param_data$max_cluster == 1, param]
            cluster2_data <- param_data[param_data$max_cluster == 2, param]
            
            pooled_sd <- sqrt(((length(cluster1_data) - 1) * var(cluster1_data, na.rm = TRUE) + 
                                 (length(cluster2_data) - 1) * var(cluster2_data, na.rm = TRUE)) / 
                                (length(cluster1_data) + length(cluster2_data) - 2))
            
            effect_size <- abs(mean(cluster2_data, na.rm = TRUE) - mean(cluster1_data, na.rm = TRUE)) / pooled_sd
          } else {
            effect_size <- NA
          }
        }
        
        # 保存结果
        if(class(test_result) != "try-error") {
          results <- rbind(results, data.frame(
            Parameter = gsub("_T0$", "", param),
            Data_Type = data_type,
            Region = region,
            Variable_Type = variable_type,
            Cluster1_Mean = ifelse("1" %in% names(means), means["1"], NA),
            Cluster2_Mean = ifelse("2" %in% names(means), means["2"], NA),
            Mean_Difference = ifelse(length(means) >= 2, means["2"] - means["1"], NA),
            P_Value = test_result$p.value,
            Effect_Size = effect_size,
            Test_Method = test_method,
            Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # 多重比较校正
  if(nrow(results) > 0) {
    results$P_Adjusted <- p.adjust(results$P_Value, method = "fdr")
    results$Significant_Adjusted <- ifelse(results$P_Adjusted < 0.05, "Yes", "No")
    results <- results %>% arrange(Data_Type, Region, P_Value)
  }
  
  return(results)
}

# ================== 7. 可视化函数（修改为logMAR版本）==================

# 修改可视化函数，针对logMAR视力优化
create_baseline_visualizations_late_recovery_logmar <- function(data, stats_results) {
  
  cat("\n===== 创建Late Recovery基线可视化 (logMAR版本) =====\n")
  
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  
  # 在函数开始处定义 analysis_params
  analysis_params <- names(data)[!names(data) %in% c("ID", "max_cluster", "max_membership")]
  
  # 获取显著差异的参数（使用矫正后的p值）
  significant_baseline <- stats_results %>% 
    filter(Significant_Adjusted == "Yes") %>%
    arrange(P_Adjusted)
  
  # 为二分类变量创建堆叠条形图
  create_binary_plot_lr <- function(data, param, param_clean, p_value, test_method) {
    
    # 计算比例
    prop_data <- data %>%
      group_by(max_cluster, !!sym(param)) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(max_cluster) %>%
      mutate(
        total = sum(count),
        percentage = round(count / total * 100, 1)
      ) %>%
      ungroup()
    
    # 为gender变量添加标签
    if(param == "gender") {
      prop_data <- prop_data %>%
        mutate(
          gender_label = case_when(
            !!sym(param) == 0 ~ "Female",
            !!sym(param) == 1 ~ "Male",
            TRUE ~ as.character(!!sym(param))
          )
        )
      fill_var <- "gender_label"
    } else {
      prop_data$category <- factor(prop_data[[param]])
      fill_var <- "category"
    }
    
    # 格式化p值
    p_text <- if(is.na(p_value)) {
      "p = N/A"
    } else if(p_value < 0.001) {
      "p(adj) < 0.001"
    } else {
      paste("p(adj) =", round(p_value, 3))
    }
    
    # 创建堆叠条形图
    p_binary <- ggplot(prop_data, aes(x = as.factor(max_cluster), y = percentage, fill = !!sym(fill_var))) +
      geom_col(position = "stack", alpha = 0.8) +
      geom_text(aes(label = paste0(count, "\n(", percentage, "%)")), 
                position = position_stack(vjust = 0.5), 
                color = "white", fontface = "bold", size = 3) +
      scale_fill_manual(
        values = if(param == "gender") c("Female" = "#FF69B4", "Male" = "#4169E1") else c("#a388bf", "#bc982f"),
        name = if(param == "gender") "Gender" else param_clean
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      ) +
      labs(
        title = paste("Pre-operative Baseline:", param_clean, "Distribution"),
        subtitle = paste("Late Recovery Cluster Comparison |", p_text),
        x = "Late Recovery Cluster",
        y = "Percentage",
        caption = paste("Statistical test:", test_method, "| Numbers show count and percentage")
      )
    
    # 保存图片
    ggsave(paste0("plots/baseline_", param, "_distribution_late_recovery_logmar.pdf"), 
           p_binary, width = 8, height = 6)
    ggsave(paste0("plots/baseline_", param, "_distribution_late_recovery_logmar.png"), 
           p_binary, width = 8, height = 6, dpi = 300)
  }
  
  # 🎯 特别优化：为logMAR视力创建箱线图
  create_logmar_vision_plot <- function(data, param, p_value, test_method) {
    
    # 格式化p值
    p_text <- if(is.na(p_value)) {
      "p = N/A"
    } else if(p_value < 0.001) {
      "p(adj) < 0.001"
    } else {
      paste("p(adj) =", round(p_value, 3))
    }
    
    # 计算描述性统计
    vision_stats <- data %>%
      group_by(max_cluster) %>%
      summarise(
        n = sum(!is.na(pre_vision)),
        mean_logmar = mean(pre_vision, na.rm = TRUE),
        median_logmar = median(pre_vision, na.rm = TRUE),
        sd_logmar = sd(pre_vision, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # 创建logMAR视力的特殊箱线图
    p_logmar <- ggplot(data, aes(x = as.factor(max_cluster), y = pre_vision, 
                                 fill = as.factor(max_cluster))) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_manual(
        values = c("1" = "#a388bf", "2" = "#bc982f"),
        name = "Late Recovery\nCluster"
      ) +
      # 添加参考线（正常视力 logMAR = 0）
      geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
      annotate("text", x = 0.5, y = 0.05, label = "Normal Vision (20/20)", 
               color = "red", size = 3, hjust = 0) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      )
    
    # 添加p值标注
    if(!is.na(p_value)) {
      y_max <- max(data$pre_vision, na.rm = TRUE)
      y_min <- min(data$pre_vision, na.rm = TRUE)
      y_range <- y_max - y_min
      y_pos <- y_max + 0.1 * y_range
      
      p_logmar <- p_logmar +
        geom_segment(aes(x = 1, xend = 2, y = y_pos, yend = y_pos), 
                     color = "black", inherit.aes = FALSE) +
        geom_segment(aes(x = 1, xend = 1, y = y_pos, yend = y_pos - 0.02 * y_range), 
                     color = "black", inherit.aes = FALSE) +
        geom_segment(aes(x = 2, xend = 2, y = y_pos, yend = y_pos - 0.02 * y_range), 
                     color = "black", inherit.aes = FALSE) +
        annotate("text", x = 1.5, y = y_pos + 0.03 * y_range, 
                 label = p_text, hjust = 0.5, size = 4, fontface = "bold",
                 color = ifelse(!is.na(p_value) && p_value < 0.05, "red", "black"))
    }
    
    # 添加标签和说明
    p_logmar <- p_logmar +
      labs(
        title = "Pre-operative logMAR Visual Acuity",
        subtitle = paste("Late Recovery Cluster Comparison |", p_text),
        x = "Late Recovery Cluster",
        y = "logMAR Visual Acuity (lower = better vision)",
        caption = paste("Statistical test:", test_method, 
                        "| Red dashed line = Normal vision (20/20)", 
                        "| Higher logMAR values indicate worse vision")
      )
    
    # 保存图片
    ggsave("plots/baseline_logmar_vision_boxplot_late_recovery.pdf", 
           p_logmar, width = 10, height = 8)
    ggsave("plots/baseline_logmar_vision_boxplot_late_recovery.png", 
           p_logmar, width = 10, height = 8, dpi = 300)
    
    # 同时生成描述性统计表
    stat_table <- vision_stats %>%
      mutate(
        mean_snellen_equiv = paste0("20/", round(10^mean_logmar * 20)),
        logmar_interpretation = case_when(
          mean_logmar <= 0 ~ "Normal/Better (≤20/20)",
          mean_logmar <= 0.3 ~ "Mild impairment (20/25-20/40)",
          mean_logmar <= 0.7 ~ "Moderate impairment (20/50-20/100)",
          mean_logmar <= 1.0 ~ "Severe impairment (20/125-20/200)",
          TRUE ~ "Profound impairment (>20/200)"
        )
      ) %>%
      dplyr::select(max_cluster, n, mean_logmar, mean_snellen_equiv, logmar_interpretation)
    
    write.csv(stat_table, "logmar_vision_statistics_by_cluster.csv", row.names = FALSE)
    
    cat("✓ logMAR视力专用图表已保存\n")
    cat("✓ 视力描述性统计已保存: logmar_vision_statistics_by_cluster.csv\n")
  }
  
  # 为连续变量创建箱线图
  create_continuous_plot_lr <- function(data, param, param_clean, p_value, test_method) {
    
    # 格式化p值
    p_text <- if(is.na(p_value)) {
      "p = N/A"
    } else if(p_value < 0.001) {
      "p(adj) < 0.001"
    } else {
      paste("p(adj) =", round(p_value, 3))
    }
    
    # 创建箱线图
    p_box <- ggplot(data, aes(x = as.factor(max_cluster), y = !!sym(param), 
                              fill = as.factor(max_cluster))) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      scale_fill_manual(
        values = c("1" = "#a388bf", "2" = "#bc982f"),
        name = "Late Recovery\nCluster"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")
      )
    
    # 添加p值标注
    if(!is.na(p_value)) {
      y_max <- max(data[[param]], na.rm = TRUE)
      y_min <- min(data[[param]], na.rm = TRUE)
      y_range <- y_max - y_min
      y_pos <- y_max + 0.1 * y_range
      
      p_box <- p_box +
        geom_segment(aes(x = 1, xend = 2, y = y_pos, yend = y_pos), 
                     color = "black", inherit.aes = FALSE) +
        geom_segment(aes(x = 1, xend = 1, y = y_pos, yend = y_pos - 0.02 * y_range), 
                     color = "black", inherit.aes = FALSE) +
        geom_segment(aes(x = 2, xend = 2, y = y_pos, yend = y_pos - 0.02 * y_range), 
                     color = "black", inherit.aes = FALSE) +
        annotate("text", x = 1.5, y = y_pos + 0.03 * y_range, 
                 label = p_text, hjust = 0.5, size = 4, fontface = "bold",
                 color = ifelse(!is.na(p_value) && p_value < 0.05, "red", "black"))
    }
    
    # 添加标签
    p_box <- p_box +
      labs(
        title = paste("Pre-operative Baseline:", param_clean),
        subtitle = paste("Late Recovery Cluster Comparison |", p_text),
        x = "Late Recovery Cluster",
        y = param_clean,
        caption = paste("Statistical test:", test_method, "| Individual points show patients")
      )
    
    ggsave(paste0("plots/baseline_", param, "_boxplot_late_recovery_logmar.pdf"), 
           p_box, width = 8, height = 6)
    ggsave(paste0("plots/baseline_", param, "_boxplot_late_recovery_logmar.png"), 
           p_box, width = 8, height = 6, dpi = 300)
  }
  
  # 处理显著差异的参数
  if(nrow(significant_baseline) > 0) {
    
    significant_params <- paste0(significant_baseline$Parameter, 
                                 ifelse(significant_baseline$Data_Type != "Baseline", "_T0", ""))
    available_sig_params <- intersect(significant_params, analysis_params)
    
    if(length(available_sig_params) > 0) {
      
      for(param in available_sig_params) {
        
        param_clean <- gsub("_T0$|_", " ", param)
        param_base <- gsub("_T0$", "", param)
        
        # 获取参数信息
        param_info <- significant_baseline %>% filter(Parameter == param_base)
        
        if(nrow(param_info) > 0) {
          p_value <- param_info$P_Adjusted[1]
          variable_type <- param_info$Variable_Type[1]
          test_method <- param_info$Test_Method[1]
          
          # 特别处理logMAR视力
          if(param == "pre_vision") {
            create_logmar_vision_plot(data, param, p_value, test_method)
          } else if(variable_type == "Binary") {
            create_binary_plot_lr(data, param, param_clean, p_value, test_method)
          } else {
            create_continuous_plot_lr(data, param, param_clean, p_value, test_method)
          }
        }
      }
    }
  }
  
  # 为主要参数创建图（即使不显著）
  main_params <- c("pre_vision", "age", "gender")
  for(param in main_params) {
    if(param %in% names(data)) {
      param_info <- stats_results %>% filter(Parameter == param)
      
      if(nrow(param_info) > 0) {
        p_value <- param_info$P_Adjusted[1]
        variable_type <- param_info$Variable_Type[1]
        test_method <- param_info$Test_Method[1]
        param_clean <- gsub("_", " ", param)
        
        # 特别处理logMAR视力
        if(param == "pre_vision") {
          create_logmar_vision_plot(data, param, p_value, test_method)
        } else if(variable_type == "Binary") {
          create_binary_plot_lr(data, param, param_clean, p_value, test_method)
        } else {
          create_continuous_plot_lr(data, param, param_clean, p_value, test_method)
        }
      }
    }
  }
  
  # OCTA基线特征分组箱线图（保持不变）
  create_grouped_octa_baseline_boxplots <- function(data, analysis_params) {
    
    # 识别OCTA参数
    octa_params <- analysis_params[grepl("SVP|ICP|DCP|Choroid|GCL|INL|Retina", analysis_params) & grepl("_T0$", analysis_params)]
    
    if(length(octa_params) == 0) {
      cat("    No OCTA baseline parameters found\n")
      return(NULL)
    }
    
    cat("    Found", length(octa_params), "OCTA baseline parameters\n")
    
    # 分组OCTA参数
    bloodflow_params <- octa_params[grepl("SVP|ICP|DCP|Choroid", octa_params)]
    vd_params <- octa_params[grepl("VD", octa_params)]
    thickness_params <- octa_params[grepl("GCL|INL|Retina", octa_params)]
    
    # 存储生成的图
    plots_list <- list()
    
    # 创建血流参数箱线图
    if(length(bloodflow_params) > 0) {
      cat("    Creating blood flow parameters boxplot...\n")
      
      bloodflow_data <- data.frame(
        ID = character(0),
        max_cluster = numeric(0),
        Parameter_Name = character(0),
        Parameter_Clean = character(0),
        Baseline_Value = numeric(0),
        Cluster = character(0),
        stringsAsFactors = FALSE
      )
      
      for(param in bloodflow_params) {
        if(param %in% names(data)) {
          param_data <- data %>%
            dplyr::select(ID, max_cluster, all_of(param)) %>%
            filter(!is.na(.data[[param]])) %>%
            mutate(
              Parameter_Name = gsub("_T0$", "", param),
              Parameter_Clean = case_when(
                grepl("Choroid", param) ~ "Choroid",
                grepl("DCP", param) ~ "DCP", 
                grepl("ICP", param) ~ "ICP",
                grepl("SVP", param) ~ "SVP",
                TRUE ~ gsub("_T0$|_0_21", "", param)
              ),
              Baseline_Value = as.numeric(.data[[param]]),
              Cluster = as.character(max_cluster)
            ) %>%
            dplyr::select(ID, max_cluster, Parameter_Name, Parameter_Clean, Baseline_Value, Cluster)
          
          bloodflow_data <- rbind.data.frame(bloodflow_data, param_data)
        }
      }
      
      if(nrow(bloodflow_data) > 0) {
        bloodflow_data$Parameter_Clean <- factor(bloodflow_data$Parameter_Clean, 
                                                 levels = c("Choroid", "DCP", "ICP", "SVP"))
        
        p_bloodflow <- ggplot(bloodflow_data, aes(x = Parameter_Clean, y = Baseline_Value, fill = Cluster)) +
          geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 1) +
          scale_fill_manual(values = c("1" = "#a388bf", "2" = "#bc982f"), 
                            name = "Late Recovery\nCluster") +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            legend.position = "right",
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10)
          ) +
          labs(
            title = "Baseline Blood Flow Parameters (T0)",
            subtitle = "Pre-operative OCTA measurements by Late Recovery Clusters",
            x = "Blood Flow Parameters",
            y = "Baseline Value (T0)"
          )
        
        ggsave("plots/baseline_bloodflow_grouped_boxplot_late_recovery_logmar.pdf", p_bloodflow, width = 8, height = 6)
        ggsave("plots/baseline_bloodflow_grouped_boxplot_late_recovery_logmar.png", p_bloodflow, width = 8, height = 6, dpi = 300)
        
        plots_list[["bloodflow"]] <- p_bloodflow
        cat("    ✓ Baseline blood flow grouped boxplot saved\n")
      }
    }
    
    # 创建厚度参数箱线图
    if(length(thickness_params) > 0) {
      cat("    Creating thickness parameters boxplot...\n")
      
      thickness_data <- data.frame(
        ID = character(0),
        max_cluster = numeric(0),
        Parameter_Clean = character(0),
        Baseline_Value = numeric(0),
        Cluster = character(0),
        stringsAsFactors = FALSE
      )
      
      for(param in thickness_params) {
        if(param %in% names(data)) {
          param_data <- data %>%
            dplyr::select(ID, max_cluster, all_of(param)) %>%
            filter(!is.na(.data[[param]])) %>%
            mutate(
              Parameter_Clean = case_when(
                grepl("GCL.IPL|GCL_IPL", param) ~ "GCL.IPL",
                grepl("INL", param) ~ "INL",
                grepl("OuterRetina", param) ~ "OuterRetina",
                grepl("Retina", param) & !grepl("Outer", param) ~ "Retina",
                TRUE ~ gsub("_T0$|_0_21|_", " ", param)
              ),
              Baseline_Value = as.numeric(.data[[param]]),
              Cluster = as.character(max_cluster)
            ) %>%
            dplyr::select(ID, max_cluster, Parameter_Clean, Baseline_Value, Cluster)
          
          thickness_data <- rbind.data.frame(thickness_data, param_data)
        }
      }
      
      if(nrow(thickness_data) > 0) {
        thickness_data$Parameter_Clean <- factor(thickness_data$Parameter_Clean, 
                                                 levels = c("GCL.IPL", "INL", "OuterRetina", "Retina"))
        
        p_thickness <- ggplot(thickness_data, aes(x = Parameter_Clean, y = Baseline_Value, fill = Cluster)) +
          geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.size = 1) +
          scale_fill_manual(values = c("1" = "#a388bf", "2" = "#bc982f"), 
                            name = "Late Recovery\nCluster") +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            legend.position = "right",
            axis.text.x = element_text(size = 10, angle = 0)
          ) +
          labs(
            title = "Baseline Retinal Thickness Parameters (T0)",
            subtitle = "Pre-operative thickness measurements by Late Recovery Clusters",
            x = "Thickness Parameters",
            y = "Baseline Value (T0, μm)"
          )
        
        ggsave("plots/baseline_thickness_grouped_boxplot_late_recovery_logmar.pdf", p_thickness, width = 9, height = 6)
        ggsave("plots/baseline_thickness_grouped_boxplot_late_recovery_logmar.png", p_thickness, width = 9, height = 6, dpi = 300)
        
        plots_list[["thickness"]] <- p_thickness
        cat("    ✓ Baseline thickness grouped boxplot saved\n")
      }
    }
    
    # 创建组合图（如果有多个图）
    if(length(plots_list) >= 2) {
      cat("    Creating combined plot...\n")
      
      if(!requireNamespace("gridExtra", quietly = TRUE)) {
        install.packages("gridExtra")
      }
      library(gridExtra)
      library(grid)
      
      plots_for_combine <- list()
      plot_names <- names(plots_list)
      
      for(i in 1:length(plots_list)) {
        if(i < length(plots_list)) {
          plots_for_combine[[i]] <- plots_list[[i]] + theme(legend.position = "none")
        } else {
          plots_for_combine[[i]] <- plots_list[[i]] + theme(legend.position = "right")
        }
      }
      
      if(length(plots_for_combine) == 2) {
        combined_plot <- grid.arrange(
          plots_for_combine[[1]], plots_for_combine[[2]],
          ncol = 2,
          top = textGrob("Pre-operative OCTA Parameters by Late Recovery Clusters (logMAR Vision)", 
                         gp = gpar(fontsize = 16, fontface = "bold")),
          widths = c(1, 1.2)
        )
      } else if(length(plots_for_combine) == 3) {
        combined_plot <- grid.arrange(
          plots_for_combine[[1]], plots_for_combine[[2]], plots_for_combine[[3]],
          ncol = 3,
          top = textGrob("Pre-operative OCTA Parameters by Late Recovery Clusters (logMAR Vision)", 
                         gp = gpar(fontsize = 16, fontface = "bold")),
          widths = c(1, 1, 1.2)
        )
      }
      
      ggsave("plots/baseline_octa_combined_grouped_boxplots_late_recovery_logmar.pdf", combined_plot, 
             width = 15, height = 5, device = "pdf")
      ggsave("plots/baseline_octa_combined_grouped_boxplots_late_recovery_logmar.png", combined_plot, 
             width = 15, height = 5, dpi = 300)
      
      cat("    ✓ Combined baseline OCTA grouped boxplots saved\n")
    }
    
    return(plots_list)
  }
  
  octa_baseline_plots <- create_grouped_octa_baseline_boxplots(data, analysis_params)
  
  # PCA分析（修改标题以反映logMAR）
  if(length(analysis_params) > 2) {
    
    pca_data <- data %>%
      dplyr::select(all_of(analysis_params)) %>%
      na.omit()
    
    if(nrow(pca_data) > 3 && ncol(pca_data) > 1) {
      
      pca_result <- prcomp(pca_data, scale. = TRUE)
      
      # 获取对应的聚类信息
      pca_indices <- as.numeric(rownames(pca_data))
      cluster_info <- data[pca_indices, c("max_cluster", "max_membership")]
      
      pca_plot_data <- data.frame(
        PC1 = pca_result$x[,1],
        PC2 = pca_result$x[,2],
        Cluster = cluster_info$max_cluster,
        Membership = cluster_info$max_membership
      )
      
      p_pca <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = as.factor(Cluster), 
                                         alpha = Membership)) +
        geom_point(size = 3) +
        stat_ellipse(aes(group = Cluster), level = 0.95) +
        scale_color_manual(
          values = c("1" = "#a388bf", "2" = "#bc982f"),
          name = "Late Recovery\nCluster"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.title = element_text(face = "bold")
        ) +
        labs(
          title = "PCA of Pre-operative Baseline Characteristics (logMAR Vision)",
          subtitle = "Based on Late Recovery Clustering | T0 OCTA (0_21 Wide-field) + logMAR Vision parameters",
          x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
          y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
          caption = "Ellipses show 95% confidence regions | Alpha indicates cluster membership confidence | Focus: 0_21 region + logMAR vision"
        )
      
      ggsave("plots/baseline_pca_late_recovery_logmar.pdf", p_pca, width = 10, height = 8)
      ggsave("plots/baseline_pca_late_recovery_logmar.png", p_pca, width = 10, height = 8, dpi = 300)
    }
  }
  
  cat("Late Recovery基线可视化（logMAR版本）完成！\n")
}

# ================== 8. 创建logMAR基线特征分组热图 ==================

create_baseline_characteristics_heatmap_lr_logmar <- function(data, stats_results) {
  
  cat("\n===== 创建Late Recovery术前基线特征分组热图 (logMAR版本) =====\n")
  
  # 确定要分析的参数
  analysis_params <- names(data)[!names(data) %in% c("ID", "max_cluster", "max_membership")]
  
  # 按类别分组参数
  baseline_params <- analysis_params[analysis_params %in% c("pre_vision", "age", "gender")]
  bloodflow_params <- analysis_params[grepl("SVP|ICP|DCP|Choroid", analysis_params) & grepl("_T0$", analysis_params)]
  thickness_params <- analysis_params[grepl("GCL|INL|Retina", analysis_params) & grepl("_T0$", analysis_params)]
  
  # 计算每个cluster的均值
  calculate_group_means <- function(param_list, category_name) {
    if(length(param_list) == 0) return(NULL)
    
    means_data <- data %>%
      group_by(max_cluster) %>%
      summarise(across(all_of(param_list), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
      pivot_longer(cols = -max_cluster, names_to = "Parameter", values_to = "Mean_Value") %>%
      mutate(
        Category = category_name,
        Parameter_Clean = case_when(
          Parameter == "pre_vision" ~ "logMAR Vision",
          TRUE ~ gsub("_T0$|_", " ", Parameter)
        )
      )
    
    return(means_data)
  }
  
  # 计算各类别的均值
  baseline_means <- calculate_group_means(baseline_params, "Baseline Characteristics (logMAR)")
  bloodflow_means <- calculate_group_means(bloodflow_params, "Blood Flow - Wide-field (0_21)")
  thickness_means <- calculate_group_means(thickness_params, "Thickness - Wide-field (0_21)")
  
  # 合并所有数据
  all_means <- bind_rows(baseline_means, bloodflow_means, thickness_means) %>%
    filter(!is.na(Mean_Value))
  
  if(nrow(all_means) == 0) {
    cat("警告：没有可用的参数数据！\n")
    return(NULL)
  }
  
  # 创建分面热图
  p_heatmap <- ggplot(all_means, aes(x = Parameter_Clean, y = as.factor(max_cluster), fill = Mean_Value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Mean_Value, 2)), color = "black", size = 3, fontface = "bold") +
    facet_wrap(~ Category, scales = "free_x", ncol = 1, strip.position = "top") +
    scale_fill_gradient2(
      low = "#542788", 
      mid = "white", 
      high = "#f1a340", 
      midpoint = 0,
      name = "Mean\nValue"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold", color = "black"),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = "Pre-operative Baseline Characteristics by Late Recovery Clusters (logMAR Vision)",
      subtitle = paste("T0 OCTA (0_21 Wide-field) + logMAR Vision | n =", nrow(data)),
      x = "",
      y = "Late Recovery Cluster",
      caption = "Values show pre-operative measurements before surgery | logMAR: lower = better vision | Focus: 0_21 region"
    )
  
  # 保存图片
  ggsave("plots/baseline_characteristics_heatmap_late_recovery_logmar.pdf", p_heatmap, 
         width = 16, height = 10, device = "pdf")
  ggsave("plots/baseline_characteristics_heatmap_late_recovery_logmar.png", p_heatmap, 
         width = 16, height = 10, dpi = 300)
  
  cat("✓ Late Recovery基线特征热图（logMAR版本）已保存\n")
  
  # 创建显著差异参数的重点热图
  if(nrow(stats_results) > 0) {
    significant_params <- stats_results %>% 
      filter(Significant == "Yes") %>% 
      pull(Parameter)
    
    if(length(significant_params) > 0) {
      significant_data <- all_means %>%
        filter(gsub("_T0$", "", Parameter) %in% significant_params)
      
      if(nrow(significant_data) > 0) {
        p_heatmap_sig <- ggplot(significant_data, aes(x = Parameter_Clean, y = as.factor(max_cluster), fill = Mean_Value)) +
          geom_tile(color = "white", size = 0.8) +
          geom_text(aes(label = round(Mean_Value, 2)), color = "black", size = 4, fontface = "bold") +
          facet_wrap(~ Category, scales = "free_x", ncol = 1) +
          scale_fill_gradient2(
            low = "#542788", 
            mid = "white", 
            high = "#f1a340", 
            midpoint = 0,
            name = "Mean\nValue"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
            axis.text.y = element_text(size = 12, face = "bold"),
            strip.text = element_text(size = 12, face = "bold"),
            panel.grid = element_blank()
          ) +
          labs(
            title = "Significantly Different Pre-operative Parameters (logMAR Vision)",
            subtitle = "Late Recovery Clusters - Only parameters with significant baseline differences",
            x = "Baseline Parameters (0_21 Region + logMAR Vision)",
            y = "Late Recovery Cluster"
          )
        
        ggsave("plots/baseline_significant_heatmap_late_recovery_logmar.pdf", p_heatmap_sig, 
               width = 12, height = 8)
        ggsave("plots/baseline_significant_heatmap_late_recovery_logmar.png", p_heatmap_sig, 
               width = 12, height = 8, dpi = 300)
        
        cat("✓ Late Recovery显著差异参数热图（logMAR版本）已保存\n")
      }
    }
  }
  
  return(list(
    all_means = all_means,
    plot = p_heatmap
  ))
}

# ================== 9. 执行Late Recovery基线差异统计分析（logMAR版本）==================

# 运行统计分析
baseline_stats_lr_logmar <- analyze_baseline_differences_late_recovery(baseline_comprehensive)

# 保存统计结果
write.csv(baseline_stats_lr_logmar, "baseline_differences_statistics_late_recovery_logmar.csv", row.names = FALSE)

cat("\n===== Late Recovery统计分析结果摘要（logMAR版本）=====\n")
cat("总分析参数:", nrow(baseline_stats_lr_logmar), "\n")
cat("统计显著参数:", sum(baseline_stats_lr_logmar$Significant == "Yes"), "\n")
cat("矫正后显著参数:", sum(baseline_stats_lr_logmar$Significant_Adjusted == "Yes"), "\n")

# 显示前几个最显著的结果
if(nrow(baseline_stats_lr_logmar) > 0) {
  top_results <- baseline_stats_lr_logmar %>% 
    filter(Significant_Adjusted == "Yes") %>% 
    arrange(P_Adjusted) %>% 
    head(5)
  
  if(nrow(top_results) > 0) {
    cat("\n前5个最显著的差异（logMAR版本）:\n")
    print(top_results %>% 
            dplyr::select(Parameter, Data_Type, Region, P_Adjusted, Effect_Size, Test_Method))
  } else {
    cat("\n没有发现统计显著的基线差异\n")
  }
}

# 特别检查logMAR视力差异
vision_result <- baseline_stats_lr_logmar %>% filter(Parameter == "pre_vision")
if(nrow(vision_result) > 0) {
  cat("\n🔍 logMAR视力差异详情:\n")
  cat("- Cluster 1 平均logMAR:", round(vision_result$Cluster1_Mean, 3), "\n")
  cat("- Cluster 2 平均logMAR:", round(vision_result$Cluster2_Mean, 3), "\n")
  cat("- 差异 (Cluster2-Cluster1):", round(vision_result$Mean_Difference, 3), "\n")
  cat("- P值:", round(vision_result$P_Value, 4), "\n")
  cat("- 调整后P值:", round(vision_result$P_Adjusted, 4), "\n")
  cat("- 效应量 (Cohen's d):", round(vision_result$Effect_Size, 3), "\n")
  cat("- 显著性:", vision_result$Significant_Adjusted, "\n")
  
  # 转换为近似Snellen等效值
  cluster1_snellen <- paste0("20/", round(10^vision_result$Cluster1_Mean * 20))
  cluster2_snellen <- paste0("20/", round(10^vision_result$Cluster2_Mean * 20))
  cat("- Cluster 1 近似Snellen:", cluster1_snellen, "\n")
  cat("- Cluster 2 近似Snellen:", cluster2_snellen, "\n")
  
  # 临床意义解释
  if(abs(vision_result$Mean_Difference) >= 0.3) {
    cat("- 临床意义: 大于1.5行视力差异（临床显著）\n")
  } else if(abs(vision_result$Mean_Difference) >= 0.2) {
    cat("- 临床意义: 约1行视力差异（临床相关）\n")
  } else {
    cat("- 临床意义: 小于1行视力差异（临床轻微）\n")
  }
}

# 创建可视化
create_baseline_visualizations_late_recovery_logmar(baseline_comprehensive, baseline_stats_lr_logmar)

# 生成基线特征热图
if(nrow(baseline_comprehensive) > 0 && exists("baseline_stats_lr_logmar")) {
  baseline_heatmap_results_lr_logmar <- create_baseline_characteristics_heatmap_lr_logmar(baseline_comprehensive, baseline_stats_lr_logmar)
  if(!is.null(baseline_heatmap_results_lr_logmar)) {
    cat("✓ Late Recovery基线特征分组热图（logMAR版本）生成完成\n")
  }
} else {
  cat("警告：无法生成热图 - 检查数据或统计分析结果\n")
}

# ================== 10. 效应量分析（logMAR版本）==================

analyze_effect_sizes_lr_logmar <- function(stats_results) {
  
  cat("\n===== Late Recovery术前基线差异效应量分析（logMAR版本）=====\n")
  
  if(nrow(stats_results) == 0) {
    cat("无统计结果可分析\n")
    return(NULL)
  }
  
  # 按效应量大小分类
  effect_summary <- stats_results %>%
    filter(!is.na(Effect_Size)) %>%
    mutate(
      Effect_Magnitude = case_when(
        Effect_Size < 0.2 ~ "Negligible (< 0.2)",
        Effect_Size < 0.5 ~ "Small (0.2-0.5)",
        Effect_Size < 0.8 ~ "Medium (0.5-0.8)",
        TRUE ~ "Large (≥ 0.8)"
      )
    ) %>%
    group_by(Data_Type, Effect_Magnitude) %>%
    summarise(
      Count = n(),
      Mean_Effect_Size = mean(Effect_Size),
      Significant_Count = sum(Significant == "Yes"),
      .groups = 'drop'
    )
  
  cat("效应量分布（logMAR版本）:\n")
  print(effect_summary)
  
  # 重要发现的解释
  large_effects <- stats_results %>%
    filter(Effect_Size >= 0.8) %>%
    arrange(desc(Effect_Size))
  
  medium_effects <- stats_results %>%
    filter(Effect_Size >= 0.5 & Effect_Size < 0.8) %>%
    arrange(desc(Effect_Size))
  
  if(nrow(large_effects) > 0) {
    cat("\n🔍 大效应量参数 (≥ 0.8, logMAR版本):\n")
    print(large_effects %>% 
            dplyr::select(Parameter, Data_Type, Region, Effect_Size, P_Value, Significant))
  }
  
  if(nrow(medium_effects) > 0) {
    cat("\n📊 中等效应量参数 (0.5-0.8, logMAR版本):\n")
    print(medium_effects %>% 
            dplyr::select(Parameter, Data_Type, Region, Effect_Size, P_Value, Significant))
  }
  
  # logMAR视力特别分析
  vision_effect <- stats_results %>% filter(Parameter == "pre_vision")
  if(nrow(vision_effect) > 0) {
    cat("\n👁️ logMAR视力效应量分析:\n")
    cat("- 效应量:", round(vision_effect$Effect_Size, 3), "\n")
    cat("- 效应量级别:", case_when(
      vision_effect$Effect_Size >= 0.8 ~ "大效应",
      vision_effect$Effect_Size >= 0.5 ~ "中等效应", 
      vision_effect$Effect_Size >= 0.2 ~ "小效应",
      TRUE ~ "可忽略效应"
    ), "\n")
    cat("- 临床解释: logMAR差异的", 
        ifelse(vision_effect$Effect_Size >= 0.5, "中等到大", "小到中等"), "临床意义\n")
  }
  
  # 总体结论
  total_significant <- sum(stats_results$Significant == "Yes")
  total_large_medium <- sum(stats_results$Effect_Size >= 0.5, na.rm = TRUE)
  
  cat("\n📋 效应量总结（logMAR版本）:\n")
  cat("- 总分析参数:", nrow(stats_results), "\n")
  cat("- 统计显著参数:", total_significant, "\n")
  cat("- 中等及以上效应量:", total_large_medium, "\n")
  
  return(list(
    effect_summary = effect_summary,
    large_effects = large_effects,
    medium_effects = medium_effects,
    total_significant = total_significant,
    total_large_medium = total_large_medium,
    vision_effect = vision_effect
  ))
}

# 执行效应量分析
effect_analysis_lr_logmar <- analyze_effect_sizes_lr_logmar(baseline_stats_lr_logmar)

# ================== 11. 临床意义解释（logMAR版本）==================

interpret_baseline_findings_lr_logmar <- function(stats_results, effect_analysis) {
  
  cat("\n===== Late Recovery术前基线发现的临床解释（logMAR版本）=====\n")
  
  # 判断患者差异的来源
  significant_count <- sum(stats_results$Significant == "Yes", na.rm = TRUE)
  large_effect_count <- sum(stats_results$Effect_Size >= 0.8, na.rm = TRUE)
  medium_effect_count <- sum(stats_results$Effect_Size >= 0.5 & stats_results$Effect_Size < 0.8, na.rm = TRUE)
  
  # 检查logMAR视力的特殊意义
  vision_significant <- any(stats_results$Parameter == "pre_vision" & stats_results$Significant == "Yes")
  vision_effect_size <- ifelse("vision_effect" %in% names(effect_analysis) && nrow(effect_analysis$vision_effect) > 0, 
                               effect_analysis$vision_effect$Effect_Size, 0)
  
  cat("🎯 Late Recovery阶段主要发现总结（logMAR视力分析）:\n")
  
  if(significant_count == 0) {
    cat("✨ 核心发现：Late Recovery期患者群体差异主要在恢复过程中体现\n")
    cat("📋 临床意义（logMAR视力分析）：\n")
    cat("  - 术前患者在logMAR视力和OCTA方面相对同质\n")
    cat("  - Late recovery阶段的差异更多来自个体恢复能力和中后期因素\n")
    cat("  - 术前logMAR视力不是Late Recovery预后的强预测因子\n")
    cat("  - 需要关注恢复中后期的干预和个性化康复\n")
    cat("  - 预后差异更多来自长期恢复响应性而非基线视力状态\n\n")
    
    interpretation <- "Late_Recovery_Acquired_logMAR"
    
  } else if(vision_significant && vision_effect_size >= 0.5) {
    cat("👁️ 核心发现：术前logMAR视力差异可预测Late Recovery模式\n")
    cat("📋 临床意义（logMAR视力重点分析）：\n")
    cat("  - 术前logMAR视力是Late Recovery的重要预测因子\n")
    cat("  - 较差的术前logMAR视力可能预示更慢的Late Recovery\n")
    cat("  - 可建立基于logMAR视力的Late Recovery风险分层\n")
    cat("  - 术前logMAR优化可能改善长期恢复结局\n")
    cat("  - 需要个性化的logMAR视力导向Late Recovery方案\n\n")
    
    interpretation <- "Vision_Driven_Late_Recovery"
    
  } else if(significant_count <= 3 && large_effect_count == 0) {
    cat("🔍 核心发现：术前存在轻微差异，但Late Recovery差异主要在恢复过程中显现\n")
    cat("📋 临床意义（logMAR综合分析）：\n")
    cat("  - 术前有某些预测因子（可能包括logMAR视力），但对Late Recovery预测能力有限\n")
    cat("  - 中后期恢复因素仍是主要决定因素\n")
    cat("  - 可进行基础的术前风险分层（包括logMAR评估）\n")
    cat("  - 重点仍应放在长期恢复管理优化\n\n")
    
    interpretation <- "Mixed_Late_Recovery_Dominant_logMAR"
    
  } else if(large_effect_count > 0) {
    cat("⚡ 核心发现：患者群体在术前即存在重要差异，可预测Late Recovery模式\n")
    cat("📋 临床意义（logMAR多因素分析）：\n")
    cat("  - 存在明确的Late Recovery术前预测因子（可能包括logMAR视力）\n")
    cat("  - 可建立有效的长期预后风险分层系统\n")
    cat("  - 'Late Recovery High-risk' vs 'Low-risk' 患者识别\n")
    cat("  - 个性化长期恢复方案和术前优化策略\n")
    cat("  - logMAR视力可能是关键的预测变量之一\n\n")
    
    interpretation <- "Late_Recovery_Pre_operative_logMAR"
    
  } else {
    cat("🎭 核心发现：术前和Late Recovery因素共同作用\n")
    cat("📋 临床意义（logMAR整合分析）：\n")
    cat("  - 多因素综合预测Late Recovery模型（包括logMAR视力）\n")
    cat("  - 术前评估 + 长期恢复监测的组合策略\n")
    cat("  - 个性化全程Late Recovery管理方案\n")
    cat("  - logMAR视力作为重要的基线评估指标\n\n")
    
    interpretation <- "Mixed_Late_Recovery_Balanced_logMAR"
  }
  
  # 详细参数解释（特别关注logMAR视力）
  if(nrow(stats_results) > 0) {
    
    cat("🔬 具体参数解读（Late Recovery聚类 - logMAR重点）:\n")
    
    # logMAR视力优先分析
    vision_params <- stats_results %>% filter(Parameter == "pre_vision")
    if(nrow(vision_params) > 0) {
      cat("logMAR视力分析:\n")
      param <- vision_params[1, ]
      sig_status <- ifelse(param$Significant == "Yes", "显著", "不显著")
      cat(sprintf("  - 术前logMAR视力: %s差异 (p=%.3f, d=%.2f)\n", 
                  sig_status, param$P_Value, param$Effect_Size))
      cat(sprintf("  - Cluster 1 平均: %.3f logMAR (约20/%d)\n", 
                  param$Cluster1_Mean, round(10^param$Cluster1_Mean * 20)))
      cat(sprintf("  - Cluster 2 平均: %.3f logMAR (约20/%d)\n", 
                  param$Cluster2_Mean, round(10^param$Cluster2_Mean * 20)))
      cat(sprintf("  - 临床意义: %.1f行视力差异\n", abs(param$Mean_Difference) / 0.1))
      cat("\n")
    }
    
    # 其他基线特征
    other_baseline <- stats_results %>% filter(Data_Type == "Baseline" & Parameter != "pre_vision")
    if(nrow(other_baseline) > 0) {
      cat("其他基线特征:\n")
      for(i in 1:nrow(other_baseline)) {
        param <- other_baseline[i, ]
        sig_status <- ifelse(param$Significant == "Yes", "显著", "不显著")
        cat(sprintf("  - %s: %s差异 (p=%.3f, d=%.2f)\n", 
                    param$Parameter, sig_status, param$P_Value, param$Effect_Size))
      }
      cat("\n")
    }
    
    # OCTA参数
    octa_params <- stats_results %>% filter(Data_Type %in% c("Blood Flow", "Thickness"))
    if(nrow(octa_params) > 0) {
      for(data_type in unique(octa_params$Data_Type)) {
        cat(data_type, "参数:\n")
        type_params <- octa_params %>% filter(Data_Type == data_type)
        
        for(region in unique(type_params$Region)) {
          region_params <- type_params %>% filter(Region == region)
          if(nrow(region_params) > 0) {
            cat("  ", region, "区域:\n")
            significant_in_region <- sum(region_params$Significant == "Yes")
            total_in_region <- nrow(region_params)
            cat(sprintf("    显著差异: %d/%d 参数\n", significant_in_region, total_in_region))
            
            if(significant_in_region > 0) {
              sig_params <- region_params %>% filter(Significant == "Yes") %>% arrange(P_Value)
              for(j in 1:min(3, nrow(sig_params))) {
                param <- sig_params[j, ]
                cat(sprintf("      %s (p=%.3f, d=%.2f)\n", 
                            param$Parameter, param$P_Value, param$Effect_Size))
              }
            }
          }
        }
        cat("\n")
      }
    }
  }
  
  # 临床建议（logMAR视力重点）
  cat("🎯 Late Recovery临床应用建议（logMAR视力导向）:\n")
  
  if(interpretation == "Vision_Driven_Late_Recovery") {
    cat("1. 建立术前logMAR视力为核心的Late Recovery预测模型\n")
    cat("2. 差视力患者的术前优化和特殊准备\n")
    cat("3. 个性化logMAR视力导向的长期恢复期望管理\n")
    cat("4. 基于术前logMAR的分层Late Recovery随访方案\n")
    
  } else if(interpretation == "Late_Recovery_Acquired_logMAR") {
    cat("1. 重点投入Late Recovery期护理和康复优化\n")
    cat("2. 建立长期恢复监测系统（术前logMAR不是主要预测因子）\n")
    cat("3. 个性化中后期康复方案\n")
    cat("4. 关注恢复plateau期的干预\n")
    
  } else if(interpretation == "Late_Recovery_Pre_operative_logMAR") {
    cat("1. 建立包含logMAR视力的综合Late Recovery风险评估模型\n")
    cat("2. 高风险患者术前优化和counseling\n")
    cat("3. 个性化长期恢复期望管理\n")
    cat("4. 分层长期随访方案\n")
    
  } else {
    cat("1. 综合术前评估（包括logMAR视力）+ 长期恢复监测\n")
    cat("2. 动态Late Recovery风险分层\n")
    cat("3. 全程个性化恢复管理\n")
    cat("4. 多时点长期预测模型\n")
  }
  
  return(list(
    interpretation = interpretation,
    significant_count = significant_count,
    large_effect_count = large_effect_count,
    medium_effect_count = medium_effect_count,
    vision_significant = vision_significant,
    vision_effect_size = vision_effect_size
  ))
}

# 执行临床解释
clinical_interpretation_lr_logmar <- interpret_baseline_findings_lr_logmar(baseline_stats_lr_logmar, effect_analysis_lr_logmar)

# ================== 12. 生成Late Recovery综合报告（logMAR版本）==================

generate_late_recovery_baseline_report_logmar <- function(data, stats_results, clinical_interpretation, effect_analysis) {
  
  report <- paste0(
    "========================================\n",
    "Late Recovery时间窗口聚类术前基线特征差异分析报告\n",
    "logMAR视力版本 - OCTA T0 + logMAR Vision Analysis\n",
    "========================================\n\n",
    
    "🎯 研究目的:\n",
    "基于Late Recovery时间窗口聚类结果，分析患者群体差异是:\n",
    "A) 术前即存在（先天差异）- 可预测Late Recovery模式\n",
    "B) Late Recovery恢复过程中体现（获得性差异）- 主要看长期恢复\n",
    "📍 专注分析：0_21区域（广角区域）+ logMAR视力参数\n",
    "🔍 特别关注：logMAR视力作为Late Recovery的预测因子\n\n",
    
    "📊 数据概况（logMAR版本）:\n",
    "- 聚类来源: Late Recovery时间窗口聚类\n",
    "- 分析患者数: ", nrow(data), "\n",
    "- 术前参数总数: ", ncol(data) - 3, "\n",
    "- logMAR视力/基线参数: ", sum(stats_results$Data_Type == "Baseline"), "\n",
    "- OCTA血流参数 (0_21): ", sum(stats_results$Data_Type == "Blood Flow"), "\n",
    "- OCTA厚度参数 (0_21): ", sum(stats_results$Data_Type == "Thickness"), "\n\n",
    
    "🔍 核心发现（logMAR重点分析）:\n",
    "- 统计显著差异参数: ", clinical_interpretation$significant_count, "\n",
    "- 大效应量参数 (d≥0.8): ", clinical_interpretation$large_effect_count, "\n",
    "- 中等效应量参数 (d≥0.5): ", clinical_interpretation$medium_effect_count, "\n",
    "- logMAR视力显著性: ", ifelse(clinical_interpretation$vision_significant, "是", "否"), "\n",
    "- logMAR视力效应量: ", round(clinical_interpretation$vision_effect_size, 3), "\n",
    "- Late Recovery差异来源: ", clinical_interpretation$interpretation, "\n\n"
  )
  
  # 添加logMAR视力的具体分析
  vision_stats <- stats_results %>% filter(Parameter == "pre_vision")
  if(nrow(vision_stats) > 0) {
    report <- paste0(report,
                     "👁️ logMAR视力详细分析:\n",
                     sprintf("- Cluster 1 平均logMAR: %.3f (约20/%d)\n", 
                             vision_stats$Cluster1_Mean, round(10^vision_stats$Cluster1_Mean * 20)),
                     sprintf("- Cluster 2 平均logMAR: %.3f (约20/%d)\n", 
                             vision_stats$Cluster2_Mean, round(10^vision_stats$Cluster2_Mean * 20)),
                     sprintf("- 差异: %.3f logMAR (约%.1f行视力差异)\n", 
                             vision_stats$Mean_Difference, abs(vision_stats$Mean_Difference) / 0.1),
                     sprintf("- 统计显著性: p=%.4f (调整后p=%.4f)\n", 
                             vision_stats$P_Value, vision_stats$P_Adjusted),
                     sprintf("- 效应量: %.3f (%s)\n", vision_stats$Effect_Size,
                             case_when(
                               vision_stats$Effect_Size >= 0.8 ~ "大效应",
                               vision_stats$Effect_Size >= 0.5 ~ "中等效应",
                               vision_stats$Effect_Size >= 0.2 ~ "小效应",
                               TRUE ~ "可忽略效应"
                             )), "\n")
  }
  
  # 添加具体发现
  if(clinical_interpretation$significant_count > 0) {
    significant_params <- stats_results %>% 
      filter(Significant == "Yes") %>% 
      arrange(P_Value)
    
    report <- paste0(report,
                     "📋 显著差异参数详情 (0_21区域, Late Recovery聚类, logMAR版本):\n")
    
    for(i in 1:min(5, nrow(significant_params))) {
      param <- significant_params[i, ]
      report <- paste0(report,
                       sprintf("  %d. %s (%s - %s)\n", i, param$Parameter, param$Data_Type, param$Region),
                       sprintf("     差异: %.3f, p=%.4f, 效应量=%.2f\n", 
                               param$Mean_Difference, param$P_Value, param$Effect_Size))
    }
    report <- paste0(report, "\n")
  }
  
  # 临床意义（特别强调logMAR）
  report <- paste0(report,
                   "🏥 Late Recovery临床意义（logMAR视力导向）:\n")
  
  if(clinical_interpretation$interpretation == "Vision_Driven_Late_Recovery") {
    report <- paste0(report,
                     "👁️ 术前logMAR视力是Late Recovery的重要预测因子\n",
                     "💡 提示：\n",
                     "  - 建立logMAR视力为核心的Late Recovery预测模型\n",
                     "  - 差视力患者需要特殊的术前准备和期望管理\n",
                     "  - 个性化logMAR导向的长期恢复方案\n",
                     "  - 基于术前视力的分层Late Recovery随访\n\n")
  } else if(clinical_interpretation$interpretation == "Late_Recovery_Acquired_logMAR") {
    report <- paste0(report,
                     "✨ 患者在术前相对同质（包括logMAR视力），Late Recovery差异主要在长期恢复过程中体现\n",
                     "💡 提示：\n",
                     "  - 术前logMAR视力不是主要的Late Recovery预测因子\n",
                     "  - 长期恢复能力和中后期因素是关键\n",
                     "  - 重点投入长期恢复管理优化\n",
                     "  - 建立Late Recovery期干预策略\n\n")
  } else if(clinical_interpretation$interpretation == "Late_Recovery_Pre_operative_logMAR") {
    report <- paste0(report,
                     "⚡ 患者在术前即存在重要差异（可能包括logMAR视力），可预测Late Recovery结局\n",
                     "💡 提示：\n",
                     "  - 建立包含logMAR视力的综合Late Recovery风险分层系统\n",
                     "  - 高风险患者术前优化和counseling\n",
                     "  - 个性化长期恢复方案选择\n",
                     "  - Late Recovery期望管理\n\n")
  } else {
    report <- paste0(report,
                     "🎭 术前和Late Recovery因素共同作用，需综合管理（logMAR视力为重要因子）\n",
                     "💡 提示：\n",
                     "  - 建立包含logMAR视力的多因素Late Recovery预测模型\n",
                     "  - 术前评估（包括logMAR） + 长期恢复监测结合\n",
                     "  - 全程个性化Late Recovery管理策略\n",
                     "  - 动态长期风险分层系统\n\n")
  }
  
  # 研究价值
  report <- paste0(report,
                   "🔬 科学价值（logMAR重点）:\n",
                   "1. 明确了Late Recovery患者差异的时间起源和logMAR视力的作用\n",
                   "2. 专注0_21区域（广角）+ logMAR视力提供Late Recovery重要洞察\n",
                   "3. 为logMAR视力导向的Late Recovery个性化医疗提供证据基础\n",
                   "4. 指导长期恢复资源分配策略\n",
                   "5. 支持包含logMAR视力的Late Recovery预测模型开发\n\n",
                   
                   "📈 下一步研究（logMAR导向）:\n",
                   "1. 扩大样本量验证logMAR视力对Late Recovery的预测能力\n",
                   "2. 纵向随访验证logMAR视力的长期预测价值\n",
                   "3. 开发包含logMAR视力的Late Recovery临床决策支持工具\n",
                   "4. 多中心logMAR视力导向的Late Recovery验证研究\n",
                   "5. 探索术前logMAR视力优化对Late Recovery结局的影响\n\n",
                   
                   "📁 输出文件（logMAR版本）:\n",
                   "- baseline_differences_statistics_late_recovery_logmar.csv: 详细统计结果\n",
                   "- logmar_vision_statistics_by_cluster.csv: logMAR视力聚类统计\n",
                   "- plots/baseline_characteristics_heatmap_late_recovery_logmar.pdf: 术前特征热图\n",
                   "- plots/baseline_logmar_vision_boxplot_late_recovery.pdf: logMAR视力专用箱线图\n",
                   "- plots/baseline_*_late_recovery_logmar.pdf: 带p值的Late Recovery箱线图\n",
                   "- plots/baseline_pca_late_recovery_logmar.pdf: 术前特征PCA分析\n\n",
                   
                   "🎯 分析重点（logMAR版本）:\n",
                   "- 专注0_21区域（广角区域）OCTA参数\n",
                   "- 结合术前logMAR视力和基线特征\n",
                   "- 基于Late Recovery时间窗口聚类\n",
                   "- 所有可视化均标注Late Recovery和logMAR\n",
                   "- 提供logMAR视力导向的Late Recovery明确临床指导\n",
                   "- logMAR视力作为核心预测因子的重点分析\n\n",
                   
                   "📊 logMAR视力解读指南:\n",
                   "- logMAR = 0: 20/20视力（正常）\n",
                   "- logMAR = 0.3: 20/40视力（轻度受损）\n",
                   "- logMAR = 0.7: 20/100视力（中度受损）\n",
                   "- logMAR = 1.0: 20/200视力（重度受损）\n",
                   "- logMAR差异0.1 ≈ 1行视力差异\n",
                   "- logMAR差异≥0.3 = 临床显著差异\n\n",
                   
                   "生成时间: ", Sys.time(), "\n",
                   "========================================\n"
  )
  
  # 保存报告
  writeLines(report, "Late_Recovery_Baseline_Analysis_Report_logMAR.txt")
  cat(report)
  
  return(report)
}

# 生成报告
late_recovery_baseline_report_logmar <- generate_late_recovery_baseline_report_logmar(
  baseline_comprehensive, baseline_stats_lr_logmar, clinical_interpretation_lr_logmar, effect_analysis_lr_logmar
)

# ================== 13. 保存所有Late Recovery分析结果（logMAR版本）==================

# 保存基线数据
write.csv(baseline_comprehensive, "baseline_comprehensive_data_late_recovery_logmar.csv", row.names = FALSE)

# 保存效应量分析
if(!is.null(effect_analysis_lr_logmar)) {
  write.csv(effect_analysis_lr_logmar$effect_summary, "effect_size_summary_late_recovery_logmar.csv", row.names = FALSE)
  
  if(nrow(effect_analysis_lr_logmar$large_effects) > 0) {
    write.csv(effect_analysis_lr_logmar$large_effects, "large_effect_parameters_late_recovery_logmar.csv", row.names = FALSE)
  }
  
  if(nrow(effect_analysis_lr_logmar$medium_effects) > 0) {
    write.csv(effect_analysis_lr_logmar$medium_effects, "medium_effect_parameters_late_recovery_logmar.csv", row.names = FALSE)
  }
  
  # 保存logMAR视力特别分析
  if("vision_effect" %in% names(effect_analysis_lr_logmar) && nrow(effect_analysis_lr_logmar$vision_effect) > 0) {
    write.csv(effect_analysis_lr_logmar$vision_effect, "logmar_vision_effect_analysis.csv", row.names = FALSE)
  }
}

# 保存临床解释结果（logMAR版本）
clinical_summary_lr_logmar <- data.frame(
  Analysis_Date = Sys.Date(),
  Clustering_Source = "Late Recovery Time Window",
  Vision_Measure = "logMAR Visual Acuity",
  Total_Patients = nrow(baseline_comprehensive),
  Total_Parameters = ncol(baseline_comprehensive) - 3,
  Significant_Parameters = clinical_interpretation_lr_logmar$significant_count,
  Large_Effects = clinical_interpretation_lr_logmar$large_effect_count,
  Medium_Effects = clinical_interpretation_lr_logmar$medium_effect_count,
  Vision_Significant = clinical_interpretation_lr_logmar$vision_significant,
  Vision_Effect_Size = clinical_interpretation_lr_logmar$vision_effect_size,
  Primary_Difference_Source = clinical_interpretation_lr_logmar$interpretation,
  Clinical_Implication = case_when(
    clinical_interpretation_lr_logmar$interpretation == "Vision_Driven_Late_Recovery" ~ "logMAR vision drives late recovery outcomes",
    clinical_interpretation_lr_logmar$interpretation == "Late_Recovery_Acquired_logMAR" ~ "Focus on late recovery care, vision less predictive",
    clinical_interpretation_lr_logmar$interpretation == "Late_Recovery_Pre_operative_logMAR" ~ "Develop comprehensive late recovery risk stratification",
    TRUE ~ "Comprehensive pre-late recovery management including logMAR"
  )
)

write.csv(clinical_summary_lr_logmar, "clinical_interpretation_summary_late_recovery_logmar.csv", row.names = FALSE)

# 创建logMAR视力对比总结表
if(nrow(baseline_stats_lr_logmar) > 0) {
  vision_comparison_summary <- baseline_comprehensive %>%
    group_by(max_cluster) %>%
    summarise(
      n = n(),
      mean_logmar = mean(pre_vision, na.rm = TRUE),
      sd_logmar = sd(pre_vision, na.rm = TRUE),
      median_logmar = median(pre_vision, na.rm = TRUE),
      q25_logmar = quantile(pre_vision, 0.25, na.rm = TRUE),
      q75_logmar = quantile(pre_vision, 0.75, na.rm = TRUE),
      min_logmar = min(pre_vision, na.rm = TRUE),
      max_logmar = max(pre_vision, na.rm = TRUE),
      # 转换为Snellen等效值
      mean_snellen_equiv = paste0("20/", round(10^mean_logmar * 20)),
      median_snellen_equiv = paste0("20/", round(10^median_logmar * 20)),
      # 视力水平分类
      normal_vision_n = sum(pre_vision <= 0, na.rm = TRUE),
      mild_impair_n = sum(pre_vision > 0 & pre_vision <= 0.3, na.rm = TRUE),
      moderate_impair_n = sum(pre_vision > 0.3 & pre_vision <= 0.7, na.rm = TRUE),
      severe_impair_n = sum(pre_vision > 0.7, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      normal_vision_pct = round(normal_vision_n / n * 100, 1),
      mild_impair_pct = round(mild_impair_n / n * 100, 1),
      moderate_impair_pct = round(moderate_impair_n / n * 100, 1),
      severe_impair_pct = round(severe_impair_n / n * 100, 1)
    )
  
  write.csv(vision_comparison_summary, "logmar_vision_comparison_by_cluster_detailed.csv", row.names = FALSE)
  cat("✓ 详细的logMAR视力聚类对比总结已保存\n")
}

# ================== 14. 最终总结（logMAR版本）==================

cat("\n🎉 Late Recovery时间窗口聚类术前基线特征分析完成！（logMAR视力版本）\n")
cat("========================================\n")

# 显示主要结论
cat("\n🎯 主要结论 (专注0_21区域 + logMAR视力, Late Recovery聚类)：\n")

# 根据分析结果给出不同的结论
if(clinical_interpretation_lr_logmar$interpretation == "Vision_Driven_Late_Recovery") {
  cat("术前logMAR视力是Late Recovery的重要预测因子\n")
  cat("→ 重点：建立logMAR视力导向的Late Recovery预测和管理系统\n")
} else if(clinical_interpretation_lr_logmar$interpretation == "Late_Recovery_Acquired_logMAR") {
  cat("Late Recovery患者群体差异主要在长期恢复过程中体现，术前logMAR视力预测能力有限\n")
  cat("→ 重点：优化Late Recovery期护理和长期恢复管理\n")
} else if(clinical_interpretation_lr_logmar$interpretation == "Late_Recovery_Pre_operative_logMAR") {
  cat("患者在术前即存在重要差异（包括logMAR视力），可预测Late Recovery模式\n")
  cat("→ 重点：建立包含logMAR视力的Late Recovery术前风险分层系统\n")
} else {
  cat("术前存在差异（logMAR视力可能相关），Late Recovery差异需综合管理\n")
  cat("→ 重点：综合术前logMAR评估和Late Recovery管理\n")
}

# 显示logMAR视力的具体发现
vision_result <- baseline_stats_lr_logmar %>% filter(Parameter == "pre_vision")
if(nrow(vision_result) > 0) {
  cat("\n👁️ logMAR视力具体发现：\n")
  cat(sprintf("- 统计显著性: %s (p=%.4f)\n", 
              ifelse(vision_result$Significant == "Yes", "显著", "不显著"),
              vision_result$P_Value))
  cat(sprintf("- 效应量: %.3f (%s)\n", 
              vision_result$Effect_Size,
              case_when(
                vision_result$Effect_Size >= 0.8 ~ "大效应",
                vision_result$Effect_Size >= 0.5 ~ "中等效应", 
                vision_result$Effect_Size >= 0.2 ~ "小效应",
                TRUE ~ "可忽略效应"
              )))
  cat(sprintf("- 视力差异: 约%.1f行 (%.3f logMAR)\n", 
              abs(vision_result$Mean_Difference) / 0.1,
              abs(vision_result$Mean_Difference)))
}

cat("\n📁 主要输出文件（logMAR版本）：\n")
output_files_lr_logmar <- c(
  "baseline_differences_statistics_late_recovery_logmar.csv",
  "baseline_comprehensive_data_late_recovery_logmar.csv", 
  "clinical_interpretation_summary_late_recovery_logmar.csv",
  "logmar_vision_comparison_by_cluster_detailed.csv",
  "logmar_vision_statistics_by_cluster.csv",
  "Late_Recovery_Baseline_Analysis_Report_logMAR.txt"
)

for(file in output_files_lr_logmar) {
  if(file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  }
}

cat("\n📊 可视化文件 (Late Recovery + logMAR标注)：\n")
viz_files_lr_logmar <- list.files("plots", pattern = "late_recovery_logmar\\.(pdf|png)$", full.names = FALSE)
for(file in viz_files_lr_logmar) {
  cat(sprintf("✓ plots/%s\n", file))
}

# 特别突出logMAR视力图
if(file.exists("plots/baseline_logmar_vision_boxplot_late_recovery.pdf")) {
  cat("\n🌟 重点图表：\n")
  cat("✓ plots/baseline_logmar_vision_boxplot_late_recovery.pdf - logMAR视力专用分析图\n")
}

cat("\n🔍 分析特色（logMAR版本）：\n")
cat("1. 所有视力数据转换为标准化logMAR格式\n")
cat("2. 专业的logMAR视力统计分析和可视化\n")
cat("3. logMAR视力作为Late Recovery预测因子的重点评估\n")
cat("4. 临床意义明确的视力差异解释（行数转换）\n")
cat("5. Snellen等效值对照，便于临床理解\n")

cat("\n🎯 临床应用价值：\n")
cat("- 为logMAR视力导向的Late Recovery管理提供科学依据\n")
cat("- 支持个性化的术前视力评估和期望管理\n")
cat("- 指导基于视力状态的长期随访策略\n")
cat("- 促进标准化的视力结局评估体系建立\n\n")

cat("📝 注意事项：\n")
cat("- logMAR值越小表示视力越好\n")
cat("- logMAR差异0.1约等于1行视力差异\n") 
cat("- logMAR差异≥0.3被认为有临床显著意义\n")
cat("- 所有统计分析基于logMAR连续变量进行\n\n")

cat("✅ 分析完成！请查看生成的报告和图表以获取详细结果。\n")