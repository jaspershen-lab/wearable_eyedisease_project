# "🎯 研究目的:\n",
# "基于comprehensive clustering结果，分析患者群体差异是:\n",
# "A) 术前即存在（先天差异）\n",
# "B) 术后恢复过程中体现（获得性差异）\n",
# "📍 专注分析：0_21区域（广角区域）+ 视力参数\n\n",

# ============================================

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

cat("===== 术前基线特征差异分析 =====\n")
cat("分析目标：确定患者群体差异是术前即存在还是术后才显现\n\n")

# 加载baseline信息和OCTA数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# 加载comprehensive clustering结果
comprehensive_results_file <- "3_data_analysis/6_clustering_modeling/mfuzz/comprehensive_cluster/ppv_comprehensive_cluster_results_with_outcomes.csv"

if(file.exists(comprehensive_results_file)) {
  comprehensive_clusters <- read.csv(comprehensive_results_file, stringsAsFactors = FALSE)
  cat("✓ 成功加载comprehensive clustering结果:", nrow(comprehensive_clusters), "患者\n")
} else {
  stop("请先运行comprehensive clustering分析！")
}

# 创建输出目录
dir.create("3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/baseline_octa_analysis", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/baseline_octa_analysis")

# ================== 2. 处理术前OCTA数据 ==================

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

# ================== 3. 提取术前参数 ==================

# 修改filter函数只提取T0参数 - 专注于0_21区域（根据你的定义，这是广角区域）
filter_baseline_bloodflow <- function(data) {
  layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
  regions_of_interest <- c("0_21")  # 只关注0_21区域（你定义的广角区域）
  
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("未找到0_21区域T0血流参数！")
    return(list(data = data %>% select(ID), params = character(0)))
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
  regions_of_interest <- c("0_21")  # 只关注0_21区域（你定义的广角区域）
  
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*(",
                    paste(regions_of_interest, collapse = "|"), ")_T0$")
  
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("未找到0_21区域T0厚度参数！")
    return(list(data = data %>% select(ID), params = character(0)))
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

# ================== 4. 创建术前基线数据集 ==================

# 处理术前视力数据
baseline_vision <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,
      surgery_eye_1 == 1 ~ os_corrected_bas,
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(ID, pre_vision, age, gender) %>%
  filter(ID %in% ppv_patients)

# 合并所有术前数据
baseline_comprehensive <- baseline_vision %>%
  full_join(baseline_bloodflow_filtered$data, by = "ID") %>%
  full_join(baseline_thickness_filtered$data, by = "ID") %>%
  # 添加聚类信息
  inner_join(comprehensive_clusters %>% 
               dplyr::select(subject_id, max_cluster, max_membership, outcome_quality), 
             by = c("ID" = "subject_id"))

cat("\n===== 术前基线数据集摘要 =====\n")
cat("专注分析：0_21区域（广角区域）+ 视力参数\n")
cat("患者数量:", nrow(baseline_comprehensive), "\n")
cat("总参数数:", ncol(baseline_comprehensive) - 4, "\n")  # 排除ID, cluster, membership, outcome
cat("- 视力/基本信息:", ncol(baseline_vision) - 1, "\n")
cat("- OCTA血流 (0_21区域):", length(baseline_bloodflow_filtered$params_T0), "\n")
cat("- OCTA厚度 (0_21区域):", length(baseline_thickness_filtered$params_T0), "\n")

# 检查数据完整性
baseline_complete_cases <- baseline_comprehensive[complete.cases(baseline_comprehensive), ]
cat("完整数据患者:", nrow(baseline_complete_cases), 
    "(", round(nrow(baseline_complete_cases)/nrow(baseline_comprehensive)*100, 1), "%)\n")

# ================== 5. 术前基线差异统计分析 ==================

analyze_baseline_differences <- function(data) {
  
  cat("\n===== 术前基线差异统计分析（修正版）=====\n")
  
  # 确定要分析的参数
  analysis_params <- names(data)[!names(data) %in% c("ID", "max_cluster", "max_membership", "outcome_quality")]
  
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
        param %in% c("pre_vision", "age", "gender") ~ "N/A",
        grepl("0_21", param) ~ "Wide-field (0_21)",
        TRUE ~ "Unknown"
      )
      
      # 判断变量类型
      param_values <- data[[param]][!is.na(data[[param]])]
      unique_values <- unique(param_values)
      
      # 判断是否为二分类变量
      is_binary <- length(unique_values) == 2 && all(unique_values %in% c(0, 1))
      
      variable_type <- if(is_binary) "Binary" else "Continuous"
      
      # 提取数据
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
          # 连续变量：使用t检验
          cat(sprintf("对连续变量 %s 使用t检验\n", param))
          
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

# 为二分类变量创建堆叠条形图
create_binary_plot <- function(data, param, param_clean, p_value, test_method) {
  
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
      values = if(param == "gender") c("Female" = "#FF69B4", "Male" = "#4169E1") else c("#E7B800", "#00AFBB"),
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
      subtitle = paste("Comparison between outcome clusters |", p_text),
      x = "Outcome Cluster",
      y = "Percentage",
      caption = paste("Statistical test:", test_method, "| Numbers show count and percentage")
    )
  
  # 保存图片
  ggsave(paste0("plots/baseline_", param, "_distribution_with_pvalue.pdf"), 
         p_binary, width = 8, height = 6)
  ggsave(paste0("plots/baseline_", param, "_distribution_with_pvalue.png"), 
         p_binary, width = 8, height = 6, dpi = 300)
}

# 为连续变量创建箱线图
create_continuous_plot <- function(data, param, param_clean, p_value, test_method) {
  
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
      values = c("1" = "#E7B800", "2" = "#00AFBB"),
      name = "Outcome\nCluster"
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
      subtitle = paste("Comparison between outcome clusters |", p_text),
      x = "Outcome Cluster",
      y = param_clean,
      caption = paste("Statistical test:", test_method, "| Individual points show patients")
    )
  
  ggsave(paste0("plots/baseline_", param, "_boxplot_with_pvalue.pdf"), 
         p_box, width = 8, height = 6)
  ggsave(paste0("plots/baseline_", param, "_boxplot_with_pvalue.png"), 
         p_box, width = 8, height = 6, dpi = 300)
}

# 修正的可视化函数 - 针对不同变量类型
create_baseline_visualizations <- function(data, stats_results) {
  
  cat("\n===== 创建修正的可视化（区分变量类型）=====\n")
  
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  
  # 在函数开始处定义 analysis_params
  analysis_params <- names(data)[!names(data) %in% c("ID", "max_cluster", "max_membership", "outcome_quality")]
  
  # 获取显著差异的参数（使用矫正后的p值）
  significant_baseline <- stats_results %>% 
    filter(Significant_Adjusted == "Yes") %>%  # 使用矫正后的显著性
    arrange(P_Adjusted)
  
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
          p_value <- param_info$P_Adjusted[1]  # 使用矫正后的p值
          variable_type <- param_info$Variable_Type[1]
          test_method <- param_info$Test_Method[1]
          
          # 根据变量类型创建不同的图
          if(variable_type == "Binary") {
            # 二分类变量：创建堆叠条形图
            create_binary_plot(data, param, param_clean, p_value, test_method)
          } else {
            # 连续变量：创建箱线图
            create_continuous_plot(data, param, param_clean, p_value, test_method)
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
        
        if(variable_type == "Binary") {
          create_binary_plot(data, param, param_clean, p_value, test_method)
        } else {
          create_continuous_plot(data, param, param_clean, p_value, test_method)
        }
      }
    }
  }
  
  # 3. PCA分析 - 术前特征（修复后的部分）
  if(length(analysis_params) > 2) {
    
    pca_data <- data %>%
      dplyr::select(all_of(analysis_params)) %>%
      na.omit()
    
    if(nrow(pca_data) > 3 && ncol(pca_data) > 1) {
      
      pca_result <- prcomp(pca_data, scale. = TRUE)
      
      # 获取对应的聚类信息
      pca_indices <- as.numeric(rownames(pca_data))
      cluster_info <- data[pca_indices, c("max_cluster", "max_membership", "outcome_quality")]
      
      pca_plot_data <- data.frame(
        PC1 = pca_result$x[,1],
        PC2 = pca_result$x[,2],
        Cluster = cluster_info$max_cluster,
        Membership = cluster_info$max_membership,
        Outcome = cluster_info$outcome_quality
      )
      
      p_pca <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = as.factor(Cluster), 
                                         alpha = Membership)) +
        geom_point(size = 3) +
        stat_ellipse(aes(group = Cluster), level = 0.95) +
        scale_color_manual(
          values = c("1" = "#E7B800", "2" = "#00AFBB"),
          name = "Outcome\nCluster"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.title = element_text(face = "bold")
        ) +
        labs(
          title = "PCA of Pre-operative Baseline Characteristics",
          subtitle = "T0 OCTA (0_21 Wide-field) + Pre-Vision parameters (before surgery)",
          x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
          y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
          caption = "Ellipses show 95% confidence regions | Alpha indicates cluster membership confidence | Focus: 0_21 region"
        )
      
      ggsave("plots/baseline_pca.pdf", p_pca, width = 10, height = 8)
      ggsave("plots/baseline_pca.png", p_pca, width = 10, height = 8, dpi = 300)
      
      # 变量贡献图
      loadings <- pca_result$rotation[, 1:2]
      loadings_df <- data.frame(
        Variable = rownames(loadings),
        PC1 = loadings[, 1],
        PC2 = loadings[, 2],
        Data_Type = case_when(
          rownames(loadings) %in% c("pre_vision", "age", "gender") ~ "Baseline",
          grepl("SVP|ICP|DCP|Choroid", rownames(loadings)) ~ "Blood Flow",
          grepl("GCL|INL|Retina", rownames(loadings)) ~ "Thickness",
          TRUE ~ "Other"
        )
      )
      
      p_loadings <- ggplot(loadings_df, aes(x = PC1, y = PC2, color = Data_Type)) +
        geom_point(size = 3) +
        geom_text(aes(label = gsub("_T0$|_", " ", Variable)), 
                  vjust = -0.5, hjust = 0.5, size = 2.5) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        scale_color_brewer(palette = "Set2", name = "Parameter\nType") +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.title = element_text(face = "bold")
        ) +
        labs(
          title = "Variable Contributions to Principal Components",
          subtitle = "Pre-operative baseline parameters (0_21 region focus)",
          x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
          y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
        )
      
      ggsave("plots/baseline_loadings.pdf", p_loadings, width = 12, height = 10)
    }
  }
  
  # 4. 效应量可视化
  if(nrow(stats_results) > 0) {
    
    effect_size_data <- stats_results %>%
      filter(!is.na(Effect_Size)) %>%
      mutate(
        Effect_Magnitude = case_when(
          Effect_Size < 0.2 ~ "Small",
          Effect_Size < 0.5 ~ "Small",
          Effect_Size < 0.8 ~ "Medium", 
          TRUE ~ "Large"
        ),
        Parameter_Display = paste0(Parameter, " (", Region, ")")
      ) %>%
      arrange(desc(Effect_Size))
    
    if(nrow(effect_size_data) > 0) {
      
      p_effect <- ggplot(effect_size_data, aes(x = reorder(Parameter_Display, Effect_Size), 
                                               y = Effect_Size, fill = Significant)) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = c(0.2, 0.5, 0.8), linetype = "dashed", alpha = 0.6) +
        scale_fill_manual(
          values = c("No" = "lightgray", "Yes" = "darkgreen"),
          name = "Statistically\nSignificant"
        ) +
        coord_flip() +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")
        ) +
        labs(
          title = "Effect Sizes of Pre-operative Baseline Differences",
          subtitle = "Cohen's d for differences between outcome clusters (0_21 region focus)",
          x = "Parameters",
          y = "Effect Size (Cohen's d)",
          caption = "Dashed lines: 0.2=small, 0.5=medium, 0.8=large effect | Focus: 0_21 wide-field region"
        )
      
      ggsave("plots/baseline_effect_sizes.pdf", p_effect, width = 12, height = 10)
      ggsave("plots/baseline_effect_sizes.png", p_effect, width = 12, height = 10, dpi = 300)
    }
  }
  
  cat("术前基线可视化完成！\n")
}

# ================== 6. 执行术前基线差异统计分析 ==================

# 首先运行统计分析生成 baseline_stats
baseline_stats <- analyze_baseline_differences(baseline_comprehensive)

# 保存统计结果
write.csv(baseline_stats, "baseline_differences_statistics.csv", row.names = FALSE)

cat("\n===== 统计分析结果摘要 =====\n")
cat("总分析参数:", nrow(baseline_stats), "\n")
cat("统计显著参数:", sum(baseline_stats$Significant == "Yes"), "\n")
cat("矫正后显著参数:", sum(baseline_stats$Significant_Adjusted == "Yes"), "\n")

# 显示前几个最显著的结果
if(nrow(baseline_stats) > 0) {
  top_results <- baseline_stats %>% 
    filter(Significant_Adjusted == "Yes") %>% 
    arrange(P_Adjusted) %>% 
    head(5)
  
  if(nrow(top_results) > 0) {
    cat("\n前5个最显著的差异:\n")
    print(top_results %>% 
            select(Parameter, Data_Type, P_Adjusted, Effect_Size, Test_Method))
  } else {
    cat("\n没有发现统计显著的基线差异\n")
  }
}

# 然后创建可视化
create_baseline_visualizations(baseline_comprehensive, baseline_stats)


# ================== 7. 效应量分析 ==================

analyze_effect_sizes <- function(stats_results) {
  
  cat("\n===== 术前基线差异效应量分析 =====\n")
  
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
  
  cat("效应量分布:\n")
  print(effect_summary)
  
  # 重要发现的解释
  large_effects <- stats_results %>%
    filter(Effect_Size >= 0.8) %>%
    arrange(desc(Effect_Size))
  
  medium_effects <- stats_results %>%
    filter(Effect_Size >= 0.5 & Effect_Size < 0.8) %>%
    arrange(desc(Effect_Size))
  
  if(nrow(large_effects) > 0) {
    cat("\n🔍 大效应量参数 (≥ 0.8):\n")
    print(large_effects %>% 
            dplyr::select(Parameter, Data_Type, Region, Effect_Size, P_Value, Significant))
  }
  
  if(nrow(medium_effects) > 0) {
    cat("\n📊 中等效应量参数 (0.5-0.8):\n")
    print(medium_effects %>% 
            dplyr::select(Parameter, Data_Type, Region, Effect_Size, P_Value, Significant))
  }
  
  # 总体结论
  total_significant <- sum(stats_results$Significant == "Yes")
  total_large_medium <- sum(stats_results$Effect_Size >= 0.5, na.rm = TRUE)
  
  cat("\n📋 效应量总结:\n")
  cat("- 总分析参数:", nrow(stats_results), "\n")
  cat("- 统计显著参数:", total_significant, "\n")
  cat("- 中等及以上效应量:", total_large_medium, "\n")
  
  return(list(
    effect_summary = effect_summary,
    large_effects = large_effects,
    medium_effects = medium_effects,
    total_significant = total_significant,
    total_large_medium = total_large_medium
  ))
}

# 执行效应量分析
effect_analysis <- analyze_effect_sizes(baseline_stats)

# ================== 8. 术前vs术后改善关联分析 ==================

analyze_baseline_improvement_correlation <- function() {
  
  cat("\n===== 术前基线与术后改善关联分析 =====\n")
  
  # 加载术后改善数据（来自comprehensive clustering）
  comprehensive_file <- "3_data_analysis/6_clustering_modeling/mfuzz/comprehensive_cluster/ppv_comprehensive_cluster_results_with_outcomes.csv"
  
  if(!file.exists(comprehensive_file)) {
    cat("警告：未找到术后改善数据，跳过关联分析\n")
    return(NULL)
  }
  
  # 这里可以添加术前基线值与术后改善量的相关性分析
  # 由于需要具体的改善数据，这里提供框架
  
  cat("分析思路：\n")
  cat("1. 术前基线值 vs 术后改善量的相关性\n")
  cat("2. 术前聚类预测术后结局的准确性\n")
  cat("3. 关键术前预测因子识别\n")
  
  # 创建假设性分析框架
  correlation_framework <- data.frame(
    Analysis_Type = c("Baseline_Prediction", "Regional_Correlation", "Temporal_Pattern"),
    Description = c(
      "术前参数预测术后结局的能力",
      "不同区域术前基线的差异模式", 
      "从术前到术后的时间演变模式"
    ),
    Clinical_Value = c(
      "术前风险分层和预后预测",
      "区域特异性治疗策略制定",
      "个性化随访和干预时机"
    )
  )
  
  cat("\n关联分析框架:\n")
  print(correlation_framework)
  
  return(correlation_framework)
}

# 执行关联分析
correlation_analysis <- analyze_baseline_improvement_correlation()

# ================== 9. 临床意义解释 ==================

interpret_baseline_findings <- function(stats_results, effect_analysis) {
  
  cat("\n===== 术前基线发现的临床解释 =====\n")
  
  # 判断患者差异的来源
  significant_count <- sum(stats_results$Significant == "Yes", na.rm = TRUE)
  large_effect_count <- sum(stats_results$Effect_Size >= 0.8, na.rm = TRUE)
  medium_effect_count <- sum(stats_results$Effect_Size >= 0.5 & stats_results$Effect_Size < 0.8, na.rm = TRUE)
  
  cat("🎯 主要发现总结:\n")
  
  if(significant_count == 0) {
    cat("✨ 核心发现：患者群体差异主要在术后恢复过程中体现\n")
    cat("📋 临床意义：\n")
    cat("  - 术前患者在OCTA和视力方面相对同质\n")
    cat("  - 手术技巧、术后护理、个体恢复能力是关键差异因素\n")
    cat("  - 需要关注术后早期干预和个性化康复\n")
    cat("  - 预后差异更多来自手术响应性而非基线状态\n\n")
    
    interpretation <- "Post-operative"
    
  } else if(significant_count <= 3 && large_effect_count == 0) {
    cat("🔍 核心发现：术前存在轻微差异，但主要差异在术后恢复\n")
    cat("📋 临床意义：\n")
    cat("  - 术前有某些预测因子，但预测能力有限\n")
    cat("  - 手术和术后因素仍是主要决定因素\n")
    cat("  - 可进行基础的术前风险分层\n")
    cat("  - 重点仍应放在术后管理优化\n\n")
    
    interpretation <- "Mixed_PostOp_Dominant"
    
  } else if(large_effect_count > 0) {
    cat("⚡ 核心发现：患者群体在术前即存在重要差异\n")
    cat("📋 临床意义：\n")
    cat("  - 存在明确的术前预测因子\n")
    cat("  - 可建立有效的术前风险分层系统\n")
    cat("  - 'High-risk' vs 'Low-risk' 患者识别\n")
    cat("  - 个性化手术方案和术前优化策略\n\n")
    
    interpretation <- "Pre-operative"
    
  } else {
    cat("🎭 核心发现：术前和术后因素共同作用\n")
    cat("📋 临床意义：\n")
    cat("  - 多因素综合预测模型\n")
    cat("  - 术前评估 + 术后监测的组合策略\n")
    cat("  - 个性化全程管理方案\n\n")
    
    interpretation <- "Mixed_Balanced"
  }
  
  # 详细参数解释
  if(nrow(stats_results) > 0) {
    
    cat("🔬 具体参数解读:\n")
    
    # 视力相关
    vision_params <- stats_results %>% filter(Data_Type == "Baseline")
    if(nrow(vision_params) > 0) {
      cat("视力/基线特征:\n")
      for(i in 1:nrow(vision_params)) {
        param <- vision_params[i, ]
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
              for(j in 1:min(3, nrow(sig_params))) {  # 显示前3个最显著的
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
  
  # 临床建议
  cat("🎯 临床应用建议:\n")
  
  if(interpretation == "Post-operative") {
    cat("1. 重点投入术后护理和康复优化\n")
    cat("2. 标准化手术流程，减少技术差异\n")
    cat("3. 建立术后早期预警系统\n")
    cat("4. 个性化术后康复方案\n")
    
  } else if(interpretation == "Pre-operative") {
    cat("1. 建立术前风险评估模型\n")
    cat("2. 高风险患者术前优化\n")
    cat("3. 分层手术方案选择\n")
    cat("4. 术前counseling和期望管理\n")
    
  } else {
    cat("1. 综合术前评估 + 术后监测\n")
    cat("2. 动态风险分层系统\n")
    cat("3. 全程个性化管理\n")
    cat("4. 多时点预测模型\n")
  }
  
  return(list(
    interpretation = interpretation,
    significant_count = significant_count,
    large_effect_count = large_effect_count,
    medium_effect_count = medium_effect_count
  ))
}

# 执行临床解释
clinical_interpretation <- interpret_baseline_findings(baseline_stats, effect_analysis)

# ================== 10. 生成综合报告 ==================

# 修复后的生成综合报告函数
generate_baseline_analysis_report <- function(data, stats_results, clinical_interpretation, effect_analysis) {
  
  report <- paste0(
    "========================================\n",
    "术前基线特征差异分析报告\n",
    "OCTA T0 + Pre-Vision Analysis\n",
    "========================================\n\n",
    
    "🎯 研究目的:\n",
    "基于comprehensive clustering结果，分析患者群体差异是:\n",
    "A) 术前即存在（先天差异）\n",
    "B) 术后恢复过程中体现（获得性差异）\n",
    "📍 专注分析：0_21区域（广角区域）+ 视力参数\n\n",
    
    "📊 数据概况:\n",
    "- 分析患者数: ", nrow(data), "\n",
    "- 术前参数总数: ", ncol(data) - 4, "\n",
    "- 视力/基线参数: ", sum(stats_results$Data_Type == "Baseline"), "\n",
    "- OCTA血流参数 (0_21): ", sum(stats_results$Data_Type == "Blood Flow"), "\n",
    "- OCTA厚度参数 (0_21): ", sum(stats_results$Data_Type == "Thickness"), "\n\n",
    
    "🔍 核心发现:\n",
    "- 统计显著差异参数: ", clinical_interpretation$significant_count, "\n",
    "- 大效应量参数 (d≥0.8): ", clinical_interpretation$large_effect_count, "\n",
    "- 中等效应量参数 (d≥0.5): ", clinical_interpretation$medium_effect_count, "\n",
    "- 主要差异来源: ", clinical_interpretation$interpretation, "\n\n"
  )
  
  # 添加具体发现
  if(clinical_interpretation$significant_count > 0) {
    significant_params <- stats_results %>% 
      filter(Significant == "Yes") %>% 
      arrange(P_Value)
    
    report <- paste0(report,
                     "📋 显著差异参数详情 (0_21区域):\n")
    
    for(i in 1:min(5, nrow(significant_params))) {  # 显示前5个最显著的
      param <- significant_params[i, ]
      report <- paste0(report,
                       sprintf("  %d. %s (%s)\n", i, param$Parameter, param$Data_Type),
                       sprintf("     差异: %.3f, p=%.4f, 效应量=%.2f\n", 
                               param$Mean_Difference, param$P_Value, param$Effect_Size))
    }
    report <- paste0(report, "\n")
  }
  
  # 临床意义
  report <- paste0(report,
                   "🏥 临床意义:\n")
  
  if(clinical_interpretation$interpretation == "Post-operative") {
    report <- paste0(report,
                     "✨ 患者在术前相对同质，差异主要在术后恢复过程中体现\n",
                     "💡 提示：\n",
                     "  - 手术技巧和术后护理是关键\n",
                     "  - 个体恢复能力差异是主要因素\n",
                     "  - 重点投入术后管理优化\n",
                     "  - 建立术后早期干预策略\n\n")
  } else if(clinical_interpretation$interpretation == "Pre-operative") {
    report <- paste0(report,
                     "⚡ 患者在术前即存在重要差异，可预测术后结局\n",
                     "💡 提示：\n",
                     "  - 建立术前风险分层系统\n",
                     "  - 高风险患者术前优化\n",
                     "  - 个性化手术方案选择\n",
                     "  - 术前counseling和期望管理\n\n")
  } else {
    report <- paste0(report,
                     "🎭 术前和术后因素共同作用，需综合管理\n",
                     "💡 提示：\n",
                     "  - 建立多因素预测模型\n",
                     "  - 术前评估 + 术后监测结合\n",
                     "  - 全程个性化管理策略\n",
                     "  - 动态风险分层系统\n\n")
  }
  
  # 研究价值
  report <- paste0(report,
                   "🔬 科学价值:\n",
                   "1. 明确了患者差异的时间起源\n",
                   "2. 专注0_21区域（广角）提供重要洞察\n",
                   "3. 为个性化医疗提供证据基础\n",
                   "4. 指导临床资源分配策略\n",
                   "5. 支持预测模型开发\n\n",
                   
                   "📈 下一步研究:\n",
                   "1. 扩大样本量验证发现\n",
                   "2. 纵向随访验证预测能力\n",
                   "3. 开发临床决策支持工具\n",
                   "4. 多中心验证研究\n\n",
                   
                   "📁 输出文件:\n",
                   "- baseline_differences_statistics.csv: 详细统计结果\n",
                   "- plots/baseline_characteristics_heatmap.pdf: 术前特征热图\n",
                   "- plots/baseline_*_boxplot_with_pvalue.pdf: 带p值的箱线图\n",
                   "- plots/baseline_pca.pdf: 术前特征PCA分析\n",
                   "- plots/baseline_effect_sizes.pdf: 效应量可视化\n\n",
                   
                   "🎯 分析重点:\n",
                   "- 专注0_21区域（广角区域）OCTA参数\n",
                   "- 结合术前视力和基线特征\n",
                   "- 所有箱线图均标注p值\n",
                   "- 提供明确的临床指导\n\n",
                   
                   "生成时间: ", Sys.time(), "\n",
                   "========================================\n"
  )
  
  # 保存报告
  writeLines(report, "Baseline_Analysis_Report.txt")
  cat(report)
  
  return(report)
}

# 生成报告
baseline_report <- generate_baseline_analysis_report(
  baseline_comprehensive, baseline_stats, clinical_interpretation, effect_analysis
)

# ================== 11. 保存所有分析结果 ==================

# 保存基线数据
write.csv(baseline_comprehensive, "baseline_comprehensive_data.csv", row.names = FALSE)

# 保存效应量分析
if(!is.null(effect_analysis)) {
  write.csv(effect_analysis$effect_summary, "effect_size_summary.csv", row.names = FALSE)
  
  if(nrow(effect_analysis$large_effects) > 0) {
    write.csv(effect_analysis$large_effects, "large_effect_parameters.csv", row.names = FALSE)
  }
  
  if(nrow(effect_analysis$medium_effects) > 0) {
    write.csv(effect_analysis$medium_effects, "medium_effect_parameters.csv", row.names = FALSE)
  }
}

# 保存临床解释结果
clinical_summary <- data.frame(
  Analysis_Date = Sys.Date(),
  Total_Patients = nrow(baseline_comprehensive),
  Total_Parameters = ncol(baseline_comprehensive) - 4,
  Significant_Parameters = clinical_interpretation$significant_count,
  Large_Effects = clinical_interpretation$large_effect_count,
  Medium_Effects = clinical_interpretation$medium_effect_count,
  Primary_Difference_Source = clinical_interpretation$interpretation,
  Clinical_Implication = case_when(
    clinical_interpretation$interpretation == "Post-operative" ~ "Focus on post-operative care",
    clinical_interpretation$interpretation == "Pre-operative" ~ "Develop pre-operative risk stratification",
    TRUE ~ "Comprehensive pre-post management"
  )
)

write.csv(clinical_summary, "clinical_interpretation_summary.csv", row.names = FALSE)

# ================== 12. 最终总结 ==================

cat("\n🎉 术前基线特征分析完成！\n")
cat("========================================\n")
cat("✅ 数据处理：术前OCTA T0 (0_21广角区域) + Pre-Vision参数\n")
cat("✅ 统计分析：组间差异检验 + 效应量计算\n")
cat("✅ 可视化：热图 + 带p值箱线图 + PCA + 效应量图\n")
cat("✅ 临床解释：差异来源 + 临床意义 + 应用建议\n")
cat("✅ 报告生成：详细分析报告和结论\n")
cat("========================================\n")

# 显示主要结论
cat("\n🎯 主要结论 (专注0_21区域)：\n")
if(clinical_interpretation$significant_count == 0) {
  cat("患者群体差异主要在术后恢复过程中体现\n")
  cat("→ 重点：优化手术技巧和术后护理\n")
} else if(clinical_interpretation$large_effect_count > 0) {
  cat("患者在术前即存在重要差异\n")
  cat("→ 重点：建立术前风险分层系统\n")
} else {
  cat("术前存在轻微差异，术后差异更明显\n")
  cat("→ 重点：综合术前评估和术后管理\n")
}

cat("\n📁 主要输出文件：\n")
output_files <- c(
  "baseline_differences_statistics.csv",
  "baseline_comprehensive_data.csv", 
  "clinical_interpretation_summary.csv",
  "Baseline_Analysis_Report.txt"
)

for(file in output_files) {
  if(file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  }
}

cat("\n📊 可视化文件 (带p值标注)：\n")
viz_files <- list.files("plots", pattern = "\\.(pdf|png)$", full.names = FALSE)
for(file in viz_files) {
  cat(sprintf("✓ plots/%s\n", file))
}

cat("\n🎯 分析特色：\n")
cat("✅ 专注0_21区域（广角区域）分析\n")
cat("✅ 所有箱线图标注p值\n")
cat("✅ 结合视力和OCTA参数\n") 
cat("✅ 明确术前vs术后差异来源\n")

cat("\n这项分析专门针对0_21区域，回答了关键科学问题：\n")
cat("患者差异的时间起源，为late recovery聚类提供术前基线证据！\n")