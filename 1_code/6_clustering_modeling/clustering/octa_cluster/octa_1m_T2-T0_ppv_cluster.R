# OCTA 聚类分析 - 专注于PPV组
# 使用Mfuzz模糊聚类区分预后较好和较差的患者

# 加载必要的库
library(tidyverse)
library(Biobase)
library(Mfuzz)
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(corrplot)
library(r4projects)

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# -------------------- 1. 加载数据 --------------------
# 加载基线信息和OCTA数据
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")

# 创建输出目录
dir.create("3_data_analysis/6_clustering_modeling/mfuzz/octa_cluster", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/mfuzz/octa_cluster")

# -------------------- 2. 处理OCTA数据 --------------------
# 函数：处理OCTA数据
process_octa_data <- function(baseline_data, octa_data, id_column = "id") {
  features <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  return(features)
}

# 处理血流和厚度数据
octa_bloodflow_features <- process_octa_data(baseline_info, octa_bloodflow)
octa_thickness_features <- process_octa_data(baseline_info, octa_thickness)

# 函数：处理每个患者数据
process_patient_data <- function(patient_data, time_points = c("T0", "T2")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术用左眼数据，右眼和双眼手术用右眼数据
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # 处理每个时间点
  for(suffix in time_points) {
    # 选择当前时间点的列
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # 选择数据并重命名列
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# 函数：处理所有患者的数据
process_all_patients <- function(features_data) {
  patient_list <- split(features_data, features_data$ID)
  processed_data <- purrr::map(patient_list, process_patient_data)
  return(bind_rows(processed_data))
}

# 分别处理血流和厚度数据
octa_bloodflow_processed <- process_all_patients(octa_bloodflow_features)
octa_thickness_processed <- process_all_patients(octa_thickness_features)

# -------------------- 3. 提取PPV组数据 --------------------
# 获取PPV患者ID (手术类型为1的患者)
ppv_patients <- baseline_info %>%
  filter(surgery_1..0.PI.1.other. == 1) %>%
  distinct(ID) %>%
  pull(ID)

# 筛选PPV患者的OCTA数据
ppv_bloodflow <- octa_bloodflow_processed %>%
  filter(ID %in% ppv_patients)

ppv_thickness <- octa_thickness_processed %>%
  filter(ID %in% ppv_patients)

# -------------------- 4. 筛选关注的参数层 --------------------
# 函数：筛选血流参数，重点关注SVP, ICP, DCP, Choroid层
filter_bloodflow_layers <- function(data) {
  # 定义要关注的血流层
  layers_of_interest <- c("SVP", "ICP", "DCP", "Choroid")
  
  # 构建匹配模式
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*0_21_T0$")
  
  # 获取T0时间点的参数
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("未找到指定血流层的T0参数！")
    return(list(data = data, params = character(0)))
  } else {
    cat("找到", length(params_T0), "个目标血流层的T0参数:\n")
    cat(paste(params_T0, collapse = "\n"), "\n")
  }
  
  # 获取对应的T2参数
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  if(length(params_T2) < length(params_T0)) {
    missing_params <- gsub("_T0$", "_T2", params_T0[!(gsub("_T0$", "_T2", params_T0) %in% params_T2)])
    warning("缺失的T2参数: ", paste(missing_params, collapse = ", "))
  }
  
  # 只保留同时有T0和T2值的参数
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  # 最终参数列表
  final_params_T0 <- paste0(valid_base_params, "_T0")
  final_params_T2 <- paste0(valid_base_params, "_T2")
  
  # 筛选数据
  filtered_data <- data %>%
    dplyr::select(ID, all_of(c(final_params_T0, final_params_T2)))
  
  return(list(
    data = filtered_data,
    params_T0 = final_params_T0,
    params_T2 = final_params_T2,
    base_params = valid_base_params
  ))
}

# 函数：筛选厚度参数，重点关注GCL.IPL, INL, Retina层
filter_thickness_layers <- function(data) {
  # 定义要关注的厚度层
  layers_of_interest <- c("GCL.IPL", "INL", "Retina")
  
  # 构建匹配模式
  pattern <- paste0("(", paste(layers_of_interest, collapse = "|"), ").*0_21_T0$")
  
  # 获取T0时间点的参数
  params_T0 <- names(data)[grep(pattern, names(data))]
  
  if(length(params_T0) == 0) {
    warning("未找到指定厚度层的T0参数！")
    return(list(data = data, params = character(0)))
  } else {
    cat("找到", length(params_T0), "个目标厚度层的T0参数:\n")
    cat(paste(params_T0, collapse = "\n"), "\n")
  }
  
  # 获取对应的T2参数
  params_T2 <- gsub("_T0$", "_T2", params_T0)
  params_T2 <- params_T2[params_T2 %in% names(data)]
  
  if(length(params_T2) < length(params_T0)) {
    missing_params <- gsub("_T0$", "_T2", params_T0[!(gsub("_T0$", "_T2", params_T0) %in% params_T2)])
    warning("缺失的T2参数: ", paste(missing_params, collapse = ", "))
  }
  
  # 只保留同时有T0和T2值的参数
  valid_base_params <- gsub("_T0$", "", params_T0[gsub("_T0$", "_T2", params_T0) %in% params_T2])
  
  # 最终参数列表
  final_params_T0 <- paste0(valid_base_params, "_T0")
  final_params_T2 <- paste0(valid_base_params, "_T2")
  
  # 筛选数据
  filtered_data <- data %>%
    dplyr::select(ID, all_of(c(final_params_T0, final_params_T2)))
  
  return(list(
    data = filtered_data,
    params_T0 = final_params_T0,
    params_T2 = final_params_T2,
    base_params = valid_base_params
  ))
}

# 应用筛选器
ppv_bloodflow_filtered <- filter_bloodflow_layers(ppv_bloodflow)
ppv_thickness_filtered <- filter_thickness_layers(ppv_thickness)

# -------------------- 5. 计算改善值(T2-T0) --------------------
# 这将作为聚类的基础

calculate_improvement <- function(data, params_T0, params_T2) {
  # 初始化结果数据框
  result <- data %>% dplyr::select(ID)
  
  # 计算每个参数的改善值
  for(i in 1:length(params_T0)) {
    t0_param <- params_T0[i]
    t2_param <- params_T2[i]
    
    # 基础参数名
    base_param <- gsub("_T0$", "", t0_param)
    
    # 计算改善值
    result[[paste0(base_param, "_improvement")]] <- data[[t2_param]] - data[[t0_param]]
  }
  
  return(result)
}

# 计算改善值
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

# -------------------- 6. 准备聚类数据 --------------------
# 使用完全没有NA值的数据进行聚类

# 保存原始ID列表，用于后续匹配
original_ids <- ppv_bloodflow_improvement$ID

# 去掉ID列用于聚类
ppv_bloodflow_for_clustering <- ppv_bloodflow_improvement %>%
  dplyr::select(-ID)

# ---------- 6.1 详细分析缺失值 ----------
# 按行计算缺失值
na_count_by_row <- rowSums(is.na(ppv_bloodflow_for_clustering))
na_count_table <- table(na_count_by_row)

# 按列计算缺失值
na_count_by_col <- colSums(is.na(ppv_bloodflow_for_clustering))

# 打印缺失值分析
cat("\n===== OCTA数据缺失值分析 =====\n")
cat("总行数(患者数):", nrow(ppv_bloodflow_for_clustering), "\n")
cat("总列数(参数数):", ncol(ppv_bloodflow_for_clustering), "\n")

cat("\n患者的缺失值分布:\n")
for(na_count in names(na_count_table)) {
  cat("有", na_count, "个缺失值的患者:", na_count_table[na_count], "人\n")
}

cat("\n各参数的缺失值情况:\n")
for(i in 1:length(na_count_by_col)) {
  cat(names(ppv_bloodflow_for_clustering)[i], ": ", 
      na_count_by_col[i], "个缺失值 (", 
      round(na_count_by_col[i]/nrow(ppv_bloodflow_for_clustering)*100, 1), "%)\n")
}

# ---------- 6.2 筛选完整数据 ----------
# 确定完全没有NA值的行
complete_rows <- complete.cases(ppv_bloodflow_for_clustering)
complete_data <- ppv_bloodflow_for_clustering[complete_rows, ]
complete_ids <- original_ids[complete_rows]

# 打印完整数据统计
cat("\n===== 完整数据统计 =====\n")
cat("原始患者数量:", nrow(ppv_bloodflow_for_clustering), "\n")
cat("完整数据患者数量:", nrow(complete_data), 
    "(", round(nrow(complete_data)/nrow(ppv_bloodflow_for_clustering)*100, 1), "%)\n")
cat("移除了", sum(!complete_rows), "位有缺失值的患者\n\n")

# 检查最终数据集是否可用于聚类
if(nrow(complete_data) < 5) {
  stop("完整数据患者数量不足5人，无法可靠聚类。建议重新考虑使用均值填充或其他策略。")
} else if(nrow(complete_data) < 10) {
  cat("⚠️ 警告: 完整数据患者数量较少(", nrow(complete_data), 
      "人)，聚类结果可靠性可能受影响\n\n")
}

# ---------- 6.3 标准化处理后的数据 ----------
# 对最终使用的完整数据进行标准化
ppv_bloodflow_std <- as.data.frame(scale(complete_data))

# 打印处理结果摘要
cat("===== 聚类数据准备完成 =====\n")
cat("使用的参数数量:", ncol(ppv_bloodflow_std), "\n")
cat("使用的参数:", paste(names(ppv_bloodflow_std), collapse=", "), "\n")
cat("使用的患者数量:", nrow(ppv_bloodflow_std), "\n")
cat("使用的患者ID:", paste(complete_ids, collapse=", "), "\n\n")


# -------------------- 7. 执行Mfuzz聚类 --------------------
# 准备用于Mfuzz的ExpressionSet对象
prep_data_for_mfuzz <- function(data, row_ids) {
  data_matrix <- as.matrix(data)
  rownames(data_matrix) <- row_ids
  
  eset <- ExpressionSet(assayData = data_matrix)
  return(eset)
}

# 创建ExpressionSet，使用完整数据行对应的ID
ppv_bloodflow_eset <- prep_data_for_mfuzz(ppv_bloodflow_std, complete_ids)

# 为Mfuzz标准化
ppv_bloodflow_eset_std <- standardise(ppv_bloodflow_eset)

# 估计最佳模糊系数m
ppv_m <- mestimate(ppv_bloodflow_eset_std)
cat("估计的最佳模糊系数 m:", ppv_m, "\n")

# 执行2个聚类的聚类
set.seed(123)
ppv_clustering <- mfuzz(ppv_bloodflow_eset_std, c = 2, m = ppv_m)

# 获取成员度值和主要聚类分配
ppv_membership <- ppv_clustering$membership
ppv_main_clusters <- apply(ppv_membership, 1, which.max)

# 创建结果数据框
ppv_clusters_result <- data.frame(
  subject_id = complete_ids,
  max_cluster = ppv_main_clusters,
  max_membership = apply(ppv_membership, 1, max)
)

# 保存聚类结果
write.csv(ppv_clusters_result, "ppv_octa_cluster_results.csv", row.names = FALSE)


# -------------------- 8. 可视化聚类结果 --------------------
# 创建包含原始数据和聚类分配的组合数据集
ppv_bloodflow_with_clusters <- ppv_bloodflow_improvement %>%
  left_join(ppv_clusters_result, by = c("ID" = "subject_id"))

# 函数：可视化聚类
visualize_clusters <- function(data, params, group_prefix) {
  # Create plots directory
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  
  # Filter out any rows with NA in max_cluster
  data <- data %>% filter(!is.na(max_cluster))
  
  # Convert cluster to factor
  data$max_cluster <- as.factor(data$max_cluster)
  
  # Calculate mean values for each cluster
  cluster_means <- data %>%
    group_by(max_cluster) %>%
    summarise(across(all_of(params), mean, na.rm = TRUE))
  
  # Reshape for plotting
  plot_data <- cluster_means %>%
    pivot_longer(
      cols = all_of(params),
      names_to = "Parameter",
      values_to = "Mean_Improvement"
    )
  
  # Create heatmap with English labels
  p_heatmap <- ggplot(plot_data, aes(x = Parameter, y = max_cluster, fill = Mean_Improvement)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste(group_prefix, "Cluster Overview"),
      x = "OCTA Parameters",
      y = "Cluster",
      fill = "Mean\nImprovement"
    )
  
  # Save heatmap
  ggsave(paste0("plots/", group_prefix, "_cluster_heatmap.pdf"), p_heatmap, width = 12, height = 6)
  ggsave(paste0("plots/", group_prefix, "_cluster_heatmap.png"), p_heatmap, width = 12, height = 6, dpi = 300)
  
  # Create boxplots for each parameter with English labels
  for(param in params) {
    # Create readable parameter name
    param_name <- gsub("_improvement", "", param)
    
    # Create boxplot
    p_box <- ggplot(data, aes(x = max_cluster, y = .data[[param]], fill = max_cluster)) +
      geom_boxplot() +
      geom_jitter(alpha = 0.3, width = 0.2) +
      scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      labs(
        title = paste(param_name, "Improvement by Cluster"),
        x = "Cluster",
        y = "Improvement (T2-T0)"
      )
    
    # Save boxplot
    ggsave(paste0("plots/", group_prefix, "_", param_name, "_boxplot.pdf"), p_box, width = 8, height = 6)
  }
  
  # Create PCA visualization with English labels
  if(length(params) > 2) {
    # Perform PCA
    pca_data <- data %>%
      dplyr::select(all_of(params))
    
    pca_result <- prcomp(pca_data, scale. = TRUE)
    
    # Prepare plot data
    pca_plot_data <- data.frame(
      PC1 = pca_result$x[,1],
      PC2 = pca_result$x[,2],
      Cluster = data$max_cluster,
      Membership = data$max_membership
    )
    
    # Create PCA plot
    p_pca <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Cluster, alpha = Membership)) +
      geom_point(size = 3) +
      stat_ellipse(aes(group = Cluster), level = 0.95) +
      scale_color_brewer(palette = "Set1") +
      theme_bw() +
      labs(
        title = paste(group_prefix, "PCA of OCTA Improvements"),
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
      )
    
    # Save PCA plot
    ggsave(paste0("plots/", group_prefix, "_pca_plot.pdf"), p_pca, width = 10, height = 8)
    ggsave(paste0("plots/", group_prefix, "_pca_plot.png"), p_pca, width = 10, height = 8, dpi = 300)
  }
}

# 获取用于可视化的参数
improvement_params <- names(ppv_bloodflow_with_clusters)[!names(ppv_bloodflow_with_clusters) %in% 
                                                           c("ID", "max_cluster", "max_membership")]

# 创建可视化
visualize_clusters(ppv_bloodflow_with_clusters, improvement_params, "PPV_BloodFlow")

# -------------------- 9. 聚类差异的统计分析 --------------------
# 函数：执行聚类间的统计分析
analyze_cluster_differences <- function(data, params) {
  results <- data.frame(
    Parameter = character(),
    Cluster1_Mean = numeric(),
    Cluster2_Mean = numeric(),
    Mean_Difference = numeric(),
    P_Value = numeric(),
    Significant = character(),
    stringsAsFactors = FALSE
  )
  
  for(param in params) {
    # 提取此参数的数据
    param_data <- data[, c("max_cluster", param)]
    
    # 计算每个聚类的均值
    means <- tapply(param_data[[param]], param_data$max_cluster, mean, na.rm = TRUE)
    
    # 执行t检验
    test_result <- t.test(reformulate("max_cluster", param), data = param_data)
    
    # 添加到结果
    results <- rbind(results, data.frame(
      Parameter = gsub("_improvement", "", param),
      Cluster1_Mean = means["1"],
      Cluster2_Mean = means["2"],
      Mean_Difference = means["2"] - means["1"],
      P_Value = test_result$p.value,
      Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
  
  # 添加校正后的p值
  results$P_Adjusted <- p.adjust(results$P_Value, method = "fdr")
  results$Significant_Adjusted <- ifelse(results$P_Adjusted < 0.05, "Yes", "No")
  
  # 按p值排序
  results <- results %>% arrange(P_Value)
  
  return(results)
}

# 执行统计分析
ppv_stats <- analyze_cluster_differences(ppv_bloodflow_with_clusters, improvement_params)

# 保存统计结果
write.csv(ppv_stats, "ppv_octa_cluster_statistics.csv", row.names = FALSE)

# 打印显著的参数
cat("\n区分聚类的显著OCTA参数:\n")
print(head(ppv_stats[ppv_stats$Significant == "Yes", ], 10))

# -------------------- 10. 解释聚类 --------------------
# 确定哪个聚类整体上有更好的结果
interpret_clusters <- function(stats_results) {
  # 统计正负差异数量
  positive_diffs <- sum(stats_results$Mean_Difference > 0)
  negative_diffs <- sum(stats_results$Mean_Difference < 0)
  
  # 统计显著的正负差异数量
  sig_positive_diffs <- sum(stats_results$Mean_Difference > 0 & stats_results$Significant == "Yes")
  sig_negative_diffs <- sum(stats_results$Mean_Difference < 0 & stats_results$Significant == "Yes")
  
  # 确定哪个聚类"更好"（改善值更高）
  if(sig_positive_diffs > sig_negative_diffs) {
    better_cluster <- 2
    worse_cluster <- 1
  } else if(sig_negative_diffs > sig_positive_diffs) {
    better_cluster <- 1
    worse_cluster <- 2
  } else {
    # 如果在显著参数上平局，使用所有参数
    if(positive_diffs > negative_diffs) {
      better_cluster <- 2
      worse_cluster <- 1
    } else {
      better_cluster <- 1
      worse_cluster <- 2
    }
  }
  
  cat("\n聚类解释:\n")
  cat("聚类", better_cluster, "显示出更好的整体结果（改善值更高）\n")
  cat("聚类", worse_cluster, "显示出较差的整体结果（改善值更低）\n")
  
  return(list(
    better_cluster = better_cluster,
    worse_cluster = worse_cluster,
    sig_positive_diffs = sig_positive_diffs,
    sig_negative_diffs = sig_negative_diffs
  ))
}

# 解释聚类
cluster_interpretation <- interpret_clusters(ppv_stats)

# 添加结果标签到聚类结果
ppv_clusters_result$outcome_quality <- ifelse(
  ppv_clusters_result$max_cluster == cluster_interpretation$better_cluster,
  "Better",
  "Worse"
)

# 保存带有结果标签的最终结果
write.csv(ppv_clusters_result, "ppv_octa_cluster_results_with_outcomes.csv", row.names = FALSE)

# -------------------- 11. 每个聚类的汇总统计 --------------------
# 计算每个聚类的汇总统计
cluster_summary <- ppv_bloodflow_with_clusters %>%
  group_by(max_cluster) %>%
  summarise(
    Count = n(),
    Mean_Membership = mean(max_membership, na.rm = TRUE),
    Min_Membership = min(max_membership, na.rm = TRUE),
    Max_Membership = max(max_membership, na.rm = TRUE)
  ) %>%
  mutate(
    Outcome_Quality = ifelse(max_cluster == cluster_interpretation$better_cluster, "Better", "Worse")
  )

# 打印聚类汇总
cat("\n聚类汇总:\n")
print(cluster_summary)

# 保存聚类汇总
write.csv(cluster_summary, "ppv_octa_cluster_summary.csv", row.names = FALSE)

# 打印最终消息
cat("\nPPV组OCTA聚类完成。\n")
cat("聚类结果已保存至'ppv_octa_cluster_results_with_outcomes.csv'\n")
cat("聚类统计已保存至'ppv_octa_cluster_statistics.csv'\n")
cat("可视化已保存至'plots'目录\n")
