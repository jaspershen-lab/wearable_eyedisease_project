# OCTA 聚类分析 - 专注于PPV组
# 使用Mfuzz模糊聚类区分预后较好和较差的患者
# 联合血流和厚度数据进行聚类

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
dir.create("3_data_analysis/6_clustering_modeling/mfuzz/octa_cluster/combined_blood_thick", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/mfuzz/octa_cluster/combined_blood_thick")

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

# -------------------- 6. 合并血流和厚度数据 --------------------
# 合并血流和厚度的改善值数据
ppv_combined_improvement <- ppv_bloodflow_improvement %>%
  full_join(ppv_thickness_improvement, by = "ID", suffix = c("_bloodflow", "_thickness"))

# 检查合并结果
cat("\n===== 合并数据统计 =====\n")
cat("血流数据患者数:", nrow(ppv_bloodflow_improvement), "\n")
cat("厚度数据患者数:", nrow(ppv_thickness_improvement), "\n")
cat("合并后患者数:", nrow(ppv_combined_improvement), "\n")
cat("血流参数数:", ncol(ppv_bloodflow_improvement) - 1, "\n")
cat("厚度参数数:", ncol(ppv_thickness_improvement) - 1, "\n")
cat("合并后总参数数:", ncol(ppv_combined_improvement) - 1, "\n\n")

# 保存原始ID列表，用于后续匹配
original_ids_combined <- ppv_combined_improvement$ID

# 去掉ID列用于聚类
ppv_combined_for_clustering <- ppv_combined_improvement %>%
  dplyr::select(-ID)

# ---------- 6.1 详细分析缺失值 ----------
# 按行计算缺失值
na_count_by_row <- rowSums(is.na(ppv_combined_for_clustering))
na_count_table <- table(na_count_by_row)

# 按列计算缺失值
na_count_by_col <- colSums(is.na(ppv_combined_for_clustering))

# 打印缺失值分析
cat("\n===== 合并OCTA数据缺失值分析 =====\n")
cat("总行数(患者数):", nrow(ppv_combined_for_clustering), "\n")
cat("总列数(参数数):", ncol(ppv_combined_for_clustering), "\n")

cat("\n患者的缺失值分布:\n")
for(na_count in names(na_count_table)) {
  cat("有", na_count, "个缺失值的患者:", na_count_table[na_count], "人\n")
}

cat("\n各参数的缺失值情况:\n")
na_summary <- data.frame(
  Parameter = names(ppv_combined_for_clustering),
  Missing_Count = na_count_by_col,
  Missing_Percent = round(na_count_by_col/nrow(ppv_combined_for_clustering)*100, 1),
  Data_Type = ifelse(grepl("bloodflow|SVP|ICP|DCP|Choroid", names(ppv_combined_for_clustering)), "BloodFlow", "Thickness")
)

# 按数据类型和缺失比例排序显示
na_summary_sorted <- na_summary %>% arrange(Data_Type, Missing_Percent)
print(na_summary_sorted)

# ---------- 6.2 筛选完整数据 ----------
# 确定完全没有NA值的行
complete_rows_combined <- complete.cases(ppv_combined_for_clustering)
complete_data_combined <- ppv_combined_for_clustering[complete_rows_combined, ]
complete_ids_combined <- original_ids_combined[complete_rows_combined]

# 打印完整数据统计
cat("\n===== 完整数据统计 =====\n")
cat("原始患者数量:", nrow(ppv_combined_for_clustering), "\n")
cat("完整数据患者数量:", nrow(complete_data_combined), 
    "(", round(nrow(complete_data_combined)/nrow(ppv_combined_for_clustering)*100, 1), "%)\n")
cat("移除了", sum(!complete_rows_combined), "位有缺失值的患者\n")

# 分别统计血流和厚度参数的数量
bloodflow_params <- names(complete_data_combined)[grepl("SVP|ICP|DCP|Choroid", names(complete_data_combined))]
thickness_params <- names(complete_data_combined)[grepl("GCL|INL|Retina", names(complete_data_combined))]

cat("其中血流参数:", length(bloodflow_params), "个\n")
cat("其中厚度参数:", length(thickness_params), "个\n\n")

# 检查最终数据集是否可用于聚类
if(nrow(complete_data_combined) < 5) {
  stop("完整数据患者数量不足5人，无法可靠聚类。建议重新考虑使用均值填充或其他策略。")
} else if(nrow(complete_data_combined) < 10) {
  cat("⚠️ 警告: 完整数据患者数量较少(", nrow(complete_data_combined), 
      "人)，聚类结果可靠性可能受影响\n\n")
}

# ---------- 6.3 标准化处理后的数据 ----------
# 对最终使用的完整数据进行标准化
ppv_combined_std <- as.data.frame(scale(complete_data_combined))

# 打印处理结果摘要
cat("===== 聚类数据准备完成 =====\n")
cat("使用的总参数数量:", ncol(ppv_combined_std), "\n")
cat("- 血流参数数量:", length(bloodflow_params), "\n")
cat("- 厚度参数数量:", length(thickness_params), "\n")
cat("使用的患者数量:", nrow(ppv_combined_std), "\n")
cat("使用的患者ID:", paste(complete_ids_combined, collapse=", "), "\n\n")

# -------------------- 7. 执行Mfuzz聚类 --------------------
# 准备用于Mfuzz的ExpressionSet对象
prep_data_for_mfuzz <- function(data, row_ids) {
  data_matrix <- as.matrix(data)
  rownames(data_matrix) <- row_ids
  
  eset <- ExpressionSet(assayData = data_matrix)
  return(eset)
}

# 创建ExpressionSet，使用完整数据行对应的ID
ppv_combined_eset <- prep_data_for_mfuzz(ppv_combined_std, complete_ids_combined)

# 为Mfuzz标准化
ppv_combined_eset_std <- standardise(ppv_combined_eset)

# 估计最佳模糊系数m
ppv_m_combined <- mestimate(ppv_combined_eset_std)
cat("估计的最佳模糊系数 m:", ppv_m_combined, "\n")

# 执行2个聚类的聚类
set.seed(123)
ppv_clustering_combined <- mfuzz(ppv_combined_eset_std, c = 2, m = ppv_m_combined)

# 获取成员度值和主要聚类分配
ppv_membership_combined <- ppv_clustering_combined$membership
ppv_main_clusters_combined <- apply(ppv_membership_combined, 1, which.max)

# 创建结果数据框
ppv_clusters_result_combined <- data.frame(
  subject_id = complete_ids_combined,
  max_cluster = ppv_main_clusters_combined,
  max_membership = apply(ppv_membership_combined, 1, max)
)

# 保存聚类结果
write.csv(ppv_clusters_result_combined, "ppv_octa_combined_cluster_results.csv", row.names = FALSE)

# -------------------- 8. 可视化聚类结果 --------------------
# 创建包含原始数据和聚类分配的组合数据集
ppv_combined_with_clusters <- ppv_combined_improvement %>%
  left_join(ppv_clusters_result_combined, by = c("ID" = "subject_id"))

# 函数：可视化聚类（增强版，区分血流和厚度）
visualize_combined_clusters <- function(data, bloodflow_params, thickness_params, group_prefix) {
  # Create plots directory
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  
  # Filter out any rows with NA in max_cluster
  data <- data %>% filter(!is.na(max_cluster))
  
  # Convert cluster to factor
  data$max_cluster <- as.factor(data$max_cluster)
  
  # All parameters
  all_params <- c(bloodflow_params, thickness_params)
  
  # Calculate mean values for each cluster
  cluster_means <- data %>%
    group_by(max_cluster) %>%
    summarise(across(all_of(all_params), mean, na.rm = TRUE))
  
  # Reshape for plotting
  plot_data <- cluster_means %>%
    pivot_longer(
      cols = all_of(all_params),
      names_to = "Parameter",
      values_to = "Mean_Improvement"
    ) %>%
    mutate(
      Data_Type = ifelse(Parameter %in% bloodflow_params, "Blood Flow", "Thickness"),
      Parameter_Clean = gsub("_improvement", "", Parameter)
    )
  
  # Create combined heatmap
  p_heatmap_combined <- ggplot(plot_data, aes(x = Parameter_Clean, y = max_cluster, fill = Mean_Improvement)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    facet_wrap(~ Data_Type, scales = "free_x", ncol = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = paste(group_prefix, "Combined OCTA Cluster Overview"),
      x = "OCTA Parameters",
      y = "Cluster",
      fill = "Mean\nImprovement"
    )
  
  # Save combined heatmap
  ggsave(paste0("plots/", group_prefix, "_combined_cluster_heatmap.pdf"), p_heatmap_combined, width = 16, height = 8)
  ggsave(paste0("plots/", group_prefix, "_combined_cluster_heatmap.png"), p_heatmap_combined, width = 16, height = 8, dpi = 300)
  
  # Create separate heatmaps for blood flow and thickness
  p_heatmap_bloodflow <- plot_data %>% 
    filter(Data_Type == "Blood Flow") %>%
    ggplot(aes(x = Parameter_Clean, y = max_cluster, fill = Mean_Improvement)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste(group_prefix, "Blood Flow Parameters by Cluster"),
      x = "Blood Flow Parameters",
      y = "Cluster",
      fill = "Mean\nImprovement"
    )
  
  p_heatmap_thickness <- plot_data %>% 
    filter(Data_Type == "Thickness") %>%
    ggplot(aes(x = Parameter_Clean, y = max_cluster, fill = Mean_Improvement)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste(group_prefix, "Thickness Parameters by Cluster"),
      x = "Thickness Parameters", 
      y = "Cluster",
      fill = "Mean\nImprovement"
    )
  
  # Save separate heatmaps
  ggsave(paste0("plots/", group_prefix, "_bloodflow_cluster_heatmap.pdf"), p_heatmap_bloodflow, width = 12, height = 6)
  ggsave(paste0("plots/", group_prefix, "_thickness_cluster_heatmap.pdf"), p_heatmap_thickness, width = 10, height = 6)
  
  # Create PCA visualization
  if(length(all_params) > 2) {
    # Perform PCA
    pca_data <- data %>%
      dplyr::select(all_of(all_params))
    
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
        title = paste(group_prefix, "PCA of Combined OCTA Improvements"),
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
      )
    
    # Save PCA plot
    ggsave(paste0("plots/", group_prefix, "_combined_pca_plot.pdf"), p_pca, width = 10, height = 8)
    ggsave(paste0("plots/", group_prefix, "_combined_pca_plot.png"), p_pca, width = 10, height = 8, dpi = 300)
    
    # Create biplot showing variable contributions
    # Extract variable loadings
    loadings <- pca_result$rotation[, 1:2]
    loadings_df <- data.frame(
      Variable = rownames(loadings),
      PC1 = loadings[, 1],
      PC2 = loadings[, 2],
      Data_Type = ifelse(rownames(loadings) %in% bloodflow_params, "Blood Flow", "Thickness")
    )
    
    # Create biplot
    p_biplot <- ggplot() +
      geom_point(data = pca_plot_data, aes(x = PC1, y = PC2, color = Cluster, alpha = Membership), size = 3) +
      geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3, linetype = Data_Type), 
                   arrow = arrow(length = unit(0.1, "cm")), alpha = 0.7) +
      geom_text(data = loadings_df, aes(x = PC1*3.2, y = PC2*3.2, label = gsub("_improvement", "", Variable)), 
                size = 2.5, alpha = 0.8) +
      scale_color_brewer(palette = "Set1") +
      theme_bw() +
      labs(
        title = paste(group_prefix, "PCA Biplot - Combined OCTA Data"),
        x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
        y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
      )
    
    # Save biplot
    ggsave(paste0("plots/", group_prefix, "_combined_biplot.pdf"), p_biplot, width = 14, height = 10)
  }
}

# 获取血流和厚度参数列表
improvement_bloodflow_params <- names(ppv_combined_with_clusters)[grepl("SVP|ICP|DCP|Choroid", names(ppv_combined_with_clusters)) & 
                                                                    grepl("improvement", names(ppv_combined_with_clusters))]
improvement_thickness_params <- names(ppv_combined_with_clusters)[grepl("GCL|INL|Retina", names(ppv_combined_with_clusters)) & 
                                                                    grepl("improvement", names(ppv_combined_with_clusters))]

# 创建可视化
visualize_combined_clusters(ppv_combined_with_clusters, improvement_bloodflow_params, improvement_thickness_params, "PPV_Combined")

# -------------------- 9. 聚类差异的统计分析 --------------------
# 函数：执行聚类间的统计分析（增强版）
analyze_combined_cluster_differences <- function(data, bloodflow_params, thickness_params) {
  all_params <- c(bloodflow_params, thickness_params)
  
  results <- data.frame(
    Parameter = character(),
    Data_Type = character(),
    Cluster1_Mean = numeric(),
    Cluster2_Mean = numeric(),
    Mean_Difference = numeric(),
    P_Value = numeric(),
    Significant = character(),
    stringsAsFactors = FALSE
  )
  
  for(param in all_params) {
    # 确定数据类型
    data_type <- ifelse(param %in% bloodflow_params, "Blood Flow", "Thickness")
    
    # 提取此参数的数据
    param_data <- data[, c("max_cluster", param)]
    
    # 计算每个聚类的均值
    means <- tapply(param_data[[param]], param_data$max_cluster, mean, na.rm = TRUE)
    
    # 执行t检验
    test_result <- t.test(reformulate("max_cluster", param), data = param_data)
    
    # 添加到结果
    results <- rbind(results, data.frame(
      Parameter = gsub("_improvement", "", param),
      Data_Type = data_type,
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
  
  # 按数据类型和p值排序
  results <- results %>% arrange(Data_Type, P_Value)
  
  return(results)
}

# 执行统计分析
ppv_combined_stats <- analyze_combined_cluster_differences(ppv_combined_with_clusters, improvement_bloodflow_params, improvement_thickness_params)

# 保存统计结果
write.csv(ppv_combined_stats, "ppv_octa_combined_cluster_statistics.csv", row.names = FALSE)

# 分别显示血流和厚度的显著参数
cat("\n===== 区分聚类的显著OCTA参数 =====\n")
cat("\n血流参数:\n")
bloodflow_sig <- ppv_combined_stats %>% 
  filter(Data_Type == "Blood Flow" & Significant == "Yes") %>%
  arrange(P_Value)
print(bloodflow_sig)

cat("\n厚度参数:\n")
thickness_sig <- ppv_combined_stats %>% 
  filter(Data_Type == "Thickness" & Significant == "Yes") %>%
  arrange(P_Value)
print(thickness_sig)

# -------------------- 10. 解释聚类 --------------------
# 确定哪个聚类整体上有更好的结果
interpret_combined_clusters <- function(stats_results) {
  # 分别分析血流和厚度数据
  bloodflow_stats <- stats_results %>% filter(Data_Type == "Blood Flow")
  thickness_stats <- stats_results %>% filter(Data_Type == "Thickness")
  
  # 统计正负差异数量
  bf_positive_diffs <- sum(bloodflow_stats$Mean_Difference > 0)
  bf_negative_diffs <- sum(bloodflow_stats$Mean_Difference < 0)
  th_positive_diffs <- sum(thickness_stats$Mean_Difference > 0)
  th_negative_diffs <- sum(thickness_stats$Mean_Difference < 0)
  
  # 统计显著的正负差异数量
  bf_sig_positive <- sum(bloodflow_stats$Mean_Difference > 0 & bloodflow_stats$Significant == "Yes")
  bf_sig_negative <- sum(bloodflow_stats$Mean_Difference < 0 & bloodflow_stats$Significant == "Yes")
  th_sig_positive <- sum(thickness_stats$Mean_Difference > 0 & thickness_stats$Significant == "Yes")
  th_sig_negative <- sum(thickness_stats$Mean_Difference < 0 & thickness_stats$Significant == "Yes")
  
  # 综合所有数据
  total_sig_positive <- bf_sig_positive + th_sig_positive
  total_sig_negative <- bf_sig_negative + th_sig_negative
  total_positive <- bf_positive_diffs + th_positive_diffs
  total_negative <- bf_negative_diffs + th_negative_diffs
  
  # 确定哪个聚类"更好"（改善值更高）
  if(total_sig_positive > total_sig_negative) {
    better_cluster <- 2
    worse_cluster <- 1
  } else if(total_sig_negative > total_sig_positive) {
    better_cluster <- 1
    worse_cluster <- 2
  } else {
    # 如果在显著参数上平局，使用所有参数
    if(total_positive > total_negative) {
      better_cluster <- 2
      worse_cluster <- 1
    } else {
      better_cluster <- 1
      worse_cluster <- 2
    }
  }
  
  cat("\n===== 聚类解释 =====\n")
  cat("聚类", better_cluster, "显示出更好的整体结果（改善值更高）\n")
  cat("聚类", worse_cluster, "显示出较差的整体结果（改善值更低）\n\n")
  
  cat("详细分析:\n")
  cat("血流参数 - 显著改善更多的聚类2:", bf_sig_positive, "个, 聚类1:", bf_sig_negative, "个\n")
  cat("厚度参数 - 显著改善更多的聚类2:", th_sig_positive, "个, 聚类1:", th_sig_negative, "个\n")
  cat("综合 - 显著改善更多的聚类2:", total_sig_positive, "个, 聚类1:", total_sig_negative, "个\n")
  
  return(list(
    better_cluster = better_cluster,
    worse_cluster = worse_cluster,
    bloodflow_sig_positive = bf_sig_positive,
    bloodflow_sig_negative = bf_sig_negative,
    thickness_sig_positive = th_sig_positive,
    thickness_sig_negative = th_sig_negative,
    total_sig_positive = total_sig_positive,
    total_sig_negative = total_sig_negative
  ))
}

# 解释聚类
combined_cluster_interpretation <- interpret_combined_clusters(ppv_combined_stats)

# 添加结果标签到聚类结果
ppv_clusters_result_combined$outcome_quality <- ifelse(
  ppv_clusters_result_combined$max_cluster == combined_cluster_interpretation$better_cluster,
  "Better",
  "Worse"
)

# 保存带有结果标签的最终结果
write.csv(ppv_clusters_result_combined, "ppv_octa_combined_cluster_results_with_outcomes.csv", row.names = FALSE)

# -------------------- 11. 每个聚类的汇总统计 --------------------
# 计算每个聚类的汇总统计
combined_cluster_summary <- ppv_combined_with_clusters %>%
  group_by(max_cluster) %>%
  summarise(
    Count = n(),
    Mean_Membership = mean(max_membership, na.rm = TRUE),
    Min_Membership = min(max_membership, na.rm = TRUE),
    Max_Membership = max(max_membership, na.rm = TRUE)
  ) %>%
  mutate(
    Outcome_Quality = ifelse(max_cluster == combined_cluster_interpretation$better_cluster, "Better", "Worse")
  )

# 打印聚类汇总
cat("\n===== 聚类汇总 =====\n")
print(combined_cluster_summary)

# 保存聚类汇总
write.csv(combined_cluster_summary, "ppv_octa_combined_cluster_summary.csv", row.names = FALSE)

# -------------------- 12. 创建参数重要性分析 --------------------
# 分析哪些参数对聚类分离贡献最大
analyze_parameter_importance <- function(stats_results) {
  # 计算参数重要性分数（基于p值和效应大小）
  importance_scores <- stats_results %>%
    mutate(
      Effect_Size = abs(Mean_Difference),
      Importance_Score = Effect_Size * (-log10(P_Value + 1e-10)) # 避免log(0)
    ) %>%
    arrange(desc(Importance_Score))
  
  # 分别获取血流和厚度的Top参数
  top_bloodflow <- importance_scores %>%
    filter(Data_Type == "Blood Flow") %>%
    head(5)
  
  top_thickness <- importance_scores %>%
    filter(Data_Type == "Thickness") %>%
    head(5)
  
  cat("\n===== 参数重要性分析 =====\n")
  cat("\nTop 5 血流参数 (按重要性分数排序):\n")
  print(top_bloodflow %>% dplyr::select(Parameter, Mean_Difference, P_Value, Importance_Score, Significant))
  
  cat("\nTop 5 厚度参数 (按重要性分数排序):\n")
  print(top_thickness %>% dplyr::select(Parameter, Mean_Difference, P_Value, Importance_Score, Significant))
  
  # 保存重要性分析结果
  write.csv(importance_scores, "ppv_octa_combined_parameter_importance.csv", row.names = FALSE)
  
  return(importance_scores)
}

# 执行参数重要性分析
parameter_importance <- analyze_parameter_importance(ppv_combined_stats)

# -------------------- 13. 创建综合报告 --------------------
# 生成分析报告
generate_combined_report <- function() {
  report <- paste0(
    "========================================\n",
    "PPV组OCTA联合聚类分析报告\n",
    "========================================\n\n",
    "1. 数据概览:\n",
    "   - 分析患者数: ", nrow(complete_data_combined), "\n",
    "   - 血流参数数: ", length(improvement_bloodflow_params), "\n",
    "   - 厚度参数数: ", length(improvement_thickness_params), "\n",
    "   - 总参数数: ", ncol(complete_data_combined), "\n\n",
    
    "2. 聚类结果:\n",
    "   - 聚类数: 2\n",
    "   - 模糊系数m: ", round(ppv_m_combined, 3), "\n",
    "   - 预后较好组 (聚类", combined_cluster_interpretation$better_cluster, "): ", 
    sum(ppv_clusters_result_combined$max_cluster == combined_cluster_interpretation$better_cluster), "人\n",
    "   - 预后较差组 (聚类", combined_cluster_interpretation$worse_cluster, "): ", 
    sum(ppv_clusters_result_combined$max_cluster == combined_cluster_interpretation$worse_cluster), "人\n\n",
    
    "3. 统计学差异:\n",
    "   - 显著血流参数: ", sum(ppv_combined_stats$Data_Type == "Blood Flow" & ppv_combined_stats$Significant == "Yes"), "\n",
    "   - 显著厚度参数: ", sum(ppv_combined_stats$Data_Type == "Thickness" & ppv_combined_stats$Significant == "Yes"), "\n",
    "   - 总显著参数: ", sum(ppv_combined_stats$Significant == "Yes"), "\n\n",
    
    "4. 主要发现:\n",
    "   - 血流参数中，预后较好组显著改善更多的参数: ", combined_cluster_interpretation$bloodflow_sig_positive, "个\n",
    "   - 厚度参数中，预后较好组显著改善更多的参数: ", combined_cluster_interpretation$thickness_sig_positive, "个\n",
    "   - 联合分析提供了比单独分析更全面的患者分层信息\n\n",
    
    "5. 输出文件:\n",
    "   - ppv_octa_combined_cluster_results_with_outcomes.csv: 聚类结果\n",
    "   - ppv_octa_combined_cluster_statistics.csv: 统计分析\n",
    "   - ppv_octa_combined_parameter_importance.csv: 参数重要性\n",
    "   - plots/: 可视化结果\n\n"
  )
  
  # 保存报告
  writeLines(report, "PPV_OCTA_Combined_Clustering_Report.txt")
  cat(report)
}

# 生成报告
generate_combined_report()

# -------------------- 14. 比较单独分析与联合分析 --------------------
# 如果需要比较，可以同时运行血流和厚度的单独聚类
compare_approaches <- function() {
  cat("\n===== 分析方法比较建议 =====\n")
  cat("1. 联合分析的优势:\n")
  cat("   - 综合考虑血流和厚度信息，更全面\n")
  cat("   - 可能发现单独分析遗漏的患者亚型\n")
  cat("   - 提供更稳健的患者分层\n\n")
  
  cat("2. 单独分析的优势:\n")
  cat("   - 可以分别评估血流和厚度的独立预测价值\n")
  cat("   - 更容易解释具体机制\n")
  cat("   - 数据完整性要求相对较低\n\n")
  
  cat("3. 建议:\n")
  cat("   - 当前联合分析使用", nrow(complete_data_combined), "例完整数据\n")
  cat("   - 建议同时进行联合分析和单独分析进行对比\n")
  cat("   - 评估哪种方法的临床预测效果更好\n\n")
}

# 执行比较分析
compare_approaches()

# 打印最终消息
cat("========================================\n")
cat("PPV组OCTA联合聚类分析完成！\n")
cat("========================================\n")
cat("聚类结果已保存至'ppv_octa_combined_cluster_results_with_outcomes.csv'\n")
cat("统计分析已保存至'ppv_octa_combined_cluster_statistics.csv'\n")
cat("参数重要性已保存至'ppv_octa_combined_parameter_importance.csv'\n")
cat("可视化结果已保存至'plots'目录\n")
cat("分析报告已保存至'PPV_OCTA_Combined_Clustering_Report.txt'\n")
cat("========================================\n")

