# OCTA 完整聚类分析 - 传统方法 + 改进方法
# 使用K-means, 层次聚类, PAM等方法，并包含异常值处理和多种改进策略

# 加载必要的库
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(cluster)
library(corrplot)
library(NbClust)
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
dir.create("3_data_analysis/6_clustering_modeling/traditional_clustering/octa_cluster", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/traditional_clustering/octa_cluster")

# -------------------- 2. 数据处理函数 --------------------
# 处理OCTA数据
process_octa_data <- function(baseline_data, octa_data, id_column = "id") {
  features <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  return(features)
}

# 处理每个患者数据
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

# 处理所有患者的数据
process_all_patients <- function(features_data) {
  patient_list <- split(features_data, features_data$ID)
  processed_data <- purrr::map(patient_list, process_patient_data)
  return(bind_rows(processed_data))
}

# -------------------- 3. 数据预处理 --------------------
# 处理血流和厚度数据
octa_bloodflow_features <- process_octa_data(baseline_info, octa_bloodflow)
octa_thickness_features <- process_octa_data(baseline_info, octa_thickness)

# 分别处理血流和厚度数据
octa_bloodflow_processed <- process_all_patients(octa_bloodflow_features)
octa_thickness_processed <- process_all_patients(octa_thickness_features)

# 提取PPV组数据
ppv_patients <- baseline_info %>%
  filter(surgery_1..0.PI.1.other. == 1) %>%
  distinct(ID) %>%
  pull(ID)

ppv_bloodflow <- octa_bloodflow_processed %>%
  filter(ID %in% ppv_patients)

ppv_thickness <- octa_thickness_processed %>%
  filter(ID %in% ppv_patients)

cat("PPV患者数量:", length(ppv_patients), "\n")
cat("PPV血流数据行数:", nrow(ppv_bloodflow), "\n")
cat("PPV厚度数据行数:", nrow(ppv_thickness), "\n")

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
    # 简化这个复杂的嵌套表达式
    t0_to_t2 <- gsub("_T0$", "_T2", params_T0)
    missing_t2 <- t0_to_t2[!(t0_to_t2 %in% params_T2)]
    warning("缺失的T2参数: ", paste(missing_t2, collapse = ", "))
  }
  
  # 只保留同时有T0和T2值的参数
  t0_to_t2_mapping <- gsub("_T0$", "_T2", params_T0)
  valid_t0 <- params_T0[t0_to_t2_mapping %in% params_T2]
  valid_base_params <- gsub("_T0$", "", valid_t0)
  
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
cat("\n===== 应用OCTA参数筛选器 =====\n")
cat("血流数据筛选:\n")
ppv_bloodflow_filtered <- filter_bloodflow_layers(ppv_bloodflow)

cat("\n厚度数据筛选:\n")
ppv_thickness_filtered <- filter_thickness_layers(ppv_thickness)

# -------------------- 5. 计算改善值(T2-T0) --------------------
# 计算改善值的函数
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

cat("\n===== 改善值计算完成 =====\n")
cat("血流参数改善值数量:", ncol(ppv_bloodflow_improvement) - 1, "\n")
cat("厚度参数改善值数量:", ncol(ppv_thickness_improvement) - 1, "\n")

# 合并血流和厚度改善值数据
ppv_all_improvements <- ppv_bloodflow_improvement %>%
  full_join(ppv_thickness_improvement, by = "ID")

cat("合并后总改善值特征数量:", ncol(ppv_all_improvements) - 1, "\n")
cat("患者数量:", nrow(ppv_all_improvements), "\n")

# 创建聚类特征
ppv_features <- ppv_all_improvements

# -------------------- 6. 数据清洗和缺失值处理 --------------------
# 分析缺失值
original_ids <- ppv_features$ID
features_for_clustering <- ppv_features %>% dplyr::select(-ID)

na_count_by_row <- rowSums(is.na(features_for_clustering))
na_count_by_col <- colSums(is.na(features_for_clustering))

cat("\n===== 缺失值分析 =====\n")
cat("总患者数:", nrow(features_for_clustering), "\n")
cat("总特征数:", ncol(features_for_clustering), "\n")

# 只使用完整病例
complete_rows <- complete.cases(features_for_clustering)
final_data <- features_for_clustering[complete_rows, ]
final_ids <- original_ids[complete_rows]

cat("完整病例数:", nrow(final_data), "\n")
cat("排除患者数:", sum(!complete_rows), "\n")

if(sum(!complete_rows) > 0) {
  excluded_ids <- original_ids[!complete_rows]
  cat("排除的患者ID:", paste(excluded_ids, collapse = ", "), "\n")
}

# 标准化数据
final_data_scaled <- as.data.frame(scale(final_data))

# -------------------- 7. 改进聚类方法函数定义 --------------------
# 异常值检测函数
detect_and_handle_outliers <- function(data, method = "iqr") {
  cat("===== 异常值检测 =====\n")
  
  if(method == "iqr") {
    outlier_scores <- apply(data, 1, function(row) {
      feature_scores <- sapply(1:ncol(data), function(i) {
        q1 <- quantile(data[,i], 0.25, na.rm = TRUE)
        q3 <- quantile(data[,i], 0.75, na.rm = TRUE)
        iqr <- q3 - q1
        lower <- q1 - 1.5 * iqr
        upper <- q3 + 1.5 * iqr
        
        if(row[i] < lower || row[i] > upper) {
          abs(row[i] - median(data[,i], na.rm = TRUE)) / iqr
        } else {
          0
        }
      })
      sum(feature_scores)
    })
  }
  
  threshold <- quantile(outlier_scores, 0.9)
  outliers <- which(outlier_scores > threshold)
  
  cat("检测到", length(outliers), "个异常值\n")
  if(length(outliers) > 0) {
    cat("异常值样本:", paste(rownames(data)[outliers], collapse = ", "), "\n")
  }
  
  return(list(
    outlier_indices = outliers,
    outlier_scores = outlier_scores,
    threshold = threshold
  ))
}

# 鲁棒聚类方法函数
perform_robust_clustering <- function(data_scaled, final_ids, method = "all") {
  results <- list()
  
  if(method %in% c("kmeans", "all")) {
    cat("\n===== 鲁棒K-means聚类 =====\n")
    
    for(k in 2:4) {
      set.seed(123)
      kmeans_result <- kmeans(data_scaled, centers = k, nstart = 50, iter.max = 100)
      sil_score <- mean(silhouette(kmeans_result$cluster, dist(data_scaled))[,3])
      
      cluster_sizes <- table(kmeans_result$cluster)
      min_cluster_size <- min(cluster_sizes)
      balance_score <- min_cluster_size / max(cluster_sizes)
      
      cat(sprintf("K=%d: Silhouette=%.3f, 最小聚类大小=%d, 平衡度=%.3f\n", 
                  k, sil_score, min_cluster_size, balance_score))
      
      results[[paste0("kmeans_k", k)]] <- list(
        method = paste0("K-means (K=", k, ")"),
        clustering = kmeans_result,
        k = k,
        silhouette = sil_score,
        balance = balance_score,
        min_size = min_cluster_size
      )
    }
  }
  
  if(method %in% c("pam", "all")) {
    cat("\n===== PAM聚类 =====\n")
    
    for(k in 2:4) {
      pam_result <- pam(data_scaled, k = k, metric = "euclidean")
      sil_score <- pam_result$silinfo$avg.width
      
      cluster_sizes <- table(pam_result$clustering)
      min_cluster_size <- min(cluster_sizes)
      balance_score <- min_cluster_size / max(cluster_sizes)
      
      cat(sprintf("PAM K=%d: Silhouette=%.3f, 最小聚类大小=%d, 平衡度=%.3f\n", 
                  k, sil_score, min_cluster_size, balance_score))
      
      results[[paste0("pam_k", k)]] <- list(
        method = paste0("PAM (K=", k, ")"),
        clustering = pam_result,
        k = k,
        silhouette = sil_score,
        balance = balance_score,
        min_size = min_cluster_size
      )
    }
  }
  
  return(results)
}

# 排除异常值后聚类函数
cluster_without_outliers <- function(data_scaled, final_ids, outlier_indices) {
  cat("\n===== 排除异常值后重新聚类 =====\n")
  
  if(length(outlier_indices) == 0) {
    cat("没有检测到异常值，跳过此方法\n")
    return(list())
  }
  
  clean_data <- data_scaled[-outlier_indices, ]
  clean_ids <- final_ids[-outlier_indices]
  
  cat("排除", length(outlier_indices), "个异常值后，剩余", nrow(clean_data), "个样本\n")
  
  if(nrow(clean_data) < 4) {
    cat("排除异常值后样本过少，无法聚类\n")
    return(list())
  }
  
  results <- list()
  
  for(k in 2:3) {
    if(nrow(clean_data) >= k * 2) {
      set.seed(123)
      kmeans_result <- kmeans(clean_data, centers = k, nstart = 25)
      sil_score <- mean(silhouette(kmeans_result$cluster, dist(clean_data))[,3])
      
      cluster_sizes <- table(kmeans_result$cluster)
      min_cluster_size <- min(cluster_sizes)
      
      cat(sprintf("Clean K-means K=%d: Silhouette=%.3f, 最小聚类大小=%d\n", 
                  k, sil_score, min_cluster_size))
      
      # 为原始数据分配聚类（异常值标记为单独类别）
      full_clusters <- rep(k+1, length(final_ids))
      full_clusters[-outlier_indices] <- kmeans_result$cluster
      
      results[[paste0("clean_kmeans_k", k)]] <- list(
        method = paste0("K-means without outliers (K=", k, ")"),
        clustering = list(cluster = full_clusters),
        k = k,
        silhouette = sil_score,
        outliers_removed = length(outlier_indices),
        clean_ids = clean_ids
      )
    }
  }
  
  return(results)
}

# 评估和选择最佳方法函数
evaluate_and_select_best <- function(all_results, min_cluster_size = 2, prefer_balance = TRUE) {
  cat("\n===== 聚类方法综合评估 =====\n")
  
  if(length(all_results) == 0) {
    cat("没有可用的聚类结果\n")
    return(NULL)
  }
  
  evaluation <- data.frame(
    Method = character(),
    K = numeric(),
    Silhouette = numeric(),
    Balance = numeric(),
    Min_Size = numeric(),
    Score = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(name in names(all_results)) {
    result <- all_results[[name]]
    
    if(!is.null(result$min_size) && result$min_size < min_cluster_size) {
      next
    }
    
    sil_score <- ifelse(is.null(result$silhouette), 0, result$silhouette)
    balance_score <- ifelse(is.null(result$balance), 0.5, result$balance)
    
    if(prefer_balance) {
      composite_score <- 0.6 * sil_score + 0.4 * balance_score
    } else {
      composite_score <- 0.8 * sil_score + 0.2 * balance_score
    }
    
    evaluation <- rbind(evaluation, data.frame(
      Method = result$method,
      K = result$k,
      Silhouette = round(sil_score, 3),
      Balance = round(balance_score, 3),
      Min_Size = ifelse(is.null(result$min_size), NA, result$min_size),
      Score = round(composite_score, 3)
    ))
  }
  
  if(nrow(evaluation) == 0) {
    cat("没有符合条件的聚类结果\n")
    return(NULL)
  }
  
  evaluation <- evaluation[order(evaluation$Score, decreasing = TRUE), ]
  
  cat("聚类方法评估结果:\n")
  print(evaluation)
  
  best_method_name <- names(all_results)[sapply(names(all_results), function(x) 
    all_results[[x]]$method == evaluation$Method[1])][1]
  
  best_result <- all_results[[best_method_name]]
  
  cat("\n推荐的最佳聚类方法:", best_result$method, "\n")
  cat("K =", best_result$k, ", Silhouette =", round(best_result$silhouette, 3), "\n")
  
  return(list(
    best_result = best_result,
    evaluation = evaluation,
    best_method_name = best_method_name
  ))
}

# -------------------- 8. 执行改进的聚类分析 --------------------
cat("\n========== 开始改进的OCTA聚类分析 ==========\n")
cat("样本数:", nrow(final_data_scaled), "\n")
cat("特征数:", ncol(final_data_scaled), "\n\n")

all_results <- list()

# 步骤1: 异常值检测
outlier_detection <- detect_and_handle_outliers(final_data_scaled, method = "iqr")

# 步骤2: 鲁棒聚类方法
robust_results <- perform_robust_clustering(final_data_scaled, final_ids, method = "all")
all_results <- c(all_results, robust_results)

# 步骤3: 排除异常值后重新聚类
if(length(outlier_detection$outlier_indices) > 0 && 
   length(outlier_detection$outlier_indices) < nrow(final_data_scaled) * 0.3) {
  clean_results <- cluster_without_outliers(final_data_scaled, final_ids, 
                                            outlier_detection$outlier_indices)
  all_results <- c(all_results, clean_results)
}

# 步骤4: 评估和选择最佳方法
best_selection <- evaluate_and_select_best(all_results, min_cluster_size = 2, prefer_balance = TRUE)

# -------------------- 9. 更新最佳聚类结果 --------------------
if(!is.null(best_selection)) {
  best_result <- best_selection$best_result
  
  # 获取聚类结果
  if("cluster" %in% names(best_result$clustering)) {
    best_clusters <- best_result$clustering$cluster
  } else if("clustering" %in% names(best_result$clustering)) {
    best_clusters <- best_result$clustering$clustering
  } else {
    best_clusters <- best_result$clustering
  }
  
  best_k <- best_result$k
  best_method_name <- best_result$method
  
  cat("\n========== 最佳聚类结果 ==========\n")
  cat("最佳方法:", best_method_name, "\n")
  cat("聚类数:", best_k, "\n")
  cat("轮廓系数:", round(best_result$silhouette, 3), "\n")
  
  # 显示聚类分布
  cluster_distribution <- table(best_clusters)
  cat("聚类分布:\n")
  print(cluster_distribution)
  
  # 创建最终聚类结果数据框
  clustering_results <- data.frame(
    subject_id = final_ids,
    cluster = best_clusters,
    method = best_method_name,
    optimal_k = best_k,
    silhouette_score = best_result$silhouette
  )
  
} else {
  cat("未找到合适的聚类结果\n")
  stop("无法完成聚类分析")
}

# -------------------- 10. 可视化聚类结果 --------------------
dir.create("plots", recursive = TRUE, showWarnings = FALSE)

# PCA可视化
pca_result <- prcomp(final_data_scaled)

pca_plot_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Cluster = factor(best_clusters),
  ID = final_ids
)

p_pca <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = ID), vjust = -0.8, size = 3) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(
    title = paste("PPV OCTA Clustering Results -", best_method_name),
    subtitle = paste("K =", best_k, "| Silhouette =", round(best_result$silhouette, 3)),
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

ggsave("plots/ppv_octa_complete_clustering_pca.pdf", p_pca, width = 10, height = 8)
ggsave("plots/ppv_octa_complete_clustering_pca.png", p_pca, width = 10, height = 8, dpi = 300)
print(p_pca)

# -------------------- 11. 保存结果 --------------------
# 保存聚类结果
write.csv(clustering_results, "ppv_octa_complete_clustering_results.csv", row.names = FALSE)

# 创建用于相关性分析的结果
final_results_for_correlation <- data.frame(
  subject_id = clustering_results$subject_id,
  octa_cluster = clustering_results$cluster,
  clustering_method = best_method_name,
  optimal_k = best_k,
  silhouette_score = best_result$silhouette
)

write.csv(final_results_for_correlation, "ppv_octa_clusters_for_correlation.csv", row.names = FALSE)

# 保存聚类总结
cluster_summary <- clustering_results %>%
  group_by(cluster) %>%
  summarise(
    Count = n(),
    Patient_IDs = paste(sort(subject_id), collapse = ", "),
    .groups = "drop"
  )

cat("\n聚类总结:\n")
print(cluster_summary)

write.csv(cluster_summary, "ppv_octa_clustering_summary.csv", row.names = FALSE)

cat("\n========== OCTA完整聚类分析完成 ==========\n")
cat("最佳聚类方法:", best_method_name, "\n")
cat("最优聚类数:", best_k, "\n")
cat("最佳轮廓系数:", round(best_result$silhouette, 3), "\n")
cat("聚类结果保存到:\n")
cat("- ppv_octa_complete_clustering_results.csv\n")
cat("- ppv_octa_clusters_for_correlation.csv\n")
cat("- ppv_octa_clustering_summary.csv\n")
cat("可视化图形保存在plots目录\n")

# 手动排除指定异常值后的OCTA聚类分析
# 排除SH047和SH073后重新进行聚类

# 使用之前已经准备好的数据：final_data_scaled 和 final_ids

# -------------------- 手动排除异常值 --------------------
cat("========== 手动排除异常值 ==========\n")

# 指定要排除的患者ID
manual_outliers <- c("SH047", "SH073")
cat("手动指定排除的患者:", paste(manual_outliers, collapse = ", "), "\n")

# 找出这些患者在数据中的索引
outlier_indices <- which(final_ids %in% manual_outliers)
cat("找到的异常值索引:", paste(outlier_indices, collapse = ", "), "\n")

if(length(outlier_indices) == 0) {
  cat("警告: 未找到指定的异常值患者\n")
  cat("可用的患者ID:", paste(final_ids, collapse = ", "), "\n")
  stop("请检查患者ID是否正确")
}

# 创建排除异常值后的数据
clean_data_scaled <- final_data_scaled[-outlier_indices, ]
clean_ids <- final_ids[-outlier_indices]

cat("原始患者数:", length(final_ids), "\n")
cat("排除后患者数:", length(clean_ids), "\n")
cat("排除的患者:", paste(final_ids[outlier_indices], collapse = ", "), "\n")
cat("剩余患者:", paste(clean_ids, collapse = ", "), "\n")

# -------------------- 对清洁数据执行多种聚类方法 --------------------
cat("\n========== 对清洁数据执行聚类 ==========\n")

# 创建存储结果的列表
clean_results <- list()

# 1. K-means聚类（K=2到4）
cat("\n===== K-means聚类 =====\n")
for(k in 2:4) {
  if(length(clean_ids) >= k * 2) {  # 确保每个聚类至少有2个样本
    set.seed(123)
    kmeans_result <- kmeans(clean_data_scaled, centers = k, nstart = 25)
    
    # 计算轮廓系数
    if(k > 1) {
      sil_score <- mean(silhouette(kmeans_result$cluster, dist(clean_data_scaled))[,3])
    } else {
      sil_score <- 0
    }
    
    # 检查聚类平衡性
    cluster_sizes <- table(kmeans_result$cluster)
    min_cluster_size <- min(cluster_sizes)
    balance_score <- min_cluster_size / max(cluster_sizes)
    
    cat(sprintf("K-means K=%d: Silhouette=%.3f, 最小聚类大小=%d, 平衡度=%.3f\n", 
                k, sil_score, min_cluster_size, balance_score))
    
    clean_results[[paste0("clean_kmeans_k", k)]] <- list(
      method = paste0("Clean K-means (K=", k, ")"),
      clustering = kmeans_result,
      k = k,
      silhouette = sil_score,
      balance = balance_score,
      min_size = min_cluster_size
    )
  }
}

# 2. PAM聚类（K=2到4）
cat("\n===== PAM聚类 =====\n")
for(k in 2:4) {
  if(length(clean_ids) >= k * 2) {
    pam_result <- pam(clean_data_scaled, k = k, metric = "euclidean")
    sil_score <- pam_result$silinfo$avg.width
    
    cluster_sizes <- table(pam_result$clustering)
    min_cluster_size <- min(cluster_sizes)
    balance_score <- min_cluster_size / max(cluster_sizes)
    
    cat(sprintf("Clean PAM K=%d: Silhouette=%.3f, 最小聚类大小=%d, 平衡度=%.3f\n", 
                k, sil_score, min_cluster_size, balance_score))
    
    clean_results[[paste0("clean_pam_k", k)]] <- list(
      method = paste0("Clean PAM (K=", k, ")"),
      clustering = pam_result,
      k = k,
      silhouette = sil_score,
      balance = balance_score,
      min_size = min_cluster_size
    )
  }
}

# 3. 层次聚类
cat("\n===== 层次聚类 =====\n")
dist_matrix <- dist(clean_data_scaled)
hc_result <- hclust(dist_matrix, method = "ward.D2")

for(k in 2:4) {
  if(length(clean_ids) >= k * 2) {
    clusters <- cutree(hc_result, k = k)
    
    if(length(unique(clusters)) == k) {
      sil_score <- mean(silhouette(clusters, dist_matrix)[,3])
      
      cluster_sizes <- table(clusters)
      min_cluster_size <- min(cluster_sizes)
      balance_score <- min_cluster_size / max(cluster_sizes)
      
      cat(sprintf("Clean Hierarchical K=%d: Silhouette=%.3f, 最小聚类大小=%d, 平衡度=%.3f\n", 
                  k, sil_score, min_cluster_size, balance_score))
      
      clean_results[[paste0("clean_hc_k", k)]] <- list(
        method = paste0("Clean Hierarchical (K=", k, ")"),
        clustering = list(cluster = clusters),
        k = k,
        silhouette = sil_score,
        balance = balance_score,
        min_size = min_cluster_size
      )
    }
  }
}

# -------------------- 选择最佳聚类方法 --------------------
cat("\n========== 选择最佳聚类方法 ==========\n")

# 评估所有方法
evaluation <- data.frame(
  Method = character(),
  K = numeric(),
  Silhouette = numeric(),
  Balance = numeric(),
  Min_Size = numeric(),
  Score = numeric(),
  stringsAsFactors = FALSE
)

for(name in names(clean_results)) {
  result <- clean_results[[name]]
  
  # 跳过聚类大小过小的结果
  if(result$min_size < 2) {
    next
  }
  
  # 计算综合得分（平衡轮廓系数和聚类平衡度）
  composite_score <- 0.7 * result$silhouette + 0.3 * result$balance
  
  evaluation <- rbind(evaluation, data.frame(
    Method = result$method,
    K = result$k,
    Silhouette = round(result$silhouette, 3),
    Balance = round(result$balance, 3),
    Min_Size = result$min_size,
    Score = round(composite_score, 3)
  ))
}

# 按综合得分排序
evaluation <- evaluation[order(evaluation$Score, decreasing = TRUE), ]

cat("清洁数据聚类方法评估结果:\n")
print(evaluation)

if(nrow(evaluation) == 0) {
  stop("没有找到合适的聚类结果")
}

# 选择最佳方法
best_method_name <- names(clean_results)[sapply(names(clean_results), function(x) 
  clean_results[[x]]$method == evaluation$Method[1])][1]

best_result <- clean_results[[best_method_name]]

cat("\n推荐的最佳聚类方法:", best_result$method, "\n")
cat("K =", best_result$k, ", Silhouette =", round(best_result$silhouette, 3), "\n")

# 获取最佳聚类结果
if("cluster" %in% names(best_result$clustering)) {
  clean_clusters <- best_result$clustering$cluster
} else if("clustering" %in% names(best_result$clustering)) {
  clean_clusters <- best_result$clustering$clustering
} else {
  clean_clusters <- best_result$clustering
}

# -------------------- 为所有患者分配聚类（包括异常值）--------------------
cat("\n========== 为所有患者分配聚类标签 ==========\n")

# 创建完整的聚类分配
# 异常值分配到单独的类别（编号为最大聚类数+1）
max_cluster <- max(clean_clusters)
full_clusters <- rep(max_cluster + 1, length(final_ids))  # 异常值类别
full_clusters[-outlier_indices] <- clean_clusters

# 显示最终聚类分布
cat("最终聚类分布:\n")
final_cluster_table <- table(full_clusters)
print(final_cluster_table)

# 显示每个聚类的患者
cat("\n各聚类的患者分布:\n")
for(cluster_id in sort(unique(full_clusters))) {
  patients_in_cluster <- final_ids[full_clusters == cluster_id]
  if(cluster_id == max_cluster + 1) {
    cat("聚类", cluster_id, "(异常值):", paste(patients_in_cluster, collapse = ", "), "\n")
  } else {
    cat("聚类", cluster_id, ":", paste(patients_in_cluster, collapse = ", "), "\n")
  }
}

# -------------------- 创建最终聚类结果 --------------------
clustering_results_manual <- data.frame(
  subject_id = final_ids,
  cluster = full_clusters,
  method = paste(best_result$method, "- Manual outlier removal"),
  optimal_k = best_result$k,
  silhouette_score = best_result$silhouette,
  is_outlier = final_ids %in% manual_outliers
)

# -------------------- 可视化结果 --------------------
cat("\n========== 创建可视化 ==========\n")

# 创建PCA图（使用所有数据，但标记异常值）
pca_result_full <- prcomp(final_data_scaled)

pca_data <- data.frame(
  PC1 = pca_result_full$x[,1],
  PC2 = pca_result_full$x[,2],
  Cluster = factor(full_clusters),
  ID = final_ids,
  Is_Outlier = final_ids %in% manual_outliers
)

# 创建PCA散点图
p_pca_manual <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(aes(size = ifelse(Is_Outlier, 6, 4), 
                 alpha = ifelse(Is_Outlier, 0.6, 0.9),
                 shape = ifelse(Is_Outlier, 4, 16))) +  # 异常值用X形状
  geom_text(aes(label = ID), vjust = -0.8, size = 3) +
  stat_ellipse(data = subset(pca_data, !Is_Outlier), 
               aes(group = Cluster), level = 0.95, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_shape_identity() +
  theme_bw() +
  labs(
    title = paste("PPV OCTA Clustering Results - Manual Outlier Removal"),
    subtitle = paste("Method:", best_result$method, "| K =", best_result$k, 
                     "| Silhouette =", round(best_result$silhouette, 3)),
    x = paste0("PC1 (", round(summary(pca_result_full)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_full)$importance[2,2]*100, 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  annotate("text", x = max(pca_data$PC1), y = max(pca_data$PC2), 
           label = "X = Manual outliers", hjust = 1, vjust = 1, size = 3, color = "gray50")

# 保存图形
dir.create("plots", recursive = TRUE, showWarnings = FALSE)
ggsave("plots/ppv_octa_manual_outlier_removal_pca.pdf", p_pca_manual, width = 12, height = 8)
ggsave("plots/ppv_octa_manual_outlier_removal_pca.png", p_pca_manual, width = 12, height = 8, dpi = 300)
print(p_pca_manual)

# 创建清洁数据的PCA图（不包含异常值）
pca_result_clean <- prcomp(clean_data_scaled)

pca_data_clean <- data.frame(
  PC1 = pca_result_clean$x[,1],
  PC2 = pca_result_clean$x[,2],
  Cluster = factor(clean_clusters),
  ID = clean_ids
)

p_pca_clean <- ggplot(pca_data_clean, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = ID), vjust = -0.8, size = 3) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(
    title = "PPV OCTA Clustering Results - Clean Data Only",
    subtitle = paste("Method:", best_result$method, "| K =", best_result$k, 
                     "| Silhouette =", round(best_result$silhouette, 3)),
    x = paste0("PC1 (", round(summary(pca_result_clean)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result_clean)$importance[2,2]*100, 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

ggsave("plots/ppv_octa_clean_data_only_pca.pdf", p_pca_clean, width = 10, height = 8)
ggsave("plots/ppv_octa_clean_data_only_pca.png", p_pca_clean, width = 10, height = 8, dpi = 300)
print(p_pca_clean)

# -------------------- 保存结果 --------------------
cat("\n========== 保存结果 ==========\n")

# 保存手动排除异常值的聚类结果
write.csv(clustering_results_manual, "ppv_octa_manual_outlier_removal_results.csv", row.names = FALSE)

# 保存用于相关性分析的结果
correlation_results <- data.frame(
  subject_id = clustering_results_manual$subject_id,
  octa_cluster = clustering_results_manual$cluster,
  clustering_method = best_result$method,
  optimal_k = best_result$k,
  silhouette_score = best_result$silhouette
)

write.csv(correlation_results, "ppv_octa_manual_clusters_for_correlation.csv", row.names = FALSE)

# 保存聚类总结
cluster_summary_manual <- clustering_results_manual %>%
  group_by(cluster) %>%
  summarise(
    Count = n(),
    Patient_IDs = paste(sort(subject_id), collapse = ", "),
    Contains_Outliers = any(is_outlier),
    .groups = "drop"
  )

write.csv(cluster_summary_manual, "ppv_octa_manual_clustering_summary.csv", row.names = FALSE)

cat("聚类总结:\n")
print(cluster_summary_manual)

# -------------------- 结果总结 --------------------
cat("\n========== 手动排除异常值聚类分析完成 ==========\n")
cat("排除的异常值:", paste(manual_outliers, collapse = ", "), "\n")
cat("最佳聚类方法:", best_result$method, "\n")
cat("聚类数:", best_result$k, " (+ 1个异常值类)\n")
cat("清洁数据轮廓系数:", round(best_result$silhouette, 3), "\n")

cat("\n最终聚类分布:\n")
for(cluster_id in sort(unique(full_clusters))) {
  cluster_count <- sum(full_clusters == cluster_id)
  if(cluster_id == max_cluster + 1) {
    cat("聚类", cluster_id, "(异常值):", cluster_count, "人\n")
  } else {
    cat("聚类", cluster_id, ":", cluster_count, "人\n")
  }
}

cat("\n保存的文件:\n")
cat("- ppv_octa_manual_outlier_removal_results.csv: 完整聚类结果\n")
cat("- ppv_octa_manual_clusters_for_correlation.csv: 用于相关性分析\n")
cat("- ppv_octa_manual_clustering_summary.csv: 聚类总结\n")
cat("- plots/ppv_octa_manual_outlier_removal_pca.png: 包含异常值的PCA图\n")
cat("- plots/ppv_octa_clean_data_only_pca.png: 清洁数据PCA图\n")

# # 读取完整结果
# full_results <- read.csv("ppv_octa_manual_outlier_removal_results.csv")

# 筛选出非异常值的患者
clean_results_only <- clustering_results_manual %>%
  filter(is_outlier == FALSE) %>%
  dplyr::select(-is_outlier)  # 移除is_outlier列

# 保存只包含清洁数据的结果
write.csv(clean_results_only, "ppv_octa_clean_only_results.csv", row.names = FALSE)

# 检查结果
print(table(clean_results_only$cluster))
