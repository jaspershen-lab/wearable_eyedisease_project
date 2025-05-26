library(tidyverse)
library(lubridate)
library(Biobase)
library(Mfuzz)
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra) # For PCA visualization
library(r4projects)
setwd(get_project_wd())
rm(list = ls())

####read data
ppv_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)
cat_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_NoD_Surg0_8h_filtered.csv", check.names = FALSE)

str(ppv_data)
str(cat_data)

dir.create("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m")

# -------------------- 1. 准备用于聚类的多指标数据 --------------------
# 定义我们要使用的指标列表
metrics <- c("mean_rhr_1", "max_rhr_1", "min_rhr_1", "sd_rhr_1", "skew_rhr_1", "cv_rhr_1", "kurt_rhr_1",
             "mean_bo","cv_bo","steps_max", "steps_total", "deep_sleep", "total_sleep")

# 检查这些指标是否在数据集中
check_metrics <- function(data, metrics, day = "day_0") {
  for (metric in metrics) {
    col_name <- paste0(day, "_", metric)
    if (!col_name %in% colnames(data)) {
      cat("Warning: Column", col_name, "not found in dataset\n")
      metrics <- metrics[metrics != metric]
    }
  }
  return(metrics)
}

# 检查第一天的数据，确认哪些指标可用
available_metrics <- check_metrics(ppv_data, metrics)
cat("可用的指标:", paste(available_metrics, collapse = ", "), "\n")

# 我们将选择3个关键指标来平衡信息量和维度
# 1. mean_rhr_1: 静息心率的核心指标
# 2. steps_total: 活动量指标
# 3. total_sleep: 睡眠质量指标
selected_metrics <- c("cv_rhr_1", "steps_max")
# selected_metrics <- c("cv_rhr_1", "cv_bo","steps_max")
selected_metrics <- intersect(selected_metrics, available_metrics)

# 如果某些选定的指标不可用，我们退回到使用可用指标中的前3个
if (length(selected_metrics) < 2 && length(available_metrics) >= 2) {
  selected_metrics <- available_metrics[1:2]
} else if (length(selected_metrics) < 2) {
  selected_metrics <- available_metrics
}

cat("最终选择的指标:", paste(selected_metrics, collapse = ", "), "\n")

# -------------------- 2. 提取多指标数据 --------------------
# 函数：提取每个选定指标的术前-4到术后30天的数据
extract_periop_data <- function(data, metric) {
  cols <- c()
  for (day in c(-4:-1, 0:30)) {  # 修改为术后30天
    if (day < 0) {
      day_str <- paste0("day_", day, "_", metric)
    } else {
      day_str <- paste0("day_", day, "_", metric)
    }
    cols <- c(cols, day_str)
  }
  
  # 检查指定的列是否存在于数据中
  existing_cols <- cols[cols %in% colnames(data)]
  
  if (length(existing_cols) == 0) {
    cat(sprintf("警告: 找不到与指标 %s 相关的任何列\n", metric))
    # 返回只有subject_id的数据框
    return(data %>% dplyr::select(subject_id))
  }
  
  # 如果只有部分列存在，打印警告
  if (length(existing_cols) < length(cols)) {
    missing_cols <- setdiff(cols, existing_cols)
    cat(sprintf("警告: 缺少以下列: %s\n", paste(missing_cols, collapse = ", ")))
  }
  
  # 选择这些列和subject_id
  result <- data %>% 
    dplyr::select(subject_id, all_of(existing_cols))
  
  return(result)
}

# -------------------- 3. NA分析和处理函数 --------------------
analyze_na <- function(data, group_name) {
  cat("\n============", group_name, "NA值分析 ============\n")
  cat("样本量:", nrow(data), "\n\n")
  
  # 移除subject_id列和surgery_date_rhr列（如果存在）
  data_sans_id <- data %>% dplyr::select(-subject_id, -matches("surgery_date"))
  
  # 统计每个列的NA值数量和比例
  na_by_col <- data_sans_id %>%
    summarise(across(everything(), ~sum(is.na(.)), .names = "{.col}_na"),
              across(everything(), ~mean(is.na(.)) * 100, .names = "{.col}_pct"))
  
  # 输出每个列的NA值情况
  na_cols <- grep("_na$", names(na_by_col), value = TRUE)
  for (col in na_cols) {
    var_name <- gsub("_na$", "", col)
    cat(sprintf("%s: %d NA (%.2f%%)\n", 
                var_name, 
                na_by_col[[col]], 
                na_by_col[[paste0(var_name, "_pct")]]))
  }
  
  # 找出数值型列的名称
  numeric_col_names <- names(data)[sapply(data, is.numeric)]
  numeric_col_names <- setdiff(numeric_col_names, "subject_id")  # 排除subject_id
  
  # 统计每个患者的NA数量
  na_by_subject <- data %>%
    mutate(na_count = rowSums(is.na(dplyr::select(., all_of(numeric_col_names)))))
  
  # 统计不同NA数量的患者分布
  na_distribution <- na_by_subject %>%
    count(na_count) %>%
    mutate(percentage = n / sum(n) * 100)
  
  cat("\n患者的NA值分布:\n")
  print(na_distribution)
  
  # 计算总体NA值比例
  total_cells <- nrow(data) * length(numeric_col_names)
  total_na <- sum(is.na(data[, numeric_col_names, drop = FALSE]))
  total_na_pct <- total_na / total_cells * 100
  
  cat(sprintf("\n总体NA值比例: %.2f%%\n", total_na_pct))
  
  return(na_by_subject)
}

# -------------------- 4. 合并和预处理多指标数据 --------------------
# 对每个组和每个指标提取数据并分析NA值
process_multi_metric_data <- function(data, group_name, metrics) {
  # 检查并过滤掉NA比例过高的指标
  filtered_metrics <- c()
  na_threshold <- 60  # 修改为更高的阈值，适应您的数据特点
  
  for (metric in metrics) {
    # 查找相关列
    metric_cols <- grep(paste0("_", metric, "$"), colnames(data), value = TRUE)
    if (length(metric_cols) == 0) {
      cat(sprintf("警告: 找不到与 %s 相关的列\n", metric))
      next
    }
    
    # 计算该指标的NA比例
    na_percent <- mean(is.na(data[, metric_cols])) * 100
    cat(sprintf("指标 %s 的NA比例: %.2f%%\n", metric, na_percent))
    
    if (na_percent <= na_threshold) {
      filtered_metrics <- c(filtered_metrics, metric)
    } else {
      cat(sprintf("警告: 指标 %s 的NA比例超过阈值 %.0f%%, 将被排除\n", metric, na_threshold))
    }
  }
  
  if (length(filtered_metrics) == 0) {
    cat("错误: 所有指标的NA值比例都太高，无法进行分析\n")
    cat("将使用原始指标中NA比例最低的两个\n")
    
    # 计算每个指标的NA比例并选择最低的两个
    na_percents <- numeric(length(metrics))
    names(na_percents) <- metrics
    
    for (i in 1:length(metrics)) {
      metric_cols <- grep(paste0("_", metrics[i], "$"), colnames(data), value = TRUE)
      if (length(metric_cols) > 0) {
        na_percents[i] <- mean(is.na(data[, metric_cols])) * 100
      } else {
        na_percents[i] <- 100
      }
    }
    
    # 排序并选择NA比例最低的两个指标
    sorted_metrics <- names(sort(na_percents))
    filtered_metrics <- sorted_metrics[1:min(2, length(sorted_metrics))]
  }
  
  cat(sprintf("最终使用的指标: %s\n", paste(filtered_metrics, collapse = ", ")))
  
  # 先提取第一个指标以得到基础数据框
  result_df <- extract_periop_data(data, filtered_metrics[1])
  colnames(result_df)[-1] <- paste0(colnames(result_df)[-1], "_", filtered_metrics[1])
  
  # 为每个额外的指标添加相应的列
  if (length(filtered_metrics) > 1) {
    for (i in 2:length(filtered_metrics)) {
      metric_df <- extract_periop_data(data, filtered_metrics[i])
      colnames(metric_df)[-1] <- paste0(colnames(metric_df)[-1], "_", filtered_metrics[i])
      result_df <- result_df %>% 
        left_join(metric_df, by = "subject_id")
    }
  }
  
  # 添加手术日期列（如果可用）
  if ("surgery_date_rhr" %in% colnames(data)) {
    result_df <- result_df %>%
      left_join(data %>% dplyr::select(subject_id, surgery_date_rhr), by = "subject_id")
  }
  
  # 分析NA值
  na_analysis <- analyze_na(result_df, group_name)
  
  # 用均值填充NA值
  for (col in colnames(result_df)[-1]) {  # 跳过subject_id列
    if (col != "surgery_date_rhr" && is.numeric(result_df[[col]])) {  # 跳过日期列和非数值列
      # 检查该列是否有足够的非NA值进行填充
      if (sum(!is.na(result_df[[col]])) > 0) {
        result_df[is.na(result_df[[col]]), col] <- mean(result_df[[col]], na.rm = TRUE)
      } else {
        result_df[[col]] <- 0  # 如果全是NA，填充为0
        cat(sprintf("警告: 列 %s 全是NA，已填充为0\n", col))
      }
    }
  }
  
  return(result_df)
}

# 处理PPV组数据
ppv_multi <- process_multi_metric_data(ppv_data, "PPV组", selected_metrics)

# 处理白内障组数据
cat_multi <- process_multi_metric_data(cat_data, "白内障组", selected_metrics)

# -------------------- 5. 准备Mfuzz聚类的数据 --------------------
# 使用Z-score标准化每个指标
# 这对于有多种不同尺度的指标特别重要
standardize_data <- function(data) {
  # 标识数值列，排除subject_id和surgery_date
  numeric_cols <- sapply(data, is.numeric)
  numeric_cols["subject_id"] <- FALSE
  if ("surgery_date_rhr" %in% names(numeric_cols)) {
    numeric_cols["surgery_date_rhr"] <- FALSE
  }
  
  # 获取要标准化的列名
  cols_to_standardize <- names(numeric_cols)[numeric_cols]
  
  # 标准化
  data_std <- data
  for (col in cols_to_standardize) {
    # 检查列是否有足够的非NA值
    if (sum(!is.na(data[[col]])) > 1) {
      data_std[[col]] <- scale(data[[col]])
    } else {
      cat(sprintf("警告: 列 %s 没有足够的非NA值进行标准化\n", col))
    }
  }
  
  return(data_std)
}

# 标准化数据
ppv_std <- standardize_data(ppv_multi)
cat_std <- standardize_data(cat_multi)

# 准备Mfuzz聚类的数据
prep_data_for_mfuzz <- function(data) {
  # 提取数值列，排除subject_id和surgery_date
  data_matrix <- data %>%
    dplyr::select(-subject_id, -surgery_date_rhr) %>%
    as.matrix()
  
  # 设置行名
  rownames(data_matrix) <- data$subject_id
  
  # 创建ExpressionSet对象
  eset <- ExpressionSet(assayData = data_matrix)
  
  return(eset)
}

# 创建ExpressionSet对象
ppv_eset <- prep_data_for_mfuzz(ppv_std)
cat_eset <- prep_data_for_mfuzz(cat_std)

# Mfuzz需要进一步标准化
ppv_eset_std <- standardise(ppv_eset)
cat_eset_std <- standardise(cat_eset)

# -------------------- 6. 降维分析（可选） --------------------
# 对于小样本量高维数据，可以考虑先进行PCA降维
# 然后对主成分进行聚类，而不是直接对原始特征聚类

# 定义PCA函数
perform_pca <- function(data, n_components = 5) {
  # 移除非数值列
  data_numeric <- data %>% dplyr::select(-subject_id, -surgery_date_rhr)
  
  # 执行PCA
  pca_result <- prcomp(data_numeric, scale. = TRUE)
  
  # 解释方差比例
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_var <- cumsum(var_explained)
  
  # 打印解释方差
  cat("\nPCA解释方差:\n")
  for (i in 1:min(10, length(var_explained))) {
    cat(sprintf("PC%d: %.2f%% (累计: %.2f%%)\n", 
                i, var_explained[i]*100, cumulative_var[i]*100))
  }
  
  # 可视化PCA结果
  fviz_eig(pca_result, addlabels = TRUE, title = "Scree Plot")
  
  # 选择主成分数量，保证至少解释70%的方差
  n_components <- min(n_components, which(cumulative_var >= 0.7)[1])
  cat(sprintf("\n选择%d个主成分，解释了%.2f%%的总方差\n", 
              n_components, cumulative_var[n_components]*100))
  
  # 提取主成分得分
  pc_scores <- as.data.frame(pca_result$x[, 1:n_components])
  
  # 添加回subject_id
  pc_scores$subject_id <- data$subject_id
  
  return(list(
    pca_result = pca_result,
    pc_scores = pc_scores,
    n_components = n_components,
    var_explained = var_explained
  ))
}

# 执行PCA
cat("\n执行PPV组PCA降维...\n")
ppv_pca <- perform_pca(ppv_multi)
cat("\n执行白内障组PCA降维...\n")
cat_pca <- perform_pca(cat_multi)

# 使用PCA结果创建新的ExpressionSet对象
ppv_pca_eset <- ExpressionSet(assayData = as.matrix(ppv_pca$pc_scores %>% dplyr::select(-subject_id)))
rownames(ppv_pca_eset) <- ppv_multi$subject_id
ppv_pca_eset_std <- standardise(ppv_pca_eset)

cat_pca_eset <- ExpressionSet(assayData = as.matrix(cat_pca$pc_scores %>% dplyr::select(-subject_id)))
rownames(cat_pca_eset) <- cat_multi$subject_id
cat_pca_eset_std <- standardise(cat_pca_eset)

# # -------------------- 7. 执行Mfuzz聚类 --------------------
# 我们可以选择使用PCA降维后的数据进行聚类
# 或者直接使用原始多指标数据进行聚类

# 判断使用哪种数据集进行聚类
# - 如果样本量足够大，可以直接使用多指标数据
# - 如果样本量小但维度高，建议使用PCA降维后的数据
use_pca_for_clustering <- FALSE  # 可以根据需要更改此设置

# 确定模糊系数m
ppv_m <- mestimate(if(use_pca_for_clustering) ppv_pca_eset_std else ppv_eset_std)
cat_m <- mestimate(if(use_pca_for_clustering) cat_pca_eset_std else cat_eset_std)

cat("PPV组估计的最佳模糊系数 m:", ppv_m, "\n")
cat("白内障组估计的最佳模糊系数 m:", cat_m, "\n")

# 确定最佳聚类数的函数
optimal_cluster <- function(eset_std, max_c = 6, m = 1.25) {
  # 计算不同聚类数的Dmin (最小质心间距离)
  dmin_values <- numeric(max_c - 1)

  for (c in 2:max_c) {
    set.seed(123)
    cl <- mfuzz(eset_std, c = c, m = m)
    centers <- cl$centers

    # 计算所有质心之间的最小距离
    min_dist <- Inf
    for (i in 1:(c-1)) {
      for (j in (i+1):c) {
        dist_ij <- sqrt(sum((centers[i,] - centers[j,])^2))
        min_dist <- min(min_dist, dist_ij)
      }
    }

    dmin_values[c-1] <- min_dist
  }

  # 绘制Dmin图
  plot(2:max_c, dmin_values, type = "b",
       xlab = "聚类数", ylab = "最小质心间距离",
       main = "最小质心间距离vs聚类数")

  # 返回建议的聚类数 (Dmin急剧下降后的点)
  return(dmin_values)
}

# 运行以确定最佳聚类数
cat("计算PPV组的最佳聚类数...\n")
ppv_dmin_values <- optimal_cluster(
  if(use_pca_for_clustering) ppv_pca_eset_std else ppv_eset_std,
  max_c = 6,
  m = ppv_m
)

cat("计算白内障组的最佳聚类数...\n")
cat_dmin_values <- optimal_cluster(
  if(use_pca_for_clustering) cat_pca_eset_std else cat_eset_std,
  max_c = 6,
  m = cat_m
)

# 执行聚类
set.seed(123)
# 考虑到小样本量，我们限制最大聚类数为3
perform_clustering <- function(eset_std, m, group_name, max_clusters = 3) {
  # 初始尝试的聚类数
  n_clusters <- min(3, nrow(eset_std) / 5)  # 确保每个聚类至少有5个样本
  n_clusters <- max(2, round(n_clusters))    # 至少使用2个聚类

  cat(sprintf("\n尝试对%s使用%d个聚类...\n", group_name, n_clusters))
  set.seed(123)
  cl <- mfuzz(eset_std, c = n_clusters, m = m)

  # 获取每个样本的主要聚类分配
  main_clusters <- apply(cl$membership, 1, which.max)
  cluster_counts <- table(main_clusters)
  cat(sprintf("%s各聚类的成员数量:\n", group_name))
  print(cluster_counts)

  # 检查是否有空聚类或成员严重不平衡
  if(any(cluster_counts < 3) || length(cluster_counts) < n_clusters) {
    cat(sprintf("\n当尝试%d个聚类时有问题，改为使用2个聚类...\n", n_clusters))
    n_clusters <- 2
    set.seed(123)
    cl <- mfuzz(eset_std, c = n_clusters, m = m)

    # 再次检查聚类分配
    main_clusters <- apply(cl$membership, 1, which.max)
    cluster_counts <- table(main_clusters)
    cat(sprintf("%s使用2个聚类的成员数量:\n", group_name))
    print(cluster_counts)
  }

  return(list(
    cl = cl,
    n_clusters = n_clusters,
    main_clusters = main_clusters
  ))
}

# 执行PPV组聚类
ppv_clustering <- perform_clustering(
  if(use_pca_for_clustering) ppv_pca_eset_std else ppv_eset_std,
  ppv_m,
  "PPV"
)

# 执行白内障组聚类
cat_clustering <- perform_clustering(
  if(use_pca_for_clustering) cat_pca_eset_std else cat_eset_std,
  cat_m,
  "Cataract"
)

# 创建可保存的聚类结果数据框
# 正确创建包含所需列的数据框
ppv_clusters_to_save <- data.frame(
  subject_id = rownames(ppv_eset_std@assayData$exprs),
  max_cluster = ppv_clustering$main_clusters
)

# 添加最大成员度值
# 从成员度矩阵中提取每个样本的最大成员度值
ppv_max_membership <- apply(ppv_clustering$cl$membership, 1, max)
ppv_clusters_to_save$max_membership <- ppv_max_membership

# 对白内障组做同样的操作
cat_clusters_to_save <- data.frame(
  subject_id = rownames(cat_eset_std@assayData$exprs),
  max_cluster = cat_clustering$main_clusters
)
cat_max_membership <- apply(cat_clustering$cl$membership, 1, max)
cat_clusters_to_save$max_membership <- cat_max_membership

# 保存结果
write.csv(ppv_clusters_to_save, "ppv_cluster_results_all_metrics.csv", row.names = FALSE)
write.csv(cat_clusters_to_save, "cataract_cluster_results_all_metrics.csv", row.names = FALSE)


# # -------------------- 7. 执行Mfuzz聚类（不降维，直接聚类） --------------------
# 
# # 估算最佳模糊系数 m
# ppv_m <- mestimate(ppv_eset_std)
# cat_m <- mestimate(cat_eset_std)
# 
# cat("PPV组估计的最佳模糊系数 m:", ppv_m, "\n")
# cat("白内障组估计的最佳模糊系数 m:", cat_m, "\n")
# 
# # 画Dmin图，辅助选择聚类数
# optimal_cluster <- function(eset_std, max_c = 6, m = 1.25) {
#   dmin_values <- numeric(max_c - 1)
#   for (c in 2:max_c) {
#     cl <- mfuzz(eset_std, c = c, m = m)
#     centers <- cl$centers
#     min_dist <- Inf
#     for (i in 1:(c-1)) {
#       for (j in (i+1):c) {
#         dist_ij <- sqrt(sum((centers[i,] - centers[j,])^2))
#         min_dist <- min(min_dist, dist_ij)
#       }
#     }
#     dmin_values[c-1] <- min_dist
#   }
#   plot(2:max_c, dmin_values, type = "b", 
#        xlab = "聚类数", ylab = "最小质心间距离",
#        main = "最小质心间距离vs聚类数")
#   return(dmin_values)
# }
# 
# cat("计算PPV组的最佳聚类数...\n")
# ppv_dmin_values <- optimal_cluster(ppv_eset_std, max_c = 6, m = ppv_m)
# cat("计算白内障组的最佳聚类数...\n")
# cat_dmin_values <- optimal_cluster(cat_eset_std, max_c = 6, m = cat_m)
# 
# # 设定聚类数（如2或3）
# ppv_cluster_num <- 2
# cat_cluster_num <- 2
# 
# # 执行聚类
# ppv_clustering <- mfuzz(ppv_eset_std, c = ppv_cluster_num, m = ppv_m)
# cat_clustering <- mfuzz(cat_eset_std, c = cat_cluster_num, m = cat_m)
# 
# # 获取每个样本的主聚类分配
# ppv_main_clusters <- apply(ppv_clustering$membership, 1, which.max)
# cat_main_clusters <- apply(cat_clustering$membership, 1, which.max)
# 
# # 创建结果数据框
# ppv_clusters_to_save <- data.frame(
#   subject_id = rownames(ppv_eset_std@assayData$exprs),
#   max_cluster = ppv_main_clusters,
#   max_membership = apply(ppv_clustering$membership, 1, max)
# )
# 
# cat_clusters_to_save <- data.frame(
#   subject_id = rownames(cat_eset_std@assayData$exprs),
#   max_cluster = cat_main_clusters,
#   max_membership = apply(cat_clustering$membership, 1, max)
# )

# -------------------- 8. 可视化聚类结果 --------------------
# 由于我们可能使用了PCA降维数据进行聚类，原始的mfuzz.plot可能不适用
# 我们需要创建自定义可视化函数来展示聚类结果

# After performing the clustering, add this code to create the results dataframes:

# Create ppv_results with clustering information
ppv_results <- ppv_multi %>%
  mutate(
    max_cluster = ppv_clustering$main_clusters,
    cluster_membership = apply(ppv_clustering$cl$membership, 1, function(row) paste(round(row, 2), collapse=","))
  )

# Create cat_results with clustering information
cat_results <- cat_multi %>%
  mutate(
    max_cluster = cat_clustering$main_clusters,
    cluster_membership = apply(cat_clustering$cl$membership, 1, function(row) paste(round(row, 2), collapse=","))
  )

# -------------------- 8.3 Create improved cluster profile visualizations --------------------

# Function to visualize cluster profiles with membership coloring
visualize_cluster_profiles <- function(results_df, clustering_result, metrics, group_prefix) {
  # Create directory for plots if it doesn't exist
  dir.create("plots/cluster_profiles", recursive = TRUE, showWarnings = FALSE)
  
  # Get number of clusters
  n_clusters <- ncol(clustering_result$cl$membership)
  
  # For each cluster, create a plot showing all metrics
  for (cluster_id in 1:n_clusters) {
    # Filter data for this cluster
    cluster_data <- results_df %>% 
      filter(max_cluster == cluster_id)
    
    if (nrow(cluster_data) == 0) {
      cat(sprintf("Warning: No data found for %s Cluster %d\n", group_prefix, cluster_id))
      next
    }
    
    # Extract membership values for this cluster
    membership_df <- data.frame(
      subject_id = rownames(clustering_result$cl$membership),
      membership = clustering_result$cl$membership[, cluster_id],
      stringsAsFactors = FALSE
    )
    
    # Create plot data for all metrics
    plot_data <- data.frame()
    
    for (metric in metrics) {
      # Find metric columns
      metric_cols <- grep(paste0("_", metric, "$"), colnames(results_df), value = TRUE)
      
      if (length(metric_cols) == 0) next
      
      # Prepare data for this metric
      metric_data <- cluster_data %>%
        dplyr::select(subject_id, all_of(metric_cols)) %>%
        pivot_longer(
          cols = all_of(metric_cols),
          names_to = "day_metric",
          values_to = "value"
        ) %>%
        mutate(
          day = as.numeric(gsub("^day_(-?\\d+)_.*$", "\\1", day_metric)),
          metric = metric
        ) %>%
        # Join with membership values
        left_join(membership_df, by = "subject_id")
      
      plot_data <- bind_rows(plot_data, metric_data)
    }
    
    # Calculate mean profile for each metric and day
    mean_profile <- plot_data %>%
      group_by(metric, day) %>%
      summarize(value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # Create plot
    p <- ggplot() +
      # Individual lines colored by membership
      geom_line(data = plot_data, 
                aes(x = day, y = value, group = subject_id, color = membership),
                size = 0.8, alpha = 0.6) +
      # Mean profile as bold black line
      geom_line(data = mean_profile,
                aes(x = day, y = value, group = metric),
                color = "black",  
                size = 1.2) +
      # Color gradient for membership values
      scale_color_gradientn(
        colors = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
                   "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027"),
        limits = c(0.2, 1.0),
        oob = scales::squish,
        breaks = seq(0.2, 1.0, by = 0.1),
        name = "Membership"
      ) +
      # Facet by metric
      facet_wrap(~ metric, scales = "free_y", ncol = 1) +
      # Add a vertical line at day 0 (surgery day)
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      # Set x-axis breaks
      scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7)) +
      # Labels and theme
      labs(
        title = paste(group_prefix, "Cluster", cluster_id, 
                      "(n=", nrow(cluster_data), ")"),
        x = "Perioperative Days",
        y = "Value"
      ) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        legend.position = "top",
        legend.key.width = unit(3, "cm"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      )
    
    # Save the plot
    ggsave(paste0("plots/cluster_profiles/", group_prefix, "_cluster_", cluster_id, "_profile.pdf"),
           p, width = 10, height = 8)
    ggsave(paste0("plots/cluster_profiles/", group_prefix, "_cluster_", cluster_id, "_profile.png"),
           p, width = 10, height = 8, dpi = 300)
    
    # Print the plot
    print(p)
  }
}

# Additionally, create correlation heatmaps as shown in the second code
create_correlation_heatmap <- function(clustering_result, group_prefix) {
  # Create directory for correlation plots
  dir.create("plots/correlations", recursive = TRUE, showWarnings = FALSE)
  
  # Extract centers from clustering result
  centers <- clustering_result$cl$centers
  rownames(centers) <- paste0("Cluster", 1:nrow(centers))
  
  # Calculate correlation matrix
  cor_matrix <- cor(t(centers))
  
  # Check if corrplot package is available, if not, use heatmap
  if (requireNamespace("corrplot", quietly = TRUE) && 
      requireNamespace("RColorBrewer", quietly = TRUE)) {
    # Use corrplot for visualization
    pdf(paste0("plots/correlations/", group_prefix, "_cluster_correlations.pdf"), 
        width = 10, height = 10)
    corrplot::corrplot(
      cor_matrix,
      type = "full",
      diag = TRUE,
      order = "hclust",
      hclust.method = "ward.D",
      col = colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
      number.cex = 0.7,
      addCoef.col = "black"
    )
    dev.off()
    
    # Also create a PNG version
    png(paste0("plots/correlations/", group_prefix, "_cluster_correlations.png"), 
        width = 1000, height = 1000, res = 100)
    corrplot::corrplot(
      cor_matrix,
      type = "full",
      diag = TRUE,
      order = "hclust",
      hclust.method = "ward.D",
      col = colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
      number.cex = 0.7,
      addCoef.col = "black"
    )
    dev.off()
    
    # Also display in the current device
    corrplot::corrplot(
      cor_matrix,
      type = "full",
      diag = TRUE,
      order = "hclust",
      hclust.method = "ward.D",
      col = colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
      number.cex = 0.7,
      addCoef.col = "black"
    )
  } else {
    # Fallback to basic heatmap
    pdf(paste0("plots/correlations/", group_prefix, "_cluster_correlations.pdf"), 
        width = 10, height = 10)
    heatmap(cor_matrix, 
            main = paste(group_prefix, "Cluster Correlations"),
            symm = TRUE,
            cexRow = 0.9, 
            cexCol = 0.9)
    dev.off()
    
    # Also create a PNG version
    png(paste0("plots/correlations/", group_prefix, "_cluster_correlations.png"), 
        width = 1000, height = 1000, res = 100)
    heatmap(cor_matrix, 
            main = paste(group_prefix, "Cluster Correlations"),
            symm = TRUE,
            cexRow = 0.9, 
            cexCol = 0.9)
    dev.off()
    
    # Also display in the current device
    heatmap(cor_matrix, 
            main = paste(group_prefix, "Cluster Correlations"),
            symm = TRUE,
            cexRow = 0.9, 
            cexCol = 0.9)
  }
  
  cat(sprintf("%s cluster correlation heatmap created\n", group_prefix))
}

# Execute the new visualization functions
cat("\nCreating improved cluster profile visualizations for PPV group...\n")
visualize_cluster_profiles(ppv_results, ppv_clustering, selected_metrics, "PPV")

cat("\nCreating improved cluster profile visualizations for Cataract group...\n")
visualize_cluster_profiles(cat_results, cat_clustering, selected_metrics, "Cataract")

# Create correlation heatmaps
cat("\nCreating correlation heatmaps...\n")
create_correlation_heatmap(ppv_clustering, "PPV")
create_correlation_heatmap(cat_clustering, "Cataract")

# -------------------- 8.1 ggplot2自定义聚类趋势图（每个指标） --------------------
# 生成聚类分组信息后，批量绘制每个指标的聚类趋势图
# 适用于ppv_results和cat_results，指标列表为selected_metrics
ppv_results <- ppv_multi %>%
  mutate(
    max_cluster = ppv_clustering$main_clusters,
    cluster_membership = apply(ppv_clustering$cl$membership, 1, function(row) paste(round(row, 2), collapse=","))
  )

# Create cat_results with clustering information
cat_results <- cat_multi %>%
  mutate(
    max_cluster = cat_clustering$main_clusters,
    cluster_membership = apply(cat_clustering$cl$membership, 1, function(row) paste(round(row, 2), collapse=","))
  )


dir.create("plots/individual_metrics", recursive = TRUE, showWarnings = FALSE)

plot_cluster_trajectories <- function(results_df, metrics, group_prefix) {
  for (metric in metrics) {
    metric_cols <- grep(paste0("_", metric, "$"), colnames(results_df), value = TRUE)
    if (length(metric_cols) == 0) next
    plot_data <- results_df %>%
      dplyr::select(subject_id, max_cluster, all_of(metric_cols)) %>%
      pivot_longer(
        cols = all_of(metric_cols),
        names_to = "day_metric",
        values_to = "value"
      ) %>%
      mutate(
        day = as.numeric(gsub("^day_(-?\\d+)_.*$", "\\1", day_metric))
      )
    p <- ggplot(plot_data, aes(x = day, y = value, group = subject_id)) +
      geom_line(aes(color = factor(max_cluster)), alpha = 0.3) +
      stat_summary(aes(group = factor(max_cluster)), fun = mean, geom = "line", size = 1.5, color = "black") +
      facet_wrap(~ max_cluster) +
      labs(title = paste(group_prefix, "Clustered Trajectories for", metric),
           x = "Day", y = metric, color = "Cluster") +
      theme_bw() +
      theme(legend.position = "bottom")
    ggsave(paste0("plots/individual_metrics/", group_prefix, "_cluster_", metric, "_trend_custom.pdf"), p, width = 10, height = 6)
  }
}

# 对PPV组和白内障组分别绘图
plot_cluster_trajectories(ppv_results, selected_metrics, "ppv")
plot_cluster_trajectories(cat_results, selected_metrics, "cat")


# 为每个指标创建聚类趋势可视化函数
visualize_metric_clusters <- function(data, metric_name, group_name) {
  # 首先检查数据是否有聚类信息
  if (!("max_cluster" %in% colnames(data))) {
    cat("错误: 数据中没有聚类信息(max_cluster列)\n")
    return(NULL)
  }
  
  # 找出与该指标相关的列
  metric_cols <- grep(paste0("_", metric_name, "$"), colnames(data), value = TRUE)
  
  # 如果没找到相关列，尝试匹配不同的模式
  if (length(metric_cols) == 0) {
    # 尝试匹配日期模式的列
    day_pattern <- "^day_[-0-9]+_"
    all_day_cols <- grep(day_pattern, colnames(data), value = TRUE)
    metric_cols <- all_day_cols[grepl(metric_name, all_day_cols)]
  }
  
  # 如果仍然没找到，返回NULL
  if (length(metric_cols) == 0) {
    cat(sprintf("警告: 找不到与指标 %s 相关的列，无法创建可视化\n", metric_name))
    return(NULL)
  }
  
  # 准备绘图数据
  plot_data <- data %>%
    group_by(max_cluster) %>%
    summarise(across(all_of(metric_cols), mean, na.rm = TRUE)) %>%
    pivot_longer(
      cols = all_of(metric_cols),
      names_to = "day_metric",
      values_to = "value"
    )
  
  # 提取天数信息（处理可能的不同列名模式）
  plot_data <- plot_data %>%
    mutate(
      # 尝试从列名中提取日期
      day = as.numeric(gsub("^day_([^_]+).*$", "\\1", day_metric)),
      # 处理负数天数
      day = ifelse(grepl("^day_-", day_metric),
                   -as.numeric(gsub("^day_-([^_]+).*$", "\\1", day_metric)),
                   day)
    )
  
  # 获取聚类数
  n_clusters <- length(unique(data$max_cluster))
  
  # 绘制趋势图
  plot <- ggplot(plot_data, aes(x = day, y = value, color = factor(max_cluster), group = max_cluster)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7)) +
    labs(
      title = paste(group_name, "Perioperative", metric_name, "Cluster Trends"),
      x = "Perioperative Days",
      y = metric_name,
      color = "Cluster"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(plot)
}

# 获取数据中的指标列表
get_metrics <- function(data) {
  # 提取 day_X_metric 格式的列
  day_cols <- grep("^day_", colnames(data), value = TRUE)
  
  # 提取指标名称
  metrics <- unique(gsub("^day_[-0-9]+_(.+)$", "\\1", day_cols))
  
  return(metrics)
}

# 提取每个组的指标
ppv_metrics <- get_metrics(ppv_results)
cat_metrics <- get_metrics(cat_results)

# 创建并保存单指标可视化
for (metric in ppv_metrics) {
  cat(sprintf("创建PPV组 %s 指标的可视化...\n", metric))
  ppv_plot <- visualize_metric_clusters(
    ppv_results,
    metric,
    "PPV Group"
  )
  
  if (!is.null(ppv_plot)) {
    # 创建目录（如果不存在）
    dir.create("plots/individual_metrics", recursive = TRUE, showWarnings = FALSE)
    ggsave(paste0("plots/individual_metrics/ppv_cluster_", metric, "_trends.pdf"),
           ppv_plot, width = 10, height = 6)
    ggsave(paste0("plots/individual_metrics/ppv_cluster_", metric, "_trends.png"),
           ppv_plot, width = 10, height = 6, dpi = 300)
    print(ppv_plot)
  }
}

for (metric in cat_metrics) {
  cat(sprintf("创建白内障组 %s 指标的可视化...\n", metric))
  cat_plot <- visualize_metric_clusters(
    cat_results,
    metric,
    "Cataract Group"
  )
  
  if (!is.null(cat_plot)) {
    ggsave(paste0("plots/individual_metrics/cataract_cluster_", metric, "_trends.pdf"),
           cat_plot, width = 10, height = 6)
    ggsave(paste0("plots/individual_metrics/cataract_cluster_", metric, "_trends.png"),
           cat_plot, width = 10, height = 6, dpi = 300)
    print(cat_plot)
  }
}


# -------------------- 8.2 创建组合指标趋势图 --------------------
# Create ppv_results with clustering information
ppv_results <- ppv_multi %>%
  mutate(
    max_cluster = ppv_clustering$main_clusters,
    cluster_membership = apply(ppv_clustering$cl$membership, 1, function(row) paste(round(row, 2), collapse=","))
  )

# Create cat_results with clustering information
cat_results <- cat_multi %>%
  mutate(
    max_cluster = cat_clustering$main_clusters,
    cluster_membership = apply(cat_clustering$cl$membership, 1, function(row) paste(round(row, 2), collapse=","))
  )

# Get the metrics for each group
get_metrics <- function(data) {
  # Extract day_X_metric format columns
  day_cols <- grep("^day_", colnames(data), value = TRUE)

  # Extract metric names
  metrics <- unique(gsub("^day_[-0-9]+_(.+)$", "\\1", day_cols))

  return(metrics)
}

# Extract metrics for each group
ppv_metrics <- get_metrics(ppv_results)
cat_metrics <- get_metrics(cat_results)

# 创建一个综合图，在一张图上显示所有指标的聚类趋势
create_combined_trends <- function(data, metrics, group_name) {
  # 首先检查数据是否有聚类信息
  if (!("max_cluster" %in% colnames(data))) {
    cat("错误: 数据中没有聚类信息(max_cluster列)\n")
    return(NULL)
  }

  plot_list <- list()
  valid_metrics <- c()

  for (metric in metrics) {
    # 找出与该指标相关的列
    metric_cols <- grep(paste0("_", metric, "$"), colnames(data), value = TRUE)

    # 如果没找到相关列，尝试不同的匹配模式
    if (length(metric_cols) == 0) {
      # 尝试匹配日期模式的列
      day_pattern <- "^day_[-0-9]+_"
      all_day_cols <- grep(day_pattern, colnames(data), value = TRUE)
      metric_cols <- all_day_cols[grepl(metric, all_day_cols)]
    }

    # 如果仍然没找到，跳过这个指标
    if (length(metric_cols) == 0) {
      cat(sprintf("警告: 找不到与指标 %s 相关的列，将在组合图中跳过\n", metric))
      next
    }

    # 添加到有效指标列表
    valid_metrics <- c(valid_metrics, metric)

    # 准备绘图数据
    plot_data <- data %>%
      group_by(max_cluster) %>%
      summarise(across(all_of(metric_cols), mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = all_of(metric_cols),
        names_to = "day_metric",
        values_to = "value"
      ) %>%
      # 提取天数信息
      mutate(
        day = as.numeric(gsub("^day_([^_]+).*$", "\\1", day_metric)),
        # 处理负数天数
        day = ifelse(grepl("^day_-", day_metric),
                     -as.numeric(gsub("^day_-([^_]+).*$", "\\1", day_metric)),
                     day),
        metric = metric
      )

    plot_list[[metric]] <- plot_data
  }

  # 如果没有有效的指标数据，返回NULL
  if (length(plot_list) == 0) {
    cat("警告: 没有找到任何可用于绘图的指标数据\n")
    return(NULL)
  }

  # 合并所有指标数据
  combined_data <- bind_rows(plot_list)

  # 标准化每个指标的值，以便在同一图上比较
  combined_data <- combined_data %>%
    group_by(metric) %>%
    mutate(
      # 标准化值为z-score
      value_std = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
    ) %>%
    ungroup()

  # 获取聚类数
  n_clusters <- length(unique(data$max_cluster))

  # 绘制趋势图
  plot <- ggplot(combined_data, aes(x = day, y = value_std, color = factor(max_cluster))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 5, 10, 15, 20, 25, 30)) +  # 修改x轴刻度
    labs(
      title = paste(group_name, "Perioperative Metrics by Cluster"),
      x = "Perioperative Days",
      y = "Standardized Value",
      color = "Cluster"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )

  # 创建非标准化的版本
  plot_raw <- ggplot(combined_data, aes(x = day, y = value, color = factor(max_cluster))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 5, 10, 15, 20, 25, 30)) +  # 修改x轴刻度
    labs(
      title = paste(group_name, "Perioperative Metrics by Cluster (Raw Values)"),
      x = "Perioperative Days",
      y = "Value",
      color = "Cluster"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )

  return(list(
    standardized = plot,
    raw = plot_raw,
    metrics = valid_metrics
  ))
}

# 创建PPV组组合趋势图
cat("\n创建PPV组组合趋势图...\n")
ppv_combined_plots <- create_combined_trends(
  ppv_results,
  ppv_metrics,
  "PPV Group"
)

if (!is.null(ppv_combined_plots)) {
  # 创建目录（如果不存在）
  dir.create("plots/combined", recursive = TRUE, showWarnings = FALSE)

  # 保存标准化和原始值图
  ggsave("plots/combined/ppv_combined_metrics_trends_standardized.pdf",
         ppv_combined_plots$standardized, width = 12, height = 8)
  ggsave("plots/combined/ppv_combined_metrics_trends_standardized.png",
         ppv_combined_plots$standardized, width = 12, height = 8, dpi = 300)

  ggsave("plots/combined/ppv_combined_metrics_trends_raw.pdf",
         ppv_combined_plots$raw, width = 12, height = 8)
  ggsave("plots/combined/ppv_combined_metrics_trends_raw.png",
         ppv_combined_plots$raw, width = 12, height = 8, dpi = 300)

  # 打印图表
  print(ppv_combined_plots$standardized)
  print(ppv_combined_plots$raw)
}

# 创建白内障组组合趋势图
cat("\n创建白内障组组合趋势图...\n")
cat_combined_plots <- create_combined_trends(
  cat_results,
  cat_metrics,
  "Cataract Group"
)

if (!is.null(cat_combined_plots)) {
  ggsave("plots/combined/cataract_combined_metrics_trends_standardized.pdf",
         cat_combined_plots$standardized, width = 12, height = 8)
  ggsave("plots/combined/cataract_combined_metrics_trends_standardized.png",
         cat_combined_plots$standardized, width = 12, height = 8, dpi = 300)

  ggsave("plots/combined/cataract_combined_metrics_trends_raw.pdf",
         cat_combined_plots$raw, width = 12, height = 8)
  ggsave("plots/combined/cataract_combined_metrics_trends_raw.png",
         cat_combined_plots$raw, width = 12, height = 8, dpi = 300)

  print(cat_combined_plots$standardized)
  print(cat_combined_plots$raw)
}

# 创建PCA散点图
# 假设ppv_pca和ppv_clustering已存在
pca_plot_data <- data.frame(
  PC1 = ppv_pca$pca_result$x[,1],
  PC2 = ppv_pca$pca_result$x[,2],
  Cluster = factor(ppv_clustering$main_clusters)
)

# 绘制PCA散点图
pca_cluster_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.5) +
  stat_ellipse(level = 0.95) +  # 添加置信椭圆
  labs(title = "PPV Group Patient Clusters",
       x = paste0("PC1 (", round(ppv_pca$var_explained[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(ppv_pca$var_explained[2]*100, 1), "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(pca_cluster_plot)
dir.create("plots/pca", recursive = TRUE, showWarnings = FALSE)
ggsave("plots/pca/ppv_pca_clusters.pdf", pca_cluster_plot, width = 10, height = 8)
ggsave("plots/pca/ppv_pca_clusters.png", pca_cluster_plot, width = 10, height = 8, dpi = 300)

cat_pca_plot_data <- data.frame(
  PC1 = cat_pca$pca_result$x[,1],
  PC2 = cat_pca$pca_result$x[,2],
  Cluster = factor(cat_clustering$main_clusters)
)

# 绘制PCA散点图
cat_pca_cluster_plot <- ggplot(cat_pca_plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.5) +
  stat_ellipse(level = 0.95) +  # 添加置信椭圆
  labs(title = "Cataract Group Patient Clusters",
       x = paste0("PC1 (", round(cat_pca$var_explained[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(cat_pca$var_explained[2]*100, 1), "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(cat_pca_cluster_plot)
ggsave("plots/pca/cataract_pca_clusters.pdf", cat_pca_cluster_plot, width = 10, height = 8)
ggsave("plots/pca/cataract_pca_clusters.png", cat_pca_cluster_plot, width = 10, height = 8, dpi = 300)
