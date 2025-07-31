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

dir.create("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/mfuzz/multi_metrics/1m/less_timepoint")

# -------------------- 1. 准备用于聚类的多指标数据 --------------------
# 定义我们要使用的指标列表
metrics <- c("mean_rhr_1", "max_rhr_1", "min_rhr_1", "sd_rhr_1", "skew_rhr_1", "cv_rhr_1", "kurt_rhr_1",
             "mean_bo","cv_bo","steps_max", "steps_total", "steps_mean","deep_sleep", "total_sleep")

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

# 我们将选择4个关键指标来平衡信息量和维度
selected_metrics <- c("cv_rhr_1","steps_max")
# selected_metrics <- c("cv_rhr_1","cv_bo", "steps_max")
# selected_metrics <- c("cv_rhr_1","cv_bo", "steps_max","deep_sleep")
selected_metrics <- intersect(selected_metrics, available_metrics)

# 如果某些选定的指标不可用，我们退回到使用可用指标中的前几个
if (length(selected_metrics) < 2 && length(available_metrics) >= 2) {
  selected_metrics <- available_metrics[1:2]
} else if (length(selected_metrics) < 2) {
  selected_metrics <- available_metrics
}

cat("最终选择的指标:", paste(selected_metrics, collapse = ", "), "\n")

# -------------------- 2. 基于时间窗口均值的数据提取（新方法） --------------------
# 定义时间窗口方案
time_windows_weekly <- list(
  baseline = list(days = -4:-1, name = "baseline"),               # Pre-surgical baseline
  acute_recovery = list(days = 0:3, name = "acute_recovery"),     # Immediate post-op recovery
  early_recovery = list(days = 4:7, name = "early_recovery"),     # First week recovery
  mid_recovery = list(days = 8:15, name = "mid_recovery"),        # Second week
  late_recovery = list(days = 16:30, name = "late_recovery")      # Later recovery
)

# 基于时间窗口提取数据的函数
extract_periop_data_windows <- function(data, metric, time_windows = time_windows_weekly) {
  result_data <- data %>% dplyr::select(subject_id)
  
  cat(sprintf("\n处理指标: %s\n", metric))
  
  for (window_name in names(time_windows)) {
    window <- time_windows[[window_name]]
    
    # 收集该时间窗口内所有可用的列
    window_cols <- c()
    for (day in window$days) {
      if (day < 0) {
        day_str <- paste0("day_", day, "_", metric)
      } else {
        day_str <- paste0("day_", day, "_", metric)
      }
      if (day_str %in% colnames(data)) {
        window_cols <- c(window_cols, day_str)
      }
    }
    
    if (length(window_cols) > 0) {
      # 计算该时间窗口的均值
      window_data <- data %>%
        dplyr::select(subject_id, all_of(window_cols))
      
      # 计算均值，但只有当至少有一半的数据不是NA时才计算
      min_valid_points <- max(1, floor(length(window_cols) / 2))
      
      window_mean <- window_data %>%
        mutate(
          valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
          window_mean = ifelse(
            valid_count >= min_valid_points,
            rowMeans(dplyr::select(., -subject_id), na.rm = TRUE),
            NA
          )
        ) %>%
        dplyr::select(subject_id, window_mean)
      
      # 重命名列
      colnames(window_mean)[2] <- paste0(window$name, "_", metric)
      
      # 合并到结果数据框
      result_data <- result_data %>%
        left_join(window_mean, by = "subject_id")
      
      cat(sprintf("  时间窗口 %s (days %d to %d): 使用了 %d 个时间点，要求至少 %d 个有效值\n", 
                  window$name, min(window$days), max(window$days), 
                  length(window_cols), min_valid_points))
    } else {
      cat(sprintf("  警告: 时间窗口 %s 中没有找到相关数据\n", window$name))
      
      # 如果没有找到数据，添加一个NA列以保持数据结构一致
      window_mean <- data.frame(
        subject_id = data$subject_id,
        temp_col = NA
      )
      colnames(window_mean)[2] <- paste0(window$name, "_", metric)
      result_data <- result_data %>%
        left_join(window_mean, by = "subject_id")
    }
  }
  
  return(result_data)
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

# -------------------- 4. 合并和预处理多指标数据（使用时间窗口） --------------------
process_multi_metric_data_windows <- function(data, group_name, metrics, time_windows = time_windows_weekly) {
  # 检查并过滤掉NA比例过高的指标
  filtered_metrics <- c()
  na_threshold <- 60  # 修改为更高的阈值，适应您的数据特点
  
  cat(sprintf("\n处理 %s 数据...\n", group_name))
  
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
  
  # 使用时间窗口方法提取第一个指标以得到基础数据框
  result_df <- extract_periop_data_windows(data, filtered_metrics[1], time_windows)
  
  # 为每个额外的指标添加相应的列
  if (length(filtered_metrics) > 1) {
    for (i in 2:length(filtered_metrics)) {
      metric_df <- extract_periop_data_windows(data, filtered_metrics[i], time_windows)
      # 移除subject_id列，然后合并
      metric_df <- metric_df %>% dplyr::select(-subject_id)
      result_df <- cbind(result_df, metric_df)
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

# 处理PPV组数据（使用时间窗口）
ppv_multi <- process_multi_metric_data_windows(ppv_data, "PPV组", selected_metrics, time_windows_weekly)

# 处理白内障组数据（使用时间窗口）
cat_multi <- process_multi_metric_data_windows(cat_data, "白内障组", selected_metrics, time_windows_weekly)

# -------------------- 5. 准备Mfuzz聚类的数据 --------------------
# 使用Z-score标准化每个指标
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
  cols_to_exclude <- c("subject_id")
  if ("surgery_date_rhr" %in% colnames(data)) {
    cols_to_exclude <- c(cols_to_exclude, "surgery_date_rhr")
  }
  
  data_matrix <- data %>%
    dplyr::select(-all_of(cols_to_exclude)) %>%
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
perform_pca <- function(data, n_components = 5) {
  # 移除非数值列
  cols_to_exclude <- c("subject_id")
  if ("surgery_date_rhr" %in% colnames(data)) {
    cols_to_exclude <- c(cols_to_exclude, "surgery_date_rhr")
  }
  
  data_numeric <- data %>% dplyr::select(-all_of(cols_to_exclude))
  
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
  
  # 由于使用了时间窗口，我们可以使用更少的主成分
  n_components <- min(min(3, n_components), which(cumulative_var >= 0.8)[1])
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

# -------------------- 7. 执行Mfuzz聚类 --------------------
# 确定模糊系数m
ppv_m <- mestimate(ppv_eset_std)
cat_m <- mestimate(cat_eset_std)

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
  
  return(dmin_values)
}

# 运行以确定最佳聚类数
cat("计算PPV组的最佳聚类数...\n")
ppv_dmin_values <- optimal_cluster(ppv_eset_std, max_c = 6, m = ppv_m)

cat("计算白内障组的最佳聚类数...\n")
cat_dmin_values <- optimal_cluster(cat_eset_std, max_c = 6, m = cat_m)

# 执行聚类
perform_clustering <- function(eset_std, m, group_name, max_clusters = 2) {
  # 由于使用了时间窗口，数据维度降低，可以尝试更多聚类
  n_clusters <- min(2, nrow(eset_std) / 5)  # 确保每个聚类至少有4个样本
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
ppv_clustering <- perform_clustering(ppv_eset_std, ppv_m, "PPV")

# 执行白内障组聚类
cat_clustering <- perform_clustering(cat_eset_std, cat_m, "Cataract")

# 创建可保存的聚类结果数据框
ppv_clusters_to_save <- data.frame(
  subject_id = rownames(ppv_eset_std@assayData$exprs),
  max_cluster = ppv_clustering$main_clusters
)

# 添加最大成员度值
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
write.csv(ppv_clusters_to_save, "ppv_cluster_results_time_windows.csv", row.names = FALSE)
write.csv(cat_clusters_to_save, "cataract_cluster_results_time_windows.csv", row.names = FALSE)

# -------------------- 8. 可视化聚类结果 --------------------
# 创建包含聚类信息的结果数据框
ppv_results <- ppv_multi %>%
  mutate(
    max_cluster = ppv_clustering$main_clusters,
    cluster_membership = apply(ppv_clustering$cl$membership, 1, function(row) paste(round(row, 2), collapse=","))
  )

cat_results <- cat_multi %>%
  mutate(
    max_cluster = cat_clustering$main_clusters,
    cluster_membership = apply(cat_clustering$cl$membership, 1, function(row) paste(round(row, 2), collapse=","))
  )

# 时间窗口聚类趋势可视化函数
visualize_time_window_clusters <- function(results_df, clustering_result, metrics, group_prefix, time_windows) {
  # Create plots directory
  dir.create("plots/cluster_profiles_time_windows", recursive = TRUE, showWarnings = FALSE)
  
  # Get number of clusters
  n_clusters <- ncol(clustering_result$cl$membership)
  
  # Create plots for each cluster
  for (cluster_id in 1:n_clusters) {
    # Filter data for this cluster
    cluster_data <- results_df %>% 
      filter(max_cluster == cluster_id)
    
    if (nrow(cluster_data) == 0) {
      cat(sprintf("Warning: No data found for %s Cluster %d\n", group_prefix, cluster_id))
      next
    }
    
    # Extract membership values
    membership_df <- data.frame(
      subject_id = rownames(clustering_result$cl$membership),
      membership = clustering_result$cl$membership[, cluster_id],
      stringsAsFactors = FALSE
    )
    
    # Create plot data for all metrics
    plot_data <- data.frame()
    
    for (metric in metrics) {
      # Find metric-related columns (time window version)
      metric_cols <- grep(paste0("_", metric, "$"), colnames(results_df), value = TRUE)
      
      if (length(metric_cols) == 0) next
      
      # Prepare data for this metric
      metric_data <- cluster_data %>%
        dplyr::select(subject_id, all_of(metric_cols)) %>%
        pivot_longer(
          cols = all_of(metric_cols),
          names_to = "window_metric",
          values_to = "value"
        ) %>%
        mutate(
          window = gsub(paste0("_", metric, "$"), "", window_metric),
          metric = metric,
          # Assign numerical values to time windows for plotting
          window_order = case_when(
            window == "baseline" ~ 1,
            window == "acute_recovery" ~ 2,
            window == "early_recovery" ~ 3,
            window == "mid_recovery" ~ 4,
            window == "late_recovery" ~ 5,
            TRUE ~ as.numeric(factor(window))
          )
        ) %>%
        # Join with membership values
        left_join(membership_df, by = "subject_id")
      
      plot_data <- bind_rows(plot_data, metric_data)
    }
    
    # Calculate mean profile for each metric and time window
    mean_profile <- plot_data %>%
      group_by(metric, window, window_order) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # Create the plot
    p <- ggplot() +
      # Individual lines colored by membership
      geom_line(data = plot_data, 
                aes(x = window_order, y = value, group = subject_id, color = membership),
                size = 0.8, alpha = 0.6) +
      # Mean trend line (thick black line)
      geom_line(data = mean_profile,
                aes(x = window_order, y = value, group = metric),
                color = "black",  
                size = 1.2) +
      # Membership color gradient
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
      # Surgery period vertical line
      geom_vline(xintercept = 2, linetype = "dashed", color = "gray40") +
      # Set x-axis ticks and labels
      scale_x_continuous(
        breaks = 1:5,
        labels = c("Baseline", "Acute Recovery", "Early Recovery", "Mid Recovery", "Late Recovery")
      ) +
      # Labels and theme
      labs(
        title = paste(group_prefix, "Cluster", cluster_id, 
                      "(n=", nrow(cluster_data), ") - Time Windows"),
        x = "Time Windows",
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
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Save plots
    ggsave(paste0("plots/cluster_profiles_time_windows/", group_prefix, "_cluster_", cluster_id, "_windows_profile.pdf"),
           p, width = 10, height = 8)
    ggsave(paste0("plots/cluster_profiles_time_windows/", group_prefix, "_cluster_", cluster_id, "_windows_profile.png"),
           p, width = 10, height = 8, dpi = 300)
    
    # Print the plot
    print(p)
  }
}

# 创建组合指标趋势图函数（包含自定义颜色）
create_combined_trends <- function(data, metrics, group_name, time_windows = time_windows_weekly) {
  # Check if data has clustering information
  if (!("max_cluster" %in% colnames(data))) {
    cat("Error: No clustering information (max_cluster column) found in data\n")
    return(NULL)
  }
  
  plot_list <- list()
  valid_metrics <- c()
  
  for (metric in metrics) {
    # Find columns related to this metric
    metric_cols <- grep(paste0("_", metric, "$"), colnames(data), value = TRUE)
    
    # Skip if no relevant columns found
    if (length(metric_cols) == 0) {
      cat(sprintf("Warning: No columns found for metric %s, skipping in combined plot\n", metric))
      next
    }
    
    # Add to valid metrics list
    valid_metrics <- c(valid_metrics, metric)
    
    # Prepare plot data
    metric_data <- data %>%
      group_by(max_cluster) %>%
      summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = "drop") %>%
      pivot_longer(
        cols = all_of(metric_cols),
        names_to = "window_metric",
        values_to = "value"
      ) %>%
      # Extract time window information
      mutate(
        window = gsub(paste0("_", metric, "$"), "", window_metric),
        metric = metric,
        # Assign numerical values to time windows for plotting
        window_order = case_when(
          window == "baseline" ~ 1,
          window == "acute_recovery" ~ 2,
          window == "early_recovery" ~ 3,
          window == "mid_recovery" ~ 4,
          window == "late_recovery" ~ 5,
          TRUE ~ as.numeric(factor(window))
        )
      )
    
    plot_list[[metric]] <- metric_data
  }
  
  # Return NULL if no valid metric data found
  if (length(plot_list) == 0) {
    cat("Warning: No usable metric data found for plotting\n")
    return(NULL)
  }
  
  # Combine all metric data
  combined_data <- bind_rows(plot_list)
  
  # Standardize values for each metric to allow comparison on same plot
  combined_data <- combined_data %>%
    group_by(metric) %>%
    mutate(
      # Standardize values to z-scores
      value_std = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Get number of clusters
  n_clusters <- length(unique(data$max_cluster))
  
  # Create standardized trend plot with custom colors
  plot_std <- ggplot(combined_data, aes(x = window_order, y = value_std, color = factor(max_cluster), group = interaction(max_cluster, metric))) +
    geom_line(size = 1.2) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "gray40") +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("1" = "#df8859", "2" = "#0fb292")) +
    scale_x_continuous(
      breaks = 1:5,
      labels = c("Baseline", "Acute Recovery", "Early Recovery", "Mid Recovery", "Late Recovery")
    ) +
    labs(
      title = paste(group_name, "Perioperative Metrics by Cluster (Standardized)"),
      x = "Time Window",
      y = "Standardized Value (Z-score)",
      color = "Cluster"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Create raw values trend plot with custom colors
  plot_raw <- ggplot(combined_data, aes(x = window_order, y = value, color = factor(max_cluster), group = interaction(max_cluster, metric))) +
    geom_line(size = 1.2) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "gray40") +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("1" = "#df8859", "2" = "#0fb292")) +
    scale_x_continuous(
      breaks = 1:5,
      labels = c("Baseline", "Acute Recovery", "Early Recovery", "Mid Recovery", "Late Recovery")
    ) +
    labs(
      title = paste(group_name, "Perioperative Metrics by Cluster (Raw Values)"),
      x = "Time Window",
      y = "Value",
      color = "Cluster"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(list(
    standardized = plot_std,
    raw = plot_raw,
    metrics = valid_metrics
  ))
}

# 创建聚类趋势可视化
cat("\nCreating time window cluster profile visualizations for PPV group...\n")
visualize_time_window_clusters(ppv_results, ppv_clustering, selected_metrics, "PPV", time_windows_weekly)

cat("\nCreating time window cluster profile visualizations for Cataract group...\n")
visualize_time_window_clusters(cat_results, cat_clustering, selected_metrics, "Cataract", time_windows_weekly)

# 创建PCA散点图
if (exists("ppv_pca") && exists("cat_pca")) {
  # PPV组PCA图
  pca_plot_data <- data.frame(
    PC1 = ppv_pca$pca_result$x[,1],
    PC2 = ppv_pca$pca_result$x[,2],
    Cluster = factor(ppv_clustering$main_clusters)
  )
  
  pca_cluster_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95) +
    labs(title = "PPV Group Patient Clusters (Time Windows)",
         x = paste0("PC1 (", round(ppv_pca$var_explained[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(ppv_pca$var_explained[2]*100, 1), "%)")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  print(pca_cluster_plot)
  dir.create("plots/pca_time_windows", recursive = TRUE, showWarnings = FALSE)
  ggsave("plots/pca_time_windows/ppv_pca_clusters_windows.pdf", pca_cluster_plot, width = 10, height = 8)
  ggsave("plots/pca_time_windows/ppv_pca_clusters_windows.png", pca_cluster_plot, width = 10, height = 8, dpi = 300)
  
  # 白内障组PCA图
  cat_pca_plot_data <- data.frame(
    PC1 = cat_pca$pca_result$x[,1],
    PC2 = cat_pca$pca_result$x[,2],
    Cluster = factor(cat_clustering$main_clusters)
  )
  
  cat_pca_cluster_plot <- ggplot(cat_pca_plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95) +
    labs(title = "Cataract Group Patient Clusters (Time Windows)",
         x = paste0("PC1 (", round(cat_pca$var_explained[1]*100, 1), "%)"),
         y = paste0("PC2 (", round(cat_pca$var_explained[2]*100, 1), "%)")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  print(cat_pca_cluster_plot)
  ggsave("plots/pca_time_windows/cataract_pca_clusters_windows.pdf", cat_pca_cluster_plot, width = 10, height = 8)
  ggsave("plots/pca_time_windows/cataract_pca_clusters_windows.png", cat_pca_cluster_plot, width = 10, height = 8, dpi = 300)
}

# -------------------- 9. 创建组合指标趋势图 --------------------
# 创建输出目录
dir.create("plots/combined_time_windows", recursive = TRUE, showWarnings = FALSE)

# 创建PPV组组合趋势图
cat("\n创建PPV组组合趋势图...\n")
ppv_combined_plots <- create_combined_trends(
  ppv_results, 
  selected_metrics,
  "PPV组"
)

if (!is.null(ppv_combined_plots)) {
  # 保存标准化和原始值图
  ggsave("plots/combined_time_windows/ppv_combined_metrics_std.pdf",
         ppv_combined_plots$standardized, width = 12, height = 10)
  ggsave("plots/combined_time_windows/ppv_combined_metrics_std.png",
         ppv_combined_plots$standardized, width = 12, height = 10, dpi = 300)
  
  ggsave("plots/combined_time_windows/ppv_combined_metrics_raw.pdf",
         ppv_combined_plots$raw, width = 12, height = 10)
  ggsave("plots/combined_time_windows/ppv_combined_metrics_raw.png",
         ppv_combined_plots$raw, width = 12, height = 10, dpi = 300)
  
  # 打印图表
  print(ppv_combined_plots$standardized)
  print(ppv_combined_plots$raw)
  
  cat("已创建PPV组合图，包含指标:", paste(ppv_combined_plots$metrics, collapse=", "), "\n")
}

# 创建白内障组组合趋势图
cat("\n创建白内障组组合趋势图...\n")
cat_combined_plots <- create_combined_trends(
  cat_results,
  selected_metrics,
  "白内障组"
)

if (!is.null(cat_combined_plots)) {
  ggsave("plots/combined_time_windows/cataract_combined_metrics_std.pdf",
         cat_combined_plots$standardized, width = 12, height = 10)
  ggsave("plots/combined_time_windows/cataract_combined_metrics_std.png",
         cat_combined_plots$standardized, width = 12, height = 10, dpi = 300)
  
  ggsave("plots/combined_time_windows/cataract_combined_metrics_raw.pdf",
         cat_combined_plots$raw, width = 12, height = 10)
  ggsave("plots/combined_time_windows/cataract_combined_metrics_raw.png",
         cat_combined_plots$raw, width = 12, height = 10, dpi = 300)
  
  print(cat_combined_plots$standardized)
  print(cat_combined_plots$raw)
  
  cat("已创建白内障组合图，包含指标:", paste(cat_combined_plots$metrics, collapse=", "), "\n")
}

cat("\n组合指标趋势图创建完成，保存在plots/combined_time_windows目录下\n")

cat("\n========== 基于时间窗口的聚类分析完成 ==========\n")
cat("使用的时间窗口：\n")
cat("1. 基线期 (day -4 to -1)\n")
cat("2. 急性恢复期 (day 0 to 3)\n")
cat("3. 早期恢复期 (day 4 to 7)\n")
cat("4. 中期恢复期 (day 8 to 15)\n")
cat("5. 晚期恢复期 (day 16 to 30)\n")
cat("聚类结果已保存到：\n")
cat("- ppv_cluster_results_time_windows.csv\n")
cat("- cataract_cluster_results_time_windows.csv\n")
cat("可视化图形保存在plots目录中\n")
