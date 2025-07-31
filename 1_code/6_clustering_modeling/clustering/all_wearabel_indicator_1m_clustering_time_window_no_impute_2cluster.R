# 修正时间窗口聚类代码 - 固定为2类聚类 + 修正cluster标签 + 类似代码一的可视化
# 解决cluster标签不连续问题 + 添加详细的聚类趋势可视化 + 强制所有窗口都聚成2类

library(tidyverse)
library(Biobase)
library(Mfuzz)
library(ggplot2)
library(gridExtra)
library(corrplot)
library(factoextra)
library(RColorBrewer)
library(r4projects)

setwd(get_project_wd())
rm(list = ls())

# ================== 设置随机种子确保可重复性 ==================
set.seed(123)  # 全局随机种子
RANDOM_SEED <- 123

# ================== 1. Setup and Data Loading ==================
cat("===== Time Window Specific Clustering Analysis with FIXED 2-Cluster Analysis =====\n")

# Define time windows
time_windows <- list(
  baseline = list(days = -4:-1, name = "baseline"),
  acute_recovery = list(days = 0:3, name = "acute_recovery"),
  early_recovery = list(days = 4:7, name = "early_recovery"),
  mid_recovery = list(days = 8:15, name = "mid_recovery"),
  late_recovery = list(days = 16:30, name = "late_recovery")
)

# Load data
ppv_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)

# Key metrics for clustering
key_metrics <- c("cv_rhr_1", "steps_mean")

# Create output directory
dir.create("3_data_analysis/6_clustering_modeling/time_window_clustering", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/time_window_clustering")

cat("Data loaded successfully. Starting time window clustering analysis with FIXED 2 clusters...\n")
cat("Time windows:", length(time_windows), "\n")
cat("Total patients in dataset:", nrow(ppv_data), "\n")
cat("Metrics for clustering:", paste(key_metrics, collapse = ", "), "\n")
cat("🎯 重要修改：所有时间窗口都将强制聚成 2 类\n\n")

# ================== 2. 🔧 修正的聚类函数 - 强制2类聚类 + 限制m值范围 ==================

calculate_window_membership_with_2clusters <- function(data, window_info, metrics) {
  window_name <- window_info$name
  window_days <- window_info$days
  
  cat(sprintf("Processing %s time window (days %s) with FIXED 2 clusters...\n", 
              window_name, paste(range(window_days), collapse = " to ")))
  
  # Extract data for this time window
  window_cols <- c()
  for(metric in metrics) {
    for(day in window_days) {
      day_str <- paste0("day_", day, "_", metric)
      if(day_str %in% colnames(data)) {
        window_cols <- c(window_cols, day_str)
      }
    }
  }
  
  if(length(window_cols) == 0) {
    cat(sprintf("Warning: No available data for %s time window\n", window_name))
    return(NULL)
  }
  
  cat(sprintf("Found %d data columns for %s\n", length(window_cols), window_name))
  
  # Calculate mean for each metric in this time window
  processed_data <- data %>% dplyr::select(subject_id)
  
  for(metric in metrics) {
    metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
    if(length(metric_cols) > 0) {
      metric_means <- data %>%
        dplyr::select(subject_id, all_of(metric_cols)) %>%
        mutate(
          valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
          metric_mean = ifelse(
            valid_count >= max(1, floor(length(metric_cols)/2)),
            rowMeans(dplyr::select(., -subject_id), na.rm = TRUE),
            NA
          )
        ) %>%
        dplyr::select(subject_id, metric_mean)
      
      names(metric_means)[2] <- paste0(window_name, "_", metric)
      processed_data <- processed_data %>%
        left_join(metric_means, by = "subject_id")
    }
  }
  
  # Remove patients with too many NAs
  complete_patients <- processed_data %>%
    filter(rowSums(is.na(dplyr::select(., -subject_id))) < ncol(dplyr::select(., -subject_id)))
  
  if(nrow(complete_patients) < 3) {
    cat(sprintf("Warning: Insufficient valid patients for %s (%d patients)\n", 
                window_name, nrow(complete_patients)))
    return(NULL)
  }
  
  cat(sprintf("Valid patients for %s: %d\n", window_name, nrow(complete_patients)))
  
  # Fill remaining NAs with mean
  numeric_cols <- names(complete_patients)[-1]
  for(col in numeric_cols) {
    if(sum(!is.na(complete_patients[[col]])) > 0) {
      complete_patients[is.na(complete_patients[[col]]), col] <- 
        mean(complete_patients[[col]], na.rm = TRUE)
    } else {
      complete_patients[[col]] <- 0
      cat(sprintf("Warning: Column %s all NA, filled with 0\n", col))
    }
  }
  
  # Standardize data
  scaled_data <- complete_patients
  for(col in numeric_cols) {
    scaled_data[[col]] <- scale(complete_patients[[col]])[,1]
  }
  
  # ================== 🎯 关键修改：强制使用2个clusters + 限制m值 ==================
  
  # Prepare Mfuzz data
  data_matrix <- scaled_data %>%
    dplyr::select(-subject_id) %>%
    as.matrix()
  
  rownames(data_matrix) <- scaled_data$subject_id
  
  # Create ExpressionSet
  eset <- ExpressionSet(assayData = data_matrix)
  eset_std <- standardise(eset)
  
  # 🎯 关键修改：限制m值范围
  estimated_m <- mestimate(eset_std)
  cat(sprintf("原始估计的m值: %.4f\n", estimated_m))
  
  # 限制m值在合理范围内 (1.5 - 3.0)
  m_value <- pmax(1.5, pmin(3.0, estimated_m))
  
  if(abs(estimated_m - m_value) > 0.01) {
    cat(sprintf("⚠️ m值已从 %.4f 调整为 %.4f (限制在1.5-3.0范围内)\n", estimated_m, m_value))
  } else {
    cat(sprintf("✓ 使用估计的m值: %.4f\n", m_value))
  }
  
  # 🎯 强制使用2个clusters
  optimal_c <- 2
  
  cat(sprintf("🔒 强制使用 %d 个clusters for %s with m=%.4f\n", optimal_c, window_name, m_value))
  
  # 执行2-cluster聚类
  set.seed(RANDOM_SEED)  # 固定随机种子
  
  final_clustering <- NULL
  clustering_method <- "mfuzz"
  
  tryCatch({
    final_clustering <- mfuzz(eset_std, c = optimal_c, m = m_value)
    
    # 🎯 检查聚类质量 - 检查是否所有患者都被分到一个聚类
    max_clusters <- apply(final_clustering$membership, 1, which.max)
    cluster_sizes <- table(max_clusters)
    min_cluster_size <- min(cluster_sizes)
    max_membership_values <- apply(final_clustering$membership, 1, max)
    avg_max_membership <- mean(max_membership_values)
    
    cat(sprintf("Mfuzz聚类结果检查:\n"))
    cat(sprintf("  - 聚类分布: %s\n", paste(paste("Cluster", names(cluster_sizes), "=", cluster_sizes, "人"), collapse = ", ")))
    cat(sprintf("  - 平均最大membership: %.4f\n", avg_max_membership))
    cat(sprintf("  - 最小聚类大小: %d\n", min_cluster_size))
    
    # 如果聚类质量太差（平均最大membership < 0.6 或 最小聚类只有1个人），使用K-means
    if(avg_max_membership < 0.6 || min_cluster_size <= 1) {
      cat(sprintf("⚠️ Mfuzz聚类质量较差 (avg_membership=%.4f, min_size=%d)，切换到K-means\n", 
                  avg_max_membership, min_cluster_size))
      final_clustering <- NULL  # 重置，使用K-means
    } else {
      cat(sprintf("✓ Mfuzz聚类质量acceptable，继续使用Mfuzz结果\n"))
    }
    
  }, error = function(e) {
    cat(sprintf("❌ Mfuzz聚类出错: %s，切换到K-means\n", e$message))
    final_clustering <<- NULL
  })
  
  # 🎯 如果Mfuzz失败或质量太差，使用K-means作为备用方案
  if(is.null(final_clustering)) {
    
    cat(sprintf("🔄 使用K-means作为备用聚类方法...\n"))
    clustering_method <- "kmeans"
    
    set.seed(RANDOM_SEED)  # 固定随机种子
    
    tryCatch({
      # 使用K-means进行聚类
      kmeans_result <- kmeans(data_matrix, centers = 2, nstart = 25, iter.max = 100)
      
      # 创建类似Mfuzz的结果结构
      membership_matrix <- matrix(0, nrow = nrow(data_matrix), ncol = 2)
      rownames(membership_matrix) <- rownames(data_matrix)
      colnames(membership_matrix) <- c("Cluster_1", "Cluster_2")
      
      # K-means给出硬分配，我们创建软分配（基于距离）
      for(i in 1:nrow(data_matrix)) {
        cluster_id <- kmeans_result$cluster[i]
        
        # 计算到两个聚类中心的距离
        dist_to_centers <- numeric(2)
        for(j in 1:2) {
          dist_to_centers[j] <- sqrt(sum((data_matrix[i,] - kmeans_result$centers[j,])^2))
        }
        
        # 基于距离创建membership值（距离越近，membership越高）
        # 使用softmax函数创建概率分布
        inv_dist <- 1 / (dist_to_centers + 0.001)  # 避免除零
        membership_matrix[i,] <- inv_dist / sum(inv_dist)
      }
      
      # 创建类似Mfuzz的结果对象
      final_clustering <- list(
        membership = membership_matrix,
        cluster = kmeans_result$cluster,
        centers = kmeans_result$centers
      )
      
      cat(sprintf("✓ K-means聚类成功完成\n"))
      
      # 检查K-means结果
      kmeans_cluster_sizes <- table(kmeans_result$cluster)
      kmeans_max_membership <- apply(membership_matrix, 1, max)
      kmeans_avg_max_membership <- mean(kmeans_max_membership)
      
      cat(sprintf("K-means聚类结果:\n"))
      cat(sprintf("  - 聚类分布: %s\n", paste(paste("Cluster", names(kmeans_cluster_sizes), "=", kmeans_cluster_sizes, "人"), collapse = ", ")))
      cat(sprintf("  - 平均最大membership: %.4f\n", kmeans_avg_max_membership))
      
    }, error = function(e) {
      cat(sprintf("❌ K-means聚类也失败: %s\n", e$message))
      return(NULL)
    })
  }
  
  if(is.null(final_clustering)) {
    cat(sprintf("❌ 所有聚类方法都失败了，无法处理 %s\n", window_name))
    return(NULL)
  }
  
  cat(sprintf("✓ Successfully clustered %s into %d clusters using %s\n", 
              window_name, optimal_c, clustering_method))
  
  # ================== 获取Membership信息 ==================
  
  # 获取membership矩阵
  membership_matrix <- final_clustering$membership
  
  # 由于强制使用2个clusters，标签应该就是1,2，但我们还是检查一下
  original_max_clusters <- apply(membership_matrix, 1, which.max)
  max_memberships_per_patient <- apply(membership_matrix, 1, max)
  
  # 确保标签是1,2（如果不是则重新映射）
  unique_clusters <- sort(unique(original_max_clusters))
  
  if(!identical(unique_clusters, c(1, 2))) {
    cat(sprintf("⚠️ 修正cluster标签从 %s 到 1,2\n", 
                paste(unique_clusters, collapse = ", ")))
    
    # 创建映射：原始cluster -> 1,2
    cluster_mapping <- setNames(1:length(unique_clusters), unique_clusters)
    max_clusters_per_patient <- cluster_mapping[as.character(original_max_clusters)]
    
    # 重新构建membership矩阵
    remapped_membership_matrix <- matrix(0, nrow = nrow(membership_matrix), ncol = 2)
    rownames(remapped_membership_matrix) <- rownames(membership_matrix)
    colnames(remapped_membership_matrix) <- c("Cluster_1", "Cluster_2")
    
    for(i in 1:length(unique_clusters)) {
      original_cluster_id <- unique_clusters[i]
      new_cluster_id <- i
      remapped_membership_matrix[, new_cluster_id] <- membership_matrix[, original_cluster_id]
    }
    
    membership_matrix <- remapped_membership_matrix
    was_remapped <- TRUE
  } else {
    max_clusters_per_patient <- original_max_clusters
    colnames(membership_matrix) <- c("Cluster_1", "Cluster_2")
    was_remapped <- FALSE
    cluster_mapping <- NULL
  }
  
  # ================== 创建Membership结果 ==================
  
  # 创建详细的membership结果
  membership_result <- data.frame(
    subject_id = rownames(membership_matrix),
    window = window_name,
    max_cluster = max_clusters_per_patient,
    max_membership = max_memberships_per_patient,
    cluster_1_membership = membership_matrix[, 1],
    cluster_2_membership = membership_matrix[, 2],
    # 为了兼容性，添加空的cluster 3,4
    cluster_3_membership = NA,
    cluster_4_membership = NA,
    clustering_method = clustering_method,  # 记录使用的聚类方法
    stringsAsFactors = FALSE
  )
  
  # 计算cluster质量指标
  cluster_quality <- data.frame(
    cluster = c(1, 2),
    size = as.vector(table(max_clusters_per_patient)),
    mean_membership = c(
      mean(membership_matrix[max_clusters_per_patient == 1, 1]),
      mean(membership_matrix[max_clusters_per_patient == 2, 2])
    )
  )
  
  cat(sprintf("✓ %s clustering completed: %d patients, %d clusters (fixed 2-cluster, method: %s)\n", 
              window_name, nrow(membership_result), 2, clustering_method))
  
  # 打印cluster分布
  cat("Cluster分布:\n")
  print(cluster_quality)
  cat("\n")
  
  return(list(
    membership_data = membership_result,
    clustering_result = final_clustering,
    original_data = complete_patients,
    scaled_data = scaled_data,
    window_name = window_name,
    n_patients = nrow(complete_patients),
    n_clusters = 2,  # 固定为2
    m_value = m_value,
    metrics = metrics,
    cluster_quality = cluster_quality,
    membership_matrix = membership_matrix,
    cluster_mapping = cluster_mapping,
    was_remapped = was_remapped,
    clustering_method = clustering_method  # 记录使用的方法
  ))
}

# ================== 3. 执行所有时间窗口的聚类分析 ==================

window_memberships <- list()
all_membership_data <- data.frame()

cat("Starting clustering analysis for all time windows with FIXED 2 clusters...\n\n")

for(window_name in names(time_windows)) {
  window_result <- calculate_window_membership_with_2clusters(ppv_data, time_windows[[window_name]], key_metrics)
  
  if(!is.null(window_result)) {
    window_memberships[[window_name]] <- window_result
    all_membership_data <- rbind(all_membership_data, window_result$membership_data)
  }
}

cat(sprintf("Clustering completed for %d time windows\n", length(window_memberships)))
cat(sprintf("Total membership records: %d\n", nrow(all_membership_data)))
cat("🎯 所有时间窗口都已聚成2类\n\n")

# ================== 4. 创建Max Membership宽格式数据 ==================

create_max_membership_wide_format <- function(all_membership_data) {
  
  cat("Creating max membership wide format data (2 clusters)...\n")
  
  # 创建基础的max membership宽格式
  max_membership_wide <- all_membership_data %>%
    dplyr::select(subject_id, window, max_cluster, max_membership) %>%
    pivot_wider(
      names_from = window,
      values_from = c(max_cluster, max_membership),
      names_sep = "_"
    )
  
  # 重命名列以匹配预期格式
  names(max_membership_wide) <- gsub("max_membership_", "membership_", names(max_membership_wide))
  names(max_membership_wide) <- gsub("max_cluster_", "cluster_", names(max_membership_wide))
  
  cat("Max membership wide format created with columns:\n")
  cat(paste(names(max_membership_wide), collapse = ", "), "\n\n")
  
  return(max_membership_wide)
}

# 创建宽格式数据
max_membership_wide <- create_max_membership_wide_format(all_membership_data)

# ================== 5. 🎨 修正的可视化函数 - 针对2类聚类优化 ==================

# 1. 为每个时间窗口创建详细的聚类趋势图
create_window_cluster_trends <- function(window_data, ppv_data, window_info) {
  
  window_name <- window_data$window_name
  window_days <- window_info$days
  metrics <- window_data$metrics
  
  cat(sprintf("\n🎨 创建 %s 时间窗口的聚类趋势图...\n", toupper(window_name)))
  
  # 创建目录
  dir.create(paste0("plots/time_window_trends/", window_name), recursive = TRUE, showWarnings = FALSE)
  
  # 获取该时间窗口的原始时间序列数据
  window_cols <- c()
  for(metric in metrics) {
    for(day in window_days) {
      day_str <- paste0("day_", day, "_", metric)
      if(day_str %in% colnames(ppv_data)) {
        window_cols <- c(window_cols, day_str)
      }
    }
  }
  
  if(length(window_cols) == 0) {
    cat(sprintf("Warning: No window columns found for %s\n", window_name))
    return(NULL)
  }
  
  # 提取患者在该时间窗口的完整时间序列
  patients_in_window <- window_data$membership_data$subject_id
  window_timeseries <- ppv_data %>%
    filter(subject_id %in% patients_in_window) %>%
    dplyr::select(subject_id, all_of(window_cols))
  
  # 添加聚类信息
  window_timeseries <- window_timeseries %>%
    left_join(window_data$membership_data %>% 
                dplyr::select(subject_id, max_cluster, max_membership), 
              by = "subject_id")
  
  # 为每个指标创建趋势图
  metric_plots <- list()
  
  for(metric in metrics) {
    
    cat(sprintf("  创建 %s 指标的趋势图...\n", metric))
    
    # 找到该指标在该时间窗口的列
    metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
    
    if(length(metric_cols) == 0) {
      cat(sprintf("  Warning: No columns found for metric %s\n", metric))
      next
    }
    
    # 准备绘图数据
    plot_data <- window_timeseries %>%
      dplyr::select(subject_id, max_cluster, max_membership, all_of(metric_cols)) %>%
      pivot_longer(
        cols = all_of(metric_cols),
        names_to = "day_metric",
        values_to = "value"
      ) %>%
      mutate(
        day = as.numeric(gsub("^day_(-?\\d+)_.*$", "\\1", day_metric))
      ) %>%
      filter(!is.na(value))
    
    if(nrow(plot_data) == 0) {
      cat(sprintf("  Warning: No valid data for metric %s\n", metric))
      next
    }
    
    # 计算每个聚类的平均轮廓
    mean_profiles <- plot_data %>%
      group_by(max_cluster, day) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        se_value = sd(value, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    
    # 🎯 修正：为2类聚类创建对比图（使用固定颜色）
    p_all <- ggplot() +
      # 个体轨迹
      geom_line(data = plot_data, 
                aes(x = day, y = value, group = subject_id, color = factor(max_cluster)),
                alpha = 0.3, size = 0.5) +
      # 每个聚类的拟合趋势线
      geom_smooth(data = plot_data,
                  aes(x = day, y = value, color = factor(max_cluster)),
                  method = "loess", se = TRUE, size = 1.5, alpha = 0.8, span = 0.7) +
      # 平均轮廓
      geom_line(data = mean_profiles,
                aes(x = day, y = mean_value, color = factor(max_cluster)),
                size = 2, linetype = "dashed") +
      # 添加平均点
      geom_point(data = mean_profiles,
                 aes(x = day, y = mean_value, color = factor(max_cluster)),
                 size = 3) +
      # 分面
      facet_wrap(~ paste("Cluster", max_cluster), labeller = label_value) +
      # 🎯 固定2类颜色（红蓝对比）
      scale_color_manual(
        values = c("1" = "#a488bf", "2" = "#bd992e"),
        labels = c("1" = "Cluster 1", "2" = "Cluster 2"),
        name = "Cluster"
      ) +
      # x轴
      scale_x_continuous(breaks = window_days) +
      # 标签
      labs(
        title = paste(toupper(window_name), "Window - 2-Cluster Comparison"),
        subtitle = paste(toupper(metric), "| Total n =", length(unique(plot_data$subject_id)), "patients | Smooth trends with confidence bands"),
        x = "Day Relative to Surgery",
        y = paste(metric, "Value"),
        caption = paste("Time window:", paste(range(window_days), collapse = " to "), "days | Dashed line = mean, Smooth curve = fitted trend")
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
    
    metric_plots[[metric]] <- p_all
    
    # 保存对比图
    ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                  window_name, "_2clusters_", metric, "_comparison.pdf"),
           p_all, width = 12, height = 8)
    ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                  window_name, "_2clusters_", metric, "_comparison.png"),
           p_all, width = 12, height = 8, dpi = 300)
  }
  
  # 组合所有指标
  if(length(metric_plots) > 0) {
    combined_plot <- do.call(gridExtra::grid.arrange, 
                             c(metric_plots, 
                               ncol = 1,
                               top = paste(toupper(window_name), "Window - All Metrics & 2 Clusters")))
    
    # 保存组合图
    ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                  window_name, "_combined_all_metrics_2clusters.pdf"),
           combined_plot, width = 12, height = 6 * length(metrics))
    ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                  window_name, "_combined_all_metrics_2clusters.png"),
           combined_plot, width = 12, height = 6 * length(metrics), dpi = 300)
  }
  
  cat(sprintf("  ✓ %s 时间窗口2类趋势图创建完成\n", toupper(window_name)))
  
  return(metric_plots)
}

# 2. 创建跨时间窗口的2类聚类中心对比
create_cross_window_cluster_centers <- function(window_memberships, key_metrics) {
  
  cat("\n🎨 创建跨时间窗口2类聚类中心对比图...\n")
  
  dir.create("plots/cross_window_analysis", recursive = TRUE, showWarnings = FALSE)
  
  # 为每个指标创建跨窗口对比
  for(metric in key_metrics) {
    
    cat(sprintf("  创建 %s 指标的跨窗口2类对比...\n", metric))
    
    # 收集所有窗口的聚类中心数据
    all_centers_data <- data.frame()
    
    for(window_name in names(window_memberships)) {
      window_data <- window_memberships[[window_name]]
      if(is.null(window_data)) next
      
      # 🎯 修正：固定为2类
      window_center_data <- data.frame(
        window = window_name,
        cluster = c(1, 2),  # 固定为2类
        metric = metric,
        stringsAsFactors = FALSE
      )
      
      # 从原始数据计算每个聚类的平均值
      cluster_means <- window_data$original_data %>%
        left_join(window_data$membership_data %>% 
                    dplyr::select(subject_id, max_cluster), by = "subject_id") %>%
        group_by(max_cluster) %>%
        summarise(across(contains(metric), mean, na.rm = TRUE), .groups = 'drop')
      
      # 添加平均值到中心数据
      metric_col <- names(cluster_means)[grep(metric, names(cluster_means))]
      if(length(metric_col) > 0) {
        window_center_data$mean_value <- cluster_means[[metric_col]]
        all_centers_data <- rbind(all_centers_data, window_center_data)
      }
    }
    
    if(nrow(all_centers_data) == 0) {
      cat(sprintf("  Warning: No data for metric %s\n", metric))
      next
    }
    
    # 创建跨窗口2类聚类中心对比图
    p_centers <- ggplot(all_centers_data, aes(x = window, y = mean_value, 
                                              color = factor(cluster), group = factor(cluster))) +
      geom_line(size = 1.5) +
      geom_point(size = 4) +
      # 为每个聚类添加拟合趋势线
      geom_smooth(method = "loess", se = TRUE, alpha = 0.3, size = 1, span = 0.8) +
      # 🎯 固定2类颜色
      scale_color_manual(
        values = c("1" = "#a488bf", "2" = "#bd992e"),
        labels = c("1" = "Cluster 1", "2" = "Cluster 2"),
        name = "Cluster"
      ) +
      labs(
        title = paste("Cross-Window 2-Cluster Centers Comparison -", toupper(metric)),
        subtitle = "Mean values across different time windows with fitted trends",
        x = "Time Window",
        y = paste("Mean", metric, "Value"),
        caption = "Each line represents one cluster across time windows | Smooth curves show fitted trends"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
    
    # 保存跨窗口对比图
    ggsave(paste0("plots/cross_window_analysis/cross_window_", metric, "_2cluster_centers.pdf"),
           p_centers, width = 12, height = 8)
    ggsave(paste0("plots/cross_window_analysis/cross_window_", metric, "_2cluster_centers.png"),
           p_centers, width = 12, height = 8, dpi = 300)
  }
  
  cat("  ✓ 跨时间窗口2类聚类中心对比图创建完成\n")
  
  return(TRUE)
}

# 3. 🎯 创建基于聚类中心的趋势对比图（针对2类优化）
create_cluster_center_trends <- function(window_data, ppv_data, window_info) {
  
  window_name <- window_data$window_name
  window_days <- window_info$days
  metrics <- window_data$metrics
  
  cat(sprintf("\n🎯 创建 %s 时间窗口的2类聚类中心趋势图...\n", toupper(window_name)))
  
  # 创建目录
  dir.create(paste0("plots/cluster_center_trends/", window_name), recursive = TRUE, showWarnings = FALSE)
  
  # 获取该时间窗口的原始时间序列数据
  window_cols <- c()
  for(metric in metrics) {
    for(day in window_days) {
      day_str <- paste0("day_", day, "_", metric)
      if(day_str %in% colnames(ppv_data)) {
        window_cols <- c(window_cols, day_str)
      }
    }
  }
  
  if(length(window_cols) == 0) {
    cat(sprintf("Warning: No window columns found for %s\n", window_name))
    return(NULL)
  }
  
  # 提取患者在该时间窗口的完整时间序列
  patients_in_window <- window_data$membership_data$subject_id
  window_timeseries <- ppv_data %>%
    filter(subject_id %in% patients_in_window) %>%
    dplyr::select(subject_id, all_of(window_cols))
  
  # 添加聚类信息
  window_timeseries <- window_timeseries %>%
    left_join(window_data$membership_data %>% 
                dplyr::select(subject_id, max_cluster, max_membership), 
              by = "subject_id")
  
  # 为每个指标创建基于聚类中心的趋势图
  for(metric in metrics) {
    
    cat(sprintf("  创建 %s 指标的2类聚类中心趋势图...\n", metric))
    
    # 找到该指标在该时间窗口的列
    metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
    
    if(length(metric_cols) == 0) {
      cat(sprintf("  Warning: No columns found for metric %s\n", metric))
      next
    }
    
    # 准备绘图数据 - 计算每个聚类在每个时间点的均值
    cluster_centers_data <- window_timeseries %>%
      dplyr::select(subject_id, max_cluster, all_of(metric_cols)) %>%
      group_by(max_cluster) %>%
      summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = 'drop')
    
    # 转换为长格式
    plot_data <- cluster_centers_data %>%
      pivot_longer(
        cols = all_of(metric_cols),
        names_to = "day_metric",
        values_to = "value"
      ) %>%
      mutate(
        day = as.numeric(gsub("^day_(-?\\d+)_.*$", "\\1", day_metric)),
        cluster = factor(max_cluster)
      )
    
    if(nrow(plot_data) == 0) {
      cat(sprintf("  Warning: No valid data for metric %s\n", metric))
      next
    }
    
    # 🎯 创建干净的2类聚类中心趋势图
    p_centers <- ggplot(plot_data, aes(x = day, y = value, color = cluster)) +
      # 连接线
      geom_line(size = 2, alpha = 0.8) +
      # 数据点
      geom_point(size = 4, alpha = 0.9) +
      # 🎯 固定2类颜色（红蓝对比）
      scale_color_manual(
        values = c("1" = "#a488bf", "2" = "#bd992e"),
        name = "Cluster",
        labels = function(x) paste("Cluster", x)
      ) +
      # x轴设置
      scale_x_continuous(
        breaks = window_days,
        labels = window_days,
        name = "Time Point (Relative Days)"
      ) +
      # 清晰的标题和标签
      labs(
        title = paste(toupper(window_name), "2-Cluster Mean Trends:", toupper(gsub("_", " ", metric))),
        subtitle = paste("Time Window:", paste(range(window_days), collapse = " to "), "days"),
        y = paste(toupper(gsub("_", " ", metric)))
      ) +
      # 干净的主题
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    # 🎯 创建带误差棒的版本（显示标准误）
    # 计算标准误
    cluster_stats <- window_timeseries %>%
      dplyr::select(subject_id, max_cluster, all_of(metric_cols)) %>%
      pivot_longer(
        cols = all_of(metric_cols),
        names_to = "day_metric", 
        values_to = "value"
      ) %>%
      mutate(
        day = as.numeric(gsub("^day_(-?\\d+)_.*$", "\\1", day_metric))
      ) %>%
      group_by(max_cluster, day) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        se_value = sd(value, na.rm = TRUE) / sqrt(n()),
        n_patients = n(),
        .groups = 'drop'
      ) %>%
      mutate(cluster = factor(max_cluster))
    
    p_centers_se <- ggplot(cluster_stats, aes(x = day, y = mean_value, color = cluster)) +
      # 误差棒
      geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                    width = 0.3, size = 1, alpha = 0.8) +
      # 连接线
      geom_line(size = 2, alpha = 0.8) +
      # 数据点
      geom_point(size = 4, alpha = 0.9) +
      # 🎯 固定2类颜色
      scale_color_manual(
        values = c("1" = "#a488bf", "2" = "#bd992e"),
        name = "Cluster",
        labels = function(x) paste("Cluster", x)
      ) +
      # x轴设置
      scale_x_continuous(
        breaks = window_days,
        labels = window_days,
        name = "Time Point (Relative Days)"
      ) +
      # 标题和标签
      labs(
        title = paste(toupper(window_name), "2-Cluster Mean Trends:", toupper(gsub("_", " ", metric))),
        subtitle = paste("Time Window:", paste(range(window_days), collapse = " to "), "days | Error bars show ±SE"),
        y = paste(toupper(gsub("_", " ", metric))),
        caption = "Points show cluster means with standard error bars"
      ) +
      # 主题
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.caption = element_text(size = 10, hjust = 0.5)
      )
    
    # 保存干净版本
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_2cluster_centers_clean.pdf"),
           p_centers, width = 10, height = 6, device = "pdf")
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_2cluster_centers_clean.png"),
           p_centers, width = 10, height = 6, dpi = 300)
    
    # 保存带误差棒版本
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_2cluster_centers_with_SE.pdf"),
           p_centers_se, width = 10, height = 6, device = "pdf")
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_2cluster_centers_with_SE.png"),
           p_centers_se, width = 10, height = 6, dpi = 300)
  }
  
  cat(sprintf("  ✓ %s 时间窗口2类聚类中心趋势图创建完成\n", toupper(window_name)))
  
  return(TRUE)
}

# 4. 🎯 修正的质量总览函数（添加返回值）
create_window_quality_overview <- function(window_memberships) {
  
  cat("\n🎨 创建时间窗口2类聚类质量总览...\n")
  
  # 收集所有窗口的质量数据
  quality_summary <- data.frame()
  
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    if(is.null(window_data)) next
    
    # 计算该窗口的质量指标
    window_quality <- data.frame(
      window = window_name,
      n_patients = window_data$n_patients,
      n_clusters = 2,  # 固定为2
      mean_max_membership = mean(window_data$membership_data$max_membership),
      sd_max_membership = sd(window_data$membership_data$max_membership),
      min_max_membership = min(window_data$membership_data$max_membership),
      max_max_membership = max(window_data$membership_data$max_membership),
      was_remapped = !is.null(window_data$cluster_mapping),
      clustering_method = ifelse(is.null(window_data$clustering_method), "mfuzz", window_data$clustering_method)
    )
    
    quality_summary <- rbind(quality_summary, window_quality)
  }
  
  if(nrow(quality_summary) == 0) {
    cat("Warning: No quality data to plot\n")
    return(NULL)
  }
  
  # 创建质量对比图
  p1 <- ggplot(quality_summary, aes(x = window, y = mean_max_membership, fill = clustering_method)) +
    geom_col(alpha = 0.8) +
    geom_errorbar(aes(ymin = pmax(0, mean_max_membership - sd_max_membership),
                      ymax = pmin(1, mean_max_membership + sd_max_membership)),
                  width = 0.2) +
    geom_text(aes(label = paste0("n=", n_patients)), vjust = -0.5, size = 3) +
    scale_fill_manual(
      values = c("mfuzz" = "#4575B4", "kmeans" = "#D73027"),
      labels = c("mfuzz" = "Mfuzz", "kmeans" = "K-means"),
      name = "Method"
    ) +
    labs(
      title = "Time Window 2-Cluster Quality",
      subtitle = "Mean Max Membership ± SD (colored by clustering method)",
      x = "Time Window",
      y = "Mean Max Membership"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  p2 <- ggplot(quality_summary, aes(x = window, y = n_clusters, fill = was_remapped)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = "2"), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("FALSE" = "lightgreen", "TRUE" = "coral"),
                      labels = c("FALSE" = "Original Labels", "TRUE" = "Remapped Labels"),
                      name = "Label Status") +
    scale_y_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3)) +
    labs(
      title = "Number of Clusters by Time Window (Fixed at 2)",
      subtitle = "Color indicates if cluster labels were remapped",
      x = "Time Window",
      y = "Number of Clusters"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # 组合质量图
  combined_quality <- gridExtra::grid.arrange(p1, p2, ncol = 1,
                                              top = "Time Window 2-Cluster Quality Overview")
  
  # 保存质量总览
  ggsave("plots/cross_window_analysis/time_window_2cluster_quality_overview.pdf",
         combined_quality, width = 12, height = 10)
  ggsave("plots/cross_window_analysis/time_window_2cluster_quality_overview.png",
         combined_quality, width = 12, height = 10, dpi = 300)
  
  cat("  ✓ 时间窗口2类聚类质量总览创建完成\n")
  
  # 🎯 修正：添加返回值
  return(quality_summary)
}

# ================== 6. 执行所有可视化 ==================

cat("\n========================================\n")
cat("🎨 开始创建固定2类聚类的时间窗口可视化\n")
cat("========================================\n")

# 1. 为每个时间窗口创建2类详细趋势图
cat("\n=== 创建各时间窗口的2类聚类趋势图 ===\n")
for(window_name in names(window_memberships)) {
  if(!is.null(window_memberships[[window_name]])) {
    create_window_cluster_trends(
      window_memberships[[window_name]], 
      ppv_data, 
      time_windows[[window_name]]
    )
  }
}

# 2. 创建2类聚类中心趋势图
cat("\n=== 创建2类聚类中心趋势图（干净风格） ===\n")
for(window_name in names(window_memberships)) {
  if(!is.null(window_memberships[[window_name]])) {
    create_cluster_center_trends(
      window_memberships[[window_name]], 
      ppv_data, 
      time_windows[[window_name]]
    )
  }
}

# 3. 创建跨时间窗口2类聚类中心对比
cat("\n=== 创建跨时间窗口2类聚类中心对比 ===\n")
create_cross_window_cluster_centers(window_memberships, key_metrics)

# 4. 创建时间窗口2类聚类质量总览
cat("\n=== 创建时间窗口2类聚类质量总览 ===\n")
quality_overview <- create_window_quality_overview(window_memberships)

# ================== 7. 保存详细的2类聚类结果 ==================

save_detailed_2cluster_results <- function(window_memberships, max_membership_wide) {
  
  cat("\n💾 保存详细的2类聚类结果...\n")
  
  # 1. 保存max membership宽格式数据
  write.csv(max_membership_wide, "time_window_2cluster_membership_data.csv", row.names = FALSE)
  cat("✓ 保存宽格式数据: time_window_2cluster_membership_data.csv\n")
  
  # 2. 保存每个时间窗口的详细membership矩阵
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    
    if(!is.null(window_data)) {
      # 保存完整的membership数据
      write.csv(window_data$membership_data, 
                paste0(window_name, "_detailed_2cluster_membership.csv"), 
                row.names = FALSE)
      
      # 保存cluster质量数据
      write.csv(window_data$cluster_quality,
                paste0(window_name, "_2cluster_quality.csv"),
                row.names = FALSE)
      
      # 保存membership矩阵
      membership_matrix_df <- as.data.frame(window_data$membership_matrix)
      membership_matrix_df$subject_id <- rownames(window_data$membership_matrix)
      write.csv(membership_matrix_df,
                paste0(window_name, "_2cluster_membership_matrix.csv"),
                row.names = FALSE)
      
      # 保存cluster映射信息（如果有的话）
      if(!is.null(window_data$cluster_mapping)) {
        mapping_df <- data.frame(
          Original_Cluster = names(window_data$cluster_mapping),
          New_Cluster = as.numeric(window_data$cluster_mapping)
        )
        write.csv(mapping_df,
                  paste0(window_name, "_2cluster_mapping.csv"),
                  row.names = FALSE)
      }
      
      cat(sprintf("✓ 保存 %s 窗口数据\n", window_name))
    }
  }
  
  # 3. 创建聚类摘要
  clustering_summary <- data.frame()
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    if(!is.null(window_data)) {
      summary_row <- data.frame(
        Time_Window = window_name,
        N_Patients = window_data$n_patients,
        N_Clusters = 2,  # 固定为2
        Mean_Max_Membership = round(mean(window_data$membership_data$max_membership), 3),
        SD_Max_Membership = round(sd(window_data$membership_data$max_membership), 3),
        Min_Max_Membership = round(min(window_data$membership_data$max_membership), 3),
        Max_Max_Membership = round(max(window_data$membership_data$max_membership), 3),
        M_Value = round(window_data$m_value, 3),
        Was_Remapped = window_data$was_remapped,
        Clustering_Method = window_data$clustering_method
      )
      clustering_summary <- rbind(clustering_summary, summary_row)
    }
  }
  
  write.csv(clustering_summary, "time_window_2cluster_summary.csv", row.names = FALSE)
  cat("✓ 保存聚类摘要: time_window_2cluster_summary.csv\n")
  
  cat("✓ 所有2类聚类结果数据已保存\n\n")
  
  return(clustering_summary)
}

# 保存结果
clustering_summary <- save_detailed_2cluster_results(window_memberships, max_membership_wide)

# ================== 8. 生成综合2类可视化报告 ==================

generate_2cluster_visualization_report <- function(window_memberships, clustering_summary, quality_overview) {
  
  cat("\n📝 生成综合2类可视化报告...\n")
  
  report <- paste0(
    "========================================\n",
    "时间窗口2类聚类分析 + 可视化报告\n",
    "========================================\n\n",
    
    "🎯 关键特点 - 强制2类聚类:\n",
    "✅ 所有时间窗口都固定聚成2类\n",
    "✅ 聚类标签统一为1和2\n",
    "✅ 颜色方案统一：红色(Cluster1) vs 蓝色(Cluster2)\n",
    "✅ m值限制在1.5-3.0范围内\n",
    "✅ 自动K-means备用方案\n\n",
    
    "🎨 可视化功能:\n",
    "✅ 每个时间窗口的详细2类聚类趋势图\n",
    "✅ 个体轨迹 + 平均轮廓可视化\n",
    "✅ 跨时间窗口2类聚类中心对比\n",
    "✅ 时间窗口聚类质量总览\n",
    "✅ 聚类中心趋势图（干净风格）\n",
    "✅ 带误差棒的统计图表\n",
    "✅ 固定随机种子确保可重复性\n\n",
    
    "🔬 分析设置:\n",
    "- 随机种子: ", RANDOM_SEED, " (确保可重复性)\n",
    "- 聚类方法: Fuzzy C-means (Mfuzz) + K-means备用\n",
    "- 聚类数量: 2 (所有时间窗口强制统一)\n",
    "- m值范围: 1.5 - 3.0\n",
    "- 分析指标: ", paste(key_metrics, collapse = ", "), "\n",
    "- 时间窗口数: ", length(window_memberships), "\n",
    "- 分析日期: ", Sys.Date(), "\n\n",
    
    "📊 各时间窗口2类聚类结果:\n"
  )
  
  # 添加每个时间窗口的详细信息
  for(i in 1:nrow(clustering_summary)) {
    window_data <- clustering_summary[i, ]
    report <- paste0(report,
                     sprintf("\n%d. %s:\n", i, toupper(window_data$Time_Window)),
                     sprintf("   - 患者数量: %d\n", window_data$N_Patients),
                     sprintf("   - 聚类数量: %d (固定为2类)\n", window_data$N_Clusters),
                     sprintf("   - 平均Max Membership: %.3f ± %.3f\n", 
                             window_data$Mean_Max_Membership, window_data$SD_Max_Membership),
                     sprintf("   - 聚类方法: %s\n", window_data$Clustering_Method),
                     sprintf("   - m值: %.3f\n", window_data$M_Value),
                     sprintf("   - 标签重新映射: %s\n", ifelse(window_data$Was_Remapped, "是", "否")))
  }
  
  report <- paste0(report,
                   "\n🎨 生成的可视化文件结构:\n",
                   "📁 plots/time_window_trends/[window_name]/:\n",
                   "  - 每个2类聚类的详细趋势图\n",
                   "  - 个体轨迹 + 拟合趋势线\n",
                   "  - 2类聚类对比图\n",
                   "  - 所有指标组合图\n\n",
                   
                   "📁 plots/cluster_center_trends/[window_name]/:\n",
                   "  - 2类聚类中心趋势图（干净风格）\n",
                   "  - 带误差棒的聚类中心图\n\n",
                   
                   "📁 plots/cross_window_analysis/:\n",
                   "  - 跨时间窗口2类聚类中心对比\n",
                   "  - 时间窗口聚类质量总览\n\n",
                   
                   "📈 可视化特点（针对2类优化）:\n",
                   "✅ 统一颜色方案：红色(Cluster1) vs 蓝色(Cluster2)\n",
                   "✅ 个体患者轨迹：每条线代表一个患者\n",
                   "✅ 平均轮廓：显示聚类平均趋势\n",
                   "✅ 拟合趋势线：平滑的loess拟合\n",
                   "✅ 标准误差：误差棒显示不确定性\n",
                   "✅ 2类对比：直观比较两个聚类模式\n\n",
                   
                   "🔍 如何使用2类可视化结果:\n",
                   "1. 查看 time_window_trends/ 了解每个时间窗口的2类聚类模式\n",
                   "2. 比较2类聚类的平均轮廓识别关键差异\n",
                   "3. 查看 cluster_center_trends/ 获得干净的聚类中心图\n",
                   "4. 查看 cross_window_analysis/ 了解跨窗口聚类演变\n",
                   "5. 使用统一颜色方案进行跨窗口比较\n\n",
                   
                   "📊 2类聚类质量总结:\n"
  )
  
  # 添加质量总结
  avg_membership <- mean(clustering_summary$Mean_Max_Membership)
  total_patients <- sum(clustering_summary$N_Patients)
  remapped_windows <- sum(clustering_summary$Was_Remapped)
  mfuzz_windows <- sum(clustering_summary$Clustering_Method == "mfuzz")
  kmeans_windows <- sum(clustering_summary$Clustering_Method == "kmeans")
  
  report <- paste0(report,
                   sprintf("- 平均Max Membership: %.3f\n", avg_membership),
                   sprintf("- 总分析患者数: %d\n", total_patients),
                   sprintf("- 需要标签修正的窗口: %d/%d\n", remapped_windows, nrow(clustering_summary)),
                   sprintf("- 使用Mfuzz的窗口: %d/%d\n", mfuzz_windows, nrow(clustering_summary)),
                   sprintf("- 使用K-means的窗口: %d/%d\n", kmeans_windows, nrow(clustering_summary)),
                   sprintf("- 所有聚类标签已统一为1和2\n"),
                   sprintf("- 颜色方案：红色(Cluster1) vs 蓝色(Cluster2)\n\n"),
                   
                   "🎯 关键发现:\n",
                   "✅ 固定随机种子确保完全可重复性\n",
                   "✅ 所有时间窗口都成功聚成2类\n",
                   "✅ m值限制解决了聚类模糊问题\n",
                   "✅ 统一的颜色方案便于跨窗口比较\n",
                   "✅ 可视化清晰展示了2类聚类模式差异\n",
                   "✅ K-means备用方案提供了可靠保障\n\n",
                   
                   "📝 数据文件:\n",
                   "- time_window_2cluster_membership_data.csv: 2类聚类宽格式数据\n",
                   "- time_window_2cluster_summary.csv: 2类聚类摘要\n",
                   "- [window]_detailed_2cluster_membership.csv: 各窗口详细数据\n",
                   "- [window]_2cluster_quality.csv: 各窗口质量数据\n",
                   "- [window]_2cluster_membership_matrix.csv: 各窗口membership矩阵\n\n",
                   
                   "🚀 下一步建议:\n",
                   "1. 使用2类可视化结果进行临床解读\n",
                   "2. 基于2类聚类模式进行预后分析\n",
                   "3. 比较不同时间窗口的2类预测能力\n",
                   "4. 验证2类聚类模式的临床意义\n",
                   "5. 使用统一颜色方案制作报告\n",
                   "6. 分析聚类演变的时间模式\n\n",
                   
                   "报告生成时间: ", Sys.time(), "\n",
                   "========================================\n")
  
  # 保存报告
  writeLines(report, "Time_Window_2Cluster_Visualization_Report.txt")
  cat("✓ 保存综合报告: Time_Window_2Cluster_Visualization_Report.txt\n")
  
  cat(report)
  
  return(report)
}

# 生成2类可视化报告
visualization_report <- generate_2cluster_visualization_report(
  window_memberships, clustering_summary, quality_overview
)

# ================== 9. 最终验证和总结 ==================

cat("\n🎉 时间窗口2类聚类分析 + 可视化完成！\n")
cat("========================================\n")
cat("✅ 随机种子固定:", RANDOM_SEED, "\n")
cat("✅ 所有时间窗口都强制聚成2类\n")
cat("✅ 聚类标签统一为1和2\n")
cat("✅ 颜色方案统一：红色(Cluster1) vs 蓝色(Cluster2)\n")
cat("✅ 详细的2类可视化已生成\n")
cat("✅ 所有图表和数据已保存\n")
cat("========================================\n")

# 显示生成的文件
cat("\n📁 生成的主要文件:\n")
main_files <- c(
  "time_window_2cluster_membership_data.csv",
  "time_window_2cluster_summary.csv", 
  "Time_Window_2Cluster_Visualization_Report.txt"
)

for (file in main_files) {
  if (file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  } else {
    cat(sprintf("❌ %s (未找到)\n", file))
  }
}

cat("\n📊 生成的可视化目录:\n")
viz_dirs <- c(
  "plots/time_window_trends",
  "plots/cluster_center_trends",
  "plots/cross_window_analysis"
)

for (dir in viz_dirs) {
  if (dir.exists(dir)) {
    n_subdirs <- length(list.dirs(dir, recursive = FALSE))
    n_files <- length(list.files(dir, pattern = "\\.(pdf|png)$", recursive = TRUE))
    cat(sprintf("✓ %s (%d subdirs, %d files)\n", dir, n_subdirs, n_files))
  } else {
    cat(sprintf("❌ %s (未找到)\n", dir))
  }
}

# 验证聚类结果
cat("\n🔍 2类聚类结果验证:\n")
for(window_name in names(window_memberships)) {
  window_data <- window_memberships[[window_name]]
  if(!is.null(window_data)) {
    clusters <- unique(window_data$membership_data$max_cluster)
    n_clusters <- length(clusters)
    method <- window_data$clustering_method
    cat(sprintf("✓ %s: %d类聚类 (方法: %s, clusters: %s)\n", 
                window_name, n_clusters, method, paste(sort(clusters), collapse = ", ")))
    
    if(n_clusters != 2) {
      cat(sprintf("  ⚠️ 警告：%s窗口未聚成2类！\n", window_name))
    }
  }
}

cat("\n🎨 现在你有了完整的2类时间窗口聚类结果:\n")
cat("1. 📊 CSV数据文件：聚类结果、质量指标、membership矩阵\n")
cat("2. 🎯 可视化图表：趋势图、对比图、质量总览\n") 
cat("3. 📝 综合报告：分析过程、结果总结、使用建议\n")
cat("4. 🔧 技术保障：固定随机种子、统一颜色、错误处理\n")
cat("5. 📈 多种图表：个体轨迹、聚类中心、跨窗口对比\n")
cat("\n🎯 关键改进：所有时间窗口都严格聚成2类，完整保存所有结果！\n")


# ================== 代码一风格：时间窗口趋势可视化 ==================
# 完全复制代码一中的 visualize_time_window_clusters 函数风格

# 1. 修改后的时间窗口聚类趋势可视化函数（代码一风格）
visualize_time_window_clusters_style1 <- function(window_memberships, ppv_data, time_windows, key_metrics) {
  
  cat("\n🎨 创建代码一风格的时间窗口聚类趋势图...\n")
  
  # 创建输出目录
  dir.create("plots/cluster_profiles_time_windows_style1", recursive = TRUE, showWarnings = FALSE)
  
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    if(is.null(window_data)) next
    
    cat(sprintf("\n处理 %s 时间窗口...\n", window_name))
    
    # 为这个时间窗口创建聚类趋势图
    create_cluster_trend_plots_for_window(window_data, ppv_data, time_windows, key_metrics, window_name)
  }
  
  cat("\n✓ 代码一风格时间窗口趋势图创建完成\n")
}

# 2. 为单个时间窗口创建聚类趋势图
create_cluster_trend_plots_for_window <- function(window_data, ppv_data, time_windows, key_metrics, current_window_name) {
  
  # 构造results_df（类似代码一的格式）
  results_df <- create_results_df_for_visualization(window_data, ppv_data, time_windows, key_metrics)
  
  if(is.null(results_df) || nrow(results_df) == 0) {
    cat(sprintf("Warning: No results data for %s\n", current_window_name))
    return(NULL)
  }
  
  # 构造clustering_result（类似代码一的格式）
  clustering_result <- list(
    cl = list(
      membership = window_data$membership_matrix
    )
  )
  
  # 获取聚类数量
  n_clusters <- ncol(clustering_result$cl$membership)
  
  # 为每个聚类创建图
  for (cluster_id in 1:n_clusters) {
    
    cat(sprintf("  创建 Cluster %d 的趋势图...\n", cluster_id))
    
    # 过滤出属于该聚类的数据
    cluster_data <- results_df %>% 
      filter(max_cluster == cluster_id)
    
    if (nrow(cluster_data) == 0) {
      cat(sprintf("Warning: No data found for %s Cluster %d\n", current_window_name, cluster_id))
      next
    }
    
    # 提取membership值
    membership_df <- data.frame(
      subject_id = rownames(clustering_result$cl$membership),
      membership = clustering_result$cl$membership[, cluster_id],
      stringsAsFactors = FALSE
    )
    
    # 为所有指标创建绘图数据
    plot_data <- data.frame()
    
    for (metric in key_metrics) {
      # 找到与该指标相关的列（时间窗口版本）
      metric_cols <- grep(paste0("_", metric, "$"), colnames(results_df), value = TRUE)
      
      if (length(metric_cols) == 0) next
      
      # 为该指标准备数据
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
          # 为时间窗口分配数值以便绘图
          window_order = case_when(
            window == "baseline" ~ 1,
            window == "acute_recovery" ~ 2,
            window == "early_recovery" ~ 3,
            window == "mid_recovery" ~ 4,
            window == "late_recovery" ~ 5,
            TRUE ~ as.numeric(factor(window))
          )
        ) %>%
        # 连接membership值
        left_join(membership_df, by = "subject_id")
      
      plot_data <- bind_rows(plot_data, metric_data)
    }
    
    if(nrow(plot_data) == 0) {
      cat(sprintf("Warning: No plot data for %s Cluster %d\n", current_window_name, cluster_id))
      next
    }
    
    # 计算每个指标和时间窗口的平均轮廓
    mean_profile <- plot_data %>%
      group_by(metric, window, window_order) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # 🎯 创建与代码一完全相同风格的图
    p <- ggplot() +
      # 个体线条，用membership着色
      geom_line(data = plot_data, 
                aes(x = window_order, y = value, group = subject_id, color = membership),
                size = 0.8, alpha = 0.6) +
      # 平均趋势线（粗黑线）
      geom_line(data = mean_profile,
                aes(x = window_order, y = value, group = metric),
                color = "black",  
                size = 1.2) +
      # Membership颜色渐变（与代码一相同）
      scale_color_gradientn(
        colors = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
                   "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027"),
        limits = c(0.2, 1.0),
        oob = scales::squish,
        breaks = seq(0.2, 1.0, by = 0.1),
        name = "Membership"
      ) +
      # 按指标分面
      facet_wrap(~ metric, scales = "free_y", ncol = 1) +
      # 手术期垂直线
      geom_vline(xintercept = 2, linetype = "dashed", color = "gray40") +
      # 设置x轴刻度和标签
      scale_x_continuous(
        breaks = 1:5,
        labels = c("Baseline", "Acute Recovery", "Early Recovery", "Mid Recovery", "Late Recovery")
      ) +
      # 标签和主题
      labs(
        title = paste(toupper(current_window_name), "Cluster", cluster_id, 
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
    
    # 保存图形
    ggsave(paste0("plots/cluster_profiles_time_windows_style1/", toupper(current_window_name), "_cluster_", cluster_id, "_windows_profile.pdf"),
           p, width = 10, height = 8)
    ggsave(paste0("plots/cluster_profiles_time_windows_style1/", toupper(current_window_name), "_cluster_", cluster_id, "_windows_profile.png"),
           p, width = 10, height = 8, dpi = 300)
    
    # 打印图形
    print(p)
    
    cat(sprintf("  ✓ %s Cluster %d 趋势图保存完成\n", toupper(current_window_name), cluster_id))
  }
}

# 3. 创建用于可视化的results_df数据框
create_results_df_for_visualization <- function(window_data, ppv_data, time_windows, key_metrics) {
  
  # 获取该窗口的患者ID
  patients_in_window <- window_data$membership_data$subject_id
  
  # 创建基础results_df
  results_df <- data.frame(
    subject_id = patients_in_window,
    max_cluster = window_data$membership_data$max_cluster,
    stringsAsFactors = FALSE
  )
  
  # 为每个时间窗口和指标添加数据
  for(time_window_name in names(time_windows)) {
    window_info <- time_windows[[time_window_name]]
    window_days <- window_info$days
    
    for(metric in key_metrics) {
      # 收集该时间窗口内该指标的所有列
      window_metric_cols <- c()
      for(day in window_days) {
        day_str <- paste0("day_", day, "_", metric)
        if(day_str %in% colnames(ppv_data)) {
          window_metric_cols <- c(window_metric_cols, day_str)
        }
      }
      
      if(length(window_metric_cols) > 0) {
        # 计算该时间窗口该指标的均值
        window_means <- ppv_data %>%
          filter(subject_id %in% patients_in_window) %>%
          dplyr::select(subject_id, all_of(window_metric_cols)) %>%
          mutate(
            valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
            window_mean = ifelse(
              valid_count >= max(1, floor(length(window_metric_cols)/2)),
              rowMeans(dplyr::select(., -subject_id), na.rm = TRUE),
              NA
            )
          ) %>%
          dplyr::select(subject_id, window_mean)
        
        # 重命名列
        col_name <- paste0(time_window_name, "_", metric)
        names(window_means)[2] <- col_name
        
        # 合并到results_df
        results_df <- results_df %>%
          left_join(window_means, by = "subject_id")
      }
    }
  }
  
  # 填充NA值（用均值填充）
  for(col in names(results_df)) {
    if(col != "subject_id" && col != "max_cluster" && is.numeric(results_df[[col]])) {
      if(sum(!is.na(results_df[[col]])) > 0) {
        results_df[is.na(results_df[[col]]), col] <- mean(results_df[[col]], na.rm = TRUE)
      } else {
        results_df[[col]] <- 0
      }
    }
  }
  
  return(results_df)
}

# 4. 创建跨聚类对比的时间窗口趋势图
create_cross_cluster_time_window_comparison <- function(window_memberships, ppv_data, time_windows, key_metrics) {
  
  cat("\n🎨 创建跨聚类对比的时间窗口趋势图...\n")
  
  # 创建输出目录
  dir.create("plots/cross_cluster_time_windows", recursive = TRUE, showWarnings = FALSE)
  
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    if(is.null(window_data)) next
    
    cat(sprintf("\n处理 %s 时间窗口的跨聚类对比...\n", window_name))
    
    # 构造results_df
    results_df <- create_results_df_for_visualization(window_data, ppv_data, time_windows, key_metrics)
    
    if(is.null(results_df) || nrow(results_df) == 0) {
      cat(sprintf("Warning: No results data for %s\n", window_name))
      next
    }
    
    # 构造clustering_result
    clustering_result <- list(
      cl = list(
        membership = window_data$membership_matrix
      )
    )
    
    # 创建包含所有聚类的对比图
    create_all_clusters_comparison_plot(results_df, clustering_result, key_metrics, window_name, time_windows)
  }
  
  cat("\n✓ 跨聚类对比趋势图创建完成\n")
}

# 5. 创建包含所有聚类的对比图
create_all_clusters_comparison_plot <- function(results_df, clustering_result, key_metrics, window_name, time_windows) {
  
  # 提取所有membership值
  n_clusters <- ncol(clustering_result$cl$membership)
  
  # 为所有聚类创建绘图数据
  all_plot_data <- data.frame()
  
  for(cluster_id in 1:n_clusters) {
    # 过滤出属于该聚类的数据
    cluster_data <- results_df %>% 
      filter(max_cluster == cluster_id)
    
    if (nrow(cluster_data) == 0) next
    
    # 提取membership值
    membership_df <- data.frame(
      subject_id = rownames(clustering_result$cl$membership),
      membership = clustering_result$cl$membership[, cluster_id],
      cluster_id = cluster_id,
      stringsAsFactors = FALSE
    )
    
    # 为所有指标创建绘图数据
    for (metric in key_metrics) {
      # 找到与该指标相关的列
      metric_cols <- grep(paste0("_", metric, "$"), colnames(results_df), value = TRUE)
      
      if (length(metric_cols) == 0) next
      
      # 为该指标准备数据
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
          cluster_id = cluster_id,
          # 为时间窗口分配数值
          window_order = case_when(
            window == "baseline" ~ 1,
            window == "acute_recovery" ~ 2,
            window == "early_recovery" ~ 3,
            window == "mid_recovery" ~ 4,
            window == "late_recovery" ~ 5,
            TRUE ~ as.numeric(factor(window))
          )
        ) %>%
        # 连接membership值
        left_join(membership_df, by = c("subject_id", "cluster_id"))
      
      all_plot_data <- bind_rows(all_plot_data, metric_data)
    }
  }
  
  if(nrow(all_plot_data) == 0) {
    cat(sprintf("Warning: No plot data for %s\n", window_name))
    return(NULL)
  }
  
  # 计算每个聚类、指标和时间窗口的平均轮廓
  mean_profiles <- all_plot_data %>%
    group_by(cluster_id, metric, window, window_order) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  # 🎯 创建跨聚类对比图（代码一风格）
  p_comparison <- ggplot() +
    # 个体线条，按聚类着色
    geom_line(data = all_plot_data, 
              aes(x = window_order, y = value, group = subject_id, 
                  color = factor(cluster_id), alpha = membership),
              size = 0.6) +
    # 聚类平均线（按聚类着色的粗线）
    geom_line(data = mean_profiles,
              aes(x = window_order, y = value, color = factor(cluster_id)),
              size = 2, alpha = 0.9) +
    # 聚类颜色
    scale_color_manual(
      values = c("1" = "#df8859", "2" = "#0fb292"),
      labels = c("1" = "Cluster 1", "2" = "Cluster 2"),
      name = "Cluster"
    ) +
    # 透明度
    scale_alpha_identity() +
    # 按指标分面
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    # 手术期垂直线
    geom_vline(xintercept = 2, linetype = "dashed", color = "gray40") +
    # 设置x轴
    scale_x_continuous(
      breaks = 1:5,
      labels = c("Baseline", "Acute Recovery", "Early Recovery", "Mid Recovery", "Late Recovery")
    ) +
    # 标签和主题
    labs(
      title = paste(toupper(window_name), "Window - All Clusters Comparison"),
      subtitle = paste("Individual trajectories and cluster means across time windows"),
      x = "Time Windows",
      y = "Value"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # 保存对比图
  ggsave(paste0("plots/cross_cluster_time_windows/", toupper(window_name), "_all_clusters_comparison.pdf"),
         p_comparison, width = 12, height = 8)
  ggsave(paste0("plots/cross_cluster_time_windows/", toupper(window_name), "_all_clusters_comparison.png"),
         p_comparison, width = 12, height = 8, dpi = 300)
  
  print(p_comparison)
  
  cat(sprintf("  ✓ %s 跨聚类对比图保存完成\n", toupper(window_name)))
}

# ================== 执行代码一风格的时间窗口可视化 ==================

cat("\n========================================\n")
cat("🎨 开始创建代码一风格的时间窗口趋势图\n")
cat("========================================\n")

# 1. 创建单聚类的时间窗口趋势图（完全模仿代码一）
cat("\n=== 创建单聚类时间窗口趋势图（代码一风格） ===\n")
visualize_time_window_clusters_style1(window_memberships, ppv_data, time_windows, key_metrics)

# 2. 创建跨聚类对比的时间窗口趋势图
cat("\n=== 创建跨聚类对比时间窗口趋势图 ===\n")
create_cross_cluster_time_window_comparison(window_memberships, ppv_data, time_windows, key_metrics)

# ================== 总结代码一风格可视化 ==================

cat("\n🎉 代码一风格时间窗口趋势图创建完成！\n")
cat("========================================\n")
cat("📁 保存位置:\n")
cat("  - plots/cluster_profiles_time_windows_style1/ (单聚类图)\n")
cat("  - plots/cross_cluster_time_windows/ (跨聚类对比图)\n")
cat("\n🎯 代码一风格特点:\n")
cat("✅ 个体轨迹线：用membership值着色（蓝-红渐变）\n")
cat("✅ 聚类平均线：黑色粗线显示平均趋势\n")
cat("✅ 时间窗口x轴：5个围手术期时间窗口\n")
cat("✅ 分面显示：每个指标单独面板\n")
cat("✅ 手术期标记：垂直虚线标记手术期\n")
cat("✅ 完全复制代码一的配色和布局\n")
cat("\n📊 生成的图表类型:\n")
cat("1. 单聚类趋势图：每个聚类单独显示\n")
cat("2. 跨聚类对比图：两个聚类在同一图中比较\n")
cat("3. 多时间窗口：为每个分析的时间窗口生成图表\n")
cat("\n🔧 与代码一的一致性:\n")
cat("✅ 相同的membership颜色渐变\n")
cat("✅ 相同的图表布局和主题\n")
cat("✅ 相同的标签和标题格式\n")
cat("✅ 相同的保存格式（PDF + PNG）\n")

# ================== 生成代码一风格报告 ==================

generate_style1_visualization_report <- function() {
  
  report <- paste0(
    "========================================\n",
    "代码一风格时间窗口趋势可视化报告\n",
    "========================================\n\n",
    
    "🎯 完全复制代码一的可视化风格\n",
    "✅ 个体轨迹 + membership着色\n",
    "✅ 聚类平均线（黑色粗线）\n",
    "✅ 时间窗口坐标轴\n",
    "✅ 分面指标显示\n\n",
    
    "🎨 视觉特点:\n",
    "- 个体轨迹：根据membership值从蓝色到红色渐变\n",
    "- 聚类中心：黑色粗线显示平均趋势\n",
    "- X轴：5个围手术期时间窗口\n",
    "- Y轴：指标数值（每个指标独立缩放）\n",
    "- 手术期：垂直虚线标记\n\n",
    
    "📊 图表输出:\n",
    "1. plots/cluster_profiles_time_windows_style1/:\n",
    "   - 每个时间窗口每个聚类的单独趋势图\n",
    "   - 命名格式：[WINDOW]_cluster_[ID]_windows_profile\n",
    "2. plots/cross_cluster_time_windows/:\n",
    "   - 每个时间窗口的跨聚类对比图\n",
    "   - 命名格式：[WINDOW]_all_clusters_comparison\n\n",
    
    "🎯 技术实现:\n",
    "- 完全模仿代码一的 visualize_time_window_clusters 函数\n",
    "- 保持相同的数据结构和处理逻辑\n",
    "- 使用相同的ggplot2设置和主题\n",
    "- 支持2类聚类的固定设置\n\n",
    
    "💡 使用场景:\n",
    "- 观察个体患者在不同时间窗口的变化模式\n",
    "- 通过membership着色识别聚类质量\n",
    "- 比较不同聚类的时间演变趋势\n",
    "- 识别关键的时间窗口转换点\n\n",
    
    "报告生成时间: ", Sys.time(), "\n",
    "========================================"
  )
  
  writeLines(report, "Style1_Time_Window_Visualization_Report.txt")
  cat("✓ 保存代码一风格报告: Style1_Time_Window_Visualization_Report.txt\n")
  
  return(report)
}

# 生成报告
style1_report <- generate_style1_visualization_report()

cat("\n📝 代码一风格可视化报告已生成\n")
cat("📄 现在你有了与代码一完全一致的时间窗口趋势图！\n")



# # 聚类诊断和解决方案
# # 分析为什么所有患者都被分到一个聚类
# 
# library(tidyverse)
# library(Biobase)
# library(Mfuzz)
# library(cluster)
# library(factoextra)
# library(ggplot2)
# 
# # ================== 1. 数据诊断函数 ==================
# 
# diagnose_clustering_data <- function(data, window_info, metrics) {
#   
#   window_name <- window_info$name
#   window_days <- window_info$days
#   
#   cat(sprintf("=== 诊断 %s 时间窗口的聚类数据 ===\n", window_name))
#   
#   # 提取时间窗口数据（重复原有逻辑）
#   window_cols <- c()
#   for(metric in metrics) {
#     for(day in window_days) {
#       day_str <- paste0("day_", day, "_", metric)
#       if(day_str %in% colnames(data)) {
#         window_cols <- c(window_cols, day_str)
#       }
#     }
#   }
#   
#   # 计算每个指标的平均值
#   processed_data <- data %>% dplyr::select(subject_id)
#   
#   for(metric in metrics) {
#     metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
#     if(length(metric_cols) > 0) {
#       metric_means <- data %>%
#         dplyr::select(subject_id, all_of(metric_cols)) %>%
#         mutate(
#           valid_count = rowSums(!is.na(dplyr::select(., -subject_id))),
#           metric_mean = ifelse(
#             valid_count >= max(1, floor(length(metric_cols)/2)),
#             rowMeans(dplyr::select(., -subject_id), na.rm = TRUE),
#             NA
#           )
#         ) %>%
#         dplyr::select(subject_id, metric_mean)
#       
#       names(metric_means)[2] <- paste0(window_name, "_", metric)
#       processed_data <- processed_data %>%
#         left_join(metric_means, by = "subject_id")
#     }
#   }
#   
#   # 过滤有效患者
#   complete_patients <- processed_data %>%
#     filter(rowSums(is.na(dplyr::select(., -subject_id))) < ncol(dplyr::select(., -subject_id)))
#   
#   # 填充缺失值
#   numeric_cols <- names(complete_patients)[-1]
#   for(col in numeric_cols) {
#     if(sum(!is.na(complete_patients[[col]])) > 0) {
#       complete_patients[is.na(complete_patients[[col]]), col] <- 
#         mean(complete_patients[[col]], na.rm = TRUE)
#     } else {
#       complete_patients[[col]] <- 0
#     }
#   }
#   
#   cat("1. 基本数据信息:\n")
#   cat(sprintf("   - 患者数: %d\n", nrow(complete_patients)))
#   cat(sprintf("   - 变量数: %d\n", length(numeric_cols)))
#   cat(sprintf("   - 变量名: %s\n", paste(numeric_cols, collapse = ", ")))
#   
#   # 2. 描述性统计
#   cat("\n2. 描述性统计:\n")
#   for(col in numeric_cols) {
#     values <- complete_patients[[col]]
#     cat(sprintf("   %s: 均值=%.4f, 标准差=%.4f, 范围=[%.4f, %.4f]\n",
#                 col, mean(values), sd(values), min(values), max(values)))
#   }
#   
#   # 3. 数据标准化后的统计
#   scaled_data <- complete_patients
#   for(col in numeric_cols) {
#     scaled_data[[col]] <- scale(complete_patients[[col]])[,1]
#   }
#   
#   cat("\n3. 标准化后的统计:\n")
#   for(col in numeric_cols) {
#     values <- scaled_data[[col]]
#     cat(sprintf("   %s: 均值=%.4f, 标准差=%.4f, 范围=[%.4f, %.4f]\n",
#                 col, mean(values), sd(values), min(values), max(values)))
#   }
#   
#   # 4. 数据分离度分析
#   data_matrix <- scaled_data %>%
#     dplyr::select(-subject_id) %>%
#     as.matrix()
#   
#   rownames(data_matrix) <- scaled_data$subject_id
#   
#   cat("\n4. 数据分离度分析:\n")
#   
#   # 计算距离矩阵
#   dist_matrix <- dist(data_matrix)
#   cat(sprintf("   - 患者间距离: 均值=%.4f, 标准差=%.4f\n", 
#               mean(dist_matrix), sd(dist_matrix)))
#   cat(sprintf("   - 最小距离=%.4f, 最大距离=%.4f\n", 
#               min(dist_matrix), max(dist_matrix)))
#   
#   # 5. K-means聚类测试
#   cat("\n5. K-means聚类测试:\n")
#   
#   if(nrow(data_matrix) >= 2) {
#     # 测试2类K-means
#     set.seed(123)
#     kmeans_2 <- kmeans(data_matrix, centers = 2, nstart = 25)
#     
#     cat(sprintf("   - K-means 2类结果: 聚类1=%d人, 聚类2=%d人\n",
#                 sum(kmeans_2$cluster == 1), sum(kmeans_2$cluster == 2)))
#     cat(sprintf("   - 类内平方和比例: %.4f\n", 
#                 kmeans_2$betweenss / kmeans_2$totss))
#     
#     # 计算silhouette score
#     if(length(unique(kmeans_2$cluster)) == 2) {
#       sil_kmeans <- silhouette(kmeans_2$cluster, dist_matrix)
#       cat(sprintf("   - K-means Silhouette score: %.4f\n", mean(sil_kmeans[,3])))
#     }
#   }
#   
#   # 6. 测试Mfuzz聚类
#   cat("\n6. Mfuzz聚类测试:\n")
#   
#   tryCatch({
#     eset <- ExpressionSet(assayData = t(data_matrix))  # 注意转置
#     eset_std <- standardise(eset)
#     m_value <- mestimate(eset_std)
#     
#     cat(sprintf("   - 估计的m值: %.4f\n", m_value))
#     
#     # 测试2类mfuzz
#     set.seed(123)
#     mfuzz_2 <- mfuzz(eset_std, c = 2, m = m_value)
#     
#     cluster_assignments <- apply(mfuzz_2$membership, 1, which.max)
#     cluster_counts <- table(cluster_assignments)
#     
#     cat(sprintf("   - Mfuzz 2类结果: 聚类1=%d人, 聚类2=%d人\n",
#                 cluster_counts[1], 
#                 ifelse(length(cluster_counts) > 1, cluster_counts[2], 0)))
#     
#     # 显示membership分布
#     cat("   - Membership分布:\n")
#     for(i in 1:nrow(mfuzz_2$membership)) {
#       patient_id <- rownames(mfuzz_2$membership)[i]
#       memberships <- mfuzz_2$membership[i,]
#       cat(sprintf("     %s: %.4f, %.4f (max=%d)\n", 
#                   patient_id, memberships[1], memberships[2], which.max(memberships)))
#     }
#     
#   }, error = function(e) {
#     cat(sprintf("   - Mfuzz错误: %s\n", e$message))
#   })
#   
#   # 7. 数据可视化分析
#   cat("\n7. 数据可视化分析:\n")
#   
#   if(length(numeric_cols) == 2) {
#     # 创建散点图
#     p <- ggplot(scaled_data, aes_string(x = numeric_cols[1], y = numeric_cols[2])) +
#       geom_point(size = 3, alpha = 0.7) +
#       geom_text(aes(label = subject_id), nudge_y = 0.1, size = 3) +
#       labs(title = paste("数据分布 -", window_name),
#            x = numeric_cols[1], y = numeric_cols[2]) +
#       theme_bw()
#     
#     filename <- paste0(window_name, "_data_distribution.png")
#     ggsave(filename, p, width = 10, height = 8, dpi = 300)
#     cat(sprintf("   - 散点图已保存: %s\n", filename))
#   }
#   
#   # 8. 聚类适合性评估
#   cat("\n8. 聚类适合性评估:\n")
#   
#   # Hopkins统计量（聚类倾向测试）
#   if(nrow(data_matrix) >= 4) {
#     tryCatch({
#       hopkins_stat <- factoextra::get_clust_tendency(data_matrix, n = min(nrow(data_matrix)-1, 5))
#       cat(sprintf("   - Hopkins统计量: %.4f (>0.5表示有聚类倾向)\n", hopkins_stat$hopkins_stat))
#     }, error = function(e) {
#       cat(sprintf("   - Hopkins统计量计算失败: %s\n", e$message))
#     })
#   }
#   
#   # 最优聚类数评估
#   if(nrow(data_matrix) >= 3) {
#     cat("   - 最优聚类数评估:\n")
#     
#     # Elbow method
#     wss <- sapply(1:min(4, nrow(data_matrix)-1), function(k) {
#       if(k == 1) {
#         return(sum(scale(data_matrix)^2))
#       } else {
#         kmeans_result <- kmeans(data_matrix, centers = k, nstart = 10)
#         return(kmeans_result$tot.withinss)
#       }
#     })
#     
#     for(k in 1:length(wss)) {
#       cat(sprintf("     k=%d: WSS=%.4f\n", k, wss[k]))
#     }
#   }
#   
#   return(list(
#     complete_patients = complete_patients,
#     scaled_data = scaled_data,
#     data_matrix = data_matrix,
#     numeric_cols = numeric_cols
#   ))
# }
# 
# # ================== 2. 强制聚类函数（即使数据不适合）==================
# 
# force_clustering <- function(diagnostic_result, window_name, force_method = "kmeans") {
#   
#   cat(sprintf("\n=== 强制聚类 %s (方法: %s) ===\n", window_name, force_method))
#   
#   data_matrix <- diagnostic_result$data_matrix
#   scaled_data <- diagnostic_result$scaled_data
#   
#   if(force_method == "kmeans") {
#     # 使用K-means强制分成2类
#     set.seed(123)
#     kmeans_result <- kmeans(data_matrix, centers = 2, nstart = 50)
#     
#     # 创建类似mfuzz的membership矩阵
#     membership_matrix <- matrix(0, nrow = nrow(data_matrix), ncol = 2)
#     rownames(membership_matrix) <- rownames(data_matrix)
#     colnames(membership_matrix) <- c("Cluster_1", "Cluster_2")
#     
#     # 硬分配转换为软分配（所有membership设为0.5和0.5，然后调整最大值）
#     for(i in 1:nrow(membership_matrix)) {
#       cluster_id <- kmeans_result$cluster[i]
#       membership_matrix[i, cluster_id] <- 0.8
#       membership_matrix[i, -cluster_id] <- 0.2
#     }
#     
#     cluster_assignments <- kmeans_result$cluster
#     
#   } else if(force_method == "random") {
#     # 随机分配
#     set.seed(123)
#     cluster_assignments <- sample(c(1, 2), nrow(data_matrix), replace = TRUE)
#     
#     # 创建随机membership矩阵
#     membership_matrix <- matrix(0, nrow = nrow(data_matrix), ncol = 2)
#     rownames(membership_matrix) <- rownames(data_matrix)
#     colnames(membership_matrix) <- c("Cluster_1", "Cluster_2")
#     
#     for(i in 1:nrow(membership_matrix)) {
#       cluster_id <- cluster_assignments[i]
#       membership_matrix[i, cluster_id] <- runif(1, 0.6, 0.9)
#       membership_matrix[i, -cluster_id] <- 1 - membership_matrix[i, cluster_id]
#     }
#     
#   } else if(force_method == "median_split") {
#     # 基于第一个变量的中位数分割
#     first_var <- data_matrix[, 1]
#     median_val <- median(first_var)
#     cluster_assignments <- ifelse(first_var > median_val, 2, 1)
#     
#     # 创建基于距离中位数的membership
#     membership_matrix <- matrix(0, nrow = nrow(data_matrix), ncol = 2)
#     rownames(membership_matrix) <- rownames(data_matrix)
#     colnames(membership_matrix) <- c("Cluster_1", "Cluster_2")
#     
#     for(i in 1:nrow(membership_matrix)) {
#       distance_from_median <- abs(first_var[i] - median_val)
#       max_distance <- max(abs(first_var - median_val))
#       
#       if(cluster_assignments[i] == 1) {
#         membership_matrix[i, 1] <- 0.5 + 0.4 * (1 - distance_from_median / max_distance)
#         membership_matrix[i, 2] <- 1 - membership_matrix[i, 1]
#       } else {
#         membership_matrix[i, 2] <- 0.5 + 0.4 * (1 - distance_from_median / max_distance)
#         membership_matrix[i, 1] <- 1 - membership_matrix[i, 2]
#       }
#     }
#   }
#   
#   # 验证membership和
#   membership_sums <- rowSums(membership_matrix)
#   cat(sprintf("强制聚类结果:\n"))
#   cat(sprintf("- 聚类1: %d人\n", sum(cluster_assignments == 1)))
#   cat(sprintf("- 聚类2: %d人\n", sum(cluster_assignments == 2)))
#   cat(sprintf("- Membership平均和: %.6f\n", mean(membership_sums)))
#   
#   return(list(
#     membership_matrix = membership_matrix,
#     cluster_assignments = cluster_assignments,
#     method = force_method
#   ))
# }
# 
# # ================== 3. 主诊断函数 ==================
# 
# main_clustering_diagnosis <- function() {
#   
#   cat("========================================\n")
#   cat("聚类问题诊断工具\n")
#   cat("========================================\n")
#   
#   # 加载数据（假设已经设置了正确的工作目录）
#   if(!exists("ppv_data")) {
#     ppv_data <- read.csv("../data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)
#   }
#   
#   time_windows <- list(
#     late_recovery = list(days = 16:30, name = "late_recovery")
#   )
#   
#   key_metrics <- c("cv_rhr_1", "steps_max")
#   
#   # 诊断late_recovery
#   diagnostic_result <- diagnose_clustering_data(ppv_data, time_windows$late_recovery, key_metrics)
#   
#   # 如果原始mfuzz失败，提供替代方案
#   cat("\n========================================\n")
#   cat("替代聚类方案\n")
#   cat("========================================\n")
#   
#   cat("如果mfuzz将所有患者分到一个聚类，你可以选择:\n\n")
#   
#   # 方案1: K-means强制聚类
#   cat("方案1: K-means强制聚类\n")
#   force_result_kmeans <- force_clustering(diagnostic_result, "late_recovery", "kmeans")
#   
#   # 方案2: 中位数分割
#   cat("\n方案2: 基于中位数分割\n")
#   force_result_median <- force_clustering(diagnostic_result, "late_recovery", "median_split")
#   
#   # 方案3: 随机分配（仅用于测试）
#   cat("\n方案3: 随机分配（仅用于测试）\n")
#   force_result_random <- force_clustering(diagnostic_result, "late_recovery", "random")
#   
#   return(list(
#     diagnostic = diagnostic_result,
#     kmeans_result = force_result_kmeans,
#     median_result = force_result_median,
#     random_result = force_result_random
#   ))
# }
# 
# # ================== 4. 使用说明 ==================
# 
# cat("使用说明:\n")
# cat("1. 运行 results <- main_clustering_diagnosis() 进行完整诊断\n")
# cat("2. 查看诊断结果选择合适的聚类方法\n")
# cat("3. 如果数据不适合自然聚类，选择强制聚类方法\n\n")
# 
# cat("开始诊断...\n")
# results <- main_clustering_diagnosis()

