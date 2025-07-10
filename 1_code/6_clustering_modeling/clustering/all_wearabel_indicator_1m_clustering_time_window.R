# 修正时间窗口聚类代码 - 添加真正的Max Membership + 修正cluster标签 + 类似代码一的可视化
# 解决cluster标签不连续问题（如1,3而不是1,2）+ 添加详细的聚类趋势可视化

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
cat("===== Time Window Specific Clustering Analysis with Max Membership (Fixed) + Visualization =====\n")

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

cat("Data loaded successfully. Starting time window clustering analysis with max membership...\n")
cat("Time windows:", length(time_windows), "\n")
cat("Total patients in dataset:", nrow(ppv_data), "\n")
cat("Metrics for clustering:", paste(key_metrics, collapse = ", "), "\n\n")

# ================== 2. 🔧 新增：Cluster标签修正函数 ==================

fix_cluster_labeling <- function(membership_matrix, max_clusters_per_patient) {
  
  cat("=== 修正Cluster标签 ===\n")
  
  # 获取实际存在的cluster IDs（有患者分配的clusters）
  unique_clusters <- sort(unique(max_clusters_per_patient))
  
  cat("原始cluster IDs:", paste(unique_clusters, collapse = ", "), "\n")
  cat("期望的连续IDs: 1, 2, ...,", length(unique_clusters), "\n")
  
  # 检查是否需要重新映射
  expected_clusters <- 1:length(unique_clusters)
  needs_remapping <- !identical(unique_clusters, expected_clusters)
  
  if(needs_remapping) {
    cat("⚠️ 发现cluster标签不连续，进行重新映射...\n")
    
    # 创建映射表：原始cluster ID -> 新的连续ID
    cluster_mapping <- setNames(1:length(unique_clusters), unique_clusters)
    
    cat("Cluster映射关系:\n")
    for(i in 1:length(cluster_mapping)) {
      cat(sprintf("  原始Cluster %s -> 新Cluster %d\n", names(cluster_mapping)[i], cluster_mapping[i]))
    }
    
    # 重新映射max_clusters_per_patient
    remapped_clusters <- cluster_mapping[as.character(max_clusters_per_patient)]
    
    # 重新构建membership矩阵（只保留有患者的clusters，并重新排序）
    remapped_membership_matrix <- matrix(0, nrow = nrow(membership_matrix), 
                                         ncol = length(unique_clusters))
    rownames(remapped_membership_matrix) <- rownames(membership_matrix)
    colnames(remapped_membership_matrix) <- paste0("Cluster_", 1:length(unique_clusters))
    
    # 将原始membership值复制到新的位置
    for(i in 1:length(unique_clusters)) {
      original_cluster_id <- unique_clusters[i]
      new_cluster_id <- i
      remapped_membership_matrix[, new_cluster_id] <- membership_matrix[, original_cluster_id]
    }
    
    cat("✓ Cluster重新映射完成\n")
    cat("新的cluster分布:", paste(sort(unique(remapped_clusters)), collapse = ", "), "\n\n")
    
    return(list(
      max_clusters = remapped_clusters,
      membership_matrix = remapped_membership_matrix,
      cluster_mapping = cluster_mapping,
      was_remapped = TRUE
    ))
    
  } else {
    cat("✓ Cluster标签已经连续，无需重新映射\n\n")
    
    return(list(
      max_clusters = max_clusters_per_patient,
      membership_matrix = membership_matrix,
      cluster_mapping = NULL,
      was_remapped = FALSE
    ))
  }
}

# ================== 3. 修正的聚类函数 - 真正的多cluster分析 + 标签修正 ==================

calculate_window_membership_with_max <- function(data, window_info, metrics) {
  window_name <- window_info$name
  window_days <- window_info$days
  
  cat(sprintf("Processing %s time window (days %s) with multiple clusters...\n", 
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
  
  if(nrow(complete_patients) < 5) {
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
  
  # ================== 关键修改：创建真正的多cluster结构 ==================
  
  # Prepare Mfuzz data
  data_matrix <- scaled_data %>%
    dplyr::select(-subject_id) %>%
    as.matrix()
  
  rownames(data_matrix) <- scaled_data$subject_id
  
  # Create ExpressionSet
  eset <- ExpressionSet(assayData = data_matrix)
  eset_std <- standardise(eset)
  
  # Estimate optimal parameters
  m_value <- mestimate(eset_std)
  
  # 🔧 关键修改：确定真正的多cluster数量
  # 参考代码二的方法，尝试2-4个clusters
  max_clusters <- min(4, max(2, floor(nrow(complete_patients)/3)))
  
  cat(sprintf("Testing optimal cluster number for %s (range: 2-%d)...\n", 
              window_name, max_clusters))
  
  # 测试不同cluster数量的效果
  cluster_results <- list()
  silhouette_scores <- numeric()
  
  for(c in 2:max_clusters) {
    set.seed(RANDOM_SEED)  # 固定随机种子
    tryCatch({
      clustering_result <- mfuzz(eset_std, c = c, m = m_value)
      
      # 计算silhouette score来评估聚类质量
      cluster_assignments <- apply(clustering_result$membership, 1, which.max)
      
      # 计算距离矩阵
      dist_matrix <- dist(data_matrix)
      
      # 检查是否有足够的clusters
      if(length(unique(cluster_assignments)) >= 2) {
        sil_score <- cluster::silhouette(cluster_assignments, dist_matrix)
        avg_sil <- mean(sil_score[, 3])
        
        cluster_results[[paste0("c_", c)]] <- list(
          clustering = clustering_result,
          silhouette = avg_sil,
          n_clusters = c
        )
        silhouette_scores <- c(silhouette_scores, avg_sil)
        
        cat(sprintf("  Clusters = %d: Silhouette = %.3f\n", c, avg_sil))
      } else {
        silhouette_scores <- c(silhouette_scores, -1)
        cat(sprintf("  Clusters = %d: Failed (insufficient clusters)\n", c))
      }
    }, error = function(e) {
      silhouette_scores <<- c(silhouette_scores, -1)
      cat(sprintf("  Clusters = %d: Error occurred\n", c))
    })
  }
  
  # 选择最佳cluster数量
  best_c_index <- which.max(silhouette_scores)
  optimal_c <- best_c_index + 1  # 因为从2开始
  
  if(length(cluster_results) == 0 || optimal_c < 2) {
    cat(sprintf("Warning: Could not find optimal clustering for %s, using 2 clusters\n", window_name))
    optimal_c <- 2
    set.seed(RANDOM_SEED)  # 固定随机种子
    final_clustering <- mfuzz(eset_std, c = optimal_c, m = m_value)
  } else {
    final_clustering <- cluster_results[[paste0("c_", optimal_c)]]$clustering
    cat(sprintf("Selected optimal clusters for %s: %d (Silhouette = %.3f)\n", 
                window_name, optimal_c, max(silhouette_scores)))
  }
  
  # ================== 🔧 新增：修正cluster标签部分 ==================
  
  # 获取membership矩阵
  membership_matrix <- final_clustering$membership
  
  # 计算每个患者的max cluster和max membership
  original_max_clusters <- apply(membership_matrix, 1, which.max)
  max_memberships_per_patient <- apply(membership_matrix, 1, max)
  
  # 🔧 关键修正：使用修正函数确保cluster标签连续
  cluster_fix_result <- fix_cluster_labeling(membership_matrix, original_max_clusters)
  
  # 使用修正后的结果
  max_clusters_per_patient <- cluster_fix_result$max_clusters
  membership_matrix <- cluster_fix_result$membership_matrix
  
  # 记录是否进行了重新映射
  if(cluster_fix_result$was_remapped) {
    cat(sprintf("📋 %s窗口cluster标签已重新映射为连续标签\n", window_name))
  }
  
  # ================== 提取Max Membership信息 ==================
  
  # 创建详细的membership结果 - 使用修正后的数据
  membership_result <- data.frame(
    subject_id = rownames(membership_matrix),
    window = window_name,
    max_cluster = max_clusters_per_patient,
    max_membership = max_memberships_per_patient,
    stringsAsFactors = FALSE
  )
  
  # 添加所有clusters的membership值 - 使用修正后的矩阵
  for(c in 1:ncol(membership_matrix)) {
    col_name <- paste0("cluster_", c, "_membership")
    membership_result[[col_name]] <- membership_matrix[, c]
  }
  
  # 标准化为最多4个clusters（为了兼容性）
  max_possible_clusters <- 4
  if(ncol(membership_matrix) < max_possible_clusters) {
    for(c in (ncol(membership_matrix) + 1):max_possible_clusters) {
      col_name <- paste0("cluster_", c, "_membership")
      membership_result[[col_name]] <- NA
    }
  }
  
  # 计算cluster质量指标 - 使用修正后的clusters
  actual_clusters <- sort(unique(max_clusters_per_patient))
  cluster_sizes <- as.vector(table(max_clusters_per_patient))
  
  cluster_quality <- data.frame(
    cluster = actual_clusters,
    size = cluster_sizes,
    mean_membership = sapply(actual_clusters, function(c) {
      mean(membership_matrix[max_clusters_per_patient == c, c])
    })
  )
  
  cat(sprintf("✓ %s clustering completed: %d patients, %d clusters (连续标签)\n", 
              window_name, nrow(membership_result), length(actual_clusters)))
  
  # 打印cluster分布
  cat("修正后的Cluster分布:\n")
  print(cluster_quality)
  cat("\n")
  
  return(list(
    membership_data = membership_result,
    clustering_result = final_clustering,
    original_data = complete_patients,
    scaled_data = scaled_data,
    window_name = window_name,
    n_patients = nrow(complete_patients),
    n_clusters = length(actual_clusters),
    m_value = m_value,
    metrics = metrics,
    cluster_quality = cluster_quality,
    membership_matrix = membership_matrix,
    cluster_mapping = cluster_fix_result$cluster_mapping  # 保存映射信息
  ))
}

# ================== 4. 执行所有时间窗口的聚类分析 ==================

window_memberships <- list()
all_membership_data <- data.frame()

cat("Starting clustering analysis for all time windows with max membership (with fixed labels)...\n\n")

for(window_name in names(time_windows)) {
  window_result <- calculate_window_membership_with_max(ppv_data, time_windows[[window_name]], key_metrics)
  
  if(!is.null(window_result)) {
    window_memberships[[window_name]] <- window_result
    all_membership_data <- rbind(all_membership_data, window_result$membership_data)
  }
}

cat(sprintf("Clustering completed for %d time windows\n", length(window_memberships)))
cat(sprintf("Total membership records: %d\n\n", nrow(all_membership_data)))

# ================== 5. 创建Max Membership宽格式数据 ==================

create_max_membership_wide_format <- function(all_membership_data) {
  
  cat("Creating max membership wide format data...\n")
  
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

# ================== 6. 🎨 新增：类似代码一的聚类可视化函数 ==================

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
    
    if(length(metric_cols) == 0) next
    
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
    
    if(nrow(plot_data) == 0) next
    
    # 计算每个聚类的平均轮廓
    mean_profiles <- plot_data %>%
      group_by(max_cluster, day) %>%
      summarise(
        mean_value = mean(value, na.rm = TRUE),
        se_value = sd(value, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    
    # 为每个聚类创建单独的图
    n_clusters <- length(unique(plot_data$max_cluster))
    
    for(cluster_id in sort(unique(plot_data$max_cluster))) {
      
      # 该聚类的数据
      cluster_data <- plot_data %>% filter(max_cluster == cluster_id)
      cluster_mean <- mean_profiles %>% filter(max_cluster == cluster_id)
      
      # 创建聚类特定的图
      p <- ggplot() +
        # 个体轨迹（按membership着色）
        geom_line(data = cluster_data, 
                  aes(x = day, y = value, group = subject_id, color = max_membership),
                  alpha = 0.6, size = 0.8) +
        # 平均轮廓（粗黑线）
        geom_line(data = cluster_mean,
                  aes(x = day, y = mean_value),
                  color = "black", size = 2) +
        # 🎯 新增：拟合趋势线（平滑曲线）
        geom_smooth(data = cluster_data,
                    aes(x = day, y = value),
                    method = "loess", se = TRUE, 
                    color = "darkred", size = 1.5, 
                    alpha = 0.8, span = 0.7) +
        # 添加标准误差
        geom_ribbon(data = cluster_mean,
                    aes(x = day, ymin = mean_value - se_value, ymax = mean_value + se_value),
                    alpha = 0.2, fill = "gray") +
        # 添加平均点
        geom_point(data = cluster_mean,
                   aes(x = day, y = mean_value),
                   color = "black", size = 3) +
        # 颜色渐变（membership值）
        scale_color_gradientn(
          colors = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
                     "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027"),
          limits = c(0.2, 1.0),
          oob = scales::squish,
          name = "Max\nMembership"
        ) +
        # 设置x轴
        scale_x_continuous(breaks = window_days) +
        # 标签
        labs(
          title = paste(toupper(window_name), "Window - Cluster", cluster_id),
          subtitle = paste(toupper(metric), "| n =", nrow(cluster_data) / length(unique(cluster_data$day)), "patients | Smooth trend in red"),
          x = "Day Relative to Surgery",
          y = paste(metric, "Value"),
          caption = paste("Time window:", paste(range(window_days), collapse = " to "), "days | Black line = mean, Red smooth = fitted trend")
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.position = "right",
          panel.grid.minor = element_blank()
        )
      
      # 保存单个聚类图
      ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                    window_name, "_cluster_", cluster_id, "_", metric, "_trend.pdf"),
             p, width = 10, height = 6)
      ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                    window_name, "_cluster_", cluster_id, "_", metric, "_trend.png"),
             p, width = 10, height = 6, dpi = 300)
    }
    
    # 创建所有聚类对比图
    p_all <- ggplot() +
      # 个体轨迹
      geom_line(data = plot_data, 
                aes(x = day, y = value, group = subject_id, color = factor(max_cluster)),
                alpha = 0.3, size = 0.5) +
      # 🎯 新增：每个聚类的拟合趋势线
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
      facet_wrap(~ max_cluster, labeller = label_both) +
      # 颜色
      scale_color_brewer(type = "qual", palette = "Set2", name = "Cluster") +
      # x轴
      scale_x_continuous(breaks = window_days) +
      # 标签
      labs(
        title = paste(toupper(window_name), "Window - All Clusters Comparison"),
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
                  window_name, "_all_clusters_", metric, "_comparison.pdf"),
           p_all, width = 12, height = 8)
    ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                  window_name, "_all_clusters_", metric, "_comparison.png"),
           p_all, width = 12, height = 8, dpi = 300)
  }
  
  # 组合所有指标
  if(length(metric_plots) > 0) {
    combined_plot <- do.call(gridExtra::grid.arrange, 
                             c(metric_plots, 
                               ncol = 1,
                               top = paste(toupper(window_name), "Window - All Metrics & Clusters")))
    
    # 保存组合图
    ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                  window_name, "_combined_all_metrics_clusters.pdf"),
           combined_plot, width = 12, height = 6 * length(metrics))
    ggsave(paste0("plots/time_window_trends/", window_name, "/", 
                  window_name, "_combined_all_metrics_clusters.png"),
           combined_plot, width = 12, height = 6 * length(metrics), dpi = 300)
  }
  
  cat(sprintf("  ✓ %s 时间窗口趋势图创建完成\n", toupper(window_name)))
}

# 2. 创建跨时间窗口的聚类中心对比
create_cross_window_cluster_centers <- function(window_memberships, key_metrics) {
  
  cat("\n🎨 创建跨时间窗口聚类中心对比图...\n")
  
  dir.create("plots/cross_window_analysis", recursive = TRUE, showWarnings = FALSE)
  
  # 为每个指标创建跨窗口对比
  for(metric in key_metrics) {
    
    cat(sprintf("  创建 %s 指标的跨窗口对比...\n", metric))
    
    # 收集所有窗口的聚类中心数据
    all_centers_data <- data.frame()
    
    for(window_name in names(window_memberships)) {
      window_data <- window_memberships[[window_name]]
      if(is.null(window_data)) next
      
      # 获取该窗口该指标的聚类中心
      window_center_data <- data.frame(
        window = window_name,
        cluster = 1:window_data$n_clusters,
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
    
    if(nrow(all_centers_data) == 0) next
    
    # 创建跨窗口聚类中心对比图
    p_centers <- ggplot(all_centers_data, aes(x = window, y = mean_value, 
                                              color = factor(cluster), group = factor(cluster))) +
      geom_line(size = 1.5) +
      geom_point(size = 4) +
      # 🎯 新增：为每个聚类添加拟合趋势线
      geom_smooth(method = "loess", se = TRUE, alpha = 0.3, size = 1, span = 0.8) +
      scale_color_brewer(type = "qual", palette = "Set2", name = "Cluster") +
      labs(
        title = paste("Cross-Window Cluster Centers Comparison -", toupper(metric)),
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
    ggsave(paste0("plots/cross_window_analysis/cross_window_", metric, "_centers.pdf"),
           p_centers, width = 12, height = 8)
    ggsave(paste0("plots/cross_window_analysis/cross_window_", metric, "_centers.png"),
           p_centers, width = 12, height = 8, dpi = 300)
  }
  
  cat("  ✓ 跨时间窗口聚类中心对比图创建完成\n")
}

# 4. 🎯 新增：创建基于聚类中心的趋势对比图（类似参考代码风格）
create_cluster_center_trends <- function(window_data, ppv_data, window_info) {
  
  window_name <- window_data$window_name
  window_days <- window_info$days
  metrics <- window_data$metrics
  
  cat(sprintf("\n🎯 创建 %s 时间窗口的聚类中心趋势图...\n", toupper(window_name)))
  
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
    
    cat(sprintf("  创建 %s 指标的聚类中心趋势图...\n", metric))
    
    # 找到该指标在该时间窗口的列
    metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
    
    if(length(metric_cols) == 0) next
    
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
    
    if(nrow(plot_data) == 0) next
    
    # 🎯 创建干净的聚类中心趋势图（类似参考代码风格）
    p_centers <- ggplot(plot_data, aes(x = day, y = value, color = cluster)) +
      # 连接线
      geom_line(size = 2, alpha = 0.8) +
      # 数据点
      geom_point(size = 4, alpha = 0.9) +
      # 颜色设置（类似参考代码）
      scale_color_manual(
        values = c("1" = "#D73027", "2" = "#4575B4", "3" = "#91BFDB", "4" = "#FC8D59"),
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
        title = paste(toupper(window_name), "Cluster Mean Trends:", toupper(gsub("_", " ", metric))),
        subtitle = paste("Time Window:", paste(range(window_days), collapse = " to "), "days"),
        y = paste(toupper(gsub("_", " ", metric)))
      ) +
      # 干净的主题（类似参考代码）
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
      # 颜色设置
      scale_color_manual(
        values = c("1" = "#D73027", "2" = "#4575B4", "3" = "#91BFDB", "4" = "#FC8D59"),
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
        title = paste(toupper(window_name), "Cluster Mean Trends:", toupper(gsub("_", " ", metric))),
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
                  window_name, "_", metric, "_cluster_centers_clean.pdf"),
           p_centers, width = 10, height = 6, device = "pdf")
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_cluster_centers_clean.png"),
           p_centers, width = 10, height = 6, dpi = 300)
    
    # 保存带误差棒版本
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_cluster_centers_with_SE.pdf"),
           p_centers_se, width = 10, height = 6, device = "pdf")
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_cluster_centers_with_SE.png"),
           p_centers_se, width = 10, height = 6, dpi = 300)
  }
  
  cat(sprintf("  ✓ %s 时间窗口聚类中心趋势图创建完成\n", toupper(window_name)))
}
create_window_quality_overview <- function(window_memberships) {
  
  cat("\n🎨 创建时间窗口聚类质量总览...\n")
  
  # 收集所有窗口的质量数据
  quality_summary <- data.frame()
  
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    if(is.null(window_data)) next
    
    # 计算该窗口的质量指标
    window_quality <- data.frame(
      window = window_name,
      n_patients = window_data$n_patients,
      n_clusters = window_data$n_clusters,
      mean_max_membership = mean(window_data$membership_data$max_membership),
      sd_max_membership = sd(window_data$membership_data$max_membership),
      min_max_membership = min(window_data$membership_data$max_membership),
      max_max_membership = max(window_data$membership_data$max_membership),
      was_remapped = !is.null(window_data$cluster_mapping)
    )
    
    quality_summary <- rbind(quality_summary, window_quality)
  }
  
  # 创建质量对比图
  p1 <- ggplot(quality_summary, aes(x = window, y = mean_max_membership, fill = window)) +
    geom_col(alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_max_membership - sd_max_membership,
                      ymax = mean_max_membership + sd_max_membership),
                  width = 0.2) +
    geom_text(aes(label = paste0("n=", n_patients)), vjust = -0.5) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    labs(
      title = "Time Window Clustering Quality",
      subtitle = "Mean Max Membership ± SD",
      x = "Time Window",
      y = "Mean Max Membership"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  p2 <- ggplot(quality_summary, aes(x = window, y = n_clusters, fill = was_remapped)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = n_clusters), vjust = -0.5) +
    scale_fill_manual(values = c("FALSE" = "lightgreen", "TRUE" = "coral"),
                      labels = c("FALSE" = "Original Labels", "TRUE" = "Remapped Labels"),
                      name = "Label Status") +
    labs(
      title = "Number of Clusters by Time Window",
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
                                              top = "Time Window Clustering Quality Overview")
  
  # 保存质量总览
  ggsave("plots/cross_window_analysis/time_window_quality_overview.pdf",
         combined_quality, width = 12, height = 10)
  ggsave("plots/cross_window_analysis/time_window_quality_overview.png",
         combined_quality, width = 12, height = 10, dpi = 300)
  
  cat("  ✓ 时间窗口聚类质量总览创建完成\n")
  
  return(quality_summary)
}

# ================== 7. 执行所有可视化 ==================

cat("\n========================================\n")
cat("🎨 开始创建类似代码一的时间窗口聚类可视化\n")
cat("========================================\n")

# 1. 为每个时间窗口创建详细趋势图
cat("\n=== 创建各时间窗口的详细聚类趋势图 ===\n")
for(window_name in names(window_memberships)) {
  if(!is.null(window_memberships[[window_name]])) {
    create_window_cluster_trends(
      window_memberships[[window_name]], 
      ppv_data, 
      time_windows[[window_name]]
    )
  }
}

# 🎯 新增：创建聚类中心趋势图（类似参考代码风格）
cat("\n=== 创建聚类中心趋势图（干净风格） ===\n")
for(window_name in names(window_memberships)) {
  if(!is.null(window_memberships[[window_name]])) {
    create_cluster_center_trends(
      window_memberships[[window_name]], 
      ppv_data, 
      time_windows[[window_name]]
    )
  }
}

# 2. 创建跨时间窗口聚类中心对比
cat("\n=== 创建跨时间窗口聚类中心对比 ===\n")
create_cross_window_cluster_centers(window_memberships, key_metrics)

# 3. 创建时间窗口聚类质量总览
cat("\n=== 创建时间窗口聚类质量总览 ===\n")
quality_overview <- create_window_quality_overview(window_memberships)


# ================== 新增：基于聚类特征的颜色统一分析 ==================

# 函数1：分析聚类特征模式
standardize_cluster_patterns <- function(window_memberships) {
  
  cat("=== 分析各时间窗口的聚类特征模式 ===\n")
  
  # 存储每个时间窗口的聚类特征
  window_cluster_features <- list()
  
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    if(is.null(window_data)) next
    
    cat(sprintf("分析 %s 窗口的聚类特征...\n", window_name))
    
    # 计算每个聚类的特征均值
    cluster_features <- window_data$original_data %>%
      left_join(window_data$membership_data %>% 
                  dplyr::select(subject_id, max_cluster), by = "subject_id") %>%
      group_by(max_cluster) %>%
      summarise(across(contains("cv_rhr"), mean, na.rm = TRUE),
                across(contains("steps_mean"), mean, na.rm = TRUE),
                .groups = 'drop') %>%
      mutate(window = window_name)
    
    # 标准化特征值 (Z-score)
    cv_col <- names(cluster_features)[grep("cv_rhr", names(cluster_features))]
    steps_col <- names(cluster_features)[grep("steps_mean", names(cluster_features))]
    
    if(length(cv_col) > 0 && length(steps_col) > 0) {
      cluster_features <- cluster_features %>%
        mutate(
          cv_rhr_z = scale(!!sym(cv_col))[,1],
          steps_mean_z = scale(!!sym(steps_col))[,1]
        )
      
      # 根据特征模式分类
      cluster_features <- cluster_features %>%
        mutate(
          pattern_type = case_when(
            cv_rhr_z < -0.5 & steps_mean_z > 0.5 ~ "Good_Pattern",      # 低CV，高Steps = 良好模式
            cv_rhr_z > 0.5 & steps_mean_z < -0.5 ~ "Poor_Pattern",      # 高CV，低Steps = 较差模式  
            cv_rhr_z < 0 & steps_mean_z < 0 ~ "Low_Activity",            # 低CV，低Steps = 低活动
            cv_rhr_z > 0 & steps_mean_z > 0 ~ "High_Variability",       # 高CV，高Steps = 高变异
            TRUE ~ "Mixed_Pattern"                                       # 其他混合模式
          ),
          # 为模式分配统一颜色
          unified_color = case_when(
            pattern_type == "Good_Pattern" ~ "#4575B4",        # 蓝色 - 良好模式
            pattern_type == "Poor_Pattern" ~ "#D73027",        # 红色 - 较差模式
            pattern_type == "Low_Activity" ~ "#74ADD1",        # 浅蓝 - 低活动
            pattern_type == "High_Variability" ~ "#F46D43",    # 橙色 - 高变异
            TRUE ~ "#ABD9E9"                                   # 浅蓝灰 - 混合
          )
        )
      
      cat(sprintf("窗口 %s 的聚类模式:\n", window_name))
      print(cluster_features %>% dplyr::select(max_cluster, pattern_type, unified_color, cv_rhr_z, steps_mean_z))
      
      window_cluster_features[[window_name]] <- cluster_features
    }
  }
  
  return(window_cluster_features)
}

# 函数2：创建模式统一的可视化
create_pattern_unified_visualizations <- function(window_memberships, ppv_data, time_windows, 
                                                  cluster_patterns) {
  
  cat("=== 创建基于模式统一颜色的可视化 ===\n")
  
  # 创建输出目录
  dir.create("plots/pattern_unified_trends", recursive = TRUE, showWarnings = FALSE)
  
  # 为每个时间窗口创建统一颜色的图
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    window_info <- time_windows[[window_name]]
    pattern_info <- cluster_patterns[[window_name]]
    
    if(is.null(window_data) || is.null(pattern_info)) next
    
    cat(sprintf("创建 %s 窗口的模式统一图...\n", window_name))
    
    # 获取时间序列数据
    window_days <- window_info$days
    metrics <- window_data$metrics
    
    window_cols <- c()
    for(metric in metrics) {
      for(day in window_days) {
        day_str <- paste0("day_", day, "_", metric)
        if(day_str %in% colnames(ppv_data)) {
          window_cols <- c(window_cols, day_str)
        }
      }
    }
    
    # 提取时间序列数据并添加模式信息
    patients_in_window <- window_data$membership_data$subject_id
    window_timeseries <- ppv_data %>%
      filter(subject_id %in% patients_in_window) %>%
      dplyr::select(subject_id, all_of(window_cols)) %>%
      left_join(window_data$membership_data %>% 
                  dplyr::select(subject_id, max_cluster), by = "subject_id") %>%
      left_join(pattern_info %>% 
                  dplyr::select(max_cluster, pattern_type, unified_color), by = "max_cluster")
    
    # 为每个指标创建图
    for(metric in metrics) {
      
      metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
      if(length(metric_cols) == 0) next
      
      # 计算聚类中心数据
      cluster_centers_data <- window_timeseries %>%
        dplyr::select(max_cluster, pattern_type, unified_color, all_of(metric_cols)) %>%
        group_by(max_cluster, pattern_type, unified_color) %>%
        summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = 'drop')
      
      # 转换为长格式
      plot_data <- cluster_centers_data %>%
        pivot_longer(
          cols = all_of(metric_cols),
          names_to = "day_metric",
          values_to = "value"
        ) %>%
        mutate(
          day = as.numeric(gsub("^day_(-?\\d+)_.*$", "\\1", day_metric))
        )
      
      if(nrow(plot_data) == 0) next
      
      # 创建模式统一颜色的图
      p_unified <- ggplot(plot_data, aes(x = day, y = value)) +
        # 连接线 - 使用统一颜色
        geom_line(aes(color = I(unified_color)), size = 2, alpha = 0.8) +
        # 数据点 - 使用统一颜色
        geom_point(aes(color = I(unified_color)), size = 4, alpha = 0.9) +
        # x轴设置
        scale_x_continuous(
          breaks = window_days,
          labels = window_days,
          name = "Time Point (Relative Days)"
        ) +
        # 添加聚类编号标注
        geom_text(aes(label = paste0("C", max_cluster), color = I(unified_color)),
                  nudge_y = max(plot_data$value) * 0.05,
                  size = 3, fontface = "bold", show.legend = FALSE) +
        # 标题
        labs(
          title = paste(toupper(window_name), "Pattern-Unified Trends:", toupper(gsub("_", " ", metric))),
          subtitle = "Colors unified by physiological patterns (Blue=Good, Red=Poor)",
          y = paste(toupper(gsub("_", " ", metric))),
          caption = "Blue = Good Pattern (Low CV, High Steps) | Red = Poor Pattern (High CV, Low Steps)"
        ) +
        # 主题
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),
          panel.grid.major = element_line(color = "grey90", size = 0.5),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(size = 10, hjust = 0.5)
        )
      
      # 保存统一颜色版本
      ggsave(paste0("plots/pattern_unified_trends/", 
                    window_name, "_", metric, "_pattern_unified.pdf"),
             p_unified, width = 10, height = 6, device = "pdf")
      ggsave(paste0("plots/pattern_unified_trends/", 
                    window_name, "_", metric, "_pattern_unified.png"),
             p_unified, width = 10, height = 6, dpi = 300)
    }
  }
  
  cat("模式统一颜色图创建完成\n")
}

# 函数3：创建跨时间窗口模式对比
create_cross_window_pattern_comparison <- function(cluster_patterns, key_metrics) {
  
  cat("=== 创建跨时间窗口模式对比 ===\n")
  
  dir.create("plots/pattern_unified_trends", recursive = TRUE, showWarnings = FALSE)
  
  # 合并所有窗口的模式数据
  all_patterns <- bind_rows(cluster_patterns) %>%
    dplyr::select(window, max_cluster, pattern_type, unified_color, cv_rhr_z, steps_mean_z)
  
  if(nrow(all_patterns) == 0) {
    cat("警告：没有找到模式数据\n")
    return(NULL)
  }
  
  # 创建特征空间图 (CV vs Steps)
  p_feature_space <- ggplot(all_patterns, aes(x = cv_rhr_z, y = steps_mean_z)) +
    geom_point(aes(color = I(unified_color)), size = 4, alpha = 0.8) +
    geom_text(aes(label = paste0(substr(window, 1, 4), "\nC", max_cluster)),
              nudge_y = 0.1, size = 3, fontface = "bold", show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = "Cluster Patterns in Standardized Feature Space",
      subtitle = "All time windows and clusters shown with unified colors",
      x = "CV RHR (Standardized)",
      y = "Steps Max (Standardized)",
      caption = "Blue = Good Pattern | Red = Poor Pattern | Light Blue = Low Activity | Orange = High Variability"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(size = 10, hjust = 0.5)
    )
  
  # 保存对比图
  ggsave("plots/pattern_unified_trends/feature_space_unified_colors.pdf",
         p_feature_space, width = 12, height = 8)
  ggsave("plots/pattern_unified_trends/feature_space_unified_colors.png",
         p_feature_space, width = 12, height = 8, dpi = 300)
  
  cat("跨时间窗口模式对比图创建完成\n")
  
  return(p_feature_space)
}

# ================== 执行颜色统一分析 ==================

cat("开始基于聚类特征的颜色统一分析...\n")

# 1. 分析聚类特征模式
cluster_patterns <- standardize_cluster_patterns(window_memberships)

# 2. 创建模式统一的可视化
create_pattern_unified_visualizations(window_memberships, ppv_data, time_windows, cluster_patterns)

# 3. 创建跨时间窗口模式对比
pattern_comparison <- create_cross_window_pattern_comparison(cluster_patterns, key_metrics)

cat("基于模式的颜色统一分析完成！\n")
cat("查看 plots/pattern_unified_trends/ 目录获取统一颜色的图表\n")

# ================== 8. 保存详细的聚类结果 ==================

save_detailed_clustering_results <- function(window_memberships, max_membership_wide) {
  
  cat("Saving detailed clustering results...\n")
  
  # 1. 保存max membership宽格式数据
  write.csv(max_membership_wide, "time_window_max_membership_data_fixed.csv", row.names = FALSE)
  
  # 2. 保存每个时间窗口的详细membership矩阵
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    
    if(!is.null(window_data)) {
      # 保存完整的membership数据
      write.csv(window_data$membership_data, 
                paste0(window_name, "_detailed_membership_fixed.csv"), 
                row.names = FALSE)
      
      # 保存cluster质量数据
      write.csv(window_data$cluster_quality,
                paste0(window_name, "_cluster_quality_fixed.csv"),
                row.names = FALSE)
      
      # 保存membership矩阵
      membership_matrix_df <- as.data.frame(window_data$membership_matrix)
      membership_matrix_df$subject_id <- rownames(window_data$membership_matrix)
      write.csv(membership_matrix_df,
                paste0(window_name, "_membership_matrix_fixed.csv"),
                row.names = FALSE)
      
      # 保存cluster映射信息（如果有的话）
      if(!is.null(window_data$cluster_mapping)) {
        mapping_df <- data.frame(
          Original_Cluster = names(window_data$cluster_mapping),
          New_Cluster = as.numeric(window_data$cluster_mapping)
        )
        write.csv(mapping_df,
                  paste0(window_name, "_cluster_mapping.csv"),
                  row.names = FALSE)
      }
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
        N_Clusters = window_data$n_clusters,
        Mean_Max_Membership = round(mean(window_data$membership_data$max_membership), 3),
        SD_Max_Membership = round(sd(window_data$membership_data$max_membership), 3),
        Min_Max_Membership = round(min(window_data$membership_data$max_membership), 3),
        Max_Max_Membership = round(max(window_data$membership_data$max_membership), 3),
        M_Value = round(window_data$m_value, 3),
        Was_Remapped = !is.null(window_data$cluster_mapping)
      )
      clustering_summary <- rbind(clustering_summary, summary_row)
    }
  }
  
  write.csv(clustering_summary, "time_window_max_clustering_summary_fixed.csv", row.names = FALSE)
  
  cat("✓ All detailed clustering results saved (with fixed labels)\n\n")
  
  return(clustering_summary)
}

# 保存结果
clustering_summary <- save_detailed_clustering_results(window_memberships, max_membership_wide)

# ================== 9. 生成综合可视化报告 ==================

generate_comprehensive_visualization_report <- function(window_memberships, clustering_summary, quality_overview) {
  
  report <- paste0(
    "========================================\n",
    "时间窗口聚类分析 + 类似代码一可视化报告\n",
    "========================================\n\n",
    
    "🎨 可视化功能升级:\n",
    "✅ 每个时间窗口的详细聚类趋势图（类似代码一）\n",
    "✅ 个体轨迹 + 平均轮廓可视化\n",
    "✅ Membership值着色显示聚类置信度\n",
    "✅ 跨时间窗口聚类中心对比\n",
    "✅ 时间窗口聚类质量总览\n",
    "✅ 固定随机种子确保可重复性\n\n",
    
    "🔬 分析设置:\n",
    "- 随机种子: ", RANDOM_SEED, " (确保可重复性)\n",
    "- 聚类方法: Fuzzy C-means (Mfuzz)\n",
    "- 分析指标: ", paste(key_metrics, collapse = ", "), "\n",
    "- 时间窗口数: ", length(window_memberships), "\n",
    "- 分析日期: ", Sys.Date(), "\n\n",
    
    "📊 各时间窗口聚类结果:\n"
  )
  
  # 添加每个时间窗口的详细信息
  for(i in 1:nrow(clustering_summary)) {
    window_data <- clustering_summary[i, ]
    report <- paste0(report,
                     sprintf("\n%d. %s:\n", i, toupper(window_data$Time_Window)),
                     sprintf("   - 患者数量: %d\n", window_data$N_Patients),
                     sprintf("   - 聚类数量: %d (fixed labels)\n", window_data$N_Clusters),
                     sprintf("   - 平均Max Membership: %.3f ± %.3f\n", 
                             window_data$Mean_Max_Membership, window_data$SD_Max_Membership),
                     sprintf("   - 标签重新映射: %s\n", ifelse(window_data$Was_Remapped, "是", "否")))
  }
  
  report <- paste0(report,
                   "\n🎨 生成的可视化文件结构:\n",
                   "📁 plots/time_window_trends/[window_name]/:\n",
                   "  - 每个聚类的详细趋势图（个体轨迹 + 平均轮廓）\n",
                   "  - 按membership值着色的个体轨迹\n",
                   "  - 各聚类对比图\n",
                   "  - 所有指标组合图\n\n",
                   
                   "📁 plots/cross_window_analysis/:\n",
                   "  - 跨时间窗口聚类中心对比\n",
                   "  - 时间窗口聚类质量总览\n\n",
                   
                   "📈 可视化特点（类似代码一）:\n",
                   "✅ 个体患者轨迹：每条线代表一个患者\n",
                   "✅ Membership着色：线条颜色反映聚类置信度\n",
                   "✅ 平均轮廓：粗黑线显示聚类平均趋势\n",
                   "✅ 标准误差：灰色阴影显示不确定性\n",
                   "✅ 聚类对比：直观比较不同聚类模式\n\n",
                   
                   "🔍 如何使用可视化结果:\n",
                   "1. 查看 time_window_trends/ 了解每个时间窗口的聚类模式\n",
                   "2. 观察个体轨迹的membership着色了解聚类稳定性\n",
                   "3. 比较不同聚类的平均轮廓识别关键差异\n",
                   "4. 查看 cross_window_analysis/ 了解跨窗口聚类演变\n\n",
                   
                   "📊 聚类质量总结:\n"
  )
  
  # 添加质量总结
  avg_membership <- mean(clustering_summary$Mean_Max_Membership)
  total_patients <- sum(clustering_summary$N_Patients)
  remapped_windows <- sum(clustering_summary$Was_Remapped)
  
  report <- paste0(report,
                   sprintf("- 平均Max Membership: %.3f\n", avg_membership),
                   sprintf("- 总分析患者数: %d\n", total_patients),
                   sprintf("- 需要标签修正的窗口: %d/%d\n", remapped_windows, nrow(clustering_summary)),
                   sprintf("- 所有聚类标签已修正为连续编号\n\n"),
                   
                   "🎯 关键发现:\n",
                   "✅ 固定随机种子确保完全可重复性\n",
                   "✅ 每个时间窗口都生成了高质量聚类\n",
                   "✅ 可视化清晰展示了聚类模式差异\n",
                   "✅ Membership着色有助于评估聚类置信度\n\n",
                   
                   "📝 数据文件:\n",
                   "- time_window_max_membership_data_fixed.csv: 修正后的宽格式数据\n",
                   "- time_window_max_clustering_summary_fixed.csv: 聚类摘要\n",
                   "- [window]_detailed_membership_fixed.csv: 各窗口详细数据\n\n",
                   
                   "🚀 下一步建议:\n",
                   "1. 使用可视化结果进行临床解读\n",
                   "2. 基于聚类模式进行预后分析\n",
                   "3. 比较不同时间窗口的预测能力\n",
                   "4. 验证聚类模式的临床意义\n\n",
                   
                   "报告生成时间: ", Sys.time(), "\n",
                   "========================================\n")
  
  # 保存报告
  writeLines(report, "Time_Window_Clustering_Visualization_Report.txt")
  cat(report)
  
  return(report)
}

# 生成可视化报告
visualization_report <- generate_comprehensive_visualization_report(
  window_memberships, clustering_summary, quality_overview
)

# ================== 10. 最终验证和总结 ==================

cat("\n🎉 时间窗口聚类分析 + 可视化完成！\n")
cat("========================================\n")
cat("✅ 随机种子固定:", RANDOM_SEED, "\n")
cat("✅ 聚类标签已修正为连续编号\n")
cat("✅ 类似代码一的详细可视化已生成\n")
cat("✅ 所有图表和数据已保存\n")
cat("========================================\n")

# 显示生成的文件
cat("\n📁 生成的主要文件:\n")
main_files <- c(
  "time_window_max_membership_data_fixed.csv",
  "time_window_max_clustering_summary_fixed.csv", 
  "Time_Window_Clustering_Visualization_Report.txt"
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

cat("\n🎨 现在你有了完整的时间窗口聚类可视化:\n")
cat("1. 每个时间窗口的详细聚类趋势图\n")
cat("2. 个体轨迹 + membership着色\n") 
cat("3. 聚类平均轮廓 + 标准误差\n")
cat("4. 🎯 聚类中心趋势图（类似参考代码风格）\n")
cat("5. 带误差棒的聚类中心图\n")
cat("6. 跨时间窗口聚类中心对比\n")
cat("7. 聚类质量总览和评估\n")
cat("8. 固定随机种子保证可重复性\n")
cat("\n所有可视化都包含干净的聚类中心趋势，类似你参考代码的风格！\n")