library(tidyverse)
library(Biobase)
library(Mfuzz)
library(ggplot2)
library(gridExtra)
library(cluster)
library(r4projects)

setwd(get_project_wd())
rm(list = ls())
set.seed(123)

# =============== 数据加载和设置 ===============
cat("开始时间窗口聚类分析...\n")

# 定义时间窗口
time_windows <- list(
  baseline = list(days = -4:-1, name = "baseline"),
  acute_recovery = list(days = 0:3, name = "acute_recovery"),
  early_recovery = list(days = 4:7, name = "early_recovery"),
  mid_recovery = list(days = 8:15, name = "mid_recovery"),
  late_recovery = list(days = 16:30, name = "late_recovery")
)

# 加载数据
ppv_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)
key_metrics <- c("cv_rhr_1", "steps_max")

# 创建输出目录
dir.create("3_data_analysis/7_figures/figure2/time_window_clustering", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/7_figures/figure2/time_window_clustering")

# =============== 🔧 修改后的聚类分析函数 ===============
calculate_window_clustering <- function(data, window_info, metrics) {
  window_name <- window_info$name
  window_days <- window_info$days
  
  cat(sprintf("处理 %s 时间窗口...\n", window_name))
  
  # 提取时间窗口数据
  window_cols <- c()
  for(metric in metrics) {
    for(day in window_days) {
      day_str <- paste0("day_", day, "_", metric)
      if(day_str %in% colnames(data)) {
        window_cols <- c(window_cols, day_str)
      }
    }
  }
  
  if(length(window_cols) == 0) return(NULL)
  
  # 计算每个指标的均值
  processed_data <- data %>% dplyr::select(subject_id)
  
  for(metric in metrics) {
    metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
    if(length(metric_cols) > 0) {
      metric_means <- data %>%
        dplyr::select(subject_id, all_of(metric_cols)) %>%
        mutate(
          metric_mean = rowMeans(dplyr::select(., -subject_id), na.rm = TRUE)
        ) %>%
        dplyr::select(subject_id, metric_mean)
      
      names(metric_means)[2] <- paste0(window_name, "_", metric)
      processed_data <- processed_data %>%
        left_join(metric_means, by = "subject_id")
    }
  }
  
  # 移除NA过多的患者
  complete_patients <- processed_data %>%
    filter(rowSums(is.na(dplyr::select(., -subject_id))) < ncol(dplyr::select(., -subject_id)))
  
  if(nrow(complete_patients) < 5) return(NULL)
  
  # 填充NA
  numeric_cols <- names(complete_patients)[-1]
  for(col in numeric_cols) {
    if(sum(!is.na(complete_patients[[col]])) > 0) {
      complete_patients[is.na(complete_patients[[col]]), col] <- 
        mean(complete_patients[[col]], na.rm = TRUE)
    } else {
      complete_patients[[col]] <- 0
    }
  }
  
  # 标准化数据
  scaled_data <- complete_patients
  for(col in numeric_cols) {
    scaled_data[[col]] <- scale(complete_patients[[col]])[,1]
  }
  
  # 创建ExpressionSet
  data_matrix <- scaled_data %>% dplyr::select(-subject_id) %>% as.matrix()
  rownames(data_matrix) <- scaled_data$subject_id
  
  eset <- ExpressionSet(assayData = data_matrix)
  eset_std <- standardise(eset)
  m_value <- mestimate(eset_std)
  
  # =============== 🎯 完全按照代码二的方式：强制多聚类 ===============
  
  # 确定真正的多cluster数量（代码二方法）
  max_clusters <- min(4, max(2, floor(nrow(complete_patients)/3)))
  
  cat(sprintf("Testing optimal cluster number for %s (range: 2-%d)...\n", 
              window_name, max_clusters))
  
  # 测试不同cluster数量的效果（代码二方法）
  cluster_results <- list()
  silhouette_scores <- numeric()
  
  for(c in 2:max_clusters) {
    set.seed(123)  # 固定随机种子
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
  
  # 选择最佳cluster数量（代码二方法）
  best_c_index <- which.max(silhouette_scores)
  optimal_c <- best_c_index + 1  # 因为从2开始
  
  if(length(cluster_results) == 0 || optimal_c < 2) {
    cat(sprintf("Warning: Could not find optimal clustering for %s, using 2 clusters\n", window_name))
    optimal_c <- 2
    set.seed(123)  # 固定随机种子
    best_clustering <- mfuzz(eset_std, c = optimal_c, m = m_value)
  } else {
    best_clustering <- cluster_results[[paste0("c_", optimal_c)]]$clustering
    cat(sprintf("Selected optimal clusters for %s: %d (Silhouette = %.3f)\n", 
                window_name, optimal_c, max(silhouette_scores)))
  }
  
  # =============== 代码二的聚类标签修正方法 ===============
  
  # 获取membership矩阵
  membership_matrix <- best_clustering$membership
  
  # 计算每个患者的max cluster和max membership
  original_max_clusters <- apply(membership_matrix, 1, which.max)
  max_memberships_per_patient <- apply(membership_matrix, 1, max)
  
  # 🔧 关键修正：使用代码二的标签修正函数
  unique_clusters <- sort(unique(original_max_clusters))
  
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
    max_clusters <- cluster_mapping[as.character(original_max_clusters)]
    
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
    
    membership_matrix <- remapped_membership_matrix
    cat("✓ Cluster重新映射完成\n")
    
  } else {
    cat("✓ Cluster标签已经连续，无需重新映射\n")
    max_clusters <- original_max_clusters
  }
  
  cat("新的cluster分布:", paste(sort(unique(max_clusters)), collapse = ", "), "\n\n")
  
  # 创建详细的membership结果 - 使用修正后的数据
  membership_result <- data.frame(
    subject_id = rownames(membership_matrix),
    window = window_name,
    max_cluster = max_clusters,
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
  actual_clusters <- sort(unique(max_clusters))
  cluster_sizes <- as.vector(table(max_clusters))
  
  cluster_quality <- data.frame(
    cluster = actual_clusters,
    size = cluster_sizes,
    mean_membership = sapply(actual_clusters, function(c) {
      mean(membership_matrix[max_clusters == c, c])
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
    clustering_result = best_clustering,
    original_data = complete_patients,
    window_name = window_name,
    metrics = metrics,
    cluster_quality = cluster_quality
  ))
}

# =============== 执行所有时间窗口聚类 ===============
window_memberships <- list()
all_membership_data <- data.frame()

for(window_name in names(time_windows)) {
  window_result <- calculate_window_clustering(ppv_data, time_windows[[window_name]], key_metrics)
  
  if(!is.null(window_result)) {
    window_memberships[[window_name]] <- window_result
    all_membership_data <- rbind(all_membership_data, window_result$membership_data)
  }
}

# 创建宽格式数据
max_membership_wide <- all_membership_data %>%
  dplyr::select(subject_id, window, max_cluster, max_membership) %>%
  pivot_wider(
    names_from = window,
    values_from = c(max_cluster, max_membership),
    names_sep = "_"
  )

names(max_membership_wide) <- gsub("max_membership_", "membership_", names(max_membership_wide))
names(max_membership_wide) <- gsub("max_cluster_", "cluster_", names(max_membership_wide))

# =============== 可视化函数（保持原样） ===============

# 1. 聚类中心趋势图
create_cluster_center_trends <- function(window_data, ppv_data, window_info) {
  window_name <- window_data$window_name
  window_days <- window_info$days
  metrics <- window_data$metrics
  
  dir.create(paste0("plots/cluster_center_trends/", window_name), recursive = TRUE, showWarnings = FALSE)
  
  window_cols <- c()
  for(metric in metrics) {
    for(day in window_days) {
      day_str <- paste0("day_", day, "_", metric)
      if(day_str %in% colnames(ppv_data)) {
        window_cols <- c(window_cols, day_str)
      }
    }
  }
  
  patients_in_window <- window_data$membership_data$subject_id
  window_timeseries <- ppv_data %>%
    filter(subject_id %in% patients_in_window) %>%
    dplyr::select(subject_id, all_of(window_cols)) %>%
    left_join(window_data$membership_data %>% 
                dplyr::select(subject_id, max_cluster), by = "subject_id")
  
  for(metric in metrics) {
    metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
    if(length(metric_cols) == 0) next
    
    cluster_centers_data <- window_timeseries %>%
      dplyr::select(subject_id, max_cluster, all_of(metric_cols)) %>%
      group_by(max_cluster) %>%
      summarise(across(all_of(metric_cols), mean, na.rm = TRUE), .groups = 'drop')
    
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
    
    p_centers <- ggplot(plot_data, aes(x = day, y = value, color = cluster)) +
      geom_line(size = 2, alpha = 0.8) +
      geom_point(size = 4, alpha = 0.9) +
      scale_color_manual(
        values = c("1" = "#b986f4", "2" = "#be9920", "3" = "#4575B4", "4" = "#D73027"),
        name = "Cluster"
      ) +
      scale_x_continuous(
        breaks = window_days,
        labels = window_days,
        name = "Time Point (Relative Days)"
      ) +
      labs(
        title = paste(toupper(window_name), "Cluster Mean Trends:", toupper(gsub("_", " ", metric))),
        subtitle = paste("Time Window:", paste(range(window_days), collapse = " to "), "days"),
        y = paste(toupper(gsub("_", " ", metric)))
      ) +
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
        panel.grid.minor = element_blank()
      )
    
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_cluster_centers.pdf"),
           p_centers, width = 10, height = 6)
    ggsave(paste0("plots/cluster_center_trends/", window_name, "/", 
                  window_name, "_", metric, "_cluster_centers.png"),
           p_centers, width = 10, height = 6, dpi = 300)
  }
}

# 2. 个体轨迹Profile图
create_window_individual_profiles <- function(window_memberships, ppv_data, time_windows) {
  
  dir.create("plots/window_individual_profiles", recursive = TRUE, showWarnings = FALSE)
  
  for(window_name in names(window_memberships)) {
    window_data <- window_memberships[[window_name]]
    window_info <- time_windows[[window_name]]
    
    if(is.null(window_data)) next
    
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
    
    if(length(window_cols) == 0) next
    
    patients_in_window <- window_data$membership_data$subject_id
    window_timeseries <- ppv_data %>%
      filter(subject_id %in% patients_in_window) %>%
      dplyr::select(subject_id, all_of(window_cols)) %>%
      left_join(window_data$membership_data %>% 
                  dplyr::select(subject_id, max_cluster, max_membership), 
                by = "subject_id")
    
    for(cluster_id in sort(unique(window_data$membership_data$max_cluster))) {
      
      cluster_data <- window_timeseries %>% filter(max_cluster == cluster_id)
      if(nrow(cluster_data) == 0) next
      
      for(metric in metrics) {
        
        metric_cols <- window_cols[grep(paste0("_", metric, "$"), window_cols)]
        if(length(metric_cols) == 0) next
        
        plot_data <- cluster_data %>%
          dplyr::select(subject_id, max_membership, all_of(metric_cols)) %>%
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
        
        mean_profile <- plot_data %>%
          group_by(day) %>%
          summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')
        
        p <- ggplot() +
          geom_line(data = plot_data, 
                    aes(x = day, y = value, group = subject_id, color = max_membership),
                    alpha = 0.7, size = 0.8) +
          geom_line(data = mean_profile,
                    aes(x = day, y = mean_value),
                    color = "black", size = 2) +
          {if(0 %in% window_days) geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 1)} +
          scale_color_gradientn(
            colors = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
                       "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027"),
            limits = c(0.2, 1.0),
            oob = scales::squish,
            breaks = seq(0.2, 1.0, by = 0.1),
            name = "Membership"
          ) +
          scale_x_continuous(
            breaks = window_days,
            labels = window_days,
            name = "Day Relative to Surgery"
          ) +
          labs(
            title = paste("Cluster", cluster_id, paste0("(n=", nrow(cluster_data), ")")),
            subtitle = toupper(gsub("_", " ", metric)),
            y = "Value"
          ) +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
            legend.position = "top",
            legend.key.width = unit(3, "cm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey95"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)
          )
        
        ggsave(paste0("plots/window_individual_profiles/", 
                      window_name, "_cluster_", cluster_id, "_", metric, "_profile.pdf"),
               p, width = 10, height = 8)
        ggsave(paste0("plots/window_individual_profiles/", 
                      window_name, "_cluster_", cluster_id, "_", metric, "_profile.png"),
               p, width = 10, height = 8, dpi = 300)
      }
    }
  }
}

# =============== 执行可视化 ===============
cat("\n创建可视化...\n")

# 1. 聚类中心趋势图
for(window_name in names(window_memberships)) {
  if(!is.null(window_memberships[[window_name]])) {
    create_cluster_center_trends(
      window_memberships[[window_name]], 
      ppv_data, 
      time_windows[[window_name]]
    )
  }
}

# 2. 个体轨迹Profile图
create_window_individual_profiles(window_memberships, ppv_data, time_windows)

# =============== 保存结果 ===============
cat("\n保存结果...\n")

write.csv(max_membership_wide, "time_window_clustering_results.csv", row.names = FALSE)
write.csv(all_membership_data, "time_window_clustering_long.csv", row.names = FALSE)

for(window_name in names(window_memberships)) {
  window_data <- window_memberships[[window_name]]
  if(!is.null(window_data)) {
    write.csv(window_data$membership_data, 
              paste0(window_name, "_clustering_results.csv"), 
              row.names = FALSE)
  }
}

cat("\n✅ 时间窗口聚类分析完成！\n")
cat("🎯 关键改进:\n")
cat("- 动态选择最佳聚类数量 (2-4个)\n")
cat("- 聚类标签自动修正为连续\n")
cat("- Late Recovery强制聚类，不再单一聚类\n")
cat("- 保持原有可视化风格\n")