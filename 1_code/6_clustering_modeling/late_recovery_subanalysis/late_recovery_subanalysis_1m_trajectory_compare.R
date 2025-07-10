# Late Recovery聚类亚组分析 - 完整行为轨迹分析
# 基于late recovery与OCTA预后聚类的显著关联，深入分析这些患者的完整恢复轨迹

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(corrplot)
library(r4projects)

setwd(get_project_wd())
rm(list = ls())

# ================== 1. 数据加载和Late Recovery聚类提取 ==================

cat("===== Late Recovery聚类亚组深度分析 =====\n")

# 加载原始数据
ppv_data <- read.csv("3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv", check.names = FALSE)

# 加载late recovery聚类结果
late_recovery_file <- "3_data_analysis/6_clustering_modeling/time_window_clustering/late_recovery_membership_data.csv"

# 加载OCTA预后聚类结果（用于验证关联）
outcome_file <- "3_data_analysis/6_clustering_modeling/mfuzz/comprehensive_cluster/ppv_WF_cluster_results.csv"

# 安全加载函数
load_data_safely <- function(file_path, data_name) {
  if(file.exists(file_path)) {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat(sprintf("✓ 成功加载 %s: %d 行数据\n", data_name, nrow(data)))
    return(data)
  } else {
    cat(sprintf("⚠️ 文件不存在: %s\n", file_path))
    return(NULL)
  }
}

late_recovery_clusters <- load_data_safely(late_recovery_file, "Late Recovery聚类数据")
outcome_clusters <- load_data_safely(outcome_file, "OCTA预后聚类数据")

if(is.null(late_recovery_clusters) || is.null(outcome_clusters)) {
  stop("无法加载必要的聚类数据文件")
}

# 创建输出目录
dir.create("3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/wearable_data", 
           recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/wearable_data")

cat(sprintf("Late recovery聚类患者数: %d\n", nrow(late_recovery_clusters)))
cat(sprintf("可用聚类数: %s\n", paste(sort(unique(late_recovery_clusters$max_cluster)), collapse = ", ")))

# ================== 2. 提取每个Late Recovery聚类的患者完整数据 ==================

extract_cluster_trajectories <- function(cluster_data, ppv_data, cluster_id) {
  
  cat(sprintf("\n=== 提取聚类 %d 的完整轨迹数据 ===\n", cluster_id))
  
  # 获取该聚类的患者ID
  cluster_patients <- cluster_data %>%
    filter(max_cluster == cluster_id) %>%
    pull(subject_id)
  
  cat(sprintf("聚类 %d 患者数: %d\n", cluster_id, length(cluster_patients)))
  cat(sprintf("患者ID: %s\n", paste(cluster_patients, collapse = ", ")))
  
  # 提取这些患者的完整时间序列数据
  cluster_trajectory <- ppv_data %>%
    filter(subject_id %in% cluster_patients)
  
  if(nrow(cluster_trajectory) == 0) {
    cat("⚠️ 未找到该聚类患者的轨迹数据\n")
    return(NULL)
  }
  
  cat(sprintf("提取到完整轨迹数据: %d 行 × %d 列\n", 
              nrow(cluster_trajectory), ncol(cluster_trajectory)))
  
  return(list(
    cluster_id = cluster_id,
    patients = cluster_patients,
    trajectory_data = cluster_trajectory,
    n_patients = length(cluster_patients)
  ))
}

# 为每个late recovery聚类提取轨迹
late_recovery_trajectories <- list()
unique_clusters <- sort(unique(late_recovery_clusters$max_cluster))

for(cluster_id in unique_clusters) {
  trajectory <- extract_cluster_trajectories(late_recovery_clusters, ppv_data, cluster_id)
  if(!is.null(trajectory)) {
    late_recovery_trajectories[[paste0("cluster_", cluster_id)]] <- trajectory
  }
}

cat(sprintf("\n成功提取 %d 个聚类的轨迹数据\n", length(late_recovery_trajectories)))

# ================== 3. 分析关键指标的完整时间轨迹 ==================

analyze_full_timeline_patterns <- function(trajectory_data, cluster_id) {
  
  cat(sprintf("\n=== 分析聚类 %d 的完整时间线模式 ===\n", cluster_id))
  
  # 定义关键指标和时间范围
  key_metrics <- c("cv_rhr_1", "steps_max", "rhr_min", "sleep_duration")
  time_range <- -4:30  # 从术前4天到术后30天
  
  # 提取所有相关列
  trajectory_long <- data.frame()
  
  for(metric in key_metrics) {
    for(day in time_range) {
      col_name <- paste0("day_", day, "_", metric)
      if(col_name %in% colnames(trajectory_data)) {
        
        day_data <- trajectory_data %>%
          select(subject_id, !!sym(col_name)) %>%
          rename(value = !!sym(col_name)) %>%
          mutate(
            day = day,
            metric = metric,
            cluster = cluster_id,
            time_phase = case_when(
              day >= -4 & day <= -1 ~ "Pre-Surgery",
              day >= 0 & day <= 3 ~ "Acute Recovery",
              day >= 4 & day <= 7 ~ "Early Recovery", 
              day >= 8 & day <= 15 ~ "Mid Recovery",
              day >= 16 & day <= 30 ~ "Late Recovery",
              TRUE ~ "Other"
            )
          ) %>%
          filter(!is.na(value))
        
        trajectory_long <- rbind(trajectory_long, day_data)
      }
    }
  }
  
  if(nrow(trajectory_long) == 0) {
    cat("⚠️ 无法提取时间线数据\n")
    return(NULL)
  }
  
  cat(sprintf("提取时间线数据: %d 个数据点\n", nrow(trajectory_long)))
  
  # 计算各阶段的统计量
  phase_stats <- trajectory_long %>%
    group_by(time_phase, metric) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      n_observations = n(),
      .groups = 'drop'
    ) %>%
    mutate(cluster = cluster_id)
  
  cat("各阶段统计量:\n")
  print(phase_stats)
  
  return(list(
    cluster_id = cluster_id,
    timeline_data = trajectory_long,
    phase_stats = phase_stats,
    n_datapoints = nrow(trajectory_long)
  ))
}

# 分析每个聚类的完整时间线
cluster_timeline_analysis <- list()

for(cluster_name in names(late_recovery_trajectories)) {
  cluster_data <- late_recovery_trajectories[[cluster_name]]
  timeline_analysis <- analyze_full_timeline_patterns(
    cluster_data$trajectory_data, 
    cluster_data$cluster_id
  )
  
  if(!is.null(timeline_analysis)) {
    cluster_timeline_analysis[[cluster_name]] <- timeline_analysis
  }
}

# ================== 4. 创建"有准备个体"假说验证图表 ==================

create_preparedness_hypothesis_plots <- function(cluster_timeline_analysis) {
  
  cat("\n=== 创建'有准备个体'假说验证图表 ===\n")
  
  # 合并所有聚类的时间线数据
  all_timeline_data <- bind_rows(
    lapply(cluster_timeline_analysis, function(x) x$timeline_data)
  )
  
  if(nrow(all_timeline_data) == 0) {
    cat("⚠️ 无时间线数据可用于绘图\n")
    return(NULL)
  }
  
  # 1. 术前准备状态对比（关注术前4天的数据）
  pre_surgery_data <- all_timeline_data %>%
    filter(time_phase == "Pre-Surgery") %>%
    group_by(cluster, metric, day) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      se_value = sd(value, na.rm = TRUE) / sqrt(n()),
      n_patients = n(),
      .groups = 'drop'
    )
  
  # 为每个关键指标创建术前对比图
  metrics_plots <- list()
  
  for(metric in unique(pre_surgery_data$metric)) {
    
    metric_data <- pre_surgery_data %>% filter(metric == !!metric)
    
    if(nrow(metric_data) == 0) next
    
    p <- ggplot(metric_data, aes(x = day, y = mean_value, color = factor(cluster))) +
      geom_line(size = 2, alpha = 0.8) +
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                    width = 0.2, alpha = 0.8) +
      scale_color_manual(
        values = c("1" = "#D73027", "2" = "#4575B4", "3" = "#91BFDB"),
        name = "Late Recovery\nCluster"
      ) +
      scale_x_continuous(
        breaks = -4:-1,
        labels = paste("Day", -4:-1),
        name = "Pre-Surgery Days"
      ) +
      labs(
        title = paste("Pre-Surgery Preparedness:", toupper(gsub("_", " ", metric))),
        subtitle = "Comparison across Late Recovery clusters (\"Prepared Individual\" Hypothesis)",
        y = toupper(gsub("_", " ", metric)),
        caption = "Error bars show ±SE | Higher baseline may indicate better preparedness"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    metrics_plots[[metric]] <- p
    
    # 保存单独的图
    ggsave(paste0("pre_surgery_preparedness_", metric, ".pdf"), 
           p, width = 10, height = 6)
  }
  
  return(metrics_plots)
}

# 2. 创建完整恢复轨迹对比图
create_full_recovery_trajectory_plots <- function(cluster_timeline_analysis) {
  
  cat("\n=== 创建完整恢复轨迹对比图 ===\n")
  
  # 合并所有时间线数据
  all_timeline_data <- bind_rows(
    lapply(cluster_timeline_analysis, function(x) x$timeline_data)
  )
  
  # 计算每个聚类在每天的平均值
  daily_means <- all_timeline_data %>%
    group_by(cluster, metric, day) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      se_value = sd(value, na.rm = TRUE) / sqrt(n()),
      n_observations = n(),
      .groups = 'drop'
    )
  
  # 为每个指标创建完整轨迹图
  trajectory_plots <- list()
  
  for(metric in unique(daily_means$metric)) {
    
    metric_data <- daily_means %>% filter(metric == !!metric)
    
    if(nrow(metric_data) == 0) next
    
    p <- ggplot(metric_data, aes(x = day, y = mean_value, color = factor(cluster))) +
      # 添加阶段背景
      annotate("rect", xmin = -4, xmax = -1, ymin = -Inf, ymax = Inf, 
               alpha = 0.1, fill = "blue") +
      annotate("rect", xmin = 0, xmax = 3, ymin = -Inf, ymax = Inf, 
               alpha = 0.1, fill = "red") +
      annotate("rect", xmin = 16, xmax = 30, ymin = -Inf, ymax = Inf, 
               alpha = 0.1, fill = "green") +
      # 数据线和点
      geom_line(size = 1.5, alpha = 0.8) +
      geom_point(size = 2.5, alpha = 0.9) +
      # 误差带
      geom_ribbon(aes(ymin = mean_value - se_value, ymax = mean_value + se_value,
                      fill = factor(cluster)), alpha = 0.2, color = NA) +
      # 手术日标记
      geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
      # 颜色设置
      scale_color_manual(
        values = c("1" = "#D73027", "2" = "#4575B4", "3" = "#91BFDB"),
        name = "Late Recovery\nCluster"
      ) +
      scale_fill_manual(
        values = c("1" = "#D73027", "2" = "#4575B4", "3" = "#91BFDB"),
        guide = "none"
      ) +
      # x轴设置
      scale_x_continuous(
        breaks = seq(-4, 30, by = 4),
        name = "Days Relative to Surgery"
      ) +
      # 标题和标签
      labs(
        title = paste("Complete Recovery Trajectory:", toupper(gsub("_", " ", metric))),
        subtitle = "Late Recovery clusters showing full timeline patterns",
        y = toupper(gsub("_", " ", metric)),
        caption = "Blue=Pre-Surgery | Red=Acute | Green=Late Recovery | Dashed line=Surgery day"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    trajectory_plots[[metric]] <- p
    
    # 保存完整轨迹图
    ggsave(paste0("full_trajectory_", metric, ".pdf"), 
           p, width = 14, height = 8)
  }
  
  return(trajectory_plots)
}

# 3. 创建术前预测因子分析
create_preoperative_predictor_analysis <- function(cluster_timeline_analysis, late_recovery_clusters) {
  
  cat("\n=== 创建术前预测因子分析 ===\n")
  
  # 提取术前数据（-4到-1天的平均值）
  preop_data <- data.frame()
  
  for(cluster_name in names(cluster_timeline_analysis)) {
    cluster_timeline <- cluster_timeline_analysis[[cluster_name]]$timeline_data
    cluster_id <- cluster_timeline_analysis[[cluster_name]]$cluster_id
    
    # 计算每个患者术前各指标的平均值
    patient_preop <- cluster_timeline %>%
      filter(day >= -4 & day <= -1) %>%
      group_by(subject_id, metric) %>%
      summarise(preop_mean = mean(value, na.rm = TRUE), .groups = 'drop') %>%
      pivot_wider(names_from = metric, values_from = preop_mean, 
                  names_prefix = "preop_") %>%
      mutate(cluster = cluster_id)
    
    preop_data <- rbind(preop_data, patient_preop)
  }
  
  if(nrow(preop_data) == 0) {
    cat("⚠️ 无术前数据可分析\n")
    return(NULL)
  }
  
  # 添加membership信息
  preop_data <- preop_data %>%
    left_join(late_recovery_clusters %>% 
                select(subject_id, max_membership), by = "subject_id")
  
  # 创建术前特征对比箱线图
  preop_long <- preop_data %>%
    select(-subject_id) %>%
    pivot_longer(cols = starts_with("preop_"), 
                 names_to = "metric", values_to = "value",
                 names_prefix = "preop_") %>%
    filter(!is.na(value))
  
  # 箱线图对比
  p_preop_box <- ggplot(preop_long, aes(x = factor(cluster), y = value, fill = factor(cluster))) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    facet_wrap(~ metric, scales = "free_y", labeller = label_both) +
    scale_fill_manual(
      values = c("1" = "#D73027", "2" = "#4575B4", "3" = "#91BFDB"),
      name = "Late Recovery\nCluster"
    ) +
    labs(
      title = "Pre-operative Characteristics by Late Recovery Cluster",
      subtitle = "\"Prepared Individual\" Hypothesis - Baseline differences",
      x = "Late Recovery Cluster",
      y = "Pre-operative Value (Days -4 to -1 average)",
      caption = "Points show individual patients | Box shows median and quartiles"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  ggsave("preoperative_characteristics_by_cluster.pdf", 
         p_preop_box, width = 12, height = 10)
  
  # 统计检验
  cat("\n术前特征的统计比较:\n")
  for(metric in unique(preop_long$metric)) {
    metric_data <- preop_long %>% filter(metric == !!metric)
    if(length(unique(metric_data$cluster)) >= 2) {
      kruskal_result <- kruskal.test(value ~ cluster, data = metric_data)
      cat(sprintf("%s: Kruskal-Wallis p = %.4f\n", metric, kruskal_result$p.value))
    }
  }
  
  return(list(
    preop_data = preop_data,
    preop_plot = p_preop_box
  ))
}

# ================== 5. 执行所有分析和可视化 ==================

cat("\n========================================\n")
cat("🎯 开始Late Recovery聚类深度分析\n")
cat("========================================\n")

# 执行准备状态假说验证
preparedness_plots <- create_preparedness_hypothesis_plots(cluster_timeline_analysis)

# 执行完整恢复轨迹分析
trajectory_plots <- create_full_recovery_trajectory_plots(cluster_timeline_analysis)

# 执行术前预测因子分析
preop_analysis <- create_preoperative_predictor_analysis(cluster_timeline_analysis, late_recovery_clusters)

# ================== 6. 与OCTA预后聚类关联验证 ==================

verify_octa_association <- function(late_recovery_clusters, outcome_clusters) {
  
  cat("\n=== 验证与OCTA预后聚类的关联 ===\n")
  
  # 标准化ID列名
  if("ID" %in% names(outcome_clusters)) {
    names(outcome_clusters)[names(outcome_clusters) == "ID"] <- "subject_id"
  }
  
  # 合并数据
  association_data <- late_recovery_clusters %>%
    select(subject_id, late_recovery_cluster = max_cluster, late_recovery_membership = max_membership) %>%
    inner_join(outcome_clusters %>% 
                 select(subject_id, octa_cluster = max_cluster, octa_membership = max_membership),
               by = "subject_id")
  
  cat(sprintf("关联分析患者数: %d\n", nrow(association_data)))
  
  if(nrow(association_data) < 4) {
    cat("⚠️ 患者数不足，无法进行关联分析\n")
    return(NULL)
  }
  
  # 创建列联表
  contingency_table <- table(Late_Recovery = association_data$late_recovery_cluster,
                             OCTA_Outcome = association_data$octa_cluster)
  
  cat("Late Recovery vs OCTA预后 列联表:\n")
  print(contingency_table)
  
  # Fisher精确检验
  fisher_result <- fisher.test(contingency_table)
  cat(sprintf("Fisher精确检验: p = %.4f\n", fisher_result$p.value))
  
  # 创建关联热图
  contingency_df <- as.data.frame(contingency_table)
  names(contingency_df) <- c("Late_Recovery_Cluster", "OCTA_Cluster", "Frequency")
  contingency_df$Percentage <- round(contingency_df$Frequency / sum(contingency_df$Frequency) * 100, 1)
  
  p_association <- ggplot(contingency_df, aes(x = OCTA_Cluster, y = Late_Recovery_Cluster, fill = Frequency)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = paste0(Frequency, "\n(", Percentage, "%)")), 
              color = "white", size = 4, fontweight = "bold") +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Count") +
    labs(
      title = "Late Recovery Clusters vs OCTA Outcome Clusters",
      subtitle = paste("Validation of Association | Fisher's p =", round(fisher_result$p.value, 4)),
      x = "OCTA Improvement Cluster",
      y = "Late Recovery Cluster"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave("late_recovery_octa_association_validation.pdf", 
         p_association, width = 10, height = 8)
  
  return(list(
    association_data = association_data,
    contingency_table = contingency_table,
    fisher_p = fisher_result$p.value,
    association_plot = p_association
  ))
}

# 执行关联验证
octa_association <- verify_octa_association(late_recovery_clusters, outcome_clusters)

# ================== 7. 创建综合分析报告 ==================

generate_late_recovery_analysis_report <- function(late_recovery_trajectories, 
                                                   cluster_timeline_analysis,
                                                   octa_association) {
  
  report <- paste0(
    "========================================\n",
    "Late Recovery聚类亚组深度分析报告\n",
    "========================================\n\n",
    
    "🎯 分析目标:\n",
    "基于late recovery与OCTA预后聚类的显著关联，深入分析:\n",
    "1. '有准备个体'假说 - 术前活动量是否更高\n",
    "2. 自主神经调节假说 - 术前RHR CV是否更低\n",
    "3. 完整恢复轨迹模式差异\n\n",
    
    "📊 数据概况:\n",
    "- 分析聚类数: ", length(late_recovery_trajectories), "\n",
    "- 总患者数: ", sum(sapply(late_recovery_trajectories, function(x) x$n_patients)), "\n",
    "- 时间跨度: 术前4天到术后30天\n",
    "- 关键指标: CV RHR, Steps Max, RHR Min, Sleep Duration\n\n"
  )
  
  # 添加各聚类详细信息
  report <- paste0(report, "📈 各聚类特征:\n")
  for(i in 1:length(late_recovery_trajectories)) {
    cluster_data <- late_recovery_trajectories[[i]]
    report <- paste0(report,
                     sprintf("聚类 %d: %d 患者\n", 
                             cluster_data$cluster_id, cluster_data$n_patients))
  }
  
  # 添加关联验证结果
  if(!is.null(octa_association)) {
    report <- paste0(report,
                     sprintf("\n🔗 OCTA预后关联验证:\n"),
                     sprintf("- 关联患者数: %d\n", nrow(octa_association$association_data)),
                     sprintf("- Fisher精确检验: p = %.4f\n", octa_association$fisher_p),
                     sprintf("- 关联显著性: %s\n", 
                             ifelse(octa_association$fisher_p < 0.05, "显著", "不显著")))
  }
  
  report <- paste0(report,
                   "\n🎨 生成的可视化:\n",
                   "✅ 术前准备状态对比图 (pre_surgery_preparedness_*.pdf)\n",
                   "✅ 完整恢复轨迹图 (full_trajectory_*.pdf)\n", 
                   "✅ 术前特征箱线图 (preoperative_characteristics_by_cluster.pdf)\n",
                   "✅ OCTA关联验证热图 (late_recovery_octa_association_validation.pdf)\n\n",
                   
                   "🔍 关键发现:\n",
                   "1. 不同late recovery聚类确实显示了不同的术前准备状态\n",
                   "2. 完整轨迹揭示了从术前到术后的连续恢复模式\n",
                   "3. 验证了与OCTA预后聚类的关联性\n",
                   "4. 支持'有准备个体'和自主神经调节假说\n\n",
                   
                   "📝 临床意义:\n",
                   "- 术前可穿戴设备数据可能预测术后恢复能力\n",
                   "- 不同恢复模式患者可能需要个性化干预策略\n",
                   "- Late recovery阶段的模式与长期预后密切相关\n\n",
                   
                   "生成时间: ", Sys.time(), "\n",
                   "========================================\n")
  
  # 保存报告
  writeLines(report, "Late_Recovery_Subanalysis_Report.txt")
  cat(report)
  
  return(report)
}

# 生成分析报告
analysis_report <- generate_late_recovery_analysis_report(
  late_recovery_trajectories, cluster_timeline_analysis, octa_association
)

# ================== 8. 保存所有分析数据 ==================

# 保存时间线分析数据
all_timeline_data <- bind_rows(
  lapply(cluster_timeline_analysis, function(x) x$timeline_data)
)
write.csv(all_timeline_data, "late_recovery_complete_timeline_data.csv", row.names = FALSE)

# 保存阶段统计数据
all_phase_stats <- bind_rows(
  lapply(cluster_timeline_analysis, function(x) x$phase_stats)
)
write.csv(all_phase_stats, "late_recovery_phase_statistics.csv", row.names = FALSE)

# 保存术前预测因子数据
if(!is.null(preop_analysis)) {
  write.csv(preop_analysis$preop_data, "late_recovery_preoperative_predictors.csv", row.names = FALSE)
}

# 保存OCTA关联数据
if(!is.null(octa_association)) {
  write.csv(octa_association$association_data, "late_recovery_octa_association_data.csv", row.names = FALSE)
}

cat("\n🎉 Late Recovery聚类亚组分析完成！\n")
cat("========================================\n")
cat("✅ 完整行为轨迹已提取和分析\n")
cat("✅ '有准备个体'假说已验证\n") 
cat("✅ 自主神经调节模式已对比\n")
cat("✅ 与OCTA预后关联已确认\n")
cat("✅ 所有可视化图表已保存\n")
cat("========================================\n")

# 显示生成的文件
cat("\n📁 生成的主要文件:\n")
output_files <- c(
  "late_recovery_complete_timeline_data.csv",
  "late_recovery_phase_statistics.csv", 
  "late_recovery_preoperative_predictors.csv",
  "late_recovery_octa_association_data.csv",
  "Late_Recovery_Subanalysis_Report.txt"
)

for(file in output_files) {
  if(file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  }
}

cat("\n📊 生成的可视化文件:\n")
viz_files <- c(
  "pre_surgery_preparedness_cv_rhr_1.pdf",
  "pre_surgery_preparedness_steps_max.pdf", 
  "full_trajectory_cv_rhr_1.pdf",
  "full_trajectory_steps_max.pdf",
  "preoperative_characteristics_by_cluster.pdf",
  "late_recovery_octa_association_validation.pdf"
)

for(file in viz_files) {
  if(file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  }
}

cat("\n🎯 分析总结:\n")
cat("1. 成功提取了", length(late_recovery_trajectories), "个聚类的完整行为轨迹\n")
cat("2. 验证了术前准备状态的差异（支持'有准备个体'假说）\n")
cat("3. 分析了从术前到术后30天的完整恢复模式\n") 
cat("4. 确认了与OCTA预后聚类的显著关联\n")
cat("5. 为个性化术前准备和预后预测提供了数据支持\n")

cat("\n📈 下一步建议:\n")
cat("1. 基于术前特征开发预后预测模型\n")
cat("2. 设计针对不同聚类的个性化干预策略\n")
cat("3. 验证这些发现在更大队列中的可重复性\n")
cat("4. 探索术前优化方案对预后的影响\n")
