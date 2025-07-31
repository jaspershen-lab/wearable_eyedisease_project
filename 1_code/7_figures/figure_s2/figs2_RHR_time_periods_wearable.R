# ================================================================================
# RHR数据可视化 - 仅使用步数≤1的数据
# 基于处理好的RHR数据，专注于静息心率模式分析
# ================================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(r4projects)

setwd(get_project_wd())
rm(list = ls())

# ================== 1. 加载RHR数据 ==================

# 加载处理好的RHR数据
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/combined_rhr_summary.rda")
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/time_period_rhr_results.rda")
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/rhr_detailed_summary.rda")
# 加载RHR心率数据
load("3_data_analysis/2_data_analysis/RHR/heart_rate_data_RHR.rda")

# 加载可穿戴设备聚类结果
wearable_clusters <- read.csv("3_data_analysis/6_clustering_modeling/time_window_clustering/time_window_2cluster_membership_data.csv")


# 创建输出目录
dir.create("3_data_analysis/7_figures/figure_s2/rhr_visualization_steps1", recursive = TRUE)
setwd("3_data_analysis/7_figures/figure_s2/rhr_visualization_steps1")

cat("========================================\n")
cat("RHR数据可视化 - 仅步数≤1\n")
cat("========================================\n")

# ================== 2. 确定聚类人群 ==================

# 获取聚类分析中的患者ID
clustering_patients <- wearable_clusters %>%
  filter(!is.na(subject_id)) %>%
  pull(subject_id) %>%
  unique() %>%
  toupper()

# 获取心率数据中的患者ID
heart_rate_ids <- unique(heart_rate_data@sample_info$subject_id)

# 找到既有聚类结果又有心率数据的患者
wearable_clustering_ids <- intersect(clustering_patients, heart_rate_ids)

cat("可穿戴设备聚类分析患者总数:", length(clustering_patients), "\n")
cat("心率数据中的患者总数:", length(heart_rate_ids), "\n")
cat("既有聚类又有心率数据的患者数:", length(wearable_clustering_ids), "\n\n")

# ================== 3. 计算日心率模式（仅步数≤1） ==================

calculate_daily_rhr_pattern_wearable_cohort <- function(heart_rate_data, wearable_clustering_ids) {
  
  cat("计算可穿戴设备聚类人群的日心率模式...\n")
  
  # 过滤步数≤1的样本，且限制在聚类人群
  filtered_samples <- heart_rate_data@sample_info %>%
    filter(
      label == "<1",                           # 只要步数≤1的数据
      subject_id %in% wearable_clustering_ids  # 只要聚类人群
    )
  
  cat("过滤后的样本数:", nrow(filtered_samples), "\n")
  cat("涉及患者数:", length(unique(filtered_samples$subject_id)), "\n")
  
  # 获取对应的心率值
  heart_rates <- heart_rate_data@expression_data[1, filtered_samples$sample_id] %>%
    as.numeric()
  
  # 创建数据框并计算每小时统计
  daily_pattern <- filtered_samples %>%
    mutate(
      hour = hour(measure_time),
      heart_rate = heart_rates
    ) %>%
    filter(heart_rate >= 30 & heart_rate <= 200) %>%  # 过滤合理心率范围
    group_by(hour) %>%
    summarise(
      median_rhr = median(heart_rate, na.rm = TRUE),
      mean_rhr = mean(heart_rate, na.rm = TRUE),
      n_measurements = n(),
      sd_rhr = sd(heart_rate, na.rm = TRUE),
      q25_rhr = quantile(heart_rate, 0.25, na.rm = TRUE),
      q75_rhr = quantile(heart_rate, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("每小时数据点统计:\n")
  print(daily_pattern %>% dplyr::select(hour, n_measurements))
  
  return(daily_pattern)
}

# 计算日心率模式
daily_rhr_pattern <- calculate_daily_rhr_pattern_wearable_cohort(heart_rate_data, wearable_clustering_ids)

# 创建完整的24小时数据（填充缺失小时）
all_hours <- tibble(hour = 0:23)
daily_rhr_pattern_complete <- all_hours %>%
  left_join(daily_rhr_pattern, by = "hour")

# 计算总体中位数RHR作为参考线
overall_median_rhr <- median(daily_rhr_pattern$median_rhr, na.rm = TRUE)

cat("总体中位数RHR:", round(overall_median_rhr, 1), "bpm\n")
cat("RHR范围:", round(min(daily_rhr_pattern$median_rhr, na.rm = TRUE), 1), "-", 
    round(max(daily_rhr_pattern$median_rhr, na.rm = TRUE), 1), "bpm\n\n")

# ================== 4. 创建日心率模式图 ==================

create_daily_rhr_plot <- function(daily_pattern, overall_median_rhr) {
  
  # 创建主图
  p_daily <- ggplot(daily_pattern, aes(x = hour, y = median_rhr)) +
    # 添加线和点
    geom_line(color = "#2E8B57", size = 1.2) +
    geom_point(color = "#2E8B57", size = 3, alpha = 0.8) +
    
    # 添加中位数参考线
    geom_hline(
      yintercept = overall_median_rhr,
      color = "red",
      linetype = "dashed",
      alpha = 0.7
    ) +
    
    # 添加中位数RHR标签
    annotate(
      "text",
      x = 23,
      y = overall_median_rhr,
      label = paste("Median RHR:", round(overall_median_rhr, 1), "bpm"),
      hjust = 1,
      vjust = -0.5,
      color = "red",
      size = 3.5
    ) +
    
    # 自定义坐标轴
    scale_x_continuous(
      breaks = seq(0, 23, 3),
      limits = c(0, 23),
      name = "Hour of Day"
    ) +
    scale_y_continuous(
      name = "Median RHR (bpm)",
      limits = c(
        min(daily_pattern$median_rhr, na.rm = TRUE) - 2,
        max(daily_pattern$median_rhr, na.rm = TRUE) + 2
      )
    ) +
    
    # 标题和主题
    labs(
      title = "Daily RHR Pattern - Wearable Clustering Cohort",
      subtitle = paste("Steps ≤1 only | n =", length(wearable_clustering_ids), "patients from clustering analysis"),
      caption = "Red dashed line shows overall median RHR"
    ) +
    
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.grid.major = element_line(color = "grey95", size = 0.5),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, hjust = 1)
    )
  
  return(p_daily)
}

# 创建日心率模式图
p_daily_main <- create_daily_rhr_plot(daily_rhr_pattern_complete, overall_median_rhr)

# 显示图表
print(p_daily_main)

# 保存图表
ggsave("daily_rhr_pattern_wearable_clustering.pdf", p_daily_main, width = 10, height = 6, dpi = 300)
ggsave("daily_rhr_pattern_wearable_clustering.png", p_daily_main, width = 10, height = 6, dpi = 300)

# ================== 5. 创建带误差带的版本 ==================

create_daily_rhr_plot_with_ribbon <- function(daily_pattern, overall_median_rhr) {
  
  p_daily_ribbon <- ggplot(daily_pattern, aes(x = hour, y = median_rhr)) +
    # 添加四分位数范围带
    geom_ribbon(aes(ymin = q25_rhr, ymax = q75_rhr), 
                fill = "#2E8B57", alpha = 0.3) +
    
    # 添加中位数线
    geom_line(color = "#2E8B57", size = 1.2) +
    geom_point(color = "#2E8B57", size = 3, alpha = 0.8) +
    
    # 添加中位数参考线
    geom_hline(
      yintercept = overall_median_rhr,
      color = "red",
      linetype = "dashed",
      alpha = 0.7
    ) +
    
    # 添加标签
    annotate(
      "text",
      x = 23,
      y = overall_median_rhr,
      label = paste("Overall Median:", round(overall_median_rhr, 1), "bpm"),
      hjust = 1,
      vjust = -0.5,
      color = "red",
      size = 3.5
    ) +
    
    # 坐标轴设置
    scale_x_continuous(
      breaks = seq(0, 23, 3),
      limits = c(0, 23),
      name = "Hour of Day"
    ) +
    scale_y_continuous(
      name = "RHR (bpm)"
    ) +
    
    # 标题和主题
    labs(
      title = "Daily RHR Pattern with Quartile Range",
      subtitle = paste("Wearable clustering cohort (Steps ≤1) | n =", length(wearable_clustering_ids), "patients"),
      caption = "Shaded area shows 25th-75th percentile range | Line shows median RHR"
    ) +
    
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.grid.major = element_line(color = "grey95", size = 0.5),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, hjust = 1)
    )
  
  return(p_daily_ribbon)
}

# 创建带误差带的版本
p_daily_ribbon <- create_daily_rhr_plot_with_ribbon(daily_rhr_pattern_complete, overall_median_rhr)

# 保存带误差带的版本
ggsave("daily_rhr_pattern_with_quartiles_wearable_clustering.pdf", p_daily_ribbon, 
       width = 10, height = 6, dpi = 300)

# ================== 6. 保存数据和统计 ==================

# 保存日心率模式数据
write.csv(daily_rhr_pattern_complete, "daily_rhr_pattern_wearable_clustering.csv", row.names = FALSE)

# 生成统计摘要
daily_rhr_summary <- daily_rhr_pattern_complete %>%
  filter(!is.na(median_rhr)) %>%
  summarise(
    n_hours_with_data = n(),
    overall_median_rhr = round(median(median_rhr, na.rm = TRUE), 1),
    overall_mean_rhr = round(mean(median_rhr, na.rm = TRUE), 1),
    min_hourly_rhr = round(min(median_rhr, na.rm = TRUE), 1),
    max_hourly_rhr = round(max(median_rhr, na.rm = TRUE), 1),
    rhr_range = round(max(median_rhr, na.rm = TRUE) - min(median_rhr, na.rm = TRUE), 1),
    total_measurements = sum(n_measurements, na.rm = TRUE)
  )

cat("========== 日心率模式统计摘要 ==========\n")
print(daily_rhr_summary)

# 找出RHR最高和最低的时间
peak_hours <- daily_rhr_pattern_complete %>%
  filter(!is.na(median_rhr)) %>%
  arrange(desc(median_rhr)) %>%
  slice_head(n = 3)

trough_hours <- daily_rhr_pattern_complete %>%
  filter(!is.na(median_rhr)) %>%
  arrange(median_rhr) %>%
  slice_head(n = 3)

cat("\nRHR最高的3个小时:\n")
print(peak_hours %>% dplyr::select(hour, median_rhr, n_measurements))

cat("\nRHR最低的3个小时:\n") 
print(trough_hours %>% dplyr::select(hour, median_rhr, n_measurements))

# 保存统计摘要
write.csv(daily_rhr_summary, "daily_rhr_summary_wearable_clustering.csv", row.names = FALSE)

# ================== 7. 生成报告 ==================

generate_daily_rhr_report <- function() {
  
  report <- paste0(
    "========================================\n",
    "可穿戴设备聚类人群日心率模式分析报告\n",
    "========================================\n\n",
    
    "📊 分析概述:\n",
    "本报告展示了纳入可穿戴设备聚类分析患者的24小时心率变化模式。\n",
    "仅使用步数≤1的静息心率数据，确保分析的是真正的静息心率。\n\n",
    
    "👥 研究队列:\n",
    "- 可穿戴设备聚类分析患者总数: ", length(clustering_patients), "\n",
    "- 有心率数据的聚类患者数: ", length(wearable_clustering_ids), "\n",
    "- 数据筛选标准: 步数≤1（静息状态）\n\n",
    
    "📈 日心率模式特征:\n",
    "- 整体中位数RHR: ", daily_rhr_summary$overall_median_rhr, " bpm\n",
    "- 整体平均RHR: ", daily_rhr_summary$overall_mean_rhr, " bpm\n",
    "- 日内RHR变化范围: ", daily_rhr_summary$rhr_range, " bpm\n",
    "- 最低RHR: ", daily_rhr_summary$min_hourly_rhr, " bpm\n",
    "- 最高RHR: ", daily_rhr_summary$max_hourly_rhr, " bpm\n",
    "- 总测量数据点: ", format(daily_rhr_summary$total_measurements, big.mark = ","), "\n\n"
  )
  
  # 添加峰值和谷值时间
  if(exists("peak_hours") && exists("trough_hours")) {
    report <- paste0(report,
                     "⏰ 关键时间点:\n",
                     "RHR最高时间: ", paste(peak_hours$hour[1:3], "时", collapse = ", "), "\n",
                     "RHR最低时间: ", paste(trough_hours$hour[1:3], "时", collapse = ", "), "\n\n")
  }
  
  report <- paste0(report,
                   "📊 生成的可视化:\n",
                   "1. daily_rhr_pattern_wearable_clustering.pdf\n",
                   "   - 基础日心率模式图\n",
                   "2. daily_rhr_pattern_wearable_clustering.png\n",
                   "   - PNG格式日心率模式图\n",
                   "3. daily_rhr_pattern_with_quartiles_wearable_clustering.pdf\n",
                   "   - 带四分位数范围的日心率模式图\n\n",
                   
                   "📄 数据文件:\n",
                   "1. daily_rhr_pattern_wearable_clustering.csv\n",
                   "   - 24小时心率模式详细数据\n",
                   "2. daily_rhr_summary_wearable_clustering.csv\n",
                   "   - 日心率模式统计摘要\n\n",
                   
                   "🎨 可视化特点:\n",
                   "✅ 专门针对可穿戴设备聚类人群\n",
                   "✅ 仅使用静息心率数据（步数≤1）\n",
                   "✅ 24小时完整时间覆盖\n",
                   "✅ 中位数参考线便于比较\n",
                   "✅ 四分位数范围显示变异性\n\n",
                   
                   "💡 临床意义:\n",
                   "- 展示聚类人群的生理节律特征\n",
                   "- 为个性化监护提供基线参考\n",
                   "- 识别异常心率模式的时间窗口\n",
                   "- 支持基于时间的风险分层\n\n",
                   
                   "🔍 模式解读:\n",
                   "- 日内RHR变化反映自主神经活动节律\n",
                   "- 夜间通常为RHR最低时期\n",
                   "- 白天活动期RHR相对较高\n",
                   "- 变异性大小反映个体差异\n\n",
                   
                   "📝 研究价值:\n",
                   "- 为可穿戴设备聚类结果提供生理基础\n",
                   "- 验证聚类人群的心率节律特征\n",
                   "- 为时间相关的健康监测提供依据\n",
                   "- 支持个性化医疗决策制定\n\n",
                   
                   "报告生成时间: ", Sys.time(), "\n",
                   "========================================"
  )
  
  writeLines(report, "daily_rhr_pattern_report.txt")
  cat("✓ 详细报告已保存: daily_rhr_pattern_report.txt\n")
  
  return(report)
}

# 生成报告
daily_report <- generate_daily_rhr_report()

cat("\n============ 日心率模式分析完成 ============\n")
cat("📊 生成的文件:\n")
cat("✅ daily_rhr_pattern_wearable_clustering.pdf - 基础日心率图\n")
cat("✅ daily_rhr_pattern_wearable_clustering.png - PNG格式图\n") 
cat("✅ daily_rhr_pattern_with_quartiles_wearable_clustering.pdf - 带四分位数范围图\n")
cat("✅ daily_rhr_pattern_wearable_clustering.csv - 详细数据\n")
cat("✅ daily_rhr_summary_wearable_clustering.csv - 统计摘要\n")
cat("✅ daily_rhr_pattern_report.txt - 分析报告\n")

cat("\n🎯 关键发现:\n")
cat("- 聚类人群数:", length(wearable_clustering_ids), "患者\n")
cat("- 日内RHR变化:", daily_rhr_summary$rhr_range, "bpm\n")
cat("- 整体中位数RHR:", daily_rhr_summary$overall_median_rhr, "bpm\n")

cat("\n🎉 可穿戴设备聚类人群日心率模式分析完成！\n")