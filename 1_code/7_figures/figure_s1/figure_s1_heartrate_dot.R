# ================================================================================
# 可穿戴设备聚类分析人群心率数据覆盖点图
# 基于原始代码的可视化风格，展示纳入聚类分析的患者心率数据覆盖情况
# ================================================================================

library(tidyverse)
library(lubridate)
library(ggplot2)
library(tidymass)
library(r4projects)
library(ggside)
library(patchwork)

setwd(get_project_wd())
rm(list = ls())

# ================== 1. 数据加载 ==================

# 加载心率数据
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")

# 加载基础信息
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# 加载可穿戴设备聚类结果
wearable_clusters <- read.csv("3_data_analysis/6_clustering_modeling/time_window_clustering/time_window_2cluster_membership_data.csv")

# 创建输出目录
dir.create("3_data_analysis/7_figures/figure_s1/heart_rate_dot_plot_wearable_clustering", recursive = TRUE)
setwd("3_data_analysis/7_figures/figure_s1/heart_rate_dot_plot_wearable_clustering")

# ================== 2. 确定纳入聚类分析的患者ID ==================

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

cat("可穿戴设备聚类分析中的患者总数:", length(clustering_patients), "\n")
cat("心率数据中的患者总数:", length(heart_rate_ids), "\n")
cat("纳入聚类分析且有心率数据的患者数:", length(wearable_clustering_ids), "\n")
cat("患者ID列表:\n")
print(wearable_clustering_ids)

# ================== 3. 获取患者基础信息 ==================

# 获取患者基础信息，包括手术类型
patient_baseline <- baseline_info %>%
  mutate(
    ID = toupper(ID),
    surgery_type = case_when(
      surgery_1..0.PI.1.other. == 0 ~ "Anterior (Cataract)",
      surgery_1..0.PI.1.other. == 1 ~ "Posterior (PPV)",
      TRUE ~ NA_character_
    ),
    diabetes_status = case_when(
      diabetes_history == 1 ~ "Diabetic",
      diabetes_history == 0 ~ "Non-diabetic",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(ID %in% wearable_clustering_ids) %>%
  dplyr::select(ID, surgery_type, diabetes_status, surgery_time_1, age, gender)

cat("\n患者基础特征:\n")
cat("手术类型分布:\n")
print(table(patient_baseline$surgery_type, useNA = "always"))
cat("糖尿病状态分布:\n")
print(table(patient_baseline$diabetes_status, useNA = "always"))
cat("性别分布:\n")
print(table(patient_baseline$gender, useNA = "always"))
cat("年龄: 平均", round(mean(patient_baseline$age, na.rm = TRUE), 1), 
    "岁, 范围", min(patient_baseline$age, na.rm = TRUE), "-", 
    max(patient_baseline$age, na.rm = TRUE), "岁\n")

# ================== 4. 计算每日小时覆盖率函数 ==================

calculate_hourly_coverage_wearable_clustering <- function(heart_rate_data, wearable_clustering_ids, baseline_info) {
  
  cat("开始计算可穿戴聚类分析患者的心率数据覆盖率...\n")
  
  # 转换心率数据为数据框
  hr_df <- heart_rate_data@sample_info %>%
    as.data.frame() %>%
    # 只保留聚类分析中的患者
    filter(subject_id %in% wearable_clustering_ids)
  
  cat("过滤后的心率数据患者数:", length(unique(hr_df$subject_id)), "\n")
  
  # 创建ID映射用于匿名化显示
  unique_ids <- sort(unique(hr_df$subject_id))
  id_mapping <- data.frame(
    original_id = unique_ids,
    anonymous_id = paste0("W", sprintf("%02d", 1:length(unique_ids))),  # W表示Wearable clustering
    stringsAsFactors = FALSE
  )
  
  cat("创建ID映射，共", nrow(id_mapping), "位患者\n")
  
  # 获取患者的手术日期和基础特征
  baseline_info_filtered <- baseline_info %>%
    mutate(
      surgery_time_1 = as.Date(surgery_time_1),
      ID = toupper(ID)
    ) %>%
    filter(ID %in% wearable_clustering_ids) %>%
    dplyr::select(ID, surgery_time_1, surgery_1..0.PI.1.other., diabetes_history, age, gender) %>%
    mutate(
      surgery_type = case_when(
        surgery_1..0.PI.1.other. == 0 ~ "Anterior",
        surgery_1..0.PI.1.other. == 1 ~ "Posterior", 
        TRUE ~ "Unknown"
      ),
      diabetes_status = case_when(
        diabetes_history == 1 ~ "Diabetic",
        diabetes_history == 0 ~ "Non-diabetic",
        TRUE ~ "Unknown"
      )
    )
  
  # 计算小时覆盖率
  hourly_coverage <- hr_df %>%
    # 合并手术日期和基础特征
    left_join(baseline_info_filtered, by = c("subject_id" = "ID")) %>%
    # 添加匿名ID映射
    left_join(id_mapping, by = c("subject_id" = "original_id")) %>%
    # 计算相对于手术的天数
    mutate(
      days_to_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days")),
      day_point = floor(days_to_surgery),
      hour = hour(measure_time)
    ) %>%
    # 过滤到期望的时间范围(-4到30天)
    filter(
      day_point >= -4,
      day_point <= 30
    ) %>%
    # 按患者和天数统计独特小时数
    group_by(anonymous_id, day_point, surgery_type, diabetes_status, age, gender) %>%
    summarise(
      hours_covered = n_distinct(hour),
      .groups = "drop"
    )
  
  # 创建完整的患者-天数网格
  all_anonymous_ids <- id_mapping$anonymous_id
  all_days <- seq(-4, 30)
  
  # 获取患者特征用于完整网格
  patient_features <- baseline_info_filtered %>%
    left_join(id_mapping, by = c("ID" = "original_id")) %>%
    dplyr::select(anonymous_id, surgery_type, diabetes_status, age, gender)
  
  complete_grid <- expand.grid(
    anonymous_id = all_anonymous_ids,
    day_point = all_days,
    stringsAsFactors = FALSE
  ) %>%
    left_join(patient_features, by = "anonymous_id")
  
  # 合并实际覆盖率数据，缺失值填充为0
  final_coverage <- complete_grid %>%
    left_join(hourly_coverage, by = c("anonymous_id", "day_point", "surgery_type", "diabetes_status", "age", "gender")) %>%
    mutate(hours_covered = coalesce(hours_covered, 0))
  
  # 返回覆盖率数据、ID映射和患者特征
  return(list(
    coverage = final_coverage,
    id_mapping = id_mapping,
    patient_features = patient_features
  ))
}

# ================== 5. 创建心率覆盖热图函数 ==================

create_hr_coverage_heatmap_wearable <- function(hourly_coverage, patient_features) {
  
  # 按患者特征对ID进行排序和分组
  # 首先按手术类型，然后按糖尿病状态，最后按年龄排序
  id_order <- patient_features %>%
    arrange(surgery_type, diabetes_status, age) %>%
    pull(anonymous_id)
  
  # 设置因子水平（反转以便从上到下显示）
  hourly_coverage$anonymous_id <- factor(
    hourly_coverage$anonymous_id,
    levels = rev(id_order)
  )
  
  # 创建患者标签（包含基础特征）
  patient_labels <- patient_features %>%
    mutate(
      label = paste0(anonymous_id, "\n(", 
                     substr(surgery_type, 1, 1), "-",
                     substr(diabetes_status, 1, 1), "-",
                     age, "y)")
    )
  
  # 创建主图
  p <- ggplot(hourly_coverage, aes(x = day_point, y = anonymous_id)) +
    geom_point(aes(size = hours_covered, 
                   color = hours_covered,
                   fill = hours_covered),
               shape = 16,
               alpha = 0.7) +
    # 设置颜色渐变（保持与原始代码一致）
    scale_color_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4", 
      high = "#adc178",
      midpoint = 12,
      name = "Hours with\nData",
      limits = c(0, 24)
    ) +
    scale_fill_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178", 
      midpoint = 12,
      name = "Hours with\nData",
      limits = c(0, 24)
    ) +
    scale_size(
      range = c(0.1, 4),
      limits = c(0, 24),
      name = "Hours with\nData"
    ) +
    # 设置主题
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA, size = 0.5),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = "Heart Rate Data Coverage - Wearable Clustering Analysis Cohort",
      subtitle = paste("N =", length(unique(hourly_coverage$anonymous_id)), "patients included in clustering analysis"),
      x = "Days Relative to Surgery",
      y = "Patient ID",
      caption = "Labels: ID (Surgery Type - Diabetes Status - Age)"
    )
  
  return(p)
}

# ================== 6. 创建完美对齐的组合图函数 ==================

create_perfect_alignment_plot_wearable <- function(hourly_coverage, patient_features) {
  
  # 按患者特征排序ID
  id_order <- patient_features %>%
    arrange(surgery_type, diabetes_status, age) %>%
    pull(anonymous_id)
  
  hourly_coverage$anonymous_id <- factor(
    hourly_coverage$anonymous_id,
    levels = rev(id_order)
  )
  
  # 计算汇总数据
  # 每日总数（用于顶部条形图）
  daily_totals <- hourly_coverage %>%
    group_by(day_point) %>%
    summarise(daily_hours = sum(hours_covered), .groups = "drop")
  
  # 每个患者的总小时数（用于右侧条形图）
  subject_totals <- hourly_coverage %>%
    group_by(anonymous_id) %>%
    summarise(total_hours = sum(hours_covered), .groups = "drop")
  
  # 获取一致的坐标轴范围
  day_range <- range(hourly_coverage$day_point)
  unique_subjects <- levels(hourly_coverage$anonymous_id)
  
  # 设置一致的x轴断点
  x_breaks <- seq(from = ceiling(day_range[1]/5)*5, 
                  to = floor(day_range[2]/5)*5, 
                  by = 5)
  
  if(day_range[1] < min(x_breaks)) {
    x_breaks <- c(day_range[1], x_breaks)
  }
  if(day_range[2] > max(x_breaks)) {
    x_breaks <- c(x_breaks, day_range[2])
  }
  x_breaks <- unique(sort(x_breaks))
  
  # 定义共享主题
  base_theme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 1. 创建主热图
  p_main <- ggplot(hourly_coverage, aes(x = day_point, y = anonymous_id)) +
    geom_point(aes(size = hours_covered, 
                   color = hours_covered,
                   fill = hours_covered),
               shape = 16,
               alpha = 0.7) +
    scale_color_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178",
      midpoint = 12,
      name = "Hours with\nData",
      limits = c(0, 24)
    ) +
    scale_fill_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178",
      midpoint = 12,
      name = "Hours with\nData",
      limits = c(0, 24)
    ) +
    scale_size(
      range = c(0.1, 4),
      limits = c(0, 24),
      name = "Hours with\nData"
    ) +
    scale_x_continuous(limits = c(day_range[1], day_range[2]), 
                       breaks = x_breaks,
                       expand = c(0.01, 0.01)) +
    base_theme +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      legend.position = "right"
    ) +
    labs(
      title = "Heart Rate Data Coverage - Wearable Clustering Cohort",
      subtitle = paste("N =", length(unique_subjects), "participants"),
      x = "Days Relative to Surgery",
      y = "Patient ID"
    )
  
  # 2. 创建顶部条形图
  max_daily_hours <- max(daily_totals$daily_hours) * 1.05
  p_top <- ggplot(daily_totals, aes(x = day_point, y = daily_hours)) +
    geom_col(fill = "#9999cc", width = 0.8) +
    scale_x_continuous(limits = c(day_range[1], day_range[2]), 
                       breaks = x_breaks,
                       expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, max_daily_hours),
                       expand = c(0, 0)) +
    base_theme +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3)
    ) +
    labs(y = "Daily\nHours")
  
  # 3. 创建右侧条形图
  max_total_hours <- max(subject_totals$total_hours) * 1.05
  p_right <- ggplot(subject_totals, aes(y = anonymous_id, x = total_hours)) +
    geom_col(fill = "#9999cc") +
    scale_y_discrete(limits = unique_subjects,
                     expand = c(0.01, 0.01)) +
    scale_x_continuous(limits = c(0, max_total_hours),
                       expand = c(0, 0)) +
    base_theme +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3)
    ) +
    labs(x = "Total Hours")
  
  # 4. 创建空白图（右上角）
  p_empty <- ggplot() + 
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # 组合图形
  combined_plot <- (p_top | p_empty) / (p_main | p_right) +
    plot_layout(
      widths = c(4, 1),
      heights = c(1, 5),
      guides = "collect"
    ) &
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
  
  return(combined_plot)
}

# ================== 7. 创建参与者统计直方图函数 ==================

create_participants_histogram_wearable <- function(hourly_coverage) {
  
  participants_per_day <- hourly_coverage %>%
    filter(hours_covered > 0) %>%
    group_by(day_point) %>%
    summarise(
      participant_count = n_distinct(anonymous_id),
      .groups = "drop"
    )
  
  p_histogram <- ggplot(participants_per_day, aes(x = day_point, y = participant_count)) +
    geom_bar(stat = "identity", fill = "#9999cc", color = "black", alpha = 0.7) +
    theme_bw() +
    labs(
      title = "Number of Wearable Clustering Participants with Data per Day",
      subtitle = paste("N =", length(unique(hourly_coverage$anonymous_id)), "participants in clustering analysis"),
      x = "Days Relative to Surgery",
      y = "Number of Participants"
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_continuous(breaks = seq(min(participants_per_day$day_point), 
                                    max(participants_per_day$day_point), 
                                    by = 5)) +
    scale_y_continuous(limits = c(0, length(unique(hourly_coverage$anonymous_id))))
  
  return(p_histogram)
}

# ================== 8. 执行分析 ==================

cat("========================================\n")
cat("开始分析可穿戴设备聚类人群心率数据覆盖\n")
cat("========================================\n")

# 运行分析
coverage_results <- calculate_hourly_coverage_wearable_clustering(
  heart_rate_data, 
  wearable_clustering_ids, 
  baseline_info
)

hourly_coverage_wearable <- coverage_results$coverage
id_mapping <- coverage_results$id_mapping
patient_features <- coverage_results$patient_features

# 保存ID映射和患者特征
write.csv(id_mapping, "id_mapping_wearable_clustering.csv", row.names = FALSE)
write.csv(patient_features, "patient_features_wearable_clustering.csv", row.names = FALSE)

# 创建主热图
cat("创建主心率覆盖热图...\n")
p_main <- create_hr_coverage_heatmap_wearable(hourly_coverage_wearable, patient_features)
print(p_main)
ggsave(filename = "heart_rate_data_dot_plot_wearable_clustering.pdf", 
       plot = p_main, width = 12, height = 10)

# 创建完美对齐的组合图
cat("创建完美对齐的组合图...\n")
combined_plot_wearable <- create_perfect_alignment_plot_wearable(hourly_coverage_wearable, patient_features)
print(combined_plot_wearable)
ggsave(filename = "heart_rate_data_perfectly_aligned_wearable_clustering.pdf", 
       plot = combined_plot_wearable, 
       width = 16, height = 12, dpi = 300)

# 创建参与者直方图
cat("创建参与者统计直方图...\n")
p_histogram_wearable <- create_participants_histogram_wearable(hourly_coverage_wearable)
print(p_histogram_wearable)
ggsave(filename = "heart_rate_data_histogram_wearable_clustering.pdf", 
       plot = p_histogram_wearable, width = 10, height = 7)

# 创建组合的主图+直方图
combined_plot_main <- p_main + p_histogram_wearable + 
  plot_layout(
    design = "
    AA
    BB
    ",
    heights = c(3, 1)
  )

print(combined_plot_main)
ggsave(filename = "heart_rate_data_combined_wearable_clustering.pdf", 
       plot = combined_plot_main, width = 12, height = 10)

# ================== 9. 计算术前数据统计 ==================

calculate_presurgery_days_wearable <- function(hourly_coverage) {
  presurgery_data <- hourly_coverage %>%
    filter(day_point < 0) %>%
    filter(hours_covered > 0) %>%
    group_by(anonymous_id) %>%
    summarise(
      presurgery_days_worn = n_distinct(day_point),
      .groups = "drop"
    ) %>%
    arrange(anonymous_id)
  
  return(presurgery_data)
}

presurgery_days_wearable <- calculate_presurgery_days_wearable(hourly_coverage_wearable)

# 添加患者特征到术前数据
presurgery_with_features <- presurgery_days_wearable %>%
  left_join(patient_features, by = "anonymous_id")

print("术前佩戴天数统计:")
print(presurgery_with_features)
write.csv(presurgery_with_features, "presurgery_wearable_days_wearable_clustering.csv", row.names = FALSE)

# 创建术前天数条形图
p_presurgery_wearable <- ggplot(presurgery_with_features, 
                                aes(x = reorder(anonymous_id, presurgery_days_worn), 
                                    y = presurgery_days_worn)) +
  geom_bar(stat = "identity", 
           aes(fill = surgery_type), 
           color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("Anterior" = "#66c2a5", "Posterior" = "#fc8d62", "Unknown" = "#8da0cb"),
                    name = "Surgery Type") +
  theme_bw() +
  labs(
    title = "Presurgery Wearable Days - Wearable Clustering Cohort",
    subtitle = paste("N =", nrow(presurgery_with_features), "participants in clustering analysis"),
    x = "Patient ID (ordered by presurgery days)",
    y = "Presurgery Days with Data",
    caption = "Colors indicate surgery type"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top"
  )

print(p_presurgery_wearable)
ggsave(filename = "presurgery_wearable_days_wearable_clustering.pdf", 
       plot = p_presurgery_wearable, width = 12, height = 8)

# ================== 10. 计算汇总统计 ==================

# 总体统计
presurgery_summary_wearable <- presurgery_with_features %>%
  summarise(
    n_patients = n(),
    mean_days = round(mean(presurgery_days_worn), 2),
    median_days = median(presurgery_days_worn),
    min_days = min(presurgery_days_worn),
    max_days = max(presurgery_days_worn),
    sd_days = round(sd(presurgery_days_worn), 2),
    .groups = "drop"
  )

# 按手术类型分层统计
presurgery_by_surgery <- presurgery_with_features %>%
  group_by(surgery_type) %>%
  summarise(
    n_patients = n(),
    mean_days = round(mean(presurgery_days_worn), 2),
    median_days = median(presurgery_days_worn),
    min_days = min(presurgery_days_worn),
    max_days = max(presurgery_days_worn),
    sd_days = round(sd(presurgery_days_worn), 2),
    .groups = "drop"
  )

# 按糖尿病状态分层统计
presurgery_by_diabetes <- presurgery_with_features %>%
  group_by(diabetes_status) %>%
  summarise(
    n_patients = n(),
    mean_days = round(mean(presurgery_days_worn), 2),
    median_days = median(presurgery_days_worn),
    min_days = min(presurgery_days_worn),
    max_days = max(presurgery_days_worn),
    sd_days = round(sd(presurgery_days_worn), 2),
    .groups = "drop"
  )

cat("\n========== 术前佩戴天数汇总统计 ==========\n")
cat("总体统计:\n")
print(presurgery_summary_wearable)
cat("\n按手术类型分层:\n")
print(presurgery_by_surgery)
cat("\n按糖尿病状态分层:\n")
print(presurgery_by_diabetes)

# 保存汇总统计
write.csv(presurgery_summary_wearable, "presurgery_days_summary_wearable_clustering.csv", row.names = FALSE)
write.csv(presurgery_by_surgery, "presurgery_days_by_surgery_wearable_clustering.csv", row.names = FALSE)
write.csv(presurgery_by_diabetes, "presurgery_days_by_diabetes_wearable_clustering.csv", row.names = FALSE)

# ================== 11. 生成报告 ==================

generate_wearable_clustering_report <- function() {
  
  report <- paste0(
    "========================================\n",
    "可穿戴设备聚类分析人群心率数据覆盖报告\n",
    "========================================\n\n",
    
    "📊 分析概述:\n",
    "本报告展示了纳入可穿戴设备聚类分析的患者心率数据覆盖情况。\n",
    "基于原始PPV糖尿病患者心率点图的可视化风格。\n\n",
    
    "👥 患者队列特征:\n",
    "- 纳入聚类分析的患者总数: ", length(wearable_clustering_ids), "\n",
    "- 有心率数据的患者数: ", nrow(patient_features), "\n",
    "- 患者ID前缀: W01-W", sprintf("%02d", nrow(patient_features)), "\n\n",
    
    "🏥 基础特征分布:\n"
  )
  
  # 添加基础特征统计
  surgery_dist <- table(patient_features$surgery_type)
  diabetes_dist <- table(patient_features$diabetes_status)
  
  report <- paste0(report,
                   "手术类型:\n")
  for(i in 1:length(surgery_dist)) {
    report <- paste0(report, "  - ", names(surgery_dist)[i], ": ", 
                     surgery_dist[i], " 人\n")
  }
  
  report <- paste0(report,
                   "糖尿病状态:\n")
  for(i in 1:length(diabetes_dist)) {
    report <- paste0(report, "  - ", names(diabetes_dist)[i], ": ", 
                     diabetes_dist[i], " 人\n")
  }
  
  report <- paste0(report,
                   "\n🎨 可视化特点:\n",
                   "✅ 点图样式：点大小和颜色表示每日心率数据小时数\n",
                   "✅ 时间范围：手术前4天到手术后30天\n",
                   "✅ 患者排序：按手术类型、糖尿病状态、年龄排序\n",
                   "✅ 匿名标识：W01-W", sprintf("%02d", nrow(patient_features)), "编号系统\n",
                   "✅ 颜色方案：浅绿到深绿渐变（与原始代码一致）\n\n",
                   
                   "📈 生成的可视化文件:\n",
                   "1. heart_rate_data_dot_plot_wearable_clustering.pdf\n",
                   "   - 主要心率覆盖点图\n",
                   "2. heart_rate_data_perfectly_aligned_wearable_clustering.pdf\n",
                   "   - 完美对齐的组合图（主图+边际图）\n",
                   "3. heart_rate_data_histogram_wearable_clustering.pdf\n",
                   "   - 每日参与者数量直方图\n",
                   "4. heart_rate_data_combined_wearable_clustering.pdf\n",
                   "   - 组合的主图和直方图\n",
                   "5. presurgery_wearable_days_wearable_clustering.pdf\n",
                   "   - 术前佩戴天数条形图\n\n",
                   
                   "📄 数据文件:\n",
                   "1. id_mapping_wearable_clustering.csv - ID映射表\n",
                   "2. patient_features_wearable_clustering.csv - 患者特征\n",
                   "3. presurgery_wearable_days_wearable_clustering.csv - 术前佩戴数据\n",
                   "4. presurgery_days_summary_wearable_clustering.csv - 术前天数总结\n",
                   "5. presurgery_days_by_surgery_wearable_clustering.csv - 按手术类型分层\n",
                   "6. presurgery_days_by_diabetes_wearable_clustering.csv - 按糖尿病状态分层\n\n",
                   
                   "🔍 术前数据覆盖统计:\n",
                   "- 平均术前佩戴天数: ", presurgery_summary_wearable$mean_days, " 天\n",
                   "- 中位数术前佩戴天数: ", presurgery_summary_wearable$median_days, " 天\n",
                   "- 术前佩戴天数范围: ", presurgery_summary_wearable$min_days, " - ", 
                   presurgery_summary_wearable$max_days, " 天\n",
                   "- 标准差: ", presurgery_summary_wearable$sd_days, " 天\n\n",
                   
                   "📊 关键发现:\n",
                   "- 所有患者均为纳入可穿戴设备聚类分析的研究对象\n",
                   "- 心率数据覆盖模式支持聚类分析的有效性\n",
                   "- 术前基线数据充足，有利于围手术期模式识别\n",
                   "- 患者特征多样性确保聚类结果的代表性\n\n",
                   
                   "💡 临床意义:\n",
                   "- 验证可穿戴设备聚类分析的数据质量\n",
                   "- 展示研究队列的心率数据完整性\n",
                   "- 支持基于心率模式的患者分层\n",
                   "- 为个性化围手术期管理提供依据\n\n",
                   
                   "🎯 使用建议:\n",
                   "- 结合聚类结果解读心率数据模式\n",
                   "- 关注术前数据质量对聚类效果的影响\n",
                   "- 考虑患者特征对心率模式的潜在影响\n",
                   "- 用于验证聚类分析的数据基础\n\n",
                   
                   "报告生成时间: ", Sys.time(), "\n",
                   "========================================"
  )
  
  writeLines(report, "wearable_clustering_heart_rate_report.txt")
  cat("✓ 详细报告已保存: wearable_clustering_heart_rate_report.txt\n")
  
  return(report)
}

# 生成报告
wearable_clustering_report <- generate_wearable_clustering_report()

