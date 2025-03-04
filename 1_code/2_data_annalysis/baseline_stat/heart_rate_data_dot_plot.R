library(tidyverse)
library(lubridate)
library(ggplot2)
library(tidymass)
library(r4projects)
library(ggside)

setwd(get_project_wd())
rm(list = ls())


load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")


#############
dir.create("3_data_analysis/2_data_analysis/baseline_stat/heart_rate_dot_plot", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/baseline_stat/heart_rate_dot_plot")

# Function to calculate daily hour coverage
calculate_hourly_coverage <- function(heart_rate_data) {
  # Convert heart_rate_data from mass_dataset to data frame for processing
  hr_df <- heart_rate_data@sample_info %>%
    as.data.frame()
  
  # Get surgery dates for each subject
  baseline_info <- baseline_info %>%
    mutate(
      surgery_time_1 = as.Date(surgery_time_1),
      ID = toupper(ID)
    ) %>%
    dplyr::select(ID, surgery_time_1)
  
  # Calculate hours coverage
  hourly_coverage <- hr_df %>%
    # Join with surgery dates
    left_join(baseline_info, by = c("subject_id" = "ID")) %>%
    # Calculate days relative to surgery
    mutate(
      days_to_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days")),
      day_point = floor(days_to_surgery),
      hour = hour(measure_time)
    ) %>%
    # Filter to desired range (-7 to 30 days)
    filter(
      day_point >= -7,
      day_point <= 30
    ) %>%
    # Count unique hours per subject per day
    group_by(subject_id, day_point) %>%
    summarise(
      hours_covered = n_distinct(hour),
      .groups = "drop"
    )
  
  # Create complete grid with all subject-day combinations
  all_subjects <- unique(hourly_coverage$subject_id)
  all_days <- seq(-7, 30)
  complete_grid <- expand.grid(
    subject_id = all_subjects,
    day_point = all_days
  )
  
  # Join with actual coverage and fill missing with 0
  final_coverage <- complete_grid %>%
    left_join(hourly_coverage, by = c("subject_id", "day_point")) %>%
    mutate(hours_covered = coalesce(hours_covered, 0))
  
  return(final_coverage)
}

########
create_hr_coverage_heatmap <- function(hourly_coverage) {
  # Sort subject IDs numerically by extracting the number and sorting
  hourly_coverage$subject_id <- factor(
    hourly_coverage$subject_id,
    levels = rev(unique(hourly_coverage$subject_id)[order(as.numeric(gsub("SH0*", "", unique(hourly_coverage$subject_id))))])
  )
  
  # Create the plot
  p <- ggplot(hourly_coverage, aes(x = day_point, y = subject_id)) +
    geom_point(aes(size = hours_covered, 
                   color = hours_covered,
                   fill = hours_covered),
               shape = 16,   # 实心圆点，没有边框
               alpha = 0.7) +
    # 分别设置color和fill的渐变
    scale_color_gradient2(
      low = "#FFF8D5",
      mid = "#EED18A",
      high = "#FAAF00",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_fill_gradient2(
      low = "#FFF8D5",
      mid = "#EED18A",
      high = "#FAAF00",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_size(
      range = c(0.1, 4),
      limits = c(0, 24),
      name = "Hours with\nData"
    ) +
    scale_alpha(
      range = c(0.3, 1),
      guide = "none"
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
      title = "Heart Rate Data Coverage Pattern",
      x = "Days Relative to Surgery",
      y = "Subject ID"
    )
  
  return(p)
}

# Run analysis
hourly_coverage <- calculate_hourly_coverage(heart_rate_data)
p <- create_hr_coverage_heatmap(hourly_coverage)
print(p)

ggsave(filename = "heart_rate_data_dot_plot.pdf",plot = p,width=8,height = 7)



# Create histogram of participants with data per day
participants_per_day <- hourly_coverage %>%
  filter(hours_covered > 0) %>%
  group_by(day_point) %>%
  summarise(
    participant_count = n_distinct(subject_id)
  )

# Create the histogram plot
p_histogram <- ggplot(participants_per_day, aes(x = day_point, y = participant_count)) +
  geom_bar(stat = "identity", fill = "#FAAF00", color = "black", alpha = 0.7) +
  theme_bw() +
  labs(
    title = "Number of Participants with Wearable Data per Day",
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
                                  by = 5))

# Print the plot
print(p_histogram)
ggsave(filename = "heart_rate_data_histogram.pdf",plot = p_histogram,width=8,height = 7)

combined_plot <- p + p_histogram + 
  plot_layout(
    design = "
    AA
    BB
    ",
    heights = c(3, 1)  # Dot plot takes up more vertical space
  )
combined_plot
ggsave(filename = "heart_rate_data_combined.pdf",plot = combined_plot,width=8,height = 7)



# 在已有代码基础上添加计算术前佩戴天数的函数
calculate_presurgery_days <- function(hourly_coverage) {
  # 筛选术前数据（day_point < 0）
  presurgery_data <- hourly_coverage %>%
    filter(day_point < 0) %>%
    # 只计算有数据的天数（hours_covered > 0）
    filter(hours_covered > 0) %>%
    # 按受试者分组
    group_by(subject_id) %>%
    # 计算每个受试者术前有数据的天数
    summarise(
      presurgery_days_worn = n_distinct(day_point),
      .groups = "drop"
    ) %>%
    # 按ID排序
    arrange(subject_id)
  
  return(presurgery_data)
}

# 使用已计算的hourly_coverage数据
presurgery_days <- calculate_presurgery_days(hourly_coverage)

# 打印结果
print(presurgery_days, n=49)

# 保存结果到CSV文件
write.csv(presurgery_days, "presurgery_wearable_days.csv", row.names = FALSE)

# 创建术前佩戴天数的柱状图
p_presurgery <- ggplot(presurgery_days, aes(x = subject_id, y = presurgery_days_worn)) +
  geom_bar(stat = "identity", fill = "#FAAF00", color = "black", alpha = 0.7) +
  theme_bw() +
  labs(
    title = "presurgery_wearable_days",
    x = "Subejct_ID",
    y = "presurgery_days"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

# 打印图表
print(p_presurgery)
ggsave(filename = "presurgery_wearable_days.pdf", plot = p_presurgery, width = 10, height = 6)

# 计算术前佩戴天数的汇总统计
presurgery_summary <- presurgery_days %>%
  summarise(
    mean_days = mean(presurgery_days_worn),
    median_days = median(presurgery_days_worn),
    min_days = min(presurgery_days_worn),
    max_days = max(presurgery_days_worn),
    sd_days = sd(presurgery_days_worn)
  )

print(presurgery_summary)
write.csv(presurgery_summary, "presurgery_days_summary.csv", row.names = FALSE)







# 基于已有的hourly_coverage数据，筛选高质量参与者

# 定义函数来评估参与者数据质量
evaluate_participant_quality <- function(hourly_coverage) {
  
  # 计算每个参与者的数据质量指标
  participant_quality <- hourly_coverage %>%
    group_by(subject_id) %>%
    summarise(
      # 1. 总体覆盖天数（有任何数据的天数）
      total_days_with_data = sum(hours_covered > 0),
      
      # 2. 术前覆盖天数
      presurgery_days = sum(day_point < 0 & hours_covered > 0),
      
      # 3. 术后覆盖天数
      postsurgery_days = sum(day_point >= 0 & hours_covered > 0),
      
      # 4. 高质量天数（每天至少12小时数据）
      high_quality_days = sum(hours_covered >= 12),
      
      # 5. 术前高质量天数
      presurgery_high_quality = sum(day_point < 0 & hours_covered >= 12),
      
      # 6. 术后高质量天数
      postsurgery_high_quality = sum(day_point >= 0 & hours_covered >= 12),
      
      # 7. 平均每天覆盖小时数
      mean_hours_per_day = mean(hours_covered),
      
      # 8. 连续性指标 - 最长连续监测天数（有一定数据，如至少6小时）
      max_consecutive_days = {
        # 找出有数据的天数（至少6小时）
        days_with_data <- day_point[hours_covered >= 6]
        calculate_max_consecutive_days(days_with_data)
      },
      
      # 9. 整体数据完整率（总记录小时数/理论最大小时数）
      data_completeness = sum(hours_covered) / (n() * 24),
      
      .groups = "drop"
    )
  
  return(participant_quality)
}

# 辅助函数：计算最长连续监测天数
calculate_max_consecutive_days <- function(days_data) {
  # 确保天数已排序
  days_data <- sort(days_data)
  
  if (length(days_data) == 0) {
    return(0)
  }
  
  # 计算连续天数
  gaps <- diff(days_data)
  # 找出所有间隔为1的位置（连续的天数）
  breaks <- which(gaps != 1)
  
  # 计算每段连续天数的长度
  if (length(breaks) == 0) {
    # 如果没有断点，则所有天数都是连续的
    max_consecutive <- length(days_data)
  } else {
    # 添加起始点和终止点
    breaks_extended <- c(0, breaks, length(days_data))
    lengths <- diff(breaks_extended)
    max_consecutive <- max(lengths)
  }
  
  return(max_consecutive)
}

# 筛选高质量参与者
select_high_quality_participants <- function(participant_quality, criteria = NULL) {
  if (is.null(criteria)) {
    # 默认筛选标准
    criteria <- list(
      min_total_days = 20,            # 至少20天有数据
      min_presurgery_days = 3,        # 术前至少3天有数据
      min_postsurgery_days = 10,      # 术后至少10天有数据
      min_high_quality_days = 15,     # 至少15天为高质量数据(>=12小时)
      min_consecutive_days = 7,       # 至少有7天连续监测
      min_data_completeness = 0.4     # 整体完整率至少40%
    )
  }
  
  # 应用筛选标准
  high_quality_participants <- participant_quality %>%
    filter(
      total_days_with_data >= criteria$min_total_days,
      presurgery_days >= criteria$min_presurgery_days,
      postsurgery_days >= criteria$min_postsurgery_days,
      high_quality_days >= criteria$min_high_quality_days,
      max_consecutive_days >= criteria$min_consecutive_days,
      data_completeness >= criteria$min_data_completeness
    ) %>%
    # 按总体质量排序
    arrange(desc(total_days_with_data * data_completeness))
  
  return(high_quality_participants)
}

# 执行筛选
participant_quality <- evaluate_participant_quality(hourly_coverage)
high_quality_participants <- select_high_quality_participants(participant_quality)

# 打印筛选结果
print(high_quality_participants)

# 将结果保存到CSV文件
write.csv(participant_quality, "participant_quality_metrics.csv", row.names = FALSE)
write.csv(high_quality_participants, "high_quality_participants.csv", row.names = FALSE)

# 创建数据质量可视化
# 1. 参与者数据完整率柱状图
p_completeness <- ggplot(participant_quality, aes(x = reorder(subject_id, -data_completeness), 
                                                  y = data_completeness * 100)) +
  geom_bar(stat = "identity", fill = "#FAAF00", color = "black", alpha = 0.7) +
  geom_hline(yintercept = 40, linetype = "dashed", color = "red") +  # 40%阈值线
  theme_bw() +
  labs(
    title = "Participant data integrity rate",
    x = "Subject ID",
    y = "Data integrity rate (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank()
  )

# 2. 高质量参与者与所有参与者对比
comparison_data <- participant_quality %>%
  mutate(is_high_quality = subject_id %in% high_quality_participants$subject_id) %>%
  group_by(is_high_quality) %>%
  summarise(
    avg_days = mean(total_days_with_data),
    avg_completeness = mean(data_completeness) * 100,
    avg_high_quality = mean(high_quality_days),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(avg_days, avg_completeness, avg_high_quality),
               names_to = "metric", values_to = "value") %>%
  mutate(
    group = ifelse(is_high_quality, "High quality participants", "All participants"),
    metric = case_when(
      metric == "avg_days" ~ "Average monitoring days",
      metric == "avg_completeness" ~ "Average integrity (%)",
      metric == "avg_high_quality" ~ "Average high quality days"
    )
  )

p_comparison <- ggplot(comparison_data, aes(x = metric, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(
    title = "High quality participants vs All participants",
    x = NULL,
    y = "Value",
    fill = NULL
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = c("High quality participants" = "#FAAF00", "All participants" = "#A9A9A9"))

# 打印图表
print(p_completeness)
print(p_comparison)

# 保存图表
ggsave("participant_data_completeness.pdf", plot = p_completeness, width = 10, height = 6)
ggsave("quality_comparison.pdf", plot = p_comparison, width = 8, height = 6)

# 创建高质量参与者热图
high_quality_coverage <- hourly_coverage %>%
  filter(subject_id %in% high_quality_participants$subject_id)

p_hq_heatmap <- create_hr_coverage_heatmap(high_quality_coverage)
print(p_hq_heatmap)
ggsave("high_quality_participants_heatmap.pdf", plot = p_hq_heatmap, width = 8, height = 7)

# 返回筛选结果摘要
summary_results <- data.frame(
  total_participants = length(unique(hourly_coverage$subject_id)),
  high_quality_count = nrow(high_quality_participants),
  percentage = nrow(high_quality_participants) / length(unique(hourly_coverage$subject_id)) * 100
)

print(summary_results)





# 函数用于计算术前时间段覆盖统计
analyze_presurgery_time_period_coverage <- function(heart_rate_data, baseline_info) {
  # 获取心率数据
  hr_df <- heart_rate_data@sample_info %>%
    as.data.frame()
  
  # 处理基线信息
  baseline_info <- baseline_info %>%
    mutate(
      surgery_time_1 = as.Date(surgery_time_1),
      ID = toupper(ID)
    ) %>%
    dplyr::select(ID, surgery_time_1)
  
  # 计算每天每个时间段的小时覆盖率
  time_period_coverage <- hr_df %>%
    # 连接手术日期
    left_join(baseline_info, by = c("subject_id" = "ID")) %>%
    # 计算相对于手术的天数
    mutate(
      days_to_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days")),
      day_point = floor(days_to_surgery),
      hour = hour(measure_time),
      # 定义时间段
      time_period = case_when(
        hour >= 6 & hour < 18 ~ "daytime",
        TRUE ~ "nighttime"  # 覆盖 18-23 和 0-5 小时
      )
    ) %>%
    # 仅筛选术前10天
    filter(
      day_point >= -10,
      day_point < 0
    ) %>%
    # 计算每个时间段的小时覆盖率
    group_by(subject_id, day_point, time_period) %>%
    summarise(
      hours_covered = n_distinct(hour),
      .groups = "drop"
    )
  
  # 创建完整的网格，包含所有主题-天-时间段组合
  all_subjects <- unique(hr_df$subject_id)
  all_days <- seq(-10, -1)
  all_periods <- c("daytime", "nighttime")
  complete_grid <- expand.grid(
    subject_id = all_subjects,
    day_point = all_days,
    time_period = all_periods
  )
  
  # 连接实际覆盖率并填充缺失值为0
  final_coverage <- complete_grid %>%
    left_join(time_period_coverage, by = c("subject_id", "day_point", "time_period")) %>%
    mutate(hours_covered = coalesce(hours_covered, 0))
  
  # 计算每个参与者每天的时间段覆盖情况
  daily_coverage_status <- final_coverage %>%
    # 将数据转换为宽格式，每个时间段有一列
    pivot_wider(
      id_cols = c(subject_id, day_point),
      names_from = time_period,
      values_from = hours_covered
    ) %>%
    # 标记是否两个时间段都满足阈值
    mutate(
      meets_daytime_threshold = daytime >= 6,
      meets_nighttime_threshold = nighttime >= 6,
      meets_both_thresholds = meets_daytime_threshold & meets_nighttime_threshold
    )
  
  # 计算每天有多少人满足条件
  daily_participant_count <- daily_coverage_status %>%
    group_by(day_point) %>%
    summarise(
      total_participants = n_distinct(subject_id),
      participants_with_daytime = sum(meets_daytime_threshold),
      participants_with_nighttime = sum(meets_nighttime_threshold),
      participants_with_both = sum(meets_both_thresholds),
      percentage_with_both = participants_with_both / total_participants * 100,
      .groups = "drop"
    ) %>%
    arrange(day_point)
  
  # 计算每个参与者有多少天满足条件
  participant_day_count <- daily_coverage_status %>%
    group_by(subject_id) %>%
    summarise(
      total_days = n(),
      days_with_daytime = sum(meets_daytime_threshold),
      days_with_nighttime = sum(meets_nighttime_threshold),
      days_with_both = sum(meets_both_thresholds),
      percentage_days_with_both = days_with_both / total_days * 100,
      .groups = "drop"
    ) %>%
    arrange(desc(days_with_both))
  
  # 返回结果列表
  return(list(
    daily_coverage = final_coverage,
    daily_participant_count = daily_participant_count,
    participant_day_count = participant_day_count,
    daily_coverage_status = daily_coverage_status
  ))
}

# 执行分析
presurgery_results <- analyze_presurgery_time_period_coverage(heart_rate_data, baseline_info)

# 打印每天满足条件的参与者数量
print(presurgery_results$daily_participant_count)

# 保存结果
write.csv(presurgery_results$daily_participant_count, "presurgery_daily_time_coverage_count.csv", row.names = FALSE)
write.csv(presurgery_results$participant_day_count, "presurgery_participant_time_coverage.csv", row.names = FALSE)

# 创建每天满足条件参与者数量的柱状图
create_daily_participant_count_plot <- function(daily_participant_count) {
  # 将数据转换为长格式用于ggplot绘图
  plot_data <- daily_participant_count %>%
    dplyr::select(day_point, participants_with_daytime, participants_with_nighttime, participants_with_both) %>%
    pivot_longer(
      cols = c(participants_with_daytime, participants_with_nighttime, participants_with_both),
      names_to = "coverage_type",
      values_to = "participant_count"
    ) %>%
    mutate(
      coverage_type = case_when(
        coverage_type == "participants_with_daytime" ~ "Daytime ≥ 6h",
        coverage_type == "participants_with_nighttime" ~ "Nighttime ≥ 6h",
        coverage_type == "participants_with_both" ~ "Both periods ≥ 6h"
      ),
      # 确保因子级别顺序正确
      coverage_type = factor(coverage_type, 
                             levels = c("Daytime ≥ 6h", "Nighttime ≥ 6h", "Both periods ≥ 6h"))
    )
  
  # 创建柱状图
  p <- ggplot(plot_data, aes(x = day_point, y = participant_count, fill = coverage_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = participant_count), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, 
              size = 3) +
    scale_fill_manual(values = c("Daytime ≥ 6h" = "#FFD580", 
                                 "Nighttime ≥ 6h" = "#4682B4", 
                                 "Both periods ≥ 6h" = "#FAAF00")) +
    theme_bw() +
    labs(
      title = "Number of Participants Meeting Time Period Coverage Thresholds",
      subtitle = "Pre-surgery days (-10 to -1)",
      x = "Days Relative to Surgery",
      y = "Number of Participants",
      fill = "Coverage Type"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0),
      legend.position = "top"
    ) +
    scale_x_continuous(breaks = seq(-10, -1)) +
    scale_y_continuous(limits = c(0, max(daily_participant_count$total_participants) * 1.1))
  
  return(p)
}

# 创建百分比柱状图
create_daily_percentage_plot <- function(daily_participant_count) {
  p <- ggplot(daily_participant_count, aes(x = day_point, y = percentage_with_both)) +
    geom_bar(stat = "identity", fill = "#FAAF00", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", percentage_with_both)), 
              vjust = -0.5, 
              size = 3) +
    theme_bw() +
    labs(
      title = "Percentage of Participants with ≥6 Hours in Both Time Periods",
      subtitle = "Pre-surgery days (-10 to -1)",
      x = "Days Relative to Surgery",
      y = "Percentage of Participants"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0)
    ) +
    scale_x_continuous(breaks = seq(-10, -1)) +
    scale_y_continuous(limits = c(0, 100))
  
  return(p)
}

# 生成并保存图表
p_count <- create_daily_participant_count_plot(presurgery_results$daily_participant_count)
p_percentage <- create_daily_percentage_plot(presurgery_results$daily_participant_count)

print(p_count)
print(p_percentage)

ggsave("presurgery_daily_participant_count.pdf", plot = p_count, width = 10, height = 6)
ggsave("presurgery_daily_percentage.pdf", plot = p_percentage, width = 10, height = 6)

# 创建热图展示每个参与者的数据覆盖情况
create_presurgery_heatmap <- function(daily_coverage_status) {
  # 准备热图数据
  heatmap_data <- daily_coverage_status %>%
    mutate(
      coverage_status = case_when(
        meets_both_thresholds ~ "Both ≥ 6h",
        meets_daytime_threshold ~ "Only Daytime ≥ 6h",
        meets_nighttime_threshold ~ "Only Nighttime ≥ 6h",
        TRUE ~ "Neither ≥ 6h"
      ),
      coverage_status = factor(coverage_status, 
                               levels = c("Both ≥ 6h", "Only Daytime ≥ 6h", 
                                          "Only Nighttime ≥ 6h", "Neither ≥ 6h"))
    )
  
  # 排序参与者ID
  heatmap_data$subject_id <- factor(
    heatmap_data$subject_id,
    levels = rev(unique(heatmap_data$subject_id)[order(as.numeric(gsub("SH0*", "", unique(heatmap_data$subject_id))))])
  )
  
  # 创建热图
  p <- ggplot(heatmap_data, aes(x = day_point, y = subject_id, fill = coverage_status)) +
    geom_tile(color = "white", size = 0.3) +
    scale_fill_manual(values = c(
      "Both ≥ 6h" = "#FAAF00",
      "Only Daytime ≥ 6h" = "#FFD580",
      "Only Nighttime ≥ 6h" = "#4682B4",
      "Neither ≥ 6h" = "#F5F5F5"
    )) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Presurgery Time Period Coverage Status by Participant",
      x = "Days Relative to Surgery",
      y = "Subject ID",
      fill = "Coverage Status"
    ) +
    scale_x_continuous(breaks = seq(-10, -1))
  
  return(p)
}

# 生成并保存热图
p_heatmap <- create_presurgery_heatmap(presurgery_results$daily_coverage_status)
print(p_heatmap)
ggsave("presurgery_coverage_heatmap.pdf", plot = p_heatmap, width = 10, height = 12)

# 总结统计
summary_stats <- presurgery_results$daily_participant_count %>%
  summarise(
    avg_participants_with_both = mean(participants_with_both),
    max_participants_with_both = max(participants_with_both),
    min_participants_with_both = min(participants_with_both),
    avg_percentage_with_both = mean(percentage_with_both),
    .groups = "drop"
  )

print(summary_stats)
