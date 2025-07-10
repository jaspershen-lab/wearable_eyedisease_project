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
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_fill_gradient2(
      low = "#f0ead2",
      mid = "#dde5b4",
      high = "#adc178",
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



# Using patchwork for perfectly aligned plots with matching axes
library(patchwork)
library(grid)
library(ggplot2)

# 修改 create_perfect_alignment_plot 函数，使上方和右侧柱状图大小一致
create_perfect_alignment_plot <- function(hourly_coverage) {
  # Sort subject IDs numerically 
  hourly_coverage$subject_id <- factor(
    hourly_coverage$subject_id,
    levels = rev(unique(hourly_coverage$subject_id)[order(as.numeric(gsub("SH0*", "", unique(hourly_coverage$subject_id))))])
  )
  
  # Calculate summary data for the additional plots
  # Daily totals for each day (for top bar chart)
  daily_totals <- hourly_coverage %>%
    group_by(day_point) %>%
    summarise(daily_hours = sum(hours_covered), .groups = "drop")
  
  # Total hours per subject across all days (for right bar chart)
  subject_totals <- hourly_coverage %>%
    group_by(subject_id) %>%
    summarise(total_hours = sum(hours_covered), .groups = "drop")
  
  # Get the range of days and unique subject IDs for consistent axis limits
  day_range <- range(hourly_coverage$day_point)
  unique_subjects <- levels(hourly_coverage$subject_id)
  
  # Define shared theme elements for consistency
  base_theme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 1. Create main heatmap with exact dimensions and fixed plot ratio
  p_main <- ggplot(hourly_coverage, aes(x = day_point, y = subject_id)) +
    geom_point(aes(size = hours_covered, 
                   color = hours_covered,
                   fill = hours_covered),
               shape = 16,
               alpha = 0.7) +
    scale_color_gradient2(
      low = "#c8d6b0",
      mid = "#769f4a", 
      high = "#5a7836",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_fill_gradient2(
      low = "#c8d6b0",
      mid = "#769f4a",
      high = "#5a7836", 
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_size(
      range = c(0.1, 4),
      limits = c(0, 24),
      name = "Hours with\nData"
    ) +
    # Set exact limits to match side plots - critical for alignment
    scale_x_continuous(limits = c(day_range[1], day_range[2]), 
                       breaks = seq(-5, 30, by = 5),
                       expand = c(0.01, 0.01)) +
    base_theme +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      legend.position = "right"
    ) +
    labs(
      title = "Heart Rate Data Coverage Pattern - PPV Diabetes Patients",
      subtitle = paste("N =", length(unique_subjects), "participants"),
      x = "Days Relative to Surgery",
      y = "Patient ID"
    )
  
  # 2. Create top bar chart with EXACTLY matching x-axis - 调整高度比例
  max_daily_hours <- max(daily_totals$daily_hours) * 1.05
  p_top <- ggplot(daily_totals, aes(x = day_point, y = daily_hours)) +
    geom_col(fill = "#d9d0b4") +
    # Match the x-axis exactly with main plot - critical for alignment
    scale_x_continuous(limits = c(day_range[1], day_range[2]), 
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
  
  # 3. Create right bar chart with EXACTLY matching y-axis - 调整宽度比例
  max_total_hours <- max(subject_totals$total_hours) * 1.05
  p_right <- ggplot(subject_totals, aes(y = subject_id, x = total_hours)) +
    geom_col(fill = "#d9d0b4") +
    # Set y-axis to exactly match main plot - critical for alignment
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
  
  # 4. Create empty plot for top-right corner
  p_empty <- ggplot() + 
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # 关键修改：调整layout的widths和heights使上方和右侧图表大小一致
  combined_plot <- (p_top | p_empty) / (p_main | p_right) +
    plot_layout(
      widths = c(4, 1),      # 保持主图和右侧图的宽度比例
      heights = c(1, 4),     # 修改：使上方图和主图的高度比例为1:4，与右侧图的宽度比例一致
      guides = "collect"     # Collect all legends
    ) &
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))  # Remove outer margins
  
  return(combined_plot)
}

# 使用修正的函数重新创建图表
combined_plot <- create_perfect_alignment_plot(hourly_coverage)
print(combined_plot)

# 保存修正的图表
ggsave(filename = "heart_rate_data_perfectly_aligned_fixed.pdf", 
       plot = combined_plot, 
       width = 12,       # 稍微调整宽度
       height = 10,      # 稍微调整高度
       dpi = 300)

# 另一种方法：如果您想要更精确的控制，可以使用相等的比例
create_perfect_alignment_plot_equal <- function(hourly_coverage) {
  # 前面的代码保持不变...
  hourly_coverage$subject_id <- factor(
    hourly_coverage$subject_id,
    levels = rev(unique(hourly_coverage$subject_id)[order(as.numeric(gsub("SH0*", "", unique(hourly_coverage$subject_id))))])
  )
  
  daily_totals <- hourly_coverage %>%
    group_by(day_point) %>%
    summarise(daily_hours = sum(hours_covered), .groups = "drop")
  
  subject_totals <- hourly_coverage %>%
    group_by(subject_id) %>%
    summarise(total_hours = sum(hours_covered), .groups = "drop")
  
  day_range <- range(hourly_coverage$day_point)
  unique_subjects <- levels(hourly_coverage$subject_id)
  
  base_theme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 主图
  p_main <- ggplot(hourly_coverage, aes(x = day_point, y = subject_id)) +
    geom_point(aes(size = hours_covered, 
                   color = hours_covered,
                   fill = hours_covered),
               shape = 16,
               alpha = 0.7) +
    scale_color_gradient2(
      low = "#c8d6b0",
      mid = "#769f4a", 
      high = "#5a7836",
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_fill_gradient2(
      low = "#c8d6b0",
      mid = "#769f4a",
      high = "#5a7836", 
      midpoint = 12,
      name = "Hours with\nData"
    ) +
    scale_size(
      range = c(0.1, 4),
      limits = c(0, 24),
      name = "Hours with\nData"
    ) +
    scale_x_continuous(limits = c(day_range[1], day_range[2]), 
                       breaks = seq(-5, 30, by = 5),
                       expand = c(0.01, 0.01)) +
    base_theme +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      legend.position = "right"
    ) +
    labs(
      title = "Heart Rate Data Coverage Pattern - PPV Diabetes Patients",
      subtitle = paste("N =", length(unique_subjects), "participants"),
      x = "Days Relative to Surgery",
      y = "Patient ID"
    )
  
  # 上方图
  max_daily_hours <- max(daily_totals$daily_hours) * 1.05
  p_top <- ggplot(daily_totals, aes(x = day_point, y = daily_hours)) +
    geom_col(fill = "#d9d0b4") +
    scale_x_continuous(limits = c(day_range[1], day_range[2]), 
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
  
  # 右侧图
  max_total_hours <- max(subject_totals$total_hours) * 1.05
  p_right <- ggplot(subject_totals, aes(y = subject_id, x = total_hours)) +
    geom_col(fill = "#d9d0b4") +
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
  
  p_empty <- ggplot() + 
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  # 使用相等的比例 - 1:1的比例使上方和右侧图表大小一致
  combined_plot <- (p_top | p_empty) / (p_main | p_right) +
    plot_layout(
      widths = c(3, 1),      # 主图:右侧图 = 3:1
      heights = c(1, 3),     # 上方图:主图 = 1:3，与宽度比例一致
      guides = "collect"
    ) &
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
  
  return(combined_plot)
}

# 创建相等比例的图表
combined_plot_equal <- create_perfect_alignment_plot_equal(hourly_coverage)
print(combined_plot_equal)

# 保存相等比例的图表
ggsave(filename = "heart_rate_data_perfectly_aligned_equal.pdf", 
       plot = combined_plot_equal, 
       width = 10,
       height = 10,       # 使用正方形比例
       dpi = 300)

cat("Generated aligned plots with consistent bar chart sizes!\n")






# Create histogram of participants with data per day
participants_per_day <- hourly_coverage %>%
  filter(hours_covered > 0) %>%
  group_by(day_point) %>%
  summarise(
    participant_count = n_distinct(subject_id)
  )

# Create the histogram plot
p_histogram <- ggplot(participants_per_day, aes(x = day_point, y = participant_count)) +
  geom_bar(stat = "identity", fill = "#d9d0b4", color = "black", alpha = 0.7) +
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
  geom_bar(stat = "identity", fill = "#d9d0b4", color = "black", alpha = 0.7) +
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







# # 基于已有的hourly_coverage数据，筛选高质量参与者
# 
# # 定义函数来评估参与者数据质量
# evaluate_participant_quality <- function(hourly_coverage) {
#   
#   # 计算每个参与者的数据质量指标
#   participant_quality <- hourly_coverage %>%
#     group_by(subject_id) %>%
#     summarise(
#       # 1. 总体覆盖天数（有任何数据的天数）
#       total_days_with_data = sum(hours_covered > 0),
#       
#       # 2. 术前覆盖天数
#       presurgery_days = sum(day_point < 0 & hours_covered > 0),
#       
#       # 3. 术后覆盖天数
#       postsurgery_days = sum(day_point >= 0 & hours_covered > 0),
#       
#       # 4. 高质量天数（每天至少12小时数据）
#       high_quality_days = sum(hours_covered >= 12),
#       
#       # 5. 术前高质量天数
#       presurgery_high_quality = sum(day_point < 0 & hours_covered >= 12),
#       
#       # 6. 术后高质量天数
#       postsurgery_high_quality = sum(day_point >= 0 & hours_covered >= 12),
#       
#       # 7. 平均每天覆盖小时数
#       mean_hours_per_day = mean(hours_covered),
#       
#       # 8. 连续性指标 - 最长连续监测天数（有一定数据，如至少6小时）
#       max_consecutive_days = {
#         # 找出有数据的天数（至少6小时）
#         days_with_data <- day_point[hours_covered >= 6]
#         calculate_max_consecutive_days(days_with_data)
#       },
#       
#       # 9. 整体数据完整率（总记录小时数/理论最大小时数）
#       data_completeness = sum(hours_covered) / (n() * 24),
#       
#       .groups = "drop"
#     )
#   
#   return(participant_quality)
# }
# 
# # 辅助函数：计算最长连续监测天数
# calculate_max_consecutive_days <- function(days_data) {
#   # 确保天数已排序
#   days_data <- sort(days_data)
#   
#   if (length(days_data) == 0) {
#     return(0)
#   }
#   
#   # 计算连续天数
#   gaps <- diff(days_data)
#   # 找出所有间隔为1的位置（连续的天数）
#   breaks <- which(gaps != 1)
#   
#   # 计算每段连续天数的长度
#   if (length(breaks) == 0) {
#     # 如果没有断点，则所有天数都是连续的
#     max_consecutive <- length(days_data)
#   } else {
#     # 添加起始点和终止点
#     breaks_extended <- c(0, breaks, length(days_data))
#     lengths <- diff(breaks_extended)
#     max_consecutive <- max(lengths)
#   }
#   
#   return(max_consecutive)
# }
# 
# # 筛选高质量参与者
# select_high_quality_participants <- function(participant_quality, criteria = NULL) {
#   if (is.null(criteria)) {
#     # 默认筛选标准
#     criteria <- list(
#       min_total_days = 20,            # 至少20天有数据
#       min_presurgery_days = 3,        # 术前至少3天有数据
#       min_postsurgery_days = 10,      # 术后至少10天有数据
#       min_high_quality_days = 15,     # 至少15天为高质量数据(>=12小时)
#       min_consecutive_days = 7,       # 至少有7天连续监测
#       min_data_completeness = 0.4     # 整体完整率至少40%
#     )
#   }
#   
#   # 应用筛选标准
#   high_quality_participants <- participant_quality %>%
#     filter(
#       total_days_with_data >= criteria$min_total_days,
#       presurgery_days >= criteria$min_presurgery_days,
#       postsurgery_days >= criteria$min_postsurgery_days,
#       high_quality_days >= criteria$min_high_quality_days,
#       max_consecutive_days >= criteria$min_consecutive_days,
#       data_completeness >= criteria$min_data_completeness
#     ) %>%
#     # 按总体质量排序
#     arrange(desc(total_days_with_data * data_completeness))
#   
#   return(high_quality_participants)
# }
# 
# # 执行筛选
# participant_quality <- evaluate_participant_quality(hourly_coverage)
# high_quality_participants <- select_high_quality_participants(participant_quality)
# 
# # 打印筛选结果
# print(high_quality_participants)
# 
# # 将结果保存到CSV文件
# write.csv(participant_quality, "participant_quality_metrics.csv", row.names = FALSE)
# write.csv(high_quality_participants, "high_quality_participants.csv", row.names = FALSE)
# 
# # 创建数据质量可视化
# # 1. 参与者数据完整率柱状图
# p_completeness <- ggplot(participant_quality, aes(x = reorder(subject_id, -data_completeness), 
#                                                   y = data_completeness * 100)) +
#   geom_bar(stat = "identity", fill = "#FAAF00", color = "black", alpha = 0.7) +
#   geom_hline(yintercept = 40, linetype = "dashed", color = "red") +  # 40%阈值线
#   theme_bw() +
#   labs(
#     title = "Participant data integrity rate",
#     x = "Subject ID",
#     y = "Data integrity rate (%)"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#     panel.grid.minor = element_blank()
#   )
# 
# # 2. 高质量参与者与所有参与者对比
# comparison_data <- participant_quality %>%
#   mutate(is_high_quality = subject_id %in% high_quality_participants$subject_id) %>%
#   group_by(is_high_quality) %>%
#   summarise(
#     avg_days = mean(total_days_with_data),
#     avg_completeness = mean(data_completeness) * 100,
#     avg_high_quality = mean(high_quality_days),
#     .groups = "drop"
#   ) %>%
#   pivot_longer(cols = c(avg_days, avg_completeness, avg_high_quality),
#                names_to = "metric", values_to = "value") %>%
#   mutate(
#     group = ifelse(is_high_quality, "High quality participants", "All participants"),
#     metric = case_when(
#       metric == "avg_days" ~ "Average monitoring days",
#       metric == "avg_completeness" ~ "Average integrity (%)",
#       metric == "avg_high_quality" ~ "Average high quality days"
#     )
#   )
# 
# p_comparison <- ggplot(comparison_data, aes(x = metric, y = value, fill = group)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_bw() +
#   labs(
#     title = "High quality participants vs All participants",
#     x = NULL,
#     y = "Value",
#     fill = NULL
#   ) +
#   theme(
#     legend.position = "top",
#     panel.grid.minor = element_blank()
#   ) +
#   scale_fill_manual(values = c("High quality participants" = "#FAAF00", "All participants" = "#A9A9A9"))
# 
# # 打印图表
# print(p_completeness)
# print(p_comparison)
# 
# # 保存图表
# ggsave("participant_data_completeness.pdf", plot = p_completeness, width = 10, height = 6)
# ggsave("quality_comparison.pdf", plot = p_comparison, width = 8, height = 6)
# 
# # 创建高质量参与者热图
# high_quality_coverage <- hourly_coverage %>%
#   filter(subject_id %in% high_quality_participants$subject_id)
# 
# p_hq_heatmap <- create_hr_coverage_heatmap(high_quality_coverage)
# print(p_hq_heatmap)
# ggsave("high_quality_participants_heatmap.pdf", plot = p_hq_heatmap, width = 8, height = 7)
# 
# # 返回筛选结果摘要
# summary_results <- data.frame(
#   total_participants = length(unique(hourly_coverage$subject_id)),
#   high_quality_count = nrow(high_quality_participants),
#   percentage = nrow(high_quality_participants) / length(unique(hourly_coverage$subject_id)) * 100
# )
# 
# print(summary_results)





# 修改函数用于计算术前7天到术后30天的时间段覆盖统计
analyze_extended_time_period_coverage <- function(heart_rate_data, baseline_info) {
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
    # 修改为术前7天到术后30天
    filter(
      day_point >= -7,
      day_point <= 30
    ) %>%
    # 计算每个时间段的小时覆盖率
    group_by(subject_id, day_point, time_period) %>%
    summarise(
      hours_covered = n_distinct(hour),
      .groups = "drop"
    )
  
  # 创建完整的网格，包含所有主题-天-时间段组合
  all_subjects <- unique(hr_df$subject_id)
  all_days <- seq(-7, 30)  # 修改为-7到30
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
extended_results <- analyze_extended_time_period_coverage(heart_rate_data, baseline_info)

# 保存结果
write.csv(extended_results$daily_participant_count, "extended_daily_time_coverage_count.csv", row.names = FALSE)
write.csv(extended_results$participant_day_count, "extended_participant_time_coverage.csv", row.names = FALSE)

# 修改创建柱状图函数以适应新的时间范围
create_extended_participant_count_plot <- function(daily_participant_count) {
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
    scale_fill_manual(values = c("Daytime ≥ 6h" = "#d6eadf", 
                                 "Nighttime ≥ 6h" = "#eac4d5", 
                                 "Both periods ≥ 6h" = "#95b8d1")) +
    theme_bw() +
    labs(
      title = "Number of Participants Meeting Time Period Coverage Thresholds",
      subtitle = "From 7 days before to 30 days after surgery",
      x = "Days Relative to Surgery",
      y = "Number of Participants",
      fill = "Coverage Type"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0),
      legend.position = "top"
    ) +
    scale_x_continuous(breaks = seq(-5, 30, by = 5)) +  # 更新刻度以匹配更大的范围
    scale_y_continuous(limits = c(0, max(daily_participant_count$total_participants) * 1.1))
  
  return(p)
}

# 修改百分比柱状图
create_extended_percentage_plot <- function(daily_participant_count) {
  p <- ggplot(daily_participant_count, aes(x = day_point, y = percentage_with_both)) +
    geom_bar(stat = "identity", fill = "#95b8d1", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", percentage_with_both)), 
              vjust = -0.5, 
              size = 2.5) +  # 减小字体大小以适应更多的条形
    theme_bw() +
    labs(
      title = "Percentage of Participants with ≥6 Hours in Both Time Periods",
      subtitle = "From 7 days before to 30 days after surgery",
      x = "Days Relative to Surgery",
      y = "Percentage of Participants"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)  # 旋转标签以适应更多数据点
    ) +
    scale_x_continuous(breaks = seq(-5, 30, by = 5)) +  # 更新刻度以匹配更大的范围
    scale_y_continuous(limits = c(0, 100))
  
  return(p)
}

# 修改热图函数以适应新的时间范围
create_extended_heatmap <- function(daily_coverage_status) {
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
      "Both ≥ 6h" = "#95b8d1",
      "Only Daytime ≥ 6h" = "#d6eadf",
      "Only Nighttime ≥ 6h" = "#eac4d5",
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
      title = "Time Period Coverage Status by Participant",
      subtitle = "From 7 days before to 30 days after surgery",
      x = "Days Relative to Surgery",
      y = "Subject ID",
      fill = "Coverage Status"
    ) +
    scale_x_continuous(breaks = seq(-5, 30, by = 5))  # 更新刻度以匹配更大的范围
  
  return(p)
}

# 生成并保存图表
p_extended_count <- create_extended_participant_count_plot(extended_results$daily_participant_count)
p_extended_percentage <- create_extended_percentage_plot(extended_results$daily_participant_count)
p_extended_heatmap <- create_extended_heatmap(extended_results$daily_coverage_status)

print(p_extended_count)
print(p_extended_percentage)
print(p_extended_heatmap)

# 保存图表
ggsave("extended_daily_participant_count.pdf", plot = p_extended_count, width = 12, height = 6)
ggsave("extended_daily_percentage.pdf", plot = p_extended_percentage, width = 12, height = 6)
ggsave("extended_coverage_heatmap.pdf", plot = p_extended_heatmap, width = 12, height = 10)

# 计算术前和术后的统计摘要
period_summary <- extended_results$daily_participant_count %>%
  mutate(period = ifelse(day_point < 0, "Pre-surgery", "Post-surgery")) %>%
  group_by(period) %>%
  summarise(
    avg_participants_with_both = mean(participants_with_both),
    max_participants_with_both = max(participants_with_both),
    min_participants_with_both = min(participants_with_both),
    avg_percentage_with_both = mean(percentage_with_both),
    days_count = n(),
    .groups = "drop"
  )

# 总体统计
overall_summary <- extended_results$daily_participant_count %>%
  summarise(
    avg_participants_with_both = mean(participants_with_both),
    max_participants_with_both = max(participants_with_both),
    min_participants_with_both = min(participants_with_both),
    avg_percentage_with_both = mean(percentage_with_both),
    days_count = n(),
    .groups = "drop"
  ) %>%
  mutate(period = "Overall")

# 合并统计信息
combined_summary <- bind_rows(period_summary, overall_summary)
print(combined_summary)
write.csv(combined_summary, "extended_time_period_summary.csv", row.names = FALSE)

# 创建术前与术后的比较柱状图
create_period_comparison_plot <- function(daily_participant_count) {
  # 计算术前和术后的平均值
  period_data <- daily_participant_count %>%
    mutate(period = ifelse(day_point < 0, "Pre-surgery (-7 to -1)", "Post-surgery (0 to 30)")) %>%
    group_by(period) %>%
    summarise(
      avg_total = mean(total_participants),
      avg_daytime = mean(participants_with_daytime),
      avg_nighttime = mean(participants_with_nighttime),
      avg_both = mean(participants_with_both),
      avg_percentage = mean(percentage_with_both),
      .groups = "drop"
    ) %>%
    # 将数据转换为长格式用于绘图
    pivot_longer(
      cols = starts_with("avg_"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    # 调整标签
    mutate(
      metric = case_when(
        metric == "avg_total" ~ "Total participants",
        metric == "avg_daytime" ~ "Daytime ≥ 6h",
        metric == "avg_nighttime" ~ "Nighttime ≥ 6h",
        metric == "avg_both" ~ "Both periods ≥ 6h",
        metric == "avg_percentage" ~ "% with both periods ≥ 6h"
      ),
      # 设置因子顺序
      metric = factor(metric, levels = c("Total participants", "Daytime ≥ 6h", 
                                         "Nighttime ≥ 6h", "Both periods ≥ 6h",
                                         "% with both periods ≥ 6h")),
      period = factor(period, levels = c("Pre-surgery (-7 to -1)", "Post-surgery (0 to 30)"))
    )
  
  # 分离百分比数据，因为它使用不同的y轴刻度
  count_data <- period_data %>% filter(metric != "% with both periods ≥ 6h")
  percentage_data <- period_data %>% filter(metric == "% with both periods ≥ 6h")
  
  # 创建计数比较图
  p_count <- ggplot(count_data, aes(x = metric, y = value, fill = period)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", value)), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, 
              size = 3) +
    scale_fill_manual(values = c("Pre-surgery (-7 to -1)" = "#809bce", 
                                 "Post-surgery (0 to 30)" = "#b8e0d4")) +
    theme_bw() +
    labs(
      title = "Average Participant Counts by Time Period",
      x = NULL,
      y = "Average Number of Participants",
      fill = "Period"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  # 创建百分比比较图
  p_percentage <- ggplot(percentage_data, aes(x = period, y = value, fill = period)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", value)), 
              vjust = -0.5, 
              size = 3) +
    scale_fill_manual(values = c("Pre-surgery (-7 to -1)" = "#809bce", 
                                 "Post-surgery (0 to 30)" = "#b8e0d4")) +
    theme_bw() +
    labs(
      title = "Average Percentage of Participants with Both Time Periods ≥ 6h",
      x = NULL,
      y = "Percentage"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0, 100))
  
  return(list(p_count = p_count, p_percentage = p_percentage))
}

# 生成比较图表
comparison_plots <- create_period_comparison_plot(extended_results$daily_participant_count)
print(comparison_plots$p_count)
print(comparison_plots$p_percentage)

# 保存比较图表
ggsave("period_count_comparison.pdf", plot = comparison_plots$p_count, width = 10, height = 6)
ggsave("period_percentage_comparison.pdf", plot = comparison_plots$p_percentage, width = 8, height = 6)
