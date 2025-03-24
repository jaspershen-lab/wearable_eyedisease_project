library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


########读取数据
load(
  "3_data_analysis/2_data_analysis/RHR/heart_rate_data_RHR.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############创建输出目录
dir.create("3_data_analysis/2_data_analysis/RHR/RHR_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/RHR/RHR_time_periods")


###############处理基线信息
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)


###############根据活动水平计算RHR
# 为指定标签计算RHR的函数
calculate_rhr_by_label <- function(data, label_value) {
  # 获取具有特定标签的样本信息
  filtered_sample_info <- data@sample_info %>%
    dplyr::filter(label == label_value)
  
  # 获取对应的心率数据
  heart_rate_values <- data@expression_data[1, filtered_sample_info$sample_id] %>%
    as.numeric()
  
  # 创建包含样本信息和心率值的数据框
  data.frame(
    sample_id = filtered_sample_info$sample_id,
    heart_rate = heart_rate_values,
    timestamp = filtered_sample_info$measure_time,
    subject_id = filtered_sample_info$subject_id
  ) %>%
    # 过滤合理的心率范围
    dplyr::filter(heart_rate >= 30 & heart_rate <= 200) %>%
    arrange(subject_id, timestamp)
}

# 为两种活动水平计算RHR
rhr_label_1 <- calculate_rhr_by_label(heart_rate_data, "<1")   # 静息状态(步数<=1)
rhr_label_50 <- calculate_rhr_by_label(heart_rate_data, "<50") # 轻度活动(步数<=50)

###############分时间窗口处理RHR数据
# 完整的RHR处理函数 - 包含所有时间窗口
process_rhr_data <- function(rhr_data, label_suffix) {
  # 计算基本时间信息
  base_data <- rhr_data %>%
    left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
    mutate(
      # 计算相对于手术时间的小时数（正值表示术后，负值表示术前）
      hours_to_surgery = as.numeric(difftime(timestamp, surgery_time_1, units = "hours"))
    )
  
  # 术前时间窗口处理 ---------------------------------
  
  # 术前3天 (hours -72 to 0)
  pre_3d <- base_data %>%
    filter(hours_to_surgery >= -72 & hours_to_surgery <= 0) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "pre_surgery_3d")
  
  # 术前3-7天 (hours -168 to -72)
  pre_3d_to_7d <- base_data %>%
    filter(hours_to_surgery >= -168 & hours_to_surgery < -72) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "pre_surgery_3d_to_7d")
  
  # 术前7天全部 (hours -168 to 0)
  pre_7d_all <- base_data %>%
    filter(hours_to_surgery >= -168 & hours_to_surgery <= 0) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "pre_surgery_7d_all")
  
  # 术前所有数据 (hours -720 to 0, 假设30天)
  pre_all <- base_data %>%
    filter(hours_to_surgery >= -720 & hours_to_surgery <= 0) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "pre_surgery_all")
  
  # 术后时间窗口处理 ---------------------------------
  
  # 术后1-3天 (hours 24 to 72)
  post_1to3d <- base_data %>%
    filter(hours_to_surgery > 24 & hours_to_surgery <= 72) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "post_surgery_1to3d")
  
  # 术后4-6天 (hours 72 to 144)
  post_4to6d <- base_data %>%
    filter(hours_to_surgery > 72 & hours_to_surgery <= 144) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "post_surgery_4to6d")
  
  # 术后6天 (hours 0 to 144) - 代替原来的7天
  post_6d <- base_data %>%
    filter(hours_to_surgery > 0 & hours_to_surgery <= 144) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "post_surgery_6d")
  
  # 术后第1周 (hours 0 to 168) - 保留原有的时间窗口
  post_1w <- base_data %>%
    filter(hours_to_surgery > 0 & hours_to_surgery <= 168) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "post_surgery_1w")
  
  # 术后7-30天 (hours 168 to 720)
  post_7d_to_30d <- base_data %>%
    filter(hours_to_surgery > 168 & hours_to_surgery <= 720) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "post_surgery_7d_to_30d")
  
  # 术后23-30天 (hours 552 to 720)
  post_day23_to_30 <- base_data %>%
    filter(hours_to_surgery > 552 & hours_to_surgery <= 720) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "post_surgery_day23_to_30")
  
  # 术后27-30天 (hours 648 to 720)
  post_day27_to_30 <- base_data %>%
    filter(hours_to_surgery > 648 & hours_to_surgery <= 720) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "post_surgery_day27_to_30")
  
  # 术后>30天 (hours > 720)
  post_over_30d <- base_data %>%
    filter(hours_to_surgery > 720) %>%
    group_by(subject_id) %>%
    summarise(
      mean_hr = mean(heart_rate),
      min_hr = min(heart_rate),
      max_hr = max(heart_rate),
      median_hr = median(heart_rate),
      q1_hr = quantile(heart_rate, 0.25),
      q3_hr = quantile(heart_rate, 0.75),
      iqr_hr = IQR(heart_rate),
      sd_hr = sd(heart_rate),
      n_measurements = n(),
      .groups = "drop"
    ) %>%
    mutate(time_period = "post_surgery_over_30d")
  
  # 合并所有时间窗口结果
  all_results <- bind_rows(
    pre_3d, 
    pre_3d_to_7d, 
    pre_7d_all, 
    pre_all,
    post_1to3d,
    post_4to6d,
    post_6d,
    post_1w,
    post_7d_to_30d,
    post_day23_to_30,
    post_day27_to_30,
    post_over_30d
  ) %>%
    mutate(label = label_suffix)
  
  return(all_results)
}

# 使用该函数处理两种活动水平的数据
rhr_summary_1 <- process_rhr_data(rhr_label_1, "steps_1")
rhr_summary_50 <- process_rhr_data(rhr_label_50, "steps_50")

# 合并两种活动水平的结果
combined_rhr_summary <- bind_rows(rhr_summary_1, rhr_summary_50) %>%
  mutate(
    label_time = paste(label, time_period, sep = "_")
  )

###############绘制RHR模式图
# 更新绘图函数，包含新添加的时间窗口
plot_rhr_patterns <- function(combined_rhr_summary) {
  # 处理数据
  plot_data <- combined_rhr_summary %>%
    mutate(
      steps = ifelse(grepl("steps_1", label_time), "1", "50"),
      time_period = case_when(
        grepl("pre_surgery_3d_to_7d", label_time) ~ "Pre 3d to 7d",
        grepl("pre_surgery_3d", label_time) ~ "Pre 3d",
        grepl("pre_surgery_7d_all", label_time) ~ "Pre 7d all",
        grepl("pre_surgery_all", label_time) ~ "Pre all",
        grepl("post_surgery_1to3d", label_time) ~ "Post 1-3d",    # 新增的时间窗
        grepl("post_surgery_4to6d", label_time) ~ "Post 4-6d",    # 新增的时间窗
        grepl("post_surgery_6d", label_time) ~ "Post 6d",         # 修改后的时间窗
        grepl("post_surgery_day27_to_30", label_time) ~ "Post day 27-30",
        grepl("post_surgery_day23_to_30", label_time) ~ "Post day 23-30",
        grepl("post_surgery_1w", label_time) ~ "Post 1w",
        grepl("post_surgery_7d_to_30d", label_time) ~ "Post 7d to 30d",
        grepl("post_surgery_over_30d", label_time) ~ "Post >30d",
        TRUE ~ "Other"
      )
    ) %>%
    # 计算每个时间段的平均统计量
    group_by(steps, time_period) %>%
    summarise(
      mean_rhr = mean(mean_hr, na.rm = TRUE),
      median_rhr = mean(median_hr, na.rm = TRUE),
      min_rhr = mean(min_hr, na.rm = TRUE),
      max_rhr = mean(max_hr, na.rm = TRUE),
      iqr_rhr = mean(iqr_hr, na.rm = TRUE),
      sd_rhr = mean(sd_hr, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    )
  
  # 设置时间窗口顺序，包括新添加的时间窗口
  plot_data$time_period <- factor(plot_data$time_period,
                                  levels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d",
                                             "Post 1-3d", "Post 4-6d", "Post 6d", "Post 1w", 
                                             "Post 7d to 30d", "Post day 23-30", "Post day 27-30", "Post >30d"))
  
  # 计算总体中位数
  median_rhr <- median(combined_rhr_summary$median_hr, na.rm = TRUE)
  
  # 创建平均RHR图
  p_mean <- ggplot(plot_data, aes(x = time_period, y = mean_rhr, 
                                  color = steps, group = steps)) +
    geom_line(size = 1) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = median_rhr, linetype = "dashed", 
               color = "red", alpha = 0.7) +
    annotate("text", x = 6, y = median_rhr, 
             label = "Overall Median RHR", hjust = 1, vjust = -0.5, 
             color = "red") +
    scale_y_continuous(limits = c(65, 85), breaks = seq(65, 85, 5)) +
    scale_color_manual(
      values = c("1" = "steelblue", "50" = "#FF9999"),
      name = "Max. steps",
      labels = c("≤1", "≤50")
    ) +
    labs(title = "Perioperative Mean RHR Patterns",
         x = "Time Period",
         y = "Mean RHR (bpm)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # 创建中位数RHR图
  p_median <- ggplot(plot_data, aes(x = time_period, y = median_rhr, 
                                    color = steps, group = steps)) +
    geom_line(size = 1) +
    geom_point(size = 3, alpha = 0.7) +
    scale_y_continuous(limits = c(65, 85), breaks = seq(65, 85, 5)) +
    scale_color_manual(
      values = c("1" = "steelblue", "50" = "#FF9999"),
      name = "Max. steps",
      labels = c("≤1", "≤50")
    ) +
    labs(title = "Perioperative Median RHR Patterns",
         x = "Time Period",
         y = "Median RHR (bpm)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # 创建包含所有指标的汇总表
  summary_table <- plot_data %>%
    dplyr::select(steps, time_period, mean_rhr, median_rhr, min_rhr, max_rhr, iqr_rhr, sd_rhr, n) %>%
    mutate(across(where(is.numeric), ~round(., 1)))
  
  return(list(mean_plot = p_mean, median_plot = p_median, summary = summary_table))
}

# 生成图表
rhr_plots <- plot_rhr_patterns(combined_rhr_summary)

# 保存平均RHR图
ggsave(
  "perioperative_mean_rhr_patterns.pdf",
  plot = rhr_plots$mean_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# 保存中位数RHR图
ggsave(
  "perioperative_median_rhr_patterns.pdf",
  plot = rhr_plots$median_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# 创建并保存汇总统计对象
rhr_detailed_summary <- rhr_plots$summary
save(rhr_detailed_summary, file = "rhr_detailed_summary.rda", compress = "xz")

# 打印汇总表
print("RHR Statistics Summary by Time Period and Activity Level:")
print(rhr_plots$summary)

###############创建宽格式数据（用于后续分析）
create_wide_format_rhr <- function(combined_rhr_summary) {
  # 创建基本数据框，包含所有唯一的subject_id
  subject_ids <- unique(combined_rhr_summary$subject_id)
  result <- data.frame(subject_id = subject_ids)
  
  # 尝试从baseline_info获取手术日期
  if(exists("baseline_info_processed")) {
    surgery_dates <- baseline_info_processed %>%
      mutate(surgery_date = format(surgery_time_1, "%Y-%m-%d")) %>%
      dplyr::select(ID, surgery_date)
    
    result <- result %>%
      left_join(surgery_dates, by = c("subject_id" = "ID"))
  } else {
    # 如果无法获取手术日期，添加空列
    result$surgery_date <- NA
  }
  
  # 遍历每个统计量和时间窗口的组合，创建宽格式数据
  stat_types <- c("mean", "median", "min", "max", "iqr", "sd", "n_measurements")
  
  for(row in 1:nrow(combined_rhr_summary)) {
    subj <- combined_rhr_summary$subject_id[row]
    time_period <- combined_rhr_summary$time_period[row]
    label <- combined_rhr_summary$label[row]
    steps_val <- sub("steps_", "", label)
    
    # 对于每个统计量
    for(stat in stat_types) {
      if(stat == "n_measurements") {
        # 直接使用n_measurements列
        val <- combined_rhr_summary$n_measurements[row]
        col_name <- paste0(time_period, "_", stat, "_steps_", steps_val)
      } else {
        # 使用hr后缀列，但在列名中添加rhr
        hr_col <- paste0(stat, "_hr")
        if(hr_col %in% names(combined_rhr_summary)) {
          val <- combined_rhr_summary[[hr_col]][row]
          # 新的列名格式，在统计类型后添加_rhr
          col_name <- paste0(time_period, "_", stat, "_rhr_steps_", steps_val)
          
          # 更新结果数据框
          result[result$subject_id == subj, col_name] <- val
        }
      }
    }
  }
  
  return(result)
}

# 创建并保存宽格式RHR数据
time_period_rhr_results <- create_wide_format_rhr(combined_rhr_summary)

# 验证列名格式
cat("新的RHR列名格式示例:\n")
print(names(time_period_rhr_results)[1:20])

# 检查新添加时间窗口的列数量
cat("\n检查新添加时间窗口的列数:\n")
cat("术后1-3天(post_surgery_1to3d)列数: ", 
    sum(grepl("post_surgery_1to3d", names(time_period_rhr_results))), "\n")
cat("术后4-6天(post_surgery_4to6d)列数: ", 
    sum(grepl("post_surgery_4to6d", names(time_period_rhr_results))), "\n")
cat("术后6天(post_surgery_6d)列数: ", 
    sum(grepl("post_surgery_6d", names(time_period_rhr_results))), "\n")

# 保存结果
save(time_period_rhr_results, 
     file = "time_period_rhr_results.rda", 
     compress = "xz")

# 保存原始汇总数据，以备后用
save(combined_rhr_summary,
     file = "combined_rhr_summary.rda", 
     compress = "xz")

# 输出确认信息
cat("\nRHR时间窗口数据处理完成，结果已保存。\n")
cat("总共处理了", length(unique(combined_rhr_summary$subject_id)), "个受试者的数据。\n")
cat("总共创建了", length(unique(combined_rhr_summary$time_period)), "个时间窗口。\n")

# 创建日心率模式图（每日心率变化）
calculate_pattern_1 <- function(data) {
  # 过滤步数<=1的样本
  filtered_samples <- data@sample_info %>%
    dplyr::filter(label == "<1")
  
  # 获取对应的心率值
  heart_rates <- data@expression_data[1, filtered_samples$sample_id] %>%
    as.numeric()
  
  # 创建数据框
  filtered_samples %>%
    mutate(
      hour = hour(measure_time) + minute(measure_time)/60,
      heart_rate = heart_rates
    ) %>%
    dplyr::filter(heart_rate >= 30 & heart_rate <= 200) %>%
    group_by(hour = floor(hour)) %>%
    summarise(
      median_rhr = median(heart_rate),
      n_measurements = n(),
      sd_rhr = sd(heart_rate),
      .groups = "drop"
    ) %>%
    mutate(label = "1")
}

calculate_pattern_50 <- function(data) {
  # 过滤步数<=50的样本
  filtered_samples <- data@sample_info %>%
    dplyr::filter(label %in% c("<1", "<50"))
  
  # 获取对应的心率值
  heart_rates <- data@expression_data[1, filtered_samples$sample_id] %>%
    as.numeric()
  
  # 创建数据框
  filtered_samples %>%
    mutate(
      hour = hour(measure_time) + minute(measure_time)/60,
      heart_rate = heart_rates
    ) %>%
    dplyr::filter(heart_rate >= 30 & heart_rate <= 200) %>%
    group_by(hour = floor(hour)) %>%
    summarise(
      median_rhr = median(heart_rate),
      n_measurements = n(),
      sd_rhr = sd(heart_rate),
      .groups = "drop"
    ) %>%
    mutate(label = "50")
}

# 计算两种活动水平的日心率模式
daily_rhr_pattern_1 <- calculate_pattern_1(heart_rate_data)
daily_rhr_pattern_50 <- calculate_pattern_50(heart_rate_data)

# 创建完整的小时序列
all_hours <- tibble(hour = 0:23)

# 对两种模式填充缺失小时
daily_rhr_pattern_1 <- all_hours %>%
  left_join(daily_rhr_pattern_1, by = "hour") %>%
  mutate(label = "1")

daily_rhr_pattern_50 <- all_hours %>%
  left_join(daily_rhr_pattern_50, by = "hour") %>%
  mutate(label = "50")

# 合并两种模式
daily_rhr_pattern <- bind_rows(daily_rhr_pattern_1, daily_rhr_pattern_50)

# 计算总体中位数RHR
overall_median_rhr <- {
  filtered_samples <- heart_rate_data@sample_info %>%
    dplyr::filter(label %in% c("<1", "<50"))
  heart_rates <- heart_rate_data@expression_data[1, filtered_samples$sample_id] %>%
    as.numeric()
  median(heart_rates)
}

# 创建图表
p_daily <- ggplot(daily_rhr_pattern, aes(x = hour, y = median_rhr, color = label, group = label)) +
  # 添加点和线
  geom_line(size = 1) +
  geom_point(size = 2, alpha = 0.6) +
  # 添加中位数参考线
  geom_hline(
    yintercept = overall_median_rhr,
    color = "red",
    linetype = "dashed"
  ) +
  # 添加中位数RHR标签
  annotate(
    "text",
    x = 23,
    y = overall_median_rhr,
    label = "Median RHR",
    hjust = 1,
    vjust = -0.5,
    color = "red"
  ) +
  # 自定义坐标轴
  scale_x_continuous(
    breaks = seq(0, 24, 4),
    limits = c(0, 24),
    name = "Time of day"
  ) +
  scale_y_continuous(
    name = "Median RHR",
    limits = c(55, 80)
  ) +
  # 自定义颜色
  scale_color_manual(
    values = c("1" = "steelblue", "50" = "#FF9999"),
    name = "Max. steps"
  ) +
  # 使用theme_bw并自定义
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.15),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank()   # 移除次网格线
  ) +
  # 添加标题
  labs(title = "Daily RHR Pattern (10-min window)")

# 保存图表
ggsave(
  "daily_rhr_pattern.pdf",
  plot = p_daily,
  width = 10,
  height = 6,
  dpi = 300
)

# 打印按小时汇总的统计数据
print("\n按小时汇总的测量数据:")
daily_rhr_pattern %>%
  group_by(label) %>%
  summarise(
    mean_rhr = mean(median_rhr, na.rm = TRUE),
    min_rhr = min(median_rhr, na.rm = TRUE),
    max_rhr = max(median_rhr, na.rm = TRUE),
    n_hours = sum(!is.na(median_rhr))
  ) %>%
  print()

# 创建心率直方图
create_heart_rate_data <- function(heart_rate_data) {
  # 获取所有心率数据
  filtered_samples <- heart_rate_data@sample_info %>%
    dplyr::filter(label %in% c("<1", "<50"))
  
  # 获取对应的心率值
  heart_rates <- heart_rate_data@expression_data[1, filtered_samples$sample_id] %>%
    as.numeric()
  
  # 创建数据框
  data.frame(heart_rate = heart_rates) %>%
    # 添加心率的合理范围过滤
    dplyr::filter(heart_rate >= 30 & heart_rate <= 200)
}

# 绘制直方图的函数
plot_heart_rate_histogram <- function(data) {
  # 创建心率数据框
  heart_rate_df <- create_heart_rate_data(data)
  
  # 计算中位数
  median_hr <- median(heart_rate_df$heart_rate, na.rm = TRUE)
  
  # 创建图形
  ggplot(heart_rate_df, aes(x = heart_rate)) +
    # 添加直方图
    geom_histogram(
      binwidth = 2,  # 每2bpm一个bin
      fill = "#a6c0d5",  # 浅蓝色填充
      color = "white",   # 白色边框
      boundary = 0       # 确保bin从整数开始
    ) +
    # 添加中位数垂直线
    geom_vline(
      xintercept = median_hr,
      color = "#16165d",  # 参考线
      size = 1
    ) +
    # 添加中位数值文本
    annotate(
      "text",
      x = median_hr,
      y = Inf,
      label = sprintf("%.0f bpm", median_hr),
      vjust = 2,
      size = 4
    ) +
    # 自定义坐标轴
    scale_x_continuous(
      name = "RHR",
      limits = c(50, 110),
      breaks = seq(50, 110, 20)
    ) +
    scale_y_continuous(
      name = "Frequency",
      expand = c(0, 0)  # 让y轴从0开始
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

# 创建并保存RHR直方图
p_hist <- plot_heart_rate_histogram(heart_rate_data)
ggsave(
  "RHR_histogram.pdf",
  plot = p_hist,
  width = 6,
  height = 4,
  dpi = 300
)
