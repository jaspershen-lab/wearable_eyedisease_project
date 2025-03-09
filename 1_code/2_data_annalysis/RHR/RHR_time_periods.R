library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


########read data
load(
  "3_data_analysis/2_data_analysis/RHR/heart_rate_data_RHR.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/RHR/RHR_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/RHR/RHR_time_periods")


###############rhr by time window
# Process baseline info
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)


# Function to calculate RHR for a specific label
calculate_rhr_by_label <- function(data, label_value) {
  # Get sample info with the specific label
  filtered_sample_info <- data@sample_info %>%
    dplyr::filter(label == label_value)
  
  # Get corresponding heart rate data
  heart_rate_values <- data@expression_data[1, filtered_sample_info$sample_id] %>%
    as.numeric()
  
  # Create data frame with sample info and heart rate values
  data.frame(
    sample_id = filtered_sample_info$sample_id,
    heart_rate = heart_rate_values,
    timestamp = filtered_sample_info$measure_time,
    subject_id = filtered_sample_info$subject_id
  ) %>%
    dplyr::filter(heart_rate >= 30 & heart_rate <= 200) %>%
    arrange(subject_id, timestamp)
}

# Calculate RHR for both labels
rhr_label_1 <- calculate_rhr_by_label(heart_rate_data, "<1")
rhr_label_50 <- calculate_rhr_by_label(heart_rate_data, "<50")

# Function to process RHR data with time periods - UPDATED with additional statistics
###############修改后的函数 process_rhr_data
# 添加了两个新的时间窗：术后第30天前3天和术后前一周
# 修改后的process_rhr_data函数
process_rhr_data <- function(rhr_data, label_suffix) {
  # Calculate basic time information
  base_data <- rhr_data %>%
    left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
    mutate(
      hours_to_surgery = as.numeric(difftime(surgery_time_1, timestamp, units = "hours"))
    )
  
  # Calculate data for pre-surgery 3d window
  pre_3d <- base_data %>%
    filter(hours_to_surgery >= 0 & hours_to_surgery < 72) %>%
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
  
  # Calculate data for pre-surgery 3d to 7d window
  pre_3d_to_7d <- base_data %>%
    filter(hours_to_surgery >= 72 & hours_to_surgery < 168) %>%
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
  
  # Calculate data for pre-surgery 7d all window
  pre_7d_all <- base_data %>%
    filter(hours_to_surgery >= 0 & hours_to_surgery < 168) %>%
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
  
  # Calculate data for all pre-surgery data
  pre_all <- base_data %>%
    filter(hours_to_surgery >= 0) %>%
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
  
  # Combine all pre-surgery data
  result <- bind_rows(pre_3d, pre_3d_to_7d, pre_7d_all, pre_all) %>%
    mutate(label = label_suffix)
  
  # 计算术后第30天前3天的时间窗 (day 27-30, hours -720 to -648)
  post_30d_3d <- base_data %>%
    filter(hours_to_surgery < 0 & hours_to_surgery >= -720 & hours_to_surgery > -792) %>%
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
  
  # 计算术后第30天前一周的时间窗 (day 23-30, hours -720 to -552)
  post_30d_1w <- base_data %>%
    filter(hours_to_surgery < 0 & hours_to_surgery >= -720 & hours_to_surgery > -720+168) %>%
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
  
  # 计算术后前一周的时间窗（前7天，hours 0 to -168）
  post_1w <- base_data %>%
    filter(hours_to_surgery < 0 & hours_to_surgery >= -168) %>%
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
  
  # 现有的术后7天、7-30天、>30天窗口计算
  post_data <- base_data %>%
    filter(hours_to_surgery < 0) %>%
    mutate(
      time_period = case_when(
        hours_to_surgery >= -168 ~ "post_surgery_7d",
        hours_to_surgery >= -720 & hours_to_surgery < -168 ~ "post_surgery_7d_to_30d",
        hours_to_surgery < -720 ~ "post_surgery_over_30d"
      )
    ) %>%
    group_by(subject_id, time_period) %>%
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
    mutate(label = label_suffix)
  
  # 合并所有结果，包括新增的时间窗
  all_results <- bind_rows(
    result, 
    post_30d_3d %>% mutate(label = label_suffix),
    post_30d_1w %>% mutate(label = label_suffix), 
    post_1w %>% mutate(label = label_suffix),
    post_data
  )
  
  return(all_results)
}

# 更新绘图函数，添加新的时间窗
plot_rhr_patterns <- function(combined_rhr_summary) {
  # Process data
  plot_data <- combined_rhr_summary %>%
    mutate(
      steps = ifelse(grepl("steps_1", label_time), "1", "50"),
      time_period = case_when(
        grepl("pre_surgery_3d_to_7d", label_time) ~ "Pre 3d to 7d",
        grepl("pre_surgery_3d", label_time) ~ "Pre 3d",
        grepl("pre_surgery_7d_all", label_time) ~ "Pre 7d all",
        grepl("pre_surgery_all", label_time) ~ "Pre all",
        grepl("post_surgery_day27_to_30", label_time) ~ "Post day 27-30", # 新增的时间窗
        grepl("post_surgery_day23_to_30", label_time) ~ "Post day 23-30", # 新增的时间窗
        grepl("post_surgery_1w", label_time) ~ "Post 1w",                # 新增的时间窗
        grepl("post_surgery_7d$", label_time) ~ "Post 7d",               # 使用正则表达式结尾匹配
        grepl("post_surgery_7d_to_30d", label_time) ~ "Post 7d to 30d",
        grepl("post_surgery_over_30d", label_time) ~ "Post >30d",
        TRUE ~ "Other"
      )
    ) %>%
    # Calculate average statistics for each time period
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
  
  # 设置时间窗口顺序，包括新增的时间窗
  plot_data$time_period <- factor(plot_data$time_period,
                                  levels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d",
                                             "Post 1w", "Post 7d", "Post 7d to 30d", "Post day 23-30", "Post day 27-30", "Post >30d"))
  
  # Calculate overall median
  median_rhr <- median(combined_rhr_summary$median_hr, na.rm = TRUE)
  
  # Create plot for mean RHR
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
  
  # Create plot for median RHR
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
  
  # Create a summary table with all metrics
  summary_table <- plot_data %>%
    dplyr::select(steps, time_period, mean_rhr, median_rhr, min_rhr, max_rhr, iqr_rhr, sd_rhr, n) %>%
    mutate(across(where(is.numeric), ~round(., 1)))
  
  return(list(mean_plot = p_mean, median_plot = p_median, summary = summary_table))
}

# 重新运行计算combined_rhr_summary的相关代码
# Process both datasets with the updated function
rhr_summary_1 <- process_rhr_data(rhr_label_1, "steps_1")
rhr_summary_50 <- process_rhr_data(rhr_label_50, "steps_50")

# Combine the results
combined_rhr_summary <- bind_rows(rhr_summary_1, rhr_summary_50) %>%
  mutate(
    label_time = paste(label, time_period, sep = "_")
  )

# 然后创建和保存图表
rhr_plots <- plot_rhr_patterns(combined_rhr_summary)
rhr_plots


# Save mean RHR plot
ggsave(
  "perioperative_mean_rhr_patterns.pdf",
  plot = rhr_plots$mean_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Save median RHR plot
ggsave(
  "perioperative_median_rhr_patterns.pdf",
  plot = rhr_plots$median_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Create summary stats object and save it
rhr_detailed_summary <- rhr_plots$summary
save(rhr_detailed_summary, file = "rhr_detailed_summary.rda", compress = "xz")

# Print summary table
print("RHR Statistics Summary by Time Period and Activity Level:")
print(rhr_plots$summary)

#######plot
# Function to calculate daily pattern for steps <= 1
calculate_pattern_1 <- function(data) {
  # Filter sample info for steps <= 1
  filtered_samples <- data@sample_info %>%
    dplyr::filter(label == "<1")
  
  # Get corresponding heart rate values
  heart_rates <- data@expression_data[1, filtered_samples$sample_id] %>%
    as.numeric()
  
  # Create data frame with filtered data
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

# Function to calculate daily pattern for steps <= 50 (including steps <= 1)
calculate_pattern_50 <- function(data) {
  # Filter sample info for steps <= 50 (includes both "<1" and "<50")
  filtered_samples <- data@sample_info %>%
    dplyr::filter(label %in% c("<1", "<50"))
  
  # Get corresponding heart rate values
  heart_rates <- data@expression_data[1, filtered_samples$sample_id] %>%
    as.numeric()
  
  # Create data frame with filtered data
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

# Calculate patterns
daily_rhr_pattern_1 <- calculate_pattern_1(heart_rate_data)
daily_rhr_pattern_50 <- calculate_pattern_50(heart_rate_data)

# Create complete hour sequence
all_hours <- tibble(hour = 0:23)

# Fill in missing hours with NA for both patterns
daily_rhr_pattern_1 <- all_hours %>%
  left_join(daily_rhr_pattern_1, by = "hour") %>%
  mutate(label = "1")

daily_rhr_pattern_50 <- all_hours %>%
  left_join(daily_rhr_pattern_50, by = "hour") %>%
  mutate(label = "50")

# Combine the patterns
daily_rhr_pattern <- bind_rows(daily_rhr_pattern_1, daily_rhr_pattern_50)

# Calculate overall median RHR
overall_median_rhr <- {
  filtered_samples <- heart_rate_data@sample_info %>%
    dplyr::filter(label %in% c("<1", "<50"))
  heart_rates <- heart_rate_data@expression_data[1, filtered_samples$sample_id] %>%
    as.numeric()
  median(heart_rates)
}

# Create the plot
ggplot(daily_rhr_pattern, aes(x = hour, y = median_rhr, color = label, group = label)) +
  # Add points and lines
  geom_line(size = 1) +
  geom_point(size = 2, alpha = 0.6) +
  # Add the median reference line
  geom_hline(
    yintercept = overall_median_rhr,
    color = "red",
    linetype = "dashed"
  ) +
  # Add median RHR label
  annotate(
    "text",
    x = 23,
    y = overall_median_rhr,
    label = "Median RHR",
    hjust = 1,
    vjust = -0.5,
    color = "red"
  ) +
  # Customize scales
  scale_x_continuous(
    breaks = seq(0, 24, 4),
    limits = c(0, 24),
    name = "Time of day"
  ) +
  scale_y_continuous(
    name = "Median RHR",
    limits = c(55, 80)
  ) +
  # Customize colors
  scale_color_manual(
    values = c("1" = "steelblue", "50" = "#FF9999"),
    name = "Max. steps"
  ) +
  # Use theme_bw and customize
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.15),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  # Add title
  labs(title = "Daily RHR Pattern (10-min window)")

# Save the plot
ggsave(
  "daily_rhr_pattern.pdf",
  width = 10,
  height = 6,
  dpi = 300
)

# Print summary statistics
print("\nMeasurements summary by hour:")
daily_rhr_pattern %>%
  group_by(label) %>%
  summarise(
    mean_rhr = mean(median_rhr, na.rm = TRUE),
    min_rhr = min(median_rhr, na.rm = TRUE),
    max_rhr = max(median_rhr, na.rm = TRUE),
    n_hours = sum(!is.na(median_rhr))
  ) %>%
  print()



# Calculate heart rate histogram
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
      fill = "#a6c0d5",  # 浅红色填充
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

p <- plot_heart_rate_histogram(heart_rate_data)
p

ggsave(
  "RHR_histogram.pdf",
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


#####列名修改
create_wide_format_rhr <- function(combined_rhr_summary) {
  # 创建基本数据框，包含所有唯一的subject_id
  subject_ids <- unique(combined_rhr_summary$subject_id)
  result <- data.frame(subject_id = subject_ids)
  
  # 尝试从baseline_info获取手术日期
  # 如果baseline_info存在且可用
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


# 查看combined_rhr_summary的结构，以便调试
str(combined_rhr_summary)

# 创建并保存宽格式RHR数据
time_period_rhr_results <- create_wide_format_rhr(combined_rhr_summary)

# 验证列名格式
cat("新的RHR列名格式示例:\n")
print(names(time_period_rhr_results)[1:10])

# 保存结果
save(time_period_rhr_results, 
     file = "time_period_rhr_results.rda", 
     compress = "xz")
