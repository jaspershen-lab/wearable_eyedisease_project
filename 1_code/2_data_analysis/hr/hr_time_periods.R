library(tidyverse)
library(tidymass)
library(r4projects)
library(moments)  # 用于计算偏度和峰度
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


######## 读取数据
load(
  "3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############创建输出目录
dir.create("3_data_analysis/2_data_analysis/hr/hr_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/hr/hr_time_periods")


###############心率(HR)时间窗口处理
# 处理基线信息
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)


# 提取心率数据的函数
extract_hr_data <- function(heart_rate_data) {
  # 获取心率值
  hr_values <- heart_rate_data@expression_data[1, ] %>%
    as.numeric()
  
  # 创建包含心率值的数据框
  data.frame(
    sample_id = colnames(heart_rate_data@expression_data),
    heart_rate = hr_values,
    timestamp = heart_rate_data@sample_info$measure_time,
    subject_id = heart_rate_data@sample_info$subject_id
  ) %>%
    dplyr::filter(heart_rate >= 30 & heart_rate <= 220) %>%  # 有效心率范围
    arrange(subject_id, timestamp)
}

# 提取心率数据
hr_data <- extract_hr_data(heart_rate_data)

# 计算连续RR间隔的RMSSD (适用于每个受试者的时间窗口)
calculate_rmssd <- function(heart_rates, timestamps) {
  # 确保按时间排序
  if(length(heart_rates) <= 1) {
    return(NA) # 需要至少两个点才能计算RMSSD
  }
  
  # 计算相邻点的差值
  hr_diffs <- diff(heart_rates)
  
  # 返回均方根
  sqrt(mean(hr_diffs^2, na.rm = TRUE))
}

# 心率数据时间窗口处理函数
process_hr_time_periods <- function(hr_data, baseline_info_processed) {
  # 计算基本时间信息
  base_data <- hr_data %>%
    left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
    mutate(
      hours_to_surgery = as.numeric(difftime(timestamp, surgery_time_1, units = "hours"))
    )
  
  # 定义时间窗口统计函数
  calculate_time_period_stats <- function(data, min_hours, max_hours, period_name) {
    period_data <- data %>%
      filter(hours_to_surgery >= min_hours & hours_to_surgery < max_hours)
    
    # 如果没有数据，返回NA行
    if(nrow(period_data) == 0) {
      return(data.frame())
    }
    
    # 按受试者分组计算RMSSD
    rmssd_by_subject <- period_data %>%
      group_by(subject_id) %>%
      arrange(timestamp) %>%
      summarize(
        hr_rmssd = calculate_rmssd(heart_rate, timestamp),
        .groups = "drop"
      )
    
    # 计算其他统计量
    main_stats <- period_data %>%
      group_by(subject_id) %>%
      summarise(
        hr_mean = mean(heart_rate),
        hr_min = min(heart_rate),
        hr_max = max(heart_rate),
        hr_median = median(heart_rate),
        hr_q1 = quantile(heart_rate, 0.25),
        hr_q3 = quantile(heart_rate, 0.75),
        hr_iqr = IQR(heart_rate),
        hr_sd = sd(heart_rate),
        hr_cv = (sd(heart_rate) / mean(heart_rate)) * 100, # CV as percentage
        hr_skew = skewness(heart_rate),  # 偏度
        hr_kurt = kurtosis(heart_rate),  # 峰度
        n_measurements = n(),
        .groups = "drop"
      )
    
    # 合并RMSSD和其他统计量
    main_stats %>%
      left_join(rmssd_by_subject, by = "subject_id") %>%
      mutate(time_period = period_name)
  }
  
  # 处理术前时间窗口
  pre_3d <- calculate_time_period_stats(base_data, -72, 0, "pre_surgery_3d")
  pre_3d_to_7d <- calculate_time_period_stats(base_data, -168, -72, "pre_surgery_3d_to_7d")
  pre_7d_all <- calculate_time_period_stats(base_data, -168, 0, "pre_surgery_7d_all")
  pre_all <- calculate_time_period_stats(base_data, -720, 0, "pre_surgery_all")
  
  # 处理术后时间窗口
  # 术后1-3天
  post_1to3d <- calculate_time_period_stats(base_data, 24, 72, "post_surgery_1to3d")
  
  # 术后4-6天
  post_4to6d <- calculate_time_period_stats(base_data, 72, 144, "post_surgery_4to6d")
  
  # 术后6天 (0-144小时)
  post_6d <- calculate_time_period_stats(base_data, 0, 144, "post_surgery_6d")
  
  # 术后1周 (0-168小时)
  post_1w <- calculate_time_period_stats(base_data, 0, 168, "post_surgery_1w")
  
  # 术后7天
  post_7d <- calculate_time_period_stats(base_data, 0, 168, "post_surgery_7d")
  
  # 术后7-30天
  post_7d_to_30d <- calculate_time_period_stats(base_data, 168, 720, "post_surgery_7d_to_30d")
  
  # 术后23-30天
  post_day23_to_30 <- calculate_time_period_stats(base_data, 552, 720, "post_surgery_day23_to_30")
  
  # 术后27-30天
  post_day27_to_30 <- calculate_time_period_stats(base_data, 648, 720, "post_surgery_day27_to_30")
  
  # 术后>30天
  post_over_30d <- calculate_time_period_stats(base_data, 720, Inf, "post_surgery_over_30d")
  
  # 合并所有时间窗口结果
  bind_rows(
    pre_3d, pre_3d_to_7d, pre_7d_all, pre_all,
    post_1to3d, post_4to6d, post_6d,
    post_1w, post_7d, post_7d_to_30d,
    post_day23_to_30, post_day27_to_30,
    post_over_30d
  )
}

# 处理心率数据
hr_time_periods <- process_hr_time_periods(hr_data, baseline_info_processed)

# 为每个统计量创建宽格式数据的函数
create_wide_format <- function(hr_time_periods) {
  # 要处理的统计量列表
  stat_cols <- c("mean", "min", "max", "median", "q1", "q3", "iqr", "sd", "cv", "skew", "kurt", "rmssd", "n_measurements")
  
  # 处理每个统计量
  stats_list <- list()
  
  for(stat in stat_cols) {
    # 对于n_measurements，直接使用
    if(stat == "n_measurements") {
      col_name <- "n_measurements"
      # 创建正确的列名，保持原样
      new_col_name <- "col_name"
    } else {
      # 对于其他统计量，使用hr_前缀
      col_name <- paste0("hr_", stat)
      # 创建正确的列名，添加hr标识
      new_col_name <- "hr_col_name"
    }
    
    # 处理每个时间窗口的统计量
    if(stat == "n_measurements") {
      # 对于测量数，保持原始列名格式
      stats_list[[stat]] <- hr_time_periods %>%
        mutate(col_name = paste0(time_period, "_", stat)) %>%
        dplyr::select(subject_id, col_name, !!sym(col_name)) %>%
        pivot_wider(
          names_from = col_name,
          values_from = !!sym(col_name)
        )
    } else {
      # 对于其他统计量，使用新的列名格式
      stats_list[[stat]] <- hr_time_periods %>%
        # 创建新的列名格式，添加hr作为统计类型的标识
        mutate(hr_col_name = paste0(time_period, "_hr_", stat)) %>%
        dplyr::select(subject_id, hr_col_name, !!sym(col_name)) %>%
        pivot_wider(
          names_from = hr_col_name,
          values_from = !!sym(col_name)
        )
    }
  }
  
  # 合并所有统计量的数据
  result <- stats_list[[1]]
  for(i in 2:length(stats_list)) {
    result <- result %>%
      left_join(stats_list[[i]], by = "subject_id")
  }
  
  return(result)
}

# 创建包含手术日期的最终结果
create_hr_time_period_results <- function(hr_time_periods, baseline_info_processed) {
  # 创建宽格式数据
  hr_wide <- create_wide_format(hr_time_periods)
  
  # 创建最终结果，包含手术日期
  tibble(subject_id = unique(hr_time_periods$subject_id)) %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = format(surgery_time_1, "%Y-%m-%d")) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    left_join(hr_wide, by = "subject_id") %>%
    # 重新排序列，将subject_id和surgery_date放在开头
    dplyr::select(subject_id, surgery_date, everything())
}

# 创建宽格式结果，使用新的列名格式
hr_time_period_results <- create_hr_time_period_results(hr_time_periods, baseline_info_processed)

# 计算汇总统计数据
calculate_hr_summary_stats <- function(hr_time_periods) {
  hr_time_periods %>%
    group_by(time_period) %>%
    summarise(
      mean_hr = mean(hr_mean, na.rm = TRUE),
      median_hr = median(hr_median, na.rm = TRUE),
      min_hr = min(hr_min, na.rm = TRUE),
      max_hr = max(hr_max, na.rm = TRUE),
      mean_sd = mean(hr_sd, na.rm = TRUE),
      mean_cv = mean(hr_cv, na.rm = TRUE),
      mean_iqr = mean(hr_iqr, na.rm = TRUE),
      mean_skew = mean(hr_skew, na.rm = TRUE),
      mean_kurt = mean(hr_kurt, na.rm = TRUE),
      mean_rmssd = mean(hr_rmssd, na.rm = TRUE),
      n_subjects = n(),
      .groups = "drop"
    )
}

# 计算汇总统计量
hr_time_period_summary_stats <- calculate_hr_summary_stats(hr_time_periods)

# 验证列名格式
cat("新的HR列名格式示例:\n")
print(names(hr_time_period_results)[1:20])

# 检查新添加时间窗口的列数量
cat("\n检查新添加时间窗口的列数:\n")
cat("术后1-3天(post_surgery_1to3d)列数: ", 
    sum(grepl("post_surgery_1to3d", names(hr_time_period_results))), "\n")
cat("术后4-6天(post_surgery_4to6d)列数: ", 
    sum(grepl("post_surgery_4to6d", names(hr_time_period_results))), "\n")
cat("术后6天(post_surgery_6d)列数: ", 
    sum(grepl("post_surgery_6d", names(hr_time_period_results))), "\n")

# 保存结果
save(hr_time_period_results, file = "hr_time_period_results.rda", compress = "xz")
save(hr_time_period_summary_stats, file = "hr_time_period_summary_stats.rda", compress = "xz")
save(hr_time_periods, file = "hr_time_periods_long.rda", compress = "xz")

# 打印最终数据验证
cat("\n最终数据验证:\n")
cat("唯一受试者数量: ", length(unique(hr_time_periods$subject_id)), "\n")
cat("时间窗口数量: ", length(unique(hr_time_periods$time_period)), "\n")
cat("宽格式中的列数: ", ncol(hr_time_period_results), "\n")

###############绘制心率模式图
# 绘图函数

# 按时间窗口绘制平均心率
plot_mean_hr_by_period <- function(summary_stats) {
  # 设置时间窗口顺序
  summary_stats$time_period <- factor(
    summary_stats$time_period,
    levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
               "post_surgery_1to3d", "post_surgery_4to6d", "post_surgery_6d",
               "post_surgery_1w", "post_surgery_7d", "post_surgery_7d_to_30d", 
               "post_surgery_day23_to_30", "post_surgery_day27_to_30", "post_surgery_over_30d"),
    labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
               "Post 1-3d", "Post 4-6d", "Post 6d",
               "Post 1w", "Post 7d", "Post 7d to 30d", 
               "Post day 23-30", "Post day 27-30", "Post >30d")
  )
  
  # 计算总体中位数
  overall_median <- median(summary_stats$median_hr, na.rm = TRUE)
  
  # 创建图表
  ggplot(summary_stats, aes(x = time_period, y = mean_hr)) +
    geom_bar(stat = "identity", fill = "#a6c0d5", alpha = 0.7, width = 0.6) +
    geom_point(size = 3, color = "#16165d") +
    geom_hline(yintercept = overall_median, linetype = "dashed", 
               color = "red", alpha = 0.7) +
    annotate("text", x = 6, y = overall_median, 
             label = "Overall Median HR", hjust = 1, vjust = -0.5, 
             color = "red") +
    # 添加误差条，使用mean_sd
    geom_errorbar(aes(ymin = mean_hr - mean_sd, ymax = mean_hr + mean_sd),
                  width = 0.2, alpha = 0.5) +
    labs(title = "Heart Rate by Time Period",
         x = "Time Period",
         y = "Mean Heart Rate (BPM)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# 心率变异性指标图 (RMSSD)
plot_hr_variability <- function(summary_stats) {
  # 设置时间窗口顺序
  summary_stats$time_period <- factor(
    summary_stats$time_period,
    levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
               "post_surgery_1to3d", "post_surgery_4to6d", "post_surgery_6d",
               "post_surgery_1w", "post_surgery_7d", "post_surgery_7d_to_30d", 
               "post_surgery_day23_to_30", "post_surgery_day27_to_30", "post_surgery_over_30d"),
    labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
               "Post 1-3d", "Post 4-6d", "Post 6d",
               "Post 1w", "Post 7d", "Post 7d to 30d", 
               "Post day 23-30", "Post day 27-30", "Post >30d")
  )
  
  # 创建图表
  ggplot(summary_stats, aes(x = time_period, y = mean_rmssd)) +
    geom_bar(stat = "identity", fill = "#d5a6a6", alpha = 0.7, width = 0.6) +
    geom_point(size = 3, color = "#5d1616") +
    labs(title = "Heart Rate Variability (RMSSD) by Time Period",
         x = "Time Period",
         y = "Mean RMSSD") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# 统计比较图
plot_statistics_comparison <- function(summary_stats) {
  # 设置时间窗口顺序
  summary_stats$time_period <- factor(
    summary_stats$time_period,
    levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
               "post_surgery_1to3d", "post_surgery_4to6d", "post_surgery_6d",
               "post_surgery_1w", "post_surgery_7d", "post_surgery_7d_to_30d", 
               "post_surgery_day23_to_30", "post_surgery_day27_to_30", "post_surgery_over_30d"),
    labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
               "Post 1-3d", "Post 4-6d", "Post 6d",
               "Post 1w", "Post 7d", "Post 7d to 30d", 
               "Post day 23-30", "Post day 27-30", "Post >30d")
  )
  
  # 创建长格式数据用于绘制多个统计量
  long_stats <- summary_stats %>%
    dplyr::select(time_period, mean_hr, median_hr, mean_cv, mean_iqr) %>%
    pivot_longer(
      cols = c(mean_hr, median_hr, mean_cv, mean_iqr),
      names_to = "statistic",
      values_to = "value"
    ) %>%
    mutate(
      statistic = factor(
        statistic,
        levels = c("mean_hr", "median_hr", "mean_cv", "mean_iqr"),
        labels = c("Mean", "Median", "CV (%)", "IQR")
      )
    )
  
  # 创建分面图
  ggplot(long_stats, aes(x = time_period, y = value, fill = statistic)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    facet_wrap(~ statistic, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Heart Rate Statistics by Time Period",
         x = "Time Period",
         y = "Value") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# 创建并保存图表
p_mean <- plot_mean_hr_by_period(hr_time_period_summary_stats)
p_variability <- plot_hr_variability(hr_time_period_summary_stats)
p_stats <- plot_statistics_comparison(hr_time_period_summary_stats)

# 显示图表
print(p_mean)
print(p_variability)
print(p_stats)

# 保存图表
ggsave("hr_by_time_period.pdf", p_mean, width = 10, height = 6, dpi = 300)
ggsave("hr_variability_by_time_period.pdf", p_variability, width = 10, height = 6, dpi = 300)
ggsave("hr_statistics_comparison.pdf", p_stats, width = 12, height = 8, dpi = 300)

# 打印汇总统计数据
print("按时间窗口的心率汇总统计:")
print(hr_time_period_summary_stats)
