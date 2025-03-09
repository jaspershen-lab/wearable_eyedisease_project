library(tidyverse)
library(tidymass)
library(r4projects)
library(moments)  # 用于计算偏度和峰度
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


######## Read data
load(
  "3_data_analysis/1_data_preparation/wearable_data/2_blood_oxygen/blood_oxygen_data.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/bo/bo_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/bo/bo_time_periods")


###############BO by time periods
# Process baseline info
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)


# Function to extract blood oxygen data
extract_bo_data <- function(blood_oxygen_data) {
  # Get blood oxygen values
  bo_values <- blood_oxygen_data@expression_data[1, ] %>%
    as.numeric()
  
  # Create data frame with blood oxygen values
  data.frame(
    sample_id = colnames(blood_oxygen_data@expression_data),
    blood_oxygen = bo_values,
    timestamp = blood_oxygen_data@sample_info$measure_time,
    subject_id = blood_oxygen_data@sample_info$subject_id
  ) %>%
    dplyr::filter(blood_oxygen >= 70 & blood_oxygen <= 100) %>%  # Valid blood oxygen range
    arrange(subject_id, timestamp)
}

# Extract blood oxygen data
bo_data <- extract_bo_data(blood_oxygen_data)

# 修改后的blood oxygen数据时间窗口处理函数
process_bo_time_periods <- function(bo_data, baseline_info_processed) {
  # Calculate basic time information
  base_data <- bo_data %>%
    left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
    mutate(
      hours_to_surgery = as.numeric(difftime(surgery_time_1, timestamp, units = "hours"))
    )
  
  # Define time periods function
  calculate_time_period_stats <- function(data, min_hours, max_hours, period_name) {
    data %>%
      filter(hours_to_surgery >= min_hours & hours_to_surgery < max_hours) %>%
      group_by(subject_id) %>%
      summarise(
        bo_mean = mean(blood_oxygen),
        bo_min = min(blood_oxygen),
        bo_max = max(blood_oxygen),
        bo_median = median(blood_oxygen),
        bo_q1 = quantile(blood_oxygen, 0.25),
        bo_q3 = quantile(blood_oxygen, 0.75),
        bo_iqr = IQR(blood_oxygen),
        bo_sd = sd(blood_oxygen),
        bo_cv = (sd(blood_oxygen) / mean(blood_oxygen)) * 100, # CV as percentage
        bo_skew = skewness(blood_oxygen),  # 偏度
        bo_kurt = kurtosis(blood_oxygen),  # 峰度
        n_measurements = n(),
        .groups = "drop"
      ) %>%
      mutate(time_period = period_name)
  }
  
  # 原有的时间窗口
  pre_3d <- calculate_time_period_stats(base_data, 0, 72, "pre_surgery_3d")
  pre_3d_to_7d <- calculate_time_period_stats(base_data, 72, 168, "pre_surgery_3d_to_7d")
  pre_7d_all <- calculate_time_period_stats(base_data, 0, 168, "pre_surgery_7d_all")
  pre_all <- calculate_time_period_stats(base_data, 0, Inf, "pre_surgery_all")
  post_7d <- calculate_time_period_stats(base_data, -168, 0, "post_surgery_7d")
  post_7d_to_30d <- calculate_time_period_stats(base_data, -720, -168, "post_surgery_7d_to_30d")
  post_over_30d <- calculate_time_period_stats(base_data, -Inf, -720, "post_surgery_over_30d")
  
  # 新增的时间窗口：术后前一周 (Post 1w)
  post_1w <- calculate_time_period_stats(base_data, -168, 0, "post_surgery_1w")
  
  # 新增的时间窗口：术后第30天前3天 (day 27-30, hours -720 to -648)
  post_day27_to_30 <- calculate_time_period_stats(base_data, -720, -648, "post_surgery_day27_to_30")
  
  # 新增的时间窗口：术后第30天前一周 (day 23-30, hours -720 to -552)
  post_day23_to_30 <- calculate_time_period_stats(base_data, -720, -552, "post_surgery_day23_to_30")
  
  # Combine all time periods，包括新增的时间窗口
  bind_rows(
    pre_3d, pre_3d_to_7d, pre_7d_all, pre_all, 
    post_1w, post_7d, post_7d_to_30d, 
    post_day23_to_30, post_day27_to_30, 
    post_over_30d
  )
}

# Process blood oxygen data
bo_time_periods <- process_bo_time_periods(bo_data, baseline_info_processed)

# Create wide format for each statistic
create_wide_format <- function(bo_time_periods) {
  # List of statistics to process
  stat_cols <- c("mean", "min", "max", "median", "q1", "q3", "iqr", "sd", "cv", "skew", "kurt", "n_measurements")
  
  # Process each statistic
  stats_list <- list()
  
  for(stat in stat_cols) {
    col_name <- paste0("bo_", stat)
    # For n_measurements, use it directly
    if(stat == "n_measurements") {
      col_name <- "n_measurements"
    }
    
    stats_list[[stat]] <- bo_time_periods %>%
      mutate(col_name = paste0(time_period, "_", stat)) %>%
      dplyr::select(subject_id, col_name, !!sym(col_name)) %>%
      pivot_wider(
        names_from = col_name,
        values_from = !!sym(col_name)
      )
  }
  
  # Combine all statistics
  result <- stats_list[[1]]
  for(i in 2:length(stats_list)) {
    result <- result %>%
      left_join(stats_list[[i]], by = "subject_id")
  }
  
  return(result)
}

# Create wide format result
bo_wide <- create_wide_format(bo_time_periods)

# Create final result with surgery date
bo_time_period_results <- tibble(subject_id = unique(bo_time_periods$subject_id)) %>%
  left_join(
    baseline_info_processed %>%
      mutate(surgery_date = format(surgery_time_1, "%Y-%m-%d")) %>%
      dplyr::select(ID, surgery_date),
    by = c("subject_id" = "ID")
  ) %>%
  left_join(bo_wide, by = "subject_id") %>%
  # Reorder columns to put subject_id and surgery_date at the beginning
  dplyr::select(subject_id, surgery_date, everything())

# Calculate summary statistics
calculate_bo_summary_stats <- function(bo_time_periods) {
  bo_time_periods %>%
    group_by(time_period) %>%
    summarise(
      mean_bo = mean(bo_mean, na.rm = TRUE),
      median_bo = median(bo_median, na.rm = TRUE),
      min_bo = min(bo_min, na.rm = TRUE),
      max_bo = max(bo_max, na.rm = TRUE),
      mean_sd = mean(bo_sd, na.rm = TRUE),
      mean_cv = mean(bo_cv, na.rm = TRUE),
      mean_iqr = mean(bo_iqr, na.rm = TRUE),
      mean_skew = mean(bo_skew, na.rm = TRUE),
      mean_kurt = mean(bo_kurt, na.rm = TRUE),
      n_subjects = n(),
      .groups = "drop"
    )
}

# Calculate summary statistics
bo_time_period_summary_stats <- calculate_bo_summary_stats(bo_time_periods)

# Save results
save(bo_time_period_results, file = "bo_time_period_results.rda", compress = "xz")
save(bo_time_period_summary_stats, file = "bo_time_period_summary_stats.rda", compress = "xz")
save(bo_time_periods, file = "bo_time_periods_long.rda", compress = "xz")

# Print summary validation
cat("\nFinal data validation:\n")
cat("Number of unique subjects: ", length(unique(bo_time_periods$subject_id)), "\n")
cat("Number of time periods: ", length(unique(bo_time_periods$time_period)), "\n")
cat("Number of columns in wide format: ", ncol(bo_time_period_results), "\n")

# 更新平均血氧绘图函数
plot_mean_bo_by_period <- function(summary_stats) {
  # 更新时间窗口顺序和标签，包括新增的时间窗口
  summary_stats$time_period <- factor(
    summary_stats$time_period,
    levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
               "post_surgery_1w", "post_surgery_7d", "post_surgery_7d_to_30d", 
               "post_surgery_day23_to_30", "post_surgery_day27_to_30", "post_surgery_over_30d"),
    labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
               "Post 1w", "Post 7d", "Post 7d to 30d", 
               "Post day 23-30", "Post day 27-30", "Post >30d")
  )
  
  # Calculate overall median
  overall_median <- median(summary_stats$median_bo, na.rm = TRUE)
  
  # Create plot
  ggplot(summary_stats, aes(x = time_period, y = mean_bo)) +
    geom_bar(stat = "identity", fill = "#a6c0d5", alpha = 0.7, width = 0.6) +
    geom_point(size = 3, color = "#16165d") +
    geom_hline(yintercept = overall_median, linetype = "dashed", 
               color = "red", alpha = 0.7) +
    annotate("text", x = 6, y = overall_median, 
             label = "Overall Median SpO2", hjust = 1, vjust = -0.5, 
             color = "red") +
    # Add error bars using mean_sd
    geom_errorbar(aes(ymin = mean_bo - mean_sd, ymax = mean_bo + mean_sd),
                  width = 0.2, alpha = 0.5) +
    scale_y_continuous(limits = c(90, 100), breaks = seq(90, 100, 2)) +
    labs(title = "Blood Oxygen Saturation (SpO2) by Time Period",
         x = "Time Period",
         y = "Mean SpO2 (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# 更新统计对比绘图函数
plot_statistics_comparison <- function(summary_stats) {
  # 更新时间窗口顺序和标签，包括新增的时间窗口
  summary_stats$time_period <- factor(
    summary_stats$time_period,
    levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
               "post_surgery_1w", "post_surgery_7d", "post_surgery_7d_to_30d", 
               "post_surgery_day23_to_30", "post_surgery_day27_to_30", "post_surgery_over_30d"),
    labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
               "Post 1w", "Post 7d", "Post 7d to 30d", 
               "Post day 23-30", "Post day 27-30", "Post >30d")
  )
  
  # Create long format for plotting multiple statistics
  long_stats <- summary_stats %>%
    dplyr::select(time_period, mean_bo, median_bo, mean_cv, mean_iqr) %>%
    pivot_longer(
      cols = c(mean_bo, median_bo, mean_cv, mean_iqr),
      names_to = "statistic",
      values_to = "value"
    ) %>%
    mutate(
      statistic = factor(
        statistic,
        levels = c("mean_bo", "median_bo", "mean_cv", "mean_iqr"),
        labels = c("Mean", "Median", "CV (%)", "IQR")
      )
    )
  
  # Create faceted plot
  ggplot(long_stats, aes(x = time_period, y = value, fill = statistic)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    facet_wrap(~ statistic, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Blood Oxygen Statistics by Time Period",
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

# Create and save plots
p_mean <- plot_mean_bo_by_period(bo_time_period_summary_stats)
p_stats <- plot_statistics_comparison(bo_time_period_summary_stats)
p_mean
p_stats

# Save plots
ggsave("bo_by_time_period.pdf", p_mean, width = 10, height = 6, dpi = 300)
ggsave("bo_statistics_comparison.pdf", p_stats, width = 12, height = 8, dpi = 300)

# Print summary stats
print("Blood Oxygen Summary Statistics by Time Period:")
print(bo_time_period_summary_stats)



# 修改后的宽格式BO数据创建函数
# 这个函数会直接生成包含"bo"的正确列名格式

create_wide_format <- function(bo_time_periods) {
  # List of statistics to process
  stat_cols <- c("mean", "min", "max", "median", "q1", "q3", "iqr", "sd", "cv", "skew", "kurt", "n_measurements")
  
  # Process each statistic
  stats_list <- list()
  
  for(stat in stat_cols) {
    # 对于n_measurements，直接使用
    if(stat == "n_measurements") {
      col_name <- "n_measurements"
      # 创建正确的列名，保持原样
      new_col_name <- "col_name"
    } else {
      # 对于其他统计量，使用bo_前缀
      col_name <- paste0("bo_", stat)
      # 创建正确的列名，添加bo标识
      new_col_name <- "bo_col_name"
    }
    
    # 处理每个时间窗口的统计量
    if(stat == "n_measurements") {
      # 对于测量数，保持原始列名格式
      stats_list[[stat]] <- bo_time_periods %>%
        mutate(col_name = paste0(time_period, "_", stat)) %>%
        dplyr::select(subject_id, col_name, !!sym(col_name)) %>%
        pivot_wider(
          names_from = col_name,
          values_from = !!sym(col_name)
        )
    } else {
      # 对于其他统计量，使用新的列名格式
      stats_list[[stat]] <- bo_time_periods %>%
        # 创建新的列名格式，添加bo作为统计类型的标识
        mutate(bo_col_name = paste0(time_period, "_bo_", stat)) %>%
        dplyr::select(subject_id, bo_col_name, !!sym(col_name)) %>%
        pivot_wider(
          names_from = bo_col_name,
          values_from = !!sym(col_name)
        )
    }
  }
  
  # Combine all statistics
  result <- stats_list[[1]]
  for(i in 2:length(stats_list)) {
    result <- result %>%
      left_join(stats_list[[i]], by = "subject_id")
  }
  
  return(result)
}

# 修改后的BO时间窗口结果生成函数
create_bo_time_period_results <- function(bo_time_periods, baseline_info_processed) {
  # 创建宽格式数据
  bo_wide <- create_wide_format(bo_time_periods)
  
  # 创建最终结果，包含手术日期
  tibble(subject_id = unique(bo_time_periods$subject_id)) %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = format(surgery_time_1, "%Y-%m-%d")) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    left_join(bo_wide, by = "subject_id") %>%
    # 重新排序列，将subject_id和surgery_date放在开头
    dplyr::select(subject_id, surgery_date, everything())
}

# 使用该函数的示例代码（不要执行，仅作为参考）

# # 处理blood oxygen数据
bo_time_periods <- process_bo_time_periods(bo_data, baseline_info_processed)

# # 创建宽格式结果，使用新的列名格式
bo_time_period_results <- create_bo_time_period_results(bo_time_periods, baseline_info_processed)

# # 验证列名格式
 cat("新的BO列名格式示例:\n")
 print(names(bo_time_period_results)[1:10])

# # 保存结果
save(bo_time_period_results, file = "bo_time_period_results.rda", compress = "xz")
