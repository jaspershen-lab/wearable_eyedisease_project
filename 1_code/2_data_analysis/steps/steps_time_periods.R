library(tidyverse)
library(tidymass)
library(r4projects)
library(moments)  # 用于计算偏度和峰度（如果需要）
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


######## Read data
load(
  "3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details/daily_workout_details_data.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/steps/steps_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/steps/steps_time_periods")


###############steps by time periods
# Process baseline info
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)


# Function to extract steps data
extract_steps_data <- function(daily_workout_details_data) {
  # Get steps data from expression data
  steps_values <- daily_workout_details_data@expression_data["steps", ] %>%
    as.numeric()
  
  # Create data frame with steps values
  data.frame(
    sample_id = colnames(daily_workout_details_data@expression_data),
    steps = steps_values,
    timestamp = daily_workout_details_data@sample_info$measure_time,
    subject_id = daily_workout_details_data@sample_info$subject_id
  ) %>%
    dplyr::filter(steps >= 0) %>% # Remove negative values if any
    arrange(subject_id, timestamp)
}

# Extract steps data
steps_data <- extract_steps_data(daily_workout_details_data)

# Process steps data with time periods
process_steps_time_periods <- function(steps_data, baseline_info_processed) {
  # Calculate basic time information
  base_data <- steps_data %>%
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
        steps_total = sum(steps, na.rm = TRUE),
        steps_mean = mean(steps, na.rm = TRUE),
        steps_median = median(steps, na.rm = TRUE),
        steps_max = max(steps, na.rm = TRUE),
        steps_min = min(steps, na.rm = TRUE),
        steps_sd = sd(steps, na.rm = TRUE),
        steps_cv = (sd(steps, na.rm = TRUE) / mean(steps, na.rm = TRUE)) * 100, # CV as percentage
        steps_skew = skewness(steps, na.rm = TRUE),  # 偏度
        steps_kurt = kurtosis(steps, na.rm = TRUE),  # 峰度
        n_measurements = n(),
        .groups = "drop"
      ) %>%
      mutate(time_period = period_name)
  }
  
  # Calculate statistics for each time period
  pre_3d <- calculate_time_period_stats(base_data, 0, 72, "pre_surgery_3d")
  pre_3d_to_7d <- calculate_time_period_stats(base_data, 72, 168, "pre_surgery_3d_to_7d")
  pre_7d_all <- calculate_time_period_stats(base_data, 0, 168, "pre_surgery_7d_all")
  pre_all <- calculate_time_period_stats(base_data, 0, Inf, "pre_surgery_all")
  post_7d <- calculate_time_period_stats(base_data, -168, 0, "post_surgery_7d")
  post_7d_to_30d <- calculate_time_period_stats(base_data, -720, -168, "post_surgery_7d_to_30d") # 修改变量名为更准确的描述
  post_over_30d <- calculate_time_period_stats(base_data, -Inf, -720, "post_surgery_over_30d")
  
  # Combine all time periods
  bind_rows(pre_3d, pre_3d_to_7d, pre_7d_all, pre_all, post_7d, post_7d_to_30d, post_over_30d) # 更新变量名
}

# Process steps data
steps_time_periods <- process_steps_time_periods(steps_data, baseline_info_processed)

# Create wide format for each statistic
create_wide_format <- function(steps_time_periods) {
  # List of statistics to process
  stat_cols <- c("total", "mean", "median", "max", "min", "sd", "cv", "skew", "kurt", "n_measurements")
  
  # Process each statistic
  stats_list <- list()
  
  for(stat in stat_cols) {
    col_name <- paste0("steps_", stat)
    # For n_measurements, use it directly
    if(stat == "n_measurements") {
      col_name <- "n_measurements"
    }
    
    stats_list[[stat]] <- steps_time_periods %>%
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
steps_wide <- create_wide_format(steps_time_periods)

# Create final result with surgery date
steps_time_period_results <- tibble(subject_id = unique(steps_time_periods$subject_id)) %>%
  left_join(
    baseline_info_processed %>%
      mutate(surgery_date = format(surgery_time_1, "%Y-%m-%d")) %>%
      dplyr::select(ID, surgery_date),
    by = c("subject_id" = "ID")
  ) %>%
  left_join(steps_wide, by = "subject_id") %>%
  # Reorder columns to put subject_id and surgery_date at the beginning
  dplyr::select(subject_id, surgery_date, everything())

# Calculate summary statistics
calculate_steps_summary_stats <- function(steps_time_periods) {
  steps_time_periods %>%
    group_by(time_period) %>%
    summarise(
      mean_steps_total = mean(steps_total, na.rm = TRUE),
      mean_steps_mean = mean(steps_mean, na.rm = TRUE),
      mean_steps_median = median(steps_median, na.rm = TRUE),
      mean_steps_max = mean(steps_max, na.rm = TRUE),
      mean_steps_sd = mean(steps_sd, na.rm = TRUE),
      mean_steps_cv = mean(steps_cv, na.rm = TRUE),
      mean_steps_skew = mean(steps_skew, na.rm = TRUE),
      mean_steps_kurt = mean(steps_kurt, na.rm = TRUE),
      mean_n_measurements = mean(n_measurements, na.rm = TRUE),
      n_subjects = n(),
      .groups = "drop"
    )
}

# Calculate summary statistics
steps_time_period_summary_stats <- calculate_steps_summary_stats(steps_time_periods)

# Save results
save(steps_time_period_results, file = "steps_time_period_results.rda", compress = "xz")
save(steps_time_period_summary_stats, file = "steps_time_period_summary_stats.rda", compress = "xz")
save(steps_time_periods, file = "steps_time_periods_long.rda", compress = "xz")

# Print summary validation
cat("\nFinal data validation:\n")
cat("Number of unique subjects: ", length(unique(steps_time_periods$subject_id)), "\n")
cat("Number of time periods: ", length(unique(steps_time_periods$time_period)), "\n")
cat("Number of columns in wide format: ", ncol(steps_time_period_results), "\n")

# Create plots

# Plot average daily steps by time period
plot_average_steps_by_period <- function(summary_stats) {
  # Set time period order
  summary_stats$time_period <- factor(
    summary_stats$time_period,
    levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
               "post_surgery_7d", "post_surgery_7d_to_30d", "post_surgery_over_30d"), # 修正标签
    labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
               "Post 7d", "Post 7d to 30d", "Post >30d") # 更新标签文字
  )
  
  # Print data to debug
  print("Data used for plotting:")
  print(summary_stats[, c("time_period", "mean_steps_mean", "mean_steps_sd")])
  
  # Create plot with improvements
  ggplot(summary_stats, aes(x = time_period, y = mean_steps_mean)) +
    geom_bar(stat = "identity", fill = "#4CAF50", alpha = 0.7, width = 0.6) +
    # Prevent negative error bars and add text labels
    geom_errorbar(aes(ymin = pmax(0, mean_steps_mean - mean_steps_sd), 
                      ymax = mean_steps_mean + mean_steps_sd),
                  width = 0.2, alpha = 0.7) +
    # Add text labels on top of bars
    geom_text(aes(label = sprintf("%.0f", mean_steps_mean)), 
              vjust = -0.5, size = 3.5) +
    # 10,000 steps reference line
    geom_hline(yintercept = 10000, linetype = "dashed", 
               color = "red", alpha = 0.7) +
    annotate("text", x = 6, y = 10000, 
             label = "10,000 steps goal", hjust = 1, vjust = -0.5, 
             color = "red") +
    # Dynamic y-axis scaling
    scale_y_continuous(
      limits = c(0, max(max(summary_stats$mean_steps_mean + summary_stats$mean_steps_sd, na.rm = TRUE), 10000) * 1.1),
      breaks = seq(0, ceiling(max(c(summary_stats$mean_steps_mean, 10000), na.rm = TRUE)/1000)*1000, by = 2000)
    ) +
    labs(title = "Average Daily Steps by Time Period",
         x = "Time Period",
         y = "Average Steps per Day") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# Plot total steps by time period
plot_total_steps_by_period <- function(summary_stats) {
  # Set time period order
  summary_stats$time_period <- factor(
    summary_stats$time_period,
    levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
               "post_surgery_7d", "post_surgery_7d_to_30d", "post_surgery_over_30d"), # 修正标签
    labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
               "Post 7d", "Post 7d to 30d", "Post >30d") # 更新标签文字
  )
  
  # Create plot
  ggplot(summary_stats, aes(x = time_period, y = mean_steps_total)) +
    geom_bar(stat = "identity", fill = "#2196F3", alpha = 0.7, width = 0.6) +
    geom_text(aes(label = sprintf("%.0f", mean_steps_total)), 
              vjust = -0.5, size = 3.5) +
    labs(title = "Total Steps by Time Period",
         x = "Time Period",
         y = "Total Steps") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# Plot statistics comparison
plot_statistics_comparison <- function(summary_stats) {
  # Set time period order
  summary_stats$time_period <- factor(
    summary_stats$time_period,
    levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
               "post_surgery_7d", "post_surgery_7d_to_30d", "post_surgery_over_30d"), # 修正标签
    labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
               "Post 7d", "Post 7d to 30d", "Post >30d") # 更新标签文字
  )
  
  # Create long format for plotting multiple statistics
  long_stats <- summary_stats %>%
    dplyr::select(time_period, mean_steps_mean, mean_steps_median, mean_steps_max) %>%
    pivot_longer(
      cols = c(mean_steps_mean, mean_steps_median, mean_steps_max),
      names_to = "statistic",
      values_to = "value"
    ) %>%
    mutate(
      statistic = factor(
        statistic,
        levels = c("mean_steps_mean", "mean_steps_median", "mean_steps_max"),
        labels = c("Mean", "Median", "Maximum")
      )
    )
  
  # Create faceted plot
  ggplot(long_stats, aes(x = time_period, y = value, fill = statistic)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    facet_wrap(~ statistic, scales = "free_y", ncol = 3) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Steps Statistics by Time Period",
         x = "Time Period",
         y = "Steps") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# Create and save plots
p_avg <- plot_average_steps_by_period(steps_time_period_summary_stats)
p_total <- plot_total_steps_by_period(steps_time_period_summary_stats)
p_stats <- plot_statistics_comparison(steps_time_period_summary_stats)
p_avg 
p_total
p_stats

# Save plots
ggsave("average_steps_by_time_period.pdf", p_avg, width = 10, height = 6, dpi = 300)
ggsave("total_steps_by_time_period.pdf", p_total, width = 10, height = 6, dpi = 300)
ggsave("steps_statistics_comparison.pdf", p_stats, width = 12, height = 6, dpi = 300)

# Print summary stats
print("Steps Summary Statistics by Time Period:")
print(steps_time_period_summary_stats)
