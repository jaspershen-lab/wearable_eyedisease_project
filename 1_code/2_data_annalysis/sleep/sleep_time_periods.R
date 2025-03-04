library(tidyverse)
library(tidymass)
library(r4projects)
library(moments)  # 用于计算偏度和峰度
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


######## Read data
load(
  "3_data_analysis/1_data_preparation/wearable_data/7_sleep/sleep_data.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/sleep/sleep_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/sleep/sleep_time_periods")


###############sleep by time periods
# Process baseline info
baseline_info_processed <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
  ) %>%
  dplyr::select(ID, surgery_time_1)


# Function to extract sleep data
extract_sleep_data <- function(sleep_data) {
  # Get all sleep metrics
  sleep_metrics <- sleep_data@expression_data %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
  
  # Add sample information
  sleep_metrics <- sleep_metrics %>%
    mutate(
      sample_id = rownames(.),
      timestamp = sleep_data@sample_info$measure_time,
      subject_id = sleep_data@sample_info$subject_id,
      sleep_start = sleep_data@sample_info$sleep_start_time,
      sleep_end = sleep_data@sample_info$sleep_end_time
    )
  
  # Rename sleep metrics columns
  colnames(sleep_metrics)[1:7] <- c(
    "light_sleep",
    "deep_sleep",
    "dream_sleep",
    "awake",
    "total_sleep",
    "daytime_sleep",
    "sleep_score"
  )
  
  # Add sleep duration calculation
  sleep_metrics %>%
    mutate(
      sleep_duration = as.numeric(difftime(sleep_end, sleep_start, units = "mins"))
    ) %>%
    arrange(subject_id, sleep_start)
}

# Extract sleep data
sleep_data_processed <- extract_sleep_data(sleep_data)

# Process sleep data with time periods
process_sleep_time_periods <- function(sleep_data_processed, baseline_info_processed) {
  # Calculate basic time information
  base_data <- sleep_data_processed %>%
    left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
    mutate(
      # Use sleep_start for calculating time to surgery
      hours_to_surgery = as.numeric(difftime(surgery_time_1, sleep_start, units = "hours"))
    )
  
  # Define time periods function
  calculate_time_period_stats <- function(data, min_hours, max_hours, period_name) {
    data %>%
      filter(hours_to_surgery >= min_hours & hours_to_surgery < max_hours) %>%
      group_by(subject_id) %>%
      summarise(
        light_sleep_total = sum(light_sleep, na.rm = TRUE),
        deep_sleep_total = sum(deep_sleep, na.rm = TRUE),
        dream_sleep_total = sum(dream_sleep, na.rm = TRUE),
        awake_total = sum(awake, na.rm = TRUE),
        total_sleep_total = sum(total_sleep, na.rm = TRUE),
        daytime_sleep_total = sum(daytime_sleep, na.rm = TRUE),
        sleep_score_mean = mean(sleep_score, na.rm = TRUE),
        sleep_duration_mean = mean(sleep_duration, na.rm = TRUE),
        sleep_duration_median = median(sleep_duration, na.rm = TRUE),
        sleep_duration_min = min(sleep_duration, na.rm = TRUE),
        sleep_duration_max = max(sleep_duration, na.rm = TRUE),
        sleep_duration_sd = sd(sleep_duration, na.rm = TRUE),
        sleep_duration_cv = (sd(sleep_duration, na.rm = TRUE) / mean(sleep_duration, na.rm = TRUE)) * 100,
        n_sleep_sessions = n(),
        .groups = "drop"
      ) %>%
      mutate(time_period = period_name)
  }
  
  # 按照步数分析的时间窗口修改 - 统一时间窗口
  pre_3d <- calculate_time_period_stats(base_data, 0, 72, "pre_surgery_3d")
  pre_3d_to_7d <- calculate_time_period_stats(base_data, 72, 168, "pre_surgery_3d_to_7d")
  pre_7d_all <- calculate_time_period_stats(base_data, 0, 168, "pre_surgery_7d_all")
  pre_all <- calculate_time_period_stats(base_data, 0, Inf, "pre_surgery_all")
  post_7d <- calculate_time_period_stats(base_data, -168, 0, "post_surgery_7d")
  post_7d_to_30d <- calculate_time_period_stats(base_data, -720, -168, "post_surgery_7d_to_30d")
  post_over_30d <- calculate_time_period_stats(base_data, -Inf, -720, "post_surgery_over_30d")
  
  # Combine all time periods
  bind_rows(pre_3d, pre_3d_to_7d, pre_7d_all, pre_all, post_7d, post_7d_to_30d, post_over_30d)
}

# Process sleep data by time periods
sleep_time_periods <- process_sleep_time_periods(sleep_data_processed, baseline_info_processed)

# Create wide format for each statistic
create_wide_format <- function(sleep_time_periods) {
  # List of metrics to process
  metrics <- c(
    "light_sleep_total", "deep_sleep_total", "dream_sleep_total", "awake_total", 
    "total_sleep_total", "daytime_sleep_total", "sleep_score_mean", 
    "sleep_duration_mean", "sleep_duration_median", "sleep_duration_min", 
    "sleep_duration_max", "sleep_duration_sd", "sleep_duration_cv", "n_sleep_sessions"
  )
  
  # Process each metric
  stats_list <- list()
  
  for(metric in metrics) {
    stats_list[[metric]] <- sleep_time_periods %>%
      mutate(col_name = paste0(time_period, "_", metric)) %>%
      dplyr::select(subject_id, col_name, !!sym(metric)) %>%
      pivot_wider(
        names_from = col_name,
        values_from = !!sym(metric)
      )
  }
  
  # Combine all metrics
  result <- stats_list[[1]]
  for(i in 2:length(stats_list)) {
    result <- result %>%
      left_join(stats_list[[i]], by = "subject_id")
  }
  
  return(result)
}

# Create wide format result
sleep_wide <- create_wide_format(sleep_time_periods)

# Create final result with surgery date
sleep_time_period_results <- tibble(subject_id = unique(sleep_time_periods$subject_id)) %>%
  left_join(
    baseline_info_processed %>%
      mutate(surgery_date = format(surgery_time_1, "%Y-%m-%d")) %>%
      dplyr::select(ID, surgery_date),
    by = c("subject_id" = "ID")
  ) %>%
  left_join(sleep_wide, by = "subject_id") %>%
  # Reorder columns to put subject_id and surgery_date at the beginning
  dplyr::select(subject_id, surgery_date, everything())

# Calculate summary statistics
calculate_sleep_summary_stats <- function(sleep_time_periods) {
  sleep_time_periods %>%
    group_by(time_period) %>%
    summarise(
      mean_light_sleep = mean(light_sleep_total, na.rm = TRUE),
      mean_deep_sleep = mean(deep_sleep_total, na.rm = TRUE),
      mean_dream_sleep = mean(dream_sleep_total, na.rm = TRUE),
      mean_total_sleep = mean(total_sleep_total, na.rm = TRUE),
      mean_sleep_score = mean(sleep_score_mean, na.rm = TRUE),
      mean_sleep_duration = mean(sleep_duration_mean, na.rm = TRUE),
      median_sleep_duration = median(sleep_duration_median, na.rm = TRUE),
      mean_sleep_duration_cv = mean(sleep_duration_cv, na.rm = TRUE),
      mean_n_sleep_sessions = mean(n_sleep_sessions, na.rm = TRUE),
      n_subjects = n(),
      .groups = "drop"
    )
}

# Calculate summary statistics
sleep_time_period_summary_stats <- calculate_sleep_summary_stats(sleep_time_periods)

# Save results
save(sleep_time_period_results, file = "sleep_time_period_results.rda", compress = "xz")
save(sleep_time_period_summary_stats, file = "sleep_time_period_summary_stats.rda", compress = "xz")
save(sleep_time_periods, file = "sleep_time_periods_long.rda", compress = "xz")

# Print summary validation
cat("\nFinal data validation:\n")
cat("Number of unique subjects: ", length(unique(sleep_time_periods$subject_id)), "\n")
cat("Number of time periods: ", length(unique(sleep_time_periods$time_period)), "\n")
cat("Number of columns in wide format: ", ncol(sleep_time_period_results), "\n")

# Create plots

# Plot total sleep time by time period
plot_total_sleep_by_period <- function(sleep_time_periods) {
  # Prepare data
  plot_data <- sleep_time_periods %>%
    mutate(time_period_factor = factor(
      time_period,
      levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
                 "post_surgery_7d", "post_surgery_7d_to_30d", "post_surgery_over_30d"),
      labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
                 "Post 7d", "Post 7d to 30d", "Post >30d")
    )) %>%
    group_by(time_period_factor) %>%
    summarise(
      mean_total_sleep = mean(total_sleep_total, na.rm = TRUE),
      sd_total_sleep = sd(total_sleep_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Convert to hours for better visualization
  plot_data <- plot_data %>%
    mutate(
      mean_total_sleep_hours = mean_total_sleep / 60,
      sd_total_sleep_hours = sd_total_sleep / 60
    )
  
  # Create plot
  ggplot(plot_data, aes(x = time_period_factor, y = mean_total_sleep_hours)) +
    geom_bar(stat = "identity", fill = "#7CB9E8", alpha = 0.7, width = 0.6) +
    geom_errorbar(aes(ymin = mean_total_sleep_hours - sd_total_sleep_hours/60, 
                      ymax = mean_total_sleep_hours + sd_total_sleep_hours/60),
                  width = 0.2, alpha = 0.7) +
    labs(title = "Total Sleep Time by Time Period",
         x = "Time Period",
         y = "Total Sleep Time (hours)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# Plot sleep composition by time period
plot_sleep_composition <- function(sleep_time_periods) {
  # Prepare data
  plot_data <- sleep_time_periods %>%
    mutate(time_period_factor = factor(
      time_period,
      levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
                 "post_surgery_7d", "post_surgery_7d_to_30d", "post_surgery_over_30d"),
      labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
                 "Post 7d", "Post 7d to 30d", "Post >30d")
    )) %>%
    group_by(time_period_factor) %>%
    summarise(
      light_sleep = mean(light_sleep_total, na.rm = TRUE) / 60,
      deep_sleep = mean(deep_sleep_total, na.rm = TRUE) / 60,
      dream_sleep = mean(dream_sleep_total, na.rm = TRUE) / 60,
      awake = mean(awake_total, na.rm = TRUE) / 60,
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(light_sleep, deep_sleep, dream_sleep, awake),
      names_to = "sleep_type",
      values_to = "hours"
    ) %>%
    mutate(sleep_type = factor(sleep_type, 
                               levels = c("light_sleep", "deep_sleep", "dream_sleep", "awake"),
                               labels = c("Light Sleep", "Deep Sleep", "REM Sleep", "Awake")
    ))
  
  # Create plot
  ggplot(plot_data, aes(x = time_period_factor, y = hours, fill = sleep_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Blues", direction = -1) +
    labs(title = "Sleep Composition by Time Period",
         x = "Time Period",
         y = "Duration (hours)",
         fill = "Sleep Phase") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
}

# Plot sleep score by time period
plot_sleep_score <- function(sleep_time_periods) {
  # Prepare data
  plot_data <- sleep_time_periods %>%
    mutate(time_period_factor = factor(
      time_period,
      levels = c("pre_surgery_all", "pre_surgery_7d_all", "pre_surgery_3d_to_7d", "pre_surgery_3d",
                 "post_surgery_7d", "post_surgery_7d_to_30d", "post_surgery_over_30d"),
      labels = c("Pre all", "Pre 7d all", "Pre 3d to 7d", "Pre 3d", 
                 "Post 7d", "Post 7d to 30d", "Post >30d")
    )) %>%
    group_by(time_period_factor) %>%
    summarise(
      mean_score = mean(sleep_score_mean, na.rm = TRUE),
      sd_score = sd(sleep_score_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create plot
  ggplot(plot_data, aes(x = time_period_factor, y = mean_score)) +
    geom_bar(stat = "identity", fill = "#4682B4", alpha = 0.7, width = 0.6) +
    geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score),
                  width = 0.2, alpha = 0.7) +
    labs(title = "Average Sleep Score by Time Period",
         x = "Time Period",
         y = "Sleep Score") +
    scale_y_continuous(limits = c(0, 100)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# Create and save plots
p_total_sleep <- plot_total_sleep_by_period(sleep_time_periods)
p_composition <- plot_sleep_composition(sleep_time_periods)
p_score <- plot_sleep_score(sleep_time_periods)
p_total_sleep 
p_composition
p_score

# Save plots
ggsave("total_sleep_by_time_period.pdf", p_total_sleep, width = 10, height = 6, dpi = 300)
ggsave("sleep_composition_by_time_period.pdf", p_composition, width = 10, height = 6, dpi = 300)
ggsave("sleep_score_by_time_period.pdf", p_score, width = 10, height = 6, dpi = 300)

# Print summary stats
print("Sleep Summary Statistics by Time Period:")
print(sleep_time_period_summary_stats)
