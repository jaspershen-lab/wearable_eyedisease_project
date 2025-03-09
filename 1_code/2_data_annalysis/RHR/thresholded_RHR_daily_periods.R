library(tidyverse)
library(tidymass)
library(r4projects)
library(moments)  # 添加moments包用于计算偏度和峰度
setwd(get_project_wd())
rm(list = ls())
library(lubridate)


######## Read data
load(
  "3_data_analysis/2_data_analysis/RHR/heart_rate_data_RHR.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/RHR/RHR_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/RHR/RHR_time_periods")

# Enhanced function with time-based thresholding
calculate_daily_rhr_with_thresholds <- function(data, baseline_info) {
  # Set fixed time ranges
  PRE_SURGERY_DAYS <- -10
  POST_SURGERY_DAYS <- 90
  
  # Get subject IDs from RHR data
  rhr_subjects <- data@sample_info %>%
    dplyr::select(subject_id) %>%
    distinct() %>%
    pull(subject_id)
  
  # Process baseline information
  baseline_info_processed <- baseline_info %>%
    dplyr::filter(ID %in% rhr_subjects) %>%
    mutate(
      surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
    ) %>%
    dplyr::select(ID, surgery_time_1)
  
  # Extract and process heart rate data with additional thresholding
  # More reliable approach to determine weekends
  process_heart_rate_data <- function(label_value) {
    filtered_sample_info <- data@sample_info %>%
      dplyr::filter(label == label_value)
    
    heart_rate_values <- data@expression_data[1, filtered_sample_info$sample_id] %>%
      as.numeric()
    
    data.frame(
      sample_id = filtered_sample_info$sample_id,
      heart_rate = heart_rate_values,
      timestamp = filtered_sample_info$measure_time,
      subject_id = filtered_sample_info$subject_id
    ) %>%
      dplyr::filter(
        heart_rate >= 30 & heart_rate <= 200,
        subject_id %in% baseline_info_processed$ID
      ) %>%
      # Add time-based information for thresholding
      mutate(
        hour_of_day = hour(timestamp),
        # Ensure weekend determination is correct
        # wday() returns 1-7 with 1 being Sunday by default
        day_of_week = wday(timestamp),
        is_weekend = wday(timestamp, week_start = 1) %in% c(6, 7), # Saturday and Sunday
        month = month(timestamp),
        # Create time periods based on hour of day
        time_period = case_when(
          # hour_of_day >= 7 & hour_of_day <= 9 ~ "morning",
          hour_of_day >= 6 & hour_of_day < 18 ~ "daytime",
          # hour_of_day >= 17 & hour_of_day <= 22 ~ "evening",
          TRUE ~ "nighttime"  # This covers hours 18-23 and 0-5
        ),
        # Create season based on month (Northern Hemisphere)
        season = case_when(
          month %in% c(12, 1, 2) ~ "winter",
          month %in% c(3, 4, 5) ~ "spring",
          month %in% c(6, 7, 8, 9) ~ "summer",
          month %in% c(10, 11) ~ "fall"
        )
      )
  }
  
  # Calculate statistics by various thresholds
  calculate_stats_by_threshold <- function(rhr_data) {
    # Join with baseline info
    rhr_data <- rhr_data %>%
      left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
      mutate(
        days_to_surgery = as.numeric(difftime(timestamp, surgery_time_1, units = "days")),
        day_point = floor(days_to_surgery)
      ) %>%
      dplyr::filter(
        day_point >= PRE_SURGERY_DAYS,
        day_point <= POST_SURGERY_DAYS
      )
    
    # Calculate statistics for each threshold combination
    stats_list <- list()
    
    # Base daily stats (all data combined)
    stats_list$daily_all <- rhr_data %>%
      group_by(subject_id, day_point) %>%
      summarise(
        daily_rhr_mean = mean(heart_rate),
        daily_rhr_min = min(heart_rate),
        daily_rhr_max = max(heart_rate),
        daily_rhr_median = median(heart_rate),
        daily_rhr_sd = sd(heart_rate),
        daily_rhr_cv = (sd(heart_rate) / mean(heart_rate)) * 100,
        daily_rhr_iqr = IQR(heart_rate),
        daily_rhr_skew = skewness(heart_rate),
        daily_rhr_kurt = kurtosis(heart_rate),
        n_measurements = n(),
        threshold = "all",
        .groups = "drop"
      )
    
    # By time period
    stats_list$time_period <- rhr_data %>%
      group_by(subject_id, day_point, time_period) %>%
      summarise(
        daily_rhr_mean = mean(heart_rate),
        daily_rhr_min = min(heart_rate),
        daily_rhr_max = max(heart_rate),
        daily_rhr_median = median(heart_rate),
        daily_rhr_sd = sd(heart_rate),
        daily_rhr_cv = (sd(heart_rate) / mean(heart_rate)) * 100,
        daily_rhr_iqr = IQR(heart_rate),
        daily_rhr_skew = skewness(heart_rate),
        daily_rhr_kurt = kurtosis(heart_rate),
        n_measurements = n(),
        threshold = paste0("time_", time_period),
        .groups = "drop"
      )
    
    # By day type (weekend/weekday)
    stats_list$day_type <- rhr_data %>%
      group_by(subject_id, day_point, is_weekend) %>%
      summarise(
        daily_rhr_mean = mean(heart_rate),
        daily_rhr_min = min(heart_rate),
        daily_rhr_max = max(heart_rate),
        daily_rhr_median = median(heart_rate),
        daily_rhr_sd = sd(heart_rate),
        daily_rhr_cv = (sd(heart_rate) / mean(heart_rate)) * 100,
        daily_rhr_iqr = IQR(heart_rate),
        daily_rhr_skew = skewness(heart_rate),
        daily_rhr_kurt = kurtosis(heart_rate),
        n_measurements = n(),
        threshold = ifelse(is_weekend, "day_weekend", "day_weekday"),
        .groups = "drop"
      ) %>%
      dplyr::select(-is_weekend)
    
    # By season
    stats_list$season <- rhr_data %>%
      group_by(subject_id, day_point, season) %>%
      summarise(
        daily_rhr_mean = mean(heart_rate),
        daily_rhr_min = min(heart_rate),
        daily_rhr_max = max(heart_rate),
        daily_rhr_median = median(heart_rate),
        daily_rhr_sd = sd(heart_rate),
        daily_rhr_cv = (sd(heart_rate) / mean(heart_rate)) * 100,
        daily_rhr_iqr = IQR(heart_rate),
        daily_rhr_skew = skewness(heart_rate),
        daily_rhr_kurt = kurtosis(heart_rate),
        n_measurements = n(),
        threshold = paste0("season_", season),
        .groups = "drop"
      )
    
    # Combined time period and day type
    stats_list$time_day <- rhr_data %>%
      group_by(subject_id, day_point, time_period, is_weekend) %>%
      summarise(
        daily_rhr_mean = mean(heart_rate),
        daily_rhr_min = min(heart_rate),
        daily_rhr_max = max(heart_rate),
        daily_rhr_median = median(heart_rate),
        daily_rhr_sd = sd(heart_rate),
        daily_rhr_cv = (sd(heart_rate) / mean(heart_rate)) * 100,
        daily_rhr_iqr = IQR(heart_rate),
        daily_rhr_skew = skewness(heart_rate),
        daily_rhr_kurt = kurtosis(heart_rate),
        n_measurements = n(),
        threshold = paste0("time_", time_period, "_", ifelse(is_weekend, "weekend", "weekday")),
        .groups = "drop"
      ) %>%
      dplyr::select(-is_weekend)
    
    # Combine all results
    result <- bind_rows(stats_list)
    
    return(result)
  }
  
  # Calculate statistics for both labels
  thresholded_stats_1 <- calculate_stats_by_threshold(process_heart_rate_data("<1"))
  thresholded_stats_50 <- calculate_stats_by_threshold(process_heart_rate_data("<50"))
  
  # Add label identifier
  thresholded_stats_1 <- thresholded_stats_1 %>% mutate(rhr_type = "rhr_1")
  thresholded_stats_50 <- thresholded_stats_50 %>% mutate(rhr_type = "rhr_50")
  
  # Combine results
  all_stats <- bind_rows(thresholded_stats_1, thresholded_stats_50)
  
  # Create wide format
  wide_result <- all_stats %>%
    mutate(
      column_name = paste0("day_", day_point, "_", threshold, "_", rhr_type)
    ) %>%
    pivot_wider(
      id_cols = subject_id,
      names_from = c(column_name),
      values_from = c(daily_rhr_mean, daily_rhr_min, daily_rhr_max, daily_rhr_median, 
                      daily_rhr_sd, daily_rhr_cv, daily_rhr_iqr, daily_rhr_skew, daily_rhr_kurt)
    )
  
  # Add surgery date
  final_result <- wide_result %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = as.Date(surgery_time_1)) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    relocate(subject_id, surgery_date)
  
  # Validation
  cat("\nFinal data validation:\n")
  cat("Number of unique subjects: ", nrow(final_result), "\n")
  cat("Number of columns: ", ncol(final_result), "\n")
  
  return(final_result)
}

# Run the enhanced function
thresholded_rhr_result <- calculate_daily_rhr_with_thresholds(heart_rate_data, baseline_info)

# Save the result
save(thresholded_rhr_result, file = "thresholded_rhr_result.rda", compress = "xz")

# Count valid entries by threshold type
threshold_count_summary <- names(thresholded_rhr_result) %>%
  tibble(column_name = .) %>%
  filter(grepl("daily_rhr", column_name)) %>%
  mutate(
    day = as.numeric(str_extract(column_name, "-?\\d+")),
    threshold = str_extract(column_name, "time_[^_]+|day_[^_]+|season_[^_]+|all|time_[^_]+_[^_]+"),
    stat_type = str_extract(column_name, "mean|min|max|median|sd|cv|iqr|skew|kurt"),
    rhr_type = str_extract(column_name, "rhr_\\d+"),
    valid_count = sapply(thresholded_rhr_result[column_name], function(x) sum(!is.na(x)))
  )

# Summarize by threshold type
threshold_summary <- threshold_count_summary %>%
  group_by(threshold, stat_type, rhr_type) %>%
  summarise(
    avg_valid_count = mean(valid_count, na.rm = TRUE),
    min_valid_count = min(valid_count, na.rm = TRUE),
    max_valid_count = max(valid_count, na.rm = TRUE),
    .groups = "drop"
  )

print(threshold_summary)