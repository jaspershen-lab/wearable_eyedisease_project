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

# Enhanced function with time-based thresholding (similar to RHR analysis)
calculate_daily_bo_with_thresholds <- function(data, baseline_info) {
  # Set fixed time ranges
  PRE_SURGERY_DAYS <- -10
  POST_SURGERY_DAYS <- 90
  
  # Get subject IDs from blood oxygen data
  bo_subjects <- data@sample_info %>%
    dplyr::select(subject_id) %>%
    distinct() %>%
    pull(subject_id)
  
  # Process baseline information
  baseline_info_processed <- baseline_info %>%
    dplyr::filter(ID %in% bo_subjects) %>%
    mutate(
      surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
    ) %>%
    dplyr::select(ID, surgery_time_1)
  
  # Extract and process blood oxygen data with additional thresholding
  process_blood_oxygen_data <- function() {
    blood_oxygen_values <- data@expression_data[1, ] %>%
      as.numeric()
    
    data.frame(
      sample_id = colnames(data@expression_data),
      blood_oxygen = blood_oxygen_values,
      timestamp = data@sample_info$measure_time,
      subject_id = data@sample_info$subject_id
    ) %>%
      dplyr::filter(
        blood_oxygen >= 70 & blood_oxygen <= 100,  # Valid blood oxygen range
        subject_id %in% baseline_info_processed$ID
      ) %>%
      # Add time-based information for thresholding (same as RHR)
      mutate(
        hour_of_day = hour(timestamp),
        day_of_week = wday(timestamp),
        is_weekend = wday(timestamp, week_start = 1) %in% c(6, 7), # Saturday and Sunday
        month = month(timestamp),
        # Create time periods based on hour of day
        time_period = case_when(
          hour_of_day >= 6 & hour_of_day < 18 ~ "daytime",
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
  calculate_stats_by_threshold <- function(bo_data) {
    # Join with baseline info
    bo_data <- bo_data %>%
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
    stats_list$daily_all <- bo_data %>%
      group_by(subject_id, day_point) %>%
      summarise(
        daily_bo_mean = mean(blood_oxygen),
        daily_bo_min = min(blood_oxygen),
        daily_bo_max = max(blood_oxygen),
        daily_bo_median = median(blood_oxygen),
        daily_bo_sd = sd(blood_oxygen),
        daily_bo_cv = (sd(blood_oxygen) / mean(blood_oxygen)) * 100,
        daily_bo_iqr = IQR(blood_oxygen),
        daily_bo_skew = skewness(blood_oxygen),
        daily_bo_kurt = kurtosis(blood_oxygen),
        n_measurements = n(),
        threshold = "all",
        .groups = "drop"
      )
    
    # By time period
    stats_list$time_period <- bo_data %>%
      group_by(subject_id, day_point, time_period) %>%
      summarise(
        daily_bo_mean = mean(blood_oxygen),
        daily_bo_min = min(blood_oxygen),
        daily_bo_max = max(blood_oxygen),
        daily_bo_median = median(blood_oxygen),
        daily_bo_sd = sd(blood_oxygen),
        daily_bo_cv = (sd(blood_oxygen) / mean(blood_oxygen)) * 100,
        daily_bo_iqr = IQR(blood_oxygen),
        daily_bo_skew = skewness(blood_oxygen),
        daily_bo_kurt = kurtosis(blood_oxygen),
        n_measurements = n(),
        threshold = paste0("time_", time_period),
        .groups = "drop"
      )
    
    # By day type (weekend/weekday)
    stats_list$day_type <- bo_data %>%
      group_by(subject_id, day_point, is_weekend) %>%
      summarise(
        daily_bo_mean = mean(blood_oxygen),
        daily_bo_min = min(blood_oxygen),
        daily_bo_max = max(blood_oxygen),
        daily_bo_median = median(blood_oxygen),
        daily_bo_sd = sd(blood_oxygen),
        daily_bo_cv = (sd(blood_oxygen) / mean(blood_oxygen)) * 100,
        daily_bo_iqr = IQR(blood_oxygen),
        daily_bo_skew = skewness(blood_oxygen),
        daily_bo_kurt = kurtosis(blood_oxygen),
        n_measurements = n(),
        threshold = ifelse(is_weekend, "day_weekend", "day_weekday"),
        .groups = "drop"
      ) %>%
      dplyr::select(-is_weekend)
    
    # By season
    stats_list$season <- bo_data %>%
      group_by(subject_id, day_point, season) %>%
      summarise(
        daily_bo_mean = mean(blood_oxygen),
        daily_bo_min = min(blood_oxygen),
        daily_bo_max = max(blood_oxygen),
        daily_bo_median = median(blood_oxygen),
        daily_bo_sd = sd(blood_oxygen),
        daily_bo_cv = (sd(blood_oxygen) / mean(blood_oxygen)) * 100,
        daily_bo_iqr = IQR(blood_oxygen),
        daily_bo_skew = skewness(blood_oxygen),
        daily_bo_kurt = kurtosis(blood_oxygen),
        n_measurements = n(),
        threshold = paste0("season_", season),
        .groups = "drop"
      )
    
    # Combined time period and day type
    stats_list$time_day <- bo_data %>%
      group_by(subject_id, day_point, time_period, is_weekend) %>%
      summarise(
        daily_bo_mean = mean(blood_oxygen),
        daily_bo_min = min(blood_oxygen),
        daily_bo_max = max(blood_oxygen),
        daily_bo_median = median(blood_oxygen),
        daily_bo_sd = sd(blood_oxygen),
        daily_bo_cv = (sd(blood_oxygen) / mean(blood_oxygen)) * 100,
        daily_bo_iqr = IQR(blood_oxygen),
        daily_bo_skew = skewness(blood_oxygen),
        daily_bo_kurt = kurtosis(blood_oxygen),
        n_measurements = n(),
        threshold = paste0("time_", time_period, "_", ifelse(is_weekend, "weekend", "weekday")),
        .groups = "drop"
      ) %>%
      dplyr::select(-is_weekend)
    
    # Combine all results
    result <- bind_rows(stats_list)
    
    return(result)
  }
  
  # Calculate statistics
  thresholded_stats <- calculate_stats_by_threshold(process_blood_oxygen_data())
  
  # Create wide format
  wide_result <- thresholded_stats %>%
    mutate(
      column_name = paste0("day_", day_point, "_", threshold, "_bo")
    ) %>%
    pivot_wider(
      id_cols = subject_id,
      names_from = c(column_name),
      values_from = c(daily_bo_mean, daily_bo_min, daily_bo_max, daily_bo_median, 
                      daily_bo_sd, daily_bo_cv, daily_bo_iqr, daily_bo_skew, daily_bo_kurt)
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

# Count valid entries by threshold type
create_threshold_summary <- function(result_data) {
  threshold_count_summary <- names(result_data) %>%
    tibble(column_name = .) %>%
    filter(grepl("daily_bo", column_name)) %>%
    mutate(
      day = as.numeric(str_extract(column_name, "-?\\d+")),
      threshold = str_extract(column_name, "time_[^_]+|day_[^_]+|season_[^_]+|all|time_[^_]+_[^_]+"),
      stat_type = str_extract(column_name, "mean|min|max|median|sd|cv|iqr|skew|kurt"),
      valid_count = sapply(result_data[column_name], function(x) sum(!is.na(x)))
    )
  
  # Summarize by threshold type
  threshold_summary <- threshold_count_summary %>%
    group_by(threshold, stat_type) %>%
    summarise(
      avg_valid_count = mean(valid_count, na.rm = TRUE),
      min_valid_count = min(valid_count, na.rm = TRUE),
      max_valid_count = max(valid_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(threshold_summary)
}

# Run the enhanced function
thresholded_bo_result <- calculate_daily_bo_with_thresholds(blood_oxygen_data, baseline_info)

# Create threshold summary
threshold_summary <- create_threshold_summary(thresholded_bo_result)

# Save the results
save(thresholded_bo_result, file = "thresholded_bo_result.rda", compress = "xz")
save(threshold_summary, file = "bo_threshold_summary.rda", compress = "xz")

# Print summary
print(threshold_summary)
