library(tidyverse)
library(tidymass)
library(r4projects)
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


calculate_daily_sleep <- function(data, baseline_info) {
  # Set fixed time ranges
  PRE_SURGERY_DAYS <- -10
  POST_SURGERY_DAYS <- 90
  
  # Get subject IDs from sleep data
  sleep_subjects <- data@sample_info %>%
    dplyr::select(subject_id) %>%
    distinct() %>%
    pull(subject_id)
  
  # Process baseline information
  baseline_info_processed <- baseline_info %>%
    dplyr::filter(ID %in% sleep_subjects) %>%
    mutate(
      surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
    ) %>%
    dplyr::select(ID, surgery_time_1)
  
  # Extract and process sleep data
  process_sleep_data <- function() {
    # Get all sleep metrics
    sleep_metrics <- data@expression_data %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()
    
    # Add sample information
    sleep_metrics <- sleep_metrics %>%
      mutate(
        sample_id = rownames(.),
        timestamp = data@sample_info$measure_time,
        subject_id = data@sample_info$subject_id,
        sleep_start = data@sample_info$sleep_start_time,
        sleep_end = data@sample_info$sleep_end_time
      )
    
    colnames(sleep_metrics)[1:7] <- c(
      "light_sleep",
      "deep_sleep",
      "dream_sleep",
      "awake",
      "total_sleep",
      "daytime_sleep",
      "sleep_score"
    )
    
    sleep_metrics %>%
      dplyr::filter(subject_id %in% baseline_info_processed$ID)
  }
  
  # Calculate daily values
  calculate_daily_stats <- function(sleep_data) {
    sleep_data %>%
      left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
      mutate(
        # sleep_start times计算天数
        days_to_surgery = as.numeric(difftime(sleep_start, surgery_time_1, units = "days")),
        day_point = floor(days_to_surgery),
        # 计算睡眠持续时间（min）
        sleep_duration = as.numeric(difftime(sleep_end, sleep_start, units = "mins"))
      ) %>%
      dplyr::filter(
        day_point >= PRE_SURGERY_DAYS,
        day_point <= POST_SURGERY_DAYS
      ) %>%
      group_by(subject_id, day_point) %>%
      summarise(
        light_sleep = sum(light_sleep, na.rm = TRUE),
        deep_sleep = sum(deep_sleep, na.rm = TRUE),
        dream_sleep = sum(dream_sleep, na.rm = TRUE),
        awake = sum(awake, na.rm = TRUE),
        total_sleep = sum(total_sleep, na.rm = TRUE),
        daytime_sleep = sum(daytime_sleep, na.rm = TRUE),
        sleep_score = mean(sleep_score, na.rm = TRUE),
        # 添加睡眠时长的统计量
        median_sleep_duration = median(sleep_duration, na.rm = TRUE),
        mean_sleep_duration = mean(sleep_duration, na.rm = TRUE),
        min_sleep_duration = min(sleep_duration, na.rm = TRUE),
        max_sleep_duration = max(sleep_duration, na.rm = TRUE),
        sd_sleep_duration = sd(sleep_duration, na.rm = TRUE),
        n_measudreaments = n(),
        .groups = "drop"
      )
  }
  
  # Calculate statistics
  daily_stats <- calculate_daily_stats(process_sleep_data())
  
  # Create wide format
  create_wide_format <- function(daily_stats) {
    metrics <- c(
      "light_sleep",
      "deep_sleep",
      "dream_sleep",
      "awake",
      "total_sleep",
      "daytime_sleep",
      "sleep_score",
      "median_sleep_duration",
      "mean_sleep_duration",
      "min_sleep_duration",
      "max_sleep_duration",
      "sd_sleep_duration"
    )
    
    result <- daily_stats %>%
      pivot_longer(
        cols = all_of(metrics),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(col_name = paste0("day_", day_point, "_", metric)) %>%
      dplyr::select(subject_id, col_name, value) %>%
      pivot_wider(
        names_from = col_name,
        values_from = value
      )
    
    return(result)
  }
  
  # Process data
  sleep_wide <- create_wide_format(daily_stats)
  
  # Create final result
  result <- tibble(subject_id = unique(baseline_info_processed$ID)) %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = as.Date(surgery_time_1)) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    left_join(sleep_wide, by = "subject_id") %>%
    distinct() %>%
    arrange(subject_id)
  
  # Sort columns
  col_order <- c("subject_id", "surgery_date")
  time_cols <- setdiff(names(result), col_order)
  sorted_time_cols <- time_cols[order(as.numeric(gsub(".*day_(.+)_.*", "\\1", time_cols)))]
  result <- result[, c(col_order, sorted_time_cols)]
  
  # Validation
  cat("\nFinal data validation:\n")
  cat("Number of unique subjects: ", nrow(result), "\n")
  cat("Number of columns: ", ncol(result), "\n")
  
  return(result)
}

# Function to calculate summary statistics
calculate_daily_sleep_summary <- function(data) {
  # Get all time columns (excluding subject_id and surgery_date)
  time_cols <- names(data)[-(1:2)]
  
  # Extract day and metric type
  summary_stats <- data.frame(
    column = time_cols,
    day = as.numeric(str_extract(time_cols, "-?\\d+")),
    metric = str_extract(time_cols, "(light_sleep|deep_sleep|dream_sleep|awake|total_sleep|daytime_sleep|sleep_score|median_sleep_duration|mean_sleep_duration|min_sleep_duration|max_sleep_duration|sd_sleep_duration)")
  ) %>%
    mutate(valid_count = sapply(data[time_cols], function(x) sum(!is.na(x))))
  
  return(summary_stats)
}

# Function to calculate daily non-NA counts
calculate_daily_sleep_counts <- function(data) {
  # Get all time columns (excluding subject_id and surgery_date)
  time_cols <- names(data)[-(1:2)]
  
  # Calculate count of non-NA values for each column
  counts <- sapply(data[time_cols], function(x) sum(!is.na(x)))
  
  # Create result dataframe
  result <- data.frame(
    day = as.numeric(str_extract(names(counts), "-?\\d+")),
    metric = str_extract(names(counts), "(light_sleep|deep_sleep|dream_sleep|awake|total_sleep|daytime_sleep|sleep_score|median_sleep_duration|mean_sleep_duration|min_sleep_duration|max_sleep_duration|sd_sleep_duration)"),
    valid_count = counts
  ) %>%
    arrange(day, metric)
  
  return(result)
}


daily_sleep_result <- calculate_daily_sleep(sleep_data, baseline_info)
summary_stats <- calculate_daily_sleep_summary(daily_sleep_result)
daily_counts <- calculate_daily_sleep_counts(daily_sleep_result)

save(daily_sleep_result, file = "daily_sleep_result.rda", compress = "xz")

