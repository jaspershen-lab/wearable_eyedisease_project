library(tidyverse)
library(tidymass)
library(r4projects)
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

daily_workout_details_data@sample_info <- daily_workout_details_data@sample_info %>%
  mutate(
    measure_time = as.POSIXct(
      str_extract(sample_id, "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )
##########
calculate_daily_steps <- function(data, baseline_info) {
  # Set fixed time ranges
  PRE_SURGERY_DAYS <- -10
  POST_SURGERY_DAYS <- 90
  
  # Get subject IDs from data
  steps_subjects <- data@sample_info %>%
    dplyr::select(subject_id) %>%
    distinct() %>%
    pull(subject_id)
  
  # Process baseline information
  baseline_info_processed <- baseline_info %>%
    dplyr::filter(ID %in% steps_subjects) %>%
    mutate(
      surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
    ) %>%
    dplyr::select(ID, surgery_time_1)
  
  # Extract and process steps data
  process_steps_data <- function() {
    # Get steps data from expression data
    steps_values <- data@expression_data["steps", ] %>%
      as.numeric()
    
    data.frame(
      sample_id = colnames(data@expression_data),
      steps = steps_values,
      timestamp = data@sample_info$measure_time,
      subject_id = data@sample_info$subject_id
    ) %>%
      dplyr::filter(
        steps >= 0, # Remove negative values if any
        subject_id %in% baseline_info_processed$ID
      )
  }
  
  # Calculate daily statistics
  calculate_daily_stats <- function(steps_data) {
    steps_data %>%
      left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
      mutate(
        days_to_surgery = as.numeric(difftime(timestamp, surgery_time_1, units = "days")),
        day_point = floor(days_to_surgery)
      ) %>%
      dplyr::filter(
        day_point >= PRE_SURGERY_DAYS,
        day_point <= POST_SURGERY_DAYS
      ) %>%
      group_by(subject_id, day_point) %>%
      summarise(
        daily_steps_total = sum(steps, na.rm = TRUE),
        daily_steps_mean = mean(steps, na.rm = TRUE),
        daily_steps_max = max(steps, na.rm = TRUE),
        daily_steps_median = median(steps, na.rm = TRUE),
        n_measurements = n(),
        .groups = "drop"
      )
  }
  
  # Calculate statistics
  daily_stats <- calculate_daily_stats(process_steps_data())
  
  # Create wide format
  create_wide_format <- function(daily_stats) {
    metrics <- c("total", "mean", "max","median")
    
    result <- daily_stats %>%
      pivot_longer(
        cols = starts_with("daily_steps_"),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        metric = str_replace(metric, "daily_steps_", ""),
        col_name = paste0("day_", day_point, "_steps_", metric)
      ) %>%
      dplyr::select(subject_id, col_name, value) %>%
      pivot_wider(
        names_from = col_name,
        values_from = value
      )
    
    return(result)
  }
  
  # Process data
  steps_wide <- create_wide_format(daily_stats)
  
  # Create final result
  result <- tibble(subject_id = unique(baseline_info_processed$ID)) %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = as.Date(surgery_time_1)) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    left_join(steps_wide, by = "subject_id") %>%
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
calculate_daily_steps_summary <- function(data) {
  # Get all time columns (excluding subject_id and surgery_date)
  time_cols <- names(data)[-(1:2)]
  
  # Extract day and metric type
  summary_stats <- data.frame(
    column = time_cols,
    day = as.numeric(str_extract(time_cols, "-?\\d+")),
    metric = str_extract(time_cols, "(total|mean|max)")
  ) %>%
    mutate(valid_count = sapply(data[time_cols], function(x) sum(!is.na(x))))
  
  return(summary_stats)
}

# Function to calculate daily non-NA counts
calculate_daily_steps_counts <- function(data) {
  # Get all time columns (excluding subject_id and surgery_date)
  time_cols <- names(data)[-(1:2)]
  
  # Calculate count of non-NA values for each column
  counts <- sapply(data[time_cols], function(x) sum(!is.na(x)))
  
  # Create result dataframe
  result <- data.frame(
    day = as.numeric(str_extract(names(counts), "-?\\d+")),
    metric = str_extract(names(counts), "(total|mean|max)"),
    valid_count = counts
  ) %>%
    arrange(day, metric)
  
  return(result)
}


daily_steps_result <- calculate_daily_steps(daily_workout_details_data, baseline_info)
summary_stats <- calculate_daily_steps_summary(daily_steps_result)
daily_counts <- calculate_daily_steps_counts(daily_steps_result)

save(daily_steps_result, file = "daily_steps_result.rda", compress = "xz")
