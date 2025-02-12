library(tidyverse)
library(tidymass)
library(r4projects)
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


########

calculate_daily_rhr <- function(data, baseline_info) {
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
  
  # Extract and process heart rate data
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
      )
  }
  
  # Calculate daily statistics
  calculate_daily_stats <- function(rhr_data) {
    rhr_data %>%
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
        daily_rhr_mean = mean(heart_rate),
        daily_rhr_min = min(heart_rate),
        daily_rhr_max = max(heart_rate),
        daily_rhr_median = median(heart_rate),
        daily_rhr_sd = sd(heart_rate),
        daily_rhr_cv = (sd(heart_rate) / mean(heart_rate)) * 100, # CV as percentage
        daily_rhr_iqr = IQR(heart_rate),
        n_measurements = n(),
        .groups = "drop"
      )
  }
  
  # Calculate statistics for both labels
  daily_stats_1 <- calculate_daily_stats(process_heart_rate_data("<1"))
  daily_stats_50 <- calculate_daily_stats(process_heart_rate_data("<50"))
  
  # Create wide format for each statistic
  create_wide_format <- function(daily_stats, suffix) {
    stats_list <- list()
    
    # List of statistics to process
    stat_cols <- c("mean", "min", "max","median", "sd", "cv", "iqr")
    
    for(stat in stat_cols) {
      col_name <- paste0("daily_rhr_", stat)
      stats_list[[stat]] <- daily_stats %>%
        mutate(col_name = paste0("day_", day_point, "_", stat, "_", suffix)) %>%
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
  
  # Process data for both labels
  rhr_wide_1 <- create_wide_format(daily_stats_1, "rhr_1")
  rhr_wide_50 <- create_wide_format(daily_stats_50, "rhr_50")
  
  # Create final result
  result <- tibble(subject_id = unique(baseline_info_processed$ID)) %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = as.Date(surgery_time_1)) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    left_join(rhr_wide_1, by = "subject_id") %>%
    left_join(rhr_wide_50, by = "subject_id") %>%
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

# Calculate summary statistics
calculate_daily_summary <- function(data) {
  # Get all time columns (excluding subject_id and surgery_date)
  time_cols <- names(data)[-(1:2)]
  
  # Extract statistics type and day
  summary_stats <- data.frame(
    column = time_cols,
    day = as.numeric(str_extract(time_cols, "-?\\d+")),
    stat_type = str_extract(time_cols, "(mean|min|max|sd|cv|iqr)"),
    rhr_type = str_extract(time_cols, "rhr_\\d+")
  ) %>%
    mutate(valid_count = sapply(data[time_cols], function(x) sum(!is.na(x))))
  
  return(summary_stats)
}

# Calculate daily RHR
daily_rhr_result <- calculate_daily_rhr(heart_rate_data, baseline_info)
summary_stats <- calculate_daily_summary(daily_rhr_result)

save(daily_rhr_result, file = "daily_rhr_result.rda", compress = "xz")

# Calculate daily non-NA counts
calculate_daily_counts <- function(data) {
  # Get all time columns (excluding subject_id and surgery_date)
  time_cols <- names(data)[-(1:2)]
  
  # Calculate count of non-NA values for each column
  counts <- sapply(data[time_cols], function(x) sum(!is.na(x)))
  
  # Create result dataframe
  result <- data.frame(
    day = as.numeric(str_extract(names(counts), "-?\\d+")),
    rhr_type = str_extract(names(counts), "rhr_\\d+"),
    valid_count = counts
  ) %>%
    arrange(day, rhr_type)
  
  return(result)
}

daily_counts <- calculate_daily_counts(daily_rhr_result)

print(daily_counts)

