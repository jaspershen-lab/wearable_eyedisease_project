library(tidyverse)
library(tidymass)
library(r4projects)
library(moments)  # 添加moments包用于计算偏度和峰度
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


###########
calculate_daily_blood_oxygen <- function(data, baseline_info) {
  # Set fixed time ranges (same as RHR analysis)
  PRE_SURGERY_DAYS <- -10
  POST_SURGERY_DAYS <- 90
  
  # Get subject IDs from blood oxygen data
  blood_oxygen_subjects <- data@sample_info %>%
    dplyr::select(subject_id) %>%
    distinct() %>%
    pull(subject_id)
  
  # Process baseline information
  baseline_info_processed <- baseline_info %>%
    dplyr::filter(ID %in% blood_oxygen_subjects) %>%
    mutate(
      surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
    ) %>%
    dplyr::select(ID, surgery_time_1)
  
  # Extract and process blood oxygen data
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
      )
  }
  
  # Calculate daily statistics
  calculate_daily_stats <- function(blood_oxygen_data) {
    blood_oxygen_data %>%
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
        daily_bo_mean = mean(blood_oxygen),
        daily_bo_min = min(blood_oxygen),
        daily_bo_max = max(blood_oxygen),
        daily_bo_median = median(blood_oxygen),
        daily_bo_sd = sd(blood_oxygen),
        daily_bo_cv = (sd(blood_oxygen) / mean(blood_oxygen)) * 100, # CV as percentage
        daily_bo_iqr = IQR(blood_oxygen),
        daily_bo_skew = skewness(blood_oxygen),  # 添加偏度计算
        daily_bo_kurt = kurtosis(blood_oxygen),  # 添加峰度计算
        n_measurements = n(),
        .groups = "drop"
      )
  }
  
  # Calculate statistics
  daily_stats <- calculate_daily_stats(process_blood_oxygen_data())
  
  # Create wide format for each statistic
  create_wide_format <- function(daily_stats) {
    stats_list <- list()
    
    # List of statistics to process
    stat_cols <- c("mean", "min", "max","median", "sd", "cv", "iqr", "skew", "kurt")
    
    for(stat in stat_cols) {
      col_name <- paste0("daily_bo_", stat)
      stats_list[[stat]] <- daily_stats %>%
        mutate(col_name = paste0("day_", day_point, "_", stat, "_bo")) %>%
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
  
  # Process data
  bo_wide <- create_wide_format(daily_stats)
  
  # Create final result
  result <- tibble(subject_id = unique(baseline_info_processed$ID)) %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = as.Date(surgery_time_1)) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    left_join(bo_wide, by = "subject_id") %>%
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
calculate_daily_bo_summary <- function(data) {
  # Get all time columns (excluding subject_id and surgery_date)
  time_cols <- names(data)[-(1:2)]
  
  # Extract statistics type and day
  summary_stats <- data.frame(
    column = time_cols,
    day = as.numeric(str_extract(time_cols, "-?\\d+")),
    stat_type = str_extract(time_cols, "(mean|min|max|sd|cv|iqr|skew|kurt)")
  ) %>%
    mutate(valid_count = sapply(data[time_cols], function(x) sum(!is.na(x))))
  
  return(summary_stats)
}

# Function to calculate daily non-NA counts
calculate_daily_bo_counts <- function(data) {
  # Get all time columns (excluding subject_id and surgery_date)
  time_cols <- names(data)[-(1:2)]
  
  # Calculate count of non-NA values for each column
  counts <- sapply(data[time_cols], function(x) sum(!is.na(x)))
  
  # Create result dataframe
  result <- data.frame(
    day = as.numeric(str_extract(names(counts), "-?\\d+")),
    stat_type = str_extract(names(counts), "(mean|min|max|median|sd|cv|iqr|skew|kurt)"),  # 添加新统计量
    valid_count = counts
  ) %>%
    arrange(day, stat_type)
  
  return(result)
}


daily_bo_result <- calculate_daily_blood_oxygen(blood_oxygen_data, baseline_info)
summary_stats <- calculate_daily_bo_summary(daily_bo_result)
daily_counts <- calculate_daily_bo_counts(daily_bo_result)

save(daily_bo_result, file = "daily_bo_result.rda", compress = "xz")
