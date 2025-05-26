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

# Load heart rate data (added from second code)
load(
  "3_data_analysis/2_data_analysis/RHR/heart_rate_data_RHR.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/steps/steps_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/steps/steps_time_periods")


# daily_workout_details_data@sample_info <- daily_workout_details_data@sample_info %>%
#   mutate(
#     measure_time = as.POSIXct(
#       str_extract(sample_id, "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"),
#       format = "%Y-%m-%d %H:%M:%S"
#     )
#   )

# Function to preprocess steps data with heart rate data (added from second code)
preprocess_steps_data <- function(steps_data, heart_rate_data) {
  # Extract heart rate data
  hr_sample_info <- heart_rate_data@sample_info
  hr_values <- as.numeric(heart_rate_data@expression_data[1, ])
  
  # Create dataframe with heart rate values
  hr_df <- data.frame(
    sample_id = colnames(heart_rate_data@expression_data),
    heart_rate = hr_values
  ) %>%
    left_join(hr_sample_info, by = "sample_id") %>%
    filter(!is.na(heart_rate)) %>%
    dplyr::select(subject_id, measure_time) %>%
    distinct()
  
  # Extract steps data
  steps_sample_info <- steps_data@sample_info
  
  # Ensure expression_data is a data.frame
  if(is.matrix(steps_data@expression_data)) {
    steps_expr <- as.data.frame(steps_data@expression_data)
  } else {
    steps_expr <- steps_data@expression_data
  }
  
  steps_values <- as.numeric(steps_expr["steps", ])
  
  # Create dataframe with steps values
  steps_df <- data.frame(
    sample_id = colnames(steps_expr),
    steps = steps_values
  ) %>%
    left_join(steps_sample_info, by = "sample_id")
  
  # Join steps data with heart rate timestamps
  merged_data <- hr_df %>%
    left_join(
      steps_df,
      by = c("subject_id", "measure_time")
    ) %>%
    mutate(
      # When heart rate exists, fill NA steps with 0
      steps = ifelse(is.na(steps), 0, steps),
      # Create or update sample_id
      sample_id = ifelse(is.na(sample_id), 
                         paste0(subject_id, "_", format(measure_time, "%Y%m%d%H%M%S")),
                         sample_id)
    )
  
  # Ensure activity and class columns exist
  if("activity" %in% names(steps_sample_info) && !("activity" %in% names(merged_data))) {
    merged_data$activity <- NA
  }
  if("activity" %in% names(merged_data)) {
    merged_data$activity[is.na(merged_data$activity)] <- "unknown"
  }
  
  if("class" %in% names(steps_sample_info) && !("class" %in% names(merged_data))) {
    merged_data$class <- NA
  }
  if("class" %in% names(merged_data)) {
    merged_data$class[is.na(merged_data$class)] <- "unknown"
  }
  
  # Create updated sample_info
  required_columns <- names(steps_sample_info)
  missing_columns <- setdiff(required_columns, names(merged_data))
  
  if(length(missing_columns) > 0) {
    merged_data[, missing_columns] <- NA
  }
  
  # Ensure steps column exists in merged_data
  if(!"steps" %in% names(merged_data)) {
    merged_data$steps <- 0
  }
  
  # Select all necessary columns and ensure uniqueness
  updated_sample_info <- merged_data %>%
    dplyr::select(all_of(required_columns), steps) %>%  # Temporarily keep steps column
    distinct(sample_id, .keep_all = TRUE)
  
  # Create updated expression_data (as data.frame)
  row_names <- rownames(steps_expr)
  
  # Initialize an empty dataframe
  new_expression_data <- as.data.frame(matrix(NA, 
                                              nrow = length(row_names), 
                                              ncol = nrow(updated_sample_info)))
  
  # Set row and column names
  rownames(new_expression_data) <- row_names
  colnames(new_expression_data) <- updated_sample_info$sample_id
  
  # Find the index of the steps row
  steps_row_index <- which(rownames(new_expression_data) == "steps")
  
  # Vectorized assignment: directly use steps column from updated_sample_info
  if(length(steps_row_index) > 0) {
    new_expression_data[steps_row_index, ] <- updated_sample_info$steps
  }
  
  # Remove temporarily kept steps column
  updated_sample_info <- updated_sample_info %>% 
    dplyr::select(all_of(required_columns))
  
  # Create updated steps data object
  updated_steps_data <- steps_data
  updated_steps_data@expression_data <- new_expression_data
  updated_steps_data@sample_info <- updated_sample_info
  
  # Print summary and validation info
  cat("Original number of steps data points:", nrow(steps_sample_info), "\n")
  cat("After filling with heart rate data:", nrow(updated_sample_info), "\n")
  cat("Number of added zero steps:", nrow(updated_sample_info) - nrow(steps_sample_info), "\n")
  cat("Number of columns in expression_data:", ncol(new_expression_data), "\n")
  cat("Number of rows in sample_info:", nrow(updated_sample_info), "\n")
  cat("Sample IDs in expression_data match sample_info:", 
      all(colnames(new_expression_data) %in% updated_sample_info$sample_id), "\n")
  
  return(updated_steps_data)
}

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
    metric = str_extract(time_cols, "(total|mean|max|median)")
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
    metric = str_extract(names(counts), "(total|mean|max|median)"),
    valid_count = counts
  ) %>%
    arrange(day, metric)
  
  return(result)
}

# Preprocess steps data with heart rate data
preprocessed_steps_data <- preprocess_steps_data(daily_workout_details_data, heart_rate_data)

# Calculate daily steps using the preprocessed data
daily_steps_result <- calculate_daily_steps(preprocessed_steps_data, baseline_info)
summary_stats <- calculate_daily_steps_summary(daily_steps_result)
daily_counts <- calculate_daily_steps_counts(daily_steps_result)

# Save the preprocessed data and results
save(preprocessed_steps_data, file = "preprocessed_steps_data.rda", compress = "xz")
save(daily_steps_result, file = "daily_steps_result_assigned.rda", compress = "xz")

# Print summary of the preprocessing
cat("Preprocessing summary:\n")
cat("Original steps data points:", nrow(daily_workout_details_data@sample_info), "\n")
cat("After filling with heart rate data:", nrow(preprocessed_steps_data@sample_info), "\n")
cat("Number of added zero steps:", 
    nrow(preprocessed_steps_data@sample_info) - nrow(daily_workout_details_data@sample_info), 
    "\n")
