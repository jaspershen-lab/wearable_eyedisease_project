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

load(
  "3_data_analysis/2_data_analysis/RHR/heart_rate_data_RHR.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/steps/steps_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/steps/steps_time_periods")


#############
calculate_daily_steps_with_thresholds <- function(data, baseline_info) {
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
  
  # Extract and process steps data with thresholding
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
      ) %>%
      # Add time-based information for thresholding
      mutate(
        hour_of_day = hour(timestamp),
        # Ensure weekend determination is correct
        is_weekend = wday(timestamp, week_start = 1) %in% c(6, 7), # Saturday and Sunday
        month = month(timestamp),
        # Create simple daytime/nighttime periods
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
  
  # Get the processed data
  processed_data <- process_steps_data()
  
  # Create thresholded results
  create_thresholded_results <- function() {
    results_list <- list()
    
    # Join with baseline info to get days relative to surgery
    processed_with_days <- processed_data %>%
      left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
      mutate(
        days_to_surgery = as.numeric(difftime(timestamp, surgery_time_1, units = "days")),
        day_point = floor(days_to_surgery)
      ) %>%
      dplyr::filter(
        day_point >= PRE_SURGERY_DAYS,
        day_point <= POST_SURGERY_DAYS
      )
    
    # 1. All data (no thresholding)
    all_data <- processed_with_days %>%
      group_by(subject_id, day_point) %>%
      summarise(
        steps_total = sum(steps, na.rm = TRUE),
        steps_mean = mean(steps, na.rm = TRUE),
        steps_max = max(steps, na.rm = TRUE),
        steps_median = median(steps, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(threshold = "all")
    
    # 2. Time period threshold (daytime/nighttime)
    time_period_data <- processed_with_days %>%
      group_by(subject_id, day_point, time_period) %>%
      summarise(
        steps_total = sum(steps, na.rm = TRUE),
        steps_mean = mean(steps, na.rm = TRUE),
        steps_max = max(steps, na.rm = TRUE),
        steps_median = median(steps, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(threshold = paste0("time_", time_period)) %>%
      dplyr::select(-time_period)
    
    # 3. Weekend/weekday threshold
    day_type_data <- processed_with_days %>%
      group_by(subject_id, day_point, is_weekend) %>%
      summarise(
        steps_total = sum(steps, na.rm = TRUE),
        steps_mean = mean(steps, na.rm = TRUE),
        steps_max = max(steps, na.rm = TRUE),
        steps_median = median(steps, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(threshold = ifelse(is_weekend, "day_weekend", "day_weekday")) %>%
      dplyr::select(-is_weekend)
    
    # 4. Season threshold
    season_data <- processed_with_days %>%
      group_by(subject_id, day_point, season) %>%
      summarise(
        steps_total = sum(steps, na.rm = TRUE),
        steps_mean = mean(steps, na.rm = TRUE),
        steps_max = max(steps, na.rm = TRUE),
        steps_median = median(steps, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(threshold = paste0("season_", season)) %>%
      dplyr::select(-season)
    
    # 5. Combined time and day type
    combined_data <- processed_with_days %>%
      group_by(subject_id, day_point, time_period, is_weekend) %>%
      summarise(
        steps_total = sum(steps, na.rm = TRUE),
        steps_mean = mean(steps, na.rm = TRUE),
        steps_max = max(steps, na.rm = TRUE),
        steps_median = median(steps, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        threshold = paste0(
          "time_", time_period, "_", 
          ifelse(is_weekend, "weekend", "weekday")
        )
      ) %>%
      dplyr::select(-time_period, -is_weekend)
    
    # Combine all results
    bind_rows(all_data, time_period_data, day_type_data, season_data, combined_data)
  }
  
  # Get all thresholded results
  all_thresholded_results <- create_thresholded_results()
  
  # Reshape to wide format - each metric separate to avoid duplicate issues
  reshape_to_wide <- function(data, metric) {
    filtered_data <- data %>% 
      dplyr::select(subject_id, day_point, threshold, !!sym(paste0("steps_", metric)))
    
    # Create a unique column name for each day-threshold combination
    filtered_data %>%
      mutate(
        column_name = paste0("day_", day_point, "_", threshold)
      ) %>%
      pivot_wider(
        id_cols = subject_id,
        names_from = column_name,
        values_from = !!sym(paste0("steps_", metric)),
        names_prefix = paste0("daily_steps_", metric, "_")
      )
  }
  
  # Create wide format for each metric
  total_wide <- reshape_to_wide(all_thresholded_results, "total")
  mean_wide <- reshape_to_wide(all_thresholded_results, "mean")
  max_wide <- reshape_to_wide(all_thresholded_results, "max")
  median_wide <- reshape_to_wide(all_thresholded_results, "median")
  
  # Join all the metrics
  result <- total_wide %>%
    left_join(mean_wide, by = "subject_id") %>%
    left_join(max_wide, by = "subject_id") %>%
    left_join(median_wide, by = "subject_id")
  
  # Add surgery date
  final_result <- result %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = as.Date(surgery_time_1)) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    relocate(subject_id, surgery_date)
  
  # Sort columns
  col_order <- c("subject_id", "surgery_date")
  time_cols <- setdiff(names(final_result), col_order)
  
  # Create a custom sorting function for time columns
  sort_time_cols <- function(cols) {
    # Extract day numbers
    day_nums <- as.numeric(str_extract(cols, "-?\\d+"))
    # Extract metric type
    metrics <- str_extract(cols, "total|mean|max|median")
    # Create a factor for metrics to control their order
    metric_factor <- factor(metrics, levels = c("total", "mean", "max", "median"))
    # Order by day, then metric
    cols[order(day_nums, metric_factor)]
  }
  
  sorted_time_cols <- sort_time_cols(time_cols)
  result <- final_result[, c(col_order, sorted_time_cols)]
  
  # Validation
  cat("\nFinal data validation:\n")
  cat("Number of unique subjects: ", nrow(result), "\n")
  cat("Number of columns: ", ncol(result), "\n")
  
  return(result)
}


# Function to supplement steps data with heart rate data
preprocess_steps_data <- function(steps_data, heart_rate_data) {
  # 提取心率数据
  hr_sample_info <- heart_rate_data@sample_info
  hr_values <- as.numeric(heart_rate_data@expression_data[1, ])
  
  # 创建包含心率值的数据框
  hr_df <- data.frame(
    sample_id = colnames(heart_rate_data@expression_data),
    heart_rate = hr_values
  ) %>%
    left_join(hr_sample_info, by = "sample_id") %>%
    filter(!is.na(heart_rate)) %>%
    dplyr::select(subject_id, measure_time) %>%
    distinct()
  
  # 提取步数数据
  steps_sample_info <- steps_data@sample_info
  
  # 确保expression_data是data.frame
  if(is.matrix(steps_data@expression_data)) {
    steps_expr <- as.data.frame(steps_data@expression_data)
  } else {
    steps_expr <- steps_data@expression_data
  }
  
  steps_values <- as.numeric(steps_expr["steps", ])
  
  # 创建包含步数值的数据框
  steps_df <- data.frame(
    sample_id = colnames(steps_expr),
    steps = steps_values
  ) %>%
    left_join(steps_sample_info, by = "sample_id")
  
  # 用步数数据连接心率时间戳
  merged_data <- hr_df %>%
    left_join(
      steps_df,
      by = c("subject_id", "measure_time")
    ) %>%
    mutate(
      # 当心率存在时，用0填充NA步数
      steps = ifelse(is.na(steps), 0, steps),
      # 创建或更新sample_id
      sample_id = ifelse(is.na(sample_id), 
                         paste0(subject_id, "_", format(measure_time, "%Y%m%d%H%M%S")),
                         sample_id)
    )
  
  # 确保有activity和class列
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
  
  # 创建更新的sample_info
  required_columns <- names(steps_sample_info)
  missing_columns <- setdiff(required_columns, names(merged_data))
  
  if(length(missing_columns) > 0) {
    merged_data[, missing_columns] <- NA
  }
  
  # 预先确保steps列存在于merged_data中
  if(!"steps" %in% names(merged_data)) {
    merged_data$steps <- 0
  }
  
  # 选择所有必要的列并确保唯一性
  updated_sample_info <- merged_data %>%
    dplyr::select(all_of(required_columns), steps) %>%  # 临时保留steps列
    distinct(sample_id, .keep_all = TRUE)
  
  # 创建更新的expression_data (作为data.frame)
  row_names <- rownames(steps_expr)
  
  # 初始化一个空的数据框
  new_expression_data <- as.data.frame(matrix(NA, 
                                              nrow = length(row_names), 
                                              ncol = nrow(updated_sample_info)))
  
  # 设置行名和列名
  rownames(new_expression_data) <- row_names
  colnames(new_expression_data) <- updated_sample_info$sample_id
  
  # 找到steps行的索引
  steps_row_index <- which(rownames(new_expression_data) == "steps")
  
  # 向量化赋值：直接使用updated_sample_info中的steps列
  if(length(steps_row_index) > 0) {
    new_expression_data[steps_row_index, ] <- updated_sample_info$steps
  }
  
  # 移除临时保留的steps列
  updated_sample_info <- updated_sample_info %>% 
    dplyr::select(all_of(required_columns))
  
  # 创建更新后的步数数据对象
  updated_steps_data <- steps_data
  updated_steps_data@expression_data <- new_expression_data
  updated_steps_data@sample_info <- updated_sample_info
  
  # 打印摘要和验证信息
  cat("Original number of steps data points:", nrow(steps_sample_info), "\n")
  cat("After filling with heart rate data:", nrow(updated_sample_info), "\n")
  cat("Number of added zero steps:", nrow(updated_sample_info) - nrow(steps_sample_info), "\n")
  cat("Number of columns in expression_data:", ncol(new_expression_data), "\n")
  cat("Number of rows in sample_info:", nrow(updated_sample_info), "\n")
  cat("Sample IDs in expression_data match sample_info:", 
      all(colnames(new_expression_data) %in% updated_sample_info$sample_id), "\n")
  
  return(updated_steps_data)
}
# Preprocess steps data before applying thresholds
preprocessed_steps_data <- preprocess_steps_data(daily_workout_details_data, heart_rate_data)

# Now apply thresholding to the preprocessed data
thresholded_steps_result <- calculate_daily_steps_with_thresholds(preprocessed_steps_data, baseline_info)

# Save the result
save(thresholded_steps_result, file = "thresholded_steps_result.rda", compress = "xz")

# Also save the preprocessed data for future reference
save(preprocessed_steps_data, file = "preprocessed_steps_data.rda", compress = "xz")

# Print summary of the preprocessing
cat("Preprocessing summary:\n")
cat("Original steps data points:", nrow(daily_workout_details_data@sample_info), "\n")
cat("After filling with heart rate data:", nrow(preprocessed_steps_data@sample_info), "\n")
cat("Number of added zero steps:", 
    nrow(preprocessed_steps_data@sample_info) - nrow(daily_workout_details_data@sample_info), 
    "\n")

# Count valid entries by threshold type
threshold_count_summary <- names(thresholded_steps_result) %>%
  tibble(column_name = .) %>%
  filter(grepl("daily_steps", column_name)) %>%
  mutate(
    day = as.numeric(str_extract(column_name, "-?\\d+")),
    threshold = str_extract(column_name, "time_[^_]+|day_[^_]+|season_[^_]+|all|time_[^_]+_[^_]+"),
    metric = str_extract(column_name, "total|mean|max|median"),
    valid_count = sapply(thresholded_steps_result[column_name], function(x) sum(!is.na(x)))
  )

# Summarize by threshold type
threshold_summary <- threshold_count_summary %>%
  group_by(threshold, metric) %>%
  summarise(
    avg_valid_count = mean(valid_count, na.rm = TRUE),
    min_valid_count = min(valid_count, na.rm = TRUE),
    max_valid_count = max(valid_count, na.rm = TRUE),
    .groups = "drop"
  )

print(threshold_summary)
