library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)
library(caret)
library(randomForest)
library(mice)



# Data Preparation for Vision Improvement Prediction
# This script prepares data for predicting post-surgery vision improvement
# using wearable device metrics, baseline disease data, and pre-surgery vision

# Step 1: Load required libraries
library(tidyverse)
library(caret)  # For modeling later

# Step 2: Load the saved data (adjust paths as needed)
# Load the thresholded result data
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/thresholded_rhr_result.rda")
load("3_data_analysis/2_data_analysis/bo/bo_time_periods/thresholded_bo_result.rda")
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/thresholded_steps_result.rda")

# Load baseline info for disease and vision data
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Step 3: Create disease and vision datasets from baseline_info
# Create disease dataset
disease_data <- baseline_info %>%
  mutate(
    cataract_2 = case_when(
      cataract == 1 ~ 0,  
      cataract %in% c(2, 3, 4) ~ 1, 
      TRUE ~ NA_real_
    ),
    dm_2 = case_when(
      diabetes_history == 1 ~ 1,  
      diabetes_history == 2 ~ 0, 
      TRUE ~ NA_real_
    ),
    hypertension_2 = case_when(
      hypertension_history == 1 ~ 1,  
      hypertension_history == 2 ~ 0, 
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(ID, cataract_2, dm_2, hypertension_2)

# Create vision dataset
vision_data <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,  # Right eye surgery
      surgery_eye_1 == 1 ~ os_corrected_bas,  # Left eye surgery
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,  # Both eyes (average)
      TRUE ~ NA_real_
    ),
    post_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1w,   # Right eye post-surgery
      surgery_eye_1 == 1 ~ os_corrected_1w,   # Left eye post-surgery
      surgery_eye_1 == 2 ~ (od_corrected_1w + os_corrected_1w)/2,   # Both eyes post-surgery
      TRUE ~ NA_real_
    ),
    vision_improvement = post_vision - pre_vision,
    vision_improved = if_else(vision_improvement >= 0, 1, 0),  # Binary improvement indicator
    vision_improved_factor = factor(vision_improved, 
                                    levels = c(0, 1), 
                                    labels = c("NoImprovement", "Improved"))  # Factor version
  ) %>%
  dplyr::select(ID, surgery_eye_1, pre_vision, post_vision, 
                vision_improvement, vision_improved, vision_improved_factor,
                age, gender)

# Step 4: Create the daily coverage status dataset
# This identifies participants with >6h of data for both day and night
analyze_time_period_coverage <- function(heart_rate_data, baseline_info) {
  # Get heart rate data from RDA file
  load(heart_rate_data)  # This loads the heart_rate_data object
  
  # Process baseline information
  baseline_info <- baseline_info %>%
    mutate(
      surgery_time_1 = as.Date(surgery_time_1)
    )
  
  # Get sample info from heart rate data
  hr_df <- heart_rate_data@sample_info %>%
    as.data.frame()
  
  # Calculate hourly coverage
  time_period_coverage <- hr_df %>%
    # Join with surgery dates
    left_join(baseline_info %>% dplyr::select(ID, surgery_time_1), by = c("subject_id" = "ID")) %>%
    # Calculate days relative to surgery
    mutate(
      days_to_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days")),
      day_point = floor(days_to_surgery),
      hour = lubridate::hour(measure_time),
      # Define time periods
      time_period = case_when(
        hour >= 6 & hour < 18 ~ "daytime",
        TRUE ~ "nighttime"  # Hours 18-23 and 0-5
      )
    ) %>%
    # Filter to our desired range (-7 to 6 days)
    filter(
      day_point >= -7,
      day_point <= 6
    ) %>%
    # Count distinct hours per subject-day-period
    group_by(subject_id, day_point, time_period) %>%
    summarise(
      hours_covered = n_distinct(hour),
      .groups = "drop"
    )
  
  # Create complete grid with all subject-day-period combinations
  all_subjects <- unique(hr_df$subject_id)
  all_days <- seq(-7, 6)
  all_periods <- c("daytime", "nighttime")
  complete_grid <- expand.grid(
    subject_id = all_subjects,
    day_point = all_days,
    time_period = all_periods,
    stringsAsFactors = FALSE
  )
  
  # Join with actual coverage and fill missing with 0
  final_coverage <- complete_grid %>%
    left_join(time_period_coverage, by = c("subject_id", "day_point", "time_period")) %>%
    mutate(hours_covered = coalesce(hours_covered, 0))
  
  # Calculate coverage status for each day
  daily_coverage_status <- final_coverage %>%
    pivot_wider(
      id_cols = c(subject_id, day_point),
      names_from = time_period,
      values_from = hours_covered
    ) %>%
    # Mark if both time periods meet threshold (>= 6 hours)
    mutate(
      meets_daytime_threshold = daytime >= 6,
      meets_nighttime_threshold = nighttime >= 6,
      meets_both_thresholds = meets_daytime_threshold & meets_nighttime_threshold
    )
  
  return(daily_coverage_status)
}

# Calculate daily coverage status
daily_coverage_status <- analyze_time_period_coverage(
  heart_rate_data = "3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda",
  baseline_info = baseline_info
)

# Step 5: Filter participants based on the >6h day and night criteria
# For each day from -7 to 6, get the IDs of participants with adequate coverage
filtered_participants_by_day <- function(daily_coverage_status, target_day) {
  daily_coverage_status %>%
    filter(day_point == target_day, meets_both_thresholds == TRUE) %>%
    pull(subject_id)
}

# Step 6: Extract metrics for each day and combine with baseline data
prepare_day_dataset <- function(day_point, 
                                rhr_data = thresholded_rhr_result, 
                                bo_data = thresholded_bo_result,
                                steps_data = thresholded_steps_result,
                                daily_coverage_status,
                                disease_data,
                                vision_data) {
  
  # Get participants with adequate coverage for this day
  valid_participants <- filtered_participants_by_day(daily_coverage_status, day_point)
  
  # Extract RHR metrics for the day
  rhr_columns <- names(rhr_data) %>%
    grep(paste0("day_", day_point, "_"), ., value = TRUE)
  
  rhr_day_data <- rhr_data %>%
    dplyr::select(subject_id, all_of(rhr_columns)) %>%
    # Rename columns to remove day prefix for clarity
    rename_with(~ gsub(paste0("day_", day_point, "_"), "", .), starts_with(paste0("day_", day_point, "_")))
  
  # Extract BO metrics for the day
  bo_columns <- names(bo_data) %>%
    grep(paste0("day_", day_point, "_"), ., value = TRUE)
  
  bo_day_data <- bo_data %>%
    dplyr::select(subject_id, all_of(bo_columns)) %>%
    # Rename columns to remove day prefix for clarity
    rename_with(~ gsub(paste0("day_", day_point, "_"), "", .), starts_with(paste0("day_", day_point, "_")))
  
  # Extract Steps metrics for the day
  steps_columns <- names(steps_data) %>%
    grep(paste0("day_", day_point, "_"), ., value = TRUE)
  
  steps_day_data <- steps_data %>%
    dplyr::select(subject_id, all_of(steps_columns)) %>%
    # Rename columns to remove day prefix for clarity
    rename_with(~ gsub(paste0("day_", day_point, "_"), "", .), starts_with(paste0("day_", day_point, "_")))
  
  # Combine all data
  combined_data <- rhr_day_data %>%
    left_join(bo_day_data, by = "subject_id") %>%
    left_join(steps_day_data, by = "subject_id") %>%
    # Join with disease data
    left_join(disease_data, by = c("subject_id" = "ID")) %>%
    # Join with vision data
    left_join(vision_data, by = c("subject_id" = "ID")) %>%
    # Filter to only include participants with adequate coverage
    filter(subject_id %in% valid_participants)
  
  # Add day information
  combined_data$day_relative_to_surgery <- day_point
  
  return(combined_data)
}

# Step 7: Generate datasets for each day
all_days <- seq(-7, 6)
day_datasets <- list()

for(day in all_days) {
  day_datasets[[as.character(day)]] <- prepare_day_dataset(
    day_point = day,
    rhr_data = thresholded_rhr_result,
    bo_data = thresholded_bo_result,
    steps_data = thresholded_steps_result,
    daily_coverage_status = daily_coverage_status,
    disease_data = disease_data,
    vision_data = vision_data
  )
}

# Step 8: Analyze sample sizes for each day
sample_sizes <- sapply(all_days, function(day) {
  nrow(day_datasets[[as.character(day)]])
})

# Create a data frame for easy viewing
day_sample_sizes <- data.frame(
  day_relative_to_surgery = all_days,
  sample_size = sample_sizes
)

# Print the sample sizes
print(day_sample_sizes)

# Step 9: Function to prepare data for modeling
prepare_for_modeling <- function(dataset, outcome_type = "continuous") {
  # Remove identifier columns
  modeling_data <- dataset %>%
    dplyr::select(-subject_id, -surgery_eye_1)
  
  # Define outcome variable based on modeling type
  if(outcome_type == "continuous") {
    # For continuous outcome (vision_improvement)
    outcome_col <- "vision_improvement"
  } else if(outcome_type == "binary") {
    # For binary classification (vision_improved)
    outcome_col <- "vision_improved"
  } else if(outcome_type == "factor") {
    # For factor outcome (vision_improved_factor)
    outcome_col <- "vision_improved_factor"
  }
  
  # Handle missing values - several approaches possible:
  # 1. Remove rows with any missing values
  # complete_data <- na.omit(modeling_data)
  
  # 2. Impute missing values - example with median imputation
  # Only impute predictors, not the outcome
  predictor_cols <- setdiff(names(modeling_data), 
                            c("vision_improvement", "vision_improved", "vision_improved_factor"))
  
  # Check if we have enough data for imputation
  if(nrow(modeling_data) > 5) {  # Arbitrary threshold
    preprocess_steps <- preProcess(modeling_data[, predictor_cols], method = c("medianImpute"))
    imputed_predictors <- predict(preprocess_steps, modeling_data[, predictor_cols])
    
    # Combine imputed predictors with outcome
    imputed_data <- cbind(imputed_predictors, modeling_data[, c(outcome_col)])
    names(imputed_data)[ncol(imputed_data)] <- outcome_col
    
    return(imputed_data)
  } else {
    # If too few observations, just return the original data
    warning("Too few observations for imputation. Returning original data.")
    return(modeling_data)
  }
}

# Step 10: Export each day's dataset to CSV for further analysis
# Create directory if it doesn't exist
output_dir <- "3_data_analysis/3_prediction_modeling/daily_prediction/daily_data/1w"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 首先检查每个day_datasets元素的结构
for(day in all_days) {
  day_str <- as.character(day)
  
  # 跳过空数据集
  if(length(day_datasets[[day_str]]) == 0) {
    cat("No data available for day", day_str, "\n\n")
    next
  }
  
  # 检查是否是数据框
  if(!is.data.frame(day_datasets[[day_str]])) {
    cat("Warning: Data for day", day_str, "is not a data frame\n")
    cat("  Class:", class(day_datasets[[day_str]]), "\n")
    
    # 如果是列表但不是数据框，尝试转换为数据框
    if(is.list(day_datasets[[day_str]])) {
      tryCatch({
        # 尝试将列表转换为数据框
        day_datasets[[day_str]] <- as.data.frame(day_datasets[[day_str]])
        cat("  Successfully converted to data frame\n")
      }, error = function(e) {
        cat("  Failed to convert to data frame:", e$message, "\n")
      })
    }
  }
  
  # 再次检查是否现在是数据框
  if(is.data.frame(day_datasets[[day_str]])) {
    # 检查数据框中是否有列表列
    list_cols <- sapply(day_datasets[[day_str]], is.list)
    if(any(list_cols)) {
      cat("  Warning: Data frame for day", day_str, "contains list columns\n")
      cat("  List columns:", names(list_cols[list_cols]), "\n")
      
      # 删除所有列表类型的列
      day_datasets[[day_str]] <- day_datasets[[day_str]][, !list_cols]
      cat("  Removed list columns for CSV export\n")
    }
    
    # 导出为CSV
    filename <- file.path(output_dir, paste0("day_", day_str, "_data.csv"))
    tryCatch({
      write.csv(day_datasets[[day_str]], file = filename, row.names = FALSE)
      cat("Exported data for day", day_str, "to", filename, "\n")
      cat("  Number of participants:", nrow(day_datasets[[day_str]]), "\n")
      cat("  Number of variables:", ncol(day_datasets[[day_str]]), "\n\n")
    }, error = function(e) {
      cat("  Failed to export data for day", day_str, ":", e$message, "\n\n")
    })
  } else {
    cat("  Unable to export data for day", day_str, "- not a proper data frame\n\n")
  }
}

# 保存样本量信息到CSV
sample_sizes_file <- file.path(output_dir, "day_sample_sizes.csv")
write.csv(day_sample_sizes, file = sample_sizes_file, row.names = FALSE)
cat("Exported sample sizes to", sample_sizes_file, "\n\n")

# 尝试合并数据 - 只合并数据框类型的元素
valid_dfs <- list()
for(day in all_days) {
  day_str <- as.character(day)
  if(length(day_datasets[[day_str]]) > 0 && is.data.frame(day_datasets[[day_str]])) {
    valid_dfs[[day_str]] <- day_datasets[[day_str]]
  }
}

if(length(valid_dfs) > 0) {
  # 尝试合并数据
  tryCatch({
    combined_data <- bind_rows(valid_dfs, .id = "day")
    combined_file <- file.path(output_dir, "all_days_combined_data.csv")
    write.csv(combined_data, file = combined_file, row.names = FALSE)
    cat("Exported combined data for all days to", combined_file, "\n")
    cat("  Total observations:", nrow(combined_data), "\n")
    cat("  Number of variables:", ncol(combined_data), "\n\n")
  }, error = function(e) {
    cat("Failed to combine datasets:", e$message, "\n\n")
  })
} else {
  cat("No valid data frames available to combine\n\n")
}

# 保存R数据对象 - 检查后的版本
save(day_datasets, file = file.path(output_dir, "vision_prediction_day_datasets.rda"))
save(day_sample_sizes, file = file.path(output_dir, "vision_prediction_sample_sizes.rda"))

cat("Data preparation complete. Datasets saved in", output_dir, "\n")



# 检查每个day_datasets元素的结构
for(day in names(day_datasets)) {
  cat("Day", day, "- Class:", class(day_datasets[[day]]), "\n")
  if(is.data.frame(day_datasets[[day]])) {
    cat("  Is data frame: YES\n")
    cat("  Dimensions:", dim(day_datasets[[day]]), "\n")
    
    # 检查是否有列表列
    list_cols <- which(sapply(day_datasets[[day]], is.list))
    if(length(list_cols) > 0) {
      cat("  Has list columns: YES -", names(day_datasets[[day]])[list_cols], "\n")
    } else {
      cat("  Has list columns: NO\n")
    }
  } else {
    cat("  Is data frame: NO\n")
  }
  cat("\n")
}


# 改进的筛选函数，只保留time_daytime、time_nighttime和all阈值的变量
filter_specific_threshold_variables <- function(datasets) {
  filtered_datasets <- list()
  
  for(day in names(datasets)) {
    # 获取当前数据集
    current_df <- datasets[[day]]
    
    # 如果是空的或不是数据框，则跳过
    if(length(current_df) == 0 || !is.data.frame(current_df)) {
      cat("跳过Day", day, "- 不是有效的数据框\n")
      next
    }
    
    # 获取所有变量名
    all_vars <- names(current_df)
    
    # 定义要保留的变量类型
    # 1. 基本标识符和结果变量
    id_vars <- c("subject_id", "surgery_eye_1", "day_relative_to_surgery", 
                 "pre_vision", "post_vision", "vision_improvement", 
                 "vision_improved", "vision_improved_factor", "age", "gender")
    
    # 2. 疾病变量
    disease_vars <- c("cataract_2", "dm_2", "hypertension_2")
    
    # 3. 更精确地筛选变量
    # 只保留以下模式的变量:
    # a. 包含"_all"但不包含"season"、"weekend"或"weekday"的变量
    all_threshold_vars <- grep("_all", all_vars, value = TRUE)
    all_threshold_vars <- all_threshold_vars[!grepl("season|weekend|weekday", all_threshold_vars)]
    
    # b. 只包含"time_daytime$"或"time_nighttime$"的变量（确保不带weekend/weekday后缀）
    daytime_vars <- grep("time_daytime$", all_vars, value = TRUE)
    nighttime_vars <- grep("time_nighttime$", all_vars, value = TRUE)
    
    # 合并所有筛选出的阈值变量
    threshold_vars <- c(all_threshold_vars, daytime_vars, nighttime_vars)
    
    # 合并所有要保留的变量
    keep_vars <- c(id_vars, disease_vars, threshold_vars)
    
    # 只保留在当前数据框中实际存在的变量
    keep_vars <- keep_vars[keep_vars %in% all_vars]
    
    # 筛选数据框
    filtered_df <- current_df[, keep_vars, drop = FALSE]
    
    # 添加到结果列表
    filtered_datasets[[day]] <- filtered_df
    
    # 打印结果
    cat("Day", day, "- 原始变量数:", length(all_vars), 
        "筛选后变量数:", ncol(filtered_df), "\n")
    
    # 检查是否仍有season或weekend变量
    remaining_season <- sum(grepl("season", names(filtered_df)))
    remaining_weekend <- sum(grepl("weekend|weekday", names(filtered_df)))
    cat("  剩余season变量:", remaining_season, "个\n")
    cat("  剩余weekend/weekday变量:", remaining_weekend, "个\n")
  }
  
  return(filtered_datasets)
}

# 应用改进的筛选函数
improved_filtered_day_datasets <- filter_specific_threshold_variables(day_datasets)

# 导出筛选后的数据集
output_dir <- "3_data_analysis/3_prediction_modeling/daily_prediction/daily_data_filtered"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 导出每个筛选后的数据集
for(day in names(improved_filtered_day_datasets)) {
  output_file <- file.path(output_dir, paste0("day_", day, "_filtered_data.csv"))
  
  tryCatch({
    write.csv(improved_filtered_day_datasets[[day]], file = output_file, row.names = FALSE)
    cat("成功导出筛选后的Day", day, "数据到", output_file, "\n")
    
    # 打印筛选后的维度
    cat("  维度:", nrow(improved_filtered_day_datasets[[day]]), "行 ×", 
        ncol(improved_filtered_day_datasets[[day]]), "列\n")
  }, error = function(e) {
    cat("导出筛选后的Day", day, "数据失败:", e$message, "\n")
  })
}

# 导出合并的筛选后数据集
tryCatch({
  improved_filtered_combined_data <- bind_rows(improved_filtered_day_datasets, .id = "day")
  improved_filtered_combined_file <- file.path(output_dir, "all_days_filtered_combined_data.csv")
  write.csv(improved_filtered_combined_data, file = improved_filtered_combined_file, row.names = FALSE)
  cat("成功导出筛选后的合并数据到", improved_filtered_combined_file, "\n")
  cat("  维度:", nrow(improved_filtered_combined_data), "行 ×", 
      ncol(improved_filtered_combined_data), "列\n")
}, error = function(e) {
  cat("合并筛选后的数据集失败:", e$message, "\n")
})
