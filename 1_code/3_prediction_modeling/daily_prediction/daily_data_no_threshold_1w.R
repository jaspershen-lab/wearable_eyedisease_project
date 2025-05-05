library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)
library(caret)
library(randomForest)
library(mice)

# Step 2: Load the saved data (adjust paths as needed)
# Load heart rate data directly
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")

# Load the thresholded result data
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/daily_rhr_result.rda")
load("3_data_analysis/2_data_analysis/bo/bo_time_periods/daily_bo_result.rda")
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/daily_steps_result_assigned.rda")
load("3_data_analysis/2_data_analysis/sleep/sleep_time_periods/daily_sleep_result.rda")

# 加载心率数据结果 - 新增行
load("3_data_analysis/2_data_analysis/hr/hr_time_periods/daily_hr_result.rda")

# Load baseline info for disease and vision data
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Load OCTA data - added from second code
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness.csv")

# Step 3: Create disease and vision datasets from baseline_info
# Create disease dataset
disease_data <- baseline_info %>%
  mutate(
    surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d"),
    month = month(surgery_time_1),
    season = case_when(
      month %in% c(12, 1, 2) ~ "winter",
      month %in% c(3, 4, 5) ~ "spring",
      month %in% c(6, 7, 8, 9) ~ "summer",
      month %in% c(10, 11) ~ "fall"
    ),
    season_factor = factor(season, levels = c("spring", "summer", "fall", "winter")),
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
  dplyr::select(ID, cataract_2, dm_2, hypertension_2, season, season_factor, month, bmi)

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

# Process OCTA blood flow data for baseline (T0)
octa_bloodflow_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_bloodflow, by = c("ID" = "id"))

# Function to process blood flow data for each patient at baseline (T0)
process_patient_bloodflow <- function(patient_data, time_points = c("T0")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # Use left eye data for left eye surgery, right eye data for right eye and both eyes surgery
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # Process each time point
  for(suffix in time_points) {
    # Select columns for current time point
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # Select data and rename columns
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# Process each patient individually for baseline (T0)
patient_list_bloodflow <- split(octa_bloodflow_features, octa_bloodflow_features$ID)
processed_data_bloodflow <- purrr::map(patient_list_bloodflow, process_patient_bloodflow)

# Combine results
octa_bloodflow_features <- bind_rows(processed_data_bloodflow)

# Create blood flow variables subset for baseline (T0)
bloodflow_var_T0 <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T0$"),  # Select all columns ending with 0_6_T0
    -matches("PA_OuterRetina_0_6_T0"),  # Exclude these columns
    -matches("PA_PED_0_6_T0")
  )

# Process OCTA thickness data for baseline (T0)
octa_thickness_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_thickness, by = c("ID" = "id"))

# Function to process thickness data for each patient at baseline (T0)
process_patient_thickness <- function(patient_data, time_points = c("T0")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # Use left eye data for left eye surgery, right eye data for right eye and both eyes surgery
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # Process each time point
  for(suffix in time_points) {
    # Select columns for current time point
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # Select data and rename columns
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# Process each patient individually for baseline (T0)
patient_list_thickness <- split(octa_thickness_features, octa_thickness_features$ID)
processed_data_thickness <- purrr::map(patient_list_thickness, process_patient_thickness)

# Combine results
octa_thickness_features <- bind_rows(processed_data_thickness)

# Create thickness variables subset for baseline (T0)
thickness_var_T0 <- octa_thickness_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T0$"),
    matches("Thickness_PED_0_6_T0")
  )

# Step 4: Create the daily coverage status dataset
# 修改后的覆盖率分析函数 - 仅要求全天总计≥8小时
analyze_time_period_coverage <- function(heart_rate_data, baseline_info) {
  # Use the already loaded heart_rate_data object
  
  # Process baseline information
  baseline_info <- baseline_info %>%
    mutate(
      surgery_time_1 = as.Date(surgery_time_1)
    )
  
  # Get sample info from heart rate data
  hr_df <- heart_rate_data@sample_info %>%
    as.data.frame()
  
  # Calculate hourly coverage by day (without separating day/night)
  time_period_coverage <- hr_df %>%
    # Join with surgery dates
    left_join(baseline_info %>% dplyr::select(ID, surgery_time_1), by = c("subject_id" = "ID")) %>%
    # Calculate days relative to surgery
    mutate(
      days_to_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days")),
      day_point = floor(days_to_surgery),
      hour = lubridate::hour(measure_time)
    ) %>%
    # Filter to our desired range (-7 to 6 days)
    filter(
      day_point >= -7,
      day_point <= 6
    ) %>%
    # Count distinct hours per subject-day
    group_by(subject_id, day_point) %>%
    summarise(
      hours_covered = n_distinct(hour),
      .groups = "drop"
    )
  
  # Create complete grid with all subject-day combinations
  all_subjects <- unique(hr_df$subject_id)
  all_days <- seq(-7, 6)
  complete_grid <- expand.grid(
    subject_id = all_subjects,
    day_point = all_days,
    stringsAsFactors = FALSE
  )
  
  # Join with actual coverage and fill missing with 0
  final_coverage <- complete_grid %>%
    left_join(time_period_coverage, by = c("subject_id", "day_point")) %>%
    mutate(hours_covered = coalesce(hours_covered, 0))
  
  # Calculate coverage status for each day (≥8 hours total)
  daily_coverage_status <- final_coverage %>%
    # Mark if day meets threshold (>= 8 hours)
    mutate(
      meets_threshold = hours_covered >= 8
    )
  
  return(daily_coverage_status)
}

# Calculate daily coverage status
# 计算每日覆盖状态 - 使用修改后的函数
daily_coverage_status <- analyze_time_period_coverage(
  heart_rate_data = heart_rate_data,
  baseline_info = baseline_info
)

# 打印每日样本量信息以进行验证
daily_coverage_count <- daily_coverage_status %>%
  group_by(day_point) %>%
  summarise(
    total_participants = n(),
    participants_with_enough_coverage = sum(meets_threshold),
    .groups = "drop"
  )

# 打印覆盖率信息
print(daily_coverage_count)

# Step 5: Filter participants based on the >6h day and night criteria
# For each day from -7 to 6, get the IDs of participants with adequate coverage
# 修改筛选参与者函数以适应新的覆盖率要求
filtered_participants_by_day <- function(daily_coverage_status, target_day) {
  daily_coverage_status %>%
    filter(day_point == target_day, meets_threshold == TRUE) %>%
    pull(subject_id)
}

# Step 6: Extract metrics for each day and combine with baseline data
# Modified to include HR data and OCTA data
prepare_day_dataset <- function(day_point, 
                                rhr_data = daily_rhr_result, 
                                bo_data = daily_bo_result,
                                steps_data = daily_steps_result,
                                sleep_data = daily_sleep_result,
                                hr_data = daily_hr_result,  # 添加心率数据参数
                                daily_coverage_status,
                                disease_data,
                                vision_data,
                                bloodflow_var_T0,
                                thickness_var_T0) {
  
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
  
  # Extract Sleep metrics for the day
  sleep_columns <- names(sleep_data) %>%
    grep(paste0("day_", day_point, "_"), ., value = TRUE)
  
  if(length(sleep_columns) > 0) {
    sleep_day_data <- sleep_data %>%
      dplyr::select(subject_id, all_of(sleep_columns)) %>%
      # Rename columns to remove day prefix for clarity
      rename_with(~ gsub(paste0("day_", day_point, "_"), "", .), starts_with(paste0("day_", day_point, "_")))
  } else {
    sleep_day_data <- data.frame(subject_id = rhr_day_data$subject_id)
  }
  
  # 提取心率指标 - 新增部分
  hr_columns <- names(hr_data) %>%
    grep(paste0("day_", day_point, "_"), ., value = TRUE)
  
  if(length(hr_columns) > 0) {
    hr_day_data <- hr_data %>%
      dplyr::select(subject_id, all_of(hr_columns)) %>%
      # 重命名列以去除日期前缀，提高清晰度
      rename_with(~ gsub(paste0("day_", day_point, "_"), "", .), starts_with(paste0("day_", day_point, "_")))
  } else {
    hr_day_data <- data.frame(subject_id = rhr_day_data$subject_id)
  }
  
  # Combine all data
  combined_data <- rhr_day_data %>%
    left_join(bo_day_data, by = "subject_id") %>%
    left_join(steps_day_data, by = "subject_id") %>%
    left_join(sleep_day_data, by = "subject_id") %>%
    left_join(hr_day_data, by = "subject_id") %>%  # 加入心率数据 - 新增行
    # Join with disease data
    left_join(disease_data, by = c("subject_id" = "ID")) %>%
    # Join with vision data
    left_join(vision_data, by = c("subject_id" = "ID")) %>%
    # Join with OCTA baseline bloodflow data (T0)
    left_join(bloodflow_var_T0, by = c("subject_id" = "ID")) %>%
    # Join with OCTA baseline thickness data (T0)
    left_join(thickness_var_T0, by = c("subject_id" = "ID")) %>%
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
    rhr_data = daily_rhr_result,
    bo_data = daily_bo_result,
    steps_data = daily_steps_result,
    sleep_data = daily_sleep_result,
    hr_data = daily_hr_result,  # 添加心率数据 - 新增行
    daily_coverage_status = daily_coverage_status,
    disease_data = disease_data,
    vision_data = vision_data,
    bloodflow_var_T0 = bloodflow_var_T0,
    thickness_var_T0 = thickness_var_T0
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
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data"
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

# 分析特征类型分布 - 新增部分
if(length(valid_dfs) > 0 && nrow(combined_data) > 0) {
  # 查看每种数据类型有多少特征
  feature_types <- sapply(combined_data, function(col) class(col)[1])
  feature_type_counts <- table(feature_types)
  cat("\n特征类型分布:\n")
  print(feature_type_counts)
  
  # 分析心率特征
  hr_features <- grep("_hr_|heart_rate", names(combined_data), value = TRUE)
  cat("\n包含的心率特征 (", length(hr_features), "):\n")
  print(hr_features[1:min(20, length(hr_features))])
  if(length(hr_features) > 20) {
    cat("...及", length(hr_features) - 20, "个更多特征\n")
  }
}

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



# 检查SH038的覆盖情况
daily_coverage_status %>% 
  filter(subject_id == "SH038") %>%
  arrange(day_point)
