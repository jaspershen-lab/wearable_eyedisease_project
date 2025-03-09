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
# Load the thresholded result data with time periods
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/time_period_rhr_results.rda")
load("3_data_analysis/2_data_analysis/bo/bo_time_periods/bo_time_period_results.rda")
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/steps_time_period_results_assigned.rda")
load("3_data_analysis/2_data_analysis/sleep/sleep_time_periods/sleep_time_period_results.rda")
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")

# Load baseline info for disease and vision data
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Load OCTA data
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
  dplyr::select(ID, cataract_2, dm_2, hypertension_2, season, season_factor, month)

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
    post_vision_1w = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1w,   # Right eye post-surgery 1 week
      surgery_eye_1 == 1 ~ os_corrected_1w,   # Left eye post-surgery 1 week
      surgery_eye_1 == 2 ~ (od_corrected_1w + os_corrected_1w)/2,   # Both eyes post-surgery 1 week
      TRUE ~ NA_real_
    ),
    post_vision_1m = case_when(
      surgery_eye_1 == 0 ~ od_corrected_1m,   # Right eye post-surgery 1 month
      surgery_eye_1 == 1 ~ os_corrected_1m,   # Left eye post-surgery 1 month
      surgery_eye_1 == 2 ~ (od_corrected_1m + os_corrected_1m)/2,   # Both eyes post-surgery 1 month
      TRUE ~ NA_real_
    ),
    vision_improvement_1w = post_vision_1w - pre_vision,
    vision_improvement_1m = post_vision_1m - pre_vision,
    vision_improved_1m = if_else(vision_improvement_1m >= 0, 1, 0),  # Binary improvement indicator
    vision_improved_factor_1m = factor(vision_improved_1m, 
                                       levels = c(0, 1), 
                                       labels = c("NoImprovement", "Improved"))  # Factor version
  ) %>%
  dplyr::select(ID, surgery_eye_1, pre_vision, post_vision_1w, post_vision_1m,
                vision_improvement_1w, vision_improvement_1m, 
                vision_improved_1m, vision_improved_factor_1m,
                age, gender)

# Process OCTA blood flow data for all time points (T0, T1, T2)
octa_bloodflow_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_bloodflow, by = c("ID" = "id"))

# Modified function to process blood flow data for each patient and time point
process_patient_bloodflow <- function(patient_data, time_points = c("T0", "T1", "T2")) {
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

# Process each patient individually for all time points
patient_list_bloodflow <- split(octa_bloodflow_features, octa_bloodflow_features$ID)
processed_data_bloodflow <- purrr::map(patient_list_bloodflow, process_patient_bloodflow)

# Combine results
octa_bloodflow_features <- bind_rows(processed_data_bloodflow)

# Create blood flow variables subset for each time point
bloodflow_var_T0 <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T0$"),  # Select all columns ending with 0_6_T0
    -matches("PA_OuterRetina_0_6_T0"),  # Exclude these columns
    -matches("PA_PED_0_6_T0")
  )

bloodflow_var_T1 <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T1$"),  # Select all columns ending with 0_6_T1
    -matches("PA_OuterRetina_0_6_T1"),  # Exclude these columns
    -matches("PA_PED_0_6_T1")
  )

bloodflow_var_T2 <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T2$"),  # Select all columns ending with 0_6_T2
    -matches("PA_OuterRetina_0_6_T2"),  # Exclude these columns
    -matches("PA_PED_0_6_T2")
  )

# Process OCTA thickness data for all time points
octa_thickness_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_thickness, by = c("ID" = "id"))

# Modified function to process thickness data for each patient and time point
process_patient_thickness <- function(patient_data, time_points = c("T0", "T1", "T2")) {
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

# Process each patient individually for all time points
patient_list_thickness <- split(octa_thickness_features, octa_thickness_features$ID)
processed_data_thickness <- purrr::map(patient_list_thickness, process_patient_thickness)

# Combine results
octa_thickness_features <- bind_rows(processed_data_thickness)

# Create thickness variables subset for each time point
thickness_var_T0 <- octa_thickness_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T0$"),
    matches("Thickness_PED_0_6_T0")
  )

thickness_var_T1 <- octa_thickness_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T1$"),
    matches("Thickness_PED_0_6_T1")
  )

thickness_var_T2 <- octa_thickness_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T2$"),
    matches("Thickness_PED_0_6_T2")
  )

# Step 4: Define time periods based on your case_when function
# 根据提供的时间窗口定义
time_periods <- list(
  "pre_surgery_3d" = c(-3, -1),            # 术前3天
  "pre_surgery_3d_to_7d" = c(-7, -4),      # 术前3-7天
  "pre_surgery_7d_all" = c(-7, -1),        # 术前全部7天
  "pre_surgery_all" = c(-30, -1),          # 术前所有数据(假设最多30天)
  "post_surgery_1w" = c(0, 7),             # 术后第1周
  "post_surgery_7d" = c(0, 7),             # 术后7天
  "post_surgery_7d_to_30d" = c(8, 30),     # 术后7-30天
  "post_surgery_day23_to_30" = c(23, 30),  # 术后23-30天
  "post_surgery_day27_to_30" = c(27, 30),  # 术后27-30天
  "post_surgery_over_30d" = c(31, 60)      # 术后>30天(假设最多60天)
)

# 函数：根据时间窗口长度计算覆盖率阈值
calculate_threshold <- function(period_range) {
  days_in_period <- period_range[2] - period_range[1] + 1
  if(days_in_period <= 3) {
    return(8)   # 1-3天的窗口：至少8个唯一小时
  } else if(days_in_period <= 7) {
    return(16)  # 4-7天的窗口：至少16个唯一小时
  } else {
    return(24)  # 8天以上的窗口：至少24个唯一小时
  }
}

# Step 5: Calculate time period coverage and determine which participants have adequate data
analyze_time_period_coverage <- function(heart_rate_data, baseline_info, time_periods) {
  # Get sample info from heart rate data
  hr_df <- heart_rate_data@sample_info %>%
    as.data.frame()
  
  # Process baseline information
  baseline_info <- baseline_info %>%
    mutate(
      surgery_time_1 = as.Date(surgery_time_1)
    )
  
  # Calculate time relative to surgery for each measurement
  time_data <- hr_df %>%
    # Join with surgery dates
    left_join(baseline_info %>% dplyr::select(ID, surgery_time_1), by = c("subject_id" = "ID")) %>%
    # Calculate days relative to surgery
    mutate(
      days_to_surgery = as.numeric(difftime(measure_time, surgery_time_1, units = "days")),
      day_point = floor(days_to_surgery),
      hour = lubridate::hour(measure_time)
    ) %>%
    # Filter to our desired range (all days covered by our time periods)
    filter(
      day_point >= min(sapply(time_periods, min)),
      day_point <= max(sapply(time_periods, max))
    )
  
  # Calculate coverage for each time period
  coverage_status <- list()
  
  for (period_name in names(time_periods)) {
    period_range <- time_periods[[period_name]]
    days_in_period <- period_range[2] - period_range[1] + 1
    
    # 计算该时间窗口的覆盖率阈值
    threshold <- calculate_threshold(period_range)
    
    # Filter data for current period
    period_data <- time_data %>%
      filter(day_point >= period_range[1], day_point <= period_range[2])
    
    # Group by subject and calculate coverage metrics
    period_coverage <- period_data %>%
      group_by(subject_id) %>%
      summarise(
        total_measurements = n(),
        unique_hours = n_distinct(paste(day_point, hour, sep = "_")),
        days_covered = n_distinct(day_point),
        time_span = max(day_point) - min(day_point) + 1,
        required_threshold = threshold,  # 添加所需阈值信息
        days_in_period = days_in_period, # 添加窗口天数信息
        .groups = "drop"
      ) %>%
      mutate(
        # 根据时间窗口长度使用不同阈值
        meets_threshold = unique_hours >= threshold
      )
    
    # Add to coverage status list
    coverage_status[[period_name]] <- period_coverage
  }
  
  # Calculate overall statistics for each time period
  coverage_summary <- lapply(names(time_periods), function(period) {
    data <- coverage_status[[period]]
    period_range <- time_periods[[period]]
    days_in_period <- period_range[2] - period_range[1] + 1
    threshold <- calculate_threshold(period_range)
    
    data.frame(
      period = period,
      days_in_period = days_in_period,
      required_threshold = threshold,
      total_subjects = nrow(data),
      subjects_with_coverage = sum(data$meets_threshold, na.rm = TRUE),
      avg_measurements = mean(data$total_measurements, na.rm = TRUE),
      avg_hours = mean(data$unique_hours, na.rm = TRUE),
      avg_days = mean(data$days_covered, na.rm = TRUE),
      coverage_rate = round(sum(data$meets_threshold, na.rm = TRUE) / nrow(data) * 100, 1)
    )
  }) %>% bind_rows()
  
  return(list(coverage_status = coverage_status, summary = coverage_summary))
}

# Calculate time period coverage
time_period_coverage <- analyze_time_period_coverage(
  heart_rate_data = heart_rate_data,
  baseline_info = baseline_info,
  time_periods = time_periods
)

# Extract coverage status and summary
time_period_coverage_status <- time_period_coverage$coverage_status
coverage_summary <- time_period_coverage$summary

# Print coverage summary
print("Coverage Summary by Time Period:")
print(coverage_summary)

# Step 6: Function to get IDs of participants with adequate coverage for a time period
filtered_participants_by_period <- function(coverage_status, period_name) {
  coverage_status[[period_name]] %>%
    filter(meets_threshold == TRUE) %>%
    pull(subject_id)
}

# Step 7: Extract time period metrics and combine with baseline data
# 函数来提取给定时间周期的列，考虑自定义的时间窗口命名
extract_period_columns <- function(data, period_name) {
  # 尝试直接匹配
  period_pattern <- paste0("^", period_name, "_")
  cols <- grep(period_pattern, names(data), value = TRUE)
  
  # 如果找不到匹配的列，尝试其他可能的命名模式
  if(length(cols) == 0) {
    # 检查名称格式差异或缩写
    simplified_pattern <- gsub("surgery_", "", period_name)
    simplified_pattern <- gsub("_to_", "_", simplified_pattern)
    cols <- grep(simplified_pattern, names(data), value = TRUE)
    
    # 检查更简单的匹配
    if(length(cols) == 0) {
      # 尝试匹配主关键词
      if(grepl("pre_surgery_3d$", period_name)) {
        cols <- grep("pre.*3d[^_]", names(data), value = TRUE)
      } else if(grepl("pre_surgery_3d_to_7d", period_name)) {
        cols <- grep("pre.*3d.*7d", names(data), value = TRUE)
      } else if(grepl("pre_surgery_7d_all", period_name)) {
        cols <- grep("pre.*7d.*all", names(data), value = TRUE)
      } else if(grepl("pre_surgery_all", period_name)) {
        cols <- grep("pre.*all", names(data), value = TRUE)
      } else if(grepl("post_surgery_1w", period_name)) {
        cols <- grep("post.*1w", names(data), value = TRUE)
      } else if(grepl("post_surgery_7d$", period_name)) {
        cols <- grep("post.*7d$", names(data), value = TRUE)
      } else if(grepl("post_surgery_7d_to_30d", period_name)) {
        cols <- grep("post.*7d.*30d", names(data), value = TRUE)
      } else if(grepl("post_surgery_day23_to_30", period_name)) {
        cols <- grep("post.*day23.*30", names(data), value = TRUE)
      } else if(grepl("post_surgery_day27_to_30", period_name)) {
        cols <- grep("post.*day27.*30", names(data), value = TRUE)
      } else if(grepl("post_surgery_over_30d", period_name)) {
        cols <- grep("post.*over_30d", names(data), value = TRUE)
      }
    }
  }
  
  return(unique(cols))
}

prepare_period_dataset <- function(period_name, 
                                   rhr_data = time_period_rhr_results, 
                                   bo_data = bo_time_period_results,
                                   steps_data = steps_time_period_results_assigned,
                                   sleep_data = sleep_time_period_results,
                                   coverage_status,
                                   disease_data,
                                   vision_data,
                                   bloodflow_var_T0,
                                   thickness_var_T0,
                                   bloodflow_var_T1,
                                   thickness_var_T1,
                                   bloodflow_var_T2,
                                   thickness_var_T2) {
  
  # Get participants with adequate coverage for this period
  valid_participants <- filtered_participants_by_period(coverage_status, period_name)
  
  # 提取各个数据源的相关列
  rhr_columns <- extract_period_columns(rhr_data, period_name)
  bo_columns <- extract_period_columns(bo_data, period_name)
  steps_columns <- extract_period_columns(steps_data, period_name)
  sleep_columns <- extract_period_columns(sleep_data, period_name)
  
  # 打印调试信息
  cat("Period:", period_name, "\n")
  cat("  RHR columns found:", length(rhr_columns), "\n")
  cat("  BO columns found:", length(bo_columns), "\n")
  cat("  Steps columns found:", length(steps_columns), "\n")
  cat("  Sleep columns found:", length(sleep_columns), "\n")
  
  # 提取RHR指标数据
  if(length(rhr_columns) > 0) {
    rhr_period_data <- rhr_data %>%
      dplyr::select(subject_id, all_of(rhr_columns))
  } else {
    # 如果没有找到匹配列，创建一个只有subject_id的数据框
    rhr_period_data <- data.frame(subject_id = unique(rhr_data$subject_id))
  }
  
  # 提取BO指标数据
  if(length(bo_columns) > 0) {
    bo_period_data <- bo_data %>%
      dplyr::select(subject_id, all_of(bo_columns))
  } else {
    bo_period_data <- data.frame(subject_id = unique(bo_data$subject_id))
  }
  
  # 提取Steps指标数据
  if(length(steps_columns) > 0) {
    steps_period_data <- steps_data %>%
      dplyr::select(subject_id, all_of(steps_columns))
  } else {
    steps_period_data <- data.frame(subject_id = unique(steps_data$subject_id))
  }
  
  # 提取Sleep指标数据
  if(length(sleep_columns) > 0) {
    sleep_period_data <- sleep_data %>%
      dplyr::select(subject_id, all_of(sleep_columns))
  } else {
    sleep_period_data <- data.frame(subject_id = unique(sleep_data$subject_id))
  }
  
  # 合并所有数据
  combined_data <- rhr_period_data %>%
    left_join(bo_period_data, by = "subject_id") %>%
    left_join(steps_period_data, by = "subject_id") %>%
    left_join(sleep_period_data, by = "subject_id") %>%
    # 与基线疾病数据合并
    left_join(disease_data, by = c("subject_id" = "ID")) %>%
    # 与视力数据合并
    left_join(vision_data, by = c("subject_id" = "ID")) %>%
    # 与OCTA基线(T0)血流数据合并
    left_join(bloodflow_var_T0, by = c("subject_id" = "ID")) %>%
    # 与OCTA基线(T0)厚度数据合并
    left_join(thickness_var_T0, by = c("subject_id" = "ID")) %>%
    # 与OCTA一周(T1)血流数据合并
    left_join(bloodflow_var_T1, by = c("subject_id" = "ID")) %>%
    # 与OCTA一周(T1)厚度数据合并
    left_join(thickness_var_T1, by = c("subject_id" = "ID")) %>%
    # 与OCTA一月(T2)血流数据合并
    left_join(bloodflow_var_T2, by = c("subject_id" = "ID")) %>%
    # 与OCTA一月(T2)厚度数据合并
    left_join(thickness_var_T2, by = c("subject_id" = "ID")) %>%
    # 只包含具有足够覆盖率的参与者
    filter(subject_id %in% valid_participants)
  
  # 添加时间周期信息
  combined_data$time_period <- period_name
  
  # 添加覆盖率相关信息
  if(period_name %in% names(coverage_status)) {
    coverage_data <- coverage_status[[period_name]] %>%
      dplyr::select(subject_id, unique_hours, days_covered, required_threshold)
    
    combined_data <- combined_data %>%
      left_join(coverage_data, by = "subject_id")
  }
  
  return(combined_data)
}

# Step 8: Generate datasets for each time period
all_periods <- names(time_periods)
period_datasets <- list()

for(period in all_periods) {
  period_datasets[[period]] <- prepare_period_dataset(
    period_name = period,
    rhr_data = time_period_rhr_results,
    bo_data = bo_time_period_results,
    steps_data = steps_time_period_results,
    sleep_data = sleep_time_period_results,
    coverage_status = time_period_coverage_status,
    disease_data = disease_data,
    vision_data = vision_data,
    bloodflow_var_T0 = bloodflow_var_T0,
    thickness_var_T0 = thickness_var_T0,
    bloodflow_var_T1 = bloodflow_var_T1,
    thickness_var_T1 = thickness_var_T1,
    bloodflow_var_T2 = bloodflow_var_T2,
    thickness_var_T2 = thickness_var_T2
  )
}

# Step 9: Analyze sample sizes for each period
sample_sizes <- sapply(all_periods, function(period) {
  nrow(period_datasets[[period]])
})

# Create a data frame for easy viewing
period_sample_sizes <- data.frame(
  time_period = all_periods,
  days_in_period = sapply(time_periods, function(x) x[2] - x[1] + 1),
  required_threshold = sapply(time_periods, calculate_threshold),
  sample_size = sample_sizes
)

# Print the sample sizes
print("Sample sizes by time period:")
print(period_sample_sizes)

# Step 10: Function to prepare data for modeling
prepare_for_modeling <- function(dataset, outcome_type = "continuous") {
  # Remove identifier columns and coverage metrics
  modeling_data <- dataset %>%
    dplyr::select(-subject_id, -surgery_eye_1, 
                  -matches("unique_hours"), -matches("days_covered"), -matches("required_threshold"))
  
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
  
  # Handle missing values with median imputation
  predictor_cols <- setdiff(names(modeling_data), 
                            c("vision_improvement", "vision_improved", "vision_improved_factor"))
  
  # Check if we have enough data for imputation
  if(nrow(modeling_data) > 5) {
    preprocess_steps <- preProcess(modeling_data[, predictor_cols], method = c("medianImpute"))
    imputed_predictors <- predict(preprocess_steps, modeling_data[, predictor_cols])
    
    # Combine imputed predictors with outcome
    imputed_data <- cbind(imputed_predictors, modeling_data[, c(outcome_col)])
    names(imputed_data)[ncol(imputed_data)] <- outcome_col
    
    return(imputed_data)
  } else {
    warning("Too few observations for imputation. Returning original data.")
    return(modeling_data)
  }
}

# Step 11: Export each period's dataset to CSV for further analysis
# Create directory if it doesn't exist
output_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Process and export each period's dataset
for(period in all_periods) {
  # Skip empty datasets
  if(length(period_datasets[[period]]) == 0) {
    cat("No data available for period", period, "\n\n")
    next
  }
  
  # Check if it's a data frame
  if(!is.data.frame(period_datasets[[period]])) {
    cat("Warning: Data for period", period, "is not a data frame\n")
    cat("  Class:", class(period_datasets[[period]]), "\n")
    
    # Try to convert list to data frame if possible
    if(is.list(period_datasets[[period]])) {
      tryCatch({
        period_datasets[[period]] <- as.data.frame(period_datasets[[period]])
        cat("  Successfully converted to data frame\n")
      }, error = function(e) {
        cat("  Failed to convert to data frame:", e$message, "\n")
      })
    }
  }
  
  # Check again if it's now a data frame
  if(is.data.frame(period_datasets[[period]])) {
    # Check for list columns that would cause problems with CSV export
    list_cols <- sapply(period_datasets[[period]], is.list)
    if(any(list_cols)) {
      cat("  Warning: Data frame for period", period, "contains list columns\n")
      cat("  List columns:", names(list_cols[list_cols]), "\n")
      
      # Remove list columns for CSV export
      period_datasets[[period]] <- period_datasets[[period]][, !list_cols]
      cat("  Removed list columns for CSV export\n")
    }
    
    # Export to CSV - 使用简化的文件名
    simplified_period <- gsub("surgery_", "", period)
    simplified_period <- gsub("_to_", "_", simplified_period)
    filename <- file.path(output_dir, paste0(simplified_period, "_data.csv"))
    
    tryCatch({
      write.csv(period_datasets[[period]], file = filename, row.names = FALSE)
      cat("Exported data for period", period, "to", filename, "\n")
      cat("  Number of participants:", nrow(period_datasets[[period]]), "\n")
      cat("  Number of variables:", ncol(period_datasets[[period]]), "\n\n")
    }, error = function(e) {
      cat("  Failed to export data for period", period, ":", e$message, "\n\n")
    })
  } else {
    cat("  Unable to export data for period", period, "- not a proper data frame\n\n")
  }
}

# Save coverage summary with additional information
coverage_summary_file <- file.path(output_dir, "time_period_coverage_summary.csv")
write.csv(coverage_summary, file = coverage_summary_file, row.names = FALSE)
cat("Exported coverage summary to", coverage_summary_file, "\n\n")

# Save sample size information
sample_sizes_file <- file.path(output_dir, "period_sample_sizes.csv")
write.csv(period_sample_sizes, file = sample_sizes_file, row.names = FALSE)
cat("Exported sample sizes to", sample_sizes_file, "\n\n")

# Combine datasets for all periods
valid_dfs <- list()
for(period in all_periods) {
  if(length(period_datasets[[period]]) > 0 && is.data.frame(period_datasets[[period]])) {
    valid_dfs[[period]] <- period_datasets[[period]]
  }
}

if(length(valid_dfs) > 0) {
  # Combine data and export
  tryCatch({
    combined_data <- bind_rows(valid_dfs, .id = "period")
    combined_file <- file.path(output_dir, "all_periods_combined_data.csv")
    write.csv(combined_data, file = combined_file, row.names = FALSE)
    cat("Exported combined data for all periods to", combined_file, "\n")
    cat("  Total observations:", nrow(combined_data), "\n")
    cat("  Number of variables:", ncol(combined_data), "\n\n")
  }, error = function(e) {
    cat("Failed to combine datasets:", e$message, "\n\n")
  })
} else {
  cat("No valid data frames available to combine\n\n")
}

# Save R data objects
save(period_datasets, file = file.path(output_dir, "vision_prediction_period_datasets.rda"))
save(period_sample_sizes, file = file.path(output_dir, "vision_prediction_period_sample_sizes.rda"))
save(time_period_coverage_status, file = file.path(output_dir, "time_period_coverage_status.rda"))
save(coverage_summary, file = file.path(output_dir, "time_period_coverage_summary.rda"))

# 创建覆盖率详细报告 - 每个参与者在每个时间窗口的覆盖情况
coverage_details <- list()
for (period_name in names(time_period_coverage_status)) {
  period_data <- time_period_coverage_status[[period_name]] %>%
    mutate(time_period = period_name)
  coverage_details[[period_name]] <- period_data
}

# 合并所有时间窗口的覆盖率数据
all_coverage_details <- bind_rows(coverage_details)

# 导出覆盖率详情
coverage_details_file <- file.path(output_dir, "participant_coverage_details.csv")
write.csv(all_coverage_details, file = coverage_details_file, row.names = FALSE)
cat("Exported participant coverage details to", coverage_details_file, "\n\n")

cat("Data preparation complete. Datasets saved in", output_dir, "\n")

# 创建阈值分析报告 - 分析不同阈值对样本量的影响
threshold_analysis <- list()
for (period_name in names(time_period_coverage_status)) {
  period_data <- time_period_coverage_status[[period_name]]
  period_range <- time_periods[[period_name]]
  days_in_period <- period_range[2] - period_range[1] + 1
  
  # 测试不同阈值
  thresholds_to_test <- c(4, 8, 12, 16, 20, 24, 28, 32)
  threshold_results <- lapply(thresholds_to_test, function(threshold) {
    subjects_meeting <- sum(period_data$unique_hours >= threshold, na.rm = TRUE)
    percentage <- round(subjects_meeting / nrow(period_data) * 100, 1)
    
    data.frame(
      time_period = period_name,
      days_in_period = days_in_period,
      threshold = threshold,
      subjects_meeting = subjects_meeting,
      total_subjects = nrow(period_data),
      percentage = percentage
    )
  }) %>% bind_rows()
  
  threshold_analysis[[period_name]] <- threshold_results
}

# 合并所有阈值分析结果
all_threshold_analysis <- bind_rows(threshold_analysis)

# 导出阈值分析报告
threshold_analysis_file <- file.path(output_dir, "threshold_impact_analysis.csv")
write.csv(all_threshold_analysis, file = threshold_analysis_file, row.names = FALSE)
cat("Exported threshold impact analysis to", threshold_analysis_file, "\n\n")

# Print structure of period datasets for verification
cat("\nStructure verification for period datasets:\n")
for(period in names(period_datasets)) {
  cat("Period", period, "- Class:", class(period_datasets[[period]]), "\n")
  if(is.data.frame(period_datasets[[period]])) {
    cat("  Is data frame: YES\n")
    cat("  Dimensions:", dim(period_datasets[[period]]), "\n")
    
    # Check for list columns
    list_cols <- which(sapply(period_datasets[[period]], is.list))
    if(length(list_cols) > 0) {
      cat("  Has list columns: YES -", names(period_datasets[[period]])[list_cols], "\n")
    } else {
      cat("  Has list columns: NO\n")
    }
  } else {
    cat("  Is data frame: NO\n")
  }
  cat("\n")
}

# 创建覆盖率统计可视化
if(require(ggplot2)) {
  # 1. 时间窗口长度与覆盖率的关系图
  coverage_by_length <- ggplot(coverage_summary, 
                               aes(x = days_in_period, y = coverage_rate)) +
    geom_point(aes(size = total_subjects), alpha = 0.7) +
    geom_text(aes(label = period), vjust = -1, size = 3) +
    geom_smooth(method = "loess", se = TRUE, color = "blue", linetype = "dashed") +
    labs(title = "Coverage Rate by Time Period Length",
         x = "Days in Period", 
         y = "Coverage Rate (%)",
         size = "Number of\nSubjects") +
    theme_minimal()
  
  # 保存图形
  coverage_plot_file <- file.path(output_dir, "coverage_by_period_length.pdf")
  ggsave(coverage_plot_file, plot = coverage_by_length, width = 8, height = 6)
  
  # 2. 阈值对样本量的影响图
  threshold_impact <- ggplot(all_threshold_analysis, 
                             aes(x = threshold, y = percentage, color = factor(days_in_period))) +
    geom_line() +
    geom_point() +
    facet_wrap(~time_period, scales = "free_y") +
    labs(title = "Impact of Different Thresholds on Sample Size",
         x = "Required Unique Hours", 
         y = "Percentage of Subjects Meeting Threshold",
         color = "Days in Period") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 保存图形
  threshold_plot_file <- file.path(output_dir, "threshold_impact_plot.pdf")
  ggsave(threshold_plot_file, plot = threshold_impact, width = 10, height = 8)
}

