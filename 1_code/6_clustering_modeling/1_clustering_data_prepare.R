# 按糖尿病状态和手术类型创建矩阵文件，并筛选每天佩戴时长超过8小时以上的人群
# 该脚本处理每日数据，创建一个矩阵，其中：
# - 行是样本（参与者）
# - 列是时间点（相对于手术的天数）
# - 数据分为两组：
#   1. 无糖尿病 + 手术类型 0/1
#   2. 有糖尿病 + 手术类型 1

library(tidyverse)
library(lubridate)
setwd(get_project_wd())
rm(list = ls())

# Load the data files
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/daily_rhr_result.rda")
load("3_data_analysis/2_data_analysis/bo/bo_time_periods/daily_bo_result.rda")
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/daily_steps_result.rda") # Note: Using daily_steps_result instead of daily_steps_result_assigned
load("3_data_analysis/2_data_analysis/sleep/sleep_time_periods/daily_sleep_result.rda")
# Load heart rate data directly for wear time calculation
load("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda")

# Print dimensions of each dataset to confirm loading
cat("Dimensions of each dataset:\n")
cat("daily_rhr_result:", dim(daily_rhr_result), "\n")
cat("daily_bo_result:", dim(daily_bo_result), "\n")
cat("daily_steps_result:", dim(daily_steps_result), "\n")
cat("daily_sleep_result:", dim(daily_sleep_result), "\n")

# Step 1: Load the baseline_info data to get diabetes status and surgery type
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# Create group information based on diabetes status and surgery type
group_info <- baseline_info %>%
  mutate(
    # Based on the conversion logic in the provided code
    dm_status = case_when(
      diabetes_history == 1 ~ "Diabetes",
      diabetes_history == 2 ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    # Use surgery_1..0.PI.1.other. as surgery type
    surgery_type = surgery_1..0.PI.1.other.
  ) %>%
  dplyr::select(ID, dm_status, surgery_type)

# Print group distribution
cat("\nDiabetes status and surgery type distribution:\n")
print(table(group_info$dm_status, group_info$surgery_type, dnn = c("Diabetes Status", "Surgery Type")))

# Step 2: Analyze daily wear time coverage using heart rate data
# 修改后的覆盖率分析函数 - 仅要求全天总计≥8小时
analyze_time_period_coverage <- function(heart_rate_data, baseline_info) {
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
      day_point >= -4,
      day_point <= 7
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
    coverage_percentage = round(sum(meets_threshold) / n() * 100, 1),
    .groups = "drop"
  )

# 打印覆盖率信息
cat("\n每日佩戴时长≥8小时的参与者数量：\n")
print(daily_coverage_count)

# Step 3: Filter participants based on 8+ hours criteria
# For each day from -7 to 6, get the IDs of participants with adequate coverage
filtered_participants_by_day <- function(daily_coverage_status, target_day) {
  daily_coverage_status %>%
    filter(day_point == target_day, meets_threshold == TRUE) %>%
    pull(subject_id)
}

# Step 4: Merge all datasets by subject_id
# All datasets have subject_id so we can use that for merging
merged_data <- daily_rhr_result %>%
  left_join(daily_bo_result, by = "subject_id") %>%
  left_join(daily_steps_result, by = "subject_id") %>%
  left_join(daily_sleep_result, by = "subject_id")

# Check if we have surgery_date columns in each dataset and rename them if they exist
merged_data <- merged_data %>%
  rename_with(~ case_when(
    . == "surgery_date.x" ~ "surgery_date_rhr",
    . == "surgery_date.y" ~ "surgery_date_bo",
    . == "surgery_date.x.x" ~ "surgery_date_steps",
    . == "surgery_date.y.y" ~ "surgery_date_sleep",
    TRUE ~ .
  ))

# Print dimensions of merged data
cat("\nMerged data dimensions before filtering:", dim(merged_data), "\n")

# Step 5: Filter the merged data to only include participants with 8+ hours of wear time
# Extract the day_point from the merged data
merged_data_with_day <- merged_data %>%
  mutate(
    # Assuming there's a measure_date column in the merged data
    # If not, adjust according to your actual data structure
    day_point = case_when(
      !is.na(day_point.x) ~ day_point.x,
      !is.na(day_point.y) ~ day_point.y, 
      !is.na(day_point.x.x) ~ day_point.x.x,
      !is.na(day_point.y.y) ~ day_point.y.y,
      TRUE ~ NA_real_
    )
  )

# Filter merged data to only include participants with sufficient coverage for each day
filtered_merged_data <- merged_data_with_day %>%
  inner_join(
    daily_coverage_status %>% 
      filter(meets_threshold == TRUE),
    by = c("subject_id", "day_point")
  )

# Print dimensions after filtering
cat("\nMerged data dimensions after filtering for 8+ hour wear time:", dim(filtered_merged_data), "\n")

# Step 6: Merge group information with the filtered health data
# Assuming subject_id in merged_data corresponds to ID in group_info
merged_grouped_data <- filtered_merged_data %>%
  left_join(group_info, by = c("subject_id" = "ID"))

# Check for missing group information
missing_group_info <- sum(is.na(merged_grouped_data$dm_status))
if (missing_group_info > 0) {
  cat("\nWarning:", missing_group_info, "rows are missing group information\n")
}

# Create output directory
output_dir <- "3_data_analysis/6_clustering_modeling/data_prepare"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Split data into the two specified groups
# 1. No Diabetes Surgery 0 group - now including Other/NA cases
no_diabetes_surgery0 <- merged_grouped_data %>%
  filter(surgery_type == 0 | is.na(dm_status) | is.na(surgery_type))

# 2. Diabetes Surgery 1 group
diabetes_surgery1 <- merged_grouped_data %>%
  filter(dm_status == "Diabetes", surgery_type == 1)

# Save the grouped data
if (nrow(no_diabetes_surgery0) > 0) {
  no_diabetes_file <- file.path(output_dir, "combined_NoD_Surg0_8h.csv")
  write.csv(no_diabetes_surgery0, no_diabetes_file, row.names = FALSE)
  cat("\nSaved No Diabetes Surgery 0 group data (8+ hour wear time) to:", no_diabetes_file, "\n")
  cat("Number of rows:", nrow(no_diabetes_surgery0), "\n")
  cat("Number of unique subject_ids:", length(unique(no_diabetes_surgery0$subject_id)), "\n")
} else {
  cat("\nNo Diabetes Surgery 0 group has no data to save\n")
}

if (nrow(diabetes_surgery1) > 0) {
  diabetes_file <- file.path(output_dir, "combined_D_Surg1_8h.csv")
  write.csv(diabetes_surgery1, diabetes_file, row.names = FALSE)
  cat("\nSaved Diabetes Surgery 1 group data (8+ hour wear time) to:", diabetes_file, "\n")
  cat("Number of rows:", nrow(diabetes_surgery1), "\n")
  cat("Number of unique subject_ids:", length(unique(diabetes_surgery1$subject_id)), "\n")
} else {
  cat("\nDiabetes Surgery 1 group has no data to save\n")
}

# Create a summary of the final merged data
cat("\nFinal data summary (8+ hour wear time):\n")
cat("Total number of rows:", nrow(merged_grouped_data), "\n")
cat("Number of rows in No Diabetes Surgery 0 group:", nrow(no_diabetes_surgery0), "\n")
cat("Number of rows in Diabetes Surgery 1 group:", nrow(diabetes_surgery1), "\n")
cat("Number of unique subject_ids in final data:", length(unique(merged_grouped_data$subject_id)), "\n")

# Save the daily coverage status for reference
coverage_file <- file.path(output_dir, "daily_coverage_status_8h.csv")
write.csv(daily_coverage_status, coverage_file, row.names = FALSE)
cat("\nSaved daily coverage status to:", coverage_file, "\n")