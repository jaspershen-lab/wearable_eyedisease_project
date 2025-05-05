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
    # Filter to our desired range (-4 to 7 days)
    filter(
      day_point >= -4,
      day_point <= 30
    ) %>%
    # Count distinct hours per subject-day
    group_by(subject_id, day_point) %>%
    summarise(
      hours_covered = n_distinct(hour),
      .groups = "drop"
    )
  
  # Create complete grid with all subject-day combinations
  all_subjects <- unique(hr_df$subject_id)
  all_days <- seq(-4, 30)
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

# 针对mfuzz聚类，保持宽格式数据结构，同时应用8小时佩戴时间筛选

# 查看合并数据的结构
cat("\n检查合并数据的列名以了解数据结构:\n")
cat("前15个列名:", paste(colnames(merged_data)[1:15], collapse=", "), "\n")

# 从daily_coverage_status中获取符合条件的subject_id和day_point组合
valid_combinations <- daily_coverage_status %>%
  filter(meets_threshold == TRUE) %>%
  select(subject_id, day_point)

# 打印筛选前的统计信息
cat("\n筛选前数据统计:\n")
cat("参与者数量:", length(unique(merged_data$subject_id)), "\n")

# 需要识别哪些列对应哪些天
day_columns <- list()
for(day in seq(-4, 30)) {
  day_prefix <- paste0("day_", day, "_")
  day_columns[[as.character(day)]] <- grep(day_prefix, colnames(merged_data), value = TRUE)
}

# 确定哪些参与者在哪些天符合8小时佩戴时间要求
valid_subject_days <- valid_combinations %>%
  group_by(subject_id) %>%
  summarise(
    valid_days = list(day_point),
    valid_days_count = n(),
    .groups = "drop"
  )

# 打印每个参与者符合条件的天数分布
cat("\n参与者符合8小时佩戴时间要求的天数分布:\n")
valid_days_distribution <- valid_subject_days %>%
  count(valid_days_count) %>%
  mutate(percentage = n / sum(n) * 100)
print(valid_days_distribution)

# 设置参与者最低要求的有效天数（例如总天数的70%）
min_required_days <- round(length(seq(-4, 30)) * 0.7)
cat("\n参与者至少需要", min_required_days, "天符合8小时佩戴时间要求 (70%的总天数)\n")

# 筛选出符合条件的参与者（至少有min_required_days天的数据）
qualified_subjects <- valid_subject_days %>%
  filter(valid_days_count >= min_required_days) %>%
  pull(subject_id)

cat("符合条件的参与者数量:", length(qualified_subjects), "\n")

# 筛选宽格式数据，只保留符合条件的参与者
filtered_wide_data <- merged_data %>%
  filter(subject_id %in% qualified_subjects)

# 为符合条件的参与者，将不符合8小时佩戴时间要求的天的数据设为NA
# 首先创建一个查找表，用于快速检查哪些(subject_id, day_point)组合是有效的
valid_lookup <- valid_combinations %>%
  mutate(is_valid = TRUE) %>%
  unite("subject_day", c(subject_id, day_point), sep = "_", remove = FALSE)

# 用于检查(subject_id, day_point)组合是否有效的函数
is_valid_combination <- function(subject_id, day) {
  lookup_key <- paste(subject_id, day, sep = "_")
  return(lookup_key %in% valid_lookup$subject_day)
}

# 循环处理每个时间点，将不符合条件的数据设为NA
for(day in seq(-4, 30)) {
  day_cols <- day_columns[[as.character(day)]]
  if(length(day_cols) > 0) {
    for(i in 1:nrow(filtered_wide_data)) {
      subject <- filtered_wide_data$subject_id[i]
      if(!is_valid_combination(subject, day)) {
        filtered_wide_data[i, day_cols] <- NA
      }
    }
  }
}

# 查看筛选后每个时间点的数据完整度
day_completeness <- data.frame(day_point = seq(-4, 30))
for(day in seq(-4, 30)) {
  day_cols <- day_columns[[as.character(day)]]
  if(length(day_cols) > 0) {
    # 选取第一列作为代表，检查该时间点有多少非NA值
    sample_col <- day_cols[1]
    valid_count <- sum(!is.na(filtered_wide_data[[sample_col]]))
    day_completeness$valid_count[day_completeness$day_point == day] <- valid_count
    day_completeness$completeness[day_completeness$day_point == day] <- 
      valid_count / nrow(filtered_wide_data) * 100
  }
}

cat("\n筛选后各时间点的数据完整度:\n")
print(day_completeness)

# 合并组信息
filtered_wide_with_groups <- filtered_wide_data %>%
  left_join(group_info, by = c("subject_id" = "ID"))

# 检查是否有缺失的组信息
missing_group_info <- sum(is.na(filtered_wide_with_groups$dm_status))
if(missing_group_info > 0) {
  cat("\n警告:", missing_group_info, "个参与者缺少组信息\n")
}

# 为两个组分别创建mfuzz格式的数据
no_diabetes_surgery0 <- filtered_wide_with_groups %>%
  filter(surgery_type == 0 | is.na(dm_status) | is.na(surgery_type))

diabetes_surgery1 <- filtered_wide_with_groups %>%
  filter(dm_status == "Diabetes", surgery_type == 1)

# 创建输出目录
output_dir <- "3_data_analysis/6_clustering_modeling/data_prepare/1m"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 保存处理后的数据
no_diabetes_file <- file.path(output_dir, "mfuzz_NoD_Surg0_8h.csv")
write.csv(no_diabetes_surgery0, no_diabetes_file, row.names = FALSE)
cat("\n保存无糖尿病+手术0组数据到:", no_diabetes_file, "\n")
cat("行数:", nrow(no_diabetes_surgery0), "\n")
cat("唯一参与者数:", length(unique(no_diabetes_surgery0$subject_id)), "\n")

diabetes_file <- file.path(output_dir, "mfuzz_D_Surg1_8h.csv")
write.csv(diabetes_surgery1, diabetes_file, row.names = FALSE)
cat("\n保存糖尿病+手术1组数据到:", diabetes_file, "\n")
cat("行数:", nrow(diabetes_surgery1), "\n")
cat("唯一参与者数:", length(unique(diabetes_surgery1$subject_id)), "\n")

# 保存参与者和时间点的完整度信息
write.csv(
  valid_subject_days %>% select(-valid_days), 
  file.path(output_dir, "participant_valid_days.csv"), 
  row.names = FALSE
)
write.csv(day_completeness, file.path(output_dir, "day_completeness.csv"), row.names = FALSE)


# 创建mfuzz专用格式的数据文件
# 首先对无糖尿病组
if(nrow(no_diabetes_surgery0) > 0) {
  # 选择需要的列，移除组信息和其他非数据列
  mfuzz_cols <- c()
  for(day in seq(-4, 30)) {
    mfuzz_cols <- c(mfuzz_cols, day_columns[[as.character(day)]])
  }
  
  # 确保所有列都存在
  mfuzz_cols <- mfuzz_cols[mfuzz_cols %in% colnames(no_diabetes_surgery0)]
  
  # 创建mfuzz格式数据
  mfuzz_data_NoD <- no_diabetes_surgery0 %>%
    select(subject_id, all_of(mfuzz_cols))
  
  # 保存mfuzz格式数据
  mfuzz_file_NoD <- file.path(output_dir, "mfuzz_format_NoD_Surg0.csv")
  write.csv(mfuzz_data_NoD, mfuzz_file_NoD, row.names = FALSE)
  cat("\n保存mfuzz格式无糖尿病组数据到:", mfuzz_file_NoD, "\n")
}

# 对糖尿病组
if(nrow(diabetes_surgery1) > 0) {
  # 选择需要的列，移除组信息和其他非数据列
  mfuzz_cols <- c()
  for(day in seq(-4, 30)) {
    mfuzz_cols <- c(mfuzz_cols, day_columns[[as.character(day)]])
  }
  
  # 确保所有列都存在
  mfuzz_cols <- mfuzz_cols[mfuzz_cols %in% colnames(diabetes_surgery1)]
  
  # 创建mfuzz格式数据
  mfuzz_data_D <- diabetes_surgery1 %>%
    select(subject_id, all_of(mfuzz_cols))
  
  # 保存mfuzz格式数据
  mfuzz_file_D <- file.path(output_dir, "mfuzz_format_D_Surg1.csv")
  write.csv(mfuzz_data_D, mfuzz_file_D, row.names = FALSE)
  cat("\n保存mfuzz格式糖尿病组数据到:", mfuzz_file_D, "\n")
}


# -----------------------------------------------------------------------------------
# 添加针对时间序列聚类的数据完整度评估部分（适用于宽格式数据）
# -----------------------------------------------------------------------------------

cat("\n开始进行时间序列聚类的数据完整度评估...\n")

# 基于我们已有的有效日期和参与者组合信息，进行宽格式数据的完整度评估
# 首先分析参与者的完整度（每个参与者有多少天满足佩戴时间要求）
participant_completeness <- valid_subject_days %>%
  select(subject_id, valid_days_count) %>%
  mutate(total_days = length(seq(-4, 30)),
         completeness_rate = valid_days_count / total_days)

# 打印参与者完整度摘要统计
cat("\n参与者数据完整度摘要:\n")
summary_stats_participants <- participant_completeness %>%
  summarise(
    min_completeness = min(completeness_rate),
    q1_completeness = quantile(completeness_rate, 0.25),
    median_completeness = median(completeness_rate),
    mean_completeness = mean(completeness_rate),
    q3_completeness = quantile(completeness_rate, 0.75),
    max_completeness = max(completeness_rate),
    participants_above_50pct = sum(completeness_rate > 0.5),
    participants_above_70pct = sum(completeness_rate > 0.7),
    participants_above_90pct = sum(completeness_rate > 0.9),
    total_participants = n()
  )
print(summary_stats_participants)

# 设定参与者完整度阈值
min_participant_completeness <- 0.7  # 参与者至少70%天数有数据

# 筛选符合条件的参与者
qualified_participants <- participant_completeness %>%
  filter(completeness_rate >= min_participant_completeness) %>%
  pull(subject_id)

cat("\n筛选结果:\n")
cat("符合条件的参与者数量 (完整度>=", min_participant_completeness, "):", length(qualified_participants), "\n")

# 检查是否需要调整阈值
if (length(qualified_participants) < 30) {
  cat("\n警告: 符合条件的参与者数量较少 (<30)。考虑降低参与者完整度阈值。\n")
  
  # 尝试不同的阈值
  thresholds <- c(0.6, 0.5, 0.4)
  for (t in thresholds) {
    alternative_participants <- participant_completeness %>%
      filter(completeness_rate >= t) %>%
      pull(subject_id)
    
    cat("- 如果将参与者完整度阈值设为", t, "，符合条件的参与者数量为:", 
        length(alternative_participants), "\n")
  }
}

# 筛选后的数据 - 只包含完整度足够高的参与者
no_diabetes_filtered <- no_diabetes_surgery0 %>%
  filter(subject_id %in% qualified_participants)

diabetes_filtered <- diabetes_surgery1 %>%
  filter(subject_id %in% qualified_participants)

# 打印筛选后的数据摘要
cat("\n筛选后的数据摘要:\n")
cat("无糖尿病组 - 筛选前行数:", nrow(no_diabetes_surgery0), 
    "筛选后行数:", nrow(no_diabetes_filtered), 
    "减少比例:", round((1 - nrow(no_diabetes_filtered)/nrow(no_diabetes_surgery0)) * 100, 1), "%\n")

cat("糖尿病组 - 筛选前行数:", nrow(diabetes_surgery1), 
    "筛选后行数:", nrow(diabetes_filtered), 
    "减少比例:", round((1 - nrow(diabetes_filtered)/nrow(diabetes_surgery1)) * 100, 1), "%\n")

# 合并筛选后的两组数据
final_filtered_data <- bind_rows(no_diabetes_filtered, diabetes_filtered)

# 保存筛选后的数据
filtered_no_diabetes_file <- file.path(output_dir, "mfuzz_NoD_Surg0_8h_filtered.csv")
write.csv(no_diabetes_filtered, filtered_no_diabetes_file, row.names = FALSE)
cat("\n保存筛选后的无糖尿病组数据到:", filtered_no_diabetes_file, "\n")

filtered_diabetes_file <- file.path(output_dir, "mfuzz_D_Surg1_8h_filtered.csv")
write.csv(diabetes_filtered, filtered_diabetes_file, row.names = FALSE)
cat("\n保存筛选后的糖尿病组数据到:", filtered_diabetes_file, "\n")

combined_filtered_file <- file.path(output_dir, "mfuzz_all_8h_filtered.csv")
write.csv(final_filtered_data, combined_filtered_file, row.names = FALSE)
cat("\n保存筛选后的合并数据到:", combined_filtered_file, "\n")

# 保存完整度信息供参考
write.csv(participant_completeness, file.path(output_dir, "participant_completeness.csv"), row.names = FALSE)
cat("\n保存参与者完整度信息到:", file.path(output_dir, "participant_completeness.csv"), "\n")
