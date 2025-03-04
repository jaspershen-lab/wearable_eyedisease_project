library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)
library(caret)
library(randomForest)

# Load the RHR data
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/daily_rhr_result.rda")
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/time_period_rhr_results.rda")
load("3_data_analysis/2_data_analysis/bo/bo_time_periods/daily_bo_result.rda")
load("3_data_analysis/2_data_analysis/sleep/sleep_time_periods/daily_sleep_result.rda")
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/daily_steps_result.rda")

#############
dir.create("3_data_analysis/2_data_analysis/combined_daily_wearbale_data/combined_daily_wearbale_data", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/combined_daily_wearbale_data/combined_daily_wearbale_data")

# 合并所有wearable数据
merge_wearable_data <- function(rhr_data, bo_data, sleep_data, steps_data) {
  # 获取所有时间点
  days <- seq(-10, 90)
  
  # 获取所有subjects
  all_subjects <- unique(c(
    rhr_data$subject_id,
    bo_data$subject_id,
    sleep_data$subject_id,
    steps_data$subject_id
  ))
  
  # 处理每个数据集
  process_dataset <- function(data, prefix) {
    # 获取该数据集的所有测量列
    measure_cols <- grep(paste0("day_.*_", prefix), names(data), value = TRUE)
    
    # 选择需要的列并重命名
    data_processed <- data %>%
      dplyr::select(subject_id, all_of(measure_cols))
    
    return(data_processed)
  }
  
  # 处理每个数据集
  rhr_processed <- process_dataset(rhr_data, "rhr")
  bo_processed <- process_dataset(bo_data, "bo")
  sleep_processed <- process_dataset(sleep_data, "sleep")
  steps_processed <- process_dataset(steps_data, "steps")
  
  # 合并所有数据
  merged_data <- all_subjects %>%
    as.data.frame() %>%
    setNames("subject_id") %>%
    left_join(rhr_processed, by = "subject_id") %>%
    left_join(bo_processed, by = "subject_id") %>%
    left_join(sleep_processed, by = "subject_id") %>%
    left_join(steps_processed, by = "subject_id")
  
  return(merged_data)
}

# 合并数据
merged_data <- merge_wearable_data(
  daily_rhr_result,
  daily_bo_result,
  daily_sleep_result,
  daily_steps_result
)

# 创建sample_info
sample_info <- data.frame(
  sample_id = names(merged_data)[-1],  # 除去subject_id列
  class = "Subject"
) 

# 创建expression_data
expression_data <- merged_data[, -1]  # 移除 subject_id 列但保持 data.frame 格式
rownames(expression_data) <- merged_data$subject_id

# 创建variable_info
variable_info <- data.frame(
  variable_id = merged_data$subject_id,
  description = "Patient ID"
)

# 创建mass_dataset
combined_daily_wearbale_data <- create_mass_dataset(
  expression_data = expression_data,
  sample_info = sample_info,
  variable_info = variable_info
)

# 保存结果
save(combined_daily_wearbale_data, file = "combined_daily_wearbale_data.rda", compress = "xz")