library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)
library(lubridate)
library(zoo)
library(rpart)
library(randomForest)
library(caret)
library(tidyverse)
library(lubridate)

###read data

load(
  "3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda"
)
load(
  "3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details/daily_workout_details_data.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")


# 修复心率数据的时间戳
heart_rate_data@sample_info <- heart_rate_data@sample_info %>%
  mutate(
    measure_time = as.POSIXct(
      str_extract(sample_id, "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )

# 修复步数数据的时间戳
daily_workout_details_data@sample_info <- daily_workout_details_data@sample_info %>%
  mutate(
    measure_time = as.POSIXct(
      str_extract(sample_id, "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )

######
dir.create("3_data_analysis/2_data_analysis/RHR/", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/RHR/")

subject_id <-
  heart_rate_data@sample_info %>%
  dplyr::pull(subject_id) %>%
  unique()


sample_info_new <-
  vector("list", length(subject_id))

names(sample_info_new) <- subject_id

for (temp_subject_id in subject_id) {
  cat(temp_subject_id, " ")
  
  temp_hr_sample_info <-
    heart_rate_data@sample_info %>%
    dplyr::filter(subject_id == temp_subject_id)
  
  temp_daily_workout_details_data <-
    daily_workout_details_data %>%
    activate_mass_dataset(what = "sample_info") %>%
    dplyr::filter(subject_id == temp_subject_id)
  
  label <-
    temp_hr_sample_info$measure_time %>%
    purrr::map(function(x) {
      cat(as.character(x), " ")
      time_diff <-
        difftime(x,
                 temp_daily_workout_details_data@sample_info$measure_time,
                 units = "mins") %>%
        as.numeric()
      temp_index <-
        which(time_diff <= 10 & time_diff >= 0)
      if (length(temp_index) == 0) {
        return("<1")
      }
      
      if (length(temp_index) > 0) {
        total_steps <-
          temp_daily_workout_details_data@expression_data[1, temp_index] %>%
          as.numeric() %>%
          sum()
        
        if (total_steps <= 1) {
          return("<1")
        }
        
        if (total_steps > 1 & total_steps <= 50) {
          return("<50")
        }
        
        if (total_steps > 50) {
          return(">50")
        }
      }
    }) %>%
    unlist()
  
  
  temp_hr_sample_info$label <- label
  
  sample_info_new[[temp_subject_id]] <- temp_hr_sample_info
  
}

sample_info_new <-
  sample_info_new %>%
  do.call(rbind, .) %>%
  as.data.frame()


heart_rate_data@sample_info <-
  heart_rate_data@sample_info %>%
  dplyr::left_join(sample_info_new[, c("sample_id", "label")], by = "sample_id")


save(heart_rate_data, file = "heart_rate_data_RHR.rda")


###

# 检查label列是否存在
if("label" %in% colnames(heart_rate_data@sample_info)) {
  # 计算各种标签的频率
  label_table <- table(heart_rate_data@sample_info$label)
  print(label_table)
  
  # 按受试者查看标签分布
  label_by_subject <- table(heart_rate_data@sample_info$subject_id, 
                            heart_rate_data@sample_info$label)
  print(label_by_subject)
} else {
  cat("数据中未找到'label'列，处理可能不完整\n")
}
