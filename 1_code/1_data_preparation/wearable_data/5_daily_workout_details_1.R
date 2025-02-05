library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)


read_dailyworkout_details_file <- function(file_path) {
  data <- read.csv(file_path,
                   stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8")
  
  # Check if the first column contains numeric timestamps or IDs
  first_col <- data[[1]]
  if (!any(grepl("-", first_col))) {
    # If the first column doesn't contain date format (no "-"), skip it
    data <- data[-1]
  }
  
  # Check number of columns and assign column names accordingly
  if(ncol(data) == 7) {
    colnames(data) <- c("数据时间", "步数(steps)", "运动距离(m)", "卡路里", "累计爬升高度", "心率", "外部ID")
    data$'活动名称' <- NA
  } else if(ncol(data) == 8) {
    colnames(data) <- c("数据时间", "活动名称", "步数(steps)", "运动距离(m)", "卡路里", "累计爬升高度", "心率", "外部ID")
  } else {
    stop(paste("Unexpected number of columns:", ncol(data)))
  }
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data_1", recursive = FALSE)
daily_workout_details_data <- purrr::map(subject_folders, function(folder) {
  subject_id <- basename(folder)
  cat(subject_id, " ")
  
  hr_files <- list.files(
    folder,
    pattern = "t_apvddepp_dailyworkoutdetail.*\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(hr_files) > 0) {
    subject_data <- map_df(hr_files, read_dailyworkout_details_file)
    
    # Add subject_id before renaming columns
    subject_data$subject_id <- subject_id
    
    colnames(subject_data) <- c("measure_time", "activity", "steps", "distance", "calorie", "climbing_height", "heartrate", "external_id", "subject_id")
    
    subject_data <- subject_data %>%
      dplyr::select(subject_id, measure_time, activity, steps, distance, calorie, climbing_height, heartrate) %>%
      dplyr::mutate(sample_id = paste(subject_id, measure_time, sep = "_")) %>%
      dplyr::mutate(measure_time = as.POSIXct(measure_time)) %>%
      dplyr::distinct(sample_id, .keep_all = TRUE)
    
    return(subject_data)
  } else {
    return(NULL)
  }
})


daily_workout_details_data <- daily_workout_details_data %>%
  do.call(rbind, .) %>%
  as.data.frame()

sample_info <- daily_workout_details_data %>%
  dplyr::select(sample_id, subject_id, measure_time, activity)

sample_info_fixed <- sample_info %>%
  mutate(
    subject_id = toupper(subject_id),  # 统一转换为大写
    subject_id = gsub("SHH|ShH", "SH", subject_id),  # 修复ShH061的问题
    measure_time = substr(measure_time, 1, 19),  # 保留到秒
    sample_id = paste(subject_id, measure_time, sep = "_")  
  )

expression_data <- daily_workout_details_data %>%
  dplyr::select(steps, distance, calorie, climbing_height, heartrate) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <- sample_info_fixed$sample_id

expression_data <- expression_data[, sample_info_fixed$sample_id]

variable_info <- data.frame(
  variable_id = c(
    "steps",
    "distance",
    "calorie",
    "climbing_height",
    "heartrate"
  )
)

sample_info_fixed$class <- "Subject"

library(massdataset)

daily_workout_details_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info_fixed,
    variable_info = variable_info
  )

dir.create("3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details", recursive = TRUE)
save(daily_workout_details_data, file = "3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details/daily_workout_details_data.rda", compress = "xz")
