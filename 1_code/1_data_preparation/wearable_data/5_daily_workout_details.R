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
  
  if(ncol(data) == 7) {
    colnames(data) <- c("数据时间", "步数(steps)", "运动距离(m)", "卡路里", "累计爬升高度", "心率", "外部ID")
    data$'活动名称' <- NA
  } else {
    colnames(data) <- c("数据时间", "活动名称", "步数(steps)", "运动距离(m)", "卡路里", "累计爬升高度", "心率", "外部ID")
  }
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data", recursive = FALSE)

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
    
    # 检查并添加缺失的activity列
    if(!"activity" %in% colnames(subject_data)) {
      subject_data$activity <- NA
    }
    
    colnames(subject_data) <- c("data_time", "activity", "steps", "distance", "calorie", "climbing_height", "heartrate", "subject_id")
    
    subject_data <- subject_data %>%
      dplyr::select(subject_id, data_time, activity, steps, distance, calorie, climbing_height, heartrate) %>%
      dplyr::mutate(sample_id = paste(subject_id, data_time, sep = "_")) %>%
      dplyr::mutate(data_time = as.POSIXct(data_time)) %>%
      dplyr::distinct(sample_id, .keep_all = TRUE)
    
    return(subject_data)
  } else {
    return(NULL)
  }
})

daily_workout_details_data[[1]] %>%
  head(1000) %>%
  ggplot(aes(data_time, steps)) +
  geom_line()


daily_workout_details_data <-
  daily_workout_details_data %>%
  do.call(rbind, .) %>%
  as.data.frame()


sample_info <-
  daily_workout_details_data %>%
  dplyr::select(sample_id, subject_id, data_time, activity)

expression_data <-
  daily_workout_details_data %>%
  dplyr::select(steps,distance,calorie,climbing_height,heartrate) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info$sample_id

variable_info <- data.frame(
  variable_id = c(
    "steps",
    "distance",
    "calorie",
    "climbing_height",
    "heartrate"
  )
)

sample_info$class <- "Subject"

library(massdataset)

daily_workout_details_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

dir.create("3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details", recursive = TRUE)
save(daily_workout_details_data, file = "3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details/daily_workout_details_data.rda", compress = "xz")
