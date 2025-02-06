library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)


read_dailyworkout_file <- function(file_path) {
  data <- read.csv(file_path,
                   stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8",
                   colClasses = "character")  
  
  # Check if the first column contains numeric timestamps or IDs
  first_col <- data[[1]]
  if (!any(grepl("-", first_col))) {
    # If the first column doesn't contain date format (no "-"), skip it
    data <- data[-1]
  }
  
  # 只保留需要的列，不管是5列还是6列的数据，都只取需要的数值列
  if(ncol(data) == 5) {
    data <- data[, c(1, 2, 3, 4, 5)]  # 时间，步数，距离，卡路里，外部ID
    colnames(data) <- c("数据时间", "步数(steps)", "运动距离(m)", "卡路里", "外部ID")
  } else {
    data <- data[, c(1, 3, 4, 5, 6)]  # 跳过activity列，只取时间，步数，距离，卡路里，外部ID
    colnames(data) <- c("数据时间", "步数(steps)", "运动距离(m)", "卡路里", "外部ID")
  }
  
  # 转换数据类型
  data$'外部ID' <- as.character(data$'外部ID')
  data$'步数(steps)' <- as.numeric(data$'步数(steps)')
  data$'运动距离(m)' <- as.numeric(data$'运动距离(m)')
  data$'卡路里' <- as.numeric(data$'卡路里')
  
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data_1", recursive = FALSE)
daily_workout_data <- purrr::map(subject_folders, function(folder) {
  cat(basename(folder), " ")
  
  hr_files <- list.files(
    folder,
    pattern = "t_apvddepp_dailyworkout.*\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(hr_files) > 0) {
    subject_data <- map_df(hr_files, read_dailyworkout_file)
    
    # Add subject_id before renaming columns
    subject_data$subject_id <- basename(folder)
    
    colnames(subject_data) <- c("measure_time", "steps", "distance", "calorie", "external_id", "subject_id")
    
    subject_data <- subject_data %>%
      dplyr::select(subject_id, measure_time, steps, distance, calorie) %>%
      dplyr::mutate(
        subject_id = as.character(subject_id)
      ) %>%
      dplyr::mutate(sample_id = paste(subject_id, measure_time, sep = "_")) %>%
      dplyr::mutate(measure_time = as.POSIXct(measure_time, format = "%Y-%m-%d %H:%M:%S")) %>%
      # 对重复数据的处理
      group_by(sample_id) %>%
      summarise(
        subject_id = unique(subject_id)[1],  # 使用unique()[1]替代first()
        measure_time = min(measure_time),          # 使用min()替代first()
        # 如果有非零值，取非零值中的最大值；如果都是零，就取0
        steps = ifelse(all(steps == 0), 0, max(steps[steps > 0], na.rm = TRUE)),
        distance = ifelse(all(distance == 0), 0, max(distance[distance > 0], na.rm = TRUE)),
        calorie = ifelse(all(calorie == 0), 0, max(calorie[calorie > 0], na.rm = TRUE)),
        .groups = "drop"
      )
    
    return(subject_data)
  } else {
    return(NULL)
  }
}) %>%
  bind_rows() %>%
  distinct(sample_id, .keep_all = TRUE)


daily_workout_data %>%
  head(1000) %>%
  ggplot(aes(measure_time, steps)) +
  geom_line()

# mass_dataset
sample_info <- daily_workout_data %>%
  dplyr::select(sample_id, subject_id, measure_time)

sample_info_fixed <- sample_info %>%
  mutate(
    subject_id = toupper(subject_id),  # 统一转换为大写
    subject_id = gsub("SHH|ShH", "SH", subject_id),  # 修复ShH061的问题
    measure_time = paste0(substr(measure_time, 1, 10), " 00:00:00"),  # 将 "+08" 替换为 "00:00:00"
    sample_id = paste(subject_id, measure_time, sep = "_")  # 重新构建 sample_id
  )

sample_info_fixed <- sample_info_fixed %>%
  distinct(sample_id, .keep_all = TRUE)

expression_data <- daily_workout_data %>%
  dplyr::select(steps, distance, calorie) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <- sample_info_fixed$sample_id

expression_data <- expression_data[, sample_info_fixed$sample_id]

variable_info <- data.frame(
  variable_id = c(
    "steps",
    "distance",
    "calorie"
  )
)

sample_info_fixed$class <- "Subject"

daily_workout_data <- create_mass_dataset(
  expression_data = expression_data,
  sample_info = sample_info_fixed,
  variable_info = variable_info
)

daily_workout_data@sample_info$measure_time <-
  as.POSIXct(daily_workout_data@sample_info$measure_time, tz = "Asia/Shanghai")

daily_workout_data <-
  daily_workout_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(subject_id, measure_time)


dir.create("3_data_analysis/1_data_preparation/wearable_data/6_daily_workout", recursive = TRUE)
save(daily_workout_data, file = "3_data_analysis/1_data_preparation/wearable_data/6_daily_workout/daily_workout_data.rda", compress = "xz")
