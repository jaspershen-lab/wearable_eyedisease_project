library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)


read_dailyheartrate_file <- function(file_path) {
  data <- read.csv(file_path,
                   stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8")
  
  # Check if the first column contains numeric timestamps or IDs
  first_col <- data[[1]]
  if (!any(grepl("-", first_col))) {
    # If the first column doesn't contain date format (no "-"), skip it
    data <- data[-1]
  }
  
  # Check number of columns and assign names accordingly
  if (ncol(data) == 4) {
    colnames(data) <- c("数据时间", "最大心率(beats/min)", "最小心率(beats/min)", "外部ID")
  } else if (ncol(data) == 3) {
    colnames(data) <- c("数据时间", "最大心率(beats/min)", "最小心率(beats/min)")
  } else {
    stop(paste("Unexpected number of columns:", ncol(data)))
  }
  
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data_1", recursive = FALSE)
daily_heart_rate_data <-
  purrr::map(subject_folders, function(folder) {
    subject_id <- basename(folder)
    cat(subject_id, " ")
    
    hr_files <- list.files(
      folder,
      pattern = "t_apvddepp_dailyheart.*\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(hr_files) > 0) {
      subject_data <- map_df(hr_files, read_dailyheartrate_file)
      
      # Adjust column names based on number of columns
      if (ncol(subject_data) == 4) {
        colnames(subject_data) <- c("measure_time", "daily_heart_rate_max", "daily_heart_rate_min", "subject_id")
        subject_data <- subject_data %>%
          dplyr::select(subject_id, measure_time, daily_heart_rate_max, daily_heart_rate_min)
      } else {
        colnames(subject_data) <- c("measure_time", "daily_heart_rate_max", "daily_heart_rate_min")
        subject_data <- subject_data %>%
          dplyr::mutate(subject_id = subject_id)  # Add subject_id column
      }
      
      subject_data <- subject_data %>%
        dplyr::mutate(sample_id = paste(subject_id, measure_time, sep = "_")) %>%
        dplyr::mutate(measure_time = as.POSIXct(measure_time, format = "%Y-%m-%d %H:%M:%S")) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      return(subject_data)
    } else {
      return(NULL)
    }
  })


daily_heart_rate_data[[2]] %>%
  head(1000) %>%
  ggplot(aes(measure_time, daily_heart_rate_max)) +
  geom_line()


daily_heart_rate_data <-
  daily_heart_rate_data %>%
  do.call(rbind, .) %>%
  as.data.frame()


sample_info <-
  daily_heart_rate_data %>%
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

expression_data <-
  daily_heart_rate_data %>%
  dplyr::select(daily_heart_rate_max,daily_heart_rate_min) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info_fixed$sample_id

expression_data <- expression_data[, sample_info_fixed$sample_id]

variable_info <- data.frame(
  variable_id = c(
    "daily_heart_rate_max",
    "daily_heart_rate_min"
  )
)

sample_info_fixed$class <- "Subject"

library(massdataset)

daily_heart_rate_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info_fixed,
    variable_info = variable_info
  )


daily_heart_rate_data@sample_info$measure_time <-
  as.POSIXct(daily_heart_rate_data@sample_info$measure_time, tz = "Asia/Shanghai")

daily_heart_rate_data <-
  daily_heart_rate_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(subject_id, measure_time)

dir.create("3_data_analysis/1_data_preparation/wearable_data/4_daily_heart_rate", recursive = TRUE)
save(daily_heart_rate_data, file = "3_data_analysis/1_data_preparation/wearable_data/4_daily_heart_rate/daily_heart_rate_data.rda", compress = "xz")
