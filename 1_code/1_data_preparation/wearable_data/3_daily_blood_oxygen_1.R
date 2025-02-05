library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)


read_dailybloodoxygen_file <- function(file_path) {
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
  if (ncol(data) == 6) {
    data <- data[, -1]
    colnames(data) <- c("数据时间", "最大血氧饱和度(%)", "最小血氧饱和度(%)", "平均血氧饱和度(%)", "外部ID")
  } else if (ncol(data) == 5) {
    colnames(data) <- c("数据时间", "最大血氧饱和度(%)", "最小血氧饱和度(%)", "平均血氧饱和度(%)", "外部ID")
  } else {
    stop(paste("Unexpected number of columns:", ncol(data)))
  }
  
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data_1", recursive = FALSE)
daily_blood_oxygen_data <-
  purrr::map(subject_folders, function(folder) {
    subject_id <- basename(folder)
    cat(subject_id, " ")
    
    hr_files <- list.files(
      folder,
      pattern = "t_apvddepp_dailybloodoxygensaturation.*\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(hr_files) > 0) {
      subject_data <- map_df(hr_files, read_dailybloodoxygen_file)
      
      # Adjust column names based on number of columns
      if (ncol(subject_data) == 5) {
        colnames(subject_data) <- c("data_time", "daily_blood_oxygen_max", "daily_blood_oxygen_min", "daily_blood_oxygen_mean", "subject_id")
        subject_data <- subject_data %>%
          dplyr::select(subject_id, data_time, daily_blood_oxygen_max, daily_blood_oxygen_min, daily_blood_oxygen_mean)
      } else {
        colnames(subject_data) <- c("data_time", "daily_blood_oxygen_max", "daily_blood_oxygen_min", "daily_blood_oxygen_mean")
        subject_data <- subject_data %>%
          dplyr::mutate(subject_id = subject_id)  # Add subject_id column
      }
      
      subject_data <- subject_data %>%
        dplyr::mutate(sample_id = paste(subject_id, data_time, sep = "_")) %>%
        dplyr::mutate(data_time = as.POSIXct(data_time, format = "%Y-%m-%d %H:%M:%S")) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      return(subject_data)
    } else {
      return(NULL)
    }
  })

daily_blood_oxygen_data[[2]] %>%
  head(1000) %>%
  ggplot(aes(data_time, daily_blood_oxygen_mean)) +
  geom_line()


daily_blood_oxygen_data <-
  daily_blood_oxygen_data %>%
  do.call(rbind, .) %>%
  as.data.frame()


sample_info <-
  daily_blood_oxygen_data %>%
  dplyr::select(sample_id, subject_id, data_time)

sample_info_fixed <- sample_info %>%
  mutate(
    subject_id = toupper(subject_id),  # 统一转换为大写
    subject_id = gsub("SHH|ShH", "SH", subject_id),  # 修复ShH061的问题
    data_time = paste0(substr(data_time, 1, 10), " 00:00:00"),  # 将 "+08" 替换为 "00:00:00"
    sample_id = paste(subject_id, data_time, sep = "_")  # 重新构建 sample_id
  )

expression_data <-
  daily_blood_oxygen_data %>%
  dplyr::select(daily_blood_oxygen_max,daily_blood_oxygen_min,daily_blood_oxygen_mean) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info_fixed$sample_id

expression_data <- expression_data[, sample_info_fixed$sample_id]

variable_info <- data.frame(
  variable_id = c(
    "daily_blood_oxygen_max",
    "daily_blood_oxygen_min", 
    "daily_blood_oxygen_mean"
  )
)

sample_info_fixed$class <- "Subject"

library(massdataset)

daily_blood_oxygen_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info_fixed,
    variable_info = variable_info
  )

dir.create("3_data_analysis/1_data_preparation/wearable_data/3_daily_blood_oxygen", recursive = TRUE)
save(daily_blood_oxygen_data, file = "3_data_analysis/1_data_preparation/wearable_data/3_daily_blood_oxygen/daily_blood_oxygen_data.rda", compress = "xz")
