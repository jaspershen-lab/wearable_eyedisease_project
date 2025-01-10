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
  colnames(data) <- c("数据时间", "最大心率(beats/min)", "最小心率(beats/min)","外部ID")
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data", recursive = FALSE)

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
      colnames(subject_data) <-
        c("data_time", "daily_heart_rate_max", "daily_heart_rate_min","subject_id")
      subject_data <-
        subject_data %>%
        dplyr::select(subject_id, data_time, daily_heart_rate_max,daily_heart_rate_min)
      
      subject_data <-
        subject_data %>%
        dplyr::mutate(sample_id = paste(subject_id,data_time, sep = "_")) %>%
        dplyr::mutate(data_time = as.POSIXct(data_time)) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      return(subject_data)
    } else{
      return(NULL)
    }
  })


daily_heart_rate_data[[2]] %>%
  head(1000) %>%
  ggplot(aes(data_time, daily_heart_rate_max)) +
  geom_line()


daily_heart_rate_data <-
  daily_heart_rate_data %>%
  do.call(rbind, .) %>%
  as.data.frame()


sample_info <-
  daily_heart_rate_data %>%
  dplyr::select(sample_id, subject_id, data_time)

expression_data <-
  daily_heart_rate_data %>%
  dplyr::select(daily_heart_rate_max,daily_heart_rate_min) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info$sample_id

variable_info <- data.frame(
  variable_id = c(
    "daily_heart_rate_max",
    "daily_heart_rate_min"
  )
)

sample_info$class <- "Subject"

library(massdataset)

daily_heart_rate_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

dir.create("3_data_analysis/1_data_preparation/wearable_data/4_daily_heart_rate", recursive = TRUE)
save(daily_heart_rate_data, file = "3_data_analysis/1_data_preparation/wearable_data/4_daily_heart_rate/daily_heart_rate_data.rda", compress = "xz")
