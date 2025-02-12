library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)

# Modified function to read sleep data files with all columns
read_sleep_file <- function(file_path) {
  data <- read.csv(file_path,
                   stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8")
  
  # Check if the first column contains numeric timestamps or IDs
  first_col <- data[[1]]
  if (!any(grepl("-", first_col))) {
    # If the first column doesn't contain date format (no "-"), skip it
    data <- data[-1]
  }
  
  colnames(data) <- c(
    "数据时间",
    "最早入睡时间",
    "最迟出睡时间",
    "浅睡时长(min)",
    "深睡时长(min)",
    "做梦时长(min)",
    "清醒时长(min)",
    "全部睡眠时长(min)",
    "白天睡眠时长",
    "睡眠得分",
    "外部ID"
  )
  return(data)
}

# Process all subjects
subject_folders <- list.dirs("2_data/wearable data_1", recursive = FALSE)
sleep_data <-
  purrr::map(subject_folders, function(folder) {
    subject_id <- basename(folder)
    cat(subject_id, " ")
    
    sleep_files <- list.files(
      folder,
      pattern = "t_apvddepp_sleep.*\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(sleep_files) > 0) {
      subject_data <- map_df(sleep_files, read_sleep_file)
      
      # Add subject_id before renaming columns
      subject_data$subject_id <- subject_id
      
      colnames(subject_data) <- c(
        "measure_time",
        "sleep_start_time",
        "sleep_end_time",
        "light_sleep_duration",
        "deep_sleep_duration",
        "dream_sleep_duration",
        "awake_duration",
        "total_sleep_duration",
        "daytime_sleep_duration",
        "sleep_score",
        "external_id",
        "subject_id"
      )
      
      subject_data <-
        subject_data %>%
        dplyr::select(
          subject_id,
          measure_time,
          sleep_start_time,
          sleep_end_time,
          light_sleep_duration,
          deep_sleep_duration,
          dream_sleep_duration,
          awake_duration,
          total_sleep_duration,
          daytime_sleep_duration,
          sleep_score
        ) %>%
        dplyr::mutate(sample_id = paste(subject_id, measure_time, sep = "_")) %>%
        dplyr::mutate(
          measure_time = as.POSIXct(measure_time, format = "%Y-%m-%d %H:%M:%S"),
          sleep_start_time = as.POSIXct(sleep_start_time, format = "%Y-%m-%d %H:%M:%S"),
          sleep_end_time = as.POSIXct(sleep_end_time, format = "%Y-%m-%d %H:%M:%S")
        ) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      return(subject_data)
    } else{
      return(NULL)
    }
  })

# Visualization example - you can modify this based on what you want to visualize
sleep_data[[2]] %>%
  head(1000) %>%
  ggplot(aes(measure_time, total_sleep_duration)) +
  geom_line() +
  labs(title = "Total Sleep Duration Over Time",
       x = "Date",
       y = "Sleep Duration (minutes)")

# Combine all data
sleep_data <-
  sleep_data %>%
  do.call(rbind, .) %>%
  as.data.frame()

# Create sample info
sample_info <-
  sleep_data %>%
  dplyr::select(sample_id, subject_id, measure_time, sleep_start_time, sleep_end_time)

# Create expression data with all numeric sleep metrics
expression_data <-
  sleep_data %>%
  dplyr::select(
    light_sleep_duration,
    deep_sleep_duration,
    dream_sleep_duration,
    awake_duration,
    total_sleep_duration,
    daytime_sleep_duration,
    sleep_score
  ) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <- sample_info$sample_id

# Create variable info for all sleep metrics
variable_info <- data.frame(
  variable_id = c(
    "light_sleep_duration",
    "deep_sleep_duration",
    "dream_sleep_duration",
    "awake_duration",
    "total_sleep_duration",
    "daytime_sleep_duration",
    "sleep_score"
  )
)

sample_info$class <- "Subject"

# Create mass dataset
library(massdataset)
sleep_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

sleep_data@sample_info$measure_time <-
  as.POSIXct(sleep_data@sample_info$measure_time, tz = "Asia/Shanghai")

sleep_data <-
  sleep_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(subject_id, measure_time)


# Create directory and save data
dir.create("3_data_analysis/1_data_preparation/wearable_data/7_sleep", recursive = TRUE)
save(sleep_data, 
     file = "3_data_analysis/1_data_preparation/wearable_data/7_sleep/sleep_data.rda", 
     compress = "xz")

