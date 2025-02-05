library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)


read_heartrate_file <- function(file_path) {
  # 读取原始数据
  data <- read.csv(file_path,
                   stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8")
  
  # 动态处理列数
  if (ncol(data) == 5) {
    # 跳过第一列（数据唯一ID）
    data <- data[, -1]
    colnames(data) <- c("数据时间", "平均心率", "测量时间", "外部ID")
  } else if (ncol(data) == 4) {
    colnames(data) <- c("数据时间", "平均心率", "测量时间", "外部ID")
  } else {
    stop("File:", file_path, "has unexpected number of columns:", ncol(data))
  }
  
  # 处理 measure_time 列
  if (is.numeric(data$"测量时间")) {
    cat("File:", file_path, "has numeric measure_time. Converting to POSIXct.\n")
    data$"测量时间" <- as.POSIXct(data$"测量时间" / 1000, origin = "1970-01-01")  # 处理毫秒级时间戳
  } else if (is.character(data$"测量时间")) {
    cat("File:", file_path, "has character measure_time. Converting to POSIXct.\n")
    data$"测量时间" <- as.POSIXct(data$"测量时间", format = "%Y-%m-%d %H:%M:%S")
  } else {
    cat("File:", file_path, "has unknown measure_time type:", class(data$"测量时间"), "\n")
  }
  
  return(data)
}



# all_subjects
subject_folders <- list.dirs("2_data/wearable data_1", recursive = FALSE)

heart_rate_data <-
  purrr::map(subject_folders, function(folder) {
    subject_id <- basename(folder)
    cat("Processing:", subject_id, "\n")
    
    hr_files <- list.files(
      folder,
      pattern = "t_apvddepp_continuousheartrate.*\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(hr_files) > 0) {
      # 添加调试信息
      cat("Found", length(hr_files), "files for", subject_id, "\n")
      
      subject_data <- purrr::map_df(hr_files, read_heartrate_file)
      cat("Read", nrow(subject_data), "rows for", subject_id, "\n")
      
      colnames(subject_data) <- c("data_time", "heart_rate", "measure_time", "subject_id")
      
      subject_data <-
        subject_data %>%
        dplyr::mutate(subject_id = subject_id) %>%  # 确保使用文件夹名作为subject_id
        dplyr::select(subject_id, measure_time, heart_rate) %>%
        dplyr::mutate(
          sample_id = paste(subject_id, measure_time, sep = "_"),
          measure_time = as.POSIXct(measure_time)
        ) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      # 添加更多调试信息
      cat("Processed", nrow(subject_data), "unique records for", subject_id, "\n")
      
      return(subject_data)
    } else {
      cat("No files found for", subject_id, "\n")
      return(NULL)
    }
  })

# 检查每个文件夹的处理结果
for(i in seq_along(heart_rate_data)) {
  folder_name <- basename(subject_folders[i])
  if(!is.null(heart_rate_data[[i]])) {
    cat(folder_name, ": ", nrow(heart_rate_data[[i]]), " rows\n")
  } else {
    cat(folder_name, ": NULL\n")
  }
}


heart_rate_data[[49]] %>%
  head(1000) %>%
  ggplot(aes(measure_time, heart_rate)) +
  geom_line()


heart_rate_data <-
  heart_rate_data %>%
  do.call(rbind, .) %>%
  as.data.frame()

sample_info <-
  heart_rate_data %>%
  dplyr::select(sample_id, subject_id, measure_time) %>%
  distinct(sample_id, .keep_all = TRUE)

# 修复sample_info
sample_info_fixed <- sample_info %>%
  mutate(
    subject_id = toupper(subject_id),  # 统一转换为大写
    subject_id = gsub("SHH|ShH", "SH", subject_id),  # 修复ShH061的问题
    measure_time = substr(measure_time, 1, 19),  # 保留到秒
    sample_id = paste(subject_id, measure_time, sep = "_")  
  )

expression_data <-
  heart_rate_data %>%
  dplyr::select(heart_rate) %>%
  t() %>%
  as.data.frame()



colnames(expression_data) <-
  sample_info_fixed$sample_id

expression_data <- expression_data[, sample_info_fixed$sample_id]


variable_info <- data.frame(variable_id = "heart_rate")


sample_info_fixed$class <- "Subject"

# mass_dataset
heart_rate_data <- create_mass_dataset(
  expression_data = expression_data,
  sample_info = sample_info_fixed,
  variable_info = variable_info
)





dir.create("3_data_analysis/1_data_preparation/wearable_data/1_heart_rate", recursive = TRUE)
save(heart_rate_data, file = "3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda", compress = "xz")
