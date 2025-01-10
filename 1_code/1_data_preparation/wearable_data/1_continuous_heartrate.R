library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)


read_heartrate_file <- function(file_path) {
  data <- read.csv(file_path,
                   stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8")
  colnames(data) <- c("数据时间", "平均心率", "测量时间", "外部ID")
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data", recursive = FALSE)

heart_rate_data <-
  purrr::map(subject_folders, function(folder) {
    subject_id <- basename(folder)
    cat(subject_id, " ")
    
    hr_files <- list.files(
      folder,
      pattern = "t_apvddepp_continuousheartrate.*\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(hr_files) > 0) {
      subject_data <- map_df(hr_files, read_heartrate_file)
      colnames(subject_data) <-
        c("data_time", "heart_rate", "measure_time", "subject_id")
      subject_data <-
        subject_data %>%
        dplyr::select(subject_id, measure_time, heart_rate)
      
      subject_data <-
        subject_data %>%
        dplyr::mutate(sample_id = paste(subject_id, measure_time, sep = "_")) %>%
        dplyr::mutate(measure_time = as.POSIXct(measure_time)) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      return(subject_data)
    } else{
      return(NULL)
    }
  })


heart_rate_data[[2]] %>%
  head(1000) %>%
  ggplot(aes(measure_time, heart_rate)) +
  geom_line()


heart_rate_data <-
  heart_rate_data %>%
  do.call(rbind, .) %>%
  as.data.frame()


sample_info <-
  heart_rate_data %>%
  dplyr::select(sample_id, subject_id, measure_time)

expression_data <-
  heart_rate_data %>%
  dplyr::select(heart_rate) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info$sample_id

variable_info <-
  data.frame(variable_id = "heart_rate")

sample_info$class <- "Subject"

library(massdataset)

heart_rate_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )


dir.create("3_data_analysis/1_data_preparation/1_heart_rate", recursive = TRUE)
save(heart_rate_data, file = "3_data_analysis/1_data_preparation/1_heart_rate/heart_rate_data.rda", compress = "xz")



# library(tidyverse)
# library(tidymass)
# library(r4projects)
# setwd(get_project_wd())
# rm(list = ls())
# source('1_code/100_tools.R')
# library(tidyr)
# library(dplyr)
#
#
#
# process_heartrate_data <- function(root_dir = "2_data/wearable data") {
#   # 读取并处理单个心率文件的函数
#   read_heartrate_file <- function(file_path) {
#     data <- read.csv(file_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
#     colnames(data) <- c("数据时间", "平均心率", "测量时间", "外部ID")
#     return(data)
#   }
#
#   # 获取所有参与者文件夹
#   subject_folders <- list.dirs(root_dir, recursive = FALSE)
#
#   # 存储所有心率数据
#   all_heartrate_data <- list()
#
#   for (folder in subject_folders) {
#     subject_id <- basename(folder)
#
#     hr_files <- list.files(
#       folder,
#       pattern = "t_apvddepp_continuousheartrate.*\\.txt$",
#       recursive = TRUE,
#       full.names = TRUE
#     )
#
#     if (length(hr_files) > 0) {
#       subject_data <- map_df(hr_files, read_heartrate_file)
#       all_heartrate_data[[subject_id]] <- subject_data
#     }
#   }
#
#   # 合并所有数据并按时间排序
#   combined_data <- bind_rows(all_heartrate_data) %>%
#     arrange(数据时间)
#
#   # 检查数据结构
#   message("Combined data columns: ", paste(colnames(combined_data), collapse = ", "))
#
#   # 创建expression_data (转置的格式)
#   expression_data <- combined_data %>%
#     dplyr::select(dplyr::all_of(c("数据时间", "平均心率", "外部ID"))) %>%
#     dplyr::distinct(数据时间, 外部ID, .keep_all = TRUE) %>%
#     tidyr::pivot_wider(
#       id_cols = 数据时间,
#       names_from = 外部ID,
#       values_from = 平均心率,
#       values_fill = NA
#     ) %>%
#     tibble::column_to_rownames("数据时间")
#
#   # 创建sample_info
#   sample_info <- data.frame(
#     sample_id = colnames(expression_data),
#     class = "subject",
#     row.names = colnames(expression_data)
#   )
#
#   # 创建variable_info
#   variable_info <- data.frame(
#     variable_id = rownames(expression_data),
#     variable_type = "continuous_heartrate",
#     row.names = rownames(expression_data)
#   )
#
#   # 打印维度信息以验证
#   message("Dimensions check:")
#   message("expression_data: ", nrow(expression_data), " rows x ", ncol(expression_data), " columns")
#   message("sample_info: ", nrow(sample_info), " rows")
#   message("variable_info: ", nrow(variable_info), " rows")
#
#   # 创建mass_dataset
#   continuous_heartrate_data <- create_mass_dataset(
#     expression_data = expression_data,
#     sample_info = sample_info,
#     variable_info = variable_info
#   )
#
#   return(continuous_heartrate_data)
# }
#
# # 使用函数处理数据
# continuous_heartrate_data <- process_heartrate_data()
