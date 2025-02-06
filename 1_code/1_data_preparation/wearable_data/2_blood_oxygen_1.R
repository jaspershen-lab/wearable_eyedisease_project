library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)


read_bloodoxygen_file <- function(file_path) {
  data <- read.csv(file_path,
                   stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8")
  
  # Check if the first column contains numeric timestamps or IDs
  first_col <- data[[1]]
  if (!any(grepl("-", first_col))) {  # Changed from "/" to "-" to match the format
    # If the first column doesn't contain date format (no "-"), skip it
    data <- data[-1]
  }
  
  # Check number of columns and assign names accordingly
  if (ncol(data) == 5) {
    # 跳过第一列（数据唯一ID）
    data <- data[, -1]
    colnames(data) <- c("数据时间", "平均血氧值(%)", "平均血氧值.测量时间", "外部ID")
  } else if (ncol(data) == 4) {
    colnames(data) <- c("数据时间", "平均血氧值(%)", "平均血氧值.测量时间", "外部ID")
  } else {
    stop(paste("Unexpected number of columns:", ncol(data)))
  }
  
  return(data)
}

# all_subjects
subject_folders <- list.dirs("2_data/wearable data_1", recursive = FALSE)
blood_oxygen_data <-
  purrr::map(subject_folders, function(folder) {
    subject_id <- basename(folder)
    cat(subject_id, " ")
    
    hr_files <- list.files(
      folder,
      pattern = "t_apvddepp_continuousbloodoxygensaturation.*\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(hr_files) > 0) {
      subject_data <- map_df(hr_files, read_bloodoxygen_file)
      
      # Adjust column names based on number of columns
      if (ncol(subject_data) == 4) {
        colnames(subject_data) <- c("data_time", "blood_oxygen", "measure_time", "subject_id")
        subject_data <- subject_data %>%
          dplyr::select(subject_id, measure_time, blood_oxygen)
      } else {
        colnames(subject_data) <- c("data_time", "blood_oxygen", "measure_time")
        subject_data <- subject_data %>%
          dplyr::mutate(subject_id = subject_id)  # Add subject_id column
      }
      
      subject_data <- subject_data %>%
        filter(!is.na(measure_time)) %>%
        dplyr::mutate(sample_id = paste(subject_id, measure_time, sep = "_")) %>%
        dplyr::mutate(measure_time = as.POSIXct(measure_time, format = "%Y-%m-%d %H:%M:%S")) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      return(subject_data)
    } else{
      return(NULL)
    }
  })


blood_oxygen_data[[2]] %>%
  head(1000) %>%
  ggplot(aes(measure_time, blood_oxygen)) +
  geom_line()


blood_oxygen_data <-
  blood_oxygen_data %>%
  do.call(rbind, .) %>%
  as.data.frame()


sample_info <-
  blood_oxygen_data %>%
  dplyr::select(sample_id, subject_id, measure_time)

sample_info_fixed <- sample_info %>%
  mutate(
    subject_id = toupper(subject_id),  # 统一转换为大写
    subject_id = gsub("SHH|ShH", "SH", subject_id),  # 修复ShH061的问题
    measure_time = substr(measure_time, 1, 19),  # 保留到秒
    sample_id = paste(subject_id, measure_time, sep = "_")  
  )

sample_info_fixed <- sample_info_fixed %>%
  filter(subject_id != "")

expression_data <-
  blood_oxygen_data %>%
  dplyr::select(blood_oxygen) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info_fixed$sample_id

expression_data <- expression_data[, sample_info_fixed$sample_id]

variable_info <-
  data.frame(variable_id = "blood_oxygen")

sample_info_fixed$class <- "Subject"

library(massdataset)

blood_oxygen_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info_fixed,
    variable_info = variable_info
  )

blood_oxygen_data@sample_info$measure_time <-
  as.POSIXct(blood_oxygen_data@sample_info$measure_time, tz = "Asia/Shanghai")

blood_oxygen_data <-
  blood_oxygen_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(subject_id, measure_time)

dir.create("3_data_analysis/1_data_preparation/wearable_data/2_blood_oxygen", recursive = TRUE)
save(blood_oxygen_data, file = "3_data_analysis/1_data_preparation/wearable_data/2_blood_oxygen/blood_oxygen_data.rda", compress = "xz")
