library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)

# 改进的血氧文件读取函数，加强时间处理
read_bloodoxygen_file <- function(file_path) {
  # 读取原始数据
  data <- tryCatch({
    read.csv(file_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  }, error = function(e) {
    stop("Error reading file:", file_path, "\nError:", e$message)
  })
  
  # 检查第一列是否包含时间戳或ID
  first_col <- data[[1]]
  if (!any(grepl("-", first_col))) {
    # 如果第一列不包含日期格式（没有"-"），跳过它
    data <- data[, -1]
  }
  
  # 根据列数分配列名
  if (ncol(data) == 5) {
    # 跳过第一列（数据唯一ID）
    data <- data[, -1]
    colnames(data) <- c("数据时间", "平均血氧值(%)", "平均血氧值.测量时间", "外部ID")
  } else if (ncol(data) == 4) {
    colnames(data) <- c("数据时间", "平均血氧值(%)", "平均血氧值.测量时间", "外部ID")
  } else {
    stop(paste("Unexpected number of columns:", ncol(data)))
  }
  
  # 处理时间列，确保正确格式化
  if ("平均血氧值.测量时间" %in% colnames(data)) {
    tryCatch({
      # 检查时间列的数据类型
      if (is.numeric(data$"平均血氧值.测量时间")) {
        # 如果是数值型，可能是Unix时间戳（毫秒）
        cat("File:", file_path, "has numeric measure_time. Converting to POSIXct.\n")
        data$"平均血氧值.测量时间" <- as.POSIXct(data$"平均血氧值.测量时间" / 1000, 
                                        origin = "1970-01-01", 
                                        tz = "Asia/Shanghai")
      } else if (is.character(data$"平均血氧值.测量时间")) {
        # 如果是字符型，尝试不同格式
        cat("File:", file_path, "has character measure_time. Converting to POSIXct.\n")
        
        # 确定字符串长度来判断可能的格式
        time_length <- unique(nchar(na.omit(data$"平均血氧值.测量时间")))
        cat("Time string lengths:", paste(time_length, collapse=", "), "\n")
        
        # 尝试几种常见格式
        if (all(time_length == 10)) {
          # 可能只有日期 YYYY-MM-DD
          data$"平均血氧值.测量时间" <- as.POSIXct(data$"平均血氧值.测量时间", 
                                          format = "%Y-%m-%d", 
                                          tz = "Asia/Shanghai")
        } else {
          # 尝试完整格式
          data$"平均血氧值.测量时间" <- as.POSIXct(data$"平均血氧值.测量时间", 
                                          format = "%Y-%m-%d %H:%M:%S", 
                                          tz = "Asia/Shanghai")
        }
      }
      
      # 验证转换结果
      na_count <- sum(is.na(data$"平均血氧值.测量时间"))
      if (na_count > 0) {
        warning("File:", file_path, "has", na_count, "NA values in measure_time after conversion")
      }
    }, error = function(e) {
      warning("Error processing measure_time in file:", file_path, "\nError:", e$message)
    })
  }
  
  return(data)
}

# 所有受试者
subject_folders <- list.dirs("2_data/wearable data_2", recursive = FALSE)

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
      # 读取并处理数据
      subject_data <- map_df(hr_files, read_bloodoxygen_file)
      
      # 根据列数调整列名
      if (ncol(subject_data) == 4) {
        colnames(subject_data) <- c("data_time", "blood_oxygen", "measure_time", "subject_id")
        subject_data <- subject_data %>%
          dplyr::select(subject_id, measure_time, blood_oxygen)
      } else {
        colnames(subject_data) <- c("data_time", "blood_oxygen", "measure_time")
        subject_data <- subject_data %>%
          dplyr::mutate(subject_id = subject_id)  # 添加受试者ID列
      }
      
      # 额外的时间格式检查和处理
      subject_data <- subject_data %>%
        filter(!is.na(measure_time)) %>%
        dplyr::mutate(
          # 确保measure_time是POSIXct类型
          measure_time = if(class(measure_time)[1] != "POSIXct") {
            as.POSIXct(measure_time, format = "%Y-%m-%d %H:%M:%S", tz = "Asia/Shanghai")
          } else {
            measure_time
          },
          # 创建样本ID
          sample_id = paste(subject_id, format(measure_time, "%Y-%m-%d %H:%M:%S"), sep = "_")
        ) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      # 验证时间格式
      time_check <- subject_data %>%
        mutate(
          has_time = format(measure_time, "%H:%M:%S") != "00:00:00",
          time_only = format(measure_time, "%H:%M:%S")
        )
      
      time_stats <- table(time_check$has_time)
      cat("\n", subject_id, "time stats:", 
          "With time:", time_stats[TRUE], 
          "Without time:", time_stats[FALSE], "\n")
      
      return(subject_data)
    } else {
      return(NULL)
    }
  })

# 检查血氧数据的时间分布
check_subject_data <- function(data_list, n = 3) {
  for (i in 1:min(n, length(data_list))) {
    if (!is.null(data_list[[i]])) {
      cat("Subject", i, "data sample:\n")
      print(head(data_list[[i]], 3))
      
      time_distribution <- data_list[[i]] %>%
        mutate(hour = hour(measure_time)) %>%
        group_by(hour) %>%
        summarise(count = n())
      
      cat("Hour distribution:\n")
      print(time_distribution)
    }
  }
}

# 检查前3个受试者的数据
check_subject_data(blood_oxygen_data)

# 合并所有受试者数据
blood_oxygen_data <-
  blood_oxygen_data %>%
  do.call(rbind, .) %>%
  as.data.frame()

# 创建样本信息
sample_info <-
  blood_oxygen_data %>%
  dplyr::select(sample_id, subject_id, measure_time)

# 修复样本信息，避免截断时间
sample_info_fixed <- sample_info %>%
  mutate(
    subject_id = toupper(subject_id),  # 统一转换为大写
    subject_id = gsub("SHH|ShH", "SH", subject_id),  # 修复ShH061的问题
    # 确保measure_time是POSIXct类型而不是截断字符串
    measure_time = as.POSIXct(measure_time, tz = "Asia/Shanghai"),
    # 重新创建样本ID以确保一致
    sample_id = paste(subject_id, format(measure_time, "%Y-%m-%d %H:%M:%S"), sep = "_")
  ) %>%
  filter(subject_id != "")

# 创建表达数据
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

# 创建mass_dataset
blood_oxygen_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info_fixed,
    variable_info = variable_info
  )

# 最终验证：确保measure_time包含完整的时间信息
final_time_check <- blood_oxygen_data@sample_info %>%
  mutate(
    has_time = format(measure_time, "%H:%M:%S") != "00:00:00",
    date_only = format(measure_time, "%Y-%m-%d"),
    time_only = format(measure_time, "%H:%M:%S")
  ) %>%
  group_by(subject_id) %>%
  summarise(
    total_records = n(),
    records_with_time = sum(has_time),
    percent_with_time = (records_with_time / total_records) * 100,
    unique_dates = n_distinct(date_only),
    unique_times = n_distinct(time_only)
  )

print("Final time check:")
print(final_time_check)

# 按受试者ID和测量时间排序
blood_oxygen_data <-
  blood_oxygen_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(subject_id, measure_time)

# 创建目录并保存数据
dir.create("3_data_analysis/1_data_preparation/wearable_data/2_blood_oxygen", recursive = TRUE)
save(blood_oxygen_data, file = "3_data_analysis/1_data_preparation/wearable_data/2_blood_oxygen/blood_oxygen_data.rda", compress = "xz")