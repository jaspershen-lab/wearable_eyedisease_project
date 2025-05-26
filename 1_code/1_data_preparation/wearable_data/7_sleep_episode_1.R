library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)

# 改进的睡眠数据文件读取函数
read_sleep_file <- function(file_path) {
  tryCatch({
    # 读取原始数据
    data <- read.csv(file_path,
                     stringsAsFactors = FALSE,
                     fileEncoding = "UTF-8")
    
    # 检查第一列是否包含时间戳或ID
    first_col <- data[[1]]
    if (!any(grepl("-", first_col))) {
      # 如果第一列不包含日期格式（没有"-"），跳过它
      data <- data[-1]
    }
    
    # 分配列名
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
    
    # 处理所有时间列
    time_columns <- c("数据时间", "最早入睡时间", "最迟出睡时间")
    
    for (col in time_columns) {
      if (col %in% colnames(data)) {
        # 检查时间列的数据类型
        if (is.numeric(data[[col]])) {
          # 可能是Unix时间戳（毫秒）
          cat("File:", file_path, "column", col, "has numeric time. Converting to POSIXct.\n")
          data[[col]] <- as.POSIXct(data[[col]] / 1000, 
                                    origin = "1970-01-01", 
                                    tz = "Asia/Shanghai")
        } else if (is.character(data[[col]])) {
          # 字符型，尝试不同格式
          cat("File:", file_path, "column", col, "has character time. Converting to POSIXct.\n")
          
          # 取非NA值的样本
          sample_values <- na.omit(data[[col]])
          if (length(sample_values) > 0) {
            sample_time <- sample_values[1]
            cat("Sample time string for", col, ":", sample_time, "\n")
            
            # 根据字符串长度尝试不同格式
            if (nchar(sample_time) == 10) {
              # 可能只有日期 YYYY-MM-DD
              data[[col]] <- as.POSIXct(data[[col]], 
                                        format = "%Y-%m-%d", 
                                        tz = "Asia/Shanghai")
            } else {
              # 尝试完整的日期时间格式
              data[[col]] <- as.POSIXct(data[[col]], 
                                        format = "%Y-%m-%d %H:%M:%S", 
                                        tz = "Asia/Shanghai")
            }
          }
        }
        
        # 验证转换结果
        na_count <- sum(is.na(data[[col]]))
        if (na_count > 0) {
          warning("File:", file_path, "column", col, "has", na_count, "NA values after conversion")
        }
      }
    }
    
    return(data)
  }, error = function(e) {
    warning("Error processing file:", file_path, "\nError:", e$message)
    return(NULL)
  })
}

# 处理所有受试者
subject_folders <- list.dirs("2_data/wearable data_2", recursive = FALSE)

# 创建进度跟踪函数
process_with_progress <- function(folders, pattern, function_name) {
  total <- length(folders)
  result_list <- vector("list", total)
  
  cat("开始处理", total, "个受试者的睡眠数据...\n")
  
  for (i in seq_along(folders)) {
    folder <- folders[i]
    subject_id <- basename(folder)
    cat("[", i, "/", total, "] 处理受试者:", subject_id, "\n")
    
    # 查找匹配的文件
    sleep_files <- list.files(
      folder,
      pattern = pattern,
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(sleep_files) > 0) {
      cat("  找到", length(sleep_files), "个睡眠数据文件\n")
      
      # 读取并合并数据
      subject_data_list <- lapply(sleep_files, read_sleep_file)
      subject_data_list <- subject_data_list[!sapply(subject_data_list, is.null)]
      
      if (length(subject_data_list) == 0) {
        cat("  所有文件处理失败，跳过该受试者\n")
        result_list[[i]] <- NULL
        next
      }
      
      subject_data <- do.call(rbind, subject_data_list)
      
      # 添加受试者ID
      subject_data$subject_id <- subject_id
      
      # 重命名列
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
      
      # 处理数据
      subject_data <- subject_data %>%
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
        dplyr::mutate(
          # 确保所有时间列都是POSIXct类型并有正确的时区
          measure_time = if(class(measure_time)[1] != "POSIXct") {
            as.POSIXct(measure_time, tz = "Asia/Shanghai")
          } else {
            measure_time
          },
          sleep_start_time = if(class(sleep_start_time)[1] != "POSIXct") {
            as.POSIXct(sleep_start_time, tz = "Asia/Shanghai")
          } else {
            sleep_start_time
          },
          sleep_end_time = if(class(sleep_end_time)[1] != "POSIXct") {
            as.POSIXct(sleep_end_time, tz = "Asia/Shanghai")
          } else {
            sleep_end_time
          },
          # 创建sample_id，使用完整的时间格式
          sample_id = paste(subject_id, format(measure_time, "%Y-%m-%d %H:%M:%S"), sep = "_")
        ) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE)
      
      # 验证时间格式
      time_formats <- subject_data %>%
        summarise(
          measure_time_complete = sum(format(measure_time, "%H:%M:%S") != "00:00:00"),
          start_time_complete = sum(format(sleep_start_time, "%H:%M:%S") != "00:00:00"),
          end_time_complete = sum(format(sleep_end_time, "%H:%M:%S") != "00:00:00"),
          total_records = n()
        )
      
      cat("  时间格式检查:\n",
          "    数据时间 - 完整时间记录:", time_formats$measure_time_complete, "/", time_formats$total_records, "\n",
          "    入睡时间 - 完整时间记录:", time_formats$start_time_complete, "/", time_formats$total_records, "\n",
          "    出睡时间 - 完整时间记录:", time_formats$end_time_complete, "/", time_formats$total_records, "\n")
      
      result_list[[i]] <- subject_data
      cat("  成功处理", nrow(subject_data), "条记录\n")
    } else {
      cat("  没有找到睡眠数据文件\n")
      result_list[[i]] <- NULL
    }
  }
  
  cat("所有受试者处理完成\n")
  return(result_list)
}

# 处理所有睡眠数据
sleep_data <- process_with_progress(
  subject_folders,
  "t_apvddepp_sleep.*\\.txt$",
  read_sleep_file
)

# 可视化示例 - 可以根据需要修改
if (!is.null(sleep_data[[2]])) {
  plot <- sleep_data[[2]] %>%
    head(1000) %>%
    ggplot(aes(measure_time, total_sleep_duration)) +
    geom_line() +
    labs(title = "Total Sleep Duration Over Time",
         x = "Date",
         y = "Sleep Duration (minutes)")
  
  print(plot)
}

# 合并所有数据
non_null_data <- sleep_data[!sapply(sleep_data, is.null)]
if (length(non_null_data) > 0) {
  sleep_data_combined <- do.call(rbind, non_null_data)
  cat("合并后的睡眠数据包含", nrow(sleep_data_combined), "条记录\n")
} else {
  stop("没有有效的睡眠数据可以合并！")
}

# 创建样本信息
sample_info <- sleep_data_combined %>%
  dplyr::select(sample_id, subject_id, measure_time, sleep_start_time, sleep_end_time)

# 修复样本信息，标准化受试者ID
sample_info_fixed <- sample_info %>%
  mutate(
    subject_id = toupper(subject_id),  # 统一转换为大写
    subject_id = gsub("SHH|ShH", "SH", subject_id),  # 修复ShH061的问题
    # 确保所有时间都是POSIXct类型并有正确的时区
    measure_time = as.POSIXct(measure_time, tz = "Asia/Shanghai"),
    sleep_start_time = as.POSIXct(sleep_start_time, tz = "Asia/Shanghai"),
    sleep_end_time = as.POSIXct(sleep_end_time, tz = "Asia/Shanghai"),
    # 使用标准化格式重新创建sample_id
    sample_id = paste(subject_id, format(measure_time, "%Y-%m-%d %H:%M:%S"), sep = "_")
  )

# 创建表达数据
expression_data <- sleep_data_combined %>%
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

# 确保列名与sample_info_fixed中的sample_id一致
colnames(expression_data) <- sample_info_fixed$sample_id

# 创建变量信息
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

# 设置样本类别
sample_info_fixed$class <- "Subject"

# 创建质谱数据集
library(massdataset)
sleep_data <- create_mass_dataset(
  expression_data = expression_data,
  sample_info = sample_info_fixed,
  variable_info = variable_info
)

# 最终验证时间格式
final_time_check <- sleep_data@sample_info %>%
  mutate(
    measure_time_complete = format(measure_time, "%H:%M:%S") != "00:00:00",
    start_time_complete = format(sleep_start_time, "%H:%M:%S") != "00:00:00",
    end_time_complete = format(sleep_end_time, "%H:%M:%S") != "00:00:00"
  ) %>%
  group_by(subject_id) %>%
  summarise(
    total_records = n(),
    measure_time_complete_count = sum(measure_time_complete),
    start_time_complete_count = sum(start_time_complete),
    end_time_complete_count = sum(end_time_complete),
    measure_time_complete_percent = (measure_time_complete_count / total_records) * 100,
    start_time_complete_percent = (start_time_complete_count / total_records) * 100,
    end_time_complete_percent = (end_time_complete_count / total_records) * 100
  )

cat("最终时间格式验证:\n")
print(final_time_check)

# 按受试者ID和测量时间排序
sleep_data <- sleep_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(subject_id, measure_time)

# 创建目录并保存数据
dir.create("3_data_analysis/1_data_preparation/wearable_data/7_sleep", recursive = TRUE)

# 保存主要数据集
save(sleep_data, 
     file = "3_data_analysis/1_data_preparation/wearable_data/7_sleep/sleep_data.rda", 
     compress = "xz")

# 保存时间验证结果
save(final_time_check,
     file = "3_data_analysis/1_data_preparation/wearable_data/7_sleep/time_validation.rda")

cat("睡眠数据处理和保存完成！\n")