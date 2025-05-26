library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyr)
library(dplyr)

# 改进的每日运动详情文件读取函数
read_dailyworkout_details_file <- function(file_path) {
  # 读取原始数据
  tryCatch({
    data <- read.csv(file_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    
    # 检查第一列是否包含时间戳或ID
    first_col <- data[[1]]
    if (!any(grepl("-", first_col))) {
      # 如果第一列不包含日期格式（没有"-"），跳过它
      data <- data[, -1]
    }
    
    # 根据列数分配列名
    if(ncol(data) == 7) {
      colnames(data) <- c("数据时间", "步数(steps)", "运动距离(m)", "卡路里", "累计爬升高度", "心率", "外部ID")
      data$'活动名称' <- NA
    } else if(ncol(data) == 8) {
      colnames(data) <- c("数据时间", "活动名称", "步数(steps)", "运动距离(m)", "卡路里", "累计爬升高度", "心率", "外部ID")
    } else {
      stop(paste("Unexpected number of columns:", ncol(data)))
    }
    
    # 处理时间列，确保正确格式化
    if ("数据时间" %in% colnames(data)) {
      # 检查时间列的数据类型和格式
      if (is.numeric(data$"数据时间")) {
        # 可能是Unix时间戳（毫秒）
        cat("File:", file_path, "has numeric time. Converting to POSIXct.\n")
        data$"数据时间" <- as.POSIXct(data$"数据时间" / 1000, 
                                  origin = "1970-01-01", 
                                  tz = "Asia/Shanghai")
      } else if (is.character(data$"数据时间")) {
        # 字符型，尝试不同格式
        cat("File:", file_path, "has character time. Converting to POSIXct.\n")
        
        # 确定可能的时间格式
        time_sample <- na.omit(data$"数据时间")[1]
        cat("Sample time string:", time_sample, "\n")
        
        # 根据字符串长度尝试不同格式
        if (nchar(time_sample) == 10) {
          # 可能只有日期 YYYY-MM-DD
          data$"数据时间" <- as.POSIXct(data$"数据时间", 
                                    format = "%Y-%m-%d", 
                                    tz = "Asia/Shanghai")
        } else {
          # 尝试完整的日期时间格式
          data$"数据时间" <- as.POSIXct(data$"数据时间", 
                                    format = "%Y-%m-%d %H:%M:%S", 
                                    tz = "Asia/Shanghai")
        }
      }
      
      # 验证转换结果
      na_count <- sum(is.na(data$"数据时间"))
      if (na_count > 0) {
        warning("File:", file_path, "has", na_count, "NA values in time column after conversion")
      }
    }
    
    return(data)
  }, error = function(e) {
    warning("Error processing file:", file_path, "\nError:", e$message)
    return(NULL)
  })
}

# 分析所有受试者
subject_folders <- list.dirs("2_data/wearable data_2", recursive = FALSE)

# 进度跟踪函数
process_with_progress <- function(folders, pattern, process_func) {
  total <- length(folders)
  cat("开始处理", total, "个受试者...\n")
  
  results <- list()
  
  for (i in seq_along(folders)) {
    folder <- folders[i]
    subject_id <- basename(folder)
    cat("[", i, "/", total, "] 处理受试者:", subject_id, "\n")
    
    files <- list.files(
      folder,
      pattern = pattern,
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(files) > 0) {
      cat("  找到", length(files), "个文件\n")
      result <- process_func(subject_id, files)
      results[[i]] <- result
      
      if (!is.null(result)) {
        cat("  成功处理了", nrow(result), "条记录\n")
      } else {
        cat("  处理结果为空\n")
      }
    } else {
      cat("  没有找到匹配的文件\n")
      results[[i]] <- NULL
    }
  }
  
  cat("所有受试者处理完成\n")
  return(results)
}

# 处理每个受试者的函数
process_subject_data <- function(subject_id, files) {
  # 读取并合并该受试者的所有文件
  subject_data_list <- lapply(files, read_dailyworkout_details_file)
  subject_data_list <- subject_data_list[!sapply(subject_data_list, is.null)]
  
  if (length(subject_data_list) == 0) {
    return(NULL)
  }
  
  subject_data <- do.call(rbind, subject_data_list)
  
  # 添加受试者ID
  subject_data$subject_id <- subject_id
  
  # 重命名列
  colnames(subject_data) <- c("measure_time", "activity", "steps", "distance", "calorie", "climbing_height", "heartrate", "external_id", "subject_id")
  
  # 处理数据
  subject_data <- subject_data %>%
    dplyr::select(subject_id, measure_time, activity, steps, distance, calorie, climbing_height, heartrate) %>%
    dplyr::mutate(
      # 确保measure_time是POSIXct类型
      measure_time = if(class(measure_time)[1] != "POSIXct") {
        as.POSIXct(measure_time, tz = "Asia/Shanghai")
      } else {
        measure_time
      },
      # 创建规范的sample_id
      sample_id = paste(subject_id, format(measure_time, "%Y-%m-%d %H:%M:%S"), sep = "_")
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  # 验证时间格式
  time_check <- table(format(subject_data$measure_time, "%H:%M:%S") != "00:00:00")
  cat("  时间格式检查: 有完整时间的记录:", time_check["TRUE"], 
      "只有日期的记录:", time_check["FALSE"], "\n")
  
  return(subject_data)
}

# 处理所有受试者数据
daily_workout_details_data <- process_with_progress(
  subject_folders,
  "t_apvddepp_dailyworkoutdetail.*\\.txt$",
  process_subject_data
)

# 合并所有受试者数据
daily_workout_details_data <- do.call(rbind, daily_workout_details_data)

# 检查合并后的数据
cat("合并后的数据行数:", nrow(daily_workout_details_data), "\n")

# 创建样本信息
sample_info <- daily_workout_details_data %>%
  dplyr::select(sample_id, subject_id, measure_time, activity)

# 修复样本信息，保留完整时间
sample_info_fixed <- sample_info %>%
  mutate(
    subject_id = toupper(subject_id),  # 统一转换为大写
    subject_id = gsub("SHH|ShH", "SH", subject_id),  # 修复ShH061的问题
    # 确保measure_time是POSIXct类型，不截断时间
    measure_time = as.POSIXct(measure_time, tz = "Asia/Shanghai"),
    # 重新创建sample_id以确保一致性
    sample_id = paste(subject_id, format(measure_time, "%Y-%m-%d %H:%M:%S"), sep = "_")
  )

# 创建表达数据
expression_data <- daily_workout_details_data %>%
  dplyr::select(steps, distance, calorie, climbing_height, heartrate) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <- sample_info_fixed$sample_id
expression_data <- expression_data[, sample_info_fixed$sample_id]

# 创建变量信息
variable_info <- data.frame(
  variable_id = c(
    "steps",
    "distance",
    "calorie",
    "climbing_height",
    "heartrate"
  )
)

# 设置样本类别
sample_info_fixed$class <- "Subject"

# 创建质谱数据集
library(massdataset)

daily_workout_details_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info_fixed,
    variable_info = variable_info
  )

# 最终验证时间格式
final_time_check <- daily_workout_details_data@sample_info %>%
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

cat("最终时间格式检查:\n")
print(final_time_check)

# 按受试者ID和测量时间排序
daily_workout_details_data <-
  daily_workout_details_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(subject_id, measure_time)

# 检查每个受试者的时间分布
time_distribution_check <- daily_workout_details_data@sample_info %>%
  mutate(hour = hour(measure_time)) %>%
  group_by(subject_id, hour) %>%
  summarise(count = n(), .groups = "drop")

# 显示部分结果
cat("时间分布检查 (部分样本):\n")
print(head(time_distribution_check, 20))

# 创建目录并保存数据
dir.create("3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details", recursive = TRUE)
save(daily_workout_details_data, file = "3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details/daily_workout_details_data.rda", compress = "xz")

# 保存时间检查结果，便于后续参考
save(final_time_check, time_distribution_check, 
     file = "3_data_analysis/1_data_preparation/wearable_data/5_daily_workout_details/time_validation.rda")

cat("数据处理和保存完成！\n")