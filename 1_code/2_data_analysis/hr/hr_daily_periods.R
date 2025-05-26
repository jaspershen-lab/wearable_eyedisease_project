library(tidyverse)
library(tidymass)
library(r4projects)
library(moments)  # 添加moments包用于计算偏度和峰度
setwd(get_project_wd())
rm(list = ls())
library(lubridate)

######## 读取数据
load(
  "3_data_analysis/1_data_preparation/wearable_data/1_heart_rate/heart_rate_data.rda"
)

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#############
dir.create("3_data_analysis/2_data_analysis/hr/hr_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/hr/hr_time_periods")

###########
calculate_daily_heart_rate <- function(data, baseline_info) {
  # 设置固定的时间范围
  PRE_SURGERY_DAYS <- -10
  POST_SURGERY_DAYS <- 90
  
  # 从心率数据获取受试者ID
  heart_rate_subjects <- data@sample_info %>%
    dplyr::select(subject_id) %>%
    distinct() %>%
    pull(subject_id)
  
  # 处理基线信息
  baseline_info_processed <- baseline_info %>%
    dplyr::filter(ID %in% heart_rate_subjects) %>%
    mutate(
      surgery_time_1 = as.POSIXct(surgery_time_1, format = "%Y-%m-%d")
    ) %>%
    dplyr::select(ID, surgery_time_1)
  
  # 提取和处理心率数据
  process_heart_rate_data <- function() {
    heart_rate_values <- data@expression_data[1, ] %>%
      as.numeric()
    
    data.frame(
      sample_id = colnames(data@expression_data),
      heart_rate = heart_rate_values,
      timestamp = data@sample_info$measure_time,
      subject_id = data@sample_info$subject_id
    ) %>%
      dplyr::filter(
        heart_rate >= 30 & heart_rate <= 220,  # 有效心率范围
        subject_id %in% baseline_info_processed$ID
      )
  }
  
  # 计算每日统计数据
  calculate_daily_stats <- function(heart_rate_data) {
    heart_rate_data %>%
      left_join(baseline_info_processed, by = c("subject_id" = "ID")) %>%
      mutate(
        days_to_surgery = as.numeric(difftime(timestamp, surgery_time_1, units = "days")),
        day_point = floor(days_to_surgery)
      ) %>%
      dplyr::filter(
        day_point >= PRE_SURGERY_DAYS,
        day_point <= POST_SURGERY_DAYS
      ) %>%
      group_by(subject_id, day_point) %>%
      summarise(
        daily_hr_mean = mean(heart_rate),
        daily_hr_min = min(heart_rate),
        daily_hr_max = max(heart_rate),
        daily_hr_median = median(heart_rate),
        daily_hr_sd = sd(heart_rate),
        daily_hr_cv = (sd(heart_rate) / mean(heart_rate)) * 100, # CV以百分比表示
        daily_hr_iqr = IQR(heart_rate),
        daily_hr_skew = skewness(heart_rate),  # 添加偏度计算
        daily_hr_kurt = kurtosis(heart_rate),  # 添加峰度计算
        daily_hr_rmssd = sqrt(mean(diff(heart_rate)^2)),  # 添加RMSSD（连续RR间隔差值的均方根）
        n_measurements = n(),
        .groups = "drop"
      )
  }
  
  # 计算统计数据
  daily_stats <- calculate_daily_stats(process_heart_rate_data())
  
  # 为每个统计数据创建宽格式
  create_wide_format <- function(daily_stats) {
    stats_list <- list()
    
    # 要处理的统计数据列表
    stat_cols <- c("mean", "min", "max", "median", "sd", "cv", "iqr", "skew", "kurt", "rmssd")
    
    for(stat in stat_cols) {
      col_name <- paste0("daily_hr_", stat)
      stats_list[[stat]] <- daily_stats %>%
        mutate(col_name = paste0("day_", day_point, "_", stat, "_hr")) %>%
        dplyr::select(subject_id, col_name, !!sym(col_name)) %>%
        pivot_wider(
          names_from = col_name,
          values_from = !!sym(col_name)
        )
    }
    
    # 合并所有统计数据
    result <- stats_list[[1]]
    for(i in 2:length(stats_list)) {
      result <- result %>%
        left_join(stats_list[[i]], by = "subject_id")
    }
    
    return(result)
  }
  
  # 处理数据
  hr_wide <- create_wide_format(daily_stats)
  
  # 创建最终结果
  result <- tibble(subject_id = unique(baseline_info_processed$ID)) %>%
    left_join(
      baseline_info_processed %>%
        mutate(surgery_date = as.Date(surgery_time_1)) %>%
        dplyr::select(ID, surgery_date),
      by = c("subject_id" = "ID")
    ) %>%
    left_join(hr_wide, by = "subject_id") %>%
    distinct() %>%
    arrange(subject_id)
  
  # 对列进行排序
  col_order <- c("subject_id", "surgery_date")
  time_cols <- setdiff(names(result), col_order)
  sorted_time_cols <- time_cols[order(as.numeric(gsub(".*day_(.+)_.*", "\\1", time_cols)))]
  result <- result[, c(col_order, sorted_time_cols)]
  
  # 验证
  cat("\n最终数据验证:\n")
  cat("唯一受试者数量: ", nrow(result), "\n")
  cat("列数: ", ncol(result), "\n")
  
  return(result)
}

# 计算总结统计数据的函数
calculate_daily_hr_summary <- function(data) {
  # 获取所有时间列（不包括subject_id和surgery_date）
  time_cols <- names(data)[-(1:2)]
  
  # 提取统计类型和天数
  summary_stats <- data.frame(
    column = time_cols,
    day = as.numeric(str_extract(time_cols, "-?\\d+")),
    stat_type = str_extract(time_cols, "(mean|min|max|median|sd|cv|iqr|skew|kurt|rmssd)")
  ) %>%
    mutate(valid_count = sapply(data[time_cols], function(x) sum(!is.na(x))))
  
  return(summary_stats)
}

# 计算每日非NA值计数的函数
calculate_daily_hr_counts <- function(data) {
  # 获取所有时间列（不包括subject_id和surgery_date）
  time_cols <- names(data)[-(1:2)]
  
  # 计算每列非NA值的数量
  counts <- sapply(data[time_cols], function(x) sum(!is.na(x)))
  
  # 创建结果数据框
  result <- data.frame(
    day = as.numeric(str_extract(names(counts), "-?\\d+")),
    stat_type = str_extract(names(counts), "(mean|min|max|median|sd|cv|iqr|skew|kurt|rmssd)"),
    valid_count = counts
  ) %>%
    arrange(day, stat_type)
  
  return(result)
}

# 执行计算
daily_hr_result <- calculate_daily_heart_rate(heart_rate_data, baseline_info)
summary_stats <- calculate_daily_hr_summary(daily_hr_result)
daily_counts <- calculate_daily_hr_counts(daily_hr_result)

# 保存结果
save(daily_hr_result, file = "daily_hr_result.rda", compress = "xz")

# 可选：创建一些基本的可视化
# 计算每日平均心率的可用数据量
mean_counts <- daily_counts %>%
  filter(stat_type == "mean") %>%
  arrange(day)

# 绘制每日平均心率的可用数据量
ggplot(mean_counts, aes(x = day, y = valid_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "每日平均心率的有效样本数",
    x = "手术相对天数",
    y = "有效样本数"
  )
ggsave("daily_hr_sample_counts.pdf", width = 10, height = 6)
