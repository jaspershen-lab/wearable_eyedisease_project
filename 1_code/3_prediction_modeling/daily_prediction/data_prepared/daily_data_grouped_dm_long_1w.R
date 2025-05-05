library(tidyverse)
library(lubridate)
setwd(get_project_wd())
rm(list = ls())

# 设置路径
input_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped/dm_group"
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/longitudinal_data/dm_group"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 获取NoD和D文件列表 (修改后的文件命名)
nod_files <- list.files(input_dir, pattern = "_NoD\\.csv$", full.names = TRUE)
d_files <- list.files(input_dir, pattern = "_D\\.csv$", full.names = TRUE)

# 为每组文件创建空列表，用于存储读取的数据框
nod_data_list <- list()
d_data_list <- list()

# 读取NoD文件并添加时间点信息
cat("处理No Diabetes组文件:\n")
for(file_path in nod_files) {
  # 提取日期信息
  day_str <- gsub(".*day_(.*)_NoD\\.csv", "\\1", basename(file_path))
  cat("  读取文件:", basename(file_path), "- 天数:", day_str, "\n")
  
  # 读取数据
  data <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    cat("  错误: 无法读取文件:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(data)) {
    # 添加time列以标识时间点
    data$time <- day_str
    
    # 标准化ID列名(如果需要)
    if("subject_id" %in% names(data) && !"ID" %in% names(data)) {
      data$ID <- data$subject_id
    }
    
    # 存储到列表中
    nod_data_list[[length(nod_data_list) + 1]] <- data
  }
}

# 读取D文件并添加时间点信息
cat("\n处理Diabetes组文件:\n")
for(file_path in d_files) {
  # 提取日期信息
  day_str <- gsub(".*day_(.*)_D\\.csv", "\\1", basename(file_path))
  cat("  读取文件:", basename(file_path), "- 天数:", day_str, "\n")
  
  # 读取数据
  data <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    cat("  错误: 无法读取文件:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(data)) {
    # 添加time列以标识时间点
    data$time <- day_str
    
    # 标准化ID列名(如果需要)
    if("subject_id" %in% names(data) && !"ID" %in% names(data)) {
      data$ID <- data$subject_id
    }
    
    # 存储到列表中
    d_data_list[[length(d_data_list) + 1]] <- data
  }
}

# 检查是否每个列表都有数据
if(length(nod_data_list) == 0) {
  cat("警告: 没有找到No Diabetes组数据!\n")
} else {
  # 合并NoD的所有数据框
  cat("\n合并No Diabetes组数据...\n")
  
  # 检查所有数据框的列是否一致
  column_lists <- lapply(nod_data_list, names)
  common_columns <- Reduce(intersect, column_lists)
  
  if(length(common_columns) == 0) {
    cat("错误: No Diabetes组数据框没有共同的列!\n")
  } else {
    # 保留共同的列，确保所有数据框有相同的结构
    nod_data_list <- lapply(nod_data_list, function(df) {
      df[, common_columns]
    })
    
    # 合并数据框
    nod_combined <- bind_rows(nod_data_list)
    
    # 确保时间点按正确顺序排序
    # 首先将time转换为数值型以正确排序
    nod_combined$time_num <- as.numeric(nod_combined$time)
    nod_combined <- nod_combined %>%
      arrange(ID, time_num) %>%
      dplyr::select(-time_num)  # 删除辅助列
    
    # 保存合并后的数据
    nod_output <- file.path(output_dir, "NoD_longitudinal.csv")
    write.csv(nod_combined, nod_output, row.names = FALSE)
    cat("  已保存合并后的No Diabetes组数据到:", basename(nod_output), "\n")
    cat("  共", nrow(nod_combined), "行,", length(unique(nod_combined$ID)), "个唯一ID\n")
  }
}

if(length(d_data_list) == 0) {
  cat("警告: 没有找到Diabetes组数据!\n")
} else {
  # 合并D的所有数据框
  cat("\n合并Diabetes组数据...\n")
  
  # 检查所有数据框的列是否一致
  column_lists <- lapply(d_data_list, names)
  common_columns <- Reduce(intersect, column_lists)
  
  if(length(common_columns) == 0) {
    cat("错误: Diabetes组数据框没有共同的列!\n")
  } else {
    # 保留共同的列，确保所有数据框有相同的结构
    d_data_list <- lapply(d_data_list, function(df) {
      df[, common_columns]
    })
    
    # 合并数据框
    d_combined <- bind_rows(d_data_list)
    
    # 确保时间点按正确顺序排序
    # 首先将time转换为数值型以正确排序
    d_combined$time_num <- as.numeric(d_combined$time)
    d_combined <- d_combined %>%
      arrange(ID, time_num) %>%
      dplyr::select(-time_num)  # 删除辅助列
    
    # 保存合并后的数据
    d_output <- file.path(output_dir, "D_longitudinal.csv")
    write.csv(d_combined, d_output, row.names = FALSE)
    cat("  已保存合并后的Diabetes组数据到:", basename(d_output), "\n")
    cat("  共", nrow(d_combined), "行,", length(unique(d_combined$ID)), "个唯一ID\n")
  }
}

# 输出汇总信息
cat("\n处理完成！\n")
cat("No Diabetes组数据:", length(nod_files), "个文件已合并\n")
cat("Diabetes组数据:", length(d_files), "个文件已合并\n")

library(lme4)

# 运行 LMM
model <- lmer(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + age + gender + bmi +time + pre_vision+season + (1 | ID), data = d_combined)

# 查看模型结果
summary(model)


# 选择特定时间点（例如术后第6天）的数据
day_prior4_data <- d_combined %>% 
  filter(time == "-4")

# 使用这一天的数据构建模型
model_day_prior4 <- lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
                         age + gender + bmi + pre_vision + season, 
                       data = day_prior4_data)

summary(model_day_prior4)
