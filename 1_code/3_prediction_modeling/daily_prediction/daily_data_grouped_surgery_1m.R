library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)
library(caret)
library(randomForest)
library(mice)

# 设置路径
input_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/daily_data"
output_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/daily_data_grouped"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 加载原始的baseline_info数据以获取糖尿病状态和手术类型
# 注意：可能需要调整路径以匹配您的实际文件位置
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# 从baseline_info中提取糖尿病状态和手术类型信息
group_info <- baseline_info %>%
  mutate(
    # 基于原始代码中的转换逻辑创建dm_status
    dm_status = case_when(
      diabetes_history == 1 ~ "Diabetes",
      diabetes_history == 2 ~ "No Diabetes",
      TRUE ~ NA_character_
    ),
    # 使用surgery_1..0.PI.1.other.作为手术类型
    surgery_type = surgery_1..0.PI.1.other.
  ) %>%
  dplyr::select(ID, dm_status, surgery_type)

# 显示分组信息的分布
cat("糖尿病状态和手术类型分布:\n")
print(table(group_info$dm_status, group_info$surgery_type, dnn = c("糖尿病状态", "手术类型")))

# 获取每日数据文件列表
daily_files <- list.files(input_dir, pattern = "^day_.*_data\\.csv$", full.names = TRUE)

# 创建一个数据框来存储每个文件的分组样本量
summary_df <- data.frame(
  day = character(),
  No_Diabetes_Surgery0_count = numeric(),
  Diabetes_Surgery1_count = numeric(),
  stringsAsFactors = FALSE
)

# 处理每个日期文件
for(file_path in daily_files) {
  # 提取日期信息
  day_str <- gsub(".*day_(.*)_data\\.csv", "\\1", basename(file_path))
  
  # 读取每日数据
  cat("\n处理文件:", basename(file_path), "\n")
  daily_data <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    cat("  错误: 无法读取文件:", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(daily_data)) next
  
  # 检查daily_data中的ID列名
  id_col <- NULL
  if("subject_id" %in% names(daily_data)) {
    id_col <- "subject_id"
  } else if("ID" %in% names(daily_data)) {
    id_col <- "ID"
  }
  
  if(is.null(id_col)) {
    cat("  错误: 文件中没有找到ID列(subject_id或ID)\n")
    next
  }
  
  # 将组信息与每日数据合并
  daily_grouped <- daily_data %>%
    left_join(group_info, by = setNames("ID", id_col))
  
  # 检查合并后的结果
  missing_group_info <- sum(is.na(daily_grouped$dm_status))
  if(missing_group_info > 0) {
    cat("  警告: ", missing_group_info, "行缺少分组信息\n")
  }
  
  # 根据dm_status和surgery_type分组
  # 1. No Diabetes Surgery 0组
  no_diabetes_surgery0 <- daily_grouped %>%
    filter(dm_status == "No Diabetes", surgery_type == 0)
  
  # 2. Diabetes Surgery 1组
  diabetes_surgery1 <- daily_grouped %>%
    filter(dm_status == "Diabetes", surgery_type == 1)
  
  # 输出每组样本量
  cat("  No Diabetes Surgery 0组: ", nrow(no_diabetes_surgery0), "人\n")
  cat("  Diabetes Surgery 1组: ", nrow(diabetes_surgery1), "人\n")
  
  # 添加到汇总数据框
  summary_df <- rbind(summary_df, data.frame(
    day = day_str,
    No_Diabetes_Surgery0_count = nrow(no_diabetes_surgery0),
    Diabetes_Surgery1_count = nrow(diabetes_surgery1)
  ))
  
  # 保存分组数据
  if(nrow(no_diabetes_surgery0) > 0) {
    out_file <- file.path(output_dir, paste0("day_", day_str, "_NoD_Surg0.csv"))
    write.csv(no_diabetes_surgery0, out_file, row.names = FALSE)
    cat("  已保存No Diabetes Surgery 0组数据到:", basename(out_file), "\n")
  }
  
  if(nrow(diabetes_surgery1) > 0) {
    out_file <- file.path(output_dir, paste0("day_", day_str, "_D_Surg1.csv"))
    write.csv(diabetes_surgery1, out_file, row.names = FALSE)
    cat("  已保存Diabetes Surgery 1组数据到:", basename(out_file), "\n")
  }
}

# 保存汇总信息
summary_file <- file.path(output_dir, "group_sample_sizes.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("\n样本量汇总信息已保存到:", summary_file, "\n")

# 输出完成信息
cat("\n处理完成! 已为", nrow(summary_df), "个日期文件生成分组数据\n")

# 可选：绘制样本量随时间变化的图表
if(requireNamespace("ggplot2", quietly = TRUE)) {
  # 转换day为数值，便于排序
  summary_df$day_num <- as.numeric(summary_df$day)
  
  # 按day_num排序
  summary_df <- summary_df[order(summary_df$day_num),]
  
  # 转换为长格式，便于绘图
  summary_long <- tidyr::pivot_longer(
    summary_df,
    cols = c("No_Diabetes_Surgery0_count", "Diabetes_Surgery1_count"),
    names_to = "group",
    values_to = "count"
  )
  
  # 美化组名
  summary_long$group <- factor(
    summary_long$group,
    levels = c("No_Diabetes_Surgery0_count", "Diabetes_Surgery1_count"),
    labels = c("No Diabetes Surgery 0", "Diabetes Surgery 1")
  )
  
  # 绘制样本量随时间变化的图表
  p <- ggplot2::ggplot(summary_long, ggplot2::aes(x = factor(day), y = count, fill = group)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(
      title = "各组在不同时间点的样本量",
      x = "相对于手术的天数",
      y = "样本量",
      fill = "组别"
    ) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # 保存图表
  ggplot2::ggsave(
    file.path(output_dir, "sample_sizes_by_day.pdf"),
    p,
    width = 10,
    height = 6
  )
  
}

###########################################
## 以下是合并数据表的附加代码
###########################################

# 设置合并数据的输出目录
combined_output_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/daily_data_grouped"

# 创建输出目录
dir.create(combined_output_dir, recursive = TRUE, showWarnings = FALSE)

# 获取所有分组后的数据文件
grouped_files <- list.files(output_dir, pattern = "^day_.*\\.(csv)$", full.names = TRUE)

# 排除汇总文件
grouped_files <- grouped_files[!grepl("group_sample_sizes\\.csv$", grouped_files)]

# 初始化空数据框用于存储合并的数据
combined_data <- data.frame()

# 处理每个分组文件
for(file_path in grouped_files) {
  # 提取文件名信息
  file_name <- basename(file_path)
  day_str <- gsub("day_(.*)_(NoD_Surg0|D_Surg1)\\.csv", "\\1", file_name)
  group_str <- gsub("day_.*_(NoD_Surg0|D_Surg1)\\.csv", "\\1", file_name)
  
  # 读取数据
  cat("处理文件:", file_name, "\n")
  group_data <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    cat("  错误: 无法读取文件:", e$message, "\n")
    return(NULL)
  })
  
  if(is.null(group_data)) next
  
  # 添加day列和group列以标识数据来源
  group_data$day <- day_str
  group_data$group <- ifelse(group_str == "NoD_Surg0", "No Diabetes Surgery 0", "Diabetes Surgery 1")
  
  # 添加到合并数据中
  combined_data <- bind_rows(combined_data, group_data)
}

# 检查合并后的数据
cat("\n合并后的数据维度:", dim(combined_data), "\n")

# 计算每个day和group组合的样本量
sample_summary <- combined_data %>%
  group_by(day, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = count, values_fill = 0)

cat("\n各组在不同时间点的样本量:\n")
print(sample_summary)

# 为day创建因子以正确排序
# 首先检查day列是否可以转换为数值
if(all(!is.na(suppressWarnings(as.numeric(combined_data$day))))) {
  combined_data$day_num <- as.numeric(combined_data$day)
  combined_data <- combined_data %>%
    arrange(day_num, group)
} else {
  # 如果不是纯数字，则使用原始字符串排序
  combined_data <- combined_data %>%
    arrange(day, group)
}

# 保存合并后的数据
combined_file <- file.path(combined_output_dir, "combined_all_days_groups.csv")
write.csv(combined_data, combined_file, row.names = FALSE)
cat("\n合并数据已保存到:", combined_file, "\n")

# 创建一个包含所有原始变量的数据集，但添加day和group标识
# 首先检查哪些列可能需要移除以避免重复
potential_duplicates <- c("subject_id", "ID", "day_num")
cols_to_remove <- intersect(names(combined_data), potential_duplicates)

combined_data_clean <- combined_data %>%
  # 检查是否有需要移除的列
  {if(length(cols_to_remove) > 0) dplyr::select(., -one_of(cols_to_remove)) else .}

# 保存清理后的合并数据
combined_clean_file <- file.path(combined_output_dir, "combined_clean.csv")
write.csv(combined_data_clean, combined_clean_file, row.names = FALSE)
cat("清理后的合并数据已保存到:", combined_clean_file, "\n")

# 额外添加：创建一个基本的统计摘要
# 首先识别数值型列
numeric_cols <- names(combined_data)[sapply(combined_data, is.numeric)]
# 排除可能的ID或索引列
exclude_cols <- c("ID", "subject_id", "day_num")
numeric_cols <- setdiff(numeric_cols, exclude_cols)

if(length(numeric_cols) > 0) {
  # 创建统计摘要
  stat_summary <- combined_data %>%
    group_by(day, group) %>%
    summarise(across(all_of(numeric_cols), 
                     list(mean = ~mean(., na.rm = TRUE),
                          sd = ~sd(., na.rm = TRUE),
                          median = ~median(., na.rm = TRUE),
                          n = ~sum(!is.na(.))),
                     .names = "{.col}_{.fn}"),
              .groups = "drop")
  
  # 保存统计摘要
  summary_file <- file.path(combined_output_dir, "combined_data_summary.csv")
  write.csv(stat_summary, summary_file, row.names = FALSE)
  cat("统计摘要已保存到:", basename(summary_file), "\n")
}

# 绘制合并数据的样本量图表
if(requireNamespace("ggplot2", quietly = TRUE)) {
  # 如果天数可以转换为数值，按数值排序
  if(all(!is.na(suppressWarnings(as.numeric(sample_summary$day))))) {
    sample_summary$day_num <- as.numeric(sample_summary$day)
    sample_summary <- sample_summary %>% arrange(day_num)
    day_levels <- sample_summary$day
  } else {
    day_levels <- sample_summary$day
  }
  
  # 转换为长格式
  sample_long <- sample_summary %>%
    pivot_longer(
      cols = c("No Diabetes Surgery 0", "Diabetes Surgery 1"),
      names_to = "group",
      values_to = "count"
    )
  
  # 设置day为因子以正确排序
  sample_long$day <- factor(sample_long$day, levels = day_levels)
  
  # 绘制图表
  p <- ggplot2::ggplot(sample_long, ggplot2::aes(x = day, y = count, fill = group)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(
      title = "合并数据中各组在不同时间点的样本量",
      x = "相对于手术的天数",
      y = "样本量",
      fill = "组别"
    ) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # 保存图表
  ggplot2::ggsave(
    file.path(combined_output_dir, "combined_sample_sizes_by_day.pdf"),
    p,
    width = 10,
    height = 6
  )
  cat("样本量图表已保存到: combined_sample_sizes_by_day.pdf\n")
}
