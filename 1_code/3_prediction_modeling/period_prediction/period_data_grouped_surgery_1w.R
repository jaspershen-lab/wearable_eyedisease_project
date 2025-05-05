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
input_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/time_period_data"
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/time_period_data_grouped"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义标准时间窗口
time_windows <- c(
  "pre_3d",
  "pre_3d_7d", 
  "pre_7d_all", 
  "pre_all", 
  "post_1to3d",
  "post_4to6d",
  "post_6d",
  "post_7d",
  "post_7d_30d",
  "post_1w",
  "post_day23_30",
  "post_day27_30",
  "post_over_30d"
)

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
period_files <- list.files(input_dir, pattern = "^(pre|post)_.*_data\\.csv$", full.names = TRUE)

# 创建一个数据框来存储每个文件的分组样本量
summary_df <- data.frame(
  original_filename = character(),
  standardized_window = character(),
  No_Diabetes_Surgery0_count = numeric(),
  Diabetes_Surgery1_count = numeric(),
  stringsAsFactors = FALSE
)

# 处理每个日期文件
for(file_path in period_files) {
  # 提取日期信息
  file_name <- basename(file_path)
  original_window <- sub("^(pre|post)_(.*)_data\\.csv$", "\\1_\\2", file_name)
  
  # 初始化标准化窗口变量
  standardized_window <- original_window
  
  # 针对特定的时间窗口进行特殊处理
  if(original_window == "pre_3d_to_7d" || 
     original_window == "pre_surgery_3d_to_7d" || 
     grepl("pre_3d_7d", original_window)) {
    standardized_window <- "pre_3d_7d"
    cat("\n处理文件:", file_name, "\n")
    cat("  特殊处理: 3d到7d时间窗口\n")
    cat("  标准化窗口:", standardized_window, "\n")
  } else {
    # 对于带"surgery"的文件名
    if(grepl("surgery", original_window)) {
      standardized_window <- gsub("_surgery_", "_", original_window)
    }
    
    # 处理包含"_to_"的文件名
    if(grepl("_to_", standardized_window)) {
      standardized_window <- gsub("_to_", "_", standardized_window)
    }
    
    # 确保标准化窗口名称不含"_data"后缀
    standardized_window <- gsub("_data$", "", standardized_window)
    
    cat("\n处理文件:", file_name, "\n")
    cat("  原始窗口:", original_window, "\n")
    cat("  标准化窗口:", standardized_window, "\n")
  }
  
  # 读取每日数据
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
    original_filename = file_name,
    standardized_window = standardized_window,
    No_Diabetes_Surgery0_count = nrow(no_diabetes_surgery0),
    Diabetes_Surgery1_count = nrow(diabetes_surgery1)
  ))
  
  # 处理文件
  if(nrow(no_diabetes_surgery0) > 0) {
    # 确保pre_3d_7d正确输出
    if(standardized_window == "pre_3d_7d") {
      out_file <- file.path(output_dir, paste0("pre_3d_7d_NoD_Surg0.csv"))
      write.csv(no_diabetes_surgery0, out_file, row.names = FALSE)
      cat("  已保存pre_3d_7d无糖尿病组数据到:", basename(out_file), "\n")
    } else {
      out_file <- file.path(output_dir, paste0(standardized_window, "_NoD_Surg0.csv"))
      write.csv(no_diabetes_surgery0, out_file, row.names = FALSE)
      cat("  已保存No Diabetes Surgery 0组数据到:", basename(out_file), "\n")
    }
  }
  
  if(nrow(diabetes_surgery1) > 0) {
    # 确保pre_3d_7d正确输出
    if(standardized_window == "pre_3d_7d") {
      out_file <- file.path(output_dir, paste0("pre_3d_7d_D_Surg1.csv"))
      write.csv(diabetes_surgery1, out_file, row.names = FALSE)
      cat("  已保存pre_3d_7d糖尿病组数据到:", basename(out_file), "\n")
    } else {
      out_file <- file.path(output_dir, paste0(standardized_window, "_D_Surg1.csv"))
      write.csv(diabetes_surgery1, out_file, row.names = FALSE)
      cat("  已保存Diabetes Surgery 1组数据到:", basename(out_file), "\n")
    }
  }
}

# 保存汇总信息
summary_file <- file.path(output_dir, "group_sample_sizes.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("\n样本量汇总信息已保存到:", summary_file, "\n")

# 创建文件名到标准时间窗口的映射表
mapping_table <- summary_df %>%
  dplyr::select(original_filename, standardized_window) %>%
  distinct()

# 保存映射表
write.csv(mapping_table, 
          file.path(output_dir, "filename_to_window_mapping.csv"), 
          row.names = FALSE)

cat("\n已创建文件名到标准时间窗口的映射表:", 
    file.path(output_dir, "filename_to_window_mapping.csv"), "\n")

# 输出完成信息
cat("\n处理完成! 已为", nrow(mapping_table), "个时间窗口文件生成分组数据\n")

# 可选：绘制样本量随时间变化的图表
if(requireNamespace("ggplot2", quietly = TRUE)) {
  # 确保时间窗口是因子，以便按预定义顺序排序
  valid_windows <- intersect(time_windows, summary_df$standardized_window)
  
  if(length(valid_windows) > 0) {
    summary_df$standardized_window <- factor(
      summary_df$standardized_window,
      levels = time_windows,
      ordered = TRUE
    )
  }
  
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
  p <- ggplot2::ggplot(summary_long, 
                       ggplot2::aes(x = standardized_window, y = count, fill = group)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(
      title = "各组在不同时间窗口的样本量",
      x = "相对于手术的时间窗口",
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
    file.path(output_dir, "sample_sizes_by_period.pdf"),
    p,
    width = 10,
    height = 6
  )
}