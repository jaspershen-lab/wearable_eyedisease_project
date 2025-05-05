library(tidyverse)
library(lme4)
library(caret)
library(ggplot2)
library(lmerTest)
library(corrplot)
library(patchwork)
library(r4projects)

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

# 输入输出目录设置
input_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data_grouped"
octa_dir <- "2_data/analysis_data" 
output_dir <- "3_data_analysis/5_presurgery_analysis/octa_change_prediction/t2_t0_selected_parameters"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义标准时间窗口
time_windows <- c(
  "pre_3d",        # 术前3天
  "pre_3d_7d",     # 术前3-7天
  "pre_7d_all",    # 术前7天及以上
  "post_7d",       # 术后7天内
  "post_7d_30d",   # 术后7-30天
  "post_day23_30", # 术后23-30天
  "post_day27_30", # 术后27-30天
  "post_1w",       # 术后1周
  "post_over_30d", # 术后超过30天
  "pre_all"        # 所有术前数据
)

#===============================================================
# 第1部分: 读取并预处理OCTA数据 - 仅选择指定的参数
#===============================================================

# 定义显著相关性的OCTA参数
selected_bloodflow_params <- c(
  "PA_Choriocapillaris",
  "SVD_NerveFiber",
  "PA_PED",
  "PA_Vitreous"
)

selected_thickness_params <- c(
  "Thickness_GCL.IPL",
  "Thickness_Retina",
  "Thickness_RNFL"
)

# 读取OCTA和基线数据
octa_bloodflow <- read.csv(file.path(octa_dir, "octa_data_bloodflow.csv"))
octa_thickness <- read.csv(file.path(octa_dir, "octa_data_thickness.csv"))
baseline_info <- read.csv(file.path(octa_dir, "baseline_info.csv"))

# 处理两个数据框的行数不匹配问题
mismatched_ids <- setdiff(octa_bloodflow$id, octa_thickness$id)
cat("在bloodflow但不在thickness中的ID:", mismatched_ids, "\n")

# 使用inner_join合并数据，确保只保留两个数据集中都有的ID
octa_combined <- octa_bloodflow %>%
  inner_join(octa_thickness, by = "id")
cat("合并后的数据行数:", nrow(octa_combined), "\n")

# 处理OCTA血流数据 - 只提取选定的参数
process_patient_bloodflow <- function(patient_data, param_patterns, time_points = c("T0", "T2")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术用左眼数据，右眼和双眼手术用右眼数据
  eye_code <- if(current_eye == 1) "OS" else "OD"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # 针对每个选定的参数创建列匹配模式
  cols_to_keep <- c()
  for(param in param_patterns) {
    for(tp in time_points) {
      # 修改模式匹配，适应实际列名格式
      pattern <- paste0(param, "_", eye_code, "_0_6_", tp, "$")
      matching_cols <- grep(pattern, names(patient_data), value = TRUE)
      cols_to_keep <- c(cols_to_keep, matching_cols)
    }
  }
  
  # 确保有匹配的列
  if(length(cols_to_keep) > 0) {
    # 选择数据并重命名列为更简单的格式
    param_data <- patient_data %>%
      dplyr::select("ID", all_of(cols_to_keep)) %>%
      rename_with(~ gsub(paste0("_", eye_code, "_0_6_"), "_", .), -ID)
    
    result <- result %>% left_join(param_data, by = "ID")
  }
  
  return(result)
}

# 处理OCTA厚度数据 - 只提取选定的参数（使用正确的版本）
process_patient_thickness <- function(patient_data, param_patterns, time_points = c("T0", "T2")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术用左眼数据，右眼和双眼手术用右眼数据
  eye_code <- if(current_eye == 1) "OS" else "OD"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # 针对每个选定的参数创建列匹配模式
  cols_to_keep <- c()
  for(param in param_patterns) {
    for(tp in time_points) {
      # 修改模式匹配，适应实际列名格式
      pattern <- paste0(param, "_", eye_code, "_0_6_", tp, "$")
      matching_cols <- grep(pattern, names(patient_data), value = TRUE)
      cols_to_keep <- c(cols_to_keep, matching_cols)
    }
  }
  
  # 确保有匹配的列
  if(length(cols_to_keep) > 0) {
    # 选择数据并重命名列为更简单的格式
    param_data <- patient_data %>%
      dplyr::select("ID", all_of(cols_to_keep)) %>%
      rename_with(~ gsub(paste0("_", eye_code, "_0_6_"), "_", .), -ID)
    
    result <- result %>% left_join(param_data, by = "ID")
  }
  
  return(result)
}

# 准备选定的OCTA血流数据
octa_bloodflow_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_bloodflow, by = c("ID" = "id"))

patient_list_bloodflow <- split(octa_bloodflow_features, octa_bloodflow_features$ID)
processed_data_bloodflow <- purrr::map(patient_list_bloodflow, function(data) {
  process_patient_bloodflow(data, selected_bloodflow_params)
})
octa_bloodflow_features <- bind_rows(processed_data_bloodflow)

# 准备选定的OCTA厚度数据（使用已更正的函数）
octa_thickness_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_thickness, by = c("ID" = "id"))

patient_list_thickness <- split(octa_thickness_features, octa_thickness_features$ID)
processed_data_thickness <- purrr::map(patient_list_thickness, function(data) {
  process_patient_thickness(data, selected_thickness_params)
})
octa_thickness_features <- bind_rows(processed_data_thickness)

# 打印提取的参数列，确认数据正确 - 修改grep模式来匹配新的列名格式
cat("提取的血流参数列：\n")
print(names(octa_bloodflow_features)[grepl("_T[0-2]$", names(octa_bloodflow_features))])

cat("\n提取的厚度参数列：\n")
print(names(octa_thickness_features)[grepl("_T[0-2]$", names(octa_thickness_features))])

# 如果上面的打印仍然为空，尝试查看所有列名
cat("\n所有血流特征列：\n")
print(names(octa_bloodflow_features))

cat("\n所有厚度特征列：\n")
print(names(octa_thickness_features))


#===============================================================
# 第2部分: 计算OCTA参数的T2-T0变化值
#===============================================================

# 修改 selected_params 以匹配处理后的列名
selected_params <- list(
  "PA_Choriocapillaris" = c(""),
  "SVD_NerveFiber" = c(""),
  "PA_PED" = c(""),
  "PA_Vitreous" = c(""),
  "Thickness_GCL.IPL" = c(""),
  "Thickness_Retina" = c(""),
  "Thickness_RNFL" = c("")
)

# 修改 calculate_changes 函数
calculate_changes <- function(bloodflow_data, thickness_data, id_col = "ID") {
  # 创建结果数据框
  result <- data.frame(ID = bloodflow_data[[id_col]])
  
  # 对血流参数计算变化
  for(param in c("PA_Choriocapillaris", "SVD_NerveFiber", "PA_PED", "PA_Vitreous")) {
    t0_col <- paste0(param, "_T0")
    t2_col <- paste0(param, "_T2")
    
    if(t0_col %in% names(bloodflow_data) && t2_col %in% names(bloodflow_data)) {
      change_col <- paste0(param, "_change")
      
      # 计算变化
      result[[change_col]] <- bloodflow_data[[t2_col]] - bloodflow_data[[t0_col]]
      
      # 检查数据质量
      n_valid <- sum(!is.na(result[[change_col]]))
      cat(sprintf("参数 %s 变化: 有效值 %d (%.1f%%)\n", 
                  param, n_valid, 100*n_valid/nrow(result)))
    } else {
      missing_cols <- setdiff(c(t0_col, t2_col), names(bloodflow_data))
      cat(sprintf("警告: 参数 %s 缺少列: %s\n", 
                  param, paste(missing_cols, collapse=", ")))
    }
  }
  
  # 对厚度参数计算变化
  for(param in c("Thickness_GCL.IPL", "Thickness_Retina", "Thickness_RNFL")) {
    t0_col <- paste0(param, "_T0")
    t2_col <- paste0(param, "_T2")
    
    if(t0_col %in% names(thickness_data) && t2_col %in% names(thickness_data)) {
      change_col <- paste0(param, "_change")
      
      # 计算变化
      result[[change_col]] <- thickness_data[[t2_col]] - thickness_data[[t0_col]]
      
      # 检查数据质量
      n_valid <- sum(!is.na(result[[change_col]]))
      cat(sprintf("参数 %s 变化: 有效值 %d (%.1f%%)\n", 
                  param, n_valid, 100*n_valid/nrow(result)))
    } else {
      missing_cols <- setdiff(c(t0_col, t2_col), names(thickness_data))
      cat(sprintf("警告: 参数 %s 缺少列: %s\n", 
                  param, paste(missing_cols, collapse=", ")))
    }
  }
  
  # 检查最终变化数据
  change_cols <- grep("_change$", names(result), value = TRUE)
  cat(sprintf("共计算 %d 个变化参数\n", length(change_cols)))
  
  if(length(change_cols) > 0) {
    # 添加变量到全局变量
    assign("octa_change_vars", change_cols, envir = .GlobalEnv)
    cat("变化参数变量已添加到全局环境: octa_change_vars\n")
  }
  
  return(result)
}

# 计算变化
all_changes <- calculate_changes(octa_bloodflow_features, octa_thickness_features)

# 保存变化数据
write.csv(all_changes, file.path(output_dir, "selected_octa_changes.csv"), row.names = FALSE)


# 为后续分析保存数据
bloodflow_changes <- all_changes
thickness_changes <- all_changes

# 保存变化参数数据
write.csv(bloodflow_changes, file.path(output_dir, "selected_bloodflow_changes.csv"), row.names = FALSE)
write.csv(thickness_changes, file.path(output_dir, "selected_thickness_changes.csv"), row.names = FALSE)

