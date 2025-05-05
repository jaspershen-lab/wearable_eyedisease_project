library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)
library(caret)
library(randomForest)
library(mice)

combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1w_prediction/daily_data/all_days_combined_data.csv")
# Load OCTA data - added from second code
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness.csv")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

dir.create("3_data_analysis/5_presurgery_analysis/basline_prediction/octa", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/basline_prediction/octa")


# Process OCTA blood flow data for baseline (T0)
octa_bloodflow_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_bloodflow, by = c("ID" = "id"))

# Function to process blood flow data for each patient at baseline (T0)
process_patient_bloodflow <- function(patient_data, time_points = c("T0")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # Use left eye data for left eye surgery, right eye data for right eye and both eyes surgery
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # Process each time point
  for(suffix in time_points) {
    # Select columns for current time point
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # Select data and rename columns
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# Process each patient individually for baseline (T0)
patient_list_bloodflow <- split(octa_bloodflow_features, octa_bloodflow_features$ID)
processed_data_bloodflow <- purrr::map(patient_list_bloodflow, process_patient_bloodflow)

# Combine results
octa_bloodflow_features <- bind_rows(processed_data_bloodflow)

# Create blood flow variables subset for baseline (T0)
bloodflow_var_T0 <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T0$"),  # Select all columns ending with 0_6_T0
    -matches("PA_OuterRetina_0_6_T0"),  # Exclude these columns
    -matches("PA_PED_0_6_T0")
  )

# Process OCTA thickness data for baseline (T0)
octa_thickness_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_thickness, by = c("ID" = "id"))

# Function to process thickness data for each patient at baseline (T0)
process_patient_thickness <- function(patient_data, time_points = c("T0")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # Use left eye data for left eye surgery, right eye data for right eye and both eyes surgery
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # Process each time point
  for(suffix in time_points) {
    # Select columns for current time point
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # Select data and rename columns
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# Process each patient individually for baseline (T0)
patient_list_thickness <- split(octa_thickness_features, octa_thickness_features$ID)
processed_data_thickness <- purrr::map(patient_list_thickness, process_patient_thickness)

# Combine results
octa_thickness_features <- bind_rows(processed_data_thickness)

# Create thickness variables subset for baseline (T0)
thickness_var_T0 <- octa_thickness_features %>%
  dplyr::select(
    ID,  # Keep ID column
    matches(".*_0_6_T0$"),
    matches("Thickness_PED_0_6_T0")
  )



# =====================================================
# 多元线性回归分析：可穿戴设备与OCTA变量
# =====================================================

# 加载必要的库
library(tidyverse)
library(broom)      # 用于整理回归结果
library(corrplot)   # 用于可视化相关矩阵
library(pheatmap)   # 用于高级热图绘制
library(RColorBrewer) # 用于配色方案
library(car)        # 用于VIF检测多重共线性
library(ggplot2)    # 用于绘图

# 设置随机种子保证结果可重复
set.seed(42)

# =====================================================
# 1. 准备数据
# =====================================================

# 定义变量组
wearable_vars <- c("mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep")
covariate_vars <- c("age", "gender", "bmi")
thickness_vars <- colnames(thickness_var_T0)[-1]  # 去掉ID列
bloodflow_vars <- colnames(bloodflow_var_T0)[-1]  # 去掉ID列

# 打印变量列表以确认
cat("可穿戴设备变量:", paste(wearable_vars, collapse=", "), "\n\n")
cat("协变量:", paste(covariate_vars, collapse=", "), "\n\n")
cat("厚度变量:", paste(thickness_vars, collapse=", "), "\n\n")
cat("血流变量:", paste(bloodflow_vars, collapse=", "), "\n\n")

# =====================================================
# 2. 多元线性回归函数
# =====================================================

# 函数：为单个OCTA变量进行多元线性回归
run_multiple_regression <- function(data, octa_var, wearable_vars, covariates) {
  # 准备数据
  model_data <- data %>% 
    dplyr::select(all_of(c(octa_var, wearable_vars, covariates))) %>%
    na.omit()
  
  # 如果数据不足，返回NULL
  if(nrow(model_data) < max(length(wearable_vars) + length(covariates) + 5, 15)) {
    warning(paste("变量", octa_var, "的有效样本量不足"))
    return(NULL)
  }
  
  # 构建回归公式
  formula_str <- paste(octa_var, "~", paste(c(wearable_vars, covariates), collapse = " + "))
  
  # 运行多元线性回归
  model <- lm(as.formula(formula_str), data = model_data)
  
  # 检查多重共线性
  vif_results <- try(car::vif(model), silent = TRUE)
  high_vif <- FALSE
  
  if(!inherits(vif_results, "try-error")) {
    high_vif <- any(vif_results > 5)
  }
  
  # 提取模型总结
  model_summary <- summary(model)
  
  # 使用broom整理结果
  tidy_results <- broom::tidy(model) %>%
    filter(term %in% wearable_vars)  # 只保留可穿戴设备变量的结果
  
  glance_results <- broom::glance(model)
  
  # 添加模型整体信息
  result <- tidy_results %>%
    mutate(
      octa_var = octa_var,
      sample_size = nrow(model_data),
      r_squared = glance_results$r.squared,
      adj_r_squared = glance_results$adj.r.squared,
      model_p_value = glance_results$p.value,
      high_vif = high_vif
    )
  
  # 返回结果和模型对象
  return(list(
    result = result,
    model = model,
    model_data = model_data,
    vif = if(!inherits(vif_results, "try-error")) vif_results else NULL
  ))
}

# 函数：为一组OCTA变量执行批量多元线性回归分析
perform_multiple_regression_analysis <- function(data, octa_vars, wearable_vars, covariate_vars, 
                                                 title = "多元线性回归分析", 
                                                 output_prefix = "multi_reg") {
  # 初始化结果数据框
  all_results <- data.frame()
  coefficients_matrix <- matrix(NA, nrow = length(wearable_vars), ncol = length(octa_vars),
                                dimnames = list(wearable_vars, octa_vars))
  p_matrix <- matrix(NA, nrow = length(wearable_vars), ncol = length(octa_vars),
                     dimnames = list(wearable_vars, octa_vars))
  
  # 存储所有模型对象
  all_models <- list()
  
  # 对每个OCTA变量进行多元回归
  for(octa_var in octa_vars) {
    # 运行模型
    model_result <- run_multiple_regression(data, octa_var, wearable_vars, covariate_vars)
    
    # 如果模型运行成功
    if(!is.null(model_result)) {
      # 添加结果到数据框
      all_results <- rbind(all_results, model_result$result)
      
      # 填充系数和p值矩阵
      for(i in 1:nrow(model_result$result)) {
        wearable_var <- gsub("term", "", model_result$result$term[i])
        wearable_idx <- which(wearable_vars == wearable_var)
        octa_idx <- which(octa_vars == octa_var)
        
        if(length(wearable_idx) > 0 && length(octa_idx) > 0) {
          coefficients_matrix[wearable_idx, octa_idx] <- model_result$result$estimate[i]
          p_matrix[wearable_idx, octa_idx] <- model_result$result$p.value[i]
        }
      }
      
      # 存储模型对象
      all_models[[octa_var]] <- model_result$model
    }
  }
  
  # 计算FDR校正的p值
  all_results$p_adj <- p.adjust(all_results$p.value, method = "BH")
  
  # 打印样本量信息
  cat("多元回归分析的样本量范围:", min(all_results$sample_size), "至", 
      max(all_results$sample_size), "\n")
  
  # 按p值排序
  all_results <- all_results %>%
    arrange(p.value)
  
  # 打印显著的结果
  cat("\n显著的多元回归结果 (p < 0.05):\n")
  sig_results <- all_results %>% filter(p.value < 0.05)
  if(nrow(sig_results) > 0) {
    print(sig_results %>% dplyr::select(octa_var, term, estimate, std.error, statistic, p.value))
  } else {
    cat("未发现显著关联\n")
  }
  
  # 打印FDR校正后显著的结果
  cat("\nFDR校正后显著的多元回归结果 (p_adj < 0.05):\n")
  sig_results_adj <- all_results %>% filter(p_adj < 0.05)
  if(nrow(sig_results_adj) > 0) {
    print(sig_results_adj %>% dplyr::select(octa_var, term, estimate, std.error, statistic, p.value, p_adj))
  } else {
    cat("FDR校正后未发现显著关联\n")
  }
  
  # 将结果保存为CSV
  write.csv(all_results, paste0(output_prefix, "_results.csv"), row.names = FALSE)
  write.csv(coefficients_matrix, paste0(output_prefix, "_coefficients.csv"))
  write.csv(p_matrix, paste0(output_prefix, "_pvalues.csv"))
  
  # 创建系数热图
  pdf(paste0(output_prefix, "_coef_heatmap.pdf"), width = 12, height = 10)
  # 复制矩阵用于绘图
  coefficients_plot <- coefficients_matrix
  coefficients_plot[is.na(coefficients_plot)] <- 0
  
  # 创建带有*标记显著性的系数矩阵
  coef_matrix_stars <- matrix("", nrow = nrow(coefficients_matrix), 
                              ncol = ncol(coefficients_matrix),
                              dimnames = dimnames(coefficients_matrix))
  
  for(i in 1:nrow(coefficients_matrix)) {
    for(j in 1:ncol(coefficients_matrix)) {
      if(is.na(coefficients_matrix[i,j]) || is.na(p_matrix[i,j])) {
        coef_matrix_stars[i,j] <- "NA"
      } else if(p_matrix[i,j] < 0.001) {
        coef_matrix_stars[i,j] <- paste0(round(coefficients_matrix[i,j], 3), "***")
      } else if(p_matrix[i,j] < 0.01) {
        coef_matrix_stars[i,j] <- paste0(round(coefficients_matrix[i,j], 3), "**")
      } else if(p_matrix[i,j] < 0.05) {
        coef_matrix_stars[i,j] <- paste0(round(coefficients_matrix[i,j], 3), "*")
      } else {
        coef_matrix_stars[i,j] <- as.character(round(coefficients_matrix[i,j], 3))
      }
    }
  }
  
  pheatmap(coefficients_plot,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           display_numbers = coef_matrix_stars,
           fontsize_number = 8,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = paste(title, "- 回归系数热图\n* p<0.05, ** p<0.01, *** p<0.001"))
  dev.off()
  
  # 返回结果
  return(list(
    all_results = all_results,
    coefficients_matrix = coefficients_matrix,
    p_matrix = p_matrix,
    models = all_models
  ))
}


# 步骤1: 准备术前7天的可穿戴数据均值
pre_surgery_data <- combined_data %>%
  filter(as.numeric(day) >= -7 & as.numeric(day) <= -1)  # 只选择术前7天数据

# 步骤2: 按ID分组计算术前均值
patient_preop_means <- pre_surgery_data %>%
  group_by(subject_id) %>%
  summarise(
    # 计算关键指标均值
    mean_hr = mean(mean_hr, na.rm = TRUE),           # 平均心率
    mean_rhr = mean(mean_rhr_1, na.rm = TRUE),       # 平均静息心率
    mean_bo = mean(mean_bo, na.rm = TRUE),           # 平均血氧
    total_steps = mean(steps_total, na.rm = TRUE),   # 平均总步数
    total_sleep = mean(total_sleep, na.rm = TRUE),   # 平均总睡眠
    
    # 保留其他人口统计学变量
    dm_2 = dplyr::first(dm_2),                  # 糖尿病状态
    age = dplyr::first(age),                    # 年龄
    gender = dplyr::first(gender),              # 性别
    bmi = dplyr::first(bmi)                     # 体重指数
  )

# 步骤3: 将可穿戴数据与OCTA厚度数据合并
wearable_thickness_data <- patient_preop_means %>%
  left_join(thickness_var_T0, by = c("subject_id" = "ID"))

# 步骤4: 将可穿戴数据与OCTA血流数据合并
wearable_bloodflow_data <- patient_preop_means %>%
  left_join(bloodflow_var_T0, by = c("subject_id" = "ID"))
# =====================================================
# 3. 执行多元线性回归分析 - 整体样本
# =====================================================

# 可穿戴设备与厚度变量 - 整体样本
cat("\n\n========= 多元线性回归分析: 可穿戴设备与厚度变量（整体样本）=========\n\n")
mr_thickness <- perform_multiple_regression_analysis(
  wearable_thickness_data,
  thickness_vars,
  wearable_vars,
  covariate_vars,
  "可穿戴设备与厚度变量",
  "mr_wearable_thickness_overall"
)

# 可穿戴设备与血流变量 - 整体样本
cat("\n\n========= 多元线性回归分析: 可穿戴设备与血流变量（整体样本）=========\n\n")
mr_bloodflow <- perform_multiple_regression_analysis(
  wearable_bloodflow_data,
  bloodflow_vars,
  wearable_vars,
  covariate_vars,
  "可穿戴设备与血流变量",
  "mr_wearable_bloodflow_overall"
)

# =====================================================
# 4. 执行多元线性回归分析 - 按糖尿病状态分层
# =====================================================

# 创建DM和非DM子集用于厚度分析
wearable_thickness_dm <- wearable_thickness_data %>% filter(dm_2 == 1)
wearable_thickness_nodm <- wearable_thickness_data %>% filter(dm_2 == 0)

# 创建DM和非DM子集用于血流分析
wearable_bloodflow_dm <- wearable_bloodflow_data %>% filter(dm_2 == 1)
wearable_bloodflow_nodm <- wearable_bloodflow_data %>% filter(dm_2 == 0)

# 打印分层分析的样本量
cat("\n分层分析的样本量：\n")
cat("厚度 - 糖尿病:", nrow(wearable_thickness_dm), "\n")
cat("厚度 - 非糖尿病:", nrow(wearable_thickness_nodm), "\n")
cat("血流 - 糖尿病:", nrow(wearable_bloodflow_dm), "\n")
cat("血流 - 非糖尿病:", nrow(wearable_bloodflow_nodm), "\n\n")

# 糖尿病组 - 可穿戴设备与厚度变量
cat("\n\n========= 多元线性回归分析: 可穿戴设备与厚度变量（糖尿病组）=========\n\n")
mr_thickness_dm <- perform_multiple_regression_analysis(
  wearable_thickness_dm,
  thickness_vars,
  wearable_vars,
  covariate_vars,
  "可穿戴设备与厚度变量（糖尿病组）",
  "mr_wearable_thickness_dm"
)

# 非糖尿病组 - 可穿戴设备与厚度变量
cat("\n\n========= 多元线性回归分析: 可穿戴设备与厚度变量（非糖尿病组）=========\n\n")
mr_thickness_nodm <- perform_multiple_regression_analysis(
  wearable_thickness_nodm,
  thickness_vars,
  wearable_vars,
  covariate_vars,
  "可穿戴设备与厚度变量（非糖尿病组）",
  "mr_wearable_thickness_nodm"
)

# 糖尿病组 - 可穿戴设备与血流变量
cat("\n\n========= 多元线性回归分析: 可穿戴设备与血流变量（糖尿病组）=========\n\n")
mr_bloodflow_dm <- perform_multiple_regression_analysis(
  wearable_bloodflow_dm,
  bloodflow_vars,
  wearable_vars,
  covariate_vars,
  "可穿戴设备与血流变量（糖尿病组）",
  "mr_wearable_bloodflow_dm"
)

# 非糖尿病组 - 可穿戴设备与血流变量
cat("\n\n========= 多元线性回归分析: 可穿戴设备与血流变量（非糖尿病组）=========\n\n")
mr_bloodflow_nodm <- perform_multiple_regression_analysis(
  wearable_bloodflow_nodm,
  bloodflow_vars,
  wearable_vars,
  covariate_vars,
  "可穿戴设备与血流变量（非糖尿病组）",
  "mr_wearable_bloodflow_nodm"
)

# =====================================================
# 5. 测试糖尿病状态的交互作用
# =====================================================

# 函数：为单个OCTA变量和可穿戴设备变量测试交互作用
test_diabetes_interaction_reg <- function(data, octa_var, wearable_var, covariates) {
  # 准备数据
  model_data <- data %>% 
    dplyr::select(all_of(c(octa_var, wearable_var, "dm_2", covariates))) %>%
    na.omit()
  
  # 如果数据不足，返回NULL
  if(nrow(model_data) < 15) {
    warning(paste("变量", octa_var, "和", wearable_var, "的有效样本量不足"))
    return(NULL)
  }
  
  # 构建带交互项的回归公式
  formula_str <- paste(octa_var, "~", wearable_var, "*dm_2 +", paste(covariates, collapse = " + "))
  
  # 运行回归模型
  model <- lm(as.formula(formula_str), data = model_data)
  
  # 使用broom整理结果
  tidy_results <- broom::tidy(model)
  
  # 查找交互项
  interaction_term <- paste0(wearable_var, ":dm_2")
  alt_interaction_term <- paste0("dm_2:", wearable_var)
  
  interaction_row <- tidy_results %>%
    filter(term == interaction_term | term == alt_interaction_term)
  
  if(nrow(interaction_row) == 0) {
    warning(paste("未找到交互项:", octa_var, "和", wearable_var))
    return(NULL)
  }
  
  # 提取交互效应结果
  result <- interaction_row %>%
    mutate(
      octa_var = octa_var,
      wearable_var = wearable_var,
      sample_size = nrow(model_data)
    )
  
  # 检查主效应
  main_effect_row <- tidy_results %>%
    filter(term == wearable_var)
  
  # 合并结果
  full_result <- data.frame(
    octa_var = octa_var,
    wearable_var = wearable_var,
    main_effect = if(nrow(main_effect_row) > 0) main_effect_row$estimate else NA,
    main_effect_p = if(nrow(main_effect_row) > 0) main_effect_row$p.value else NA,
    interaction_effect = result$estimate,
    interaction_se = result$std.error,
    interaction_t = result$statistic,
    interaction_p = result$p.value,
    sample_size = nrow(model_data),
    stringsAsFactors = FALSE
  )
  
  return(list(
    result = full_result,
    model = model
  ))
}

# 函数：批量测试糖尿病交互作用
test_all_diabetes_interactions_reg <- function(data, octa_vars, wearable_vars, covariate_vars,
                                               title = "糖尿病交互作用分析",
                                               output_prefix = "dm_int_reg") {
  # 初始化结果数据框
  all_interactions <- data.frame()
  interaction_matrix <- matrix(NA, nrow = length(wearable_vars), ncol = length(octa_vars),
                               dimnames = list(wearable_vars, octa_vars))
  p_matrix <- matrix(NA, nrow = length(wearable_vars), ncol = length(octa_vars),
                     dimnames = list(wearable_vars, octa_vars))
  
  # 对每个OCTA变量和可穿戴设备变量组合测试交互作用
  for(octa_var in octa_vars) {
    for(wearable_var in wearable_vars) {
      # 测试交互作用
      interaction_result <- test_diabetes_interaction_reg(data, octa_var, wearable_var, covariate_vars)
      
      # 如果模型运行成功
      if(!is.null(interaction_result)) {
        # 添加结果到数据框
        all_interactions <- rbind(all_interactions, interaction_result$result)
        
        # 填充系数和p值矩阵
        wearable_idx <- which(wearable_vars == wearable_var)
        octa_idx <- which(octa_vars == octa_var)
        interaction_matrix[wearable_idx, octa_idx] <- interaction_result$result$interaction_effect
        p_matrix[wearable_idx, octa_idx] <- interaction_result$result$interaction_p
      }
    }
  }
  
  # 计算FDR校正的p值
  all_interactions$p_adj <- p.adjust(all_interactions$interaction_p, method = "BH")
  
  # 按p值排序
  all_interactions <- all_interactions %>%
    arrange(interaction_p)
  
  # 打印显著的结果
  cat("\n显著的糖尿病交互作用 (p < 0.05):\n")
  sig_interactions <- all_interactions %>% filter(interaction_p < 0.05)
  if(nrow(sig_interactions) > 0) {
    print(sig_interactions %>% dplyr::select(octa_var, wearable_var, main_effect, main_effect_p, 
                                             interaction_effect, interaction_p))
  } else {
    cat("未发现显著交互作用\n")
  }
  
  # 打印FDR校正后显著的结果
  cat("\nFDR校正后显著的糖尿病交互作用 (p_adj < 0.05):\n")
  sig_interactions_adj <- all_interactions %>% filter(p_adj < 0.05)
  if(nrow(sig_interactions_adj) > 0) {
    print(sig_interactions_adj %>% dplyr::select(octa_var, wearable_var, main_effect, main_effect_p, 
                                                 interaction_effect, interaction_p, p_adj))
  } else {
    cat("FDR校正后未发现显著交互作用\n")
  }
  
  # 将结果保存为CSV
  write.csv(all_interactions, paste0(output_prefix, "_results.csv"), row.names = FALSE)
  write.csv(interaction_matrix, paste0(output_prefix, "_coefficients.csv"))
  write.csv(p_matrix, paste0(output_prefix, "_pvalues.csv"))
  
  # 创建交互系数热图
  pdf(paste0(output_prefix, "_interaction_heatmap.pdf"), width = 12, height = 10)
  # 复制矩阵用于绘图
  interaction_plot <- interaction_matrix
  interaction_plot[is.na(interaction_plot)] <- 0
  
  # 创建带有*标记显著性的系数矩阵
  interaction_matrix_stars <- matrix("", nrow = nrow(interaction_matrix), 
                                     ncol = ncol(interaction_matrix),
                                     dimnames = dimnames(interaction_matrix))
  
  for(i in 1:nrow(interaction_matrix)) {
    for(j in 1:ncol(interaction_matrix)) {
      if(is.na(interaction_matrix[i,j]) || is.na(p_matrix[i,j])) {
        interaction_matrix_stars[i,j] <- "NA"
      } else if(p_matrix[i,j] < 0.001) {
        interaction_matrix_stars[i,j] <- paste0(round(interaction_matrix[i,j], 3), "***")
      } else if(p_matrix[i,j] < 0.01) {
        interaction_matrix_stars[i,j] <- paste0(round(interaction_matrix[i,j], 3), "**")
      } else if(p_matrix[i,j] < 0.05) {
        interaction_matrix_stars[i,j] <- paste0(round(interaction_matrix[i,j], 3), "*")
      } else {
        interaction_matrix_stars[i,j] <- as.character(round(interaction_matrix[i,j], 3))
      }
    }
  }
  
  pheatmap(interaction_plot,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           display_numbers = interaction_matrix_stars,
           fontsize_number = 8,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = paste(title, "- 糖尿病交互作用\n* p<0.05, ** p<0.01, *** p<0.001"))
  dev.off()
  
  # 返回结果
  return(list(
    all_interactions = all_interactions,
    interaction_matrix = interaction_matrix,
    p_matrix = p_matrix
  ))
}

# 测试厚度变量的糖尿病交互作用
cat("\n\n========= 糖尿病交互作用分析: 可穿戴设备与厚度变量 =========\n\n")
dm_int_reg_thickness <- test_all_diabetes_interactions_reg(
  wearable_thickness_data,
  thickness_vars,
  wearable_vars,
  covariate_vars,
  "可穿戴设备与厚度变量",
  "dm_int_reg_wearable_thickness"
)

# 测试血流变量的糖尿病交互作用
cat("\n\n========= 糖尿病交互作用分析: 可穿戴设备与血流变量 =========\n\n")
dm_int_reg_bloodflow <- test_all_diabetes_interactions_reg(
  wearable_bloodflow_data,
  bloodflow_vars,
  wearable_vars,
  covariate_vars,
  "可穿戴设备与血流变量",
  "dm_int_reg_wearable_bloodflow"
)

# =====================================================
# 6. 汇总所有显著的结果
# =====================================================

# 创建汇总显著结果的函数
summarize_significant_regression <- function(result_list) {
  # 初始化空数据框
  all_sig <- data.frame()
  
  # 遍历结果列表
  for(name in names(result_list)) {
    if(!is.null(result_list[[name]])) {
      # 提取显著结果
      sig_results <- result_list[[name]]$all_results %>%
        filter(p.value < 0.05) %>%
        mutate(analysis = name)
      
      # 添加到汇总数据框
      all_sig <- rbind(all_sig, sig_results)
    }
  }
  
  # 如果有显著结果，按p值排序
  if(nrow(all_sig) > 0) {
    all_sig <- all_sig %>%
      arrange(p.value)
  }
  
  return(all_sig)
}

# 创建结果列表
all_mr_results <- list(
  "整体_厚度" = mr_thickness,
  "整体_血流" = mr_bloodflow,
  "糖尿病_厚度" = mr_thickness_dm,
  "糖尿病_血流" = mr_bloodflow_dm,
  "非糖尿病_厚度" = mr_thickness_nodm,
  "非糖尿病_血流" = mr_bloodflow_nodm
)

# 汇总所有显著结果
all_significant_mr <- summarize_significant_regression(all_mr_results)

# 打印并保存所有显著结果
cat("\n\n========= 所有显著的多元回归结果汇总 (p < 0.05) =========\n\n")
if(nrow(all_significant_mr) > 0) {
  print(all_significant_mr %>% dplyr::select(octa_var, term, estimate, p.value, analysis))
  write.csv(all_significant_mr, "all_significant_multiple_regression.csv", row.names = FALSE)
} else {
  cat("所有分析中均未发现显著关联\n")
}

# =====================================================
# 7. 可视化显著结果
# =====================================================

# 函数：可视化显著关联
plot_significant_regression <- function(data, octa_var, wearable_var, covariates = NULL,
                                        by_dm = FALSE, title = NULL) {
  # 准备数据
  plot_data <- data %>% 
    dplyr::select(all_of(c(octa_var, wearable_var, "dm_2", covariates))) %>%
    na.omit() %>%
    mutate(dm_status = factor(dm_2, levels = c(0, 1), 
                              labels = c("非糖尿病", "糖尿病")))
  
  # 创建散点图
  p <- ggplot(plot_data, aes_string(x = wearable_var, y = octa_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_minimal() +
    labs(
      title = ifelse(is.null(title), 
                     paste(octa_var, "与", wearable_var, "的关系"),
                     title),
      x = wearable_var,
      y = octa_var
    )
  
  # 如果需要按糖尿病状态分组
  if(by_dm) {
    p <- ggplot(plot_data, aes_string(x = wearable_var, y = octa_var, color = "dm_status")) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE) +
      theme_minimal() +
      scale_color_brewer(palette = "Set1", name = "糖尿病状态") +
      labs(
        title = ifelse(is.null(title), 
                       paste(octa_var, "与", wearable_var, "的关系 (按糖尿病状态)"),
                       title),
        x = wearable_var,
        y = octa_var
      )
  }
  
  return(p)
}

# 创建结果文件夹
dir.create("plots", showWarnings = FALSE)

# 可视化总体分析中的显著结果
if(nrow(all_significant_mr) > 0) {
  # 限制为top 10显著结果
  top_results <- head(all_significant_mr, 10)
  
  for(i in 1:nrow(top_results)) {
    result <- top_results[i, ]
    
    # 确定使用哪个数据集
    if(grepl("厚度", result$analysis)) {
      plot_data <- wearable_thickness_data
    } else {
      plot_data <- wearable_bloodflow_data
    }
    
    # 确定是否按糖尿病状态分组
    by_dm <- grepl("糖尿病", result$analysis) | grepl("非糖尿病", result$analysis)
    
    # 创建图
    p <- plot_significant_regression(
      plot_data, 
      result$octa_var, 
      gsub("term", "", result$term), 
      covariate_vars,
      by_dm = TRUE,  # 总是显示糖尿病分组以便进行比较
      title = paste0(
        result$octa_var, " 与 ", result$term,
        "\np = ", format(result$p.value, digits = 3),
        ", 分析: ", result$analysis
      )
    )
    
    # 保存图
    filename <- paste0("plots/", result$octa_var, "_", result$term, "_", 
                       gsub(" ", "_", result$analysis), ".pdf")
    ggsave(filename, p, width = 8, height = 6)
  }
  
  cat("\n保存了", min(nrow(top_results), 10), "个显著关联的可视化结果到plots文件夹\n")
} else {
  cat("\n没有显著结果可以可视化\n")
}


# =====================================================
# 8. OCTA变量的模型性能评估：糖尿病组与非糖尿病组比较
# =====================================================

# 创建模型性能评估函数
evaluate_model_performance <- function(models_list, dataset, octa_vars, group_name) {
  # 初始化结果数据框
  performance_metrics <- data.frame(
    OCTA_Variable = character(),
    R_Squared = numeric(),
    Adjusted_R_Squared = numeric(),
    RMSE = numeric(),
    MAE = numeric(),
    AIC = numeric(),
    BIC = numeric(),
    Significant_Predictors = character(),
    stringsAsFactors = FALSE
  )
  
  # 遍历所有OCTA变量模型
  for(octa_var in octa_vars) {
    # 获取模型
    model <- models_list[[octa_var]]
    
    # 如果模型存在
    if(!is.null(model)) {
      # 提取模型摘要
      model_summary <- summary(model)
      
      # 计算性能指标
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      
      # 计算残差指标
      actual_values <- dataset[[octa_var]]
      predicted_values <- predict(model, newdata = dataset)
      rmse <- sqrt(mean((actual_values - predicted_values)^2, na.rm = TRUE))
      mae <- mean(abs(actual_values - predicted_values), na.rm = TRUE)
      
      # 计算信息准则
      aic <- AIC(model)
      bic <- BIC(model)
      
      # 找出显著的预测变量 (p < 0.05)
      coef_table <- coef(model_summary)
      if(is.matrix(coef_table)) {
        significant_vars <- rownames(coef_table)[coef_table[, "Pr(>|t|)"] < 0.05 & 
                                                   rownames(coef_table) != "(Intercept)"]
        significant_vars_str <- paste(significant_vars, collapse = ", ")
      } else {
        significant_vars_str <- "No significant predictors"
      }
      
      # 添加到结果数据框
      performance_metrics <- rbind(performance_metrics, data.frame(
        OCTA_Variable = octa_var,
        R_Squared = r_squared,
        Adjusted_R_Squared = adj_r_squared,
        RMSE = rmse,
        MAE = mae,
        AIC = aic,
        BIC = bic,
        Significant_Predictors = significant_vars_str,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 按R²降序排序
  performance_metrics <- performance_metrics %>%
    arrange(desc(R_Squared))
  
  # 保存结果
  write.csv(performance_metrics, 
            paste0("octa_model_performance_", tolower(gsub(" ", "_", group_name)), ".csv"), 
            row.names = FALSE)
  
  # 返回结果
  return(performance_metrics)
}

# 创建模型性能可视化函数
plot_model_performance <- function(dm_performance, nodm_performance, metric = "R_Squared") {
  # 合并两组性能数据
  combined_performance <- rbind(
    dm_performance %>% mutate(Group = "Diabetic"),
    nodm_performance %>% mutate(Group = "Non-diabetic")
  )
  
  # 确定要使用的度量标准和标签
  metric_label <- case_when(
    metric == "R_Squared" ~ "R²",
    metric == "Adjusted_R_Squared" ~ "Adjusted R²",
    metric == "RMSE" ~ "RMSE",
    metric == "MAE" ~ "MAE",
    metric == "AIC" ~ "AIC",
    metric == "BIC" ~ "BIC",
    TRUE ~ metric
  )
  
  # 为OCTA变量创建更具可读性的标签
  combined_performance <- combined_performance %>%
    mutate(
      OCTA_Label = case_when(
        # 根据实际OCTA变量名添加适当的映射
        grepl("Thickness", OCTA_Variable) ~ paste("Thickness:", gsub("Thickness_", "", OCTA_Variable)),
        grepl("VD", OCTA_Variable) ~ paste("Vessel Density:", gsub("VD_", "", OCTA_Variable)),
        grepl("PA", OCTA_Variable) ~ paste("Perfusion Area:", gsub("PA_", "", OCTA_Variable)),
        grepl("SVD", OCTA_Variable) ~ paste("SVD:", gsub("SVD_", "", OCTA_Variable)),
        TRUE ~ OCTA_Variable
      )
    )
  
  # 创建条形图比较模型性能
  plot <- ggplot(combined_performance, aes(x = reorder(OCTA_Label, !!sym(metric)), 
                                           y = !!sym(metric), 
                                           fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Diabetic" = "#D6604D", "Non-diabetic" = "#4393C3")) +
    labs(
      title = paste("Model Performance Comparison:", metric_label),
      subtitle = "Diabetic vs Non-diabetic Groups",
      x = "OCTA Variables",
      y = metric_label,
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # 保存图表
  ggsave(paste0("octa_model_performance_", metric, "_comparison.pdf"), 
         plot, width = 12, height = 8)
  
  return(plot)
}

# 创建预测值与实际值对比图函数
create_prediction_plot <- function(model, data, octa_var, group_name) {
  # 计算预测值
  data$predicted <- predict(model, newdata = data)
  
  # 创建散点图
  plot <- ggplot(data, aes(x = predicted, y = .data[[octa_var]])) +
    geom_point(alpha = 0.7, color = ifelse(group_name == "Diabetic", "#D6604D", "#4393C3")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = paste("Predicted vs Actual", octa_var),
      subtitle = paste(group_name, "Group"),
      x = "Predicted Values",
      y = "Actual Values"
    ) +
    annotate("text", x = min(data$predicted, na.rm = TRUE), 
             y = max(data[[octa_var]], na.rm = TRUE),
             label = sprintf("R² = %.3f\nRMSE = %.3f", 
                             summary(model)$r.squared,
                             sqrt(mean(resid(model)^2))),
             hjust = 0, vjust = 1) +
    theme_minimal()
  
  # 保存图表
  ggsave(paste0("octa_", octa_var, "_", tolower(gsub(" ", "_", group_name)), "_prediction.pdf"), 
         plot, width = 8, height = 6)
  
  return(plot)
}

# 执行性能评估 - 糖尿病组
cat("\n\n========= OCTA模型性能评估：糖尿病组 =========\n\n")
dm_performance <- evaluate_model_performance(
  mr_thickness_dm$models,  # 假设您的模型存储在这个结构中
  wearable_thickness_dm,   # 糖尿病组厚度数据
  thickness_vars,          # 厚度变量列表
  "Diabetic Thickness"     # 组名
)

# 输出糖尿病组模型性能摘要
print(dm_performance %>% dplyr::select(OCTA_Variable, R_Squared, Adjusted_R_Squared, RMSE, MAE))

# 执行性能评估 - 非糖尿病组
cat("\n\n========= OCTA模型性能评估：非糖尿病组 =========\n\n")
nodm_performance <- evaluate_model_performance(
  mr_thickness_nodm$models,  # 非糖尿病组模型
  wearable_thickness_nodm,   # 非糖尿病组厚度数据
  thickness_vars,            # 厚度变量列表
  "Non-diabetic Thickness"   # 组名
)

# 输出非糖尿病组模型性能摘要
print(nodm_performance %>% dplyr::select(OCTA_Variable, R_Squared, Adjusted_R_Squared, RMSE, MAE))

# 对血流变量重复同样的操作
cat("\n\n========= OCTA模型性能评估：糖尿病组(血流) =========\n\n")
dm_bloodflow_performance <- evaluate_model_performance(
  mr_bloodflow_dm$models,  # 糖尿病组血流模型
  wearable_bloodflow_dm,   # 糖尿病组血流数据
  bloodflow_vars,          # 血流变量列表
  "Diabetic Bloodflow"     # 组名
)

# 输出糖尿病组血流模型性能摘要
print(dm_bloodflow_performance %>% dplyr::select(OCTA_Variable, R_Squared, Adjusted_R_Squared, RMSE, MAE))

cat("\n\n========= OCTA模型性能评估：非糖尿病组(血流) =========\n\n")
nodm_bloodflow_performance <- evaluate_model_performance(
  mr_bloodflow_nodm$models,  # 非糖尿病组血流模型
  wearable_bloodflow_nodm,   # 非糖尿病组血流数据
  bloodflow_vars,            # 血流变量列表
  "Non-diabetic Bloodflow"   # 组名
)

# 输出非糖尿病组血流模型性能摘要
print(nodm_bloodflow_performance %>% dplyr::select(OCTA_Variable, R_Squared, Adjusted_R_Squared, RMSE, MAE))

# 创建性能比较可视化
# R-squared比较
r2_comparison_thickness <- plot_model_performance(dm_performance, nodm_performance, "R_Squared")
print(r2_comparison_thickness)

# RMSE比较
rmse_comparison_thickness <- plot_model_performance(dm_performance, nodm_performance, "RMSE")
print(rmse_comparison_thickness)

# 血流变量的R-squared比较
r2_comparison_bloodflow <- plot_model_performance(dm_bloodflow_performance, nodm_bloodflow_performance, "R_Squared")
print(r2_comparison_bloodflow)

# 血流变量的RMSE比较
rmse_comparison_bloodflow <- plot_model_performance(dm_bloodflow_performance, nodm_bloodflow_performance, "RMSE")
print(rmse_comparison_bloodflow)

# 为表现最好的模型创建预测图
# 厚度变量
best_thickness_dm_var <- dm_performance$OCTA_Variable[1]
best_thickness_nodm_var <- nodm_performance$OCTA_Variable[1]

best_prediction_plot_dm_thickness <- create_prediction_plot(
  mr_thickness_dm$models[[best_thickness_dm_var]],
  wearable_thickness_dm,
  best_thickness_dm_var,
  "Diabetic"
)
print(best_prediction_plot_dm_thickness)

best_prediction_plot_nodm_thickness <- create_prediction_plot(
  mr_thickness_nodm$models[[best_thickness_nodm_var]],
  wearable_thickness_nodm,
  best_thickness_nodm_var,
  "Non-diabetic"
)
print(best_prediction_plot_nodm_thickness)

# 血流变量
best_bloodflow_dm_var <- dm_bloodflow_performance$OCTA_Variable[1]
best_bloodflow_nodm_var <- nodm_bloodflow_performance$OCTA_Variable[1]

best_prediction_plot_dm_bloodflow <- create_prediction_plot(
  mr_bloodflow_dm$models[[best_bloodflow_dm_var]],
  wearable_bloodflow_dm,
  best_bloodflow_dm_var,
  "Diabetic"
)
print(best_prediction_plot_dm_bloodflow)

best_prediction_plot_nodm_bloodflow <- create_prediction_plot(
  mr_bloodflow_nodm$models[[best_bloodflow_nodm_var]],
  wearable_bloodflow_nodm,
  best_bloodflow_nodm_var,
  "Non-diabetic"
)
print(best_prediction_plot_nodm_bloodflow)

# 模型比较表格
# 合并所有性能数据
all_performance <- bind_rows(
  dm_performance %>% mutate(Category = "Thickness", Group = "Diabetic"),
  nodm_performance %>% mutate(Category = "Thickness", Group = "Non-diabetic"),
  dm_bloodflow_performance %>% mutate(Category = "Bloodflow", Group = "Diabetic"),
  nodm_bloodflow_performance %>% mutate(Category = "Bloodflow", Group = "Non-diabetic")
)

# 找出每个类别/组合最好的模型
best_models <- all_performance %>%
  group_by(Category, Group) %>%
  top_n(1, R_Squared) %>%
  ungroup() %>%
  dplyr::select(OCTA_Variable, Category, Group, R_Squared, Adjusted_R_Squared, RMSE, 
                MAE, Significant_Predictors)

# 保存最佳模型表格
write.csv(best_models, "octa_best_models_summary.csv", row.names = FALSE)

# 创建最佳模型性能比较表格
best_models_table <- best_models %>%
  mutate(
    Group_Category = paste(Group, Category, sep = " - "),
    Model_Details = paste0(
      "R² = ", round(R_Squared, 3),
      ", RMSE = ", round(RMSE, 3),
      "\nSignificant: ", Significant_Predictors
    )
  ) %>%
  dplyr::select(OCTA_Variable, Group_Category, Model_Details)

# 打印最佳模型表格
print(best_models_table)

# 额外：检验两组模型性能差异的显著性
# 这里我们计算每组模型的平均R²
group_performance_summary <- all_performance %>%
  group_by(Category, Group) %>%
  summarise(
    Mean_R2 = mean(R_Squared, na.rm = TRUE),
    Median_R2 = median(R_Squared, na.rm = TRUE),
    SD_R2 = sd(R_Squared, na.rm = TRUE),
    Mean_RMSE = mean(RMSE, na.rm = TRUE),
    Median_RMSE = median(RMSE, na.rm = TRUE),
    SD_RMSE = sd(RMSE, na.rm = TRUE),
    n_models = n(),
    .groups = "drop"
  )

# 打印组间性能摘要
print(group_performance_summary)

# 将结果写入文件
write.csv(group_performance_summary, "octa_group_performance_summary.csv", row.names = FALSE)

# 创建比较表的函数
create_comparison_table <- function(var_name, dm_model, nodm_model) {
  # 提取每个模型的系数
  dm_coefs <- coef(summary(dm_model))
  nodm_coefs <- coef(summary(nodm_model))
  
  # 创建包含两个模型的系数表
  coefs_df <- data.frame(
    Variable = rownames(dm_coefs),
    DM_Estimate = dm_coefs[, "Estimate"],
    DM_P_Value = dm_coefs[, "Pr(>|t|)"],
    NonDM_Estimate = nodm_coefs[, "Estimate"],
    NonDM_P_Value = nodm_coefs[, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ) %>%
    filter(Variable != "(Intercept)")
  
  # 添加显著性标记
  coefs_df <- coefs_df %>%
    mutate(
      DM_Significant = ifelse(DM_P_Value < 0.05, "*", ""),
      NonDM_Significant = ifelse(NonDM_P_Value < 0.05, "*", "")
    )
  
  # 保存表格
  write.csv(coefs_df, paste0("octa_", var_name, "_coefficient_comparison.csv"), row.names = FALSE)
  
  return(coefs_df)
}

# 对最佳模型创建系数比较表
best_thickness_comparison <- create_comparison_table(
  best_thickness_dm_var,
  mr_thickness_dm$models[[best_thickness_dm_var]],
  mr_thickness_nodm$models[[best_thickness_nodm_var]]
)

best_bloodflow_comparison <- create_comparison_table(
  best_bloodflow_dm_var,
  mr_bloodflow_dm$models[[best_bloodflow_dm_var]],
  mr_bloodflow_nodm$models[[best_bloodflow_nodm_var]]
)

# 打印比较表
print(best_thickness_comparison)
print(best_bloodflow_comparison)


# =====================================================
# 8. OCTA变量的模型性能评估：糖尿病组与非糖尿病组比较
# =====================================================

# 创建模型性能评估函数
evaluate_model_performance <- function(models_list, dataset, octa_vars, group_name) {
  # 初始化结果数据框
  performance_metrics <- data.frame(
    OCTA_Variable = character(),
    R_Squared = numeric(),
    Adjusted_R_Squared = numeric(),
    RMSE = numeric(),
    MAE = numeric(),
    AIC = numeric(),
    BIC = numeric(),
    Significant_Predictors = character(),
    stringsAsFactors = FALSE
  )
  
  # 遍历所有OCTA变量模型
  for(octa_var in octa_vars) {
    # 获取模型
    model <- models_list[[octa_var]]
    
    # 如果模型存在
    if(!is.null(model)) {
      # 提取模型摘要
      model_summary <- summary(model)
      
      # 计算性能指标
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      
      # 计算残差指标
      actual_values <- dataset[[octa_var]]
      predicted_values <- predict(model, newdata = dataset)
      rmse <- sqrt(mean((actual_values - predicted_values)^2, na.rm = TRUE))
      mae <- mean(abs(actual_values - predicted_values), na.rm = TRUE)
      
      # 计算信息准则
      aic <- AIC(model)
      bic <- BIC(model)
      
      # 找出显著的预测变量 (p < 0.05)
      coef_table <- coef(model_summary)
      if(is.matrix(coef_table)) {
        significant_vars <- rownames(coef_table)[coef_table[, "Pr(>|t|)"] < 0.05 & 
                                                   rownames(coef_table) != "(Intercept)"]
        significant_vars_str <- paste(significant_vars, collapse = ", ")
      } else {
        significant_vars_str <- "No significant predictors"
      }
      
      # 添加到结果数据框
      performance_metrics <- rbind(performance_metrics, data.frame(
        OCTA_Variable = octa_var,
        R_Squared = r_squared,
        Adjusted_R_Squared = adj_r_squared,
        RMSE = rmse,
        MAE = mae,
        AIC = aic,
        BIC = bic,
        Significant_Predictors = significant_vars_str,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 按R²降序排序
  performance_metrics <- performance_metrics %>%
    arrange(desc(R_Squared))
  
  # 保存结果
  write.csv(performance_metrics, 
            paste0("octa_model_performance_", tolower(gsub(" ", "_", group_name)), ".csv"), 
            row.names = FALSE)
  
  # 返回结果
  return(performance_metrics)
}

# 创建模型性能可视化函数
plot_model_performance <- function(dm_performance, nodm_performance, metric = "R_Squared") {
  # 合并两组性能数据
  combined_performance <- rbind(
    dm_performance %>% mutate(Group = "Diabetic"),
    nodm_performance %>% mutate(Group = "Non-diabetic")
  )
  
  # 确定要使用的度量标准和标签
  metric_label <- case_when(
    metric == "R_Squared" ~ "R²",
    metric == "Adjusted_R_Squared" ~ "Adjusted R²",
    metric == "RMSE" ~ "RMSE",
    metric == "MAE" ~ "MAE",
    metric == "AIC" ~ "AIC",
    metric == "BIC" ~ "BIC",
    TRUE ~ metric
  )
  
  # 为OCTA变量创建更具可读性的标签
  combined_performance <- combined_performance %>%
    mutate(
      OCTA_Label = case_when(
        # 根据实际OCTA变量名添加适当的映射
        grepl("Thickness", OCTA_Variable) ~ paste("Thickness:", gsub("Thickness_", "", OCTA_Variable)),
        grepl("VD", OCTA_Variable) ~ paste("Vessel Density:", gsub("VD_", "", OCTA_Variable)),
        grepl("PA", OCTA_Variable) ~ paste("Perfusion Area:", gsub("PA_", "", OCTA_Variable)),
        grepl("SVD", OCTA_Variable) ~ paste("SVD:", gsub("SVD_", "", OCTA_Variable)),
        TRUE ~ OCTA_Variable
      )
    )
  
  # 创建条形图比较模型性能
  plot <- ggplot(combined_performance, aes(x = reorder(OCTA_Label, !!sym(metric)), 
                                           y = !!sym(metric), 
                                           fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Diabetic" = "#D6604D", "Non-diabetic" = "#4393C3")) +
    labs(
      title = paste("Model Performance Comparison:", metric_label),
      subtitle = "Diabetic vs Non-diabetic Groups",
      x = "OCTA Variables",
      y = metric_label,
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # 保存图表
  ggsave(paste0("octa_model_performance_", metric, "_comparison.pdf"), 
         plot, width = 12, height = 8)
  
  return(plot)
}

# 创建预测值与实际值对比图函数
create_prediction_plot <- function(model, data, octa_var, group_name) {
  # 计算预测值
  data$predicted <- predict(model, newdata = data)
  
  # 创建散点图
  plot <- ggplot(data, aes(x = predicted, y = .data[[octa_var]])) +
    geom_point(alpha = 0.7, color = ifelse(group_name == "Diabetic", "#D6604D", "#4393C3")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = paste("Predicted vs Actual", octa_var),
      subtitle = paste(group_name, "Group"),
      x = "Predicted Values",
      y = "Actual Values"
    ) +
    annotate("text", x = min(data$predicted, na.rm = TRUE), 
             y = max(data[[octa_var]], na.rm = TRUE),
             label = sprintf("R² = %.3f\nRMSE = %.3f", 
                             summary(model)$r.squared,
                             sqrt(mean(resid(model)^2))),
             hjust = 0, vjust = 1) +
    theme_minimal()
  
  # 保存图表
  ggsave(paste0("octa_", octa_var, "_", tolower(gsub(" ", "_", group_name)), "_prediction.pdf"), 
         plot, width = 8, height = 6)
  
  return(plot)
}

# 执行性能评估 - 糖尿病组
cat("\n\n========= OCTA模型性能评估：糖尿病组 =========\n\n")
dm_performance <- evaluate_model_performance(
  mr_thickness_dm$models,  # 假设您的模型存储在这个结构中
  wearable_thickness_dm,   # 糖尿病组厚度数据
  thickness_vars,          # 厚度变量列表
  "Diabetic Thickness"     # 组名
)

# 输出糖尿病组模型性能摘要
print(dm_performance %>% dplyr::select(OCTA_Variable, R_Squared, Adjusted_R_Squared, RMSE, MAE))

# 执行性能评估 - 非糖尿病组
cat("\n\n========= OCTA模型性能评估：非糖尿病组 =========\n\n")
nodm_performance <- evaluate_model_performance(
  mr_thickness_nodm$models,  # 非糖尿病组模型
  wearable_thickness_nodm,   # 非糖尿病组厚度数据
  thickness_vars,            # 厚度变量列表
  "Non-diabetic Thickness"   # 组名
)

# 输出非糖尿病组模型性能摘要
print(nodm_performance %>% dplyr::select(OCTA_Variable, R_Squared, Adjusted_R_Squared, RMSE, MAE))

# 对血流变量重复同样的操作
cat("\n\n========= OCTA模型性能评估：糖尿病组(血流) =========\n\n")
dm_bloodflow_performance <- evaluate_model_performance(
  mr_bloodflow_dm$models,  # 糖尿病组血流模型
  wearable_bloodflow_dm,   # 糖尿病组血流数据
  bloodflow_vars,          # 血流变量列表
  "Diabetic Bloodflow"     # 组名
)

# 输出糖尿病组血流模型性能摘要
print(dm_bloodflow_performance %>% dplyr::select(OCTA_Variable, R_Squared, Adjusted_R_Squared, RMSE, MAE))

cat("\n\n========= OCTA模型性能评估：非糖尿病组(血流) =========\n\n")
nodm_bloodflow_performance <- evaluate_model_performance(
  mr_bloodflow_nodm$models,  # 非糖尿病组血流模型
  wearable_bloodflow_nodm,   # 非糖尿病组血流数据
  bloodflow_vars,            # 血流变量列表
  "Non-diabetic Bloodflow"   # 组名
)

