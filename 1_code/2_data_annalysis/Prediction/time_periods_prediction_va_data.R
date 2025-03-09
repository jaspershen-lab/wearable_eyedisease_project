library(tidyverse)
library(tidymass)
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(lubridate)
library(caret)
library(randomForest)
library(mice)



# 加载数据
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/time_period_rhr_results.rda")
load("3_data_analysis/2_data_analysis/bo/bo_time_periods/bo_time_period_results.rda")
load("3_data_analysis/2_data_analysis/sleep/sleep_time_periods/sleep_time_period_results.rda")
load("3_data_analysis/2_data_analysis/steps/steps_time_periods/steps_time_period_results.rda")

baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness.csv")

# 创建输出目录
dir.create("3_data_analysis/2_data_analysis/Prediction/model_data_final_time_periods", recursive = TRUE)
setwd("3_data_analysis/2_data_analysis/Prediction/model_data_final_time_periods")

# 定义时间窗口 - 这是我们要循环处理的三个时间窗口
time_periods <- c("pre_surgery_3d", "pre_surgery_3d_to_7d", "pre_surgery_7d_all")

# 加载和处理不依赖于时间窗口的数据（只需要处理一次）

####disease data
disease_data <- baseline_info %>%
  mutate(
    cataract_2 = case_when(
      cataract == 1 ~ 0,  
      cataract %in% c(2, 3, 4) ~ 1, 
      TRUE ~ NA_real_
    ),
    dm_2 = case_when(
      diabetes_history == 1 ~ 1,  
      diabetes_history == 2 ~ 0, 
      TRUE ~ NA_real_
    ),
    hypertension_2 = case_when(
      hypertension_history == 1 ~ 1,  
      hypertension_history == 2 ~ 0, 
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(ID, cataract_2, dm_2, hypertension_2)

# 处理视力数据
vision_data <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,  # 右眼手术
      surgery_eye_1 == 1 ~ os_corrected_bas,  # 左眼手术
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,  # 双眼手术取平均
      TRUE ~ NA_real_
    ),
    post_vision= case_when(
      surgery_eye_1 == 0 ~ od_corrected_1w,   # 右眼手术后
      surgery_eye_1 == 1 ~ os_corrected_1w,   # 左眼手术后
      surgery_eye_1 == 2 ~ (od_corrected_1w + os_corrected_1w)/2,   # 双眼手术后取平均
      TRUE ~ NA_real_
    ),
    vision_improvement = post_vision - pre_vision,
    vision_improved = if_else(vision_improvement >= 0, 1, 0),  # 添加二分类指标
    vision_improved_factor = factor(vision_improved, 
                                    levels = c(0, 1), 
                                    labels = c("NoImprovement", "Improved"))  # 添加因子版本
  ) %>%
  dplyr::select(ID, surgery_eye_1, pre_vision, post_vision, 
                vision_improvement, vision_improved, vision_improved_factor,
                age, gender)

# 处理OCTA bloodflow数据
octa_bloodflow_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_bloodflow, by = c("ID" = "id"))

# 处理单个患者的数据
process_patient <- function(patient_data) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术选择左眼数据，右眼和双眼手术选择右眼数据
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  # 选择列
  cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
  cols_to_keep <- cols_to_keep[grep("T0$", cols_to_keep)]
  
  # 选择数据并重命名
  result <- patient_data %>%
    dplyr::select("ID", all_of(cols_to_keep)) %>%
    rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
  
  return(result)
}

# 对每个患者分别处理
patient_list <- split(octa_bloodflow_features, octa_bloodflow_features$ID)
processed_data <- purrr::map(patient_list, process_patient)

# 合并结果
octa_bloodflow_features <- bind_rows(processed_data)

bloodflow_var <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # 保留ID列
    matches(".*_0_6_T0$"),  # 选择所有以 0_6_T0 结尾的列
    -"PA_OuterRetina_0_6_T0",  # 直接用列名排除
    -"PA_PED_0_6_T0"
  )

# 处理OCTA thickness数据
octa_thickness_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_thickness, by = c("ID" = "id"))

# 对每个患者分别处理
patient_list <- split(octa_thickness_features, octa_thickness_features$ID)
processed_data <- purrr::map(patient_list, process_patient)

# 合并结果
octa_thickness_features <- bind_rows(processed_data)

thickness_var <- octa_thickness_features %>%
  dplyr::select(
    ID,  # 保留ID列
    matches(".*_0_6_T0$"),
    "Thickness_PED_0_6_T0"
  )

# 定义处理可穿戴数据异常值的函数
detect_and_handle_outliers <- function(data, patterns, z_threshold = 3, replace_with_na = TRUE) {
  all_outliers <- list()
  
  for(col in colnames(data)) {
    # 检查列名是否匹配模式
    is_matched <- any(sapply(patterns, function(pattern) grepl(pattern, col)))
    
    if(is_matched && is.numeric(data[[col]])) {
      # 获取非NA值
      valid_values <- data[[col]][!is.na(data[[col]])]
      
      if(length(valid_values) > 5) {  # 确保有足够的数据计算z分数
        # 计算z-scores
        valid_mean <- mean(valid_values, na.rm = TRUE)
        valid_sd <- sd(valid_values, na.rm = TRUE)
        
        if(valid_sd > 0) {  # 避免除以0
          z_scores <- (data[[col]] - valid_mean) / valid_sd
          
          # 找出超过阈值的值
          outliers <- which(abs(z_scores) > z_threshold & !is.na(z_scores))
          
          # 存储和处理异常值信息
          if(length(outliers) > 0) {
            all_outliers[[col]] <- list(
              indices = outliers,
              values = data[[col]][outliers],
              z_scores = z_scores[outliers]
            )
            
            # 根据参数决定是否替换为NA
            if(replace_with_na) {
              data[outliers, col] <- NA
            }
          }
        }
      }
    }
  }
  
  return(list(data = data, outliers = all_outliers))
}

# 定义PCA函数
perform_pca <- function(data, type) {
  # 移除全为NA的列
  valid_cols <- colSums(is.na(data)) < nrow(data)
  if(any(!valid_cols)) {
    data <- data[, valid_cols, drop = FALSE]
  }
  
  # 处理剩余的NA值（均值填充）
  for(col in colnames(data)) {
    if(any(is.na(data[[col]]))) {
      col_mean <- mean(data[[col]], na.rm = TRUE)
      data[[col]][is.na(data[[col]])] <- col_mean
    }
  }
  
  # 标准化数据
  scaled_data <- scale(data)
  
  # 检查无限值
  if(any(!is.finite(scaled_data))) {
    scaled_data[!is.finite(scaled_data)] <- 0
  }
  
  # 执行PCA
  pca_result <- prcomp(scaled_data, scale. = FALSE)
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  n_components <- which(cumsum(var_explained) >= 0.8)[1]
  
  # 提取PC分数
  pc_scores <- pca_result$x[, 1:n_components]
  if(type == "Bloodflow") {
    colnames(pc_scores) <- paste0("OCTA_Flow_PC", 1:n_components)
  } else {
    colnames(pc_scores) <- paste0("OCTA_Thick_PC", 1:n_components)
  }
  
  return(list(
    pca = pca_result,
    var_explained = var_explained,
    n_components = n_components,
    pc_scores = pc_scores
  ))
}

# 为每个时间窗口创建数据集的函数
create_dataset_for_time_period <- function(time_period) {
  
  cat("\n\n=============== 处理时间窗口:", time_period, "===============\n")
  
  # 1. 提取RHR特征
  if(time_period == "pre_surgery_3d") {
    rhr_column_suffix <- "pre_surgery_3d"
  } else if(time_period == "pre_surgery_3d_to_7d") {
    rhr_column_suffix <- "pre_surgery_3d_to_7d"
  } else { # pre_surgery_7d_all
    rhr_column_suffix <- "pre_surgery_7d_all"
  }
  
  rhr_features <- time_period_rhr_results %>%
    dplyr::select(
      subject_id,
      rhr_mean = paste0("mean_hr_rhr_steps_1_", rhr_column_suffix),
      rhr_min = paste0("min_hr_rhr_steps_1_", rhr_column_suffix),
      rhr_max = paste0("max_hr_rhr_steps_1_", rhr_column_suffix),
      rhr_median = paste0("median_hr_rhr_steps_1_", rhr_column_suffix),
      rhr_sd = paste0("sd_hr_rhr_steps_1_", rhr_column_suffix),
      rhr_iqr = paste0("iqr_hr_rhr_steps_1_", rhr_column_suffix)
    )
  
  # 2. 提取血氧特征
  bo_features <- bo_time_period_results %>%
    dplyr::select(
      subject_id,
      bo_mean = paste0(time_period, "_mean"),
      bo_min = paste0(time_period, "_min"),
      bo_max = paste0(time_period, "_max"),
      bo_median = paste0(time_period, "_median"),
      bo_sd = paste0(time_period, "_sd"),
      bo_iqr = paste0(time_period, "_iqr")
    )
  
  # 3. 提取睡眠特征
  sleep_features <- sleep_time_period_results %>%
    dplyr::select(
      subject_id,
      deep_sleep = paste0(time_period, "_deep_sleep_total"),
      light_sleep = paste0(time_period, "_light_sleep_total"),
      dream_sleep = paste0(time_period, "_dream_sleep_total"),
      total_sleep = paste0(time_period, "_total_sleep_total"),
      awake = paste0(time_period, "_awake_total"),
      daytime_sleep = paste0(time_period, "_daytime_sleep_total")
    )
  
  # 4. 提取步数特征
  steps_features <- steps_time_period_results %>%
    dplyr::select(
      subject_id,
      steps_total = paste0(time_period, "_total"),
      steps_mean = paste0(time_period, "_mean"),
      steps_median = paste0(time_period, "_median"),
      steps_max = paste0(time_period, "_max")
    )
  
  # 5. 合并所有特征
  model_data <- rhr_features %>%
    inner_join(disease_data, by = c("subject_id" = "ID")) %>%
    inner_join(vision_data, by = c("subject_id" = "ID")) %>%
    inner_join(bo_features, by = "subject_id") %>%
    inner_join(sleep_features, by = "subject_id") %>%
    inner_join(steps_features, by = "subject_id") %>%
    inner_join(bloodflow_var, by = c("subject_id" = "ID")) %>%
    inner_join(thickness_var, by = c("subject_id" = "ID"))
  
  # 6. 将分类变量转换为因子
  model_data$gender <- as.factor(model_data$gender)
  model_data$cataract_2 <- as.factor(model_data$cataract_2)
  model_data$dm_2 <- as.factor(model_data$dm_2)
  model_data$hypertension_2 <- as.factor(model_data$hypertension_2)
  
  # 7. 移除vision_improvement为NA的行（目标变量）
  model_data_clean <- model_data[!is.na(model_data$vision_improvement), ]
  
  # 8. 选择用于填补的特征
  features_for_imputation <- model_data_clean %>%
    dplyr::select(
      # RHR特征
      starts_with("rhr_"),
      # 血氧特征
      starts_with("bo_"),
      # 睡眠特征
      deep_sleep, light_sleep, total_sleep, dream_sleep,
      awake, daytime_sleep,
      # 步数特征
      starts_with("steps_"),
      # 人口统计学和病史
      age, gender, cataract_2, dm_2, hypertension_2, pre_vision,
      # OCTA特征
      matches("(SVD|PA|VD).*_0_6_T0"),
      matches("Thickness.*_0_6_T0"),
      # 目标变量
      vision_improvement
    )
  
  # 9. 处理可穿戴数据异常值
  wearable_patterns <- c("^rhr_", "^bo_", "^steps_", "^deep_sleep$", "^light_sleep$", 
                         "^total_sleep$", "^dream_sleep$", "^awake$", "^daytime_sleep$")
  
  wearable_result <- detect_and_handle_outliers(
    features_for_imputation, 
    wearable_patterns, 
    z_threshold = 3, 
    replace_with_na = TRUE
  )
  features_for_imputation <- wearable_result$data
  
  # 10. 处理OCTA数据异常值
  octa_patterns <- c("(SVD|PA|VD).*_0_6_T0$", "Thickness.*_0_6_T0$")
  
  octa_result <- detect_and_handle_outliers(
    features_for_imputation, 
    octa_patterns, 
    z_threshold = 4.5,
    replace_with_na = TRUE
  )
  features_for_imputation <- octa_result$data
  
  # 11. 执行多重填补
  set.seed(123)  # 为了可重复性
  imp <- mice(features_for_imputation, 
              m = 5,
              maxit = 10,
              ridge = 0.1,
              printFlag = FALSE)
  
  # 使用第一个填补数据集
  complete_data <- complete(imp, 1)
  
  # 确保vision_improved是因子型
  complete_data$vision_improved <- factor(ifelse(complete_data$vision_improvement >= 0, 1, 0),
                                          levels = c(0, 1))
  
  # 12. OCTA PCA
  # 获取血流变量
  bloodflow_vars <- grep("(SVD|PA|VD).*_0_6_T0$", colnames(complete_data), value = TRUE)
  bloodflow_data <- complete_data[, bloodflow_vars]
  
  # 获取厚度变量
  thickness_vars <- grep("Thickness.*_0_6_T0$", colnames(complete_data), value = TRUE)
  thickness_data <- complete_data[, thickness_vars]
  
  # 执行PCA
  bloodflow_pca <- perform_pca(bloodflow_data, "Bloodflow")
  thickness_pca <- perform_pca(thickness_data, "Thickness")
  
  # 13. 创建最终数据集
  model_data_final <- cbind(
    dplyr::select(complete_data, -all_of(c(bloodflow_vars, thickness_vars))),
    bloodflow_pca$pc_scores,
    thickness_pca$pc_scores
  )
  
  # 14. 输出信息
  cat("\n时间窗口:", time_period)
  cat("\n数据集行数:", nrow(model_data_final))
  cat("\n数据集列数:", ncol(model_data_final))
  
  # 15. 保存结果
  filename_suffix <- gsub("pre_surgery_", "", time_period)
  save(complete_data, file = paste0("complete_data_", filename_suffix, ".rda"), compress = "xz")
  save(model_data_final, file = paste0("model_data_final_", filename_suffix, ".rda"), compress = "xz")
  
  return(model_data_final)
}

# 为每个时间窗口创建数据集
results <- list()
for (period in time_periods) {
  results[[period]] <- create_dataset_for_time_period(period)
}

# 打印每个数据集的维度
for (period in time_periods) {
  cat("\n", period, "数据集维度:", dim(results[[period]]))
}

# 在循环结束后添加此代码
for (period in time_periods) {
  if (!is.null(results[[period]])) {
    cat("\n", period, "数据集行数:", nrow(results[[period]]))
  } else {
    cat("\n", period, "数据集不可用")
  }
}


# 由于现有数据集没有ID，我们可以直接使用行号来匹配相应的观测
comparison_data <- data.frame(
  row_id = 1:nrow(results[["pre_surgery_3d"]]),  # 使用行号作为ID
  vision_imp_3d = results[["pre_surgery_3d"]]$vision_improvement,
  vision_imp_3d_to_7d = results[["pre_surgery_3d_to_7d"]]$vision_improvement,
  vision_imp_7d_all = results[["pre_surgery_7d_all"]]$vision_improvement
)

# 输出比较数据的前几行和维度
cat("\n比较数据维度:", dim(comparison_data))
cat("\n比较数据前三行:\n")
print(head(comparison_data, 3))

# 最后创建一个示例比较三个时间窗口的视力结果的可视化
library(tidyr)
library(ggplot2)

# 转换为长格式
long_data <- pivot_longer(comparison_data, 
                          cols = starts_with("vision_imp"), 
                          names_to = "time_window", 
                          values_to = "improvement")

# 创建箱线图比较
ggplot(long_data, aes(x = time_window, y = improvement)) +
  geom_boxplot() +
  labs(title = "视力改善比较 - 不同时间窗口",
       x = "时间窗口",
       y = "视力改善值") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 保存比较数据
save(comparison_data, file = "time_periods_comparison_data.rda", compress = "xz")
