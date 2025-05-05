library(tidyverse)
library(meta)
library(metafor)
library(gridExtra)
library(gtsummary)
library(dplyr)

# 设置工作目录
setwd(get_project_wd())  # 取决于你的项目设置
rm(list = ls())

# 加载数据
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1w_prediction/daily_data/all_days_combined_data.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness.csv")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# 创建输出目录
dir.create("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_analysis/results_csv")

# 处理OCTA血流数据（基线T0）
octa_bloodflow_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_bloodflow, by = c("ID" = "id"))

# 处理每个患者血流数据的函数
process_patient_bloodflow <- function(patient_data, time_points = c("T0")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术用左眼数据，右眼和双眼手术用右眼数据
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # 处理每个时间点
  for(suffix in time_points) {
    # 选择当前时间点的列
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # 选择数据并重命名列
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# 单独处理每个患者的基线(T0)数据
patient_list_bloodflow <- split(octa_bloodflow_features, octa_bloodflow_features$ID)
processed_data_bloodflow <- purrr::map(patient_list_bloodflow, process_patient_bloodflow)

# 合并结果
octa_bloodflow_features <- bind_rows(processed_data_bloodflow)

# 创建基线(T0)的血流变量子集
bloodflow_var_T0 <- octa_bloodflow_features %>%
  dplyr::select(
    ID,  # 保留ID列
    matches(".*_0_6_T0$"),  # 选择所有以0_6_T0结尾的列
    -matches("PA_OuterRetina_0_6_T0"),  # 排除这些列
    -matches("PA_PED_0_6_T0")
  )

# 处理OCTA厚度数据（基线T0）
octa_thickness_features <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  left_join(octa_thickness, by = c("ID" = "id"))

# 处理每个患者厚度数据的函数
process_patient_thickness <- function(patient_data, time_points = c("T0")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术用左眼数据，右眼和双眼手术用右眼数据
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # 处理每个时间点
  for(suffix in time_points) {
    # 选择当前时间点的列
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # 选择数据并重命名列
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# 单独处理每个患者的基线(T0)数据
patient_list_thickness <- split(octa_thickness_features, octa_thickness_features$ID)
processed_data_thickness <- purrr::map(patient_list_thickness, process_patient_thickness)

# 合并结果
octa_thickness_features <- bind_rows(processed_data_thickness)

# 创建基线(T0)的厚度变量子集
thickness_var_T0 <- octa_thickness_features %>%
  dplyr::select(
    ID,  # 保留ID列
    matches(".*_0_6_T0$"),
    matches("Thickness_PED_0_6_T0")
  )

# 定义OCTA变量用于分析
thickness_vars <- colnames(thickness_var_T0)[-1]  # 移除ID列
bloodflow_vars <- colnames(bloodflow_var_T0)[-1]  # 移除ID列

# 第1步：筛选术前数据（-7天至-1天）
pre_surgery_data <- combined_data %>%
  filter(as.numeric(day) >= -7 & as.numeric(day) <= -1)  # 仅术前7天数据

# 第2步：按受试者ID分组并计算术前均值
patient_preop_means <- pre_surgery_data %>%
  group_by(subject_id) %>%
  summarise(
    # 选择关键均值指标
    mean_hr = mean(mean_hr, na.rm = TRUE),           # 平均心率
    mean_rhr = mean(mean_rhr_1, na.rm = TRUE),       # 平均静息心率
    mean_bo = mean(mean_bo, na.rm = TRUE),           # 平均血氧饱和度
    total_steps = mean(steps_total, na.rm = TRUE),   # 平均总步数
    total_sleep = mean(total_sleep, na.rm = TRUE),   # 平均总睡眠时间
    
    # 保留其他人口统计变量
    dm_2 = dplyr::first(dm_2),                    # 糖尿病状态
    age = dplyr::first(age),                      # 年龄
    gender = dplyr::first(gender),                # 性别
    bmi = dplyr::first(bmi),                      # BMI
    hypertension_2 = dplyr::first(hypertension_2) # 高血压状态
  )

# 第3步：处理缺失值
# 计算缺失值比例
missing_values <- colSums(is.na(patient_preop_means))/nrow(patient_preop_means)
print(missing_values)

# 执行中位数填补
data_imputed <- patient_preop_means
for(col in names(data_imputed)[2:(ncol(data_imputed)-5)]) {  # 跳过subject_id和人口统计变量
  if(sum(is.na(data_imputed[[col]])) > 0) {
    data_imputed[[col]][is.na(data_imputed[[col]])] <- median(data_imputed[[col]], na.rm = TRUE)
  }
}

# 第4步：合并可穿戴设备数据与OCTA数据
wearable_thickness_data <- data_imputed %>%
  left_join(thickness_var_T0, by = c("subject_id" = "ID"))

wearable_bloodflow_data <- data_imputed %>%
  left_join(bloodflow_var_T0, by = c("subject_id" = "ID"))

# 第5步：将数据分为糖尿病组和非糖尿病组
wearable_thickness_dm <- wearable_thickness_data %>% filter(dm_2 == 1)
wearable_thickness_nodm <- wearable_thickness_data %>% filter(dm_2 == 0)

wearable_bloodflow_dm <- wearable_bloodflow_data %>% filter(dm_2 == 1)
wearable_bloodflow_nodm <- wearable_bloodflow_data %>% filter(dm_2 == 0)

# 打印样本容量
cat("糖尿病患者样本量(厚度):", nrow(wearable_thickness_dm), "人\n")
cat("非糖尿病患者样本量(厚度):", nrow(wearable_thickness_nodm), "人\n")
cat("糖尿病患者样本量(血流):", nrow(wearable_bloodflow_dm), "人\n")
cat("非糖尿病患者样本量(血流):", nrow(wearable_bloodflow_nodm), "人\n")

# 第6步：比较糖尿病组与非糖尿病组之间的可穿戴参数
# 创建比较表
comparison_table <- data_imputed %>%
  mutate(
    dm_status = ifelse(dm_2 == 1, "Diabetes", "No Diabetes"),
    gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female"))
  ) %>%
  dplyr::select(mean_hr, mean_rhr, mean_bo, total_steps, total_sleep, 
                age, gender, bmi, dm_status) %>%
  tbl_summary(
    by = dm_status,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    label = list(
      mean_hr ~ "Mean Heart Rate",
      mean_rhr ~ "Mean Resting Heart Rate",
      mean_bo ~ "Mean Blood Oxygen",
      total_steps ~ "Mean Total Steps",
      total_sleep ~ "Mean Total Sleep Time",
      age ~ "Age",
      gender ~ "Gender",
      bmi ~ "BMI"
    )
  ) %>%
  add_p() %>%
  add_n() %>%
  modify_header(label = "**Parameter**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Comparison by Diabetes Status**") %>%
  modify_footnote(
    all_stat_cols() ~ "Mean (SD) or n (%)",
    p.value ~ "Wilcoxon rank-sum test; Chi-square test"
  )

# 打印比较表
print(comparison_table)

# 保存比较表为HTML（如果gt包可用）
if(requireNamespace("gt", quietly = TRUE)) {
  comparison_table %>%
    as_gt() %>%
    gt::gtsave("diabetic_nondiabetic_comparison_table.html")
}

# 保存比较数据为CSV以供可视化使用
comparison_data <- data_imputed %>%
  mutate(dm_status = ifelse(dm_2 == 1, "Diabetes", "No Diabetes"),
         gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female")))
write.csv(comparison_data, "comparison_data.csv", row.names = FALSE)

# 定义可穿戴变量和协变量用于分析
wearable_vars <- c("mean_hr", "mean_rhr", "mean_bo", "total_steps", "total_sleep")
covariate_vars <- c("age", "gender", "bmi")

# 第7步：相关性分析 - OCTA变量与可穿戴设备变量

# 创建函数计算一个OCTA变量的相关性
calculate_correlation_octa <- function(data, octa_var, wearable_var) {
  # 移除含有NA值的行
  complete_data <- data %>%
    dplyr::select(all_of(c(octa_var, wearable_var))) %>%
    na.omit()
  
  if(nrow(complete_data) < 5) {  # 至少需要5个观测值
    return(data.frame(
      OCTA_Variable = octa_var,
      Wearable_Variable = wearable_var,
      Correlation = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      P_value = NA,
      N = nrow(complete_data)
    ))
  }
  
  # Pearson相关系数
  correlation <- cor(complete_data[[octa_var]], complete_data[[wearable_var]], 
                     use = "complete.obs")
  
  # Fisher z变换计算置信区间
  n <- nrow(complete_data)
  z <- 0.5 * log((1 + correlation) / (1 - correlation))
  se <- 1 / sqrt(n - 3)
  ci_lower <- tanh(z - 1.96 * se)
  ci_upper <- tanh(z + 1.96 * se)
  
  # 计算p值
  t_stat <- correlation * sqrt((n - 2) / (1 - correlation^2))
  p_value <- 2 * pt(-abs(t_stat), df = n - 2)
  
  return(data.frame(
    OCTA_Variable = octa_var,
    Wearable_Variable = wearable_var,
    Correlation = correlation,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    P_value = p_value,
    N = n
  ))
}

# 计算一组OCTA变量相关性的函数
calculate_all_correlations <- function(data, octa_vars, wearable_vars, group_name) {
  all_correlations <- data.frame()
  
  for(octa_var in octa_vars) {
    for(wearable_var in wearable_vars) {
      correlation_result <- calculate_correlation_octa(data, octa_var, wearable_var)
      correlation_result$Group <- group_name
      all_correlations <- rbind(all_correlations, correlation_result)
    }
  }
  
  # 添加FDR校正处理多重检验
  all_correlations$P_adjusted <- p.adjust(all_correlations$P_value, method = "BH")
  
  return(all_correlations)
}

# 计算厚度变量的相关性
thickness_dm_correlations <- calculate_all_correlations(
  wearable_thickness_dm, thickness_vars, wearable_vars, "Diabetes")

thickness_nodm_correlations <- calculate_all_correlations(
  wearable_thickness_nodm, thickness_vars, wearable_vars, "No Diabetes")

# 计算血流变量的相关性
bloodflow_dm_correlations <- calculate_all_correlations(
  wearable_bloodflow_dm, bloodflow_vars, wearable_vars, "Diabetes")

bloodflow_nodm_correlations <- calculate_all_correlations(
  wearable_bloodflow_nodm, bloodflow_vars, wearable_vars, "No Diabetes")

# 合并所有相关性
all_thickness_correlations <- rbind(thickness_dm_correlations, thickness_nodm_correlations)
all_bloodflow_correlations <- rbind(bloodflow_dm_correlations, bloodflow_nodm_correlations)

# 保存相关性结果
write.csv(all_thickness_correlations, "thickness_correlations.csv", row.names = FALSE)
write.csv(all_bloodflow_correlations, "bloodflow_correlations.csv", row.names = FALSE)

# 打印显著相关性
cat("\n显著厚度相关性 (p < 0.05):\n")
print(all_thickness_correlations %>% 
        filter(P_value < 0.05) %>% 
        arrange(P_value))

cat("\n显著血流相关性 (p < 0.05):\n")
print(all_bloodflow_correlations %>% 
        filter(P_value < 0.05) %>% 
        arrange(P_value))


# 第8步：为顶级OCTA变量建立多元回归模型的函数

# 选择顶级OCTA变量的函数
select_top_octa_vars <- function(correlations, n = 5) {
  # 如果不是数据框，转换为数据框
  corr_df <- as.data.frame(correlations)

  # 过滤掉NA值
  corr_df <- corr_df[!is.na(corr_df$Correlation), ]

  # 添加绝对相关性列
  corr_df$abs_corr <- abs(corr_df$Correlation)

  # 为每个OCTA变量计算最大绝对相关性
  result <- aggregate(abs_corr ~ OCTA_Variable, data = corr_df, FUN = max)

  # 按降序排列相关性
  result <- result[order(-result$abs_corr), ]

  # 选择前n个变量
  if(nrow(result) > n) {
    top_vars <- result$OCTA_Variable[1:n]
  } else {
    top_vars <- result$OCTA_Variable
  }

  return(top_vars)
}

# 建立回归模型的函数
build_regression_models <- function(data, octa_var, wearable_vars, covariates, group_name) {
  # 准备公式
  formula_str <- paste(octa_var, "~", paste(c(wearable_vars, covariates), collapse = " + "))

  # 创建回归模型
  model <- tryCatch({
    lm(as.formula(formula_str), data = data)
  }, error = function(e) {
    cat("Error in creating model for", octa_var, "in", group_name, "group:", e$message, "\n")
    return(NULL)
  })

  # 如果模型创建成功
  if(!is.null(model)) {
    # 提取系数和置信区间
    coef_table <- summary(model)$coefficients
    conf_int <- confint(model)

    # 合并到数据框
    result <- data.frame(
      OCTA_Variable = octa_var,
      Variable = rownames(coef_table),
      Coefficient = coef_table[, "Estimate"],
      Std_Error = coef_table[, "Std. Error"],
      t_value = coef_table[, "t value"],
      P_value = coef_table[, "Pr(>|t|)"],
      Lower_CI = conf_int[, "2.5 %"],
      Upper_CI = conf_int[, "97.5 %"],
      Group = group_name,
      stringsAsFactors = FALSE
    ) %>% filter(Variable != "(Intercept)")

    # 添加模型性能指标
    model_summary <- summary(model)
    r_squared <- model_summary$r.squared
    adj_r_squared <- model_summary$adj.r.squared

    # 计算RMSE
    rmse <- sqrt(mean(model_summary$residuals^2))

    # 返回结果和模型
    return(list(
      coef_data = result,
      model = model,
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      rmse = rmse
    ))
  } else {
    return(NULL)
  }
}

# 选择顶级变量
top_thickness_vars <- select_top_octa_vars(all_thickness_correlations)
top_bloodflow_vars <- select_top_octa_vars(all_bloodflow_correlations)

# 将顶级变量保存到CSV文件
write.csv(data.frame(Variable = top_thickness_vars, Type = "Thickness"),
          "top_thickness_vars.csv", row.names = FALSE)
write.csv(data.frame(Variable = top_bloodflow_vars, Type = "Blood Flow"),
          "top_bloodflow_vars.csv", row.names = FALSE)

# 为顶级厚度变量建立模型
thickness_dm_models <- list()
thickness_nodm_models <- list()
thickness_dm_results <- data.frame()
thickness_nodm_results <- data.frame()
thickness_model_performance <- data.frame()

for(var in top_thickness_vars) {
  # 糖尿病组
  dm_result <- build_regression_models(
    wearable_thickness_dm, var, wearable_vars, covariate_vars, "Diabetes"
  )

  if(!is.null(dm_result)) {
    thickness_dm_models[[var]] <- dm_result$model
    thickness_dm_results <- rbind(thickness_dm_results, dm_result$coef_data)

    # 添加性能指标
    thickness_model_performance <- rbind(
      thickness_model_performance,
      data.frame(
        OCTA_Variable = var,
        Group = "Diabetes",
        R_squared = dm_result$r_squared,
        Adj_R_squared = dm_result$adj_r_squared,
        RMSE = dm_result$rmse,
        stringsAsFactors = FALSE
      )
    )
  }

  # 非糖尿病组
  nodm_result <- build_regression_models(
    wearable_thickness_nodm, var, wearable_vars, covariate_vars, "No Diabetes"
  )

  if(!is.null(nodm_result)) {
    thickness_nodm_models[[var]] <- nodm_result$model
    thickness_nodm_results <- rbind(thickness_nodm_results, nodm_result$coef_data)

    # 添加性能指标
    thickness_model_performance <- rbind(
      thickness_model_performance,
      data.frame(
        OCTA_Variable = var,
        Group = "No Diabetes",
        R_squared = nodm_result$r_squared,
        Adj_R_squared = nodm_result$adj_r_squared,
        RMSE = nodm_result$rmse,
        stringsAsFactors = FALSE
      )
    )
  }
}

# 为顶级血流变量建立模型
bloodflow_dm_models <- list()
bloodflow_nodm_models <- list()
bloodflow_dm_results <- data.frame()
bloodflow_nodm_results <- data.frame()
bloodflow_model_performance <- data.frame()

for(var in top_bloodflow_vars) {
  # 糖尿病组
  dm_result <- build_regression_models(
    wearable_bloodflow_dm, var, wearable_vars, covariate_vars, "Diabetes"
  )

  if(!is.null(dm_result)) {
    bloodflow_dm_models[[var]] <- dm_result$model
    bloodflow_dm_results <- rbind(bloodflow_dm_results, dm_result$coef_data)

    # 添加性能指标
    bloodflow_model_performance <- rbind(
      bloodflow_model_performance,
      data.frame(
        OCTA_Variable = var,
        Group = "Diabetes",
        R_squared = dm_result$r_squared,
        Adj_R_squared = dm_result$adj_r_squared,
        RMSE = dm_result$rmse,
        stringsAsFactors = FALSE
      )
    )
  }

  # 非糖尿病组
  nodm_result <- build_regression_models(
    wearable_bloodflow_nodm, var, wearable_vars, covariate_vars, "No Diabetes"
  )

  if(!is.null(nodm_result)) {
    bloodflow_nodm_models[[var]] <- nodm_result$model
    bloodflow_nodm_results <- rbind(bloodflow_nodm_results, nodm_result$coef_data)

    # 添加性能指标
    bloodflow_model_performance <- rbind(
      bloodflow_model_performance,
      data.frame(
        OCTA_Variable = var,
        Group = "No Diabetes",
        R_squared = nodm_result$r_squared,
        Adj_R_squared = nodm_result$adj_r_squared,
        RMSE = nodm_result$rmse,
        stringsAsFactors = FALSE
      )
    )
  }
}

# 合并所有回归结果
all_thickness_reg_results <- rbind(thickness_dm_results, thickness_nodm_results)
all_bloodflow_reg_results <- rbind(bloodflow_dm_results, bloodflow_nodm_results)

# 保存回归结果
write.csv(all_thickness_reg_results, "thickness_regression_results.csv", row.names = FALSE)
write.csv(all_bloodflow_reg_results, "bloodflow_regression_results.csv", row.names = FALSE)
write.csv(rbind(thickness_model_performance, bloodflow_model_performance),
          "model_performance.csv", row.names = FALSE)

# 最后，保存最佳模型的预测数据
save_predictions <- function(model_list, data, result_file) {
  predictions <- data.frame(subject_id = data$subject_id)

  for(var_name in names(model_list)) {
    if(!is.null(model_list[[var_name]])) {
      # 预测值
      pred_values <- predict(model_list[[var_name]], newdata = data)
      predictions[[paste0(var_name, "_pred")]] <- pred_values

      # 实际值
      predictions[[paste0(var_name, "_actual")]] <- data[[var_name]]
    }
  }

  # 添加糖尿病状态
  predictions$dm_status <- ifelse(data$dm_2 == 1, "Diabetes", "No Diabetes")

  # 保存到CSV
  write.csv(predictions, result_file, row.names = FALSE)
}

# 保存模型预测结果
save_predictions(thickness_dm_models, wearable_thickness_dm, "thickness_dm_predictions.csv")
save_predictions(thickness_nodm_models, wearable_thickness_nodm, "thickness_nodm_predictions.csv")
save_predictions(bloodflow_dm_models, wearable_bloodflow_dm, "bloodflow_dm_predictions.csv")
save_predictions(bloodflow_nodm_models, wearable_bloodflow_nodm, "bloodflow_nodm_predictions.csv")

# 整合上面两个保存预测结果的函数，为可视化准备输入
prepare_top_model_data <- function(model_performance, dm_models, nodm_models,
                                   dm_data, nodm_data, type = "thickness") {
  # 找到每组R²最高的模型
  top_dm_var <- model_performance %>%
    filter(Group == "Diabetes") %>%
    arrange(desc(R_squared)) %>%
    slice(1) %>%
    pull(OCTA_Variable)

  top_nodm_var <- model_performance %>%
    filter(Group == "No Diabetes") %>%
    arrange(desc(R_squared)) %>%
    slice(1) %>%
    pull(OCTA_Variable)

  # 创建结果数据框
  result <- data.frame()

  # 糖尿病组
  if(!is.null(dm_models[[top_dm_var]])) {
    dm_pred <- data.frame(
      subject_id = dm_data$subject_id,
      octa_var = top_dm_var,
      actual = dm_data[[top_dm_var]],
      predicted = predict(dm_models[[top_dm_var]], newdata = dm_data),
      Group = "Diabetes",
      r_squared = summary(dm_models[[top_dm_var]])$r.squared,
      rmse = sqrt(mean(resid(dm_models[[top_dm_var]])^2)),
      type = type
    )
    result <- rbind(result, dm_pred)
  }

  # 非糖尿病组
  if(!is.null(nodm_models[[top_nodm_var]])) {
    nodm_pred <- data.frame(
      subject_id = nodm_data$subject_id,
      octa_var = top_nodm_var,
      actual = nodm_data[[top_nodm_var]],
      predicted = predict(nodm_models[[top_nodm_var]], newdata = nodm_data),
      Group = "No Diabetes",
      r_squared = summary(nodm_models[[top_nodm_var]])$r.squared,
      rmse = sqrt(mean(resid(nodm_models[[top_nodm_var]])^2)),
      type = type
    )
    result <- rbind(result, nodm_pred)
  }

  return(result)
}

# 准备顶级模型数据
top_thickness_model_data <- prepare_top_model_data(
  thickness_model_performance,
  thickness_dm_models,
  thickness_nodm_models,
  wearable_thickness_dm,
  wearable_thickness_nodm,
  "thickness"
)

top_bloodflow_model_data <- prepare_top_model_data(
  bloodflow_model_performance,
  bloodflow_dm_models,
  bloodflow_nodm_models,
  wearable_bloodflow_dm,
  wearable_bloodflow_nodm,
  "bloodflow"
)

# 合并顶级模型数据并保存
all_top_model_data <- rbind(top_thickness_model_data, top_bloodflow_model_data)
write.csv(all_top_model_data, "top_model_predictions.csv", row.names = FALSE)

# 第9步：结果摘要

# 显著相关性
cat("\n\n=== 结果摘要 ===\n\n")

cat("厚度与可穿戴参数的显著相关性：\n")
significant_thickness <- all_thickness_correlations %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value) %>%
  head(10) %>%
  dplyr::select(Group, OCTA_Variable, Wearable_Variable, Correlation, P_value)
print(significant_thickness)
write.csv(significant_thickness, "significant_thickness_correlations.csv", row.names = FALSE)

cat("\n血流与可穿戴参数的显著相关性：\n")
significant_bloodflow <- all_bloodflow_correlations %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value) %>%
  head(10) %>%
  dplyr::select(Group, OCTA_Variable, Wearable_Variable, Correlation, P_value)
print(significant_bloodflow)
write.csv(significant_bloodflow, "significant_bloodflow_correlations.csv", row.names = FALSE)

# 显著回归系数
cat("\n厚度模型的显著回归系数：\n")
significant_thickness_reg <- all_thickness_reg_results %>% 
  filter(P_value < 0.05, Variable %in% wearable_vars) %>% 
  arrange(P_value) %>%
  head(10) %>%
  dplyr::select(Group, OCTA_Variable, Variable, Coefficient, P_value)
print(significant_thickness_reg)
write.csv(significant_thickness_reg, "significant_thickness_regression.csv", row.names = FALSE)

cat("\n血流模型的显著回归系数：\n")
significant_bloodflow_reg <- all_bloodflow_reg_results %>% 
  filter(P_value < 0.05, Variable %in% wearable_vars) %>% 
  arrange(P_value) %>%
  head(10) %>%
  dplyr::select(Group, OCTA_Variable, Variable, Coefficient, P_value)
print(significant_bloodflow_reg)
write.csv(significant_bloodflow_reg, "significant_bloodflow_regression.csv", row.names = FALSE)

# 模型性能摘要
cat("\n厚度最优模型：\n")
best_thickness_models <- thickness_model_performance %>% 
  group_by(Group) %>%
  arrange(desc(R_squared)) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(Group, OCTA_Variable, R_squared, RMSE)
print(best_thickness_models)
write.csv(best_thickness_models, "best_thickness_models.csv", row.names = FALSE)

cat("\n血流最优模型：\n")
best_bloodflow_models <- bloodflow_model_performance %>% 
  group_by(Group) %>%
  arrange(desc(R_squared)) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(Group, OCTA_Variable, R_squared, RMSE)
print(best_bloodflow_models)
write.csv(best_bloodflow_models, "best_bloodflow_models.csv", row.names = FALSE)

# 为长格式数据创建一个文件，用于箱线图
long_data <- data_imputed %>%
  pivot_longer(
    cols = c(mean_hr, mean_rhr, mean_bo, total_steps, total_sleep),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    dm_status = ifelse(dm_2 == 1, "Diabetes", "No Diabetes"),
    parameter_label = case_when(
      parameter == "mean_hr" ~ "Mean Heart Rate",
      parameter == "mean_rhr" ~ "Mean Resting Heart Rate",
      parameter == "mean_bo" ~ "Mean Blood Oxygen",
      parameter == "total_steps" ~ "Mean Total Steps",
      parameter == "total_sleep" ~ "Mean Total Sleep"
    )
  )
write.csv(long_data, "long_format_data.csv", row.names = FALSE)

# 保存原始数据供可视化参考
write.csv(wearable_thickness_data, "wearable_thickness_data.csv", row.names = FALSE)
write.csv(wearable_bloodflow_data, "wearable_bloodflow_data.csv", row.names = FALSE)