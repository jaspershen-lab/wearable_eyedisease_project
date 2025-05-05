library(tidyverse)
library(meta)
library(metafor)
library(gridExtra)
library(gtsummary)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# 设置工作目录
setwd(get_project_wd())  # 取决于你的项目设置
rm(list = ls())

# 加载数据
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1w_prediction/daily_data/all_days_combined_data.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow_1.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness_1.csv")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

# 创建输出目录
dir.create("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_compare/results_csv", recursive = TRUE, showWarnings = FALSE)
setwd("3_data_analysis/5_presurgery_analysis/basline_prediction/octa_compare/results_csv")

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
    matches(".*_0_21_T0$"),  # 选择所有以0_21_T0结尾的列
    # -matches("PA_OuterRetina_0_21_T0"),  # 排除这些列
    # -matches("PA_PED_0_21_T0")
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
    matches(".*_0_21_T0$"),
    matches("Thickness_PED_0_21_T0")
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




##########analysis
# OCTA与可穿戴设备数据的两步分析流程
# 第一步：未校正的pearson相关分析
# 第二步：对显著关联进行多元线性回归分析（控制年龄、性别和BMI）

# 定义可穿戴变量和协变量用于分析
wearable_vars <- c("mean_rhr", "mean_bo", "total_steps", "total_sleep")
covariate_vars <- c("age", "gender", "bmi")

# 第一步：创建未校正的pearson相关分析函数
calculate_pearson_correlation <- function(data, octa_var, wearable_var) {
  # 移除含有NA值的行
  complete_data <- data %>%
    dplyr::select(all_of(c(octa_var, wearable_var))) %>%
    na.omit()
  
  if(nrow(complete_data) < 5) {  # 至少需要5个观测值
    return(data.frame(
      OCTA_Variable = octa_var,
      Wearable_Variable = wearable_var,
      Correlation = NA,
      P_value = NA,
      N = nrow(complete_data)
    ))
  }
  
  # pearson相关系数
  pearson_test <- cor.test(complete_data[[octa_var]], complete_data[[wearable_var]], 
                            method = "pearson", exact = FALSE)
  
  correlation <- pearson_test$estimate
  p_value <- pearson_test$p.value
  n <- nrow(complete_data)
  
  return(data.frame(
    OCTA_Variable = octa_var,
    Wearable_Variable = wearable_var,
    Correlation = correlation,
    P_value = p_value,
    N = n
  ))
}

# 计算一组OCTA变量未校正相关性的函数
calculate_all_correlations <- function(data, octa_vars, wearable_vars, group_name) {
  all_correlations <- data.frame()
  
  for(octa_var in octa_vars) {
    for(wearable_var in wearable_vars) {
      correlation_result <- calculate_pearson_correlation(data, octa_var, wearable_var)
      correlation_result$Group <- group_name
      all_correlations <- rbind(all_correlations, correlation_result)
    }
  }
  
  return(all_correlations)
}

# 计算厚度变量的相关性（未校正）
thickness_dm_correlations <- calculate_all_correlations(
  wearable_thickness_dm, thickness_vars, wearable_vars, "DR")
thickness_nodm_correlations <- calculate_all_correlations(
  wearable_thickness_nodm, thickness_vars, wearable_vars, "Cataract")

# 计算血流变量的相关性（未校正）
bloodflow_dm_correlations <- calculate_all_correlations(
  wearable_bloodflow_dm, bloodflow_vars, wearable_vars, "DR")
bloodflow_nodm_correlations <- calculate_all_correlations(
  wearable_bloodflow_nodm, bloodflow_vars, wearable_vars, "Cataract")

# 合并所有相关性
all_thickness_correlations <- rbind(thickness_dm_correlations, thickness_nodm_correlations)
all_bloodflow_correlations <- rbind(bloodflow_dm_correlations, bloodflow_nodm_correlations)


# 筛选显著的相关性 (p < 0.05)
significant_thickness_correlations <- all_thickness_correlations %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value)

significant_bloodflow_correlations <- all_bloodflow_correlations %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value)

# 打印显著相关性
print(significant_thickness_correlations)
print(significant_bloodflow_correlations)


# 保存相关性结果
write.csv(all_thickness_correlations, "thickness_correlations_pearson.csv", row.names = FALSE)
write.csv(all_bloodflow_correlations, "bloodflow_correlations_pearson.csv", row.names = FALSE)



# 为厚度和血流数据分别创建热图

# 1. 厚度相关性热图
# 处理数据以便于绘图
heatmap_thickness_data <- all_thickness_correlations %>%
  # 移除缺失值
  filter(!is.na(Correlation)) %>%
  # 保留需要的列
  dplyr::select(OCTA_Variable, Wearable_Variable, Correlation, P_value, Group) %>%
  # 改进变量名称以便于显示
  mutate(
    OCTA_Variable = gsub("Thickness_", "", OCTA_Variable),
    OCTA_Variable = gsub("_0_21_T0", "", OCTA_Variable),
    # 添加星号标记显著性
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**", 
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 为糖尿病和非糖尿病组分别创建热图
create_heatmap <- function(data, group_name, title) {
  group_data <- data %>% filter(Group == group_name)
  
  # 创建热图
  ggplot(group_data, aes(x = Wearable_Variable, y = OCTA_Variable, fill = Correlation)) +
    geom_tile(color = "white") +
    # 在每个单元格中添加相关系数值和显著性标记
    geom_text(aes(label = paste0(sprintf("%.2f", Correlation), Significance)), 
              color = "black", size = 3) +
    # 使用红蓝色阶，中间值为白色
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", 
      midpoint = 0, limit = c(-1, 1), space = "Lab", 
      name = "Pearson\nCorrelation"
    ) +
    # 美化图表
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    ) +
    labs(
      title = title,
      subtitle = paste("Group:", group_name),
      x = "Wearable Device Variables",
      y = "OCTA Thickness Variables"
    )
}

# 创建糖尿病组热图
thickness_dm_heatmap <- create_heatmap(
  heatmap_thickness_data, 
  "DR", 
  "Correlation between OCTA Thickness and Wearable Variables"
)

# 创建非糖尿病组热图
thickness_nodm_heatmap <- create_heatmap(
  heatmap_thickness_data, 
  "Cataract", 
  "Correlation between OCTA Thickness and Wearable Variables"
)


# 2. 血流相关性热图
# 处理数据以便于绘图
heatmap_bloodflow_data <- all_bloodflow_correlations %>%
  # 移除缺失值
  filter(!is.na(Correlation)) %>%
  # 保留需要的列
  dplyr::select(OCTA_Variable, Wearable_Variable, Correlation, P_value, Group) %>%
  # 改进变量名称以便于显示
  mutate(
    OCTA_Variable = gsub("PA_", "", OCTA_Variable),
    OCTA_Variable = gsub("_0_21_T0", "", OCTA_Variable),
    # 添加星号标记显著性
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**", 
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 创建糖尿病组血流热图
bloodflow_dm_heatmap <- create_heatmap(
  heatmap_bloodflow_data, 
  "DR", 
  "Correlation between OCTA Blood Flow and Wearable Variables"
)

# 创建非糖尿病组血流热图
bloodflow_nodm_heatmap <- create_heatmap(
  heatmap_bloodflow_data, 
  "Cataract", 
  "Correlation between OCTA Blood Flow and Wearable Variables"
)

# 显示创建的热图
print(thickness_dm_heatmap)
print(thickness_nodm_heatmap)
print(bloodflow_dm_heatmap)
print(bloodflow_nodm_heatmap)

# 组合并保存所有热图
# 厚度热图组合
ggsave("thickness_correlation_heatmap_DR.pdf", 
       thickness_dm_heatmap, width = 10, height = 12, dpi = 300)

ggsave("thickness_correlation_heatmap_nonDR.pdf", 
       thickness_nodm_heatmap, width = 10, height = 12, dpi = 300)

# 血流热图组合
ggsave("bloodflow_correlation_heatmap_DR.pdf", 
       bloodflow_dm_heatmap, width = 10, height = 12, dpi = 300)

ggsave("bloodflow_correlation_heatmap_nonDR.pdf", 
       bloodflow_nodm_heatmap, width = 10, height = 12, dpi = 300)






# 第二步：对显著相关性进行多元线性回归分析（控制年龄、性别和BMI）

# 创建回归分析函数
run_linear_regression <- function(data, octa_var, wearable_var, covariates) {
  # 创建回归公式
  formula_str <- paste(octa_var, "~", wearable_var, "+", paste(covariates, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  # 运行回归模型
  model <- lm(formula_obj, data = data)
  
  # 获取回归摘要
  model_summary <- summary(model)
  
  # 提取关注的可穿戴设备变量的系数信息
  coef_index <- which(names(coef(model)) == wearable_var)
  if(length(coef_index) > 0) {
    coefficient <- coef(model)[coef_index]
    std_error <- model_summary$coefficients[coef_index, 2]
    t_value <- model_summary$coefficients[coef_index, 3]
    p_value <- model_summary$coefficients[coef_index, 4]
    
    # 计算95%置信区间
    ci <- confint(model, wearable_var, level = 0.95)
    ci_lower <- ci[1]
    ci_upper <- ci[2]
    
    # 模型整体评估
    r_squared <- model_summary$r.squared
    adj_r_squared <- model_summary$adj.r.squared
    f_statistic <- model_summary$fstatistic[1]
    model_p_value <- pf(
      model_summary$fstatistic[1], 
      model_summary$fstatistic[2], 
      model_summary$fstatistic[3], 
      lower.tail = FALSE
    )
    
    return(data.frame(
      OCTA_Variable = octa_var,
      Wearable_Variable = wearable_var,
      Coefficient = coefficient,
      Std_Error = std_error,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      T_value = t_value,
      P_value = p_value,
      R_squared = r_squared,
      Adj_R_squared = adj_r_squared,
      F_statistic = f_statistic,
      Model_P_value = model_p_value
    ))
  } else {
    return(NULL)
  }
}

# 对显著的相关性进行回归分析
regression_results_thickness <- data.frame()
regression_results_bloodflow <- data.frame()

# 处理厚度变量的回归分析
if(nrow(significant_thickness_correlations) > 0) {
  for(i in 1:nrow(significant_thickness_correlations)) {
    row <- significant_thickness_correlations[i, ]
    
    # 选择适当的数据集
    dataset <- if(row$Group == "DR") wearable_thickness_dm else wearable_thickness_nodm
    
    # 运行回归分析
    reg_result <- run_linear_regression(
      dataset, 
      row$OCTA_Variable, 
      row$Wearable_Variable, 
      covariate_vars
    )
    
    if(!is.null(reg_result)) {
      reg_result$Group <- row$Group
      reg_result$Original_Correlation <- row$Correlation
      reg_result$Original_P_value <- row$P_value
      regression_results_thickness <- rbind(regression_results_thickness, reg_result)
    }
  }
}

# 处理血流变量的回归分析
if(nrow(significant_bloodflow_correlations) > 0) {
  for(i in 1:nrow(significant_bloodflow_correlations)) {
    row <- significant_bloodflow_correlations[i, ]
    
    # 选择适当的数据集
    dataset <- if(row$Group == "DR") wearable_bloodflow_dm else wearable_bloodflow_nodm
    
    # 运行回归分析
    reg_result <- run_linear_regression(
      dataset, 
      row$OCTA_Variable, 
      row$Wearable_Variable, 
      covariate_vars
    )
    
    if(!is.null(reg_result)) {
      reg_result$Group <- row$Group
      reg_result$Original_Correlation <- row$Correlation
      reg_result$Original_P_value <- row$P_value
      regression_results_bloodflow <- rbind(regression_results_bloodflow, reg_result)
    }
  }
}


# 筛选回归后仍显著的结果
significant_reg_thickness <- regression_results_thickness %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value)

significant_reg_bloodflow <- regression_results_bloodflow %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value)

# 打印回归后仍显著的结果
print(significant_reg_thickness)
print(significant_reg_bloodflow)


# 保存回归结果
write.csv(regression_results_thickness, "thickness_regression_results.csv", row.names = FALSE)
write.csv(regression_results_bloodflow, "bloodflow_regression_results.csv", row.names = FALSE)

# 提供一些图形可视化
# 创建散点图函数
plot_significant_association <- function(data, octa_var, wearable_var, group) {
  # 选择适当的数据集
  dataset <- if(group == "DR") {
    if(octa_var %in% thickness_vars) wearable_thickness_dm else wearable_bloodflow_dm
  } else {
    if(octa_var %in% thickness_vars) wearable_thickness_nodm else wearable_bloodflow_nodm
  }
  
  # 计算相关系数和p值
  # 首先过滤掉NA值
  complete_data <- dataset %>%
    dplyr::select(all_of(c(octa_var, wearable_var))) %>%
    na.omit()
  
  # 计算相关性
  cor_test <- cor.test(complete_data[[octa_var]], complete_data[[wearable_var]], 
                       method = "pearson", exact = FALSE)
  r_value <- cor_test$estimate
  p_value <- cor_test$p.value
  
  # 创建显示文本
  if(p_value < 0.001) {
    p_text <- "p < 0.001"
  } else {
    p_text <- paste0("p = ", sprintf("%.3f", p_value))
  }
  stat_text <- paste0("r = ", sprintf("%.2f", r_value), "\n", p_text)
  
  # 创建散点图
  p <- ggplot(dataset, aes_string(x = wearable_var, y = octa_var)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    # 添加相关系数和p值
    annotate("text", 
             x = max(dataset[[wearable_var]], na.rm = TRUE) - 
               (max(dataset[[wearable_var]], na.rm = TRUE) - 
                  min(dataset[[wearable_var]], na.rm = TRUE)) * 0.3,
             y = max(dataset[[octa_var]], na.rm = TRUE),
             label = stat_text, 
             hjust = 0, vjust = 1, 
             size = 5, 
             fontface = "bold") +
    theme_bw() +
    labs(
      title = paste("Association between", wearable_var, "and", octa_var),
      subtitle = paste("Group:", group),
      x = wearable_var,
      y = octa_var
    )
  
  # 保存图形
  ggsave(
    paste0("plot_", group, "_", octa_var, "_", wearable_var, ".pdf"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

# 为显著的关联创建图形
plots_thickness <- list()
if(nrow(significant_reg_thickness) > 0) {
  for(i in 1:nrow(significant_reg_thickness)) {
    row <- significant_reg_thickness[i, ]
    plots_thickness[[i]] <- plot_significant_association(
      data = NULL,  # 在函数内选择数据集
      octa_var = row$OCTA_Variable,
      wearable_var = row$Wearable_Variable,
      group = row$Group
    )
  }
}

plots_bloodflow <- list()
if(nrow(significant_reg_bloodflow) > 0) {
  for(i in 1:nrow(significant_reg_bloodflow)) {
    row <- significant_reg_bloodflow[i, ]
    plots_bloodflow[[i]] <- plot_significant_association(
      data = NULL,  # 在函数内选择数据集
      octa_var = row$OCTA_Variable,
      wearable_var = row$Wearable_Variable,
      group = row$Group
    )
  }
}




##########全人群analysis
# OCTA与可穿戴设备数据的两步分析流程
# 第一步：未校正的pearson相关分析
# 第二步：对显著关联进行多元线性回归分析（控制年龄、性别和BMI）

# 定义可穿戴变量和协变量用于分析
wearable_vars <- c("mean_rhr", "mean_bo", "total_steps", "total_sleep")
covariate_vars <- c("age", "gender", "bmi")

# 第4步：合并可穿戴设备数据与OCTA数据 (对全部人群)
wearable_thickness_data <- data_imputed %>%
  left_join(thickness_var_T0, by = c("subject_id" = "ID"))

wearable_bloodflow_data <- data_imputed %>%
  left_join(bloodflow_var_T0, by = c("subject_id" = "ID"))

# 第一步：创建未校正的pearson相关分析函数
calculate_pearson_correlation <- function(data, octa_var, wearable_var) {
  # 移除含有NA值的行
  complete_data <- data %>%
    dplyr::select(all_of(c(octa_var, wearable_var))) %>%
    na.omit()
  
  if(nrow(complete_data) < 5) {  # 至少需要5个观测值
    return(data.frame(
      OCTA_Variable = octa_var,
      Wearable_Variable = wearable_var,
      Correlation = NA,
      P_value = NA,
      N = nrow(complete_data)
    ))
  }
  
  # pearson相关系数
  pearson_test <- cor.test(complete_data[[octa_var]], complete_data[[wearable_var]], 
                           method = "pearson", exact = FALSE)
  
  correlation <- pearson_test$estimate
  p_value <- pearson_test$p.value
  n <- nrow(complete_data)
  
  return(data.frame(
    OCTA_Variable = octa_var,
    Wearable_Variable = wearable_var,
    Correlation = correlation,
    P_value = p_value,
    N = n
  ))
}

# 计算一组OCTA变量未校正相关性的函数
calculate_all_correlations <- function(data, octa_vars, wearable_vars) {
  all_correlations <- data.frame()
  
  for(octa_var in octa_vars) {
    for(wearable_var in wearable_vars) {
      correlation_result <- calculate_pearson_correlation(data, octa_var, wearable_var)
      all_correlations <- rbind(all_correlations, correlation_result)
    }
  }
  
  return(all_correlations)
}

# 计算厚度变量的相关性（未校正）- 全部人群
thickness_correlations <- calculate_all_correlations(
  wearable_thickness_data, thickness_vars, wearable_vars)

# 计算血流变量的相关性（未校正）- 全部人群
bloodflow_correlations <- calculate_all_correlations(
  wearable_bloodflow_data, bloodflow_vars, wearable_vars)


# 筛选显著的相关性 (p < 0.05)
significant_thickness_correlations <- thickness_correlations %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value)

significant_bloodflow_correlations <- bloodflow_correlations %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value)

# 打印显著相关性
print(significant_thickness_correlations)
print(significant_bloodflow_correlations)

# 保存相关性结果
write.csv(thickness_correlations, "thickness_correlations_pearson_all.csv", row.names = FALSE)
write.csv(bloodflow_correlations, "bloodflow_correlations_pearson_all.csv", row.names = FALSE)


# 为厚度和血流数据分别创建热图

# 1. 厚度相关性热图
# 处理数据以便于绘图
heatmap_thickness_data <- thickness_correlations %>%
  # 移除缺失值
  filter(!is.na(Correlation)) %>%
  # 保留需要的列
  dplyr::select(OCTA_Variable, Wearable_Variable, Correlation, P_value) %>%
  # 改进变量名称以便于显示
  mutate(
    OCTA_Variable = gsub("Thickness_", "", OCTA_Variable),
    OCTA_Variable = gsub("_0_21_T0", "", OCTA_Variable),
    # 添加星号标记显著性
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**", 
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 创建热图函数
create_heatmap <- function(data, title) {
  # 创建热图
  ggplot(data, aes(x = Wearable_Variable, y = OCTA_Variable, fill = Correlation)) +
    geom_tile(color = "white") +
    # 在每个单元格中添加相关系数值和显著性标记
    geom_text(aes(label = paste0(sprintf("%.2f", Correlation), Significance)), 
              color = "black", size = 3) +
    # 使用红蓝色阶，中间值为白色
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", 
      midpoint = 0, limit = c(-1, 1), space = "Lab", 
      name = "Pearson\nCorrelation"
    ) +
    # 美化图表
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    ) +
    labs(
      title = title,
      subtitle = "All Participants",
      x = "Wearable Device Variables",
      y = "OCTA Variables"
    )
}

# 创建厚度热图
thickness_heatmap <- create_heatmap(
  heatmap_thickness_data, 
  "Correlation between OCTA Thickness and Wearable Variables"
)

# 2. 血流相关性热图
# 处理数据以便于绘图
heatmap_bloodflow_data <- bloodflow_correlations %>%
  # 移除缺失值
  filter(!is.na(Correlation)) %>%
  # 保留需要的列
  dplyr::select(OCTA_Variable, Wearable_Variable, Correlation, P_value) %>%
  # 改进变量名称以便于显示
  mutate(
    OCTA_Variable = gsub("PA_", "", OCTA_Variable),
    OCTA_Variable = gsub("_0_21_T0", "", OCTA_Variable),
    # 添加星号标记显著性
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**", 
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 创建血流热图
bloodflow_heatmap <- create_heatmap(
  heatmap_bloodflow_data, 
  "Correlation between OCTA Blood Flow and Wearable Variables"
)

# 显示创建的热图
print(thickness_heatmap)
print(bloodflow_heatmap)

# 保存热图
ggsave("thickness_correlation_heatmap_all.pdf", 
       thickness_heatmap, width = 10, height = 12, dpi = 300)

ggsave("bloodflow_correlation_heatmap_all.pdf", 
       bloodflow_heatmap, width = 10, height = 12, dpi = 300)

# 第二步：对显著相关性进行多元线性回归分析（控制年龄、性别和BMI）

# 创建回归分析函数
run_linear_regression <- function(data, octa_var, wearable_var, covariates) {
  # 创建回归公式
  formula_str <- paste(octa_var, "~", wearable_var, "+ dm_2 +", paste(covariates, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  # 运行回归模型
  model <- lm(formula_obj, data = data)
  
  # 获取回归摘要
  model_summary <- summary(model)
  
  # 提取关注的可穿戴设备变量的系数信息
  coef_index <- which(names(coef(model)) == wearable_var)
  if(length(coef_index) > 0) {
    coefficient <- coef(model)[coef_index]
    std_error <- model_summary$coefficients[coef_index, 2]
    t_value <- model_summary$coefficients[coef_index, 3]
    p_value <- model_summary$coefficients[coef_index, 4]
    
    # 计算95%置信区间
    ci <- confint(model, wearable_var, level = 0.95)
    ci_lower <- ci[1]
    ci_upper <- ci[2]
    
    # 模型整体评估
    r_squared <- model_summary$r.squared
    adj_r_squared <- model_summary$adj.r.squared
    f_statistic <- model_summary$fstatistic[1]
    model_p_value <- pf(
      model_summary$fstatistic[1], 
      model_summary$fstatistic[2], 
      model_summary$fstatistic[3], 
      lower.tail = FALSE
    )
    
    return(data.frame(
      OCTA_Variable = octa_var,
      Wearable_Variable = wearable_var,
      Coefficient = coefficient,
      Std_Error = std_error,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      T_value = t_value,
      P_value = p_value,
      R_squared = r_squared,
      Adj_R_squared = adj_r_squared,
      F_statistic = f_statistic,
      Model_P_value = model_p_value
    ))
  } else {
    return(NULL)
  }
}

# 对显著的相关性进行回归分析
regression_results_thickness <- data.frame()
regression_results_bloodflow <- data.frame()

# 处理厚度变量的回归分析
if(nrow(significant_thickness_correlations) > 0) {
  for(i in 1:nrow(significant_thickness_correlations)) {
    row <- significant_thickness_correlations[i, ]
    
    # 运行回归分析
    reg_result <- run_linear_regression(
      wearable_thickness_data, 
      row$OCTA_Variable, 
      row$Wearable_Variable, 
      covariate_vars
    )
    
    if(!is.null(reg_result)) {
      reg_result$Original_Correlation <- row$Correlation
      reg_result$Original_P_value <- row$P_value
      regression_results_thickness <- rbind(regression_results_thickness, reg_result)
    }
  }
}

# 处理血流变量的回归分析
if(nrow(significant_bloodflow_correlations) > 0) {
  for(i in 1:nrow(significant_bloodflow_correlations)) {
    row <- significant_bloodflow_correlations[i, ]
    
    # 运行回归分析
    reg_result <- run_linear_regression(
      wearable_bloodflow_data, 
      row$OCTA_Variable, 
      row$Wearable_Variable, 
      covariate_vars
    )
    
    if(!is.null(reg_result)) {
      reg_result$Original_Correlation <- row$Correlation
      reg_result$Original_P_value <- row$P_value
      regression_results_bloodflow <- rbind(regression_results_bloodflow, reg_result)
    }
  }
}

# 筛选回归后仍显著的结果
significant_reg_thickness <- regression_results_thickness %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value)

significant_reg_bloodflow <- regression_results_bloodflow %>% 
  filter(P_value < 0.05) %>% 
  arrange(P_value)

# 打印回归后仍显著的结果
print(significant_reg_thickness)
print(significant_reg_bloodflow)

# 保存回归结果
write.csv(regression_results_thickness, "thickness_regression_results_all.csv", row.names = FALSE)
write.csv(regression_results_bloodflow, "bloodflow_regression_results_all.csv", row.names = FALSE)



# 提供一些图形可视化
# 创建散点图函数
plot_significant_association <- function(data, octa_var, wearable_var) {
  # 获取x轴和y轴数据，以便更好地调整坐标轴范围
  x_data <- data[[wearable_var]]
  y_data <- data[[octa_var]]
  
  # 过滤出配对完整的数据（同时有x和y的值）
  complete_idx <- !is.na(x_data) & !is.na(y_data)
  x_complete <- x_data[complete_idx]
  y_complete <- y_data[complete_idx]
  
  # 当数据点太少时处理
  if(length(x_complete) < 5) {
    warning(paste("数据点太少，无法创建图表:", octa_var, "vs", wearable_var))
    return(NULL)
  }
  
  # 计算相关系数和p值
  cor_test <- cor.test(x_complete, y_complete, method = "pearson")
  r_value <- cor_test$estimate
  p_value <- cor_test$p.value
  
  # 创建p值的显示文本
  if(p_value < 0.001) {
    p_text <- "p < 0.001"
  } else if(p_value < 0.01) {
    p_text <- paste0("p = ", sprintf("%.3f", p_value))
  } else {
    p_text <- paste0("p = ", sprintf("%.2f", p_value))
  }
  
  # 创建注释文本
  stat_text <- paste0("r = ", sprintf("%.2f", r_value), "\n", p_text)
  
  # 计算坐标轴范围
  x_min <- min(x_complete, na.rm = TRUE)
  x_max <- max(x_complete, na.rm = TRUE)
  y_min <- min(y_complete, na.rm = TRUE)
  y_max <- max(y_complete, na.rm = TRUE)
  
  # 计算适当的坐标轴范围，避免多余空间
  x_range <- x_max - x_min
  y_range <- y_max - y_min
  
  x_min_plot <- x_min - x_range * 0.05  # 增加一点边距避免裁剪
  x_max_plot <- x_max + x_range * 0.05
  
  if(grepl("PA_", octa_var)) {
    # 对于血流数据，通常值域为0-1之间，适当调整
    y_min_plot <- max(0, y_min * 0.95)
    y_max_plot <- min(1.5, y_max * 1.1)  # 增加上限以容纳注释文本
  } else {
    # 对于厚度数据，使用数据范围的比例调整
    y_min_plot <- y_min - y_range * 0.05
    y_max_plot <- y_max + y_range * 0.1  # 增加上限以容纳注释文本
  }
  
  # 计算统计信息的位置 - 放在图的右上角，避免与数据点和回归线重叠
  text_x <- x_min + x_range * 0.7
  text_y <- y_max - y_range * 0.1
  
  # 设置更好的变量名显示
  y_label <- octa_var
  if(grepl("PA_", octa_var)) {
    y_label <- gsub("PA_", "", y_label)
    y_label <- gsub("_0_21_T0", "", y_label)
  } else if(grepl("Thickness_", octa_var)) {
    y_label <- gsub("Thickness_", "", y_label)
    y_label <- gsub("_0_21_T0", "", y_label)
  }
  
  x_label <- wearable_var
  if(x_label == "mean_rhr") x_label <- "Mean Resting Heart Rate"
  if(x_label == "mean_bo") x_label <- "Mean Blood Oxygen"
  if(x_label == "total_steps") x_label <- "Total Steps"
  if(x_label == "total_sleep") x_label <- "Total Sleep (min)"
  
  # 创建散点图，使用过滤后的数据框
  filtered_data <- data[complete_idx, ]
  
  p <- ggplot(filtered_data, aes_string(x = wearable_var, y = octa_var)) +
    geom_point(size = 2) +  # 增加点的大小
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    # 添加相关系数和p值文本
    annotate("text", x = text_x, y = text_y, 
             label = stat_text, 
             hjust = 0, vjust = 1, 
             size = 5, 
             fontface = "bold") +
    theme_bw() +
    # 设置x轴和y轴范围以避免多余空间
    coord_cartesian(xlim = c(x_min_plot, x_max_plot), 
                    ylim = c(y_min_plot, y_max_plot)) +  # 使用coord_cartesian代替scale_x/y_continuous
    labs(
      title = paste("Association between", x_label, "and", y_label),
      subtitle = "All Participants",
      x = x_label,
      y = y_label
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()  # 移除次要网格线以减少视觉干扰
    )
  
  # 保存图形
  ggsave(
    paste0("plot_all_", octa_var, "_", wearable_var, ".pdf"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

# 为显著的关联创建图形
plots_thickness <- list()
if(nrow(significant_reg_thickness) > 0) {
  for(i in 1:nrow(significant_reg_thickness)) {
    row <- significant_reg_thickness[i, ]
    plots_thickness[[i]] <- plot_significant_association(
      data = wearable_thickness_data,
      octa_var = row$OCTA_Variable,
      wearable_var = row$Wearable_Variable
    )
  }
}

plots_bloodflow <- list()
if(nrow(significant_reg_bloodflow) > 0) {
  for(i in 1:nrow(significant_reg_bloodflow)) {
    row <- significant_reg_bloodflow[i, ]
    plots_bloodflow[[i]] <- plot_significant_association(
      data = wearable_bloodflow_data,
      octa_var = row$OCTA_Variable,
      wearable_var = row$Wearable_Variable
    )
  }
}




##################################
# 白内障与糖尿病患者OCTA参数比较 #
##################################
#############################################
# 厚度参数比较 #
#############################################

# 创建比较函数
compare_groups <- function(data1, data2, vars, group1_name, group2_name) {
  results <- data.frame()
  
  for(var in vars) {
    # 提取两组数据
    var_data1 <- data1[[var]]
    var_data1 <- var_data1[!is.na(var_data1)]
    
    var_data2 <- data2[[var]]
    var_data2 <- var_data2[!is.na(var_data2)]
    
    # 计算基本统计量
    n1 <- length(var_data1)
    n2 <- length(var_data2)
    mean1 <- mean(var_data1, na.rm = TRUE)
    mean2 <- mean(var_data2, na.rm = TRUE)
    sd1 <- sd(var_data1, na.rm = TRUE)
    sd2 <- sd(var_data2, na.rm = TRUE)
    
    # 执行统计检验
    if(n1 < 3 || n2 < 3) {
      p_value <- NA
      test_type <- "Insufficient data"
    } else {
      # 检查正态性
      shapiro_test1 <- try(shapiro.test(var_data1), silent = TRUE)
      shapiro_test2 <- try(shapiro.test(var_data2), silent = TRUE)
      
      if(!inherits(shapiro_test1, "try-error") && !inherits(shapiro_test2, "try-error") &&
         shapiro_test1$p.value > 0.05 && shapiro_test2$p.value > 0.05) {
        # 正态分布用t检验
        t_test_result <- t.test(var_data1, var_data2, var.equal = FALSE)
        p_value <- t_test_result$p.value
        test_type <- "t-test"
      } else {
        # 非正态分布用Wilcoxon检验
        wilcox_test_result <- try(wilcox.test(var_data1, var_data2, exact = FALSE), silent = TRUE)
        if(!inherits(wilcox_test_result, "try-error")) {
          p_value <- wilcox_test_result$p.value
          test_type <- "Wilcoxon"
        } else {
          p_value <- NA
          test_type <- "Test failed"
        }
      }
    }
    
    # 计算差异
    mean_diff <- mean1 - mean2
    percent_diff <- (mean_diff / mean2) * 100
    
    # 创建结果
    comparison <- data.frame(
      Variable = var,
      Group1 = group1_name,
      Group2 = group2_name,
      N1 = n1,
      N2 = n2,
      Mean1 = mean1,
      SD1 = sd1,
      Mean2 = mean2,
      SD2 = sd2,
      Mean_Difference = mean_diff,
      Percent_Difference = percent_diff,
      Test_Type = test_type,
      P_value = p_value
    )
    
    results <- rbind(results, comparison)
  }
  
  # 添加多重比较校正
  results$P_adjusted <- p.adjust(results$P_value, method = "BH")
  
  # 添加显著性标记
  results$Significance <- ""
  results$Significance[results$P_adjusted < 0.05] <- "*"
  results$Significance[results$P_adjusted < 0.01] <- "**"
  results$Significance[results$P_adjusted < 0.001] <- "***"
  
  return(results)
}

# 比较厚度参数
thickness_comparison <- compare_groups(
  wearable_thickness_nodm, wearable_thickness_dm, 
  thickness_vars, "Cataract", "DR"
)

# 比较血流参数
bloodflow_comparison <- compare_groups(
  wearable_bloodflow_nodm, wearable_bloodflow_dm, 
  bloodflow_vars, "Cataract", "DR"
)


## 创建美化的厚度比较表格
improved_thickness_comparison <- thickness_comparison %>%
  mutate(
    Variable_Display = gsub("Thickness_", "", Variable),
    Variable_Display = gsub("_0_21_T0", "", Variable_Display),
    Mean1_Display = sprintf("%.2f ± %.2f", Mean1, SD1),
    Mean2_Display = sprintf("%.2f ± %.2f", Mean2, SD2),
    Diff_Display = sprintf("%.2f (%.1f%%)", Mean_Difference, Percent_Difference),
    P_Display = case_when(
      P_adjusted < 0.001 ~ paste0("<0.001", Significance),
      TRUE ~ paste0(sprintf("%.3f", P_adjusted), Significance)
    )
  ) %>%
  dplyr::select(
    Variable_Display, N1, N2, Mean1_Display, Mean2_Display, 
    Diff_Display, Test_Type, P_Display
  )

# 单独执行排序
improved_thickness_comparison <- improved_thickness_comparison %>% 
  arrange(thickness_comparison$P_adjusted)


# 创建美化的血流比较表格
improved_bloodflow_comparison <- bloodflow_comparison %>%
  mutate(
    Variable_Display = gsub("PA_", "", Variable),
    Variable_Display = gsub("_0_21_T0", "", Variable_Display),
    Mean1_Display = sprintf("%.3f ± %.3f", Mean1, SD1),
    Mean2_Display = sprintf("%.3f ± %.3f", Mean2, SD2),
    Diff_Display = sprintf("%.3f (%.1f%%)", Mean_Difference, Percent_Difference),
    P_Display = case_when(
      P_adjusted < 0.001 ~ paste0("<0.001", Significance),
      TRUE ~ paste0(sprintf("%.3f", P_adjusted), Significance)
    )
  ) %>%
  dplyr::select(
    Variable_Display, N1, N2, Mean1_Display, Mean2_Display, 
    Diff_Display, Test_Type, P_Display
  )

# 单独执行排序
improved_bloodflow_comparison <- improved_bloodflow_comparison %>% 
  arrange(bloodflow_comparison$P_adjusted)

# 保存表格
write.csv(improved_thickness_comparison, "thickness_comparison_table.csv", row.names = FALSE)
write.csv(improved_bloodflow_comparison, "bloodflow_comparison_table.csv", row.names = FALSE)

# 输出显著差异的参数数量
cat("厚度参数中有显著差异的数量：", sum(thickness_comparison$P_adjusted < 0.05), "\n")
cat("血流参数中有显著差异的数量：", sum(bloodflow_comparison$P_adjusted < 0.05), "\n")

# 打印前5个最显著差异的厚度参数
cat("\n前5个最显著差异的厚度参数：\n")
print(improved_thickness_comparison %>% head(5) %>% dplyr::select(Variable_Display, Mean1_Display, Mean2_Display, P_Display))

# 打印前5个最显著差异的血流参数
cat("\n前5个最显著差异的血流参数：\n")
print(improved_bloodflow_comparison %>% head(10) %>% dplyr::select(Variable_Display, Mean1_Display, Mean2_Display, P_Display))

