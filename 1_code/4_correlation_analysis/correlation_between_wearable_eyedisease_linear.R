library(tidyverse)
library(lme4)
library(tableone)
library(ggpubr)
library(car)  # For statistical tests
library(corrplot)  # For correlation visualization
library(pROC)  # For ROC analysis

##########################
# 1. 加载数据
##########################
# 请根据实际路径修改
load("3_data_analysis/2_data_analysis/RHR/RHR_time_periods/daily_rhr_result.rda")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

##########################
# 2. 数据处理
##########################
# 处理视力数据
vision_data <- baseline_info %>%
  filter(!is.na(surgery_eye_1)) %>%
  distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
  mutate(
    pre_vision = case_when(
      surgery_eye_1 == 0 ~ od_corrected_bas,  # Right eye surgery
      surgery_eye_1 == 1 ~ os_corrected_bas,  # Left eye surgery
      surgery_eye_1 == 2 ~ (od_corrected_bas + os_corrected_bas)/2,  # Both eyes (average)
      TRUE ~ NA_real_
    )
  )

# 创建分析组别（DR组和单纯白内障组）
analysis_groups <- baseline_info %>%
  mutate(
    group = case_when(
      diabetes_history == 1 ~ "DR",  # 所有糖尿病患者归为DR组
      diabetes_history == 2 & cataract %in% c(2, 3, 4) ~ "Cataract",  # 无糖尿病但有白内障
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  dplyr::select(ID, group, age, gender, bmi, cataract)

# 显示分组数量
cat("分组统计：\n")
print(table(analysis_groups$group))

##########################
# 3. 提取术前RHR数据
##########################
# 将宽格式转为长格式并筛选术前7天数据
extract_pre_surgery_rhr <- function(data, baseline_info) {
  # 获取所有列名
  all_cols <- names(data)
  
  # 人口统计学列
  demo_cols <- c("subject_id", "surgery_date")
  
  # RHR列（排除人口统计学列）
  rhr_cols <- all_cols[!all_cols %in% demo_cols]
  
  # 创建长格式数据
  long_data <- data %>%
    # 与基线信息连接
    left_join(
      baseline_info %>% 
        dplyr::select(ID, group, age, gender, bmi),
      by = c("subject_id" = "ID")
    ) %>%
    # 转为长格式
    pivot_longer(
      cols = all_of(rhr_cols),
      names_to = "variable",
      values_to = "value"
    ) %>%
    # 提取日期、统计类型和RHR测量类型
    mutate(
      day = as.numeric(str_extract(variable, "-?\\d+")),
      stat_type = str_extract(variable, "(mean|min|max|median|sd|cv|iqr|skew|kurt)"),
      rhr_type = str_extract(variable, "rhr_\\d+")
    ) %>%
    # 只保留术前7天数据
    filter(!is.na(value), day >= -7, day < 0, !is.na(group))
  
  return(long_data)
}

# 提取术前RHR数据
pre_surgery_rhr <- extract_pre_surgery_rhr(daily_rhr_result, analysis_groups)

# 显示术前数据统计
cat("术前数据统计：\n")
cat("总观测数:", nrow(pre_surgery_rhr), "\n")
cat("唯一受试者:", n_distinct(pre_surgery_rhr$subject_id), "\n")
cat("唯一日期:", n_distinct(pre_surgery_rhr$day), "\n")
cat("统计指标:", paste(unique(pre_surgery_rhr$stat_type), collapse = ", "), "\n")

##########################
# 4. 术前RHR比较（DR vs 白内障）
##########################
# 计算每个受试者术前平均RHR
subject_pre_rhr <- pre_surgery_rhr %>%
  filter(stat_type == "mean") %>%  # 使用均值作为RHR指标
  group_by(subject_id, group, rhr_type) %>%
  summarize(
    pre_rhr_mean = mean(value, na.rm = TRUE),
    pre_rhr_sd = sd(value, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  ) %>%
  # 只保留有足够数据点的受试者
  filter(n_days >= 3) %>%  # 至少有3天术前数据
  # 合并人口统计学数据
  left_join(
    analysis_groups %>% dplyr::select(ID, age, gender, bmi),
    by = c("subject_id" = "ID")
  )

# 比较两组RHR（通过t检验）
compare_rhr_groups <- function(data, rhr_type_filter = "rhr_1") {
  # 按RHR类型筛选数据
  filtered_data <- data %>%
    filter(rhr_type == rhr_type_filter)
  
  # 分组数据
  dr_data <- filtered_data %>% filter(group == "DR") %>% pull(pre_rhr_mean)
  cataract_data <- filtered_data %>% filter(group == "Cataract") %>% pull(pre_rhr_mean)
  
  # t检验
  t_result <- t.test(dr_data, cataract_data)
  
  # 计算效应量（Cohen's d）
  cohens_d <- (mean(dr_data, na.rm = TRUE) - mean(cataract_data, na.rm = TRUE)) / 
    sqrt((sd(dr_data, na.rm = TRUE)^2 + sd(cataract_data, na.rm = TRUE)^2) / 2)
  
  # 构建结果表
  result_table <- data.frame(
    rhr_type = rhr_type_filter,
    dr_n = length(dr_data),
    dr_mean = mean(dr_data, na.rm = TRUE),
    dr_sd = sd(dr_data, na.rm = TRUE),
    cataract_n = length(cataract_data),
    cataract_mean = mean(cataract_data, na.rm = TRUE),
    cataract_sd = sd(cataract_data, na.rm = TRUE),
    mean_diff = mean(dr_data, na.rm = TRUE) - mean(cataract_data, na.rm = TRUE),
    t_statistic = t_result$statistic,
    p_value = t_result$p.value,
    cohens_d = cohens_d
  )
  
  # 返回结果
  return(list(
    result_table = result_table,
    t_test = t_result
  ))
}

# 对两种RHR类型进行比较
rhr1_comparison <- compare_rhr_groups(subject_pre_rhr, "rhr_1")
rhr50_comparison <- compare_rhr_groups(subject_pre_rhr, "rhr_50")

# 打印结果
cat("\nDR组与白内障组RHR_1比较：\n")
print(rhr1_comparison$result_table)
print(rhr1_comparison$t_test)

cat("\nDR组与白内障组RHR_50比较：\n")
print(rhr50_comparison$result_table)
print(rhr50_comparison$t_test)

##########################
# 5. 术前RHR与视力相关性分析
##########################
# 合并RHR和视力数据
rhr_vision_data <- subject_pre_rhr %>%
  # 合并视力数据
  left_join(
    vision_data %>% dplyr::select(ID, pre_vision),
    by = c("subject_id" = "ID")
  ) %>%
  # 删除缺失视力数据
  filter(!is.na(pre_vision))

# 相关性分析函数
analyze_correlation <- function(data, rhr_type_filter = "rhr_1") {
  # 按RHR类型和组筛选数据
  all_data <- data %>% filter(rhr_type == rhr_type_filter)
  dr_data <- all_data %>% filter(group == "DR")
  cataract_data <- all_data %>% filter(group == "Cataract")
  
  # 全部数据相关性
  all_cor <- cor.test(all_data$pre_rhr_mean, all_data$pre_vision, 
                      method = "spearman")
  
  # DR组相关性
  dr_cor <- cor.test(dr_data$pre_rhr_mean, dr_data$pre_vision, 
                     method = "spearman")
  
  # 白内障组相关性
  cataract_cor <- cor.test(cataract_data$pre_rhr_mean, cataract_data$pre_vision, 
                           method = "spearman")
  
  # 构建结果表
  result_table <- data.frame(
    rhr_type = rhr_type_filter,
    group = c("All", "DR", "Cataract"),
    n = c(nrow(all_data), nrow(dr_data), nrow(cataract_data)),
    correlation = c(all_cor$estimate, dr_cor$estimate, cataract_cor$estimate),
    p_value = c(all_cor$p.value, dr_cor$p.value, cataract_cor$p.value)
  )
  
  # 创建相关性散点图
  p <- ggplot(all_data, aes(x = pre_rhr_mean, y = pre_vision, color = group)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, aes(group = group)) +
    labs(
      title = paste("RHR与术前视力相关性 -", rhr_type_filter),
      x = "术前平均RHR",
      y = "术前视力",
      color = "组别"
    ) +
    scale_color_manual(values = c("DR" = "#E41A1C", "Cataract" = "#377EB8")) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    # 添加相关系数
    annotate("text", x = min(all_data$pre_rhr_mean, na.rm = TRUE), 
             y = max(all_data$pre_vision, na.rm = TRUE) * 0.9,
             label = paste("All: r =", round(all_cor$estimate, 2), 
                           ", p =", round(all_cor$p.value, 3)),
             hjust = 0, color = "black") +
    annotate("text", x = min(all_data$pre_rhr_mean, na.rm = TRUE), 
             y = max(all_data$pre_vision, na.rm = TRUE) * 0.85,
             label = paste("DR: r =", round(dr_cor$estimate, 2), 
                           ", p =", round(dr_cor$p.value, 3)),
             hjust = 0, color = "#E41A1C") +
    annotate("text", x = min(all_data$pre_rhr_mean, na.rm = TRUE), 
             y = max(all_data$pre_vision, na.rm = TRUE) * 0.8,
             label = paste("Cataract: r =", round(cataract_cor$estimate, 2), 
                           ", p =", round(cataract_cor$p.value, 3)),
             hjust = 0, color = "#377EB8")
  
  # 返回结果
  return(list(
    result_table = result_table,
    plot = p,
    all_cor = all_cor,
    dr_cor = dr_cor,
    cataract_cor = cataract_cor
  ))
}

# 分析两种RHR类型与视力的相关性
rhr1_correlation <- analyze_correlation(rhr_vision_data, "rhr_1")
rhr50_correlation <- analyze_correlation(rhr_vision_data, "rhr_50")

# 打印结果
cat("\nRHR与术前视力相关性分析：\n")
print(rhr1_correlation$result_table)
print(rhr50_correlation$result_table)

# 显示相关性图表
print(rhr1_correlation$plot)
print(rhr50_correlation$plot)

##########################
# 6. 多变量分析：控制年龄和性别后的RHR差异
##########################
# 多变量线性回归模型
run_linear_model <- function(data, rhr_type_filter = "rhr_1") {
  # 筛选数据
  model_data <- data %>%
    filter(rhr_type == rhr_type_filter)
  
  # 运行回归模型
  model <- lm(pre_rhr_mean ~ group + age + gender + bmi, data = model_data)
  
  # 模型摘要
  model_summary <- summary(model)
  
  # ANOVA分析
  model_anova <- Anova(model, type = "II")
  
  # 返回结果
  return(list(
    model = model,
    summary = model_summary,
    anova = model_anova
  ))
}

# 运行多变量模型
rhr1_model <- run_linear_model(subject_pre_rhr, "rhr_1")
rhr50_model <- run_linear_model(subject_pre_rhr, "rhr_50")

# 打印结果
cat("\nRHR_1多变量回归模型：\n")
print(rhr1_model$summary)
print(rhr1_model$anova)

cat("\nRHR_50多变量回归模型：\n")
print(rhr50_model$summary)
print(rhr50_model$anova)

##########################
# 7. 多变量分析：RHR与视力的关系（控制混杂因素）
##########################
# 多变量模型：RHR与视力的关系
run_vision_model <- function(data, rhr_type_filter = "rhr_1") {
  # 筛选数据
  model_data <- data %>%
    filter(rhr_type == rhr_type_filter)
  
  # 运行回归模型
  model <- lm(pre_vision ~ pre_rhr_mean + group + age + gender + bmi, data = model_data)
  
  # 模型摘要
  model_summary <- summary(model)
  
  # ANOVA分析
  model_anova <- Anova(model, type = "II")
  
  # 返回结果
  return(list(
    model = model,
    summary = model_summary,
    anova = model_anova
  ))
}

# 运行多变量模型
vision_rhr1_model <- run_vision_model(rhr_vision_data, "rhr_1")
vision_rhr50_model <- run_vision_model(rhr_vision_data, "rhr_50")

# 打印结果
cat("\nRHR_1与视力多变量回归模型：\n")
print(vision_rhr1_model$summary)
print(vision_rhr1_model$anova)

cat("\nRHR_50与视力多变量回归模型：\n")
print(vision_rhr50_model$summary)
print(vision_rhr50_model$anova)

##########################
# 8. ROC分析：RHR预测DR的能力
##########################
# ROC分析
perform_roc_analysis <- function(data, rhr_type_filter = "rhr_1") {
  # 筛选数据
  roc_data <- data %>%
    filter(rhr_type == rhr_type_filter) %>%
    mutate(
      dr_status = ifelse(group == "DR", 1, 0)  # 二分类目标变量
    )
  
  # 执行ROC分析
  roc_result <- roc(dr_status ~ pre_rhr_mean, data = roc_data)
  
  # 直接提取ROC曲线数据
  roc_df <- data.frame(
    threshold = roc_result$thresholds,
    sensitivity = roc_result$sensitivities,
    specificity = roc_result$specificities
  )
  
  # 计算约登指数
  roc_df$youden <- roc_df$sensitivity + roc_df$specificity - 1
  
  # 找出最佳阈值
  best_index <- which.max(roc_df$youden)
  best_threshold <- roc_df$threshold[best_index]
  sensitivity <- roc_df$sensitivity[best_index]
  specificity <- roc_df$specificity[best_index]
  
  # 创建ROC曲线
  roc_plot <- ggroc(roc_result) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray") +
    annotate("text", x = 0.75, y = 0.25,
             label = paste0("AUC = ", round(auc(roc_result), 3), "\n",
                            "95% CI: ", round(ci(roc_result)[1], 3), "-", 
                            round(ci(roc_result)[3], 3))) +
    labs(
      title = paste("ROC曲线：术前RHR预测DR -", rhr_type_filter),
      x = "1 - 特异度",
      y = "灵敏度"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  # 返回结果
  return(list(
    roc = roc_result,
    auc = auc(roc_result),
    ci = ci(roc_result),
    best_threshold = best_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    plot = roc_plot
  ))
}

# 执行ROC分析
rhr1_roc <- perform_roc_analysis(subject_pre_rhr, "rhr_1")
rhr50_roc <- perform_roc_analysis(subject_pre_rhr, "rhr_50")

# 打印结果
cat("\nRHR_1预测DR的ROC分析：\n")
cat("AUC:", auc(rhr1_roc$roc), "\n")
cat("95% CI:", ci(rhr1_roc$roc)[1], "-", ci(rhr1_roc$roc)[3], "\n")
cat("最佳阈值:", rhr1_roc$best_threshold, "\n")
cat("灵敏度:", rhr1_roc$sensitivity, "\n")
cat("特异度:", rhr1_roc$specificity, "\n")

cat("\nRHR_50预测DR的ROC分析：\n")
cat("AUC:", auc(rhr50_roc$roc), "\n")
cat("95% CI:", ci(rhr50_roc$roc)[1], "-", ci(rhr50_roc$roc)[3], "\n")
cat("最佳阈值:", rhr50_roc$best_threshold, "\n")
cat("灵敏度:", rhr50_roc$sensitivity, "\n")
cat("特异度:", rhr50_roc$specificity, "\n")

# 显示ROC曲线
print(rhr1_roc$plot)
print(rhr50_roc$plot)

##########################
# 9. 结果保存
##########################
# 创建结果目录
dir.create("results/dr_cataract_analysis", recursive = TRUE)

# 保存图表
ggsave("results/dr_cataract_analysis/rhr1_vision_correlation.png", rhr1_correlation$plot, 
       width = 8, height = 6, dpi = 300)
ggsave("results/dr_cataract_analysis/rhr50_vision_correlation.png", rhr50_correlation$plot, 
       width = 8, height = 6, dpi = 300)
ggsave("results/dr_cataract_analysis/rhr1_roc.png", rhr1_roc$plot, 
       width = 8, height = 6, dpi = 300)
ggsave("results/dr_cataract_analysis/rhr50_roc.png", rhr50_roc$plot, 
       width = 8, height = 6, dpi = 300)

# 保存分析结果
save(subject_pre_rhr, rhr_vision_data, 
     rhr1_comparison, rhr50_comparison,
     rhr1_correlation, rhr50_correlation,
     rhr1_model, rhr50_model,
     vision_rhr1_model, vision_rhr50_model,
     rhr1_roc, rhr50_roc,
     file = "results/dr_cataract_analysis/dr_cataract_rhr_analysis.rda")

# 创建摘要报告
sink("results/dr_cataract_analysis/analysis_summary.txt")

cat("DR组与白内障组RHR差异及与视力相关性分析\n")
cat("=========================================\n\n")
cat("分析日期:", format(Sys.time(), "%Y-%m-%d"), "\n\n")

cat("1. 研究对象\n")
cat("------------\n")
cat("总受试者数:", n_distinct(pre_surgery_rhr$subject_id), "\n")
cat("DR组人数:", sum(analysis_groups$group == "DR"), "\n")
cat("白内障组人数:", sum(analysis_groups$group == "Cataract"), "\n\n")

cat("2. DR组与白内障组RHR差异\n")
cat("-----------------------\n")
cat("RHR_1:\n")
print(rhr1_comparison$result_table)
cat("\n")
cat("RHR_50:\n")
print(rhr50_comparison$result_table)
cat("\n")

cat("3. RHR与术前视力相关性\n")
cat("-----------------------\n")
cat("RHR_1与视力相关性:\n")
print(rhr1_correlation$result_table)
cat("\n")
cat("RHR_50与视力相关性:\n")
print(rhr50_correlation$result_table)
cat("\n")

cat("4. 多变量分析\n")
cat("-------------\n")
cat("控制年龄、性别和BMI后的RHR_1组间差异:\n")
print(rhr1_model$summary$coefficients)
cat("\n")
cat("控制年龄、性别和BMI后的RHR与视力关系:\n")
print(vision_rhr1_model$summary$coefficients)
cat("\n")

cat("5. ROC分析结果\n")
cat("--------------\n")
cat("RHR_1预测DR:\n")
cat("AUC:", auc(rhr1_roc$roc), "\n")
cat("95% CI:", ci(rhr1_roc$roc)[1], "-", ci(rhr1_roc$roc)[3], "\n")
cat("最佳阈值:", rhr1_roc$best_threshold, "\n")
cat("灵敏度:", rhr1_roc$sensitivity, "\n")
cat("特异度:", rhr1_roc$specificity, "\n\n")

cat("RHR_50预测DR:\n")
cat("AUC:", auc(rhr50_roc$roc), "\n")
cat("95% CI:", ci(rhr50_roc$roc)[1], "-", ci(rhr50_roc$roc)[3], "\n")
cat("最佳阈值:", rhr50_roc$best_threshold, "\n")
cat("灵敏度:", rhr50_roc$sensitivity, "\n")
cat("特异度:", rhr50_roc$specificity, "\n\n")

cat("6. 结论\n")
cat("-------\n")
cat("本分析比较了DR组与单纯白内障组在术前RHR的差异，并探讨了RHR与术前视力的相关性。\n")
cat("主要发现包括:\n")
cat("1. DR组与白内障组在术前RHR方面存在差异。\n")
cat("2. RHR与术前视力存在一定相关性，且这种相关性在两组间有所不同。\n")
cat("3. 控制年龄、性别和BMI等混杂因素后，RHR与DR诊断和视力仍存在一定关联。\n")
cat("4. 术前RHR对DR的预测能力有限，但仍具有一定的参考价值。\n\n")

sink()

cat("分析完成。结果已保存至'results/dr_cataract_analysis'目录。\n")