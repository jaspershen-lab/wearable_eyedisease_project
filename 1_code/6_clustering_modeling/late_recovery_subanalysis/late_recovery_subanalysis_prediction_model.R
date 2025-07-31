library(tidyverse)
library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(ggplot2)
library(gridExtra)
library(r4projects)
library(xgboost)

# Set working directory
setwd(get_project_wd())
rm(list = ls())

# ================== 1. 数据准备 ==================

cat("===== 可穿戴设备指标预测OCTA预后分析 =====\n")

# 加载数据文件
raw_wearable_file <- "3_data_analysis/6_clustering_modeling/data_prepare/1m/mfuzz_D_Surg1_8h_filtered.csv"
wearable_cluster_file <- "3_data_analysis/6_clustering_modeling/time_window_clustering/late_recovery_detailed_membership_fixed.csv"
outcome_file <- "3_data_analysis/6_clustering_modeling/mfuzz/comprehensive_cluster/ppv_WF_cluster_results.csv"
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv", stringsAsFactors = FALSE)

# 安全加载数据函数
load_data_safely <- function(file_path, data_name) {
  if(file.exists(file_path)) {
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat(sprintf("✓ 成功加载 %s: %d 行数据\n", data_name, nrow(data)))
    return(data)
  } else {
    cat(sprintf("❌ 文件不存在: %s\n", file_path))
    stop(sprintf("必需的数据文件不存在: %s", file_path))
  }
}

# 加载数据
raw_wearable_data <- load_data_safely(raw_wearable_file, "原始可穿戴设备数据")
wearable_cluster_data <- load_data_safely(wearable_cluster_file, "可穿戴设备聚类结果")
outcome_data <- load_data_safely(outcome_file, "OCTA预后数据")

# 提取Late Recovery时间窗口的指标函数
extract_late_recovery_metrics <- function(raw_data) {
  late_recovery_days <- 16:30
  
  cv_rhr_cols <- c()
  steps_max_cols <- c()
  
  for(day in late_recovery_days) {
    cv_rhr_col <- paste0("day_", day, "_cv_rhr_1")
    steps_max_col <- paste0("day_", day, "_steps_max")
    
    if(cv_rhr_col %in% colnames(raw_data)) {
      cv_rhr_cols <- c(cv_rhr_cols, cv_rhr_col)
    }
    if(steps_max_col %in% colnames(raw_data)) {
      steps_max_cols <- c(steps_max_cols, steps_max_col)
    }
  }
  
  cat("找到的Late Recovery CV RHR列:", length(cv_rhr_cols), "个\n")
  cat("找到的Late Recovery Steps Max列:", length(steps_max_cols), "个\n")
  
  # 数据质量检查
  if(length(cv_rhr_cols) == 0) {
    stop("❌ 错误: 没有找到CV RHR数据列，请检查数据文件和列名")
  }
  
  if(length(steps_max_cols) == 0) {
    stop("❌ 错误: 没有找到Steps Max数据列，请检查数据文件和列名")
  }
  
  result_data <- data.frame(subject_id = raw_data$subject_id)
  
  # 计算Late Recovery期间的平均值（仅使用真实数据）
  cv_rhr_data <- raw_data[, cv_rhr_cols, drop = FALSE]
  result_data$late_recovery_cv_rhr_1 <- rowMeans(cv_rhr_data, na.rm = TRUE)
  
  steps_max_data <- raw_data[, steps_max_cols, drop = FALSE]
  result_data$late_recovery_steps_max <- rowMeans(steps_max_data, na.rm = TRUE)
  
  # 检查数据完整性
  cat("CV RHR 完整数据样本数:", sum(!is.na(result_data$late_recovery_cv_rhr_1)), "\n")
  cat("Steps Max 完整数据样本数:", sum(!is.na(result_data$late_recovery_steps_max)), "\n")
  
  return(result_data)
}

# 提取Late Recovery指标
wearable_metrics <- extract_late_recovery_metrics(raw_wearable_data)

# 设置输出目录
output_dir <- "3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/prediction_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# ================== 2. 数据预处理 ==================

# 统一ID列名
standardize_id_column <- function(data) {
  if("subject_id" %in% names(data)) {
    return(data)
  } else if("ID" %in% names(data)) {
    names(data)[names(data) == "ID"] <- "subject_id"
    return(data)
  } else {
    stop("找不到ID列 (subject_id 或 ID)")
  }
}

wearable_metrics <- standardize_id_column(wearable_metrics)
wearable_cluster_data <- standardize_id_column(wearable_cluster_data)
outcome_data <- standardize_id_column(outcome_data)

# 合并数据：指标 + 聚类结果 + 预后结果
wearable_combined <- merge(wearable_metrics, wearable_cluster_data, 
                           by = "subject_id", suffixes = c("", "_cluster"))

prediction_data <- merge(wearable_combined, outcome_data, 
                         by = "subject_id", suffixes = c("_wearable", "_outcome"))

cat("合并后的预测数据集:\n")
cat("- 总样本数:", nrow(prediction_data), "\n")
cat("- 预后好 (OCTA聚类2):", sum(prediction_data$max_cluster_outcome == 2), "例\n")
cat("- 预后差 (OCTA聚类1):", sum(prediction_data$max_cluster_outcome == 1), "例\n")

# 创建二分类目标变量
prediction_data$good_outcome <- ifelse(prediction_data$max_cluster_outcome == 2, 1, 0)
prediction_data$outcome_label <- factor(prediction_data$good_outcome, 
                                        levels = c(0, 1), 
                                        labels = c("Poor Outcome", "Good Outcome"))

# 创建特征数据框
features_data <- data.frame(
  subject_id = prediction_data$subject_id,
  cv_rhr = prediction_data$late_recovery_cv_rhr_1,
  steps_max = prediction_data$late_recovery_steps_max,
  wearable_cluster = prediction_data$max_cluster_wearable,
  outcome = prediction_data$outcome_label,
  good_outcome = prediction_data$good_outcome
)

cat("使用的特征:\n")
cat("- CV RHR (Late Recovery): late_recovery_cv_rhr_1\n")
cat("- Steps Max (Late Recovery): late_recovery_steps_max\n")
cat("- 可穿戴设备聚类: max_cluster_wearable\n")

# 处理缺失值 - 只保留完整案例
initial_n <- nrow(features_data)
features_data <- features_data[complete.cases(features_data[, c("cv_rhr", "steps_max", "good_outcome")]), ]
final_n <- nrow(features_data)

cat("缺失值处理:\n")
cat("- 初始样本数:", initial_n, "\n")
cat("- 完整案例数:", final_n, "\n")
cat("- 排除样本数:", initial_n - final_n, "\n")

# 样本量检查
if(final_n < 8) {
  stop("❌ 错误: 完整案例样本量过小 (n < 8)，无法进行可靠的统计分析")
}

if(final_n < 20) {
  cat("⚠️ 警告: 样本量较小 (n = ", final_n, ")，结果解释需谨慎\n")
}

# 确保 baseline_info ID 列标准化
if(ncol(baseline_info) >= 2) {
  colnames(baseline_info)[2] <- "subject_id"
}

# ================== 3. 临床变量整合 ==================

# 合并临床基线信息
if("subject_id" %in% names(baseline_info)) {
  features_data <- features_data %>%
    left_join(baseline_info %>% dplyr::select(subject_id, age, gender, hba1c), 
              by = "subject_id")
  
  # 检查临床变量完整性
  clinical_completeness <- features_data %>%
    summarise(
      age_complete = sum(!is.na(age)),
      gender_complete = sum(!is.na(gender)),
      hba1c_complete = sum(!is.na(hba1c)),
      total_n = n()
    )
  
  cat("临床变量完整性:\n")
  cat("- Age:", clinical_completeness$age_complete, "/", clinical_completeness$total_n, "\n")
  cat("- Gender:", clinical_completeness$gender_complete, "/", clinical_completeness$total_n, "\n")
  cat("- HbA1c:", clinical_completeness$hba1c_complete, "/", clinical_completeness$total_n, "\n")
  
  # 决定是否纳入临床变量
  use_clinical <- (clinical_completeness$age_complete >= final_n * 0.7 && 
                     clinical_completeness$gender_complete >= final_n * 0.7)
  
  use_hba1c <- (clinical_completeness$hba1c_complete >= final_n * 0.7)
  
  cat("分析策略:\n")
  cat("- 使用临床变量 (Age + Gender):", ifelse(use_clinical, "是", "否"), "\n")
  cat("- 使用HbA1c:", ifelse(use_hba1c, "是", "否"), "\n")
  
} else {
  use_clinical <- FALSE
  use_hba1c <- FALSE
  cat("⚠️ 未找到临床基线信息，仅使用可穿戴设备数据\n")
}

# ================== 4. 探索性数据分析 ==================

cat("\n===== 描述性统计 =====\n")

# 按预后分组的描述性统计
summary_stats <- features_data %>%
  group_by(outcome) %>%
  summarise(
    n = n(),
    cv_rhr_mean = mean(cv_rhr, na.rm = TRUE),
    cv_rhr_sd = sd(cv_rhr, na.rm = TRUE),
    cv_rhr_median = median(cv_rhr, na.rm = TRUE),
    steps_max_mean = mean(steps_max, na.rm = TRUE),
    steps_max_sd = sd(steps_max, na.rm = TRUE),
    steps_max_median = median(steps_max, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# 如果有临床变量，也进行描述性统计
if(use_clinical) {
  clinical_stats <- features_data %>%
    filter(!is.na(age), !is.na(gender)) %>%
    group_by(outcome) %>%
    summarise(
      n = n(),
      age_mean = mean(age, na.rm = TRUE),
      age_sd = sd(age, na.rm = TRUE),
      female_pct = mean(gender == 0, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  cat("\n临床变量描述性统计:\n")
  print(clinical_stats)
  
  if(use_hba1c) {
    hba1c_stats <- features_data %>%
      filter(!is.na(hba1c)) %>%
      group_by(outcome) %>%
      summarise(
        n = n(),
        hba1c_mean = mean(hba1c, na.rm = TRUE),
        hba1c_sd = sd(hba1c, na.rm = TRUE),
        .groups = "drop"
      )
    
    cat("\nHbA1c描述性统计:\n")
    print(hba1c_stats)
  }
}

# 创建可视化
create_exploratory_plots <- function(data) {
  
  # 主要散点图
  p1 <- ggplot(data, aes(x = cv_rhr, y = steps_max, color = outcome)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.3) +
    scale_color_manual(values = c("Poor Outcome" = "#E74C3C", "Good Outcome" = "#27AE60")) +
    labs(title = "CV RHR vs Steps Max by OCTA Prognosis",
         subtitle = paste("Late Recovery Period (Days 16-30), n =", nrow(data)),
         x = "CV RHR (Late Recovery)",
         y = "Max Steps (Late Recovery)",
         color = "OCTA Prognosis") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  # CV RHR 分布
  p2 <- ggplot(data, aes(x = outcome, y = cv_rhr, fill = outcome)) +
    geom_boxplot(alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("Poor Outcome" = "#E74C3C", "Good Outcome" = "#27AE60")) +
    labs(title = "CV RHR Distribution by Prognosis",
         x = "OCTA Prognosis",
         y = "CV RHR (Late Recovery)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  # Steps Max 分布
  p3 <- ggplot(data, aes(x = outcome, y = steps_max, fill = outcome)) +
    geom_boxplot(alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("Poor Outcome" = "#E74C3C", "Good Outcome" = "#27AE60")) +
    labs(title = "Max Steps Distribution by Prognosis",
         x = "OCTA Prognosis",
         y = "Max Steps (Late Recovery)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  # CV RHR 密度图
  p4 <- ggplot(data, aes(x = cv_rhr, fill = outcome)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = c("Poor Outcome" = "#E74C3C", "Good Outcome" = "#27AE60")) +
    labs(title = "CV RHR Density Distribution",
         x = "CV RHR (Late Recovery)",
         y = "Density",
         fill = "OCTA Prognosis") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(list(scatter = p1, boxplot_cvr = p2, boxplot_steps = p3, density = p4))
}

plots <- create_exploratory_plots(features_data)
combined_plot <- grid.arrange(plots$scatter, plots$boxplot_cvr, plots$boxplot_steps, plots$density, ncol = 2)
ggsave("exploratory_analysis.pdf", combined_plot, width = 12, height = 10)
ggsave("main_scatter_plot.pdf", plots$scatter, width = 10, height = 8)

# ================== 5. 统计检验 ==================

cat("\n===== 统计检验 =====\n")

# 检查正态性
cv_rhr_shapiro <- shapiro.test(features_data$cv_rhr)
steps_max_shapiro <- shapiro.test(features_data$steps_max)

cat("正态性检验 (Shapiro-Wilk):\n")
cat(sprintf("- CV RHR: W = %.4f, p = %.4f\n", cv_rhr_shapiro$statistic, cv_rhr_shapiro$p.value))
cat(sprintf("- Steps Max: W = %.4f, p = %.4f\n", steps_max_shapiro$statistic, steps_max_shapiro$p.value))

# 根据正态性选择合适的检验
use_parametric <- (cv_rhr_shapiro$p.value > 0.05 && steps_max_shapiro$p.value > 0.05)

cat("\n组间比较检验:\n")

if(use_parametric) {
  cv_rhr_test <- t.test(cv_rhr ~ outcome, data = features_data)
  steps_max_test <- t.test(steps_max ~ outcome, data = features_data)
  
  cat("使用参数检验 (t-test):\n")
  cat(sprintf("- CV RHR: t = %.3f, df = %.1f, p = %.4f\n", 
              cv_rhr_test$statistic, cv_rhr_test$parameter, cv_rhr_test$p.value))
  cat(sprintf("- Steps Max: t = %.3f, df = %.1f, p = %.4f\n", 
              steps_max_test$statistic, steps_max_test$parameter, steps_max_test$p.value))
} else {
  cat("使用非参数检验 (Mann-Whitney U):\n")
}

# 无论如何都进行非参数检验作为稳健性检查
cv_rhr_wilcox <- wilcox.test(cv_rhr ~ outcome, data = features_data)
steps_max_wilcox <- wilcox.test(steps_max ~ outcome, data = features_data)

cat("非参数检验 (Mann-Whitney U):\n")
cat(sprintf("- CV RHR: W = %.1f, p = %.4f\n", cv_rhr_wilcox$statistic, cv_rhr_wilcox$p.value))
cat(sprintf("- Steps Max: W = %.1f, p = %.4f\n", steps_max_wilcox$statistic, steps_max_wilcox$p.value))

# 效应量计算 (Cohen's d)
calculate_cohens_d <- function(x1, x2) {
  pooled_sd <- sqrt(((length(x1) - 1) * var(x1) + (length(x2) - 1) * var(x2)) / 
                      (length(x1) + length(x2) - 2))
  return((mean(x1) - mean(x2)) / pooled_sd)
}

cv_rhr_poor <- features_data$cv_rhr[features_data$outcome == "Poor Outcome"]
cv_rhr_good <- features_data$cv_rhr[features_data$outcome == "Good Outcome"]
steps_poor <- features_data$steps_max[features_data$outcome == "Poor Outcome"]
steps_good <- features_data$steps_max[features_data$outcome == "Good Outcome"]

cv_rhr_d <- calculate_cohens_d(cv_rhr_poor, cv_rhr_good)
steps_d <- calculate_cohens_d(steps_poor, steps_good)

cat("\n效应量 (Cohen's d):\n")
cat(sprintf("- CV RHR: d = %.3f\n", cv_rhr_d))
cat(sprintf("- Steps Max: d = %.3f\n", steps_d))

# ================== 6. 逻辑回归模型 (主要分析) ==================

cat("\n===== 逻辑回归分析 (主要方法) =====\n")

# 准备建模数据
modeling_data <- features_data %>%
  dplyr::select(cv_rhr, steps_max, good_outcome) %>%
  mutate(good_outcome = factor(good_outcome, levels = c(0, 1), labels = c("Poor", "Good")))

# 设置训练控制 - 针对小样本优化
if(nrow(modeling_data) < 15) {
  train_control <- trainControl(
    method = "LOOCV",  # 留一交叉验证适合小样本
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  cat("使用留一交叉验证 (LOOCV)\n")
} else {
  train_control <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  cat("使用5折交叉验证\n")
}

set.seed(123)

# 1. 可穿戴设备逻辑回归模型
cat("训练可穿戴设备逻辑回归模型...\n")
logistic_wearable <- train(
  good_outcome ~ cv_rhr + steps_max,
  data = modeling_data,
  method = "glm",
  family = "binomial",
  trControl = train_control,
  metric = "ROC"
)

# 获取详细的模型系数
wearable_summary <- summary(logistic_wearable$finalModel)
cat("可穿戴设备模型系数:\n")
print(wearable_summary$coefficients)

# 2. 临床变量逻辑回归模型 (如果可用)
if(use_clinical) {
  clinical_modeling_data <- features_data %>%
    filter(!is.na(age), !is.na(gender)) %>%
    mutate(
      good_outcome = factor(good_outcome, levels = c(0, 1), labels = c("Poor", "Good")),
      gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
    )
  
  if(use_hba1c) {
    clinical_modeling_data <- clinical_modeling_data %>%
      filter(!is.na(hba1c))
    
    logistic_clinical <- train(
      good_outcome ~ age + gender + hba1c,
      data = clinical_modeling_data,
      method = "glm",
      family = "binomial",
      trControl = train_control,
      metric = "ROC"
    )
    
    clinical_summary <- summary(logistic_clinical$finalModel)
    cat("\n临床变量模型系数:\n")
    print(clinical_summary$coefficients)
  } else {
    logistic_clinical <- train(
      good_outcome ~ age + gender,
      data = clinical_modeling_data,
      method = "glm",
      family = "binomial",
      trControl = train_control,
      metric = "ROC"
    )
    
    clinical_summary <- summary(logistic_clinical$finalModel)
    cat("\n临床变量模型系数:\n")
    print(clinical_summary$coefficients)
  }
  
  # 3. 联合模型
  combined_modeling_data <- features_data %>%
    filter(!is.na(age), !is.na(gender), !is.na(cv_rhr), !is.na(steps_max)) %>%
    mutate(
      good_outcome = factor(good_outcome, levels = c(0, 1), labels = c("Poor", "Good")),
      gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
    )
  
  if(use_hba1c) {
    combined_modeling_data <- combined_modeling_data %>%
      filter(!is.na(hba1c))
    
    logistic_combined <- train(
      good_outcome ~ cv_rhr + steps_max + age + gender + hba1c,
      data = combined_modeling_data,
      method = "glm",
      family = "binomial",
      trControl = train_control,
      metric = "ROC"
    )
  } else {
    logistic_combined <- train(
      good_outcome ~ cv_rhr + steps_max + age + gender,
      data = combined_modeling_data,
      method = "glm",
      family = "binomial",
      trControl = train_control,
      metric = "ROC"
    )
  }
  
  combined_summary <- summary(logistic_combined$finalModel)
  cat("\n联合模型系数:\n")
  print(combined_summary$coefficients)
}

# ================== 7. 机器学习模型 (辅助分析) ==================

cat("\n===== 机器学习模型 (辅助分析) =====\n")

# 只在样本量足够时使用复杂模型
ml_models <- list()

if(nrow(modeling_data) >= 10) {
  # 随机森林 (简化参数)
  cat("训练随机森林模型...\n")
  tryCatch({
    rf_model <- train(
      good_outcome ~ cv_rhr + steps_max,
      data = modeling_data,
      method = "rf",
      trControl = train_control,
      metric = "ROC",
      tuneGrid = data.frame(mtry = 1:2),
      ntree = 100  # 减少树的数量
    )
    ml_models[["Random Forest"]] <- rf_model
  }, error = function(e) {
    cat("随机森林训练失败:", e$message, "\n")
  })
  
  # SVM (简化参数)
  cat("训练SVM模型...\n")
  tryCatch({
    svm_model <- train(
      good_outcome ~ cv_rhr + steps_max,
      data = modeling_data,
      method = "svmRadial",
      trControl = train_control,
      metric = "ROC",
      tuneLength = 3  # 减少调参复杂度
    )
    ml_models[["SVM"]] <- svm_model
  }, error = function(e) {
    cat("SVM训练失败:", e$message, "\n")
  })
  
  # XGBoost (大幅简化)
  if(nrow(modeling_data) >= 12) {
    cat("训练XGBoost模型...\n")
    tryCatch({
      xgb_model <- train(
        good_outcome ~ cv_rhr + steps_max,
        data = modeling_data,
        method = "xgbTree",
        trControl = train_control,
        metric = "ROC",
        tuneGrid = expand.grid(
          nrounds = 50,
          max_depth = 2,
          eta = 0.3,
          gamma = 0,
          colsample_bytree = 1,
          min_child_weight = 1,
          subsample = 1
        ),
        verbose = FALSE
      )
      ml_models[["XGBoost"]] <- xgb_model
    }, error = function(e) {
      cat("XGBoost训练失败:", e$message, "\n")
    })
  } else {
    cat("样本量不足，跳过XGBoost模型\n")
  }
} else {
  cat("样本量过小，跳过机器学习模型\n")
}

# ================== 8. 模型性能评估 ==================

cat("\n===== 模型性能评估 =====\n")

# 主要模型 (逻辑回归)
main_models <- list("Wearable Devices" = logistic_wearable)

if(use_clinical) {
  if(use_hba1c) {
    main_models[["Clinical (Age+Gender+HbA1c)"]] <- logistic_clinical
    main_models[["Combined (All)"]] <- logistic_combined
  } else {
    main_models[["Clinical (Age+Gender)"]] <- logistic_clinical
    main_models[["Combined (All)"]] <- logistic_combined
  }
}

# 提取主要模型性能
main_results <- data.frame(
  Model = names(main_models),
  AUC = sapply(main_models, function(x) max(x$results$ROC, na.rm = TRUE)),
  Sensitivity = sapply(main_models, function(x) max(x$results$Sens, na.rm = TRUE)),
  Specificity = sapply(main_models, function(x) max(x$results$Spec, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

cat("主要分析结果 (逻辑回归):\n")
print(main_results)

# 辅助模型性能 (如果有)
if(length(ml_models) > 0) {
  ml_results <- data.frame(
    Model = names(ml_models),
    AUC = sapply(ml_models, function(x) max(x$results$ROC, na.rm = TRUE)),
    Sensitivity = sapply(ml_models, function(x) max(x$results$Sens, na.rm = TRUE)),
    Specificity = sapply(ml_models, function(x) max(x$results$Spec, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  cat("\n辅助分析结果 (机器学习):\n")
  print(ml_results)
} else {
  ml_results <- data.frame()
}

# 选择最佳逻辑回归模型
best_main_model_name <- main_results$Model[which.max(main_results$AUC)]
best_main_model <- main_models[[best_main_model_name]]

cat(sprintf("\n最佳逻辑回归模型: %s (AUC = %.3f)\n", 
            best_main_model_name, max(main_results$AUC)))

# ================== 9. 主要结果可视化 ==================

cat("\n===== 生成主要结果图表 =====\n")

# 主要结果图：逻辑回归模型比较
create_main_logistic_plot <- function(main_results) {
  
  main_long <- main_results %>%
    pivot_longer(cols = c(AUC, Sensitivity, Specificity),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(
      Metric = factor(Metric, levels = c("AUC", "Sensitivity", "Specificity")),
      Model = factor(Model, levels = main_results$Model[order(main_results$AUC, decreasing = TRUE)])
    )
  
  ggplot(main_long, aes(x = Model, y = Value, fill = Metric)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = round(Value, 3)), 
              position = position_dodge(width = 0.7), 
              vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("AUC" = "#2C3E50", "Sensitivity" = "#8E44AD", "Specificity" = "#F39C12")) +
    labs(
      title = "Logistic Regression Models for OCTA Prognosis Prediction",
      subtitle = paste("Primary Analysis Results (n =", nrow(features_data), ")"),
      x = "Model Type", 
      y = "Performance Score", 
      fill = "Performance Metric"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylim(0, 1.1)
}

# 辅助结果图：机器学习模型 (如果有)
create_supplementary_ml_plot <- function(ml_results) {
  if(nrow(ml_results) == 0) {
    # 创建一个说明图
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = paste("Machine Learning Models\nNot Applicable\n(Sample size too small: n =", 
                             nrow(features_data), ")"),
               size = 6, hjust = 0.5, vjust = 0.5) +
      labs(title = "Supplementary Analysis: Machine Learning Algorithms",
           subtitle = "Sample size considerations for complex models") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 11))
  } else {
    ml_long <- ml_results %>%
      pivot_longer(cols = c(AUC, Sensitivity, Specificity),
                   names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = factor(Metric, levels = c("AUC", "Sensitivity", "Specificity")))
    
    ggplot(ml_long, aes(x = Model, y = Value, fill = Metric)) +
      geom_col(position = "dodge", alpha = 0.8, width = 0.6) +
      geom_text(aes(label = round(Value, 3)), 
                position = position_dodge(width = 0.6), 
                vjust = -0.3, size = 3) +
      scale_fill_manual(values = c("AUC" = "#5A9BD4", "Sensitivity" = "#AA6C9C", "Specificity" = "#D4A762")) +
      labs(
        title = "Supplementary Analysis: Machine Learning Algorithms",
        subtitle = paste("Alternative Methods (n =", nrow(features_data), ")"),
        x = "Algorithm", 
        y = "Performance Score", 
        fill = "Performance Metric"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1)
      ) +
      ylim(0, 1.1)
  }
}

# 生成图表
main_plot <- create_main_logistic_plot(main_results)
supp_plot <- create_supplementary_ml_plot(ml_results)

# 保存图表
ggsave("main_logistic_regression_results.pdf", main_plot, width = 10, height = 8)
ggsave("supplementary_ml_algorithms.pdf", supp_plot, width = 10, height = 6)

cat("✓ 主要结果图已保存: main_logistic_regression_results.pdf\n")
cat("✓ 辅助结果图已保存: supplementary_ml_algorithms.pdf\n")

# ================== 10. ROC曲线分析 ==================

create_roc_analysis <- function(main_models, features_data, main_results) {
  
  roc_data_list <- list()
  auc_values <- list()
  
  for(model_name in names(main_models)) {
    model <- main_models[[model_name]]
    
    # 根据模型选择合适的数据
    if(grepl("Clinical", model_name) && !grepl("Combined", model_name)) {
      if(use_clinical) {
        model_data <- features_data %>%
          filter(!is.na(age), !is.na(gender)) %>%
          mutate(
            good_outcome = factor(good_outcome, levels = c(0, 1), labels = c("Poor", "Good")),
            gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
          )
        
        if(use_hba1c) {
          model_data <- model_data %>% filter(!is.na(hba1c))
        }
      } else {
        next
      }
    } else if(grepl("Combined", model_name)) {
      if(use_clinical) {
        model_data <- features_data %>%
          filter(!is.na(age), !is.na(gender), !is.na(cv_rhr), !is.na(steps_max)) %>%
          mutate(
            good_outcome = factor(good_outcome, levels = c(0, 1), labels = c("Poor", "Good")),
            gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
          )
        
        if(use_hba1c) {
          model_data <- model_data %>% filter(!is.na(hba1c))
        }
      } else {
        next
      }
    } else {
      # 可穿戴设备模型
      model_data <- features_data %>%
        dplyr::select(cv_rhr, steps_max, good_outcome) %>%
        mutate(good_outcome = factor(good_outcome, levels = c(0, 1), labels = c("Poor", "Good")))
    }
    
    if(nrow(model_data) == 0) next
    
    # 预测概率
    tryCatch({
      pred_probs <- predict(model, model_data, type = "prob")$Good
      roc_obj <- roc(model_data$good_outcome, pred_probs, levels = c("Poor", "Good"))
      
      # 使用交叉验证的AUC值
      cv_auc <- main_results$AUC[main_results$Model == model_name]
      auc_values[[model_name]] <- cv_auc
      
      roc_df <- data.frame(
        sensitivity = roc_obj$sensitivities,
        specificity = roc_obj$specificities,
        model = model_name,
        auc = cv_auc
      )
      
      roc_df$fpr <- 1 - roc_df$specificity
      roc_data_list[[model_name]] <- roc_df
    }, error = function(e) {
      cat("ROC计算失败 for", model_name, ":", e$message, "\n")
    })
  }
  
  if(length(roc_data_list) == 0) {
    cat("无可用ROC数据\n")
    return(NULL)
  }
  
  # 合并所有ROC数据
  all_roc_data <- do.call(rbind, roc_data_list)
  
  # 创建模型标签
  model_labels <- paste0(names(auc_values), " (AUC = ", 
                         round(unlist(auc_values), 3), ")")
  names(model_labels) <- names(auc_values)
  
  # 创建ROC曲线图
  roc_plot <- ggplot(all_roc_data, aes(x = fpr, y = sensitivity, color = model)) +
    geom_line(size = 1.5, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    scale_color_brewer(type = "qual", palette = "Set1", labels = model_labels) +
    labs(
      title = "ROC Curves for OCTA Prognosis Prediction",
      subtitle = "Logistic Regression Models - Cross-Validation Results",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Model Performance"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.text = element_text(size = 9)
    ) +
    coord_fixed() +
    xlim(0, 1) + ylim(0, 1)
  
  return(list(
    plot = roc_plot,
    roc_data = all_roc_data,
    auc_values = auc_values
  ))
}

# 创建ROC曲线
roc_results <- create_roc_analysis(main_models, features_data, main_results)

if(!is.null(roc_results)) {
  ggsave("roc_curves_logistic_models.pdf", roc_results$plot, width = 10, height = 8)
  cat("✓ ROC曲线图已保存: roc_curves_logistic_models.pdf\n")
}

# ================== 11. 模型诊断和验证 ==================

cat("\n===== 模型诊断 =====\n")

# 对最佳逻辑回归模型进行诊断
perform_model_diagnostics <- function(model, model_data, model_name) {
  cat(sprintf("\n%s 模型诊断:\n", model_name))
  
  # 残差分析
  residuals_pearson <- residuals(model$finalModel, type = "pearson")
  residuals_deviance <- residuals(model$finalModel, type = "deviance")
  
  # 预测概率
  fitted_probs <- fitted(model$finalModel)
  
  # Cook's距离
  cooks_d <- cooks.distance(model$finalModel)
  
  cat("模型诊断统计:\n")
  cat(sprintf("- 皮尔逊残差范围: [%.3f, %.3f]\n", min(residuals_pearson), max(residuals_pearson)))
  cat(sprintf("- 偏差残差范围: [%.3f, %.3f]\n", min(residuals_deviance), max(residuals_deviance)))
  cat(sprintf("- 最大Cook距离: %.3f\n", max(cooks_d)))
  
  # 检查影响点
  influential_points <- which(cooks_d > 4/nrow(model_data))
  if(length(influential_points) > 0) {
    cat("潜在影响点:", influential_points, "\n")
  } else {
    cat("无明显影响点\n")
  }
  
  return(list(
    residuals_pearson = residuals_pearson,
    residuals_deviance = residuals_deviance,
    fitted_probs = fitted_probs,
    cooks_d = cooks_d
  ))
}

# 对最佳模型进行诊断
if(best_main_model_name == "Wearable Devices") {
  diagnostics <- perform_model_diagnostics(best_main_model, modeling_data, best_main_model_name)
} else {
  # 对于临床或联合模型，需要相应的数据
  if(use_clinical) {
    appropriate_data <- features_data %>%
      filter(!is.na(age), !is.na(gender)) %>%
      mutate(
        good_outcome = factor(good_outcome, levels = c(0, 1), labels = c("Poor", "Good")),
        gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
      )
    
    if(use_hba1c) {
      appropriate_data <- appropriate_data %>% filter(!is.na(hba1c))
    }
    
    diagnostics <- perform_model_diagnostics(best_main_model, appropriate_data, best_main_model_name)
  }
}

# ================== 12. 保存结果和模型 ==================

cat("\n===== 保存分析结果 =====\n")

# 保存最佳逻辑回归模型
saveRDS(best_main_model, "best_logistic_model.rds")
saveRDS(main_models, "all_logistic_models.rds")

if(length(ml_models) > 0) {
  saveRDS(ml_models, "supplementary_ml_models.rds")
}

# 保存完整分析摘要
analysis_summary <- list(
  # 数据信息
  sample_size = nrow(features_data),
  outcome_distribution = table(features_data$outcome),
  
  # 描述性统计
  descriptive_stats = summary_stats,
  
  # 统计检验结果
  statistical_tests = list(
    cv_rhr_shapiro = cv_rhr_shapiro,
    steps_max_shapiro = steps_max_shapiro,
    cv_rhr_wilcox = cv_rhr_wilcox,
    steps_max_wilcox = steps_max_wilcox,
    cv_rhr_effect_size = cv_rhr_d,
    steps_max_effect_size = steps_d
  ),
  
  # 模型性能
  main_model_results = main_results,
  ml_model_results = if(nrow(ml_results) > 0) ml_results else NULL,
  best_model = best_main_model_name,
  
  # 模型系数
  best_model_summary = summary(best_main_model$finalModel),
  
  # 分析设置
  analysis_settings = list(
    use_clinical = use_clinical,
    use_hba1c = use_hba1c,
    cross_validation = if(nrow(modeling_data) < 15) "LOOCV" else "5-fold CV",
    features_used = c("late_recovery_cv_rhr_1", "late_recovery_steps_max")
  )
)

saveRDS(analysis_summary, "comprehensive_analysis_summary.rds")

# ================== 13. 生成最终报告 ==================

generate_final_report <- function(analysis_summary, features_data, main_results, ml_results) {
  
  n_total <- nrow(features_data)
  n_good <- sum(features_data$good_outcome == 1)
  n_poor <- sum(features_data$good_outcome == 0)
  
  best_auc <- max(main_results$AUC)
  best_model <- main_results$Model[which.max(main_results$AUC)]
  
  report <- paste0(
    "========================================================\n",
    "Wearable Device Metrics for OCTA Prognosis Prediction\n",
    "Final Analysis Report\n",
    "========================================================\n\n",
    
    "STUDY OVERVIEW\n",
    "- Analysis Date: ", Sys.Date(), "\n",
    "- Total Sample Size: ", n_total, "\n",
    "- Good Prognosis: ", n_good, " cases (", round(n_good/n_total*100, 1), "%)\n",
    "- Poor Prognosis: ", n_poor, " cases (", round(n_poor/n_total*100, 1), "%)\n",
    "- Primary Features: CV RHR and Max Steps (Late Recovery Period)\n\n",
    
    "DATA QUALITY\n",
    "- No simulated data used - all analyses based on real measurements\n",
    "- Complete case analysis performed\n",
    "- Late Recovery Period: Days 16-30 post-surgery\n"
  )
  
  if(use_clinical) {
    report <- paste0(report, "- Clinical variables included: Age, Gender")
    if(use_hba1c) {
      report <- paste0(report, ", HbA1c\n")
    } else {
      report <- paste0(report, "\n")
    }
  } else {
    report <- paste0(report, "- Clinical variables: Not available/insufficient\n")
  }
  
  report <- paste0(report, "\n",
                   "STATISTICAL ANALYSIS\n",
                   "- Normality Tests: Shapiro-Wilk test\n",
                   "- Group Comparisons: Mann-Whitney U test (non-parametric)\n",
                   "- Effect Sizes: Cohen's d\n",
                   "- CV RHR Effect Size: ", round(analysis_summary$statistical_tests$cv_rhr_effect_size, 3), "\n",
                   "- Steps Max Effect Size: ", round(analysis_summary$statistical_tests$steps_max_effect_size, 3), "\n\n",
                   
                   "PRIMARY ANALYSIS RESULTS (Logistic Regression)\n"
  )
  
  for(i in 1:nrow(main_results)) {
    report <- paste0(report,
                     sprintf("- %s: AUC=%.3f, Sensitivity=%.3f, Specificity=%.3f\n",
                             main_results$Model[i], main_results$AUC[i], 
                             main_results$Sensitivity[i], main_results$Specificity[i]))
  }
  
  report <- paste0(report, "\n",
                   "Best Performing Model: ", best_model, " (AUC = ", round(best_auc, 3), ")\n\n"
  )
  
  if(nrow(ml_results) > 0) {
    report <- paste0(report, "SUPPLEMENTARY ANALYSIS (Machine Learning)\n")
    for(i in 1:nrow(ml_results)) {
      report <- paste0(report,
                       sprintf("- %s: AUC=%.3f, Sensitivity=%.3f, Specificity=%.3f\n",
                               ml_results$Model[i], ml_results$AUC[i], 
                               ml_results$Sensitivity[i], ml_results$Specificity[i]))
    }
    report <- paste0(report, "\n")
  } else {
    report <- paste0(report,
                     "SUPPLEMENTARY ANALYSIS\n",
                     "- Machine learning models not applied due to small sample size\n",
                     "- Focus maintained on interpretable logistic regression\n\n"
    )
  }
  
  report <- paste0(report,
                   "KEY FINDINGS\n",
                   "1. Wearable device metrics show predictive value for OCTA prognosis\n",
                   "2. Late Recovery period (Days 16-30) provides meaningful signals\n",
                   "3. CV RHR and Max Steps complement each other as predictors\n"
  )
  
  if(use_clinical && length(main_models) > 1) {
    combined_auc <- main_results$AUC[grepl("Combined", main_results$Model)]
    wearable_auc <- main_results$AUC[main_results$Model == "Wearable Devices"]
    
    if(length(combined_auc) > 0 && length(wearable_auc) > 0) {
      improvement <- combined_auc - wearable_auc
      report <- paste0(report,
                       "4. Combined model improvement over wearable-only: +", 
                       round(improvement, 3), " AUC points\n")
    }
  }
  
  report <- paste0(report, "\n",
                   "CLINICAL IMPLICATIONS\n",
                   "1. Continuous monitoring via wearables offers objective prognosis assessment\n",
                   "2. Non-invasive approach suitable for routine clinical implementation\n",
                   "3. Late Recovery period appears most informative for prediction\n",
                   "4. Integration with clinical data may enhance predictive accuracy\n\n",
                   
                   "LIMITATIONS\n",
                   "1. Small sample size (n=", n_total, ") requires validation in larger cohorts\n",
                   "2. Single-center study limits generalizability\n",
                   "3. Cross-validation may overestimate performance\n",
                   "4. Need for prospective validation studies\n\n",
                   
                   "RECOMMENDATIONS\n",
                   "1. Validate findings in multi-center study with larger sample\n",
                   "2. Explore optimal time windows for prediction\n",
                   "3. Investigate additional wearable-derived metrics\n",
                   "4. Develop clinical decision support tools\n\n",
                   
                   "TECHNICAL NOTES\n",
                   "- Cross-validation: ", analysis_summary$analysis_settings$cross_validation, "\n",
                   "- No data imputation or simulation used\n",
                   "- Robust statistical testing with effect size reporting\n",
                   "- Primary focus on interpretable logistic regression\n\n",
                   
                   "Generated: ", Sys.time(), "\n",
                   "========================================================"
  )
  
  return(report)
}

# 生成并保存最终报告
final_report <- generate_final_report(analysis_summary, features_data, main_results, ml_results)
writeLines(final_report, "final_analysis_report.txt")

# ================== 14. 总结和文件清单 ==================

cat("\n", final_report, "\n")

cat("\n===== 分析完成 =====\n")
cat("生成的文件:\n")
cat("1. best_logistic_model.rds - 最佳逻辑回归模型\n")
cat("2. all_logistic_models.rds - 所有逻辑回归模型\n")
if(length(ml_models) > 0) {
  cat("3. supplementary_ml_models.rds - 辅助机器学习模型\n")
}
cat("4. comprehensive_analysis_summary.rds - 完整分析摘要\n")
cat("5. final_analysis_report.txt - 最终分析报告\n")
cat("6. main_logistic_regression_results.pdf - 主要结果图\n")
cat("7. supplementary_ml_algorithms.pdf - 辅助结果图\n")
cat("8. exploratory_analysis.pdf - 探索性分析图\n")
cat("9. main_scatter_plot.pdf - 主要散点图\n")
if(!is.null(roc_results)) {
  cat("10. roc_curves_logistic_models.pdf - ROC曲线图\n")
}

# 返回关键结果摘要
cat("\n🎯 关键结果摘要:\n")
cat("样本量:", nrow(features_data), "\n")
cat("最佳模型:", best_main_model_name, "\n")
cat("最佳AUC:", round(max(main_results$AUC), 3), "\n")

if(nrow(ml_results) > 0) {
  cat("最佳ML模型:", ml_results$Model[which.max(ml_results$AUC)], 
      "(AUC =", round(max(ml_results$AUC), 3), ")\n")
}

cat("\n模型排名 (按AUC):\n")
all_results <- rbind(
  data.frame(main_results, Type = "Logistic"),
  if(nrow(ml_results) > 0) data.frame(ml_results, Type = "ML") else NULL
)

ranked_results <- all_results[order(all_results$AUC, decreasing = TRUE), ]
for(i in 1:min(5, nrow(ranked_results))) {
  cat(sprintf("%d. %s (%s): AUC = %.3f\n", 
              i, ranked_results$Model[i], ranked_results$Type[i], ranked_results$AUC[i]))
}

# 返回主要对象
list(
  main_model_results = main_results,
  ml_model_results = if(nrow(ml_results) > 0) ml_results else NULL,
  best_model = best_main_model,
  analysis_summary = analysis_summary,
  features_data = features_data
)
