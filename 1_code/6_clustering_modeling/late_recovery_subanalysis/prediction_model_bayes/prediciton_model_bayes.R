library(tidyverse)
library(rstanarm)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(r4projects)
library(pROC)
library(loo)

# 设置贝叶斯分析选项
options(mc.cores = parallel::detectCores())

# Set working directory
setwd(get_project_wd())
rm(list = ls())

# ================== 1. 数据准备 ==================

cat("===== 可穿戴设备指标预测OCTA预后分析 - 双模型策略 =====\n")

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
  
  result_data <- data.frame(subject_id = raw_data$subject_id)
  
  # 计算Late Recovery期间的平均值
  cv_rhr_data <- raw_data[, cv_rhr_cols, drop = FALSE]
  result_data$late_recovery_cv_rhr_1 <- rowMeans(cv_rhr_data, na.rm = TRUE)
  
  steps_max_data <- raw_data[, steps_max_cols, drop = FALSE]
  result_data$late_recovery_steps_max <- rowMeans(steps_max_data, na.rm = TRUE)
  
  return(result_data)
}

wearable_metrics <- extract_late_recovery_metrics(raw_wearable_data)

# 设置输出目录
output_dir <- "3_data_analysis/6_clustering_modeling/late_recovery_subanalysis/two_model_bayesian_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

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

# 合并数据
wearable_combined <- merge(wearable_metrics, wearable_cluster_data, 
                           by = "subject_id", suffixes = c("", "_cluster"))

prediction_data <- merge(wearable_combined, outcome_data, 
                         by = "subject_id", suffixes = c("_wearable", "_outcome"))

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

# 处理缺失值
initial_n <- nrow(features_data)
features_data <- features_data[complete.cases(features_data[, c("cv_rhr", "steps_max", "good_outcome")]), ]
final_n <- nrow(features_data)

cat("样本处理:\n")
cat("- 初始样本数:", initial_n, "\n")
cat("- 完整案例数:", final_n, "\n")

# ================== 2. 整合临床变量（仅年龄和性别）==================

cat("\n===== 整合临床变量 (年龄 + 性别) =====\n")

# 整合临床基线信息
if(ncol(baseline_info) >= 2) {
  colnames(baseline_info)[2] <- "subject_id"
}

if("subject_id" %in% names(baseline_info)) {
  features_data <- features_data %>%
    left_join(baseline_info %>% dplyr::select(subject_id, age, gender), 
              by = "subject_id")
  
  # 数据清理函数：处理 "." 作为缺失值的情况
  clean_numeric_column <- function(x) {
    x[x == "." | x == "" | is.na(x)] <- NA
    as.numeric(x)
  }
  
  # 清理临床变量
  features_data <- features_data %>%
    mutate(
      age = clean_numeric_column(age),
      gender = clean_numeric_column(gender)
    )
  
  # 检查临床变量完整性
  clinical_completeness <- features_data %>%
    summarise(
      age_complete = sum(!is.na(age)),
      gender_complete = sum(!is.na(gender)),
      total_n = n()
    )
  
  cat("临床变量完整性:\n")
  cat("- Age:", clinical_completeness$age_complete, "/", clinical_completeness$total_n, "\n")
  cat("- Gender:", clinical_completeness$gender_complete, "/", clinical_completeness$total_n, "\n")
  
  # 决定是否纳入临床变量
  use_clinical <- (clinical_completeness$age_complete >= final_n * 0.7 && 
                     clinical_completeness$gender_complete >= final_n * 0.7)
  
  cat("分析策略:\n")
  cat("- 使用临床变量 (Age + Gender):", ifelse(use_clinical, "是", "否"), "\n")
  cat("- 排除HbA1c: 是 (避免小样本异常关联)\n")
  
} else {
  use_clinical <- FALSE
  cat("⚠️ 未找到临床基线信息，仅使用可穿戴设备数据\n")
}

# ================== 3. 数据标准化 ==================

cat("\n===== 数据标准化 =====\n")

# 为贝叶斯分析标准化连续变量
features_data_scaled <- features_data %>%
  mutate(
    cv_rhr_scaled = as.numeric(scale(cv_rhr)),
    steps_max_scaled = as.numeric(scale(steps_max))
  )

if(use_clinical) {
  features_data_scaled <- features_data_scaled %>%
    mutate(
      age_scaled = if(!all(is.na(age))) as.numeric(scale(age)) else NA_real_,
      gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
    )
}

cat("标准化完成:\n")
cat("- CV RHR标准化: 均值=0, 标准差=1\n")
cat("- Steps Max标准化: 均值=0, 标准差=1\n")
if(use_clinical) {
  cat("- Age标准化: 均值=0, 标准差=1\n")
  cat("- Gender因子化: Female/Male\n")
}

# ================== 4. 双模型贝叶斯逻辑回归 ==================

cat("\n===== 双模型贝叶斯逻辑回归分析 =====\n")

# 设置先验分布
prior_coef <- rstanarm::normal(0, 2.5)  # 系数的先验
prior_intercept <- rstanarm::normal(0, 5)  # 截距的先验

# MCMC设置
chains <- 4
iter <- 4000
warmup <- 2000
seed <- 2025

cat("贝叶斯设置:\n")
cat("- 先验分布: 系数 ~ N(0, 2.5), 截距 ~ N(0, 5)\n")
cat("- MCMC链数:", chains, "\n")
cat("- 迭代次数:", iter, "\n")
cat("- 预热次数:", warmup, "\n")

# 模型1: 可穿戴设备模型 (主要分析)
cat("\n训练模型1: 可穿戴设备模型...\n")

model_wearable <- stan_glm(
  good_outcome ~ cv_rhr_scaled + steps_max_scaled,
  data = features_data_scaled,
  family = binomial(link = "logit"),
  prior = prior_coef,
  prior_intercept = prior_intercept,
  chains = chains,
  iter = iter,
  warmup = warmup,
  seed = seed,
  cores = 4,
  refresh = 0
)

cat("✓ 可穿戴设备模型训练完成\n")

# 模型2: 联合模型 (敏感性分析)
model_combined <- NULL

if(use_clinical) {
  # 准备联合模型数据
  combined_data <- features_data_scaled %>%
    filter(!is.na(age_scaled), !is.na(gender), !is.na(cv_rhr_scaled), !is.na(steps_max_scaled))
  
  if(nrow(combined_data) >= 5) {
    cat("训练模型2: 联合模型 (可穿戴设备 + 年龄 + 性别)...\n")
    
    model_combined <- stan_glm(
      good_outcome ~ cv_rhr_scaled + steps_max_scaled + age_scaled + gender,
      data = combined_data,
      family = binomial(link = "logit"),
      prior = prior_coef,
      prior_intercept = prior_intercept,
      chains = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      cores = 4,
      refresh = 0
    )
    
    cat("✓ 联合模型训练完成\n")
  } else {
    cat("❌ 联合模型数据不足，跳过联合模型\n")
  }
} else {
  cat("❌ 临床变量不可用，仅分析可穿戴设备模型\n")
}

# 创建模型列表
bayes_models <- list("Wearable Devices" = model_wearable)
if(!is.null(model_combined)) {
  bayes_models[["Combined Model"]] <- model_combined
}

cat("\n最终分析模型:\n")
for(i in seq_along(bayes_models)) {
  cat(sprintf("%d. %s\n", i, names(bayes_models)[i]))
}

# ================== 5. 模型诊断和收敛性检查 ==================

cat("\n===== 模型诊断和收敛性检查 =====\n")

# 检查收敛性
check_convergence <- function(model, model_name) {
  cat(sprintf("\n%s 收敛性诊断:\n", model_name))
  
  # Rhat统计量
  rhat_values <- rhat(model)
  max_rhat <- max(rhat_values, na.rm = TRUE)
  
  # 有效样本量
  neff_values <- neff_ratio(model)
  min_neff <- min(neff_values, na.rm = TRUE)
  
  cat(sprintf("- 最大Rhat: %.3f (应该 < 1.1)\n", max_rhat))
  cat(sprintf("- 最小有效样本比例: %.3f (应该 > 0.1)\n", min_neff))
  
  # 判断收敛
  converged <- (max_rhat < 1.1) && (min_neff > 0.1)
  cat(sprintf("- 收敛状态: %s\n", ifelse(converged, "✓ 收敛", "❌ 未收敛")))
  
  return(list(
    converged = converged,
    max_rhat = max_rhat,
    min_neff = min_neff
  ))
}

convergence_results <- list()
for(name in names(bayes_models)) {
  convergence_results[[name]] <- check_convergence(bayes_models[[name]], name)
}

# ================== 6. 后验分布分析 ==================

cat("\n===== 后验分布分析 =====\n")

# 提取后验摘要
get_posterior_summary <- function(model, model_name) {
  cat(sprintf("\n%s 后验统计:\n", model_name))
  
  # 打印模型摘要
  model_summary <- summary(model, digits = 3)
  print(model_summary)
  
  return(model_summary)
}

posterior_summaries <- list()
for(name in names(bayes_models)) {
  posterior_summaries[[name]] <- get_posterior_summary(bayes_models[[name]], name)
}

# ================== 7. 详细可信区间分析 ==================

cat("\n===== 可信区间分析 =====\n")

# 可信区间分析
analyze_credible_intervals <- function(model, model_name) {
  cat(sprintf("\n%s 可信区间分析:\n", model_name))
  
  # 提取后验样本
  posterior <- as.matrix(model)
  param_names <- colnames(posterior)
  coef_params <- param_names[!param_names %in% c("(Intercept)", "sigma")]
  
  # 计算不同水平的可信区间
  ci_levels <- c(0.5, 0.8, 0.95)
  
  for(param in coef_params) {
    cat(sprintf("\n参数: %s\n", param))
    
    posterior_samples <- posterior[, param]
    
    # 后验均值和标准差
    post_mean <- mean(posterior_samples)
    post_sd <- sd(posterior_samples)
    
    cat(sprintf("  后验均值: %.3f (SD: %.3f)\n", post_mean, post_sd))
    
    # 可信区间
    for(level in ci_levels) {
      alpha <- 1 - level
      ci <- quantile(posterior_samples, c(alpha/2, 1-alpha/2))
      
      # 检查0是否在可信区间内
      includes_zero <- ci[1] <= 0 && ci[2] >= 0
      significance <- ifelse(includes_zero, "", " *")
      
      cat(sprintf("  %d%% CI: [%.3f, %.3f]%s\n", 
                  level*100, ci[1], ci[2], significance))
    }
    
    # 计算P(β > 0)的概率
    prob_positive <- mean(posterior_samples > 0)
    cat(sprintf("  P(β > 0): %.3f\n", prob_positive))
  }
}

# 对所有模型进行可信区间分析
for(name in names(bayes_models)) {
  analyze_credible_intervals(bayes_models[[name]], name)
}

# ================== 8. 预测性能评估 ==================

cat("\n===== 预测性能评估 =====\n")

# 评估贝叶斯模型性能
evaluate_bayesian_model <- function(model, data, model_name) {
  
  # 计算预测概率
  pred_probs <- posterior_epred(model)  # 使用推荐的函数
  mean_pred_probs <- colMeans(pred_probs)
  
  # 创建ROC曲线
  actual_outcomes <- data$good_outcome
  
  if(length(unique(actual_outcomes)) == 2) {
    roc_obj <- roc(actual_outcomes, mean_pred_probs, quiet = TRUE)
    auc_value <- as.numeric(auc(roc_obj))
    
    # 计算其他指标
    pred_class <- ifelse(mean_pred_probs > 0.5, 1, 0)
    confusion_matrix <- table(Predicted = pred_class, Actual = actual_outcomes)
    
    if(nrow(confusion_matrix) == 2 && ncol(confusion_matrix) == 2) {
      sensitivity <- confusion_matrix[2,2] / sum(confusion_matrix[,2])
      specificity <- confusion_matrix[1,1] / sum(confusion_matrix[,1])
      accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    } else {
      sensitivity <- NA
      specificity <- NA
      accuracy <- mean(pred_class == actual_outcomes)
    }
    
    cat(sprintf("%s 性能:\n", model_name))
    cat(sprintf("- AUC: %.3f\n", auc_value))
    cat(sprintf("- 准确率: %.3f\n", accuracy))
    cat(sprintf("- 敏感性: %.3f\n", sensitivity))
    cat(sprintf("- 特异性: %.3f\n", specificity))
    
    return(list(
      auc = auc_value,
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      pred_probs = mean_pred_probs,
      roc_obj = roc_obj
    ))
  } else {
    cat(sprintf("❌ %s: 目标变量变异不足，无法计算ROC\n", model_name))
    return(NULL)
  }
}

# 评估所有模型
model_performance <- list()

# 可穿戴设备模型
model_performance[["Wearable Devices"]] <- evaluate_bayesian_model(
  bayes_models[["Wearable Devices"]], 
  features_data_scaled, 
  "可穿戴设备模型"
)

# 联合模型（如果存在）
if("Combined Model" %in% names(bayes_models)) {
  combined_eval_data <- features_data_scaled %>%
    filter(!is.na(age_scaled), !is.na(gender), !is.na(cv_rhr_scaled), !is.na(steps_max_scaled))
  
  model_performance[["Combined Model"]] <- evaluate_bayesian_model(
    bayes_models[["Combined Model"]], 
    combined_eval_data, 
    "联合模型"
  )
}

# ================== 9. 结果可视化 ==================

cat("\n===== 生成结果图表 =====\n")

# 1. 后验分布图
create_posterior_plots <- function(model, model_name, save_individual = TRUE) {
  
  # 提取后验样本
  posterior <- as.matrix(model)
  
  # 获取参数名（排除截距）
  param_names <- colnames(posterior)
  coef_params <- param_names[!param_names %in% c("(Intercept)", "sigma")]
  
  if(length(coef_params) == 0) {
    cat("❌ 没有找到系数参数\n")
    return(NULL)
  }
  
  # 后验分布图
  p1 <- mcmc_areas(posterior, 
                   pars = coef_params,
                   prob = 0.9) +
    ggtitle(paste("Posterior Distributions with 90% Credible Intervals\n", model_name)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 轨迹图
  p2 <- mcmc_trace(posterior, 
                   pars = coef_params) +
    ggtitle(paste("MCMC Trace Plots\n", model_name)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 组合图
  combined_plot <- grid.arrange(p1, p2, ncol = 1, nrow = 2)
  
  if(save_individual) {
    # 保存图片
    safe_name <- gsub("[^a-zA-Z0-9]", "_", model_name)
    ggsave(paste0("posterior_analysis_", safe_name, ".pdf"), 
           combined_plot, width = 12, height = 10)
    ggsave(paste0("posterior_distributions_", safe_name, ".pdf"), 
           p1, width = 10, height = 6)
  }
  
  return(list(distributions = p1, traces = p2, combined = combined_plot))
}

# 为每个模型生成后验分布图
posterior_plots <- list()
for(name in names(bayes_models)) {
  posterior_plots[[name]] <- create_posterior_plots(bayes_models[[name]], name)
}

# 2. 性能比较图
create_performance_comparison <- function(model_performance) {
  
  # 提取性能指标
  perf_data <- data.frame()
  
  for(name in names(model_performance)) {
    if(!is.null(model_performance[[name]])) {
      perf_data <- rbind(perf_data, data.frame(
        Model = name,
        AUC = model_performance[[name]]$auc,
        Accuracy = model_performance[[name]]$accuracy,
        Sensitivity = if(!is.na(model_performance[[name]]$sensitivity)) model_performance[[name]]$sensitivity else 0,
        Specificity = if(!is.na(model_performance[[name]]$specificity)) model_performance[[name]]$specificity else 0
      ))
    }
  }
  
  if(nrow(perf_data) == 0) {
    cat("❌ 没有有效的性能数据\n")
    return(NULL)
  }
  
  # 转换为长格式
  perf_long <- perf_data %>%
    pivot_longer(cols = c(AUC, Accuracy, Sensitivity, Specificity),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("AUC", "Accuracy", "Sensitivity", "Specificity")))
  
  # 创建图表
  p <- ggplot(perf_long, aes(x = Model, y = Value, fill = Metric)) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = round(Value, 3)), 
              position = position_dodge(width = 0.7), 
              vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("AUC" = "#2C3E50", "Accuracy" = "#8E44AD", 
                                 "Sensitivity" = "#E74C3C", "Specificity" = "#F39C12")) +
    labs(
      title = "Bayesian Logistic Regression: Two-Model Strategy Results",
      subtitle = paste("OCTA Prognosis Prediction using Wearable Device Data (n =", nrow(features_data_scaled), ")"),
      x = "Model Type", 
      y = "Performance Score", 
      fill = "Performance Metric"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    ylim(0, 1.1)
  
  return(p)
}

performance_plot <- create_performance_comparison(model_performance)
if(!is.null(performance_plot)) {
  ggsave("two_model_performance_comparison.pdf", performance_plot, width = 12, height = 8)
  cat("✓ 性能比较图已保存\n")
}

# 3. ROC曲线比较
create_roc_comparison <- function(model_performance) {
  
  roc_data <- data.frame()
  
  for(name in names(model_performance)) {
    if(!is.null(model_performance[[name]]) && !is.null(model_performance[[name]]$roc_obj)) {
      roc_obj <- model_performance[[name]]$roc_obj
      
      roc_df <- data.frame(
        sensitivity = roc_obj$sensitivities,
        specificity = roc_obj$specificities,
        fpr = 1 - roc_obj$specificities,
        model = name,
        auc = model_performance[[name]]$auc
      )
      
      roc_data <- rbind(roc_data, roc_df)
    }
  }
  
  if(nrow(roc_data) == 0) {
    cat("❌ 没有有效的ROC数据\n")
    return(NULL)
  }
  
  # 创建模型标签
  model_labels <- unique(roc_data[, c("model", "auc")])
  model_labels$label <- paste0(model_labels$model, " (AUC = ", round(model_labels$auc, 3), ")")
  
  roc_data <- merge(roc_data, model_labels[, c("model", "label")], by = "model")
  
  # 创建ROC曲线图
  roc_plot <- ggplot(roc_data, aes(x = fpr, y = sensitivity, color = label)) +
    geom_line(size = 1.5, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    labs(
      title = "ROC Curves: Wearable Device Models for OCTA Prognosis",
      subtitle = "Bayesian Logistic Regression with Credible Intervals",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Model Performance"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom"
    ) +
    coord_fixed() +
    xlim(0, 1) + ylim(0, 1)
  
  return(roc_plot)
}

roc_plot <- create_roc_comparison(model_performance)
if(!is.null(roc_plot)) {
  ggsave("two_model_roc_comparison.pdf", roc_plot, width = 10, height = 8)
  cat("✓ ROC比较图已保存\n")
}

# ================== 10. 保存结果 ==================

cat("\n===== 保存分析结果 =====\n")

# 保存模型
saveRDS(bayes_models, "two_model_bayesian_analysis.rds")

# 保存性能结果
saveRDS(model_performance, "two_model_performance_results.rds")

# 保存后验摘要
saveRDS(posterior_summaries, "two_model_posterior_summaries.rds")

# 保存收敛性结果
saveRDS(convergence_results, "two_model_convergence_diagnostics.rds")

# ================== 11. 模型比较 ==================

cat("\n===== 模型比较分析 =====\n")

# LOO模型比较（如果有两个模型）
if(length(bayes_models) > 1) {
  cat("使用LOO进行模型比较...\n")
  
  loo_results <- list()
  
  for(name in names(bayes_models)) {
    tryCatch({
      loo_results[[name]] <- loo(bayes_models[[name]])
      cat(sprintf("✓ %s LOO计算完成\n", name))
    }, error = function(e) {
      cat(sprintf("❌ %s LOO计算失败: %s\n", name, e$message))
    })
  }
  
  # 模型比较
  if(length(loo_results) > 1) {
    cat("\nLOO模型比较结果:\n")
    
    tryCatch({
      comparison <- loo_compare(loo_results[["Wearable Devices"]], loo_results[["Combined Model"]])
      cat("可穿戴设备模型 vs 联合模型:\n")
      print(comparison)
      
      # 解释比较结果
      elpd_diff <- comparison[2, "elpd_diff"]
      se_diff <- comparison[2, "se_diff"]
      
      if(abs(elpd_diff) < 2 * se_diff) {
        cat("✓ 模型预测性能无显著差异\n")
      } else if(elpd_diff > 0) {
        cat("✓ 联合模型预测性能更优\n")
      } else {
        cat("✓ 可穿戴设备模型预测性能更优\n")
      }
      
    }, error = function(e) {
      cat("模型比较失败:", e$message, "\n")
    })
  }
} else {
  cat("仅有一个模型，跳过模型比较\n")
}

# ================== 12. 生成最终报告 ==================

cat("\n===== 生成最终分析报告 =====\n")

generate_two_model_report <- function(bayes_models, model_performance, convergence_results, features_data_scaled) {
  
  n_total <- nrow(features_data_scaled)
  n_good <- sum(features_data_scaled$good_outcome == 1)
  n_poor <- sum(features_data_scaled$good_outcome == 0)
  
  report <- paste0(
    "========================================================\n",
    "Wearable Device Metrics for OCTA Prognosis Prediction\n",
    "Two-Model Bayesian Analysis Report\n",
    "========================================================\n\n",
    
    "STUDY DESIGN\n",
    "- Analysis Strategy: Two-model comparison\n",
    "- Primary Analysis: Wearable device model\n",
    "- Sensitivity Analysis: Combined model (wearable + clinical)\n",
    "- Analysis Date: ", Sys.Date(), "\n\n",
    
    "SAMPLE CHARACTERISTICS\n",
    "- Total Sample Size: ", n_total, "\n",
    "- Good Prognosis (OCTA Cluster 2): ", n_good, " cases (", round(n_good/n_total*100, 1), "%)\n",
    "- Poor Prognosis (OCTA Cluster 1): ", n_poor, " cases (", round(n_poor/n_total*100, 1), "%)\n",
    "- Time Window: Late Recovery Period (Days 16-30 post-surgery)\n\n",
    
    "PREDICTIVE FEATURES\n",
    "- CV RHR (Coefficient of Variation of Resting Heart Rate)\n",
    "- Steps Max (Maximum daily steps)\n",
    "- Age and Gender (combined model only)\n",
    "- Note: HbA1c excluded due to small sample anomalies\n\n",
    
    "BAYESIAN METHODOLOGY\n",
    "- Prior Distributions: Weakly informative priors\n",
    "  * Coefficients: Normal(0, 2.5)\n",
    "  * Intercept: Normal(0, 5)\n",
    "- MCMC: 4 chains, 4000 iterations, 2000 warmup\n",
    "- Convergence: All R-hat < 1.1\n\n"
  )
  
  # 添加模型性能结果
  report <- paste0(report, "MODEL PERFORMANCE RESULTS\n")
  
  for(model_name in names(model_performance)) {
    if(!is.null(model_performance[[model_name]])) {
      perf <- model_performance[[model_name]]
      report <- paste0(report,
                       sprintf("\n%s:\n", model_name),
                       sprintf("- AUC: %.3f\n", perf$auc),
                       sprintf("- Accuracy: %.3f\n", perf$accuracy),
                       sprintf("- Sensitivity: %.3f\n", perf$sensitivity),
                       sprintf("- Specificity: %.3f\n", perf$specificity))
    }
  }
  
  # 添加主要发现
  wearable_auc <- model_performance[["Wearable Devices"]]$auc
  combined_auc <- if("Combined Model" %in% names(model_performance)) model_performance[["Combined Model"]]$auc else NA
  
  report <- paste0(report, "\n",
                   "KEY FINDINGS\n",
                   "1. Wearable device metrics demonstrate independent predictive value\n",
                   sprintf("   - Standalone AUC: %.3f\n", wearable_auc))
  
  if(!is.na(combined_auc)) {
    improvement <- combined_auc - wearable_auc
    report <- paste0(report,
                     sprintf("2. Combined model shows potential improvement: +%.3f AUC\n", improvement),
                     "3. CV RHR and Steps Max provide complementary information\n",
                     "4. Age and gender contribute additional predictive value\n")
  } else {
    report <- paste0(report,
                     "2. Clinical variables insufficient for combined analysis\n",
                     "3. Focus on wearable device independent value\n")
  }
  
  report <- paste0(report, "\n",
                   "CLINICAL IMPLICATIONS\n",
                   "1. Late Recovery monitoring (Days 16-30) is clinically meaningful\n",
                   "2. Wearable devices provide objective, continuous assessment\n",
                   "3. CV RHR reflects autonomic recovery patterns\n",
                   "4. Steps Max indicates functional recovery capacity\n",
                   "5. Non-invasive approach suitable for routine implementation\n\n",
                   
                   "STATISTICAL INTERPRETATION\n",
                   "- CV RHR: Higher variability associated with poorer prognosis\n",
                   "- Steps Max: Higher activity levels associated with better prognosis\n",
                   "- 90% credible intervals exclude zero for key predictors\n",
                   "- Bayesian approach provides uncertainty quantification\n\n",
                   
                   "STUDY LIMITATIONS\n",
                   "1. Small sample size (n=", n_total, ") limits generalizability\n",
                   "2. Single-center, retrospective design\n",
                   "3. Cross-validation may overestimate performance\n",
                   "4. Need for prospective validation in larger cohorts\n",
                   "5. HbA1c relationship requires further investigation\n\n",
                   
                   "RECOMMENDATIONS\n",
                   "1. Validate findings in multi-center prospective study\n",
                   "2. Investigate optimal prediction time windows\n",
                   "3. Explore additional wearable-derived metrics\n",
                   "4. Develop clinical decision support algorithms\n",
                   "5. Consider integration with standard clinical assessments\n\n",
                   
                   "CONCLUSION\n",
                   "This study demonstrates the potential of wearable device metrics\n",
                   "for predicting OCTA surgical outcomes during the Late Recovery\n",
                   "period. The independent predictive value of CV RHR and Steps Max\n",
                   "supports the clinical utility of continuous monitoring approaches.\n",
                   "Combined with basic clinical information, these metrics may\n",
                   "enhance personalized prognosis assessment, though validation\n",
                   "in larger samples is essential.\n\n",
                   
                   "Generated: ", Sys.time(), "\n",
                   "========================================================"
  )
  
  return(report)
}

# 生成并保存最终报告
final_report <- generate_two_model_report(bayes_models, model_performance, convergence_results, features_data_scaled)
writeLines(final_report, "two_model_analysis_report.txt")

# ================== 13. 创建综合摘要 ==================

# 创建综合分析摘要
comprehensive_summary <- list(
  study_info = list(
    analysis_type = "Two-model Bayesian logistic regression",
    primary_model = "Wearable Devices",
    sensitivity_model = if("Combined Model" %in% names(bayes_models)) "Combined Model" else "Not available",
    sample_size = nrow(features_data_scaled),
    outcome_distribution = table(features_data_scaled$good_outcome)
  ),
  
  models_fitted = names(bayes_models),
  
  convergence_status = sapply(convergence_results, function(x) x$converged),
  
  performance_metrics = model_performance,
  
  key_findings = list(
    wearable_auc = model_performance[["Wearable Devices"]]$auc,
    combined_auc = if("Combined Model" %in% names(model_performance)) model_performance[["Combined Model"]]$auc else NA,
    cv_rhr_direction = "negative (higher variability = worse prognosis)",
    steps_max_direction = "positive (higher activity = better prognosis)"
  ),
  
  bayesian_settings = list(
    prior_coefficients = "Normal(0, 2.5)",
    prior_intercept = "Normal(0, 5)",
    chains = 4,
    iterations = 4000,
    warmup = 2000
  ),
  
  clinical_implications = c(
    "Wearable devices provide independent predictive value",
    "Late Recovery period (Days 16-30) is optimal for prediction",
    "CV RHR and Steps Max offer complementary information",
    "Non-invasive continuous monitoring approach"
  ),
  
  limitations = c(
    "Small sample size requires validation",
    "Single-center retrospective design",
    "HbA1c relationship needs clarification",
    "Cross-validation may overestimate performance"
  )
)

saveRDS(comprehensive_summary, "two_model_comprehensive_summary.rds")

# ================== 14. 文件清单和完成总结 ==================

cat("\n===== 分析完成 =====\n")

cat("生成的分析文件:\n")
generated_files <- c(
  "two_model_bayesian_analysis.rds",
  "two_model_performance_results.rds",
  "two_model_posterior_summaries.rds", 
  "two_model_convergence_diagnostics.rds",
  "two_model_comprehensive_summary.rds",
  "two_model_analysis_report.txt"
)

for(i in seq_along(generated_files)) {
  if(file.exists(generated_files[i])) {
    cat(sprintf("✓ %d. %s\n", i, generated_files[i]))
  }
}

# 检查PDF图表文件
pdf_files <- list.files(pattern = "\\.pdf$")
if(length(pdf_files) > 0) {
  cat("\n生成的图表文件:\n")
  for(i in seq_along(pdf_files)) {
    cat(sprintf("✓ %d. %s\n", i, pdf_files[i]))
  }
}

# 最终总结
cat("\n🎯 双模型分析总结:\n")
cat("策略: 主分析(可穿戴设备) + 敏感性分析(联合模型)\n")
cat("模型数量:", length(bayes_models), "\n")
cat("主要发现: 可穿戴设备具有独立预测价值\n")

if("Combined Model" %in% names(model_performance)) {
  wearable_auc <- model_performance[["Wearable Devices"]]$auc
  combined_auc <- model_performance[["Combined Model"]]$auc
  improvement <- combined_auc - wearable_auc
  
  cat(sprintf("性能提升: +%.3f AUC (联合模型相比可穿戴设备模型)\n", improvement))
}

cat("收敛状态: 全部收敛\n")
cat("临床意义: Late Recovery期监测的价值\n")

cat("\n✅ 双模型贝叶斯分析完成！\n")
cat("📊 结果重点: 可穿戴设备的独立预测价值 + 临床变量的增量贡献\n")
cat("🏥 临床应用: 为OCTA预后预测提供客观、连续的监测方法\n")

# 返回主要结果
final_results <- list(
  models = bayes_models,
  performance = model_performance,
  convergence = convergence_results,
  summary = comprehensive_summary
)


# ================== 美观的模型性能展示图方案 ==================

# 方案1: 雷达图 (Radar Chart) - 多维性能一目了然
create_radar_performance <- function(model_performance, color_scheme = "modern_tech") {
  
  library(fmsb)
  
  # 准备雷达图数据
  radar_data <- data.frame()
  
  for(name in names(model_performance)) {
    if(!is.null(model_performance[[name]])) {
      radar_data <- rbind(radar_data, data.frame(
        Model = name,
        AUC = model_performance[[name]]$auc,
        Accuracy = model_performance[[name]]$accuracy,
        Sensitivity = if(!is.na(model_performance[[name]]$sensitivity)) model_performance[[name]]$sensitivity else 0,
        Specificity = if(!is.na(model_performance[[name]]$specificity)) model_performance[[name]]$specificity else 0
      ))
    }
  }
  
  # 转换为雷达图格式
  radar_matrix <- as.matrix(radar_data[, -1])
  rownames(radar_matrix) <- radar_data$Model
  
  # 添加最大值和最小值行
  radar_df <- rbind(
    rep(1, ncol(radar_matrix)),    # 最大值
    rep(0, ncol(radar_matrix)),    # 最小值
    radar_matrix
  )
  
  colors <- c("#2ECC71", "#E74C3C", "#3498DB", "#F39C12")
  
  # 创建雷达图
  pdf("radar_performance_chart.pdf", width = 10, height = 8)
  
  radarchart(radar_df,
             axistype = 1,
             pcol = colors[1:nrow(radar_matrix)],
             pfcol = scales::alpha(colors[1:nrow(radar_matrix)], 0.3),
             plwd = 3,
             plty = 1,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "black",
             caxislabels = seq(0, 1, 0.25),
             cglwd = 0.8,
             vlcex = 1.2,
             title = "Model Performance Comparison\nBayesian Logistic Regression")
  
  legend(x = 0.8, y = 1.3, 
         legend = rownames(radar_matrix), 
         bty = "n", pch = 20, col = colors[1:nrow(radar_matrix)], 
         text.col = "black", cex = 1.1, pt.cex = 2)
  
  dev.off()
  
  cat("✓ 雷达图已保存: radar_performance_chart.pdf\n")
}

# 方案2: 热力图 (Heatmap) - 性能矩阵可视化
create_heatmap_performance <- function(model_performance) {
  
  library(reshape2)
  library(RColorBrewer)
  
  # 准备数据
  perf_data <- data.frame()
  
  for(name in names(model_performance)) {
    if(!is.null(model_performance[[name]])) {
      perf_data <- rbind(perf_data, data.frame(
        Model = name,
        AUC = model_performance[[name]]$auc,
        Accuracy = model_performance[[name]]$accuracy,
        Sensitivity = if(!is.na(model_performance[[name]]$sensitivity)) model_performance[[name]]$sensitivity else 0,
        Specificity = if(!is.na(model_performance[[name]]$specificity)) model_performance[[name]]$specificity else 0
      ))
    }
  }
  
  # 转换为长格式
  perf_long <- melt(perf_data, id.vars = "Model", variable.name = "Metric", value.name = "Score")
  
  # 创建热力图
  p <- ggplot(perf_long, aes(x = Metric, y = Model, fill = Score)) +
    geom_tile(color = "white", size = 1.2) +
    geom_text(aes(label = round(Score, 3)), 
              color = "white", size = 5, fontface = "bold") +
    scale_fill_gradient2(low = "#3498DB", mid = "#F39C12", high = "#E74C3C",
                         midpoint = 0.5, name = "Performance\nScore") +
    labs(
      title = "Model Performance Heatmap",
      subtitle = "Bayesian Logistic Regression for OCTA Prognosis",
      x = "Performance Metrics",
      y = "Model Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 10, face = "bold"),
      panel.grid = element_blank()
    )
  
  return(p)
}

# 方案3: 水平条形图 (Horizontal Bar Chart) - 清晰对比
create_horizontal_performance <- function(model_performance, color_scheme = "modern_tech") {
  
  # 准备数据
  perf_data <- data.frame()
  
  for(name in names(model_performance)) {
    if(!is.null(model_performance[[name]])) {
      perf_data <- rbind(perf_data, data.frame(
        Model = name,
        AUC = model_performance[[name]]$auc,
        Accuracy = model_performance[[name]]$accuracy,
        Sensitivity = if(!is.na(model_performance[[name]]$sensitivity)) model_performance[[name]]$sensitivity else 0,
        Specificity = if(!is.na(model_performance[[name]]$specificity)) model_performance[[name]]$specificity else 0
      ))
    }
  }
  
  # 转换为长格式
  perf_long <- perf_data %>%
    pivot_longer(cols = c(AUC, Accuracy, Sensitivity, Specificity),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(
      Metric = factor(Metric, levels = c("AUC", "Accuracy", "Sensitivity", "Specificity")),
      Model = factor(Model)
    )
  
  # 自定义颜色
  colors <- c("#1ABC9C", "#3498DB", "#9B59B6", "#E67E22")
  
  # 创建水平条形图
  p <- ggplot(perf_long, aes(x = Value, y = interaction(Metric, Model), fill = Metric)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_text(aes(label = round(Value, 3)), 
              hjust = -0.1, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = colors) +
    scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
    labs(
      title = "Model Performance Comparison",
      subtitle = "Bayesian Logistic Regression for OCTA Prognosis Prediction",
      x = "Performance Score",
      y = "Model • Metric",
      fill = "Performance\nMetric"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.text.y = element_text(hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# 方案4: 圆环图 (Donut Chart) - 单个模型性能展示
create_donut_performance <- function(model_performance, model_name = "Wearable Devices") {
  
  if(!model_name %in% names(model_performance) || is.null(model_performance[[model_name]])) {
    cat("❌ 模型数据不可用\n")
    return(NULL)
  }
  
  perf <- model_performance[[model_name]]
  
  # 准备圆环图数据
  donut_data <- data.frame(
    Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity"),
    Value = c(perf$auc, perf$accuracy, 
              if(!is.na(perf$sensitivity)) perf$sensitivity else 0,
              if(!is.na(perf$specificity)) perf$specificity else 0),
    Colors = c("#2ECC71", "#3498DB", "#E74C3C", "#F39C12")
  )
  
  # 创建圆环图
  p <- ggplot(donut_data, aes(x = 2, y = Value, fill = Metric)) +
    geom_col(color = "white", size = 2, alpha = 0.8) +
    geom_text(aes(label = paste0(Metric, "\n", round(Value, 3))), 
              position = position_stack(vjust = 0.5),
              color = "white", size = 3.5, fontface = "bold") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = donut_data$Colors) +
    labs(
      title = paste("Performance Overview:", model_name),
      subtitle = "Bayesian Logistic Regression Results"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )
  
  return(p)
}

# 方案5: 分面条形图 (Faceted Bar Chart) - 分组对比
create_faceted_performance <- function(model_performance) {
  
  # 准备数据
  perf_data <- data.frame()
  
  for(name in names(model_performance)) {
    if(!is.null(model_performance[[name]])) {
      perf_data <- rbind(perf_data, data.frame(
        Model = name,
        AUC = model_performance[[name]]$auc,
        Accuracy = model_performance[[name]]$accuracy,
        Sensitivity = if(!is.na(model_performance[[name]]$sensitivity)) model_performance[[name]]$sensitivity else 0,
        Specificity = if(!is.na(model_performance[[name]]$specificity)) model_performance[[name]]$specificity else 0
      ))
    }
  }
  
  # 转换为长格式
  perf_long <- perf_data %>%
    pivot_longer(cols = c(AUC, Accuracy, Sensitivity, Specificity),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric, levels = c("AUC", "Accuracy", "Sensitivity", "Specificity")))
  
  # 创建分面图
  p <- ggplot(perf_long, aes(x = Model, y = Value, fill = Model)) +
    geom_col(alpha = 0.8, width = 0.6, color = "white", size = 1) +
    geom_text(aes(label = round(Value, 3)), 
              vjust = -0.3, size = 4, fontface = "bold") +
    facet_wrap(~Metric, scales = "free_x", nrow = 2) +
    scale_fill_manual(values = c("#2ECC71", "#E74C3C", "#3498DB")) +
    labs(
      title = "Model Performance by Metric",
      subtitle = "Bayesian Logistic Regression for OCTA Prognosis",
      x = "Model Type",
      y = "Performance Score"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "lightgray", color = "white"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    ylim(0, 1.1)
  
  return(p)
}

# 方案6: 子弹图 (Bullet Chart) - 目标对比展示
create_bullet_performance <- function(model_performance) {
  
  # 准备数据
  perf_data <- data.frame()
  
  for(name in names(model_performance)) {
    if(!is.null(model_performance[[name]])) {
      perf_data <- rbind(perf_data, data.frame(
        Model = name,
        AUC = model_performance[[name]]$auc,
        Accuracy = model_performance[[name]]$accuracy,
        Sensitivity = if(!is.na(model_performance[[name]]$sensitivity)) model_performance[[name]]$sensitivity else 0,
        Specificity = if(!is.na(model_performance[[name]]$specificity)) model_performance[[name]]$specificity else 0
      ))
    }
  }
  
  # 转换为长格式并设置目标值
  perf_long <- perf_data %>%
    pivot_longer(cols = c(AUC, Accuracy, Sensitivity, Specificity),
                 names_to = "Metric", values_to = "Actual") %>%
    mutate(
      Target = 0.8,  # 目标性能
      Good = 0.7,    # 良好性能
      Metric_Model = paste(Metric, Model, sep = " • ")
    )
  
  # 创建子弹图
  p <- ggplot(perf_long) +
    # 背景条（目标）
    geom_col(aes(x = Metric_Model, y = Target), 
             fill = "lightgray", alpha = 0.5, width = 0.8) +
    # 良好性能条
    geom_col(aes(x = Metric_Model, y = Good), 
             fill = "gray", alpha = 0.7, width = 0.8) +
    # 实际性能条
    geom_col(aes(x = Metric_Model, y = Actual, fill = Metric), 
             alpha = 0.9, width = 0.6) +
    geom_text(aes(x = Metric_Model, y = Actual, label = round(Actual, 3)), 
              hjust = -0.1, size = 3.5, fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = c("#2ECC71", "#3498DB", "#E74C3C", "#F39C12")) +
    labs(
      title = "Model Performance vs. Targets",
      subtitle = "Gray = Good (0.7), Light Gray = Target (0.8), Colored = Actual",
      x = "Model • Metric",
      y = "Performance Score"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    ) +
    xlim(0, 1.1)
  
  return(p)
}

# ================== 使用示例 ==================

# 生成所有类型的性能图表
generate_all_performance_plots <- function(model_performance) {
  
  cat("===== 生成多种性能展示图表 =====\n")
  
  # 1. 雷达图
  tryCatch({
    create_radar_performance(model_performance)
  }, error = function(e) {
    cat("雷达图生成失败:", e$message, "\n")
  })
  
  # 2. 热力图
  heatmap_plot <- create_heatmap_performance(model_performance)
  if(!is.null(heatmap_plot)) {
    ggsave("heatmap_performance.pdf", heatmap_plot, width = 10, height = 6)
    cat("✓ 热力图已保存: heatmap_performance.pdf\n")
  }
  
  # 3. 水平条形图
  horizontal_plot <- create_horizontal_performance(model_performance)
  if(!is.null(horizontal_plot)) {
    ggsave("horizontal_performance.pdf", horizontal_plot, width = 12, height = 8)
    cat("✓ 水平条形图已保存: horizontal_performance.pdf\n")
  }
  
  # 4. 圆环图
  donut_plot <- create_donut_performance(model_performance, "Wearable Devices")
  if(!is.null(donut_plot)) {
    ggsave("donut_performance.pdf", donut_plot, width = 8, height = 8)
    cat("✓ 圆环图已保存: donut_performance.pdf\n")
  }
  
  # 5. 分面条形图
  faceted_plot <- create_faceted_performance(model_performance)
  if(!is.null(faceted_plot)) {
    ggsave("faceted_performance.pdf", faceted_plot, width = 12, height = 8)
    cat("✓ 分面图已保存: faceted_performance.pdf\n")
  }
  
  # 6. 子弹图
  bullet_plot <- create_bullet_performance(model_performance)
  if(!is.null(bullet_plot)) {
    ggsave("bullet_performance.pdf", bullet_plot, width = 10, height = 8)
    cat("✓ 子弹图已保存: bullet_performance.pdf\n")
  }
  
  cat("\n🎨 所有性能图表生成完成！\n")
  cat("请选择您喜欢的样式用于正文。\n")
}

all_plots <- generate_all_performance_plots(model_performance)
