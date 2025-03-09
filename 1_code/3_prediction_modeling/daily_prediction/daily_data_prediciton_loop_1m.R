# Comprehensive Analysis for Vision Improvement Prediction
# This script analyzes which day's data best predicts post-operative vision improvement at 1 month

library(glmnet)
library(caret)
library(randomForest)
library(xgboost)
library(tidyverse)
library(pROC)
library(ggplot2)
library(mice)
library(VIM)
setwd(get_project_wd())
rm(list = ls())

# Function to process data from a specific day
analyze_day_data <- function(day_data, day_label) {
  # Remove rows where vision_improvement_1m is NA (our target variable)
  day_data <- day_data[!is.na(day_data$vision_improvement_1m), ]
  
  # Select features for imputation
  features_for_imputation <- day_data %>%
    dplyr::select(
      # RHR features
      mean_rhr_1, min_rhr_1, max_rhr_1, median_rhr_1, sd_rhr_1, iqr_rhr_1, skew_rhr_1, kurt_rhr_1,
      # BO features
      mean_bo, min_bo, max_bo, median_bo, sd_bo, iqr_bo, skew_bo, kurt_bo,
      # Steps features
      steps_total, steps_mean, steps_max,
      # Sleep features
      total_sleep,deep_sleep,
      # One-week post-surgery vision
      post_vision_1w, vision_improvement_1w,
      # Demographics and medical history
      age, gender, cataract_2, dm_2, hypertension_2, pre_vision,season,
      # Target variable
      vision_improvement_1m
    )
  
  # First, check for missing data
  missing_summary <- sapply(features_for_imputation, function(x) sum(is.na(x)))
  
  # Only perform imputation if there are missing values
  if(sum(missing_summary) > 0) {
    # Set the seed for reproducibility
    set.seed(1234)
    
    # Configure the imputation method
    imputation_methods <- make.method(features_for_imputation)
    
    # Create 5 imputed datasets
    imp <- mice(features_for_imputation, m=5, method=imputation_methods, 
                maxit=50, seed=1234, printFlag=FALSE)
    
    # Create a complete dataset using the first imputation
    imputed_data <- complete(imp, 1)
    
    # Replace the missing values in the original dataset
    day_data_imputed <- day_data
    for(col in names(features_for_imputation)) {
      if(sum(is.na(day_data[[col]])) > 0) {
        day_data_imputed[[col]] <- imputed_data[[col]]
      }
    }
  } else {
    day_data_imputed <- day_data
  }
  
  # Define selected features
  selected_features <- c(
    # RHR features
    "mean_rhr_1", "min_rhr_1", "max_rhr_1", "median_rhr_1", "sd_rhr_1", "iqr_rhr_1", 
    # BO features
    "mean_bo", "min_bo", "max_bo", "median_bo", "sd_bo", "iqr_bo", 
    # Steps features
    "steps_total", "steps_mean", "steps_max",
    #sleep
    " total_sleep","deep_sleep",
    # One-week post-surgery vision data
    "post_vision_1w", "vision_improvement_1w",
    # Demographics and medical history
    "age", "gender", "cataract_2", "dm_2", "hypertension_2", "pre_vision","season"
  )
  
  # Prepare data for LASSO
  # 确保所有选定的特征都存在于数据集中
  available_features <- intersect(selected_features, colnames(day_data_imputed))
  
  # 检查是否有特征缺失，如果有则输出警告
  if(length(available_features) < length(selected_features)) {
    missing_features <- setdiff(selected_features, colnames(day_data_imputed))
    warning(paste("以下特征在", day_label, "的数据中不存在:", paste(missing_features, collapse=", ")))
  }
  
  # 为分类变量创建模型矩阵
  categorical_vars <- c("gender", "dm_2", "cataract_2", "hypertension_2","season")
  categorical_in_model <- intersect(categorical_vars, available_features)
  
  # Add debugging code here, right before the model matrix creation
  cat("Day:", day_label, "- Data columns:", paste(colnames(day_data_imputed), collapse=", "), "\n")
  cat("Available features:", paste(available_features, collapse=", "), "\n")
  cat("Rows after filtering:", nrow(day_data_imputed), "\n")
  cat("Categorical variables in model:", paste(categorical_in_model, collapse=", "), "\n")
  
  # Modify the model matrix creation code
  if(length(categorical_in_model) > 0) {
    # Check if there are any categorical variables left after filtering
    # Build formula string with remaining categorical variables
    formula_str <- paste("~ 0 +", paste(available_features, collapse = " + "))
    x_formula <- as.formula(formula_str)
    
    # Print debugging information
    cat("Creating model matrix with formula:", formula_str, "\n")
    
    # Use tryCatch to handle potential errors in model.matrix
    x <- tryCatch({
      model.matrix(x_formula, data = day_data_imputed)
    }, error = function(e) {
      cat("Error in model.matrix:", e$message, "\n")
      cat("Trying alternative approach without model.matrix...\n")
      # If model.matrix fails, use direct matrix conversion
      as.matrix(day_data_imputed[, available_features])
    })
  } else {
    # No categorical variables, use simple matrix conversion
    cat("No categorical variables in model, using simple matrix conversion\n")
    x <- as.matrix(day_data_imputed[, available_features])
  }
  
  y <- as.numeric(day_data_imputed$vision_improvement_1m)
  
  # Fit LASSO model with cross-validation
  set.seed(1234)
  cv_fit <- cv.glmnet(x, y, type.measure="deviance", alpha=1, nfolds=5)
  
  # Extract LASSO-selected features
  feature_all <- as.data.frame(as.matrix(coef(cv_fit, s=cv_fit$lambda.min)))
  colnames(feature_all) <- "coff"
  feature_opt <- feature_all %>% filter(abs(coff) > 0)
  lasso_selected_features <- rownames(feature_opt)
  
  # Remove intercept from selected features if present
  lasso_selected_features <- lasso_selected_features[lasso_selected_features != "(Intercept)"]
  
  # 准备用于训练的特征集
  # 需要将model.matrix产生的哑变量名转换回原始特征名
  categorical_vars <- c("gender", "dm_2", "cataract_2", "hypertension_2","season")
  
  # 添加检查：如果LASSO没有选出任何特征，则使用所有可用特征
  if(length(lasso_selected_features) == 0) {
    cat("LASSO没有选出任何特征，使用所有可用特征\n")
    original_selected_features <- available_features
  } else {
    # 从LASSO选择的特征中提取原始变量名
    original_selected_features <- c()
    
    for(feature in lasso_selected_features) {
      # 处理哑变量名称，如gender1或genderFemale
      for(cat_var in categorical_vars) {
        if(grepl(paste0("^", cat_var), feature)) {
          if(!cat_var %in% original_selected_features) {
            original_selected_features <- c(original_selected_features, cat_var)
          }
          next
        }
      }
      
      # 如果不是哑变量，则是连续变量，直接添加
      if(!any(sapply(categorical_vars, function(cv) grepl(paste0("^", cv), feature)))) {
        original_selected_features <- c(original_selected_features, feature)
      }
    }
  }
  # 确保所有选择的特征都在数据集中
  available_original_features <- intersect(original_selected_features, colnames(day_data_imputed))
  
  # 准备用于训练的数据集
  model_data <- day_data_imputed[, c(available_original_features, "vision_improvement_1m")]
  
  # 设置5折交叉验证
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    savePredictions = "final",
    summaryFunction = defaultSummary
  )
  
  # 检查数据集大小，如果样本量太小，则减少折数
  if(nrow(model_data) < 10) {
    warning(paste(day_label, "数据集样本量太小 (n =", nrow(model_data), ")，改用留一法交叉验证"))
    ctrl <- trainControl(
      method = "LOOCV",
      savePredictions = "final",
      summaryFunction = defaultSummary
    )
  }
  
  # Train models
  set.seed(123)
  
  # Linear Regression
  lm_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "lm",
    trControl = ctrl,
    metric = "RMSE"
  )
  
  # Random Forest
  rf_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "rf",
    trControl = ctrl,
    metric = "RMSE"
  )
  
  # XGBoost
  xgb_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "xgbTree",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 5
  )
  
  # LASSO
  lasso_grid <- expand.grid(
    alpha = 1,
    lambda = seq(0.001, 0.1, length.out = 10)
  )
  
  lasso_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid = lasso_grid,
    metric = "RMSE"
  )
  
  # Elastic Net
  elastic_net_grid <- expand.grid(
    alpha = c(0.2, 0.5, 0.8),
    lambda = seq(0.01, 0.1, length.out = 5)
  )
  
  enet_model <- train(
    vision_improvement_1m ~ .,
    data = model_data,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid = elastic_net_grid,
    metric = "RMSE"
  )
  
  # Calculate performance metrics for each model
  calculate_performance <- function(model) {
    model$pred %>%
      group_by(Resample) %>%
      summarise(
        RMSE = sqrt(mean((obs - pred)^2)),
        MAE = mean(abs(obs - pred)),
        R2 = cor(obs, pred)^2
      ) %>%
      summarise(
        RMSE_mean = mean(RMSE),
        RMSE_sd = sd(RMSE),
        MAE_mean = mean(MAE),
        MAE_sd = sd(MAE),
        R2_mean = mean(R2),
        R2_sd = sd(R2)
      )
  }
  
  # Get performance for each model
  lm_performance <- calculate_performance(lm_model)
  rf_performance <- calculate_performance(rf_model)
  xgb_performance <- calculate_performance(xgb_model)
  lasso_performance <- calculate_performance(lasso_model)
  enet_performance <- calculate_performance(enet_model)
  
  # Combine results
  results <- bind_rows(
    lm_performance %>% mutate(Model = "Linear Regression"),
    rf_performance %>% mutate(Model = "Random Forest"),
    xgb_performance %>% mutate(Model = "XGBoost"),
    lasso_performance %>% mutate(Model = "LASSO"),
    enet_performance %>% mutate(Model = "Elastic Net")
  ) %>%
    mutate(Day = day_label)
  
  # Optional: Create classification models as well if needed
  # (code omitted for brevity, but could follow the same pattern as in your original script)
  
  return(results)
}

# Main analysis
# Load all day datasets - extended to match the new day range (-7 to 29)
day_files <- list(
  list(file = "day_-7_data.csv", label = "Day -7"),
  list(file = "day_-6_data.csv", label = "Day -6"),
  list(file = "day_-5_data.csv", label = "Day -5"),
  list(file = "day_-4_data.csv", label = "Day -4"),
  list(file = "day_-3_data.csv", label = "Day -3"),
  list(file = "day_-2_data.csv", label = "Day -2"),
  list(file = "day_-1_data.csv", label = "Day -1"),
  list(file = "day_0_data.csv", label = "Day 0"),  # 添加手术当天
  list(file = "day_1_data.csv", label = "Day 1"),
  list(file = "day_2_data.csv", label = "Day 2"),
  list(file = "day_3_data.csv", label = "Day 3"),
  list(file = "day_4_data.csv", label = "Day 4"),
  list(file = "day_5_data.csv", label = "Day 5"),
  list(file = "day_6_data.csv", label = "Day 6"),
  list(file = "day_7_data.csv", label = "Day 7"),
  list(file = "day_8_data.csv", label = "Day 8"),
  list(file = "day_9_data.csv", label = "Day 9"),
  list(file = "day_10_data.csv", label = "Day 10"),
  list(file = "day_11_data.csv", label = "Day 11"),
  list(file = "day_12_data.csv", label = "Day 12"),
  list(file = "day_13_data.csv", label = "Day 13"),
  list(file = "day_14_data.csv", label = "Day 14"),
  list(file = "day_15_data.csv", label = "Day 15"),
  list(file = "day_16_data.csv", label = "Day 16"),
  list(file = "day_17_data.csv", label = "Day 17"),
  list(file = "day_18_data.csv", label = "Day 18"),
  list(file = "day_19_data.csv", label = "Day 19"),
  list(file = "day_20_data.csv", label = "Day 20"),
  list(file = "day_21_data.csv", label = "Day 21"),
  list(file = "day_22_data.csv", label = "Day 22"),
  list(file = "day_23_data.csv", label = "Day 23"),
  list(file = "day_24_data.csv", label = "Day 24"),
  list(file = "day_25_data.csv", label = "Day 25"),
  list(file = "day_26_data.csv", label = "Day 26"),
  list(file = "day_27_data.csv", label = "Day 27"),
  list(file = "day_28_data.csv", label = "Day 28"),
  list(file = "day_29_data.csv", label = "Day 29")
)

# Process each day's data
all_results <- data.frame()

for (day_info in day_files) {
  cat("Processing", day_info$label, "...\n")
  
  # Load data for this day
  day_data <- NULL
  
  tryCatch({
    day_data <- read_csv(paste0("3_data_analysis/3_prediction_modeling/1m_prediction/daily_data/", day_info$file), 
                         show_col_types = FALSE)
    
    # 处理分类变量
    # 明确指定哪些是分类变量
    categorical_vars <- c("gender", "dm_2", "cataract_2", "hypertension_2", "vision_improved_1m","season")
    
    # 将分类变量转换为因子
    for(col in intersect(categorical_vars, names(day_data))) {
      day_data[[col]] <- as.factor(day_data[[col]])
    }
    
    # 处理其他可能的字符型列
    char_cols <- sapply(day_data, is.character)
    if(any(char_cols)) {
      for(col in names(day_data)[char_cols]) {
        if(col != "subject_id" && col != "vision_improved_factor_1m" && 
           !(col %in% categorical_vars)) {
          # 尝试转换为数值
          day_data[[col]] <- as.numeric(day_data[[col]])
        }
      }
    }
    
  }, error = function(e) {
    cat("无法加载", day_info$file, ":", e$message, "\n")
    return(NULL)
  })
  
  # 检查是否成功加载数据
  if(is.null(day_data)) {
    cat("跳过", day_info$label, "因为无法加载数据\n")
    next
  }
  
  # 添加这段检查代码 - 这是缺少的部分
  if(!"vision_improvement_1m" %in% names(day_data)) {
    cat("跳过", day_info$label, "因为缺少目标变量 vision_improvement_1m\n")
    next
  }
  
  # 检查数据集大小
  if(nrow(day_data) < 5) {
    cat("跳过", day_info$label, "因为样本量太小 (n =", nrow(day_data), ")\n")
    next
  }
  
  # 计算有效样本量(非NA的vision_improvement_1m)
  valid_samples <- sum(!is.na(day_data$vision_improvement_1m))
  if(valid_samples < 5) {
    cat("跳过", day_info$label, "因为有效样本量太小 (有效样本量 =", valid_samples, ")\n")
    next
  }
  
  # 分析数据
  tryCatch({
    day_results <- analyze_day_data(day_data, day_info$label)
    
    # 加入总体结果
    if(!is.null(day_results) && nrow(day_results) > 0) {
      all_results <- bind_rows(all_results, day_results)
      cat(day_info$label, "处理完成\n")
    } else {
      cat(day_info$label, "没有产生有效结果\n")
    }
  }, error = function(e) {
    cat("分析", day_info$label, "时出错:", e$message, "\n")
  })
}


#####save results
# 设置输出目录
output_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/daily_data/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 首先定义日期顺序
day_order <- c(
  "Day -7", "Day -6", "Day -5", "Day -4", "Day -3", "Day -2", "Day -1", "Day 0",
  "Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7",
  "Day 8", "Day 9", "Day 10", "Day 11", "Day 12", "Day 13", "Day 14",
  "Day 15", "Day 16", "Day 17", "Day 18", "Day 19", "Day 20", "Day 21",
  "Day 22", "Day 23", "Day 24", "Day 25", "Day 26", "Day 27", "Day 28", "Day 29"
)

# 确保Day是有序因子
all_results$Day <- factor(all_results$Day, levels = day_order, ordered = TRUE)

# 保存全部结果数据
write.csv(all_results, file.path(output_dir, "all_model_results.csv"), row.names = FALSE)

# 创建并保存最佳天数结果
best_days <- all_results %>%
  group_by(Model) %>%
  slice_max(order_by = R2_mean, n = 1) %>%
  dplyr::select(Model, Day, R2_mean, R2_sd, RMSE_mean) %>%
  arrange(desc(R2_mean))

write.csv(best_days, file.path(output_dir, "best_days_by_model.csv"), row.names = FALSE)

# 保存总体模型表现
best_models <- all_results %>%
  group_by(Model) %>%
  summarise(
    Avg_R2 = mean(R2_mean),
    Max_R2 = max(R2_mean),
    Min_R2 = min(R2_mean),
    Std_R2 = sd(R2_mean)
  ) %>%
  arrange(desc(Avg_R2))

write.csv(best_models, file.path(output_dir, "overall_model_performance.csv"), row.names = FALSE)

# 创建并保存线图
line_plot <- ggplot(all_results, aes(x = Day, y = R2_mean, color = Model, group = Model)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = R2_mean - R2_sd, ymax = R2_mean + R2_sd), width = 0.2) +
  labs(
    title = "1-Month Vision Improvement Prediction Performance by Day",
    subtitle = "Comparing different models across pre- and post-operative days",
    x = "Days Relative to Surgery",
    y = "R² (coefficient of determination)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次网格线
  )
line_plot
# 保存线图
ggsave(file.path(output_dir, "prediction_performance_line.pdf"), line_plot, width = 10, height = 7)

# 创建并保存热图
heatmap_data <- all_results %>%
  dplyr::select(Day, Model, R2_mean) %>%
  pivot_wider(names_from = Model, values_from = R2_mean)

heatmap_plot <- ggplot(heatmap_data, aes(x = Day)) +
  geom_tile(aes(y = "Linear Regression", fill = `Linear Regression`)) +
  geom_tile(aes(y = "Random Forest", fill = `Random Forest`)) +
  geom_tile(aes(y = "XGBoost", fill = `XGBoost`)) +
  geom_tile(aes(y = "LASSO", fill = `LASSO`)) +
  geom_tile(aes(y = "Elastic Net", fill = `Elastic Net`)) +
  geom_text(aes(y = "Linear Regression", label = round(`Linear Regression`, 3)), color = "white") +
  geom_text(aes(y = "Random Forest", label = round(`Random Forest`, 3)), color = "white") +
  geom_text(aes(y = "XGBoost", label = round(`XGBoost`, 3)), color = "white") +
  geom_text(aes(y = "LASSO", label = round(`LASSO`, 3)), color = "white") +
  geom_text(aes(y = "Elastic Net", label = round(`Elastic Net`, 3)), color = "white") +
  scale_fill_gradient(low = "pink", high = "darkred") +
  labs(
    title = "1-Month Prediction R² Heatmap by Day and Model",
    x = "Days Relative to Surgery",
    y = "Model",
    fill = "R²"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

heatmap_plot 

# 保存热图
ggsave(file.path(output_dir, "r2_heatmap.pdf"), heatmap_plot, width = 12, height = 6)

# 创建并保存摘要条形图
summary_plot <- ggplot(best_days, aes(x = reorder(Model, -R2_mean), y = R2_mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = R2_mean - R2_sd, ymax = R2_mean + R2_sd), width = 0.2) +
  geom_text(aes(label = paste0(round(R2_mean, 3), " (", Day, ")")), vjust = -0.5) +
  labs(
    title = "Best 1-Month Prediction Performance by Model",
    subtitle = "Showing highest R² value and corresponding day",
    x = "Model",
    y = "Maximum R²"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_plot

# 保存摘要图
ggsave(file.path(output_dir, "best_model_summary.pdf"), summary_plot, width = 8, height = 6)

# 创建基于1周视力预测1月视力的基准模型
# 收集所有包含post_vision_1w作为特征的有效天数数据
days_with_1w <- all_results %>%
  filter(R2_mean > 0) %>%  # 只选择有效的天数
  dplyr::select(Day, Model, R2_mean) %>%
  group_by(Day) %>%
  summarise(max_R2 = max(R2_mean)) %>%
  arrange(desc(max_R2))

# 导出1周视力对1月预测的贡献分析
write.csv(days_with_1w, file.path(output_dir, "one_week_vision_contribution.csv"), row.names = FALSE)
