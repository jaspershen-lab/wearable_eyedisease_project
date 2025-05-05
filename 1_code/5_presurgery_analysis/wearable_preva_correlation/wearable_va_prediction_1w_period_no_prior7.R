# 术后视力改善预测模型 - 按时间窗口分析
# 使用不同时间窗口的可穿戴设备数据预测术后1周视力改善值

library(tidyverse)
library(caret)
library(gridExtra)
library(pROC)
library(glmnet)
library(ggplot2)

# 设置工作目录
setwd(get_project_wd())
rm(list = ls())

#-------------------------------
# 1. 定义标准时间窗口
#-------------------------------

# 定义标准时间窗口
time_windows <- c(
  "pre_7d_all",   # 术前7天及以前
  "pre_3d_7d",    # 术前3-7天
  "pre_3d",       # 术前3天
  "post_3d",      # 术后3天
  "post_4d_6d",   # 术后4-6天
  "post_6d_all"   # 术后7天
)

# 创建输出目录
output_dir <- "3_data_analysis/5_presurgery_analysis/vision_prediction/results_by_timewindow/1w"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------
# 2. 加载数据
#-------------------------------

# 加载基础数据 - 使用日期信息
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1m_prediction/daily_data/all_days_combined_data.csv")

# 检查combined_data是否包含vision_improvement_1w
if(!"vision_improvement_1w" %in% colnames(combined_data)) {
  stop("vision_improvement_1w变量不在combined_data中，请检查数据")
}

# 检查combined_data是否包含pre_vision
if(!"pre_vision" %in% colnames(combined_data)) {
  stop("pre_vision变量不在combined_data中，请检查数据")
}

#-------------------------------
# 3. 按时间窗口准备可穿戴设备预测变量
#-------------------------------

# 创建一个函数，根据时间窗口筛选数据并计算均值
prepare_timewindow_data <- function(data, time_window) {
  # 根据时间窗口定义日期范围
  date_range <- switch(time_window,
                       "pre_7d_all" = c(-6, -1),      # 术前7天以前(假设最多考虑术前30天)
                       "pre_3d_7d" = c(-6, -4),        # 术前3-7天
                       "pre_3d" = c(-3, -1),           # 术前3天
                       "post_3d" = c(1, 3),            # 术后3天
                       "post_4d_6d" = c(4, 6),       # 术后4-6天
                       "post_6d_all" = c(1, 6),    # 术后1-6天
                       c(-30, 30)                       # 默认范围，应该不会用到
  )
  
  # 筛选指定时间窗口的数据
  filtered_data <- data %>%
    filter(as.numeric(day) >= date_range[1] & as.numeric(day) <= date_range[2])
  
  # 按受试者ID分组并计算均值
  aggregated_data <- filtered_data %>%
    group_by(subject_id) %>%
    summarise(
      # 选择关键均值指标
      mean_rhr_1 = mean(mean_rhr_1, na.rm = TRUE),     # 平均静息心率
      mean_bo = mean(mean_bo, na.rm = TRUE),           # 平均血氧饱和度
      steps_total = mean(steps_total, na.rm = TRUE),   # 平均总步数
      total_sleep = mean(total_sleep, na.rm = TRUE),   # 平均总睡眠时间
      
      # 保留其他人口统计变量
      dm_2 = dplyr::first(dm_2),                    # 糖尿病状态
      age = dplyr::first(age),                      # 年龄
      gender = dplyr::first(gender),                # 性别
      bmi = dplyr::first(bmi),                      # BMI
      hypertension_2 = dplyr::first(hypertension_2), # 高血压状态
      
      # 保留视力相关变量
      pre_vision = dplyr::first(pre_vision),            # 术前视力
      vision_improvement_1w = dplyr::first(vision_improvement_1w)  # 术后1周视力改善（目标变量）
    )
  
  # 中位数填补缺失值
  for(col in c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep")) {
    if(sum(is.na(aggregated_data[[col]])) > 0) {
      aggregated_data[[col]][is.na(aggregated_data[[col]])] <- median(aggregated_data[[col]], na.rm = TRUE)
    }
  }
  
  return(aggregated_data)
}

#-------------------------------
# 4. 定义要分析的变量
#-------------------------------

# 目标变量 - 术后1周视力改善
target_var <- "vision_improvement_1w"

# 定义可穿戴设备的预测变量
wearable_vars <- list(
  "mean_rhr_1",    # 静息心率
  "mean_bo",       # 血氧
  "steps_total",   # 步数
  "total_sleep"    # 睡眠时间
)

# 定义其他预测变量
additional_vars <- list(
  "pre_vision"      # 术前视力
)

#-------------------------------
# 5. 创建一个函数来存储每个时间窗口的模型性能
#-------------------------------

# 预测模型构建函数 - 适用于视力改善
build_prediction_models <- function(model_data, target_var, additional_vars, time_window = "unknown") {
  # 创建数据框来存储每个预测器组合的模型性能指标
  performance_df <- data.frame(
    time_window = character(),
    target_var = character(),
    predictor = character(),
    r_squared = numeric(),
    adj_r_squared = numeric(),
    cv_r_squared = numeric(),
    cv_rmse = numeric(),
    sample_size = numeric(),
    predictor_coef = numeric(),
    predictor_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 检查目标变量是否存在
  if(!(target_var %in% names(model_data))) {
    cat("目标变量不存在:", target_var, "\n")
    return(performance_df)
  }
  
  cat("\n===== 构建视力改善预测模型:", target_var, "- 时间窗口:", time_window, "=====\n")
  
  # 准备基本模型数据 - 移除含有NA的行
  base_vars <- c(target_var, unlist(wearable_vars), unlist(additional_vars), "age", "gender", "bmi", "dm_2")
  base_data <- model_data %>%
    dplyr::select(all_of(base_vars)) %>%
    na.omit()
  
  if(nrow(base_data) < 10) {  # 确保有足够的样本
    cat("样本量不足:", target_var, "\n")
    return(performance_df)
  }
  
  cat("样本量:", nrow(base_data), "\n")
  
  # 定义交叉验证函数
  cv_model_performance <- function(model_formula, model_data) {
    # 创建训练控制对象
    train_control <- trainControl(
      method = "cv",
      number = 5,
      savePredictions = "final"
    )
    
    # 使用caret的train函数进行交叉验证
    set.seed(123)
    cv_model <- tryCatch({
      train(
        model_formula,
        data = model_data,
        method = "lm",
        trControl = train_control
      )
    }, error = function(e) {
      cat("交叉验证错误:", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(cv_model)) {
      # 如果交叉验证失败，使用常规线性模型
      lm_model <- lm(model_formula, data = model_data)
      return(list(
        r_squared = summary(lm_model)$r.squared,
        rmse = sqrt(mean(lm_model$residuals^2))
      ))
    }
    
    # 返回交叉验证性能指标
    return(list(
      r_squared = cv_model$results$Rsquared,
      rmse = cv_model$results$RMSE
    ))
  }
  
  # 合并可穿戴设备变量为公式字符串
  wearable_formula_str <- paste(unlist(wearable_vars), collapse = " + ")
  additional_formula_str <- paste(unlist(additional_vars), collapse = " + ")
  
  # 1. 所有预测变量模型（包括pre_vision）
  all_formula <- as.formula(paste(
    target_var, "~", 
    wearable_formula_str, " + ",
    additional_formula_str, " + ",
    "age + gender + bmi + dm_2"
  ))
  
  all_model <- tryCatch({
    lm(all_formula, data = base_data)
  }, error = function(e) {
    cat("构建全模型错误:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(all_model)) {
    all_summary <- summary(all_model)
    all_cv <- cv_model_performance(all_formula, base_data)
    
    # 添加到性能数据框
    performance_df <- rbind(performance_df, data.frame(
      time_window = time_window,
      target_var = target_var,
      predictor = "All Predictors",
      r_squared = all_summary$r.squared,
      adj_r_squared = all_summary$adj.r.squared,
      cv_r_squared = all_cv$r_squared,
      cv_rmse = all_cv$rmse,
      sample_size = nrow(base_data),
      predictor_coef = NA,  # 对于全部预测变量模型，不记录单个系数
      predictor_p_value = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # 2. 对每个可穿戴设备变量构建单独的模型
  for(wearable_var in wearable_vars) {
    # 创建模型公式，包含pre_vision
    var_formula <- as.formula(paste(
      target_var, "~", 
      wearable_var, " + ",
      additional_formula_str, " + ",
      "age + gender + bmi + dm_2"
    ))
    
    # 构建模型
    var_model <- tryCatch({
      lm(var_formula, data = base_data)
    }, error = function(e) {
      cat("构建单变量模型错误:", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(var_model)) {
      var_summary <- summary(var_model)
      var_cv <- cv_model_performance(var_formula, base_data)
      
      # 获取预测变量的系数和p值
      coef_data <- coef(var_summary)
      var_coef <- coef_data[wearable_var, 1]  # 系数
      var_p_value <- coef_data[wearable_var, 4]  # p值
      
      # 添加到性能数据框
      performance_df <- rbind(performance_df, data.frame(
        time_window = time_window,
        target_var = target_var,
        predictor = wearable_var,
        r_squared = var_summary$r.squared,
        adj_r_squared = var_summary$adj.r.squared,
        cv_r_squared = var_cv$r_squared,
        cv_rmse = var_cv$rmse,
        sample_size = nrow(base_data),
        predictor_coef = var_coef,
        predictor_p_value = var_p_value,
        stringsAsFactors = FALSE
      ))
      
      # 打印单个预测变量模型信息
      cat("\n--- 单一预测变量:", wearable_var, "---\n")
      cat("R²:", var_summary$r.squared, "\n")
      cat("调整后R²:", var_summary$adj.r.squared, "\n")
      cat("交叉验证R²:", var_cv$r_squared, "\n")
      cat("交叉验证RMSE:", var_cv$rmse, "\n")
      cat(wearable_var, "系数:", var_coef, "p值:", var_p_value, "\n")
    }
  }
  
  # 3. 只使用pre_vision变量构建模型（没有可穿戴设备变量）
  pre_vision_formula <- as.formula(paste(
    target_var, "~", 
    additional_formula_str, " + ",
    "age + gender + bmi + dm_2"
  ))
  
  pre_vision_model <- tryCatch({
    lm(pre_vision_formula, data = base_data)
  }, error = function(e) {
    cat("构建术前视力模型错误:", e$message, "\n")
    return(NULL)
  })
  
  if(!is.null(pre_vision_model)) {
    pre_vision_summary <- summary(pre_vision_model)
    pre_vision_cv <- cv_model_performance(pre_vision_formula, base_data)
    
    # 获取pre_vision的系数和p值
    coef_data <- coef(pre_vision_summary)
    pre_vision_coef <- coef_data["pre_vision", 1]  # 系数
    pre_vision_p_value <- coef_data["pre_vision", 4]  # p值
    
    # 添加到性能数据框
    performance_df <- rbind(performance_df, data.frame(
      time_window = time_window,
      target_var = target_var,
      predictor = "pre_vision",
      r_squared = pre_vision_summary$r.squared,
      adj_r_squared = pre_vision_summary$adj.r.squared,
      cv_r_squared = pre_vision_cv$r_squared,
      cv_rmse = pre_vision_cv$rmse,
      sample_size = nrow(base_data),
      predictor_coef = pre_vision_coef,
      predictor_p_value = pre_vision_p_value,
      stringsAsFactors = FALSE
    ))
    
    # 打印术前视力模型信息
    cat("\n--- 术前视力预测变量: pre_vision ---\n")
    cat("R²:", pre_vision_summary$r.squared, "\n")
    cat("调整后R²:", pre_vision_summary$adj.r.squared, "\n")
    cat("交叉验证R²:", pre_vision_cv$r_squared, "\n")
    cat("交叉验证RMSE:", pre_vision_cv$rmse, "\n")
    cat("pre_vision 系数:", pre_vision_coef, "p值:", pre_vision_p_value, "\n")
  }
  
  return(performance_df)
}

#-------------------------------
# 6. 创建存储所有模型性能的数据框
#-------------------------------

# 为视力改善数据创建存储性能的数据框
all_vision_performance <- data.frame()

#-------------------------------
# 7. 对每个时间窗口运行预测模型
#-------------------------------

cat("\n开始为每个时间窗口构建视力改善预测模型...\n")

for(time_window in time_windows) {
  cat("\n=== 处理时间窗口:", time_window, "===\n")
  
  # 获取当前时间窗口的可穿戴设备数据
  current_window_data <- prepare_timewindow_data(combined_data, time_window)
  
  # 构建视力改善预测模型
  cat("\n开始构建视力改善预测模型 - 时间窗口:", time_window, "...\n")
  vision_performance <- tryCatch({
    build_prediction_models(current_window_data, target_var, additional_vars, time_window)
  }, error = function(e) {
    cat("Error in vision model:", e$message, "\n")
    return(data.frame())
  })
  
  # 将当前窗口的性能结果追加到总结果
  all_vision_performance <- rbind(all_vision_performance, vision_performance)
}

#-------------------------------
# 8. 保存模型性能结果
#-------------------------------

# 保存到CSV文件
write.csv(all_vision_performance, file.path(output_dir, "vision_improvement_prediction_by_timewindow.csv"), row.names = FALSE)
cat("\n所有模型性能结果已保存到:", file.path(output_dir, "vision_improvement_prediction_by_timewindow.csv"), "\n")

#-------------------------------
# 9. 创建性能变化趋势图表 (折线图)
#-------------------------------

# 按预测变量绘制不同时间窗口的模型性能折线图
create_trend_plots <- function(performance_data) {
  # 准备绘图数据 - 按预测变量分组
  plot_data <- performance_data %>%
    select(time_window, predictor, cv_r_squared) %>%
    # 确保时间窗口按照定义的顺序排序
    mutate(time_window = factor(time_window, levels = time_windows))
  
  # 为不同预测变量定义颜色
  color_values <- c(
    "All Predictors" = "#513a52", 
    "mean_rhr_1" = "#82a1bf", 
    "mean_bo" = "#faaa93", 
    "steps_total" = "#feefc4",
    "total_sleep" = "#7fb77e"
  )
  
  # 创建折线图
  p <- ggplot(plot_data, aes(x = time_window, y = cv_r_squared, color = predictor, group = predictor)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_manual(
      values = color_values,
      labels = c(
        "All Predictors" = "All Predictors", 
        "mean_rhr_1" = " RestingHeart Rate", 
        "mean_bo" = "Blood Oxygen", 
        "steps_total" = "Total Steps",
        "total_sleep" = "Sleep Duration"
      )
    ) +
    labs(
      title = "Prediction of Vision Improvement (1 Week)",
      x = "Time Window",
      y = "Cross-validated R²",
      color = "Predictors"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # 保存图表
  plot_filename <- "vision_improvement_performance_trend.pdf"
  ggsave(file.path(output_dir, plot_filename), p, width = 10, height = 6)
  
  # 同时保存PNG格式
  plot_filename_png <- "vision_improvement_performance_trend.png"
  ggsave(file.path(output_dir, plot_filename_png), p, width = 10, height = 6, dpi = 300)
  
  cat("已保存视力改善预测趋势图\n")
}

# 绘制视力改善预测模型性能趋势图
create_trend_plots(all_vision_performance)

#-------------------------------
# 10. 创建热图显示不同时间窗口和预测变量的性能
#-------------------------------

create_heatmap <- function(performance_data) {
  # 创建热图数据 - 排除"All Predictors"和"pre_vision"
  heatmap_data <- performance_data %>%
    filter(!predictor %in% c("All Predictors", "pre_vision")) %>%
    select(time_window, predictor, cv_r_squared) %>%
    # 确保时间窗口按照定义的顺序排序
    mutate(
      time_window = factor(time_window, levels = time_windows),
      predictor = factor(predictor, levels = c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep"),
                         labels = c("Resting Heart Rate", "Blood Oxygen", "Total Steps", "Sleep Duration"))
    )
  
  # 创建热图
  p_heatmap <- ggplot(heatmap_data, aes(x = time_window, y = predictor, fill = cv_r_squared)) +
    geom_tile() +
    scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                         midpoint = median(heatmap_data$cv_r_squared, na.rm = TRUE)) +
    labs(
      title = "Predictive Power for Vision Improvement (1 Week)",
      x = "Time Window",
      y = "Predictor",
      fill = "Cross-validated R²"
    ) + 
    scale_y_discrete(
      labels = c(
        "Heart Rate" = "Resting Heart Rate", 
        "Blood Oxygen" = "Blood Oxygen", 
        "Total Steps" = "Total Steps",
        "Sleep Duration" = "Sleep Duration"
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # 保存热图
  heatmap_filename <- "vision_improvement_heatmap.pdf"
  ggsave(file.path(output_dir, heatmap_filename), p_heatmap, width = 10, height = 6)
  
  # 同时保存PNG格式
  heatmap_filename_png <- "vision_improvement_heatmap.png"
  ggsave(file.path(output_dir, heatmap_filename_png), p_heatmap, width = 10, height = 6, dpi = 300)
  
  cat("已保存视力改善预测热图\n")
}

# 绘制视力改善预测模型热图
create_heatmap(all_vision_performance)

#-------------------------------
# 11. 创建预测变量重要性图表
#-------------------------------

# 分析不同预测变量的重要性（基于p值和R²）
create_importance_plot <- function(performance_data) {
  # 筛选单个预测变量的数据
  var_data <- performance_data %>%
    filter(!predictor %in% c("All Predictors"))
  
  # 对每个预测变量，计算显著结果（p<0.05）的比例
  significance_data <- var_data %>%
    group_by(predictor) %>%
    summarise(
      total_windows = n(),
      significant_windows = sum(predictor_p_value < 0.05, na.rm = TRUE),
      significance_rate = significant_windows / total_windows * 100,
      mean_rsquared = mean(cv_r_squared, na.rm = TRUE),
      mean_coef_abs = mean(abs(predictor_coef), na.rm = TRUE)
    ) %>%
    arrange(desc(significance_rate))
  
  # 创建水平条形图 - 显著性比例
  p1 <- ggplot(significance_data, aes(x = reorder(predictor, significance_rate), y = significance_rate)) +
    geom_bar(stat = "identity", fill = "#4292c6") +
    geom_text(aes(label = sprintf("%.1f%%", significance_rate)), hjust = -0.2) +
    labs(
      title = "Predictor significance for 1-week vision improvement",
      subtitle = "Percentage of time windows with statistically significant effect (p<0.05)",
      x = "",
      y = "Percentage of significant windows (%)"
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 12),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.2)))
  
  # 创建水平条形图 - 平均R²
  p2 <- ggplot(significance_data, aes(x = reorder(predictor, mean_rsquared), y = mean_rsquared)) +
    geom_bar(stat = "identity", fill = "#2ca25f") +
    geom_text(aes(label = sprintf("%.3f", mean_rsquared)), hjust = -0.2) +
    labs(
      title = "Average predictive power for 1-week vision improvement",
      subtitle = "Mean cross-validated R² across all time windows",
      x = "",
      y = "Mean cross-validated R²"
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 12),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  # 保存重要性图表
  ggsave(file.path(output_dir, "predictor_significance_rate.pdf"), p1, width = 8, height = 6)
  ggsave(file.path(output_dir, "predictor_mean_rsquared.pdf"), p2, width = 8, height = 6)
  
  # 同时保存PNG格式
  ggsave(file.path(output_dir, "predictor_significance_rate.png"), p1, width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "predictor_mean_rsquared.png"), p2, width = 8, height = 6, dpi = 300)
  
  # 保存重要性数据
  write.csv(significance_data, file.path(output_dir, "predictor_importance.csv"), row.names = FALSE)
  
  cat("已保存预测变量重要性图表和数据\n")
}

# 创建预测变量重要性图表
create_importance_plot(all_vision_performance)
