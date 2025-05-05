# 每日数据视力改善预测模型
# 使用每天的可穿戴设备数据预测术后1周视力改善值

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
# 1. 定义每日时间点
#-------------------------------

# 定义每日时间点（从术前7天到术后6天）
daily_timepoints <- c(
   "-6", "-5", "-4", "-3", "-2", "-1",  # 术前7天到术前1天
  "1", "2", "3", "4", "5", "6"               # 术后1天到术后6天
)

# 创建输出目录
output_dir <- "3_data_analysis/5_presurgery_analysis/vision_prediction/results_by_day/1w/no_prior7"
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

#-------------------------------
# 3. 按每日时间点准备可穿戴设备预测变量
#-------------------------------

# 创建一个函数，根据特定日期筛选数据
prepare_daily_data <- function(data, day_point) {
  # 筛选指定日期的数据
  filtered_data <- data %>%
    filter(as.numeric(day) == as.numeric(day_point))
  
  # 如果筛选后数据量很少，不继续处理
  if(nrow(filtered_data) < 5) {
    cat("日期", day_point, "的样本数量不足:", nrow(filtered_data), "\n")
    return(NULL)
  }
  
  # 数据预处理
  processed_data <- filtered_data %>%
    select(
      subject_id,
      # 选择关键指标
      mean_rhr_1,     # 平均静息心率
      mean_bo,        # 平均血氧饱和度
      steps_total,    # 平均总步数
      total_sleep,    # 平均总睡眠时间
      
      # 保留其他人口统计变量
      dm_2,           # 糖尿病状态
      age,            # 年龄
      gender,         # 性别
      bmi,            # BMI
      hypertension_2, # 高血压状态
      pre_vision,
      
      # 保留视力改善值
      vision_improvement_1w  # 术后1周视力改善值
    )
  
  # 中位数填补缺失值
  for(col in c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep")) {
    if(sum(is.na(processed_data[[col]])) > 0) {
      processed_data[[col]][is.na(processed_data[[col]])] <- median(processed_data[[col]], na.rm = TRUE)
    }
  }
  
  return(processed_data)
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

#-------------------------------
# 5. 创建一个函数来存储每个日期的模型性能
#-------------------------------

# 预测模型构建函数 - 适用于视力改善
build_prediction_models <- function(model_data, target_var, day_point = "unknown") {
  # 创建数据框来存储每个预测器组合的模型性能指标
  performance_df <- data.frame(
    day_point = character(),
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
  
  cat("\n===== 构建视力改善预测模型:", target_var, "- 日期:", day_point, "=====\n")
  
  # 准备基本模型数据 - 移除含有NA的行
  base_vars <- c(target_var, unlist(wearable_vars), "age", "gender", "bmi", "dm_2","pre_vision")
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
  
  # 合并前三个变量为公式字符串
  wearable_formula_str <- paste(unlist(wearable_vars), collapse = " + ")
  
  # 1. 所有预测变量模型
  all_formula <- as.formula(paste(
    target_var, "~", 
    wearable_formula_str, " + ",
    "age + gender + bmi + dm_2+pre_vision"
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
      day_point = day_point,
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
    # 创建模型公式
    var_formula <- as.formula(paste(
      target_var, "~", 
      wearable_var, "+", 
      "age + gender + bmi + dm_2+pre_vision"
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
        day_point = day_point,
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
  
  return(performance_df)
}

#-------------------------------
# 6. 创建存储所有模型性能的数据框
#-------------------------------

# 为视力改善数据创建存储性能的数据框
all_vision_performance <- data.frame()

# 创建存储每个日期样本量的数据框
sample_size_data <- data.frame(
  day_point = character(),
  sample_size = numeric(),
  stringsAsFactors = FALSE
)

#-------------------------------
# 7. 对每个日期运行预测模型
#-------------------------------

cat("\n开始为每个日期构建视力改善预测模型...\n")

for(day_point in daily_timepoints) {
  cat("\n=== 处理日期:", day_point, "===\n")
  
  # 获取当前日期的可穿戴设备数据
  current_day_data <- prepare_daily_data(combined_data, day_point)
  
  # 如果该日期的数据不足，跳过
  if(is.null(current_day_data)) {
    cat("跳过日期", day_point, "，因为数据不足\n")
    # 记录样本量为0
    sample_size_data <- rbind(sample_size_data, data.frame(
      day_point = day_point,
      sample_size = 0,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  # 记录该日期的样本量
  sample_size_data <- rbind(sample_size_data, data.frame(
    day_point = day_point,
    sample_size = nrow(current_day_data),
    stringsAsFactors = FALSE
  ))
  
  # 构建视力改善预测模型
  cat("\n开始构建视力改善预测模型 - 日期:", day_point, "...\n")
  vision_performance <- tryCatch({
    build_prediction_models(current_day_data, target_var, day_point)
  }, error = function(e) {
    cat("Error in vision model:", e$message, "\n")
    return(data.frame())
  })
  
  # 将当前日期的性能结果追加到总结果
  all_vision_performance <- rbind(all_vision_performance, vision_performance)
}

#-------------------------------
# 8. 保存模型性能结果
#-------------------------------

# 保存到CSV文件
write.csv(all_vision_performance, file.path(output_dir, "vision_improvement_prediction_by_day.csv"), row.names = FALSE)
cat("\n所有模型性能结果已保存到:", file.path(output_dir, "vision_improvement_prediction_by_day.csv"), "\n")

# 保存每日样本量数据
write.csv(sample_size_data, file.path(output_dir, "daily_sample_size.csv"), row.names = FALSE)
cat("\n每日样本量数据已保存到:", file.path(output_dir, "daily_sample_size.csv"), "\n")

#-------------------------------
# 9. 创建性能变化趋势图表 (折线图)
#-------------------------------

# 按预测变量绘制不同日期的模型性能折线图
create_trend_plots <- function(performance_data) {
  # 准备绘图数据 - 按预测变量分组
  plot_data <- performance_data %>%
    select(day_point, predictor, cv_r_squared) %>%
    # 确保日期按照定义的顺序排序
    mutate(day_point = factor(day_point, levels = daily_timepoints))
  
  # 为不同预测变量定义颜色
  color_values <- c(
    "All Predictors" = "#513a52", 
    "mean_rhr_1" = "#82a1bf", 
    "mean_bo" = "#faaa93", 
    "steps_total" = "#feefc4",
    "total_sleep" = "#7fb77e"  # 添加睡眠变量的颜色
  )
  
  # 创建折线图
  p <- ggplot(plot_data, aes(x = day_point, y = cv_r_squared, color = predictor, group = predictor)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_manual(
      values = color_values,
      labels = c(
        "All Predictors" = "All Predictors", 
        "mean_rhr_1" = "Resting Heart Rate", 
        "mean_bo" = "Blood Oxygen", 
        "steps_total" = "Total Steps",
        "total_sleep" = "Sleep Duration"
      )
    ) +
    labs(
      title = "Prediction of Vision Improvement (1 Week)",
      x = "Day Relative to Surgery",
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
  plot_filename <- "vision_improvement_daily_performance_trend.pdf"
  ggsave(file.path(output_dir, plot_filename), p, width = 10, height = 6)
  
  # 同时保存PNG格式
  plot_filename_png <- "vision_improvement_daily_performance_trend.png"
  ggsave(file.path(output_dir, plot_filename_png), p, width = 10, height = 6, dpi = 300)
  
  cat("已保存视力改善预测趋势图\n")
}

# 绘制视力改善预测模型性能趋势图
create_trend_plots(all_vision_performance)

#-------------------------------
# 10. 创建热图显示不同日期和预测变量的性能
#-------------------------------

create_heatmap <- function(performance_data) {
  # 创建热图数据 - 排除"All Predictors"
  heatmap_data <- performance_data %>%
    filter(predictor != "All Predictors") %>%
    select(day_point, predictor, cv_r_squared) %>%
    # 确保日期按照定义的顺序排序
    mutate(
      day_point = factor(day_point, levels = daily_timepoints),
      predictor = factor(predictor, levels = c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep"),
                         labels = c("Resting Heart Rate", "Blood Oxygen", "Total Steps", "Sleep Duration"))
    )
  
  # 创建热图
  p_heatmap <- ggplot(heatmap_data, aes(x = day_point, y = predictor, fill = cv_r_squared)) +
    geom_tile() +
    scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                         midpoint = median(heatmap_data$cv_r_squared, na.rm = TRUE)) +
    labs(
      title = "Predictive Power for Vision Improvement (1 Week)",
      x = "Day Relative to Surgery",
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
  heatmap_filename <- "vision_improvement_daily_heatmap.pdf"
  ggsave(file.path(output_dir, heatmap_filename), p_heatmap, width = 10, height = 6)
  
  # 同时保存PNG格式
  heatmap_filename_png <- "vision_improvement_daily_heatmap.png"
  ggsave(file.path(output_dir, heatmap_filename_png), p_heatmap, width = 10, height = 6, dpi = 300)
  
  cat("已保存视力改善预测热图\n")
}

# 绘制视力改善预测模型热图
create_heatmap(all_vision_performance)

#-------------------------------
# 11. 创建每日样本量图表
#-------------------------------
# 绘制每日样本量柱状图
plot_sample_size <- function(sample_data) {
  # 准备绘图数据
  plot_data <- sample_data %>%
    mutate(
      day_point = factor(day_point, levels = daily_timepoints),
      day_label = ifelse(as.numeric(day_point) < 0, 
                         paste0("Day ", day_point), 
                         paste0("Day ", day_point))
    )
  
  # 创建柱状图
  p <- ggplot(plot_data, aes(x = day_point, y = sample_size)) +
    geom_bar(stat = "identity", width = 0.7, fill = "#4292c6") +
    geom_text(aes(label = sample_size), vjust = -0.5, size = 3.5) +
    labs(
      title = "Daily sample size distribution",
      x = "Relative date to surgery",
      y = "Sample size"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) # 为文本标签留出空间
  
  return(p)
}

# 绘制样本量图表
p_sample_size <- plot_sample_size(sample_size_data)

# 保存样本量图表
sample_size_filename <- "daily_sample_size.pdf"
ggsave(file.path(output_dir, sample_size_filename), p_sample_size, width = 10, height = 6)
