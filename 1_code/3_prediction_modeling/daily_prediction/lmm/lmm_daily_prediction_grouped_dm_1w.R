library(tidyverse)
library(lme4)
library(caret)
library(ggplot2)
setwd(get_project_wd())
rm(list = ls())

# 设置工作路径（请根据需要调整）
input_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped/dm_group"
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped/lmm_model_performance/dm_group"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 获取所有分组的日期文件列表
d_files <- list.files(input_dir, pattern = "_D\\.csv$", full.names = TRUE)
nod_files <- list.files(input_dir, pattern = "_NoD\\.csv$", full.names = TRUE)

# 输出找到的文件信息
cat("找到", length(d_files), "个糖尿病组文件和", length(nod_files), "个非糖尿病组文件\n")

# 将文件名转换为时间点列表
get_time_point <- function(file_path) {
  time_str <- gsub(".*day_(.*)_[DN].*\\.csv$", "\\1", basename(file_path))
  return(time_str)
}

time_points_d <- sapply(d_files, get_time_point)
time_points_nod <- sapply(nod_files, get_time_point)

# 获取两组共有的时间点
time_points <- intersect(time_points_d, time_points_nod)

# 对时间点进行排序（从-7到6）
time_points <- sort(as.numeric(time_points))
time_points_str <- as.character(time_points)

cat("分析将使用共有的", length(time_points), "个时间点:", paste(time_points, collapse=", "), "\n")

# 创建数据框来存储每个时间点糖尿病组的性能指标
performance_df_d <- data.frame(
  time_point = character(),
  r_squared = numeric(),
  adj_r_squared = numeric(),
  rmse = numeric(),
  hr_only_r2 = numeric(),
  bo_only_r2 = numeric(),
  steps_only_r2 = numeric(),
  all_r2 = numeric(),
  sample_size = numeric(),
  stringsAsFactors = FALSE
)

# 创建数据框来存储每个时间点非糖尿病组的性能指标
performance_df_nod <- data.frame(
  time_point = character(),
  r_squared = numeric(),
  adj_r_squared = numeric(),
  rmse = numeric(),
  hr_only_r2 = numeric(),
  bo_only_r2 = numeric(),
  steps_only_r2 = numeric(),
  all_r2 = numeric(),
  sample_size = numeric(),
  stringsAsFactors = FALSE
)

# 函数计算交叉验证的R²
cv_r2 <- function(model_formula, data) {
  # 使用trainControl进行交叉验证以避免因子水平问题
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
      data = data,
      method = "lm",
      trControl = train_control
    )
  }, error = function(e) {
    cat(paste0("    交叉验证错误: ", e$message, "\n"))
    return(NULL)
  })
  
  if(is.null(cv_model)) return(NA)
  
  # 返回交叉验证R²
  return(cv_model$results$Rsquared)
}

# 对每个时间点构建糖尿病组的模型
cat("\n开始为每个时间点构建糖尿病组预测模型...\n")

for (i in 1:length(time_points)) {
  tp <- time_points_str[i]
  cat(paste0("处理时间点: ", tp, " (糖尿病组)\n"))
  
  # 读取该时间点的糖尿病组数据
  file_path <- grep(paste0("day_", tp, "_D\\.csv$"), d_files, value = TRUE)
  
  if (length(file_path) == 0) {
    cat(paste0("  警告: 没有找到时间点 ", tp, " 的糖尿病组数据文件\n"))
    next
  }
  
  # 读取数据
  day_data <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    cat(paste0("  错误: 无法读取文件: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(day_data)) next
  
  # 检查数据是否包含需要的列
  required_cols <- c("vision_improvement", "mean_rhr_1", "mean_bo", "steps_total", 
                     "age", "gender", "bmi", "pre_vision", "season")
  missing_cols <- required_cols[!required_cols %in% names(day_data)]
  
  if (length(missing_cols) > 0) {
    cat(paste0("  警告: 数据缺少以下列: ", paste(missing_cols, collapse = ", "), "\n"))
    next
  }
  
  # 移除缺失值
  day_data <- day_data[complete.cases(day_data[, required_cols]), ]
  
  if (nrow(day_data) < 10) {
    cat(paste0("  警告: 时间点 ", tp, " 的糖尿病组有效数据少于10行，可能导致模型不可靠\n"))
    next
  }
  
  # 显示样本量
  cat(paste0("  糖尿病组样本量: ", nrow(day_data), "\n"))
  
  # 预处理分类变量
  if("season" %in% names(day_data)) {
    day_data$season <- factor(day_data$season)
  }
  if("gender" %in% names(day_data)) {
    day_data$gender <- factor(day_data$gender)
  }
  
  # 设置随机种子以确保可重复性
  set.seed(123)
  
  # 构建完整模型
  full_model <- tryCatch({
    lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
         age + gender + bmi + pre_vision + season, data = day_data)
  }, error = function(e) {
    cat(paste0("  错误: 构建完整模型失败: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(full_model)) next
  
  # 计算各模型的交叉验证R²
  cat("    计算交叉验证R²值...\n")
  
  # 处理可能的错误并提供回退方案
  tryCatch({
    # 完整模型
    full_cv_r2 <- cv_r2(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
                          age + gender + bmi + pre_vision + season, day_data)
    if(is.na(full_cv_r2)) {
      # 如果交叉验证失败，使用简单的R²作为后备
      full_cv_r2 <- summary(full_model)$r.squared
      cat("    警告: 完整模型交叉验证失败，使用R²值替代\n")
    }
    
    # HR模型
    hr_cv_r2 <- cv_r2(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, day_data)
    if(is.na(hr_cv_r2)) {
      hr_model <- lm(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, data = day_data)
      hr_cv_r2 <- summary(hr_model)$r.squared
      cat("    警告: HR模型交叉验证失败，使用R²值替代\n")
    }
    
    # BO模型
    bo_cv_r2 <- cv_r2(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, day_data)
    if(is.na(bo_cv_r2)) {
      bo_model <- lm(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, data = day_data)
      bo_cv_r2 <- summary(bo_model)$r.squared
      cat("    警告: BO模型交叉验证失败，使用R²值替代\n")
    }
    
    # 步数模型
    steps_cv_r2 <- cv_r2(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, day_data)
    if(is.na(steps_cv_r2)) {
      steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, data = day_data)
      steps_cv_r2 <- summary(steps_model)$r.squared
      cat("    警告: 步数模型交叉验证失败，使用R²值替代\n")
    }
  }, error = function(e) {
    cat(paste0("    错误: 计算交叉验证R²失败: ", e$message, "\n"))
    # 如果交叉验证完全失败，使用普通R²
    full_model <- lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
                       age + gender + bmi + pre_vision + season, data = day_data)
    full_cv_r2 <<- summary(full_model)$r.squared
    
    hr_model <- lm(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, data = day_data)
    hr_cv_r2 <<- summary(hr_model)$r.squared
    
    bo_model <- lm(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, data = day_data)
    bo_cv_r2 <<- summary(bo_model)$r.squared
    
    steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, data = day_data)
    steps_cv_r2 <<- summary(steps_model)$r.squared
    
    cat("    使用常规R²值作为替代\n")
  })
  
  # 将信息添加到糖尿病组性能数据框
  performance_df_d <- rbind(performance_df_d, data.frame(
    time_point = tp,
    r_squared = summary(full_model)$r.squared,
    adj_r_squared = summary(full_model)$adj.r.squared,
    rmse = sqrt(mean(full_model$residuals^2)),
    hr_only_r2 = hr_cv_r2,
    bo_only_r2 = bo_cv_r2,
    steps_only_r2 = steps_cv_r2,
    all_r2 = full_cv_r2,
    sample_size = nrow(day_data),
    stringsAsFactors = FALSE
  ))
  
  cat(paste0("  完成时间点 ", tp, " 的糖尿病组模型构建。 交叉验证R²: ", round(full_cv_r2, 4), "\n"))
}

# 对每个时间点构建非糖尿病组的模型
cat("\n开始为每个时间点构建非糖尿病组预测模型...\n")

for (i in 1:length(time_points)) {
  tp <- time_points_str[i]
  cat(paste0("处理时间点: ", tp, " (非糖尿病组)\n"))
  
  # 读取该时间点的非糖尿病组数据
  file_path <- grep(paste0("day_", tp, "_NoD\\.csv$"), nod_files, value = TRUE)
  
  if (length(file_path) == 0) {
    cat(paste0("  警告: 没有找到时间点 ", tp, " 的非糖尿病组数据文件\n"))
    next
  }
  
  # 读取数据
  day_data <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    cat(paste0("  错误: 无法读取文件: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(day_data)) next
  
  # 检查数据是否包含需要的列
  required_cols <- c("vision_improvement", "mean_rhr_1", "mean_bo", "steps_total", 
                     "age", "gender", "bmi", "pre_vision", "season")
  missing_cols <- required_cols[!required_cols %in% names(day_data)]
  
  if (length(missing_cols) > 0) {
    cat(paste0("  警告: 数据缺少以下列: ", paste(missing_cols, collapse = ", "), "\n"))
    next
  }
  
  # 移除缺失值
  day_data <- day_data[complete.cases(day_data[, required_cols]), ]
  
  if (nrow(day_data) < 10) {
    cat(paste0("  警告: 时间点 ", tp, " 的非糖尿病组有效数据少于10行，可能导致模型不可靠\n"))
    next
  }
  
  # 显示样本量
  cat(paste0("  非糖尿病组样本量: ", nrow(day_data), "\n"))
  
  # 预处理分类变量
  if("season" %in% names(day_data)) {
    day_data$season <- factor(day_data$season)
  }
  if("gender" %in% names(day_data)) {
    day_data$gender <- factor(day_data$gender)
  }
  
  # 设置随机种子以确保可重复性
  set.seed(123)
  
  # 构建完整模型
  full_model <- tryCatch({
    lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
         age + gender + bmi + pre_vision + season, data = day_data)
  }, error = function(e) {
    cat(paste0("  错误: 构建完整模型失败: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(full_model)) next
  
  # 计算各模型的交叉验证R²
  cat("    计算交叉验证R²值...\n")
  
  # 处理可能的错误并提供回退方案
  tryCatch({
    # 完整模型
    full_cv_r2 <- cv_r2(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
                          age + gender + bmi + pre_vision + season, day_data)
    if(is.na(full_cv_r2)) {
      # 如果交叉验证失败，使用简单的R²作为后备
      full_cv_r2 <- summary(full_model)$r.squared
      cat("    警告: 完整模型交叉验证失败，使用R²值替代\n")
    }
    
    # HR模型
    hr_cv_r2 <- cv_r2(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, day_data)
    if(is.na(hr_cv_r2)) {
      hr_model <- lm(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, data = day_data)
      hr_cv_r2 <- summary(hr_model)$r.squared
      cat("    警告: HR模型交叉验证失败，使用R²值替代\n")
    }
    
    # BO模型
    bo_cv_r2 <- cv_r2(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, day_data)
    if(is.na(bo_cv_r2)) {
      bo_model <- lm(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, data = day_data)
      bo_cv_r2 <- summary(bo_model)$r.squared
      cat("    警告: BO模型交叉验证失败，使用R²值替代\n")
    }
    
    # 步数模型
    steps_cv_r2 <- cv_r2(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, day_data)
    if(is.na(steps_cv_r2)) {
      steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, data = day_data)
      steps_cv_r2 <- summary(steps_model)$r.squared
      cat("    警告: 步数模型交叉验证失败，使用R²值替代\n")
    }
  }, error = function(e) {
    cat(paste0("    错误: 计算交叉验证R²失败: ", e$message, "\n"))
    # 如果交叉验证完全失败，使用普通R²
    full_model <- lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
                       age + gender + bmi + pre_vision + season, data = day_data)
    full_cv_r2 <<- summary(full_model)$r.squared
    
    hr_model <- lm(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, data = day_data)
    hr_cv_r2 <<- summary(hr_model)$r.squared
    
    bo_model <- lm(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, data = day_data)
    bo_cv_r2 <<- summary(bo_model)$r.squared
    
    steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, data = day_data)
    steps_cv_r2 <<- summary(steps_model)$r.squared
    
    cat("    使用常规R²值作为替代\n")
  })
  
  # 将信息添加到非糖尿病组性能数据框
  performance_df_nod <- rbind(performance_df_nod, data.frame(
    time_point = tp,
    r_squared = summary(full_model)$r.squared,
    adj_r_squared = summary(full_model)$adj.r.squared,
    rmse = sqrt(mean(full_model$residuals^2)),
    hr_only_r2 = hr_cv_r2,
    bo_only_r2 = bo_cv_r2,
    steps_only_r2 = steps_cv_r2,
    all_r2 = full_cv_r2,
    sample_size = nrow(day_data),
    stringsAsFactors = FALSE
  ))
  
  cat(paste0("  完成时间点 ", tp, " 的非糖尿病组模型构建。 交叉验证R²: ", round(full_cv_r2, 4), "\n"))
}

# 处理时间点标签，确保正确排序
performance_df_d$time_point_num <- as.numeric(performance_df_d$time_point)
performance_df_d <- performance_df_d[order(performance_df_d$time_point_num), ]

performance_df_nod$time_point_num <- as.numeric(performance_df_nod$time_point)
performance_df_nod <- performance_df_nod[order(performance_df_nod$time_point_num), ]

# 将结果保存到CSV文件
write.csv(performance_df_d, file.path(output_dir, "model_performance_by_day_diabetes.csv"), row.names = FALSE)
write.csv(performance_df_nod, file.path(output_dir, "model_performance_by_day_nodiabetes.csv"), row.names = FALSE)


# 绘制糖尿病组性能变化趋势图
plot_data_d <- data.frame(
  time_point = rep(performance_df_d$time_point, 4),
  time_point_num = rep(performance_df_d$time_point_num, 4),
  predictor = c(
    rep("All Predictors", nrow(performance_df_d)),
    rep("HR", nrow(performance_df_d)),
    rep("BO", nrow(performance_df_d)),
    rep("Total Steps", nrow(performance_df_d))
  ),
  r2 = c(
    performance_df_d$all_r2,
    performance_df_d$hr_only_r2,
    performance_df_d$bo_only_r2,
    performance_df_d$steps_only_r2
  ),
  group = rep("Diabetes", nrow(performance_df_d) * 4)
)

# 绘制非糖尿病组性能变化趋势图
plot_data_nod <- data.frame(
  time_point = rep(performance_df_nod$time_point, 4),
  time_point_num = rep(performance_df_nod$time_point_num, 4),
  predictor = c(
    rep("All Predictors", nrow(performance_df_nod)),
    rep("HR", nrow(performance_df_nod)),
    rep("BO", nrow(performance_df_nod)),
    rep("Total Steps", nrow(performance_df_nod))
  ),
  r2 = c(
    performance_df_nod$all_r2,
    performance_df_nod$hr_only_r2,
    performance_df_nod$bo_only_r2,
    performance_df_nod$steps_only_r2
  ),
  group = rep("No Diabetes", nrow(performance_df_nod) * 4)
)

# 合并两组数据用于绘图
plot_data_combined <- rbind(plot_data_d, plot_data_nod)

# 绘制图形 - 糖尿病组
p_d <- ggplot(plot_data_d, aes(x = time_point_num, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for Vision Improvement - Diabetes Group",
    x = "Days relative to surgery",
    y = "Cross-validated R²",
    color = "Predictors"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# 显示糖尿病组图形
print(p_d)

# 绘制图形 - 非糖尿病组
p_nod <- ggplot(plot_data_nod, aes(x = time_point_num, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for Vision Improvement - No Diabetes Group",
    x = "Days relative to surgery",
    y = "Cross-validated R²",
    color = "Predictors"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# 显示非糖尿病组图形
print(p_nod)

# 绘制组合图形 - 分面
p_combined <- ggplot(plot_data_combined, aes(x = time_point_num, y = r2, color = predictor, group = interaction(predictor, group))) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~ group, ncol = 1) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for Vision Improvement by Diabetes Status",
    x = "Days relative to surgery",
    y = "Cross-validated R²",
    color = "Predictors"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.text = element_text(size = 12, face = "bold")
  )

# 显示组合图形
print(p_combined)

# 保存图形
ggsave(file.path(output_dir, "predictor_performance_by_day_diabetes.pdf"), p_d, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_by_day_diabetes.png"), p_d, width = 10, height = 6, dpi = 300)

ggsave(file.path(output_dir, "predictor_performance_by_day_nodiabetes.pdf"), p_nod, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_by_day_nodiabetes.png"), p_nod, width = 10, height = 6, dpi = 300)

ggsave(file.path(output_dir, "predictor_performance_by_day_combined.pdf"), p_combined, width = 10, height = 10)
ggsave(file.path(output_dir, "predictor_performance_by_day_combined.png"), p_combined, width = 10, height = 10, dpi = 300)

# 创建交互式版本 (如果安装了plotly)
if (requireNamespace("plotly", quietly = TRUE)) {
  p_d_interactive <- plotly::ggplotly(p_d)
  htmlwidgets::saveWidget(p_d_interactive, file.path(output_dir, "predictor_performance_diabetes_interactive.html"))
  
  p_nod_interactive <- plotly::ggplotly(p_nod)
  htmlwidgets::saveWidget(p_nod_interactive, file.path(output_dir, "predictor_performance_nodiabetes_interactive.html"))
  
  p_combined_interactive <- plotly::ggplotly(p_combined)
  htmlwidgets::saveWidget(p_combined_interactive, file.path(output_dir, "predictor_performance_combined_interactive.html"))
  
  cat("交互式图形已保存到:", file.path(output_dir), "目录\n")
}

# 输出糖尿病组模型系数汇总
coef_summary_d <- data.frame(
  time_point = character(),
  predictor = character(),
  estimate = numeric(),
  std_error = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  group = character(),
  stringsAsFactors = FALSE
)

# 输出非糖尿病组模型系数汇总
coef_summary_nod <- data.frame(
  time_point = character(),
  predictor = character(),
  estimate = numeric(),
  std_error = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  group = character(),
  stringsAsFactors = FALSE
)

# 重新运行每个时间点的糖尿病组模型并提取系数
for (tp in time_points_str) {
  file_path <- grep(paste0("day_", tp, "_D\\.csv$"), d_files, value = TRUE)
  
  if (length(file_path) == 0) next
  
  day_data <- tryCatch(read.csv(file_path), error = function(e) NULL)
  if (is.null(day_data)) next
  
  # 检查并移除缺失值
  required_cols <- c("vision_improvement", "mean_rhr_1", "mean_bo", "steps_total", 
                     "age", "gender", "bmi", "pre_vision", "season")
  if (!all(required_cols %in% names(day_data))) next
  
  day_data <- day_data[complete.cases(day_data[, required_cols]), ]
  if (nrow(day_data) < 10) next
  
  # 预处理分类变量
  if("season" %in% names(day_data)) {
    day_data$season <- factor(day_data$season)
  }
  if("gender" %in% names(day_data)) {
    day_data$gender <- factor(day_data$gender)
  }
  
  # 构建模型
  model <- tryCatch({
    lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
         age + gender + bmi + pre_vision + season, data = day_data)
  }, error = function(e) NULL)
  
  if (is.null(model)) next
  
  # 提取系数信息
  model_summary <- summary(model)
  coefs <- model_summary$coefficients
  
  for (i in 1:nrow(coefs)) {
    predictor_name <- rownames(coefs)[i]
    
    # 只关注生理指标
    if (predictor_name %in% c("mean_rhr_1", "mean_bo", "steps_total")) {
      coef_summary_d <- rbind(coef_summary_d, data.frame(
        time_point = tp,
        predictor = predictor_name,
        estimate = coefs[i, 1],
        std_error = coefs[i, 2],
        t_value = coefs[i, 3],
        p_value = coefs[i, 4],
        group = "Diabetes",
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 重新运行每个时间点的非糖尿病组模型并提取系数
for (tp in time_points_str) {
  file_path <- grep(paste0("day_", tp, "_NoD\\.csv$"), nod_files, value = TRUE)
  
  if (length(file_path) == 0) next
  
  day_data <- tryCatch(read.csv(file_path), error = function(e) NULL)
  if (is.null(day_data)) next
  
  # 检查并移除缺失值
  required_cols <- c("vision_improvement", "mean_rhr_1", "mean_bo", "steps_total", 
                     "age", "gender", "bmi", "pre_vision", "season")
  if (!all(required_cols %in% names(day_data))) next
  
  day_data <- day_data[complete.cases(day_data[, required_cols]), ]
  if (nrow(day_data) < 10) next
  
  # 预处理分类变量
  if("season" %in% names(day_data)) {
    day_data$season <- factor(day_data$season)
  }
  if("gender" %in% names(day_data)) {
    day_data$gender <- factor(day_data$gender)
  }
  
  # 构建模型
  model <- tryCatch({
    lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
         age + gender + bmi + pre_vision + season, data = day_data)
  }, error = function(e) NULL)
  
  if (is.null(model)) next
  
  # 提取系数信息
  model_summary <- summary(model)
  coefs <- model_summary$coefficients
  
  for (i in 1:nrow(coefs)) {
    predictor_name <- rownames(coefs)[i]
    
    # 只关注生理指标
    if (predictor_name %in% c("mean_rhr_1", "mean_bo", "steps_total")) {
      coef_summary_nod <- rbind(coef_summary_nod, data.frame(
        time_point = tp,
        predictor = predictor_name,
        estimate = coefs[i, 1],
        std_error = coefs[i, 2],
        t_value = coefs[i, 3],
        p_value = coefs[i, 4],
        group = "No Diabetes",
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 合并两组系数数据
coef_summary_combined <- rbind(coef_summary_d, coef_summary_nod)

# 保存系数汇总
write.csv(coef_summary_d, file.path(output_dir, "coefficient_summary_by_day_diabetes.csv"), row.names = FALSE)
write.csv(coef_summary_nod, file.path(output_dir, "coefficient_summary_by_day_nodiabetes.csv"), row.names = FALSE)
write.csv(coef_summary_combined, file.path(output_dir, "coefficient_summary_by_day_combined.csv"), row.names = FALSE)

cat("糖尿病组系数汇总已保存到:", file.path(output_dir, "coefficient_summary_by_day_diabetes.csv"), "\n")
cat("非糖尿病组系数汇总已保存到:", file.path(output_dir, "coefficient_summary_by_day_nodiabetes.csv"), "\n")
cat("合并系数汇总已保存到:", file.path(output_dir, "coefficient_summary_by_day_combined.csv"), "\n")

# 绘制不同时间点的系数变化图表 - 糖尿病组
if (nrow(coef_summary_d) > 0) {
  # 处理时间点，确保正确排序
  coef_summary_d$time_point_num <- as.numeric(coef_summary_d$time_point)
  coef_summary_d <- coef_summary_d[order(coef_summary_d$time_point_num), ]
  
  # 创建美化版的预测变量名称
  coef_summary_d$predictor_label <- factor(
    coef_summary_d$predictor,
    levels = c("mean_rhr_1", "mean_bo", "steps_total"),
    labels = c("Heart Rate", "Blood Oxygen", "Total Steps")
  )
  
  # 添加显著性标记
  coef_summary_d$significance <- ifelse(coef_summary_d$p_value < 0.05, "Significant", "Not Significant")
  
  # 绘制系数随时间变化的图表 - 糖尿病组
  coef_plot_d <- ggplot(coef_summary_d, 
                        aes(x = time_point_num, y = estimate, 
                            color = predictor_label, group = predictor_label,
                            alpha = significance)) +
    geom_line(size = 1) +
    geom_point(size = 3, aes(shape = significance)) +
    geom_errorbar(aes(ymin = estimate - std_error, ymax = estimate + std_error), 
                  width = 0.2, alpha = 0.7) +
    scale_color_manual(values = c("Heart Rate" = "#82a1bf", 
                                  "Blood Oxygen" = "#faaa93", 
                                  "Total Steps" = "#feefc4")) +
    scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.5)) +
    scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1)) +
    labs(
      title = "Predictor Coefficients for Vision Improvement - Diabetes Group",
      x = "Days relative to surgery",
      y = "Coefficient Estimate",
      color = "Predictor",
      alpha = "Significance",
      shape = "Significance"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # 显示糖尿病组系数图
  print(coef_plot_d)
  
  # 保存糖尿病组系数图
  ggsave(file.path(output_dir, "coefficient_changes_by_day_diabetes.pdf"), coef_plot_d, width = 10, height = 6)
  ggsave(file.path(output_dir, "coefficient_changes_by_day_diabetes.png"), coef_plot_d, width = 10, height = 6, dpi = 300)
}

# 绘制不同时间点的系数变化图表 - 非糖尿病组
if (nrow(coef_summary_nod) > 0) {
  # 处理时间点，确保正确排序
  coef_summary_nod$time_point_num <- as.numeric(coef_summary_nod$time_point)
  coef_summary_nod <- coef_summary_nod[order(coef_summary_nod$time_point_num), ]
  
  # 创建美化版的预测变量名称
  coef_summary_nod$predictor_label <- factor(
    coef_summary_nod$predictor,
    levels = c("mean_rhr_1", "mean_bo", "steps_total"),
    labels = c("Heart Rate", "Blood Oxygen", "Total Steps")
  )
  
  # 添加显著性标记
  coef_summary_nod$significance <- ifelse(coef_summary_nod$p_value < 0.05, "Significant", "Not Significant")
  
  # 绘制系数随时间变化的图表 - 非糖尿病组
  coef_plot_nod <- ggplot(coef_summary_nod, 
                          aes(x = time_point_num, y = estimate, 
                              color = predictor_label, group = predictor_label,
                              alpha = significance)) +
    geom_line(size = 1) +
    geom_point(size = 3, aes(shape = significance)) +
    geom_errorbar(aes(ymin = estimate - std_error, ymax = estimate + std_error), 
                  width = 0.2, alpha = 0.7) +
    scale_color_manual(values = c("Heart Rate" = "#82a1bf", 
                                  "Blood Oxygen" = "#faaa93", 
                                  "Total Steps" = "#feefc4")) +
    scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.5)) +
    scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1)) +
    labs(
      title = "Predictor Coefficients for Vision Improvement - No Diabetes Group",
      x = "Days relative to surgery",
      y = "Coefficient Estimate",
      color = "Predictor",
      alpha = "Significance",
      shape = "Significance"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # 显示非糖尿病组系数图
  print(coef_plot_nod)
  
  # 保存非糖尿病组系数图
  ggsave(file.path(output_dir, "coefficient_changes_by_day_nodiabetes.pdf"), coef_plot_nod, width = 10, height = 6)
  ggsave(file.path(output_dir, "coefficient_changes_by_day_nodiabetes.png"), coef_plot_nod, width = 10, height = 6, dpi = 300)
}

# 绘制合并的系数变化图表 - 分面显示两组
if (nrow(coef_summary_combined) > 0) {
  # 处理时间点，确保正确排序
  coef_summary_combined$time_point_num <- as.numeric(coef_summary_combined$time_point)
  
  # 创建美化版的预测变量名称
  coef_summary_combined$predictor_label <- factor(
    coef_summary_combined$predictor,
    levels = c("mean_rhr_1", "mean_bo", "steps_total"),
    labels = c("Heart Rate", "Blood Oxygen", "Total Steps")
  )
  
  # 添加显著性标记
  coef_summary_combined$significance <- ifelse(coef_summary_combined$p_value < 0.05, "Significant", "Not Significant")
  
  # 绘制系数随时间变化的图表 - 分面显示两组
  coef_plot_combined <- ggplot(coef_summary_combined, 
                               aes(x = time_point_num, y = estimate, 
                                   color = predictor_label, group = predictor_label,
                                   alpha = significance)) +
    geom_line(size = 1) +
    geom_point(size = 3, aes(shape = significance)) +
    geom_errorbar(aes(ymin = estimate - std_error, ymax = estimate + std_error), 
                  width = 0.2, alpha = 0.7) +
    facet_wrap(~ group, ncol = 1) +
    scale_color_manual(values = c("Heart Rate" = "#82a1bf", 
                                  "Blood Oxygen" = "#faaa93", 
                                  "Total Steps" = "#feefc4")) +
    scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.5)) +
    scale_shape_manual(values = c("Significant" = 16, "Not Significant" = 1)) +
    labs(
      title = "Predictor Coefficients for Vision Improvement by Diabetes Status",
      x = "Days relative to surgery",
      y = "Coefficient Estimate",
      color = "Predictor",
      alpha = "Significance",
      shape = "Significance"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  # 显示合并系数图
  print(coef_plot_combined)
  
  # 保存合并系数图
  ggsave(file.path(output_dir, "coefficient_changes_by_day_combined.pdf"), coef_plot_combined, width = 10, height = 10)
  ggsave(file.path(output_dir, "coefficient_changes_by_day_combined.png"), coef_plot_combined, width = 10, height = 10, dpi = 300)
}

# 创建每个预测变量的单独比较图
# Heart Rate比较
if(sum(coef_summary_d$predictor == "mean_rhr_1") > 0 && sum(coef_summary_nod$predictor == "mean_rhr_1") > 0) {
  hr_data <- coef_summary_combined[coef_summary_combined$predictor == "mean_rhr_1", ]
  
  p_hr_compare <- ggplot(hr_data, aes(x = time_point_num, y = estimate, color = group, group = group)) +
    geom_line(size = 1.2) +
    geom_point(size = 3, aes(shape = significance)) +
    geom_errorbar(aes(ymin = estimate - std_error, ymax = estimate + std_error), 
                  width = 0.2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Heart Rate Coefficient Comparison by Diabetes Status",
      x = "Days relative to surgery",
      y = "Coefficient Estimate",
      color = "Group",
      shape = "Significance"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  print(p_hr_compare)
  ggsave(file.path(output_dir, "hr_coefficient_comparison.pdf"), p_hr_compare, width = 10, height = 6)
  ggsave(file.path(output_dir, "hr_coefficient_comparison.png"), p_hr_compare, width = 10, height = 6, dpi = 300)
}

# Blood Oxygen比较
if(sum(coef_summary_d$predictor == "mean_bo") > 0 && sum(coef_summary_nod$predictor == "mean_bo") > 0) {
  bo_data <- coef_summary_combined[coef_summary_combined$predictor == "mean_bo", ]
  
  p_bo_compare <- ggplot(bo_data, aes(x = time_point_num, y = estimate, color = group, group = group)) +
    geom_line(size = 1.2) +
    geom_point(size = 3, aes(shape = significance)) +
    geom_errorbar(aes(ymin = estimate - std_error, ymax = estimate + std_error), 
                  width = 0.2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Blood Oxygen Coefficient Comparison by Diabetes Status",
      x = "Days relative to surgery",
      y = "Coefficient Estimate",
      color = "Group",
      shape = "Significance"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  print(p_bo_compare)
  ggsave(file.path(output_dir, "bo_coefficient_comparison.pdf"), p_bo_compare, width = 10, height = 6)
  ggsave(file.path(output_dir, "bo_coefficient_comparison.png"), p_bo_compare, width = 10, height = 6, dpi = 300)
}

# Steps比较
if(sum(coef_summary_d$predictor == "steps_total") > 0 && sum(coef_summary_nod$predictor == "steps_total") > 0) {
  steps_data <- coef_summary_combined[coef_summary_combined$predictor == "steps_total", ]
  
  p_steps_compare <- ggplot(steps_data, aes(x = time_point_num, y = estimate, color = group, group = group)) +
    geom_line(size = 1.2) +
    geom_point(size = 3, aes(shape = significance)) +
    geom_errorbar(aes(ymin = estimate - std_error, ymax = estimate + std_error), 
                  width = 0.2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Steps Coefficient Comparison by Diabetes Status",
      x = "Days relative to surgery",
      y = "Coefficient Estimate",
      color = "Group",
      shape = "Significance"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  print(p_steps_compare)
  ggsave(file.path(output_dir, "steps_coefficient_comparison.pdf"), p_steps_compare, width = 10, height = 6)
  ggsave(file.path(output_dir, "steps_coefficient_comparison.png"), p_steps_compare, width = 10, height = 6, dpi = 300)
}

# 总结分析结果
cat("\n=== 糖尿病组与非糖尿病组每日预测模型性能分析摘要 ===\n")

# 找出各组中最佳时间点
best_day_d <- performance_df_d[which.max(performance_df_d$all_r2), ]
best_day_nod <- performance_df_nod[which.max(performance_df_nod$all_r2), ]

# 找出各组中每个预测变量的最佳时间点
best_hr_d <- performance_df_d[which.max(performance_df_d$hr_only_r2), ]
best_bo_d <- performance_df_d[which.max(performance_df_d$bo_only_r2), ]
best_steps_d <- performance_df_d[which.max(performance_df_d$steps_only_r2), ]

best_hr_nod <- performance_df_nod[which.max(performance_df_nod$hr_only_r2), ]
best_bo_nod <- performance_df_nod[which.max(performance_df_nod$bo_only_r2), ]
best_steps_nod <- performance_df_nod[which.max(performance_df_nod$steps_only_r2), ]

# 输出糖尿病组摘要
cat("\n糖尿病组:\n")
cat("- 预测性能最好的天数: Day ", best_day_d$time_point, 
    " (R² = ", round(best_day_d$all_r2, 4), ")\n", sep = "")
cat("- HR最佳预测天数: Day ", best_hr_d$time_point, 
    " (R² = ", round(best_hr_d$hr_only_r2, 4), ")\n", sep = "")
cat("- BO最佳预测天数: Day ", best_bo_d$time_point, 
    " (R² = ", round(best_bo_d$bo_only_r2, 4), ")\n", sep = "")
cat("- Steps最佳预测天数: Day ", best_steps_d$time_point, 
    " (R² = ", round(best_steps_d$steps_only_r2, 4), ")\n", sep = "")

# 检查显著性
sig_hr_d <- coef_summary_d[coef_summary_d$predictor == "mean_rhr_1" & coef_summary_d$p_value < 0.05, ]
sig_bo_d <- coef_summary_d[coef_summary_d$predictor == "mean_bo" & coef_summary_d$p_value < 0.05, ]
sig_steps_d <- coef_summary_d[coef_summary_d$predictor == "steps_total" & coef_summary_d$p_value < 0.05, ]

if(nrow(sig_hr_d) > 0) {
  cat("- 心率在以下天数显著: ", paste(sig_hr_d$time_point, collapse=", "), "\n")
}
if(nrow(sig_bo_d) > 0) {
  cat("- 血氧在以下天数显著: ", paste(sig_bo_d$time_point, collapse=", "), "\n")
}
if(nrow(sig_steps_d) > 0) {
  cat("- 步数在以下天数显著: ", paste(sig_steps_d$time_point, collapse=", "), "\n")
}

# 输出非糖尿病组摘要
cat("\n非糖尿病组:\n")
cat("- 预测性能最好的天数: Day ", best_day_nod$time_point, 
    " (R² = ", round(best_day_nod$all_r2, 4), ")\n", sep = "")
cat("- HR最佳预测天数: Day ", best_hr_nod$time_point, 
    " (R² = ", round(best_hr_nod$hr_only_r2, 4), ")\n", sep = "")
cat("- BO最佳预测天数: Day ", best_bo_nod$time_point, 
    " (R² = ", round(best_bo_nod$bo_only_r2, 4), ")\n", sep = "")
cat("- Steps最佳预测天数: Day ", best_steps_nod$time_point, 
    " (R² = ", round(best_steps_nod$steps_only_r2, 4), ")\n", sep = "")

# 检查显著性
sig_hr_nod <- coef_summary_nod[coef_summary_nod$predictor == "mean_rhr_1" & coef_summary_nod$p_value < 0.05, ]
sig_bo_nod <- coef_summary_nod[coef_summary_nod$predictor == "mean_bo" & coef_summary_nod$p_value < 0.05, ]
sig_steps_nod <- coef_summary_nod[coef_summary_nod$predictor == "steps_total" & coef_summary_nod$p_value < 0.05, ]

if(nrow(sig_hr_nod) > 0) {
  cat("- 心率在以下天数显著: ", paste(sig_hr_nod$time_point, collapse=", "), "\n")
}
if(nrow(sig_bo_nod) > 0) {
  cat("- 血氧在以下天数显著: ", paste(sig_bo_nod$time_point, collapse=", "), "\n")
}
if(nrow(sig_steps_nod) > 0) {
  cat("- 步数在以下天数显著: ", paste(sig_steps_nod$time_point, collapse=", "), "\n")
}

# 比较两组
cat("\n两组比较:\n")
if(best_day_d$all_r2 > best_day_nod$all_r2) {
  cat("- 糖尿病组的预测模型整体表现优于非糖尿病组\n")
} else {
  cat("- 非糖尿病组的预测模型整体表现优于糖尿病组\n")
}

# 计算平均性能
avg_hr_d <- mean(performance_df_d$hr_only_r2, na.rm = TRUE)
avg_bo_d <- mean(performance_df_d$bo_only_r2, na.rm = TRUE)
avg_steps_d <- mean(performance_df_d$steps_only_r2, na.rm = TRUE)

avg_hr_nod <- mean(performance_df_nod$hr_only_r2, na.rm = TRUE)
avg_bo_nod <- mean(performance_df_nod$bo_only_r2, na.rm = TRUE)
avg_steps_nod <- mean(performance_df_nod$steps_only_r2, na.rm = TRUE)

cat("- 糖尿病组平均预测性能: HR (R² = ", round(avg_hr_d, 4), 
    "), BO (R² = ", round(avg_bo_d, 4), 
    "), Steps (R² = ", round(avg_steps_d, 4), ")\n", sep = "")

cat("- 非糖尿病组平均预测性能: HR (R² = ", round(avg_hr_nod, 4), 
    "), BO (R² = ", round(avg_bo_nod, 4), 
    "), Steps (R² = ", round(avg_steps_nod, 4), ")\n", sep = "")

# 计算术前和术后天数的平均预测力
pre_days_d <- performance_df_d[as.numeric(performance_df_d$time_point) < 0, ]
post_days_d <- performance_df_d[as.numeric(performance_df_d$time_point) >= 0, ]

pre_days_nod <- performance_df_nod[as.numeric(performance_df_nod$time_point) < 0, ]
post_days_nod <- performance_df_nod[as.numeric(performance_df_nod$time_point) >= 0, ]

if(nrow(pre_days_d) > 0 && nrow(post_days_d) > 0) {
  pre_avg_d <- mean(pre_days_d$all_r2, na.rm = TRUE)
  post_avg_d <- mean(post_days_d$all_r2, na.rm = TRUE)
  
  cat("- 糖尿病组: 术前天数平均R² = ", round(pre_avg_d, 4), 
      ", 术后天数平均R² = ", round(post_avg_d, 4), 
      ", 差异 = ", round(post_avg_d - pre_avg_d, 4), 
      " (", ifelse(post_avg_d > pre_avg_d, "术后更好", "术前更好"), ")\n", sep = "")
}

if(nrow(pre_days_nod) > 0 && nrow(post_days_nod) > 0) {
  pre_avg_nod <- mean(pre_days_nod$all_r2, na.rm = TRUE)
  post_avg_nod <- mean(post_days_nod$all_r2, na.rm = TRUE)
  
  cat("- 非糖尿病组: 术前天数平均R² = ", round(pre_avg_nod, 4), 
      ", 术后天数平均R² = ", round(post_avg_nod, 4), 
      ", 差异 = ", round(post_avg_nod - pre_avg_nod, 4), 
      " (", ifelse(post_avg_nod > pre_avg_nod, "术后更好", "术前更好"), ")\n", sep = "")
}

# 创建详细的对比表格
comparison_table <- data.frame(
  Metric = character(),
  Diabetes = numeric(),
  No_Diabetes = numeric(),
  Difference = numeric(),
  Better_Group = character(),
  stringsAsFactors = FALSE
)

# 添加整体性能指标
comparison_table <- rbind(comparison_table, data.frame(
  Metric = "平均整体预测能力 (全部天数)",
  Diabetes = mean(performance_df_d$all_r2, na.rm = TRUE),
  No_Diabetes = mean(performance_df_nod$all_r2, na.rm = TRUE),
  Difference = abs(mean(performance_df_d$all_r2, na.rm = TRUE) - mean(performance_df_nod$all_r2, na.rm = TRUE)),
  Better_Group = ifelse(mean(performance_df_d$all_r2, na.rm = TRUE) > mean(performance_df_nod$all_r2, na.rm = TRUE), 
                        "糖尿病组", "非糖尿病组"),
  stringsAsFactors = FALSE
))

# 添加术前和术后性能
if(nrow(pre_days_d) > 0 && nrow(pre_days_nod) > 0) {
  comparison_table <- rbind(comparison_table, data.frame(
    Metric = "术前天数平均预测能力",
    Diabetes = pre_avg_d,
    No_Diabetes = pre_avg_nod,
    Difference = abs(pre_avg_d - pre_avg_nod),
    Better_Group = ifelse(pre_avg_d > pre_avg_nod, "糖尿病组", "非糖尿病组"),
    stringsAsFactors = FALSE
  ))
}

if(nrow(post_days_d) > 0 && nrow(post_days_nod) > 0) {
  comparison_table <- rbind(comparison_table, data.frame(
    Metric = "术后天数平均预测能力",
    Diabetes = post_avg_d,
    No_Diabetes = post_avg_nod,
    Difference = abs(post_avg_d - post_avg_nod),
    Better_Group = ifelse(post_avg_d > post_avg_nod, "糖尿病组", "非糖尿病组"),
    stringsAsFactors = FALSE
  ))
}

# 添加各预测变量的性能
comparison_table <- rbind(comparison_table, data.frame(
  Metric = "HR平均预测能力",
  Diabetes = avg_hr_d,
  No_Diabetes = avg_hr_nod,
  Difference = abs(avg_hr_d - avg_hr_nod),
  Better_Group = ifelse(avg_hr_d > avg_hr_nod, "糖尿病组", "非糖尿病组"),
  stringsAsFactors = FALSE
))

comparison_table <- rbind(comparison_table, data.frame(
  Metric = "BO平均预测能力",
  Diabetes = avg_bo_d,
  No_Diabetes = avg_bo_nod,
  Difference = abs(avg_bo_d - avg_bo_nod),
  Better_Group = ifelse(avg_bo_d > avg_bo_nod, "糖尿病组", "非糖尿病组"),
  stringsAsFactors = FALSE
))

comparison_table <- rbind(comparison_table, data.frame(
  Metric = "Steps平均预测能力",
  Diabetes = avg_steps_d,
  No_Diabetes = avg_steps_nod,
  Difference = abs(avg_steps_d - avg_steps_nod),
  Better_Group = ifelse(avg_steps_d > avg_steps_nod, "糖尿病组", "非糖尿病组"),
  stringsAsFactors = FALSE
))

# 添加最佳天数性能
comparison_table <- rbind(comparison_table, data.frame(
  Metric = "最佳天数预测能力",
  Diabetes = max(performance_df_d$all_r2, na.rm = TRUE),
  No_Diabetes = max(performance_df_nod$all_r2, na.rm = TRUE),
  Difference = abs(max(performance_df_d$all_r2, na.rm = TRUE) - max(performance_df_nod$all_r2, na.rm = TRUE)),
  Better_Group = ifelse(max(performance_df_d$all_r2, na.rm = TRUE) > max(performance_df_nod$all_r2, na.rm = TRUE), 
                        "糖尿病组", "非糖尿病组"),
  stringsAsFactors = FALSE
))

# 输出对比表格
print(comparison_table)

# 保存对比表格
write.csv(comparison_table, file.path(output_dir, "diabetes_vs_nodiabetes_comparison.csv"), row.names = FALSE)
cat("糖尿病与非糖尿病组对比表已保存到:", file.path(output_dir, "diabetes_vs_nodiabetes_comparison.csv"), "\n")

# 找出每组最佳预测变量
best_var_d <- which.max(c(avg_hr_d, avg_bo_d, avg_steps_d))
best_var_names_d <- c("心率", "血氧", "步数")
best_var_d_name <- best_var_names_d[best_var_d]

best_var_nod <- which.max(c(avg_hr_nod, avg_bo_nod, avg_steps_nod))
best_var_names_nod <- c("心率", "血氧", "步数")
best_var_nod_name <- best_var_names_nod[best_var_nod]

# 时间趋势分析 - 糖尿病组
result_trend_d <- data.frame(
  time_point = performance_df_d$time_point_num,
  r2 = performance_df_d$all_r2
)

# 检查是否有足够的数据点进行线性回归
if(nrow(result_trend_d) >= 5) {
  trend_model_d <- lm(r2 ~ time_point, data = result_trend_d)
  trend_p_d <- summary(trend_model_d)$coefficients[2, 4]
  trend_dir_d <- ifelse(coef(trend_model_d)[2] > 0, "升高", "降低")
  
  cat("\n糖尿病组随时间的预测力变化:\n")
  cat("- 线性趋势: 随着时间的推移，预测力", trend_dir_d, 
      " (系数 = ", round(coef(trend_model_d)[2], 4), 
      ", p值 = ", round(trend_p_d, 4), ")\n", sep = "")
  
  if(trend_p_d < 0.05) {
    cat("- 该趋势统计学显著\n")
  } else {
    cat("- 该趋势统计学不显著\n")
  }
}

# 时间趋势分析 - 非糖尿病组
result_trend_nod <- data.frame(
  time_point = performance_df_nod$time_point_num,
  r2 = performance_df_nod$all_r2
)

# 检查是否有足够的数据点进行线性回归
if(nrow(result_trend_nod) >= 5) {
  trend_model_nod <- lm(r2 ~ time_point, data = result_trend_nod)
  trend_p_nod <- summary(trend_model_nod)$coefficients[2, 4]
  trend_dir_nod <- ifelse(coef(trend_model_nod)[2] > 0, "升高", "降低")
  
  cat("\n非糖尿病组随时间的预测力变化:\n")
  cat("- 线性趋势: 随着时间的推移，预测力", trend_dir_nod, 
      " (系数 = ", round(coef(trend_model_nod)[2], 4), 
      ", p值 = ", round(trend_p_nod, 4), ")\n", sep = "")
  
  if(trend_p_nod < 0.05) {
    cat("- 该趋势统计学显著\n")
  } else {
    cat("- 该趋势统计学不显著\n")
  }
}

# 总体结论
cat("\n总体结论:\n")
cat("- 糖尿病患者组最佳预测变量是", best_var_d_name, "\n", sep = "")
cat("- 非糖尿病患者组最佳预测变量是", best_var_nod_name, "\n", sep = "")

if(!is.null(pre_avg_d) && !is.null(post_avg_d) && !is.null(pre_avg_nod) && !is.null(post_avg_nod)) {
  if(post_avg_d > pre_avg_d && post_avg_nod > pre_avg_nod) {
    cat("- 对两组患者而言，术后数据的预测能力都优于术前数据\n")
  } else if(pre_avg_d > post_avg_d && pre_avg_nod > post_avg_nod) {
    cat("- 对两组患者而言，术前数据的预测能力都优于术后数据\n")
  } else if(pre_avg_d > post_avg_d && post_avg_nod > pre_avg_nod) {
    cat("- 对糖尿病患者，术前数据预测能力更强；而对非糖尿病患者，术后数据预测能力更强\n")
  } else if(post_avg_d > pre_avg_d && pre_avg_nod > post_avg_nod) {
    cat("- 对糖尿病患者，术后数据预测能力更强；而对非糖尿病患者，术前数据预测能力更强\n")
  }
}

# 预测变量的组间差异
if(best_var_d_name == best_var_nod_name) {
  cat("- 两组患者的最佳预测变量相同，均为", best_var_d_name, "\n", sep = "")
} else {
  cat("- 两组患者的最佳预测变量不同：糖尿病组为", best_var_d_name, 
      "，非糖尿病组为", best_var_nod_name, "\n", sep = "")
}

# 创建额外的热图可视化
if(requireNamespace("reshape2", quietly = TRUE)) {
  # 糖尿病组热图数据
  heatmap_data_d <- reshape2::melt(performance_df_d, 
                                   id.vars = c("time_point", "time_point_num"), 
                                   measure.vars = c("hr_only_r2", "bo_only_r2", "steps_only_r2"),
                                   variable.name = "predictor", 
                                   value.name = "r2")
  
  heatmap_data_d$predictor <- factor(heatmap_data_d$predictor,
                                     levels = c("hr_only_r2", "bo_only_r2", "steps_only_r2"),
                                     labels = c("Heart Rate", "Blood Oxygen", "Total Steps"))
  
  # 非糖尿病组热图数据
  heatmap_data_nod <- reshape2::melt(performance_df_nod, 
                                     id.vars = c("time_point", "time_point_num"), 
                                     measure.vars = c("hr_only_r2", "bo_only_r2", "steps_only_r2"),
                                     variable.name = "predictor", 
                                     value.name = "r2")
  
  heatmap_data_nod$predictor <- factor(heatmap_data_nod$predictor,
                                       levels = c("hr_only_r2", "bo_only_r2", "steps_only_r2"),
                                       labels = c("Heart Rate", "Blood Oxygen", "Total Steps"))
  
  # 糖尿病组热图
  p_heatmap_d <- ggplot(heatmap_data_d, aes(x = time_point_num, y = predictor, fill = r2)) +
    geom_tile() +
    scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                         midpoint = median(heatmap_data_d$r2, na.rm = TRUE)) +
    labs(
      title = "Predictive Power Across Days - Diabetes Group",
      x = "Days relative to surgery",
      y = "Predictor",
      fill = "R²"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # 非糖尿病组热图
  p_heatmap_nod <- ggplot(heatmap_data_nod, aes(x = time_point_num, y = predictor, fill = r2)) +
    geom_tile() +
    scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                         midpoint = median(heatmap_data_nod$r2, na.rm = TRUE)) +
    labs(
      title = "Predictive Power Across Days - No Diabetes Group",
      x = "Days relative to surgery",
      y = "Predictor",
      fill = "R²"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # 保存热图
  ggsave(file.path(output_dir, "heatmap_predictive_power_diabetes.pdf"), p_heatmap_d, width = 10, height = 6)
  ggsave(file.path(output_dir, "heatmap_predictive_power_diabetes.png"), p_heatmap_d, width = 10, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, "heatmap_predictive_power_nodiabetes.pdf"), p_heatmap_nod, width = 10, height = 6)
  ggsave(file.path(output_dir, "heatmap_predictive_power_nodiabetes.png"), p_heatmap_nod, width = 10, height = 6, dpi = 300)
  
  cat("预测力热图已保存到输出目录\n")
}

