library(tidyverse)
library(lme4)
library(caret)
library(ggplot2)
setwd(get_project_wd())
rm(list = ls())

# 设置工作路径（请根据需要调整）
input_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped"
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped/lmm_model_performance"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# 获取所有分组的日期文件列表
d_surg1_files <- list.files(input_dir, pattern = "D_Surg1\\.csv$", full.names = TRUE)
nod_surg0_files <- list.files(input_dir, pattern = "NoD_Surg0\\.csv$", full.names = TRUE)

# 将文件名转换为时间点列表
get_time_point <- function(file_path) {
  time_str <- gsub(".*day_(.*)_[DN].*\\.csv$", "\\1", basename(file_path))
  return(time_str)
}

time_points_d <- sapply(d_surg1_files, get_time_point)
time_points_nod <- sapply(nod_surg0_files, get_time_point)

# 获取两组共有的时间点
time_points <- intersect(time_points_d, time_points_nod)

# 对时间点进行排序（从-7到6）
time_points <- sort(as.numeric(time_points))
time_points_str <- as.character(time_points)

# 创建数据框来存储每个时间点模型的性能指标
performance_df <- data.frame(
  time_point = character(),
  r_squared = numeric(),
  adj_r_squared = numeric(),
  rmse = numeric(),
  hr_only_r2 = numeric(),
  bo_only_r2 = numeric(),
  steps_only_r2 = numeric(),
  stringsAsFactors = FALSE
)

# 存储每个预测变量单独的性能
hr_performance <- numeric(length(time_points))
bo_performance <- numeric(length(time_points))
steps_performance <- numeric(length(time_points))
all_performance <- numeric(length(time_points))

# 对每个时间点构建模型
cat("开始为每个时间点构建预测模型...\n")

for (i in 1:length(time_points)) {
  tp <- time_points_str[i]
  cat(paste0("处理时间点: ", tp, "\n"))
  
  # 读取该时间点的数据
  file_path <- grep(paste0("day_", tp, "_D_Surg1\\.csv$"), d_surg1_files, value = TRUE)
  
  if (length(file_path) == 0) {
    cat(paste0("  警告: 没有找到时间点 ", tp, " 的数据文件\n"))
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
    cat(paste0("  警告: 时间点 ", tp, " 的有效数据少于10行，可能导致模型不可靠\n"))
    next
  }
  
  # 构建完整模型
  full_model <- tryCatch({
    lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
         age + gender + bmi + pre_vision + season, data = day_data)
  }, error = function(e) {
    cat(paste0("  错误: 构建完整模型失败: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(full_model)) next
  
  # 构建单一预测变量模型
  hr_model <- lm(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, 
                 data = day_data)
  bo_model <- lm(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, 
                 data = day_data)
  steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, 
                    data = day_data)
  
  # 预处理分类变量
  # 确保所有分类变量都被转换为因子，并且在所有分析中使用相同的水平
  if("season" %in% names(day_data)) {
    day_data$season <- factor(day_data$season)
    cat(paste0("    季节水平: ", paste(levels(day_data$season), collapse=", "), "\n"))
  }
  if("gender" %in% names(day_data)) {
    day_data$gender <- factor(day_data$gender)
  }
  
  # 设置随机种子以确保可重复性
  set.seed(123)
  
  # 确保所有分类变量在交叉验证前都被正确处理
  if("season" %in% names(day_data)) {
    day_data$season <- factor(day_data$season)  # 确保是因子
  }
  if("gender" %in% names(day_data)) {
    day_data$gender <- factor(day_data$gender)  # 确保是因子
  }
  
  # 函数计算交叉验证的R²
  cv_r2 <- function(model_formula) {
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
        data = day_data,
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
  
  # 计算各模型的交叉验证R²
  cat("    计算交叉验证R²值...\n")
  
  # 处理可能的错误并提供回退方案
  tryCatch({
    # 完整模型
    full_cv_r2 <- cv_r2(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + 
                          age + gender + bmi + pre_vision + season)
    if(is.na(full_cv_r2)) {
      # 如果交叉验证失败，使用简单的R²作为后备
      full_cv_r2 <- summary(full_model)$r.squared
      cat("    警告: 完整模型交叉验证失败，使用R²值替代\n")
    }
    
    # HR模型
    hr_cv_r2 <- cv_r2(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season)
    if(is.na(hr_cv_r2)) {
      hr_model <- lm(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, data = day_data)
      hr_cv_r2 <- summary(hr_model)$r.squared
      cat("    警告: HR模型交叉验证失败，使用R²值替代\n")
    }
    
    # BO模型
    bo_cv_r2 <- cv_r2(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season)
    if(is.na(bo_cv_r2)) {
      bo_model <- lm(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, data = day_data)
      bo_cv_r2 <- summary(bo_model)$r.squared
      cat("    警告: BO模型交叉验证失败，使用R²值替代\n")
    }
    
    # 步数模型
    steps_cv_r2 <- cv_r2(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season)
    if(is.na(steps_cv_r2)) {
      steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, data = day_data)
      steps_cv_r2 <- summary(steps_model)$r.squared
      cat("    警告: 步数模型交叉验证失败，使用R²值替代\n")
    }
  }, error = function(e) {
    cat(paste0("    错误: 计算交叉验证R²失败: ", e$message, "\n"))
    # 如果交叉验证完全失败，使用普通R²
    full_cv_r2 <<- summary(full_model)$r.squared
    
    hr_model <- lm(vision_improvement ~ mean_rhr_1 + age + gender + bmi + pre_vision + season, data = day_data)
    hr_cv_r2 <<- summary(hr_model)$r.squared
    
    bo_model <- lm(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season, data = day_data)
    bo_cv_r2 <<- summary(bo_model)$r.squared
    
    steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, data = day_data)
    steps_cv_r2 <<- summary(steps_model)$r.squared
    
    cat("    使用常规R²值作为替代\n")
  })
  
  # 存储性能指标
  hr_performance[i] <- hr_cv_r2
  bo_performance[i] <- bo_cv_r2
  steps_performance[i] <- steps_cv_r2
  all_performance[i] <- full_cv_r2
  
  # 将信息添加到性能数据框
  performance_df <- rbind(performance_df, data.frame(
    time_point = tp,
    r_squared = summary(full_model)$r.squared,
    adj_r_squared = summary(full_model)$adj.r.squared,
    rmse = sqrt(mean(full_model$residuals^2)),
    hr_only_r2 = hr_cv_r2,
    bo_only_r2 = bo_cv_r2,
    steps_only_r2 = steps_cv_r2,
    all_r2 = full_cv_r2,
    stringsAsFactors = FALSE
  ))
  
  cat(paste0("  完成时间点 ", tp, " 的模型构建。 交叉验证R²: ", round(full_cv_r2, 4), "\n"))
}

# 处理时间点标签
performance_df$time_point_num <- as.numeric(performance_df$time_point)
performance_df <- performance_df[order(performance_df$time_point_num), ]

# 将结果保存到CSV文件
write.csv(performance_df, file.path(output_dir, "model_performance_by_day.csv"), row.names = FALSE)
cat(paste0("结果已保存到: ", file.path(output_dir, "model_performance_by_day.csv"), "\n"))

# 绘制性能变化趋势图
# 将数据转换为长格式，便于ggplot绘图
plot_data <- data.frame(
  time_point = rep(performance_df$time_point, 4),
  time_point_num = rep(performance_df$time_point_num, 4),
  predictor = c(
    rep("All Predictors", nrow(performance_df)),
    rep("HR", nrow(performance_df)),
    rep("BO", nrow(performance_df)),
    rep("Total Steps", nrow(performance_df))
  ),
  r2 = c(
    performance_df$all_r2,
    performance_df$hr_only_r2,
    performance_df$bo_only_r2,
    performance_df$steps_only_r2
  )
)

# 绘制图形
p <- ggplot(plot_data, aes(x = time_point_num, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for Vision Improvement",
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
p
# 保存图形
ggsave(file.path(output_dir, "predictor_performance_by_day.pdf"), p, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_by_day.png"), p, width = 10, height = 6, dpi = 300)


# 创建交互式版本 (如果安装了plotly)
if (requireNamespace("plotly", quietly = TRUE)) {
  p_interactive <- plotly::ggplotly(p)
  htmlwidgets::saveWidget(p_interactive, file.path(output_dir, "predictor_performance_interactive.html"))
  cat("交互式图形已保存到:", file.path(output_dir, "predictor_performance_interactive.html"), "\n")
}

# 输出模型系数汇总
coef_summary <- data.frame(
  time_point = character(),
  predictor = character(),
  estimate = numeric(),
  std_error = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# 重新运行每个时间点的模型并提取系数
for (tp in time_points_str) {
  file_path <- grep(paste0("day_", tp, "_D_Surg1\\.csv$"), d_surg1_files, value = TRUE)
  
  if (length(file_path) == 0) next
  
  day_data <- tryCatch(read.csv(file_path), error = function(e) NULL)
  if (is.null(day_data)) next
  
  # 检查并移除缺失值
  required_cols <- c("vision_improvement", "mean_rhr_1", "mean_bo", "steps_total", 
                     "age", "gender", "bmi", "pre_vision", "season")
  if (!all(required_cols %in% names(day_data))) next
  
  day_data <- day_data[complete.cases(day_data[, required_cols]), ]
  if (nrow(day_data) < 10) next
  
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
      coef_summary <- rbind(coef_summary, data.frame(
        time_point = tp,
        predictor = predictor_name,
        estimate = coefs[i, 1],
        std_error = coefs[i, 2],
        t_value = coefs[i, 3],
        p_value = coefs[i, 4],
        stringsAsFactors = FALSE
      ))
    }
  }
}

# 保存系数汇总
write.csv(coef_summary, file.path(output_dir, "coefficient_summary_by_day.csv"), row.names = FALSE)
cat("系数汇总已保存到:", file.path(output_dir, "coefficient_summary_by_day.csv"), "\n")
