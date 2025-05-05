library(tidyverse)
library(lme4)
library(caret)
library(ggplot2)
setwd(get_project_wd())
rm(list = ls())

# 设置工作路径
input_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data_grouped/dm_group"
output_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data_grouped/dm_group/lmm_model_performance/"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义标准时间窗口
time_windows <- c(
  "pre_3d",
  "pre_3d_7d", 
  "pre_7d_all", 
  "post_7d",
  "post_7d_30d",
  "post_day23_30",
  "post_day27_30"
  
)

# 定义时间窗口映射
time_window_mapping <- list(
  "pre_3d" = "pre_3d",
  "pre_3d_7d" = "pre_3d_7d", 
  "pre_7d_all" = "pre_7d_all",
  "post_day23_30" = "post_day23_30",
  "post_day27_30" = "post_day27_30",
  "post_7d_30d" = "post_7d_30d"
)

# 获取所有CSV文件
all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# 将文件分类为糖尿病组和非糖尿病组
d_files <- grep("_D\\.csv$", all_files, value = TRUE)
nod_files <- grep("_NoD\\.csv$", all_files, value = TRUE)

# 改进后的get_time_window函数
get_time_window <- function(file_path) {
  file_name <- basename(file_path)
  
  # 检查特定模式，从最具体到最一般
  if (grepl("^pre_3d_7d_[DN]", file_name)) {
    return("pre_3d_7d")
  } else if (grepl("^pre_7d_all_[DN]", file_name)) {
    return("pre_7d_all")
  } else if (grepl("^pre_3d_[DN]", file_name)) { 
    return("pre_3d")
  } else if (grepl("^post_7d_[DN]", file_name)) {
    return("post_7d")
  } else if (grepl("^post_7d_30d_[DN]", file_name)) {
    return("post_7d_30d")
  } else if (grepl("^post_day23_30_[DN]", file_name)) {
    return("post_day23_30")
  } else if (grepl("^post_day27_30_[DN]", file_name)) {
    return("post_day27_30")
  }
  
  # 如果没有匹配到预定义窗口，返回NA
  cat("警告: 无法识别文件", file_name, "的时间窗口\n")
  return(NA)
}

# 创建数据框来存储每个时间窗口模型的性能指标
performance_df <- data.frame(
  time_window = character(),
  r_squared = numeric(),
  adj_r_squared = numeric(),
  rmse = numeric(),
  hr_only_r2 = numeric(),
  bo_only_r2 = numeric(),
  steps_only_r2 = numeric(),
  improve_1w_only_r2 = numeric(),
  all_r2 = numeric(),
  sample_size = numeric(),
  stringsAsFactors = FALSE
)

# 存储每个预测变量单独的性能
hr_performance <- numeric(length(time_windows))
bo_performance <- numeric(length(time_windows))
steps_performance <- numeric(length(time_windows))
improve_1w_performance <- numeric(length(time_windows))
all_performance <- numeric(length(time_windows))

# 基于文件名跟踪我们已经处理过的窗口
processed_windows <- character(0)

# 对每个时间窗口构建模型
cat("开始为每个时间窗口构建预测模型...\n")

# 首先打印一下所有找到的文件，便于调试
cat("找到的糖尿病组(D)文件:\n")
for (file in d_files) {
  cat("  ", basename(file), "\n")
}

for (i in 1:length(d_files)) {
  file_path <- d_files[i]
  file_name <- basename(file_path)
  current_window <- get_time_window(file_path)
  
  cat("处理文件:", file_name, "\n")
  cat("  识别的时间窗口:", if(is.na(current_window)) "未识别" else current_window, "\n")
  
  if (is.na(current_window) || !current_window %in% time_windows) {
    cat("  跳过: 不在指定的时间窗口列表中\n")
    next
  }
  
  # 检查是否已处理过该窗口
  if (current_window %in% processed_windows) {
    cat("  跳过: 已经处理过", current_window, "窗口\n")
    next
  }
  
  # 记录已处理的窗口
  processed_windows <- c(processed_windows, current_window)
  
  cat(paste0("处理时间窗口: ", current_window, "\n"))
  
  # 读取数据
  period_data <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    cat(paste0("  错误: 无法读取文件: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(period_data)) next
  
  # 检查数据是否包含需要的列
  required_cols <- c("vision_improvement_1m", "min_rhr_steps_1", "bo_mean", "steps_total", 
                     "age", "gender", "bmi", "pre_vision", "season", "vision_improvement_1w")
  missing_cols <- required_cols[!required_cols %in% names(period_data)]
  
  if (length(missing_cols) > 0) {
    cat(paste0("  警告: 数据缺少以下列: ", paste(missing_cols, collapse = ", "), "\n"))
    next
  }
  
  # 移除缺失值
  period_data <- period_data[complete.cases(period_data[, required_cols]), ]
  
  if (nrow(period_data) < 10) {
    cat(paste0("  警告: 时间窗口 ", current_window, " 的有效数据少于10行，可能导致模型不可靠\n"))
    next
  }
  
  # 显示样本量
  cat(paste0("  样本量: ", nrow(period_data), "\n"))
  
  # 构建完整模型 - 加入vision_improvement_1w作为预测变量
  full_model <- tryCatch({
    lm(vision_improvement_1m ~ min_rhr_steps_1 + bo_mean + steps_total + vision_improvement_1w + 
         age + gender + bmi + pre_vision + season, data = period_data)
  }, error = function(e) {
    cat(paste0("  错误: 构建完整模型失败: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(full_model)) next
  
  # 预处理分类变量
  if("season" %in% names(period_data)) {
    period_data$season <- factor(period_data$season)
    cat(paste0("    季节水平: ", paste(levels(period_data$season), collapse=", "), "\n"))
  }
  if("gender" %in% names(period_data)) {
    period_data$gender <- factor(period_data$gender)
  }
  
  # 设置随机种子以确保可重复性
  set.seed(123)
  
  # 函数计算交叉验证的R²
  cv_r2 <- function(model_formula) {
    # 使用trainControl进行交叉验证
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
        data = period_data,
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
    full_cv_r2 <- cv_r2(vision_improvement_1m ~ min_rhr_steps_1 + bo_mean + steps_total + vision_improvement_1w + 
                          age + gender + bmi + pre_vision + season)
    if(is.na(full_cv_r2)) {
      # 如果交叉验证失败，使用简单的R²作为后备
      full_cv_r2 <- summary(full_model)$r.squared
      cat("    警告: 完整模型交叉验证失败，使用R²值替代\n")
    }
    
    # HR模型
    hr_cv_r2 <- cv_r2(vision_improvement_1m ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season)
    if(is.na(hr_cv_r2)) {
      hr_model <- lm(vision_improvement_1m ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season, data = period_data)
      hr_cv_r2 <- summary(hr_model)$r.squared
      cat("    警告: HR模型交叉验证失败，使用R²值替代\n")
    }
    
    # BO模型
    bo_cv_r2 <- cv_r2(vision_improvement_1m ~ bo_mean + age + gender + bmi + pre_vision + season)
    if(is.na(bo_cv_r2)) {
      bo_model <- lm(vision_improvement_1m ~ bo_mean + age + gender + bmi + pre_vision + season, data = period_data)
      bo_cv_r2 <- summary(bo_model)$r.squared
      cat("    警告: BO模型交叉验证失败，使用R²值替代\n")
    }
    
    # 步数模型
    steps_cv_r2 <- cv_r2(vision_improvement_1m ~ steps_total + age + gender + bmi + pre_vision + season)
    if(is.na(steps_cv_r2)) {
      steps_model <- lm(vision_improvement_1m ~ steps_total + age + gender + bmi + pre_vision + season, data = period_data)
      steps_cv_r2 <- summary(steps_model)$r.squared
      cat("    警告: 步数模型交叉验证失败，使用R²值替代\n")
    }
    
    # 1周视力改善模型
    improve_1w_cv_r2 <- cv_r2(vision_improvement_1m ~ vision_improvement_1w + age + gender + bmi + pre_vision + season)
    if(is.na(improve_1w_cv_r2)) {
      improve_1w_model <- lm(vision_improvement_1m ~ vision_improvement_1w + age + gender + bmi + pre_vision + season, data = period_data)
      improve_1w_cv_r2 <- summary(improve_1w_model)$r.squared
      cat("    警告: 1周视力改善模型交叉验证失败，使用R²值替代\n")
    }
    
  }, error = function(e) {
    cat(paste0("    错误: 计算交叉验证R²失败: ", e$message, "\n"))
    # 如果交叉验证完全失败，使用普通R²
    full_cv_r2 <<- summary(full_model)$r.squared
    
    hr_model <- lm(vision_improvement_1m ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season, data = period_data)
    hr_cv_r2 <<- summary(hr_model)$r.squared
    
    bo_model <- lm(vision_improvement_1m ~ bo_mean + age + gender + bmi + pre_vision + season, data = period_data)
    bo_cv_r2 <<- summary(bo_model)$r.squared
    
    steps_model <- lm(vision_improvement_1m ~ steps_total + age + gender + bmi + pre_vision + season, data = period_data)
    steps_cv_r2 <<- summary(steps_model)$r.squared
    
    improve_1w_model <- lm(vision_improvement_1m ~ vision_improvement_1w + age + gender + bmi + pre_vision + season, data = period_data)
    improve_1w_cv_r2 <<- summary(improve_1w_model)$r.squared
    
    cat("    使用常规R²值作为替代\n")
  })
  
  # 获取当前窗口在time_windows中的索引
  window_index <- match(current_window, time_windows)
  if (!is.na(window_index)) {
    hr_performance[window_index] <- hr_cv_r2
    bo_performance[window_index] <- bo_cv_r2
    steps_performance[window_index] <- steps_cv_r2
    improve_1w_performance[window_index] <- improve_1w_cv_r2
    all_performance[window_index] <- full_cv_r2
  }
  
  # 将信息添加到性能数据框
  performance_df <- rbind(performance_df, data.frame(
    time_window = current_window,
    r_squared = summary(full_model)$r.squared,
    adj_r_squared = summary(full_model)$adj.r.squared,
    rmse = sqrt(mean(full_model$residuals^2)),
    hr_only_r2 = hr_cv_r2,
    bo_only_r2 = bo_cv_r2,
    steps_only_r2 = steps_cv_r2,
    improve_1w_only_r2 = improve_1w_cv_r2,
    all_r2 = full_cv_r2,
    sample_size = nrow(period_data),
    stringsAsFactors = FALSE
  ))
  
  cat(paste0("  完成时间窗口 ", current_window, " 的模型构建。 交叉验证R²: ", round(full_cv_r2, 4), "\n"))
}

# 排序性能数据框，按照time_windows中的顺序
performance_df$time_window <- factor(performance_df$time_window, levels = time_windows)
performance_df <- performance_df[order(performance_df$time_window), ]

# 将结果保存到CSV文件
write.csv(performance_df, file.path(output_dir, "model_performance_by_period.csv"), row.names = FALSE)
cat(paste0("结果已保存到: ", file.path(output_dir, "model_performance_by_period.csv"), "\n"))

# 创建非糖尿病组的同样分析结果
# 重置处理过的窗口
processed_windows <- character(0)

# 创建新的性能数据框用于非糖尿病组
performance_df_nod <- data.frame(
  time_window = character(),
  r_squared = numeric(),
  adj_r_squared = numeric(),
  rmse = numeric(),
  hr_only_r2 = numeric(),
  bo_only_r2 = numeric(),
  steps_only_r2 = numeric(),
  improve_1w_only_r2 = numeric(),
  all_r2 = numeric(),
  sample_size = numeric(),
  stringsAsFactors = FALSE
)

cat("\n开始为非糖尿病组分析每个时间窗口...\n")

# 打印找到的非糖尿病组文件
cat("找到的非糖尿病组(NoD)文件:\n")
for (file in nod_files) {
  cat("  ", basename(file), "\n")
}

# 非糖尿病组的分析循环 - 与糖尿病组类似
for (i in 1:length(nod_files)) {
  file_path <- nod_files[i]
  file_name <- basename(file_path)
  current_window <- get_time_window(file_path)
  
  cat("处理文件:", file_name, "\n")
  cat("  识别的时间窗口:", if(is.na(current_window)) "未识别" else current_window, "\n")
  
  if (is.na(current_window) || !current_window %in% time_windows) {
    cat("  跳过: 不在指定的时间窗口列表中\n")
    next
  }
  
  # 检查是否已处理过该窗口
  if (current_window %in% processed_windows) {
    cat("  跳过: 已经处理过", current_window, "窗口\n")
    next
  }
  
  # 记录已处理的窗口
  processed_windows <- c(processed_windows, current_window)
  
  cat(paste0("处理时间窗口: ", current_window, "\n"))
  
  # 读取数据
  period_data <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    cat(paste0("  错误: 无法读取文件: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(period_data)) next
  
  # 检查数据是否包含需要的列
  required_cols <- c("vision_improvement_1m", "min_rhr_steps_1", "bo_mean", "steps_total", 
                     "age", "gender", "bmi", "pre_vision", "season", "vision_improvement_1w")
  missing_cols <- required_cols[!required_cols %in% names(period_data)]
  
  if (length(missing_cols) > 0) {
    cat(paste0("  警告: 数据缺少以下列: ", paste(missing_cols, collapse = ", "), "\n"))
    next
  }
  
  # 移除缺失值
  period_data <- period_data[complete.cases(period_data[, required_cols]), ]
  
  if (nrow(period_data) < 10) {
    cat(paste0("  警告: 时间窗口 ", current_window, " 的有效数据少于10行，可能导致模型不可靠\n"))
    next
  }
  
  # 显示样本量
  cat(paste0("  样本量: ", nrow(period_data), "\n"))
  
  # 构建完整模型
  full_model <- tryCatch({
    lm(vision_improvement_1m ~ min_rhr_steps_1 + bo_mean + steps_total + vision_improvement_1w + 
         age + gender + bmi + pre_vision + season, data = period_data)
  }, error = function(e) {
    cat(paste0("  错误: 构建完整模型失败: ", e$message, "\n"))
    return(NULL)
  })
  
  if (is.null(full_model)) next
  
  # 预处理分类变量
  if("season" %in% names(period_data)) {
    period_data$season <- factor(period_data$season)
    cat(paste0("    季节水平: ", paste(levels(period_data$season), collapse=", "), "\n"))
  }
  if("gender" %in% names(period_data)) {
    period_data$gender <- factor(period_data$gender)
  }
  
  # 设置随机种子以确保可重复性
  set.seed(123)
  
  # 函数计算交叉验证的R²
  cv_r2 <- function(model_formula) {
    # 使用trainControl进行交叉验证
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
        data = period_data,
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
    full_cv_r2 <- cv_r2(vision_improvement_1m ~ min_rhr_steps_1 + bo_mean + steps_total + vision_improvement_1w + 
                          age + gender + bmi + pre_vision + season)
    if(is.na(full_cv_r2)) {
      # 如果交叉验证失败，使用简单的R²作为后备
      full_cv_r2 <- summary(full_model)$r.squared
      cat("    警告: 完整模型交叉验证失败，使用R²值替代\n")
    }
    
    # HR模型
    hr_cv_r2 <- cv_r2(vision_improvement_1m ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season)
    if(is.na(hr_cv_r2)) {
      hr_model <- lm(vision_improvement_1m ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season, data = period_data)
      hr_cv_r2 <- summary(hr_model)$r.squared
      cat("    警告: HR模型交叉验证失败，使用R²值替代\n")
    }
    
    # BO模型
    bo_cv_r2 <- cv_r2(vision_improvement_1m ~ bo_mean + age + gender + bmi + pre_vision + season)
    if(is.na(bo_cv_r2)) {
      bo_model <- lm(vision_improvement_1m ~ bo_mean + age + gender + bmi + pre_vision + season, data = period_data)
      bo_cv_r2 <- summary(bo_model)$r.squared
      cat("    警告: BO模型交叉验证失败，使用R²值替代\n")
    }
    
    # 步数模型
    steps_cv_r2 <- cv_r2(vision_improvement_1m ~ steps_total + age + gender + bmi + pre_vision + season)
    if(is.na(steps_cv_r2)) {
      steps_model <- lm(vision_improvement_1m ~ steps_total + age + gender + bmi + pre_vision + season, data = period_data)
      steps_cv_r2 <- summary(steps_model)$r.squared
      cat("    警告: 步数模型交叉验证失败，使用R²值替代\n")
    }
    
    # 1周视力改善模型
    improve_1w_cv_r2 <- cv_r2(vision_improvement_1m ~ vision_improvement_1w + age + gender + bmi + pre_vision + season)
    if(is.na(improve_1w_cv_r2)) {
      improve_1w_model <- lm(vision_improvement_1m ~ vision_improvement_1w + age + gender + bmi + pre_vision + season, data = period_data)
      improve_1w_cv_r2 <- summary(improve_1w_model)$r.squared
      cat("    警告: 1周视力改善模型交叉验证失败，使用R²值替代\n")
    }
    
  }, error = function(e) {
    cat(paste0("    错误: 计算交叉验证R²失败: ", e$message, "\n"))
    # 如果交叉验证完全失败，使用普通R²
    full_cv_r2 <<- summary(full_model)$r.squared
    
    hr_model <- lm(vision_improvement_1m ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season, data = period_data)
    hr_cv_r2 <<- summary(hr_model)$r.squared
    
    bo_model <- lm(vision_improvement_1m ~ bo_mean + age + gender + bmi + pre_vision + season, data = period_data)
    bo_cv_r2 <<- summary(bo_model)$r.squared
    
    steps_model <- lm(vision_improvement_1m ~ steps_total + age + gender + bmi + pre_vision + season, data = period_data)
    steps_cv_r2 <<- summary(steps_model)$r.squared
    
    improve_1w_model <- lm(vision_improvement_1m ~ vision_improvement_1w + age + gender + bmi + pre_vision + season, data = period_data)
    improve_1w_cv_r2 <<- summary(improve_1w_model)$r.squared
    
    cat("    使用常规R²值作为替代\n")
  })
  
  # 将信息添加到非糖尿病组性能数据框
  performance_df_nod <- rbind(performance_df_nod, data.frame(
    time_window = current_window,
    r_squared = summary(full_model)$r.squared,
    adj_r_squared = summary(full_model)$adj.r.squared,
    rmse = sqrt(mean(full_model$residuals^2)),
    hr_only_r2 = hr_cv_r2,
    bo_only_r2 = bo_cv_r2,
    steps_only_r2 = steps_cv_r2,
    improve_1w_only_r2 = improve_1w_cv_r2,
    all_r2 = full_cv_r2,
    sample_size = nrow(period_data),
    stringsAsFactors = FALSE
  ))
  
  cat(paste0("  完成时间窗口 ", current_window, " 的模型构建。 交叉验证R²: ", round(full_cv_r2, 4), "\n"))
}

# 排序非糖尿病组性能数据框，按照time_windows中的顺序
performance_df_nod$time_window <- factor(performance_df_nod$time_window, levels = time_windows)
performance_df_nod <- performance_df_nod[order(performance_df_nod$time_window), ]

# 将非糖尿病组结果保存到CSV文件
write.csv(performance_df_nod, file.path(output_dir, "model_performance_by_period_NoD.csv"), row.names = FALSE)
cat(paste0("非糖尿病组结果已保存到: ", file.path(output_dir, "model_performance_by_period_NoD.csv"), "\n"))

# 绘制性能变化趋势图 - 糖尿病组 (不包括1周视力改善)
plot_data_d <- data.frame(
  time_window = rep(performance_df$time_window, 4),
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
  ),
  group = rep("Diabetes", nrow(performance_df) * 4)
)

# 绘制性能变化趋势图 - 非糖尿病组 (不包括1周视力改善)
plot_data_nod <- data.frame(
  time_window = rep(performance_df_nod$time_window, 4),
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

# 为两组分别绘制图形
# 糖尿病组
p_d <- ggplot(plot_data_d, aes(x = time_window, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for 1-Month Vision Improvement - Diabetes Group",
    x = "Time window relative to surgery",
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
print(p_d)
ggsave(file.path(output_dir, "predictor_performance_1m_diabetes.pdf"), p_d, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_1m_diabetes.png"), p_d, width = 10, height = 6, dpi = 300)

# 非糖尿病组
p_nod <- ggplot(plot_data_nod, aes(x = time_window, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for 1-Month Vision Improvement - No Diabetes Group",
    x = "Time window relative to surgery",
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
print(p_nod)
ggsave(file.path(output_dir, "predictor_performance_1m_nodiabetes.pdf"), p_nod, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_1m_nodiabetes.png"), p_nod, width = 10, height = 6, dpi = 300)

# 合并两组的图 - 分面
p_combined <- ggplot(plot_data_combined, aes(x = time_window, y = r2, color = predictor, group = interaction(predictor, group))) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~ group, ncol = 1) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for 1-Month Vision Improvement by Diabetes Status",
    x = "Time window relative to surgery",
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
    panel.grid.minor = element_line(color = "gray95"),
    strip.text = element_text(size = 12, face = "bold")
  )
print(p_combined)
ggsave(file.path(output_dir, "predictor_performance_1m_combined.pdf"), p_combined, width = 10, height = 10)
ggsave(file.path(output_dir, "predictor_performance_1m_combined.png"), p_combined, width = 10, height = 10, dpi = 300)

# 创建交互式版本 (如果安装了plotly)
if (requireNamespace("plotly", quietly = TRUE)) {
  # 糖尿病组交互式图
  p_d_interactive <- plotly::ggplotly(p_d)
  htmlwidgets::saveWidget(p_d_interactive, file.path(output_dir, "predictor_performance_1m_diabetes_interactive.html"))
  
  # 非糖尿病组交互式图
  p_nod_interactive <- plotly::ggplotly(p_nod)
  htmlwidgets::saveWidget(p_nod_interactive, file.path(output_dir, "predictor_performance_1m_nodiabetes_interactive.html"))
  
  # 合并版本的交互式图
  p_combined_interactive <- plotly::ggplotly(p_combined)
  htmlwidgets::saveWidget(p_combined_interactive, file.path(output_dir, "predictor_performance_1m_combined_interactive.html"))
  
  cat("交互式图形已保存到:", file.path(output_dir), "目录\n")
}


# # 为单个预测变量创建单独的比较图
# # 把注意力集中在1周视力改善的预测力上
# p_1w_compare <- ggplot(plot_data_combined[plot_data_combined$predictor == "1W Improvement",], 
#                       aes(x = time_window, y = r2, color = group, group = group)) +
#   geom_line(size = 1.2) +
#   geom_point(size = 3.5) +
#   labs(
#     title = "1-Week Vision Improvement as Predictor for 1-Month Outcome",
#     x = "Time window relative to surgery",
#     y = "Cross-validated R²",
#     color = "Group"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     axis.title = element_text(size = 12),
#     axis.text = element_text(size = 10),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.title = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     panel.grid.major = element_line(color = "gray90"),
#     panel.grid.minor = element_line(color = "gray95")
#   )
# print(p_1w_compare)
# ggsave(file.path(output_dir, "1w_improvement_as_predictor.pdf"), p_1w_compare, width = 10, height = 6)
# ggsave(file.path(output_dir, "1w_improvement_as_predictor.png"), p_1w_compare, width = 10, height = 6, dpi = 300)

# 创建一个热图来显示不同时间窗口和预测变量的预测能力
# 糖尿病组热图
heatmap_data_d <- reshape2::melt(performance_df, 
                                 id.vars = "time_window", 
                                 measure.vars = c("hr_only_r2", "bo_only_r2", "steps_only_r2"),
                                 variable.name = "predictor", 
                                 value.name = "r2")

# 修改预测变量名称以便更好地显示
heatmap_data_d$predictor <- factor(heatmap_data_d$predictor,
                                   levels = c("hr_only_r2", "bo_only_r2", "steps_only_r2"),
                                   labels = c("Heart Rate", "Blood Oxygen", "Total Steps"))

p_heatmap_d <- ggplot(heatmap_data_d, aes(x = time_window, y = predictor, fill = r2)) +
  geom_tile() +
  scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                       midpoint = median(heatmap_data_d$r2, na.rm = TRUE)) +
  labs(
    title = "Predictive Power for 1-Month Vision Outcome - Diabetes Group",
    x = "Time Window",
    y = "Predictor",
    fill = "R²"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p_heatmap_d)
ggsave(file.path(output_dir, "predictor_heatmap_diabetes.pdf"), p_heatmap_d, width = 10, height = 6)

# 非糖尿病组热图
heatmap_data_nod <- reshape2::melt(performance_df_nod, 
                                   id.vars = "time_window", 
                                   measure.vars = c("hr_only_r2", "bo_only_r2", "steps_only_r2"),
                                   variable.name = "predictor", 
                                   value.name = "r2")

heatmap_data_nod$predictor <- factor(heatmap_data_nod$predictor,
                                     levels = c("hr_only_r2", "bo_only_r2", "steps_only_r2"),
                                     labels = c("Heart Rate", "Blood Oxygen", "Total Steps"))

p_heatmap_nod <- ggplot(heatmap_data_nod, aes(x = time_window, y = predictor, fill = r2)) +
  geom_tile() +
  scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                       midpoint = median(heatmap_data_nod$r2, na.rm = TRUE)) +
  labs(
    title = "Predictive Power for 1-Month Vision Outcome - No Diabetes Group",
    x = "Time Window",
    y = "Predictor",
    fill = "R²"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p_heatmap_nod)
ggsave(file.path(output_dir, "predictor_heatmap_nodiabetes.pdf"), p_heatmap_nod, width = 10, height = 6)

# 创建组合热图
heatmap_data_d$group <- "Diabetes"
heatmap_data_nod$group <- "No Diabetes"
heatmap_data_combined <- rbind(heatmap_data_d, heatmap_data_nod)

p_heatmap_combined <- ggplot(heatmap_data_combined, aes(x = time_window, y = predictor, fill = r2)) +
  geom_tile() +
  facet_wrap(~ group) +
  scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                       midpoint = median(heatmap_data_combined$r2, na.rm = TRUE)) +
  labs(
    title = "Predictive Power for 1-Month Vision Outcome by Diabetes Status",
    x = "Time Window",
    y = "Predictor",
    fill = "R²"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )
print(p_heatmap_combined)
ggsave(file.path(output_dir, "predictor_heatmap_combined.pdf"), p_heatmap_combined, width = 12, height = 6)
ggsave(file.path(output_dir, "predictor_heatmap_combined.png"), p_heatmap_combined, width = 12, height = 6, dpi = 300)

# 总结分析结果
cat("\n=== 1个月视力预测模型性能分析摘要 ===\n")

# 找出每组中性能最好的时间窗口
best_window_d <- performance_df[which.max(performance_df$all_r2), ]
best_window_nod <- performance_df_nod[which.max(performance_df_nod$all_r2), ]

cat("\n糖尿病组:\n")
cat("- 预测性能最好的时间窗口: ", as.character(best_window_d$time_window), 
    " (R² = ", round(best_window_d$all_r2, 4), ")\n", sep = "")

# 找出糖尿病组中最佳的单一预测变量 (不包括1周视力改善)
best_physio_vars_d <- c(performance_df$hr_only_r2, performance_df$bo_only_r2, performance_df$steps_only_r2)
best_physio_idx_d <- which.max(best_physio_vars_d)
best_physio_names_d <- c("HR", "BO", "Total Steps")
best_physio_windows_d <- performance_df$time_window[ceiling(best_physio_idx_d/3)]
cat("- 最佳生理预测变量: ", best_physio_names_d[best_physio_idx_d %% 3 + 1], 
    " (R² = ", round(best_physio_vars_d[best_physio_idx_d], 4), 
    ", 时间窗口: ", as.character(best_physio_windows_d), ")\n", sep = "")

cat("\n非糖尿病组:\n")
cat("- 预测性能最好的时间窗口: ", as.character(best_window_nod$time_window), 
    " (R² = ", round(best_window_nod$all_r2, 4), ")\n", sep = "")

# 找出非糖尿病组中最佳的单一预测变量 (不包括1周视力改善)
best_physio_vars_nod <- c(performance_df_nod$hr_only_r2, performance_df_nod$bo_only_r2, performance_df_nod$steps_only_r2)
best_physio_idx_nod <- which.max(best_physio_vars_nod)
best_physio_names_nod <- c("HR", "BO", "Total Steps")
best_physio_windows_nod <- performance_df_nod$time_window[ceiling(best_physio_idx_nod/3)]
cat("- 最佳生理预测变量: ", best_physio_names_nod[best_physio_idx_nod %% 3 + 1], 
    " (R² = ", round(best_physio_vars_nod[best_physio_idx_nod], 4), 
    ", 时间窗口: ", as.character(best_physio_windows_nod), ")\n", sep = "")

# 总体结论
cat("\n总体结论:\n")
if (best_window_d$all_r2 > best_window_nod$all_r2) {
  cat("- 糖尿病组的1个月视力预测模型整体表现优于非糖尿病组\n")
} else {
  cat("- 非糖尿病组的1个月视力预测模型整体表现优于糖尿病组\n")
}

# 比较生理指标的平均预测性能
avg_physio_d <- mean(c(performance_df$hr_only_r2, performance_df$bo_only_r2, performance_df$steps_only_r2), na.rm = TRUE)
avg_physio_nod <- mean(c(performance_df_nod$hr_only_r2, performance_df_nod$bo_only_r2, performance_df_nod$steps_only_r2), na.rm = TRUE)

cat("- 糖尿病组生理指标平均预测力: R² = ", round(avg_physio_d, 4), "\n", sep = "")
cat("- 非糖尿病组生理指标平均预测力: R² = ", round(avg_physio_nod, 4), "\n", sep = "")

