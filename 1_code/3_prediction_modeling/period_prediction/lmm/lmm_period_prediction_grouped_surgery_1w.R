library(tidyverse)
library(lme4)
library(caret)
library(ggplot2)
setwd(get_project_wd())
rm(list = ls())

# 设置工作路径
input_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/time_period_data_grouped"
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/time_period_data_grouped/lmm_model_performance/surgery_group"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义标准时间窗口映射
time_window_mapping <- list(
  "pre_7d_all" = "pre_7d_all",
  "pre_3d_7d" = "pre_3d_7d", 
  "pre_3d" = "pre_3d",
  "post_1to3d" = "post_1to3d",
  "post_4to6d" = "post_4to6d", 
  "post_6d" = "post_6d"
)

# 定义时间窗口
time_windows <- names(time_window_mapping)

# 获取所有CSV文件
all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# 将文件分类为糖尿病组和非糖尿病组
d_surg1_files <- grep("D_Surg1\\.csv$", all_files, value = TRUE)
nod_surg0_files <- grep("NoD_Surg0\\.csv$", all_files, value = TRUE)

cat("找到的D_Surg1文件数量:", length(d_surg1_files), "\n")
cat("找到的NoD_Surg0文件数量:", length(nod_surg0_files), "\n")

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
  } else if (grepl("^post_1to3d_[DN]", file_name)) {
    return("post_1to3d")
  } else if (grepl("^post_4to6d_[DN]", file_name)) {
    return("post_4to6d")
  } else if (grepl("^post_6d_[DN]", file_name)) {
    return("post_6d")
  }
  
  # 如果没有匹配到预定义窗口，返回NA
  cat("警告: 无法识别文件", file_name, "的时间窗口\n")
  return(NA)
}

# 定义一个函数来处理每个组
process_group <- function(files, group_name) {
  # 创建数据框来存储每个时间窗口模型的性能指标
  performance_df <- data.frame(
    time_window = character(),
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
  
  # 基于文件名跟踪我们已经处理过的窗口
  processed_windows <- character(0)
  
  # 首先打印一下所有找到的文件，便于调试
  cat(paste0("找到的", group_name, "文件:\n"))
  for (file in files) {
    cat("  ", basename(file), "\n")
  }
  
  # 对每个时间窗口构建模型
  cat(paste0("开始为", group_name, "的每个时间窗口构建预测模型...\n"))
  
  for (i in 1:length(files)) {
    file_path <- files[i]
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
    required_cols <- c("vision_improvement", "min_rhr_steps_1", "bo_mean", "steps_total", 
                       "age", "gender", "bmi", "pre_vision", "season")
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
      lm(vision_improvement ~ min_rhr_steps_1 + bo_mean + steps_total + 
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
      full_cv_r2 <- cv_r2(vision_improvement ~ min_rhr_steps_1 + bo_mean + steps_total + 
                            age + gender + bmi + pre_vision + season)
      if(is.na(full_cv_r2)) {
        # 如果交叉验证失败，使用简单的R²作为后备
        full_cv_r2 <- summary(full_model)$r.squared
        cat("    警告: 完整模型交叉验证失败，使用R²值替代\n")
      }
      
      # HR模型
      hr_cv_r2 <- cv_r2(vision_improvement ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season)
      if(is.na(hr_cv_r2)) {
        hr_model <- lm(vision_improvement ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season, data = period_data)
        hr_cv_r2 <- summary(hr_model)$r.squared
        cat("    警告: HR模型交叉验证失败，使用R²值替代\n")
      }
      
      # BO模型
      bo_cv_r2 <- cv_r2(vision_improvement ~ bo_mean + age + gender + bmi + pre_vision + season)
      if(is.na(bo_cv_r2)) {
        bo_model <- lm(vision_improvement ~ bo_mean + age + gender + bmi + pre_vision + season, data = period_data)
        bo_cv_r2 <- summary(bo_model)$r.squared
        cat("    警告: BO模型交叉验证失败，使用R²值替代\n")
      }
      
      # 步数模型
      steps_cv_r2 <- cv_r2(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season)
      if(is.na(steps_cv_r2)) {
        steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, data = period_data)
        steps_cv_r2 <- summary(steps_model)$r.squared
        cat("    警告: 步数模型交叉验证失败，使用R²值替代\n")
      }
    }, error = function(e) {
      cat(paste0("    错误: 计算交叉验证R²失败: ", e$message, "\n"))
      # 如果交叉验证完全失败，使用普通R²
      full_cv_r2 <<- summary(full_model)$r.squared
      
      hr_model <- lm(vision_improvement ~ min_rhr_steps_1 + age + gender + bmi + pre_vision + season, data = period_data)
      hr_cv_r2 <<- summary(hr_model)$r.squared
      
      bo_model <- lm(vision_improvement ~ bo_mean + age + gender + bmi + pre_vision + season, data = period_data)
      bo_cv_r2 <<- summary(bo_model)$r.squared
      
      steps_model <- lm(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season, data = period_data)
      steps_cv_r2 <<- summary(steps_model)$r.squared
      
      cat("    使用常规R²值作为替代\n")
    })
    
    # 将信息添加到性能数据框
    performance_df <- rbind(performance_df, data.frame(
      time_window = current_window,
      r_squared = summary(full_model)$r.squared,
      adj_r_squared = summary(full_model)$adj.r.squared,
      rmse = sqrt(mean(full_model$residuals^2)),
      hr_only_r2 = hr_cv_r2,
      bo_only_r2 = bo_cv_r2,
      steps_only_r2 = steps_cv_r2,
      all_r2 = full_cv_r2,
      sample_size = nrow(period_data),
      stringsAsFactors = FALSE
    ))
    
    cat(paste0("  完成时间窗口 ", current_window, " 的模型构建。 交叉验证R²: ", round(full_cv_r2, 4), "\n"))
  }
  
  # 排序性能数据框，按照time_windows中的顺序
  performance_df$time_window <- factor(performance_df$time_window, levels = time_windows)
  performance_df <- performance_df[order(performance_df$time_window), ]
  
  return(performance_df)
}

# 处理D_Surg1组
performance_df_d <- process_group(d_surg1_files, "D_Surg1")

# 处理NoD_Surg0组
performance_df_nod <- process_group(nod_surg0_files, "NoD_Surg0")

# 将结果保存到CSV文件
write.csv(performance_df_d, file.path(output_dir, "model_performance_by_period_d_surg1.csv"), row.names = FALSE)
write.csv(performance_df_nod, file.path(output_dir, "model_performance_by_period_nod_surg0.csv"), row.names = FALSE)

cat(paste0("D_Surg1组结果已保存到: ", file.path(output_dir, "model_performance_by_period_d_surg1.csv"), "\n"))
cat(paste0("NoD_Surg0组结果已保存到: ", file.path(output_dir, "model_performance_by_period_nod_surg0.csv"), "\n"))

# 为D_Surg1组创建绘图数据
plot_data_d <- data.frame(
  time_window = rep(performance_df_d$time_window, 4),
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
  )
)

# 为NoD_Surg0组创建绘图数据
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
  )
)

# 绘制D_Surg1组图形
p_d <- ggplot(plot_data_d, aes(x = time_window, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for Vision Improvement - D_Surg1 Group",
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

# 保存D_Surg1组图形
print(p_d)
ggsave(file.path(output_dir, "predictor_performance_by_period_d_surg1.pdf"), p_d, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_by_period_d_surg1.png"), p_d, width = 10, height = 6, dpi = 300)

# 绘制NoD_Surg0组图形
p_nod <- ggplot(plot_data_nod, aes(x = time_window, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for Vision Improvement - NoD_Surg0 Group",
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

# 保存NoD_Surg0组图形
print(p_nod)
ggsave(file.path(output_dir, "predictor_performance_by_period_nod_surg0.pdf"), p_nod, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_by_period_nod_surg0.png"), p_nod, width = 10, height = 6, dpi = 300)
