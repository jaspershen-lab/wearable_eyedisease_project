library(tidyverse)
library(lme4)
library(caret)
library(ggplot2)
setwd(get_project_wd())
rm(list = ls())

# 设置工作路径（请根据需要调整）
input_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped"
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped/lmm_model_performance/surgery_group"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 获取所有分组的日期文件列表
d_surg1_files <- list.files(input_dir, pattern = "D_Surg1\\.csv$", full.names = TRUE)
nod_surg0_files <- list.files(input_dir, pattern = "NoD_Surg0\\.csv$", full.names = TRUE)

# 调试: 检查找到的文件
cat("找到的D_Surg1文件数量:", length(d_surg1_files), "\n")
cat("找到的NoD_Surg0文件数量:", length(nod_surg0_files), "\n")

# 将文件名转换为时间点列表
get_time_point <- function(file_path) {
  time_str <- gsub(".*day_(.*)_[DN].*\\.csv$", "\\1", basename(file_path))
  return(time_str)
}

time_points_d <- sapply(d_surg1_files, get_time_point)
time_points_nod <- sapply(nod_surg0_files, get_time_point)

# 获取两组共有的时间点
time_points <- intersect(time_points_d, time_points_nod)

# 调试: 显示找到的时间点
cat("找到的时间点数量:", length(time_points), "\n")
if(length(time_points) > 0) {
  cat("前5个时间点:", head(time_points, 5), "\n")
}

# 对时间点进行排序（从-7到6）
time_points_num <- as.numeric(time_points)
time_points_str <- as.character(time_points[order(time_points_num)])

# 调试: 显示排序后的时间点
cat("排序后时间点:", head(time_points_str, 10), "...\n")

# 定义一个函数来处理一组数据
analyze_group <- function(files, group_name, time_points_str) {
  # 创建数据框来存储每个时间点模型的性能指标
  performance_df <- data.frame(
    time_point = character(),
    group = character(),
    r_squared = numeric(),
    adj_r_squared = numeric(),
    rmse = numeric(),
    hr_only_r2 = numeric(),
    bo_only_r2 = numeric(),
    steps_only_r2 = numeric(),
    all_r2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 存储每个预测变量单独的性能
  hr_performance <- numeric(length(time_points_str))
  bo_performance <- numeric(length(time_points_str))
  steps_performance <- numeric(length(time_points_str))
  all_performance <- numeric(length(time_points_str))
  
  # 存储模型系数
  coef_summary <- data.frame(
    time_point = character(),
    group = character(),
    predictor = character(),
    estimate = numeric(),
    std_error = numeric(),
    t_value = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 对每个时间点构建模型
  cat(paste0("开始为", group_name, "组的每个时间点构建预测模型...\n"))
  
  for (i in 1:length(time_points_str)) {
    tp <- time_points_str[i]
    cat(paste0("处理", group_name, "组的时间点: ", tp, "\n"))
    
    # 确定文件模式
    file_pattern <- if(group_name == "D_Surg1") {
      paste0("day_", tp, "_D_Surg1\\.csv$")
    } else {
      paste0("day_", tp, "_NoD_Surg0\\.csv$")
    }
    
    file_path <- grep(file_pattern, files, value = TRUE)
    
    if (length(file_path) == 0) {
      cat(paste0("  警告: 没有找到", group_name, "组时间点 ", tp, " 的数据文件 (使用模式: ", file_pattern, ")\n"))
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
      cat(paste0("  警告: ", group_name, "组时间点 ", tp, " 的有效数据少于10行，可能导致模型不可靠\n"))
      next
    }
    
    # 确保所有分类变量都被转换为因子
    if("season" %in% names(day_data)) {
      day_data$season <- factor(day_data$season)
      cat(paste0("    季节水平: ", paste(levels(day_data$season), collapse=", "), "\n"))
    }
    if("gender" %in% names(day_data)) {
      day_data$gender <- factor(day_data$gender)
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
    
    # 设置随机种子以确保可重复性
    set.seed(123)
    
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
        hr_cv_r2 <- summary(hr_model)$r.squared
        cat("    警告: HR模型交叉验证失败，使用R²值替代\n")
      }
      
      # BO模型
      bo_cv_r2 <- cv_r2(vision_improvement ~ mean_bo + age + gender + bmi + pre_vision + season)
      if(is.na(bo_cv_r2)) {
        bo_cv_r2 <- summary(bo_model)$r.squared
        cat("    警告: BO模型交叉验证失败，使用R²值替代\n")
      }
      
      # 步数模型
      steps_cv_r2 <- cv_r2(vision_improvement ~ steps_total + age + gender + bmi + pre_vision + season)
      if(is.na(steps_cv_r2)) {
        steps_cv_r2 <- summary(steps_model)$r.squared
        cat("    警告: 步数模型交叉验证失败，使用R²值替代\n")
      }
    }, error = function(e) {
      cat(paste0("    错误: 计算交叉验证R²失败: ", e$message, "\n"))
      # 如果交叉验证完全失败，使用普通R²
      full_cv_r2 <<- summary(full_model)$r.squared
      hr_cv_r2 <<- summary(hr_model)$r.squared
      bo_cv_r2 <<- summary(bo_model)$r.squared
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
      group = group_name,
      r_squared = summary(full_model)$r.squared,
      adj_r_squared = summary(full_model)$adj.r.squared,
      rmse = sqrt(mean(full_model$residuals^2)),
      hr_only_r2 = hr_cv_r2,
      bo_only_r2 = bo_cv_r2,
      steps_only_r2 = steps_cv_r2,
      all_r2 = full_cv_r2,
      stringsAsFactors = FALSE
    ))
    
    # 提取系数信息
    model_summary <- summary(full_model)
    coefs <- model_summary$coefficients
    
    for (j in 1:nrow(coefs)) {
      predictor_name <- rownames(coefs)[j]
      
      # 只关注生理指标
      if (predictor_name %in% c("mean_rhr_1", "mean_bo", "steps_total")) {
        coef_summary <- rbind(coef_summary, data.frame(
          time_point = tp,
          group = group_name,
          predictor = predictor_name,
          estimate = coefs[j, 1],
          std_error = coefs[j, 2],
          t_value = coefs[j, 3],
          p_value = coefs[j, 4],
          stringsAsFactors = FALSE
        ))
      }
    }
    
    cat(paste0("  完成", group_name, "组时间点 ", tp, " 的模型构建。 交叉验证R²: ", round(full_cv_r2, 4), "\n"))
  }
  
  # 处理时间点标签
  performance_df$time_point_num <- as.numeric(performance_df$time_point)
  performance_df <- performance_df[order(performance_df$time_point_num), ]
  
  # 处理系数时间点
  coef_summary$time_point_num <- as.numeric(coef_summary$time_point)
  coef_summary <- coef_summary[order(coef_summary$time_point_num), ]
  
  # 返回两个数据框
  return(list(performance = performance_df, coefficients = coef_summary))
}

# 分别分析D_Surg1组和NoD_Surg0组
d_surg1_results <- analyze_group(d_surg1_files, "D_Surg1", time_points_str)
nod_surg0_results <- analyze_group(nod_surg0_files, "NoD_Surg0", time_points_str)

# 将结果保存到CSV文件
write.csv(d_surg1_results$performance, file.path(output_dir, "model_performance_d_surg1.csv"), row.names = FALSE)
write.csv(nod_surg0_results$performance, file.path(output_dir, "model_performance_nod_surg0.csv"), row.names = FALSE)

cat(paste0("D_Surg1组结果已保存到: ", file.path(output_dir, "model_performance_d_surg1.csv"), "\n"))
cat(paste0("NoD_Surg0组结果已保存到: ", file.path(output_dir, "model_performance_nod_surg0.csv"), "\n"))

# 为D_Surg1组创建绘图数据
plot_data_d <- data.frame(
  time_point = rep(d_surg1_results$performance$time_point, 4),
  time_point_num = rep(d_surg1_results$performance$time_point_num, 4),
  predictor = c(
    rep("All Predictors", nrow(d_surg1_results$performance)),
    rep("HR", nrow(d_surg1_results$performance)),
    rep("BO", nrow(d_surg1_results$performance)),
    rep("Total Steps", nrow(d_surg1_results$performance))
  ),
  r2 = c(
    d_surg1_results$performance$all_r2,
    d_surg1_results$performance$hr_only_r2,
    d_surg1_results$performance$bo_only_r2,
    d_surg1_results$performance$steps_only_r2
  )
)

# 为NoD_Surg0组创建绘图数据
plot_data_nod <- data.frame(
  time_point = rep(nod_surg0_results$performance$time_point, 4),
  time_point_num = rep(nod_surg0_results$performance$time_point_num, 4),
  predictor = c(
    rep("All Predictors", nrow(nod_surg0_results$performance)),
    rep("HR", nrow(nod_surg0_results$performance)),
    rep("BO", nrow(nod_surg0_results$performance)),
    rep("Total Steps", nrow(nod_surg0_results$performance))
  ),
  r2 = c(
    nod_surg0_results$performance$all_r2,
    nod_surg0_results$performance$hr_only_r2,
    nod_surg0_results$performance$bo_only_r2,
    nod_surg0_results$performance$steps_only_r2
  )
)

# 绘制D_Surg1组性能图
p_d <- ggplot(plot_data_d, aes(x = time_point_num, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for Vision Improvement - D_Surg1 Group",
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

# 绘制NoD_Surg0组性能图
p_nod <- ggplot(plot_data_nod, aes(x = time_point_num, y = r2, color = predictor, group = predictor)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("All Predictors" = "#513a52", 
                                "HR" = "#82a1bf", 
                                "BO" = "#faaa93", 
                                "Total Steps" = "#feefc4")) +
  labs(
    title = "Predictor Performance for Vision Improvement - NoD_Surg0 Group",
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

# 显示图形
print(p_d)
print(p_nod)

# 保存图形
ggsave(file.path(output_dir, "predictor_performance_d_surg1.pdf"), p_d, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_d_surg1.png"), p_d, width = 10, height = 6, dpi = 300)

ggsave(file.path(output_dir, "predictor_performance_nod_surg0.pdf"), p_nod, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_nod_surg0.png"), p_nod, width = 10, height = 6, dpi = 300)

# 创建交互式版本 (如果安装了plotly)
if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
  p_d_interactive <- plotly::ggplotly(p_d)
  htmlwidgets::saveWidget(p_d_interactive, file.path(output_dir, "predictor_performance_d_surg1_interactive.html"))
  
  p_nod_interactive <- plotly::ggplotly(p_nod)
  htmlwidgets::saveWidget(p_nod_interactive, file.path(output_dir, "predictor_performance_nod_surg0_interactive.html"))
  
  cat("交互式图形已保存到:", file.path(output_dir), "目录\n")
}
