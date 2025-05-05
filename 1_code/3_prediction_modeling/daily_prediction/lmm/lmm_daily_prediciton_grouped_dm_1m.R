library(tidyverse)
library(lme4)
library(caret)
library(ggplot2)
setwd(get_project_wd())
rm(list = ls())

# 设置工作路径（请根据需要调整）
input_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/daily_data/dm_group"
output_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/daily_data/dm_group/lmm_model_performance"

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 获取所有分组的日期文件列表 - 根据新的文件命名模式更新
d_files <- list.files(input_dir, pattern = "_D\\.csv$", full.names = TRUE)
nod_files <- list.files(input_dir, pattern = "_NoD\\.csv$", full.names = TRUE)

# 调试: 检查找到的文件
cat("找到的D文件数量:", length(d_files), "\n")
cat("找到的NoD文件数量:", length(nod_files), "\n")

# 更精确的函数来提取时间点
get_time_point <- function(file_path) {
  # 从文件名中提取日期部分
  file_name <- basename(file_path)
  # 使用更精确的模式来捕获正数和负数 - 更新为新的文件命名模式
  matches <- regexpr("day_(-?[0-9]+)_[DN]", file_name)
  if(matches != -1) {
    start_pos <- matches + 4  # "day_"的长度
    length_match <- attr(matches, "match.length") - 4 - 2  # 减去"day_"和"_D"或"_N"
    time_str <- substr(file_name, start_pos, start_pos + length_match - 1)
    return(time_str)
  } else {
    warning("无法从文件名提取时间点: ", file_name)
    return(NA)
  }
}

# 测试函数是否正常工作
if(length(d_files) > 0) {
  cat("测试时间点提取:\n")
  for(i in 1:min(5, length(d_files))) {
    cat("  文件:", basename(d_files[i]), "-> 提取时间点:", get_time_point(d_files[i]), "\n")
  }
}

time_points_d <- sapply(d_files, get_time_point)
time_points_nod <- sapply(nod_files, get_time_point)

# 移除可能的NA值
time_points_d <- time_points_d[!is.na(time_points_d)]
time_points_nod <- time_points_nod[!is.na(time_points_nod)]

# 获取两组共有的时间点
time_points <- intersect(time_points_d, time_points_nod)

# 调试: 显示找到的时间点
cat("找到的时间点数量:", length(time_points), "\n")
if(length(time_points) > 0) {
  cat("前5个时间点:", head(time_points, 5), "\n")
}

# 对时间点进行排序（从-7到29）
time_points_num <- as.numeric(time_points)
time_points_str <- time_points[order(time_points_num)]

# 调试: 显示排序后的时间点
cat("排序后时间点:", head(time_points_str, 10), "...\n")

# 定义一个函数来处理一组数据（D或NoD）
analyze_group <- function(files, group_name, time_points_str) {
  # Create dataframe to store model performance metrics for each time point
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
    has_surgery_covariate = logical(),
    stringsAsFactors = FALSE
  )
  
  # Store individual predictor performance
  hr_performance <- numeric(length(time_points_str))
  bo_performance <- numeric(length(time_points_str))
  steps_performance <- numeric(length(time_points_str))
  all_performance <- numeric(length(time_points_str))
  
  # Store model coefficients
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
  
  # For each time point, build models
  cat(paste0("开始为", group_name, "组的每个时间点构建预测模型...\n"))
  
  for (i in 1:length(time_points_str)) {
    tp <- time_points_str[i]
    cat(paste0("处理", group_name, "组的时间点: ", tp, "\n"))
    
    # Match file pattern based on time point
    if(startsWith(tp, "-")) {
      file_pattern <- paste0("day_\\", tp, "_", ifelse(group_name == "D", "D", "NoD"), "\\.csv$")
    } else {
      file_pattern <- paste0("day_", tp, "_", ifelse(group_name == "D", "D", "NoD"), "\\.csv$")
    }
    
    file_path <- grep(file_pattern, files, value = TRUE)
    
    if (length(file_path) == 0) {
      cat(paste0("  警告: 没有找到", group_name, "组时间点 ", tp, " 的数据文件 (使用模式: ", file_pattern, ")\n"))
      next
    }
    
    # Read data
    day_data <- tryCatch({
      read.csv(file_path)
    }, error = function(e) {
      cat(paste0("  错误: 无法读取文件: ", e$message, "\n"))
      return(NULL)
    })
    
    if (is.null(day_data)) next
    
    # Check required columns
    required_cols <- c("vision_improvement_1m", "mean_rhr_1", "mean_bo", "steps_total", 
                       "age", "gender", "bmi", "pre_vision", "season", "vision_improvement_1w")
    
    # Only add surgery_type to required columns if we're processing D group
    # This avoids the issue with NoD group having only one level in surgery_type
    if(group_name == "D") {
      required_cols <- c(required_cols, "surgery_type")
    }
    
    # Handle alternative surgery_type column name if needed
    if(!("surgery_type" %in% names(day_data)) && ("surgery_1..0.PI.1.other." %in% names(day_data))) {
      day_data$surgery_type <- day_data$surgery_1..0.PI.1.other.
      cat("  使用 surgery_1..0.PI.1.other. 列作为 surgery_type\n")
    }
    
    missing_cols <- required_cols[!required_cols %in% names(day_data)]
    
    if (length(missing_cols) > 0) {
      cat(paste0("  警告: 数据缺少以下列: ", paste(missing_cols, collapse = ", "), "\n"))
      # If only surgery_type is missing, continue but without it in the model
      if(length(missing_cols) == 1 && missing_cols == "surgery_type") {
        cat("  继续处理，但模型将不包含手术类型作为协变量\n")
        required_cols <- required_cols[required_cols != "surgery_type"]
      } else {
        next
      }
    }
    
    # Remove rows with missing values
    day_data <- day_data[complete.cases(day_data[, required_cols]), ]
    
    if (nrow(day_data) < 10) {
      cat(paste0("  警告: ", group_name, "组时间点 ", tp, " 的有效数据少于10行，可能导致模型不可靠\n"))
      next
    }
    
    # Convert categorical variables to factors
    if("season" %in% names(day_data)) {
      day_data$season <- factor(day_data$season)
      cat(paste0("    季节水平: ", paste(levels(day_data$season), collapse=", "), "\n"))
    }
    if("gender" %in% names(day_data)) {
      day_data$gender <- factor(day_data$gender)
    }
    
    # Only use surgery_type as a factor for D group or if it has multiple levels
    use_surgery_as_factor <- FALSE
    if("surgery_type" %in% names(day_data)) {
      day_data$surgery_type <- factor(day_data$surgery_type)
      cat(paste0("    手术类型水平: ", paste(levels(day_data$surgery_type), collapse=", "), "\n"))
      
      # Only use surgery_type in model if it has multiple levels
      if(length(levels(day_data$surgery_type)) > 1) {
        use_surgery_as_factor <- TRUE
      } else {
        cat("    警告: 手术类型仅有一个水平，将从模型中排除\n")
      }
    }
    
    # Build full model with or without surgery_type based on factor levels
    if(use_surgery_as_factor) {
      full_model_formula <- vision_improvement_1m ~ mean_rhr_1 + mean_bo + steps_total + 
        age + gender + bmi + pre_vision + season + vision_improvement_1w + surgery_type
    } else {
      full_model_formula <- vision_improvement_1m ~ mean_rhr_1 + mean_bo + steps_total + 
        age + gender + bmi + pre_vision + season + vision_improvement_1w
    }
    
    full_model <- tryCatch({
      lm(full_model_formula, data = day_data)
    }, error = function(e) {
      cat(paste0("  错误: 构建完整模型失败: ", e$message, "\n"))
      return(NULL)
    })
    
    if (is.null(full_model)) next
    
    # Build single predictor models (with or without surgery_type)
    if(use_surgery_as_factor) {
      hr_model_formula <- vision_improvement_1m ~ mean_rhr_1 + age + gender + bmi + pre_vision + season + vision_improvement_1w + surgery_type
      bo_model_formula <- vision_improvement_1m ~ mean_bo + age + gender + bmi + pre_vision + season + vision_improvement_1w + surgery_type
      steps_model_formula <- vision_improvement_1m ~ steps_total + age + gender + bmi + pre_vision + season + vision_improvement_1w + surgery_type
    } else {
      hr_model_formula <- vision_improvement_1m ~ mean_rhr_1 + age + gender + bmi + pre_vision + season + vision_improvement_1w
      bo_model_formula <- vision_improvement_1m ~ mean_bo + age + gender + bmi + pre_vision + season + vision_improvement_1w
      steps_model_formula <- vision_improvement_1m ~ steps_total + age + gender + bmi + pre_vision + season + vision_improvement_1w
    }
    
    hr_model <- lm(hr_model_formula, data = day_data)
    bo_model <- lm(bo_model_formula, data = day_data)
    steps_model <- lm(steps_model_formula, data = day_data)
    
    # Set random seed for reproducibility
    set.seed(123)
    
    # Function to calculate cross-validated R²
    cv_r2 <- function(model_formula) {
      # Use trainControl for cross-validation
      train_control <- trainControl(
        method = "cv",
        number = 5,
        savePredictions = "final"
      )
      
      # Use caret's train function for cross-validation
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
      
      # Return cross-validated R²
      return(cv_model$results$Rsquared)
    }
    
    # Calculate cross-validated R² for each model
    cat("    计算交叉验证R²值...\n")
    
    # Handle potential errors and provide fallback
    tryCatch({
      # Full model
      full_cv_r2 <- cv_r2(full_model_formula)
      if(is.na(full_cv_r2)) {
        # If cross-validation fails, use regular R² as fallback
        full_cv_r2 <- summary(full_model)$r.squared
        cat("    警告: 完整模型交叉验证失败，使用R²值替代\n")
      }
      
      # HR model
      hr_cv_r2 <- cv_r2(hr_model_formula)
      if(is.na(hr_cv_r2)) {
        hr_cv_r2 <- summary(hr_model)$r.squared
        cat("    警告: HR模型交叉验证失败，使用R²值替代\n")
      }
      
      # BO model
      bo_cv_r2 <- cv_r2(bo_model_formula)
      if(is.na(bo_cv_r2)) {
        bo_cv_r2 <- summary(bo_model)$r.squared
        cat("    警告: BO模型交叉验证失败，使用R²值替代\n")
      }
      
      # Steps model
      steps_cv_r2 <- cv_r2(steps_model_formula)
      if(is.na(steps_cv_r2)) {
        steps_cv_r2 <- summary(steps_model)$r.squared
        cat("    警告: 步数模型交叉验证失败，使用R²值替代\n")
      }
    }, error = function(e) {
      cat(paste0("    错误: 计算交叉验证R²失败: ", e$message, "\n"))
      # If cross-validation completely fails, use regular R²
      full_cv_r2 <<- summary(full_model)$r.squared
      hr_cv_r2 <<- summary(hr_model)$r.squared
      bo_cv_r2 <<- summary(bo_model)$r.squared
      steps_cv_r2 <<- summary(steps_model)$r.squared
      
      cat("    使用常规R²值作为替代\n")
    })
    
    # Store performance metrics
    hr_performance[i] <- hr_cv_r2
    bo_performance[i] <- bo_cv_r2
    steps_performance[i] <- steps_cv_r2
    all_performance[i] <- full_cv_r2
    
    # Add information to performance dataframe
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
      has_surgery_covariate = use_surgery_as_factor,
      stringsAsFactors = FALSE
    ))
    
    # Extract coefficient information
    model_summary <- summary(full_model)
    coefs <- model_summary$coefficients
    
    for (j in 1:nrow(coefs)) {
      predictor_name <- rownames(coefs)[j]
      
      # Focus only on physiological indicators and surgery_type
      if (predictor_name %in% c("mean_rhr_1", "mean_bo", "steps_total") || 
          grepl("surgery_type", predictor_name)) {
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
  
  # Process time point labels
  performance_df$time_point_num <- as.numeric(performance_df$time_point)
  performance_df <- performance_df[order(performance_df$time_point_num), ]
  
  # Process coefficient time points
  coef_summary$time_point_num <- as.numeric(coef_summary$time_point)
  coef_summary <- coef_summary[order(coef_summary$time_point_num), ]
  
  # Return both dataframes
  return(list(performance = performance_df, coefficients = coef_summary))
}

# 分别分析D组和NoD组
d_results <- analyze_group(d_files, "D", time_points_str)
nod_results <- analyze_group(nod_files, "NoD", time_points_str)

# 合并结果
all_performance <- rbind(d_results$performance, nod_results$performance)
all_coefficients <- rbind(d_results$coefficients, nod_results$coefficients)

# 将结果保存到CSV文件
write.csv(all_performance, file.path(output_dir, "model_performance_both_groups.csv"), row.names = FALSE)
write.csv(all_coefficients, file.path(output_dir, "coefficient_summary_both_groups.csv"), row.names = FALSE)

cat(paste0("结果已保存到: ", file.path(output_dir, "model_performance_both_groups.csv"), "\n"))
cat(paste0("系数已保存到: ", file.path(output_dir, "coefficient_summary_both_groups.csv"), "\n"))

# 为D组创建绘图数据
plot_data_d <- data.frame(
  time_point = rep(d_results$performance$time_point, 4),
  time_point_num = rep(d_results$performance$time_point_num, 4),
  predictor = c(
    rep("All Predictors", nrow(d_results$performance)),
    rep("HR", nrow(d_results$performance)),
    rep("BO", nrow(d_results$performance)),
    rep("Total Steps", nrow(d_results$performance))
  ),
  r2 = c(
    d_results$performance$all_r2,
    d_results$performance$hr_only_r2,
    d_results$performance$bo_only_r2,
    d_results$performance$steps_only_r2
  )
)

# 为NoD组创建绘图数据
plot_data_nod <- data.frame(
  time_point = rep(nod_results$performance$time_point, 4),
  time_point_num = rep(nod_results$performance$time_point_num, 4),
  predictor = c(
    rep("All Predictors", nrow(nod_results$performance)),
    rep("HR", nrow(nod_results$performance)),
    rep("BO", nrow(nod_results$performance)),
    rep("Total Steps", nrow(nod_results$performance))
  ),
  r2 = c(
    nod_results$performance$all_r2,
    nod_results$performance$hr_only_r2,
    nod_results$performance$bo_only_r2,
    nod_results$performance$steps_only_r2
  )
)

# 创建D组性能线图
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
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )
p_d
# 保存D组性能图
ggsave(file.path(output_dir, "predictor_performance_d.pdf"), p_d, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_d.png"), p_d, width = 10, height = 6, dpi = 300)

# 创建NoD组性能线图
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
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )
p_nod
# 保存NoD组性能图
ggsave(file.path(output_dir, "predictor_performance_nod.pdf"), p_nod, width = 10, height = 6)
ggsave(file.path(output_dir, "predictor_performance_nod.png"), p_nod, width = 10, height = 6, dpi = 300)

cat("已保存D组和NoD组的性能线图\n")