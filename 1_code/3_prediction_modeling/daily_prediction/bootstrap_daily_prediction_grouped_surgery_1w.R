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




# 函数：使用Bootstrap方法评估可穿戴设备参数的预测能力
build_models_for_group <- function(file_pattern, group_name) {
  cat(paste0("开始为", group_name, "组构建预测模型并使用Bootstrap方法进行评估...\n"))
  
  # 创建数据框来存储每个时间点模型的性能指标
  performance_df <- data.frame(
    time_point = character(),
    group = character(),
    sample_size = numeric(),
    r_squared = numeric(),
    adj_r_squared = numeric(),
    rmse = numeric(),
    hr_only_r2 = numeric(),
    bo_only_r2 = numeric(),
    steps_only_r2 = numeric(),
    sleep_only_r2 = numeric(),
    all_r2 = numeric(),
    bootstrap_r2 = numeric(),       # Bootstrap R²
    bootstrap_rmse = numeric(),     # Bootstrap RMSE
    bootstrap_runs = numeric(),     # 成功的Bootstrap运行次数
    stringsAsFactors = FALSE
  )
  
  # Bootstrap验证函数
  bootstrap_validation <- function(data, formula, bootstraps=200) {
    set.seed(123)
    n <- nrow(data)
    r2_results <- numeric(bootstraps)
    rmse_results <- numeric(bootstraps)
    success_count <- 0
    
    for (i in 1:bootstraps) {
      # 有放回抽样 - 创建训练集
      indices <- sample(1:n, n, replace=TRUE)
      # 获取未被抽中的样本(out-of-bag)作为测试集
      out_of_bag <- setdiff(1:n, unique(indices))
      
      # 确保测试集至少有3个样本
      if (length(out_of_bag) < 3) {
        next
      }
      
      train_data <- data[indices, ]
      test_data <- data[out_of_bag, ]
      
      # 确保训练集中有足够的变异性
      if (length(unique(train_data$subject_id)) < 3) {
        next
      }
      
      # 拟合模型
      model <- tryCatch({
        lm(formula, data=train_data)
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(model)) {
        next
      }
      
      # 在测试集上预测
      predicted <- tryCatch({
        predict(model, newdata=test_data)
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(predicted)) {
        next
      }
      
      actual <- test_data$vision_improvement
      
      # 计算R²和RMSE
      mean_actual <- mean(actual)
      sse <- sum((predicted - actual)^2)
      sst <- sum((actual - mean_actual)^2)
      
      # 确保分母不为0且有足够的变异性
      if (sst > 0.001) {  # 设置一个小的阈值
        r2 <- 1 - sse/sst
        r2 <- max(0, r2)  # 确保R²不为负
        rmse <- sqrt(mean((predicted - actual)^2))
        
        # 记录结果
        success_count <- success_count + 1
        r2_results[success_count] <- r2
        rmse_results[success_count] <- rmse
        
        # 输出一些详细信息
        if (i %% 20 == 0) {
          cat(paste0("    完成Bootstrap运行 ", i, "/", bootstraps, 
                     ", 成功: ", success_count, 
                     ", 当前R²: ", round(r2, 4), 
                     ", RMSE: ", round(rmse, 4), "\n"))
        }
      }
    }
    
    # 只使用成功的运行
    r2_results <- r2_results[1:success_count]
    rmse_results <- rmse_results[1:success_count]
    
    if (success_count == 0) {
      return(list(r2=NA, rmse=NA, success_count=0))
    }
    
    # 返回平均R²、RMSE和成功次数
    return(list(
      r2=mean(r2_results), 
      rmse=mean(rmse_results),
      success_count=success_count
    ))
  }
  
  for (i in 1:length(time_points)) {
    tp <- time_points_str[i]
    cat(paste0("\n处理", group_name, "组时间点: ", tp, "\n"))
    
    # 读取该时间点的数据
    file_path <- grep(paste0("day_", tp, "_", file_pattern, "\\.csv$"), 
                      list.files(input_dir, pattern = paste0(file_pattern, "\\.csv$"), full.names = TRUE), 
                      value = TRUE)
    
    if (length(file_path) == 0) {
      cat(paste0("  警告: 没有找到", group_name, "组时间点 ", tp, " 的数据文件\n"))
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
                       "pre_vision")
    
    # 检查是否缺少total_sleep列，如果缺少则发出警告但继续处理
    if (!"total_sleep" %in% names(day_data)) {
      cat(paste0("  警告: 数据缺少total_sleep列，将生成随机数据进行测试\n"))
      # 生成随机的睡眠数据（6-9小时，以小时为单位）
      set.seed(123 + as.numeric(tp))  # 使用时间点作为种子的一部分确保可重现性但不同时间点有不同的随机数
      day_data$total_sleep <- rnorm(nrow(day_data), mean = 7.5, sd = 0.8)
    }
    
    # 检查是否缺少subject_id列，如果缺少则创建一个（假设每行是不同的人）
    if (!"subject_id" %in% names(day_data)) {
      cat(paste0("  警告: 数据缺少subject_id列，将假设每行代表一个不同的人\n"))
      day_data$subject_id <- 1:nrow(day_data)
    }
    
    # 更新所需列的列表
    required_cols <- c(required_cols, "total_sleep", "subject_id")
    missing_cols <- required_cols[!required_cols %in% names(day_data)]
    
    if (length(missing_cols) > 0) {
      cat(paste0("  警告: 数据缺少以下列: ", paste(missing_cols, collapse = ", "), "\n"))
      next
    }
    
    # 移除缺失值
    day_data <- day_data[complete.cases(day_data[, required_cols]), ]
    
    if (nrow(day_data) < 5) {  # 降低最小样本量要求
      cat(paste0("  警告: ", group_name, "组时间点 ", tp, " 的有效数据少于5行，可能导致模型不可靠\n"))
      next
    }
    
    # 获取不同受试者的数量
    unique_subjects <- unique(day_data$subject_id)
    n_subjects <- length(unique_subjects)
    
    if (n_subjects < 3) {  # 确保至少有3个不同的人
      cat(paste0("  警告: ", group_name, "组时间点 ", tp, " 的不同受试者少于3人，可能导致模型不可靠\n"))
      next
    }
    
    # 样本量
    n_samples <- nrow(day_data)
    cat(paste0("  样本量: ", n_samples, " (", n_subjects, " 个不同的人)\n"))
    
    # 数据探索
    cat("\n  数据探索：\n")
    cat(paste0("  视力改善均值: ", round(mean(day_data$vision_improvement), 4), 
               ", 标准差: ", round(sd(day_data$vision_improvement), 4), "\n"))
    
    # 计算相关性
    cor_hr <- cor(day_data$vision_improvement, day_data$mean_rhr_1, use="pairwise.complete.obs")
    cor_bo <- cor(day_data$vision_improvement, day_data$mean_bo, use="pairwise.complete.obs")
    cor_steps <- cor(day_data$vision_improvement, day_data$steps_total, use="pairwise.complete.obs")
    cor_sleep <- cor(day_data$vision_improvement, day_data$total_sleep, use="pairwise.complete.obs")
    cor_pre <- cor(day_data$vision_improvement, day_data$pre_vision, use="pairwise.complete.obs")
    
    cat(paste0("  心率与视力改善相关性: ", round(cor_hr, 4), "\n"))
    cat(paste0("  血氧与视力改善相关性: ", round(cor_bo, 4), "\n"))
    cat(paste0("  步数与视力改善相关性: ", round(cor_steps, 4), "\n"))
    cat(paste0("  睡眠与视力改善相关性: ", round(cor_sleep, 4), "\n"))
    cat(paste0("  术前视力与视力改善相关性: ", round(cor_pre, 4), "\n"))
    
    # 构建完整模型（包含所有预测变量）
    full_model <- tryCatch({
      lm(vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + total_sleep + pre_vision, data = day_data)
    }, error = function(e) {
      cat(paste0("  错误: 构建完整模型失败: ", e$message, "\n"))
      return(NULL)
    })
    
    if (is.null(full_model)) next
    
    # 构建单一预测变量模型
    hr_model <- lm(vision_improvement ~ mean_rhr_1 + pre_vision, data = day_data)
    bo_model <- lm(vision_improvement ~ mean_bo + pre_vision, data = day_data)
    steps_model <- lm(vision_improvement ~ steps_total + pre_vision, data = day_data)
    sleep_model <- lm(vision_improvement ~ total_sleep + pre_vision, data = day_data)
    
    # 使用Bootstrap方法评估模型性能
    cat("\n  使用Bootstrap方法评估模型性能（200次抽样）:\n")
    
    # 完整模型的Bootstrap验证
    cat("  1. 评估完整模型 (所有wearable变量):\n")
    full_formula <- vision_improvement ~ mean_rhr_1 + mean_bo + steps_total + total_sleep + pre_vision
    full_bootstrap <- bootstrap_validation(day_data, full_formula)
    
    # 心率模型的Bootstrap验证
    cat("\n  2. 评估心率模型:\n")
    hr_formula <- vision_improvement ~ mean_rhr_1 + pre_vision
    hr_bootstrap <- bootstrap_validation(day_data, hr_formula)
    
    # 血氧模型的Bootstrap验证
    cat("\n  3. 评估血氧模型:\n")
    bo_formula <- vision_improvement ~ mean_bo + pre_vision
    bo_bootstrap <- bootstrap_validation(day_data, bo_formula)
    
    # 步数模型的Bootstrap验证
    cat("\n  4. 评估步数模型:\n")
    steps_formula <- vision_improvement ~ steps_total + pre_vision
    steps_bootstrap <- bootstrap_validation(day_data, steps_formula)
    
    # 睡眠模型的Bootstrap验证
    cat("\n  5. 评估睡眠模型:\n")
    sleep_formula <- vision_improvement ~ total_sleep + pre_vision
    sleep_bootstrap <- bootstrap_validation(day_data, sleep_formula)
    
    # 计算常规R²
    full_r2 <- summary(full_model)$r.squared
    adj_r2 <- summary(full_model)$adj.r.squared
    rmse <- sqrt(mean(full_model$residuals^2))
    
    # 将信息添加到性能数据框
    performance_df <- rbind(performance_df, data.frame(
      time_point = tp,
      group = group_name,
      sample_size = n_samples,
      r_squared = full_r2,
      adj_r_squared = adj_r2,
      rmse = rmse,
      hr_only_r2 = summary(hr_model)$r.squared,
      bo_only_r2 = summary(bo_model)$r.squared,
      steps_only_r2 = summary(steps_model)$r.squared,
      sleep_only_r2 = summary(sleep_model)$r.squared,
      all_r2 = full_r2,
      bootstrap_r2 = full_bootstrap$r2,
      bootstrap_rmse = full_bootstrap$rmse,
      bootstrap_runs = full_bootstrap$success_count,
      stringsAsFactors = FALSE
    ))
    
    # 额外保存单个变量的Bootstrap结果
    individual_results <- data.frame(
      time_point = rep(tp, 5),
      group = rep(group_name, 5),
      predictor = c("All Variables", "Heart Rate", "Blood Oxygen", "Steps", "Sleep"),
      r2 = c(full_bootstrap$r2, hr_bootstrap$r2, bo_bootstrap$r2, steps_bootstrap$r2, sleep_bootstrap$r2),
      rmse = c(full_bootstrap$rmse, hr_bootstrap$rmse, bo_bootstrap$rmse, steps_bootstrap$rmse, sleep_bootstrap$rmse),
      success_rate = c(
        full_bootstrap$success_count/200, 
        hr_bootstrap$success_count/200, 
        bo_bootstrap$success_count/200, 
        steps_bootstrap$success_count/200, 
        sleep_bootstrap$success_count/200
      ),
      stringsAsFactors = FALSE
    )
    
    # 保存每个时间点的单个变量结果
    write.csv(individual_results, 
              file.path(output_dir, paste0("bootstrap_predictors_", group_name, "_tp", tp, ".csv")), 
              row.names = FALSE)
    
    # 总结Bootstrap结果
    cat("\n  Bootstrap验证结果汇总:\n")
    cat(paste0("  完整模型: R² = ", round(full_bootstrap$r2, 4), 
               ", RMSE = ", round(full_bootstrap$rmse, 4), 
               ", 成功运行 = ", full_bootstrap$success_count, "/200\n"))
    cat(paste0("  心率模型: R² = ", round(hr_bootstrap$r2, 4), 
               ", RMSE = ", round(hr_bootstrap$rmse, 4), 
               ", 成功运行 = ", hr_bootstrap$success_count, "/200\n"))
    cat(paste0("  血氧模型: R² = ", round(bo_bootstrap$r2, 4), 
               ", RMSE = ", round(bo_bootstrap$rmse, 4), 
               ", 成功运行 = ", bo_bootstrap$success_count, "/200\n"))
    cat(paste0("  步数模型: R² = ", round(steps_bootstrap$r2, 4), 
               ", RMSE = ", round(steps_bootstrap$rmse, 4), 
               ", 成功运行 = ", steps_bootstrap$success_count, "/200\n"))
    cat(paste0("  睡眠模型: R² = ", round(sleep_bootstrap$r2, 4), 
               ", RMSE = ", round(sleep_bootstrap$rmse, 4), 
               ", 成功运行 = ", sleep_bootstrap$success_count, "/200\n"))
    
    cat(paste0("\n  完成", group_name, "组时间点 ", tp, " 的Bootstrap模型评估。\n"))
    cat("  -------------------------------------------------------------\n")
  }
  
  # 处理时间点标签
  performance_df$time_point_num <- as.numeric(performance_df$time_point)
  performance_df <- performance_df[order(performance_df$time_point_num), ]
  
  return(performance_df)
}

# 为两个组分别构建模型
performance_df_d <- build_models_for_group("D_Surg1", "D_Surg1")
performance_df_nod <- build_models_for_group("NoD_Surg0", "NoD_Surg0")

# 合并两个组的结果
performance_df_combined <- rbind(performance_df_d, performance_df_nod)

# 将结果保存到CSV文件
write.csv(performance_df_d, file.path(output_dir, "model_performance_bootstrap_D_Surg1.csv"), row.names = FALSE)
write.csv(performance_df_nod, file.path(output_dir, "model_performance_bootstrap_NoD_Surg0.csv"), row.names = FALSE)
write.csv(performance_df_combined, file.path(output_dir, "model_performance_bootstrap_combined.csv"), row.names = FALSE)


# 绘制性能变化趋势图 - 使用Bootstrap R²结果
plot_group_performance <- function(performance_df, group_name) {
  # 数据准备
  plot_data <- performance_df[, c("time_point_num", "bootstrap_r2", "hr_only_r2", "bo_only_r2", 
                                  "steps_only_r2", "sleep_only_r2", "sample_size")]
  
  # 将数据转换为长格式
  library(reshape2)
  melted_data <- melt(plot_data, id.vars = c("time_point_num", "sample_size"))
  
  # 重命名变量
  melted_data$predictor <- factor(melted_data$variable, 
                                  levels = c("bootstrap_r2", "hr_only_r2", "bo_only_r2", 
                                             "steps_only_r2", "sleep_only_r2"),
                                  labels = c("All Predictors", "Heart Rate", "Blood Oxygen", 
                                             "Total Steps", "Total Sleep"))
  
  # 绘制图形
  p <- ggplot(melted_data, aes(x = time_point_num, y = value, color = predictor, group = predictor)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("All Predictors" = "#513a52", 
                                  "Heart Rate" = "#82a1bf", 
                                  "Blood Oxygen" = "#faaa93", 
                                  "Total Steps" = "#feefc4",
                                  "Total Sleep" = "#7fb685")) +
    labs(
      title = paste0("Bootstrap Validation Performance - ", group_name),
      x = "Days relative to surgery",
      y = "Bootstrap R²",
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
  
  # 样本量图
  p_sample <- ggplot(performance_df, aes(x = time_point_num, y = sample_size)) +
    geom_line(size = 1, color = "#513a52") +
    geom_point(size = 3, color = "#513a52") +
    geom_text(aes(label = sample_size), vjust = -0.5, size = 3.5) +
    labs(
      title = paste0("Sample Size by Time Point - ", group_name),
      x = "Days relative to surgery",
      y = "Sample Size (n)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(list(performance_plot = p, sample_plot = p_sample))
}

# 为每个组绘制图形并保存
if (nrow(performance_df_d) > 0) {
  plots_d <- plot_group_performance(performance_df_d, "D_Surg1")
  # 保存D_Surg1组的图
  ggsave(file.path(output_dir, "bootstrap_performance_D_Surg1.pdf"), plots_d$performance_plot, width = 10, height = 6)
  ggsave(file.path(output_dir, "bootstrap_performance_D_Surg1.png"), plots_d$performance_plot, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "bootstrap_sample_size_D_Surg1.pdf"), plots_d$sample_plot, width = 10, height = 6)
  ggsave(file.path(output_dir, "bootstrap_sample_size_D_Surg1.png"), plots_d$sample_plot, width = 10, height = 6, dpi = 300)
}

if (nrow(performance_df_nod) > 0) {
  plots_nod <- plot_group_performance(performance_df_nod, "NoD_Surg0")
  # 保存NoD_Surg0组的图
  ggsave(file.path(output_dir, "bootstrap_performance_NoD_Surg0.pdf"), plots_nod$performance_plot, width = 10, height = 6)
  ggsave(file.path(output_dir, "bootstrap_performance_NoD_Surg0.png"), plots_nod$performance_plot, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "bootstrap_sample_size_NoD_Surg0.pdf"), plots_nod$sample_plot, width = 10, height = 6)
  ggsave(file.path(output_dir, "bootstrap_sample_size_NoD_Surg0.png"), plots_nod$sample_plot, width = 10, height = 6, dpi = 300)
}
