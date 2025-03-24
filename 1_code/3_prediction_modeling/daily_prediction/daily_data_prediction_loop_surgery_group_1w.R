# Comprehensive Analysis for Vision Improvement Prediction by Group
# 修正版：适用于分组分文件的小样本预测分析
# 增强版：添加Random Forest, XGBoost和Elastic Net模型

library(glmnet)
library(caret)
library(tidyverse)
library(mice)
library(VIM)
library(randomForest)  # 添加Random Forest包
library(xgboost)       # 添加XGBoost包
setwd(get_project_wd())
rm(list = ls())

# 定义数据目录 - 修改为直接指向包含CSV文件的目录
data_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped/"  # 假设CSV文件与脚本在同一目录
output_dir <- "3_data_analysis/3_prediction_modeling/1w_prediction/daily_data_grouped/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 定义计算性能指标的函数
calculate_performance <- function(model) {
  # 首先检查model$pred是否存在并且是一个数据框
  if(!is.null(model$pred) && is.data.frame(model$pred)) {
    # 检查是否存在需要的列
    if(all(c("obs", "pred") %in% names(model$pred))) {
      # 直接在整个数据集上计算指标
      RMSE <- sqrt(mean((model$pred$obs - model$pred$pred)^2))
      MAE <- mean(abs(model$pred$obs - model$pred$pred))
      R2 <- cor(model$pred$obs, model$pred$pred)^2
      
      # 返回整体指标
      result <- data.frame(
        RMSE_mean = RMSE,
        RMSE_sd = NA_real_,
        MAE_mean = MAE,
        MAE_sd = NA_real_,
        R2_mean = R2,
        R2_sd = NA_real_
      )
      return(result)
    }
  }
  
  # 如果无法计算性能，返回NA值
  warning("无法计算模型性能，返回NA值")
  result <- data.frame(
    RMSE_mean = NA_real_,
    RMSE_sd = NA_real_,
    MAE_mean = NA_real_,
    MAE_sd = NA_real_,
    R2_mean = NA_real_,
    R2_sd = NA_real_
  )
  return(result)
}

# 分析单个文件数据的函数
analyze_file_data <- function(file_path) {
  # 提取文件名获取日期和组别信息
  file_name <- basename(file_path)
  parts <- strsplit(file_name, "\\.")[[1]][1]  # 去掉扩展名
  
  # 解析文件名获取天数和组别
  # 假设格式为 day_{天数}_{组别}.csv
  name_parts <- strsplit(parts, "_")[[1]]
  day_str <- name_parts[2]
  group_str <- paste(name_parts[3:length(name_parts)], collapse="_")
  
  # 设置日期标签
  day_label <- paste("Day", day_str)
  group_label <- ifelse(grepl("NoD", group_str), "No Diabetes Surgery 0", "Diabetes Surgery 1")
  
  cat("\n分析文件:", file_name, "\n")
  cat("天数:", day_label, "组别:", group_label, "\n")
  
  # 加载数据
  day_data <- tryCatch({
    read_csv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    cat("无法加载文件:", e$message, "\n")
    return(NULL)
  })
  
  # 检查是否成功加载
  if(is.null(day_data) || nrow(day_data) == 0) {
    cat("文件为空或加载失败\n")
    return(NULL)
  }
  
  cat("数据集大小:", nrow(day_data), "行\n")
  
  # 检查目标变量
  if(!"vision_improvement" %in% names(day_data)) {
    cat("缺少目标变量vision_improvement\n")
    return(NULL)
  }
  
  # 去除缺失目标变量的行
  day_data <- day_data[!is.na(day_data$vision_improvement), ]
  cat("去除缺失目标变量后的样本量:", nrow(day_data), "\n")
  
  # 如果样本量太小，返回NULL
  if(nrow(day_data) < 5) {
    cat("样本量太小（<5），跳过分析\n")
    return(NULL)
  }
  
  # 为分类变量转换为因子
  if("gender" %in% names(day_data)) {
    day_data$gender <- as.factor(day_data$gender)
  }
  
  # 小样本情况下选择关键特征
  key_features <- intersect(c(
    # RHR features - 仅选择均值和标准差
    "mean_rhr_1", "sd_rhr_1",
    # BO features - 仅选择均值和标准差
    "mean_bo", "sd_bo",
    # Steps features - 仅选择总步数
    "steps_total","steps_max",
    # sleep feature - 仅选择总睡眠
    "total_sleep","deep_sleep",
    # 关键人口统计学和医疗史特征
    "age", "gender", "pre_vision","bmi","season"
  ), names(day_data))
  
  # 打印可用特征
  cat("可用特征:", paste(key_features, collapse=", "), "\n")
  
  # 检查是否有足够的特征
  if(length(key_features) < 2) {
    cat("可用特征不足，跳过分析\n")
    return(NULL)
  }
  
  # 准备模型数据 - 只包含关键特征和目标变量
  model_data <- day_data[, c(key_features, "vision_improvement")]
  
  # 去除有缺失值的行
  model_data <- model_data[complete.cases(model_data), ]
  cat("去除缺失值后的样本量:", nrow(model_data), "\n")
  
  # 如果样本量太小，返回NULL
  if(nrow(model_data) < 5) {
    cat("完整案例样本量不足，跳过分析\n")
    return(NULL)
  }
  
  # 设置留一法交叉验证
  cat("使用留一法交叉验证\n")
  ctrl <- trainControl(
    method = "LOOCV",
    savePredictions = "final",
    summaryFunction = defaultSummary
  )
  
  # 训练模型
  set.seed(123)
  
  results_list <- list()
  
  # 构建公式对象
  formula_obj <- as.formula(paste("vision_improvement ~", paste(key_features, collapse = " + ")))
  
  # 线性回归
  lm_result <- tryCatch({
    model <- train(
      formula_obj,
      data = model_data,
      method = "lm",
      trControl = ctrl,
      metric = "RMSE"
    )
    perf <- calculate_performance(model)
    perf$Model <- "Linear Regression"
    perf
  }, error = function(e) {
    cat("线性回归模型训练失败:", e$message, "\n")
    NULL
  })
  
  if(!is.null(lm_result)) {
    results_list$lm <- lm_result
    cat("线性回归模型训练成功, R² =", round(lm_result$R2_mean, 3), "\n")
  }
  
  # LASSO
  lasso_result <- tryCatch({
    lasso_grid <- expand.grid(
      alpha = 1,
      lambda = seq(0.001, 0.1, length.out = 5)
    )
    
    model <- train(
      formula_obj,
      data = model_data,
      method = "glmnet",
      tuneGrid = lasso_grid,
      trControl = ctrl,
      metric = "RMSE"
    )
    perf <- calculate_performance(model)
    perf$Model <- "LASSO"
    perf
  }, error = function(e) {
    cat("LASSO模型训练失败:", e$message, "\n")
    NULL
  })
  
  if(!is.null(lasso_result)) {
    results_list$lasso <- lasso_result
    cat("LASSO模型训练成功, R² =", round(lasso_result$R2_mean, 3), "\n")
  }
  
  # 添加Elastic Net模型
  enet_result <- tryCatch({
    enet_grid <- expand.grid(
      alpha = seq(0.2, 0.8, by = 0.2),  # 混合参数
      lambda = seq(0.001, 0.1, length.out = 3)  # 正则化参数
    )
    
    model <- train(
      formula_obj,
      data = model_data,
      method = "glmnet",
      tuneGrid = enet_grid,
      trControl = ctrl,
      metric = "RMSE"
    )
    perf <- calculate_performance(model)
    perf$Model <- "Elastic Net"
    perf
  }, error = function(e) {
    cat("Elastic Net模型训练失败:", e$message, "\n")
    NULL
  })
  
  if(!is.null(enet_result)) {
    results_list$enet <- enet_result
    cat("Elastic Net模型训练成功, R² =", round(enet_result$R2_mean, 3), "\n")
  }
  
  # 添加Random Forest模型 - 适合小样本
  rf_result <- tryCatch({
    # 对于小样本，使用较小的mtry和ntree值
    rf_grid <- expand.grid(
      mtry = c(2, max(2, floor(sqrt(length(key_features)))))
    )
    
    model <- train(
      formula_obj,
      data = model_data,
      method = "rf",
      tuneGrid = rf_grid,
      trControl = ctrl,
      metric = "RMSE",
      ntree = 100  # 减少树的数量避免过拟合
    )
    perf <- calculate_performance(model)
    perf$Model <- "Random Forest"
    perf
  }, error = function(e) {
    cat("Random Forest模型训练失败:", e$message, "\n")
    NULL
  })
  
  if(!is.null(rf_result)) {
    results_list$rf <- rf_result
    cat("Random Forest模型训练成功, R² =", round(rf_result$R2_mean, 3), "\n")
  }
  
  # 添加XGBoost模型
  xgb_result <- tryCatch({
    # 对于小样本的XGBoost参数
    xgb_grid <- expand.grid(
      nrounds = c(20, 50),            # 减少迭代次数避免过拟合
      max_depth = c(2, 3),            # 较小的树深度
      eta = c(0.1, 0.3),              # 学习率
      gamma = 0,                      # 最小损失减少
      colsample_bytree = 0.8,         # 特征采样比例
      min_child_weight = 1,           # 最小子节点权重
      subsample = 0.8                 # 样本采样比例
    )
    
    model <- train(
      formula_obj,
      data = model_data,
      method = "xgbTree",
      tuneGrid = xgb_grid,
      trControl = ctrl,
      metric = "RMSE"
    )
    perf <- calculate_performance(model)
    perf$Model <- "XGBoost"
    perf
  }, error = function(e) {
    cat("XGBoost模型训练失败:", e$message, "\n")
    NULL
  })
  
  if(!is.null(xgb_result)) {
    results_list$xgb <- xgb_result
    cat("XGBoost模型训练成功, R² =", round(xgb_result$R2_mean, 3), "\n")
  }
  
  # 过滤掉NULL结果
  results_list <- results_list[!sapply(results_list, is.null)]
  
  # 检查是否有任何结果
  if(length(results_list) == 0) {
    cat("没有模型成功训练\n")
    return(NULL)
  }
  
  # 合并所有结果
  combined_results <- bind_rows(results_list)
  
  # 添加日期和组别信息
  if(nrow(combined_results) > 0) {
    combined_results$Day <- day_label
    combined_results$Group <- group_label
    combined_results$SampleSize <- nrow(model_data)
    
    cat("分析完成，成功训练了", nrow(combined_results), "个模型\n")
    return(combined_results)
  } else {
    cat("没有产生任何有效结果\n")
    return(NULL)
  }
}

# 寻找所有CSV文件
csv_files <- list.files(data_dir, pattern = "day_.+\\.csv$", full.names = TRUE)

if(length(csv_files) == 0) {
  stop("未找到任何符合格式的CSV文件")
}

cat("找到", length(csv_files), "个CSV文件\n")

# 分析所有文件
all_results <- data.frame()

for(file_path in csv_files) {
  file_results <- analyze_file_data(file_path)
  
  if(!is.null(file_results) && nrow(file_results) > 0) {
    all_results <- bind_rows(all_results, file_results)
  }
}



# First, define all the plotting helper functions that are needed
# 辅助函数：保存图形为PDF和PNG
save_plots <- function(plot, output_dir, base_name, width = 10, height = 7) {
  ggsave(
    file.path(output_dir, paste0(base_name, ".pdf")),
    plot, width = width, height = height
  )
  
  ggsave(
    file.path(output_dir, paste0(base_name, ".png")),
    plot, width = width, height = height, dpi = 100
  )
}

# 辅助函数：创建组间比较热图
create_heatmap <- function(all_results, model_name, output_dir) {
  # 提取数据
  model_data <- all_results %>% 
    filter(Model == model_name) %>%
    dplyr::select(Day, Group, R2_mean, SampleSize)
  
  if(nrow(model_data) == 0) return(NULL)
  
  # 创建热图
  heatmap_plot <- ggplot(model_data, aes(x = Day, y = Group, fill = R2_mean)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f\nn=%d", R2_mean, SampleSize)), color = "white", size = 3) +
    scale_fill_gradient(low = "pink", high = "darkred", limits = c(0, 1), na.value = "grey90") +
    labs(
      title = paste(model_name, "- R² Comparison Between Groups"),
      x = "Days Relative to Surgery",
      y = "Group",
      fill = "R²"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  # 保存图形
  save_plots(heatmap_plot, output_dir, paste0("heatmap_", gsub(" ", "_", model_name)), width = 12, height = 6)
}

# 辅助函数：创建模型比较热图
create_model_heatmap <- function(all_results, group_name, output_dir) {
  # 提取数据
  group_data <- all_results %>% 
    filter(Group == group_name) %>%
    dplyr::select(Day, Model, R2_mean, SampleSize)
  
  if(nrow(group_data) == 0) return(NULL)
  
  # 创建热图
  heatmap_plot <- ggplot(group_data, aes(x = Day, y = Model, fill = R2_mean)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f\nn=%d", R2_mean, SampleSize)), color = "white", size = 3) +
    scale_fill_gradient(low = "pink", high = "darkred", limits = c(0, 1), na.value = "grey90") +
    labs(
      title = paste(group_name, "- R² Across Days"),
      x = "Days Relative to Surgery",
      y = "Model",
      fill = "R²"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  # 保存图形
  save_plots(heatmap_plot, output_dir, paste0("day_heatmap_", gsub(" ", "_", group_name)), width = 12, height = 6)
}

# Now define the main plotting function
# 创建绘图函数
create_all_plots <- function(all_results, output_dir) {
  # 确保输出目录存在
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 添加样本量标签，避免重复代码
  all_results <- all_results %>%
    mutate(SizeLbl = paste0("n=", SampleSize))
  
  # 定义通用的ggplot主题
  my_theme <- theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # 获取组和模型列表，避免多次重复提取
  groups <- unique(all_results$Group)
  models <- unique(all_results$Model)
  
  # 1. 为每个组创建线图 (R²随天数变化)
  message("创建组别模型表现线图...")
  
  # 使用map函数代替for循环
  purrr::walk(groups, function(group_name) {
    group_data <- all_results %>% filter(Group == group_name)
    
    if(nrow(group_data) == 0) return(NULL)
    
    # 线图
    line_plot <- ggplot(group_data, aes(x = Day, y = R2_mean, color = Model, group = Model)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      geom_text(aes(label = SizeLbl), vjust = -1.5, size = 2.5, show.legend = FALSE) +
      labs(
        title = paste("Vision Improvement Prediction -", group_name),
        subtitle = "Comparing Different Machine Learning Models with LOOCV",
        x = "Days Relative to Surgery",
        y = "R² (coefficient of determination)"
      ) +
      ylim(0, min(1, max(group_data$R2_mean) * 1.2)) +
      my_theme
    
    # 保存图形（PDF和PNG）
    save_plots(line_plot, output_dir, paste0("prediction_performance_", gsub(" ", "_", group_name)))
  })
  
  # 2. 只有当有多个组时才创建组间比较图
  if(length(groups) > 1) {
    message("创建模型组间比较图...")
    
    # 使用map函数代替for循环
    purrr::walk(models, function(model_name) {
      model_data <- all_results %>% filter(Model == model_name)
      
      if(nrow(model_data) == 0) return(NULL)
      
      # 线图比较
      comparison_plot <- ggplot(model_data, aes(x = Day, y = R2_mean, color = Group, group = Group)) +
        geom_line(size = 1) +
        geom_point(size = 3) +
        geom_text(aes(label = SizeLbl), vjust = -1.5, size = 2.5, show.legend = FALSE) +
        labs(
          title = paste(model_name, "- Group Comparison"),
          subtitle = "Comparing prediction performance between groups using LOOCV",
          x = "Days Relative to Surgery",
          y = "R²"
        ) +
        ylim(0, min(1, max(model_data$R2_mean) * 1.2)) +
        my_theme
      
      # 保存图形
      save_plots(comparison_plot, output_dir, paste0("group_comparison_", gsub(" ", "_", model_name)))
      
      # 热图比较
      create_heatmap(all_results, model_name, output_dir)
    })
  }
  
  # 3. 为每个组创建模型比较热图
  message("创建模型热图...")
  purrr::walk(groups, function(group_name) {
    create_model_heatmap(all_results, group_name, output_dir)
  })
  
  message("所有图表已创建并保存到", output_dir)
}

# After defining all functions, now we can use them

# Your main code block should be like this:
if(nrow(all_results) == 0) {
  cat("\n\n===== 警告 =====\n")
  cat("没有产生任何有效结果。这可能是由于样本量太小或数据质量问题。\n")
  
  # 创建一个空的结果文件，表明代码已运行
  write.csv(data.frame(Note="No valid results produced"), 
            file.path(output_dir, "no_results_produced.csv"), row.names = FALSE)
} else {
  cat("\n\n分析完成，共生成", nrow(all_results), "条结果记录\n")
  
  # 保存所有结果
  write.csv(all_results, file.path(output_dir, "all_model_results.csv"), row.names = FALSE)
  
  # 定义日期顺序
  day_order <- c(paste("Day", -7:-1), paste("Day", 0:29))
  
  # 确保Day是有序因子
  all_results$Day <- factor(all_results$Day, levels = day_order, ordered = TRUE)
  
  # 为每个组别创建最佳天数结果
  best_days_by_group <- all_results %>%
    group_by(Group, Model) %>%
    slice_max(order_by = R2_mean, n = 1) %>%
    dplyr::select(Group, Model, Day, R2_mean, RMSE_mean, SampleSize) %>%
    arrange(Group, desc(R2_mean))
  
  write.csv(best_days_by_group, file.path(output_dir, "best_days_by_group_model.csv"), row.names = FALSE)
  
  # 为每个组别创建总体模型表现总结
  model_performance_by_group <- all_results %>%
    group_by(Group, Model) %>%
    summarise(
      Avg_R2 = mean(R2_mean, na.rm=TRUE),
      Max_R2 = max(R2_mean, na.rm=TRUE),
      Min_R2 = min(R2_mean, na.rm=TRUE),
      Avg_RMSE = mean(RMSE_mean, na.rm=TRUE),
      Avg_SampleSize = mean(SampleSize),
      .groups = "drop"
    ) %>%
    arrange(Group, desc(Avg_R2))
  
  write.csv(model_performance_by_group, file.path(output_dir, "model_performance_by_group.csv"), row.names = FALSE)
  
  # Now we can call this function because it's already defined above
  create_all_plots(all_results, output_dir)
}

