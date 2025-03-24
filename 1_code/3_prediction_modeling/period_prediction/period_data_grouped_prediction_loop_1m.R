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

# 定义时间窗口 - 添加自第一个脚本
time_windows <- c(
  "pre_3d",
  "pre_3d_7d", 
  "pre_7d_all", 
  "post_7d",
  "post_7d_30d",
  "post_day23_30",
  "post_day27_30"
)

# 定义数据目录 - 修改为直接指向包含CSV文件的目录
data_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data_grouped/"  
output_dir <- "3_data_analysis/3_prediction_modeling/1m_prediction/time_period_data_grouped/results"
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

# 改进的文件名解析函数
parse_filename <- function(file_name) {
  # 移除扩展名
  base_name <- strsplit(file_name, "\\.")[[1]][1]
  
  # 识别组别
  if(grepl("D_Surg1", base_name)) {
    group_label <- "Diabetes Surgery 1"
  } else if(grepl("NoD_Surg0", base_name)) {
    group_label <- "No Diabetes Surgery 0"
  } else {
    return(list(success = FALSE, message = paste("无法识别组别:", file_name)))
  }
  
  # 提取时间信息
  # 首先判断是pre还是post
  if(grepl("^pre_", base_name)) {
    prefix <- "pre"
  } else if(grepl("^post_", base_name)) {
    prefix <- "post"
  } else {
    return(list(success = FALSE, message = paste("无法识别时间前缀:", file_name)))
  }
  
  # 移除前缀和组别，保留中间的时间信息
  time_info <- gsub(paste0("^", prefix, "_"), "", base_name)
  time_info <- gsub("_(NoD_Surg0|D_Surg1)$", "", time_info)
  
  # 构建完整的时间段标签
  day_label <- paste(prefix, time_info, sep="_")
  
  # 返回结果
  return(list(
    success = TRUE, 
    day_label = day_label, 
    group_label = group_label,
    prefix = prefix,
    time_info = time_info
  ))
}

# 改进的特征前缀识别函数，更灵活地处理各种时间格式
get_modified_feature_prefix <- function(day_label) {
  # 先尝试使用原始的映射表
  window_prefix_map <- list(
    "pre_3d" = "pre_surgery_3d",
    "pre_3d_7d" = "pre_surgery_3d_to_7d", 
    "pre_7d_all" = "pre_surgery_7d_all", 
    "post_4to6d" = "post_surgery_4to6d",
    "post_7d" = "post_surgery_7d",
    "post_7d_30d" = "post_surgery_7d_to_30d",
    "post_day23_30" = "post_surgery_day23_to_30",
    "post_day27_30" = "post_surgery_day27_to_30"
  )
  
  # 直接查找是否有精确匹配
  if(day_label %in% names(window_prefix_map)) {
    return(window_prefix_map[[day_label]])
  }
  
  # 如果没有精确匹配，基于更通用的模式识别
  if(grepl("^pre_", day_label)) {
    # 处理前缀为pre的各种情况
    if(grepl("3d_7d|3d_to_7d", day_label)) {
      return("pre_surgery_3d_to_7d")
    } else if(grepl("3d", day_label) && !grepl("_7d|_all", day_label)) {
      return("pre_surgery_3d")
    } else if(grepl("7d_all", day_label)) {
      return("pre_surgery_7d_all")
    } else if(grepl("all", day_label)) {
      return("pre_surgery_all")
    } else {
      # 默认前缀
      return("pre_surgery_general")
    }
  } else if(grepl("^post_", day_label)) {
    # 处理前缀为post的各种情况
    if(grepl("4to6d", day_label)) {
      return("post_surgery_4to6d")
    } else if(grepl("6d", day_label) && !grepl("4to6d", day_label)) {
      return("post_surgery_6d")
    } else if(grepl("1w", day_label)) {
      return("post_surgery_1w")
    } else if(grepl("1to3d", day_label)) {
      return("post_surgery_1to3d")
    } else if(grepl("7d_30d|7d_to_30d", day_label)) {
      return("post_surgery_7d_to_30d")
    } else if(grepl("7d", day_label) && !grepl("30d|_all", day_label)) {
      return("post_surgery_7d")
    } else if(grepl("day23_30|day27_30", day_label)) {
      return("post_surgery_30d")
    } else if(grepl("over_30d", day_label)) {
      return("post_surgery_over_30d")
    } else {
      # 默认前缀
      return("post_surgery_general")
    }
  } else {
    # 无法识别的情况
    warning("无法识别的时间段标签: ", day_label, "，使用通用前缀")
    return("general")
  }
}

# 分析单个文件数据的函数
analyze_file_data <- function(file_path) {
  # 提取文件名
  file_name <- basename(file_path)
  
  # 解析文件名
  file_info <- parse_filename(file_name)
  
  # 检查解析是否成功
  if(!file_info$success) {
    cat(file_info$message, "\n")
    return(NULL)
  }
  
  # 检查时间段是否在允许的范围内
  day_label <- file_info$day_label
  if(!(day_label %in% time_windows)) {
    cat("时间段", day_label, "不在分析范围内，跳过\n")
    return(NULL)
  }
  
  # 提取解析结果
  group_label <- file_info$group_label
  
  cat("\n分析文件:", file_name, "\n")
  cat("时间段:", day_label, "组别:", group_label, "\n")
  
  # 加载数据 (其余代码保持不变)
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
  if(!"vision_improvement_1m" %in% names(day_data)) {
    cat("缺少目标变量vision_improvement_1m\n")
    return(NULL)
  }
  
  # 去除缺失目标变量的行
  day_data <- day_data[!is.na(day_data$vision_improvement_1m), ]
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
  
  # 转换season为因子
  if("season" %in% names(day_data) && !is.factor(day_data$season)) {
    day_data$season <- as.factor(day_data$season)
  }
  
  # 使用改进的特征前缀识别函数
  prefix <- get_modified_feature_prefix(day_label)
  cat("使用特征前缀:", prefix, "\n")
  
  # 使用matches函数选择特征
  # 定义一个辅助函数，类似于matches但使用grep
  match_features <- function(pattern, names_vector) {
    grep(pattern, names_vector, value = TRUE)
  }
  
  # 选择关键特征，使用与第一个脚本相同的特征逻辑
  key_features <- c(
    # RHR features
    match_features(paste0("^", prefix, "_mean_rhr_steps_1$"), names(day_data)),
    match_features(paste0("^", prefix, "_min_rhr_steps_1$"), names(day_data)),
    match_features(paste0("^", prefix, "_max_rhr_steps_1$"), names(day_data)),
    # match_features(paste0("^", prefix, "_median_rhr_steps_1$"), names(day_data)),
    match_features(paste0("^", prefix, "_sd_rhr_steps_1$"), names(day_data)),
    # match_features(paste0("^", prefix, "_iqr_rhr_steps_1$"), names(day_data)),
    # BO features
    match_features(paste0("^", prefix, "_bo_mean$"), names(day_data)),
    # match_features(paste0("^", prefix, "_bo_min$"), names(day_data)),
    # match_features(paste0("^", prefix, "_bo_max$"), names(day_data)),
    # match_features(paste0("^", prefix, "_bo_median$"), names(day_data)),
    match_features(paste0("^", prefix, "_bo_sd$"), names(day_data)),
    # match_features(paste0("^", prefix, "_bo_iqr$"), names(day_data)),
    # Sleep features
    match_features(paste0("^", prefix, "_deep_sleep_total$"), names(day_data)),
    match_features(paste0("^", prefix, "_total_sleep_total$"), names(day_data)),
    # Steps features
    match_features(paste0("^", prefix, "_steps_total$"), names(day_data)),
    # match_features(paste0("^", prefix, "_steps_mean$"), names(day_data)),
    match_features(paste0("^", prefix, "_steps_max$"), names(day_data))
  )
  
  # 添加人口统计学和医疗史特征
  demographic_features <- c("age", "gender", "pre_vision", "season", "bmi","vision_improvement_1w")
  key_features <- c(key_features, intersect(demographic_features, names(day_data)))
  
  # 确保所有特征都在数据集中存在
  key_features <- intersect(key_features, names(day_data))
  
  # 打印可用特征
  cat("可用特征:", paste(key_features, collapse=", "), "\n")
  
  # 检查是否有足够的特征
  if(length(key_features) < 2) {
    cat("可用特征不足，跳过分析\n")
    return(NULL)
  }
  
  # 准备模型数据 - 只包含关键特征和目标变量
  model_data <- day_data[, c(key_features, "vision_improvement_1m")]
  
  # 检查缺失值
  missing_summary <- sapply(model_data, function(x) sum(is.na(x)))
  cat("缺失值数量:", paste(names(missing_summary)[missing_summary > 0], 
                      "=", missing_summary[missing_summary > 0], 
                      collapse=", "), "\n")
  
  # 如果有缺失值且样本量足够，执行缺失值插补
  if(sum(missing_summary) > 0 && nrow(model_data) >= 5) {
    # 设置随机种子
    set.seed(1234)
    
    # 配置插补方法
    imputation_methods <- make.method(model_data)
    
    # 创建5个插补数据集
    imp <- mice(model_data, m=5, method=imputation_methods, 
                maxit=50, seed=1234, printFlag=FALSE)
    
    # 使用第一个插补集创建完整数据集
    model_data <- complete(imp, 1)
    cat("缺失值插补完成\n")
  } else if(sum(missing_summary) > 0) {
    # 如果有缺失值但样本量不足，则去除有缺失值的行
    model_data <- model_data[complete.cases(model_data), ]
    cat("样本量太小无法插补，移除缺失值的行，剩余样本量:", nrow(model_data), "\n")
  }
  
  # 如果样本量太小，返回NULL
  if(nrow(model_data) < 5) {
    cat("完整案例样本量不足，跳过分析\n")
    return(NULL)
  }
  
  # 修复LASSO特征选择部分代码
  # 将这段代码替换到analyze_file_data函数中相应位置
  
  # 使用LASSO进行特征选择
  cat("使用LASSO进行特征选择...\n")
  
  # 检查是否包含分类变量
  categorical_vars <- names(model_data)[sapply(model_data, is.factor)]
  cat("检测到的分类变量:", paste(categorical_vars, collapse=", "), "\n")
  
  # 使用模型矩阵来处理分类变量，自动创建哑变量
  tryCatch({
    # 创建公式用于model.matrix
    formula_text <- paste("~", paste(key_features, collapse=" + "))
    formula_obj <- as.formula(formula_text)
    
    # 使用model.matrix创建特征矩阵，包括处理分类变量
    x <- model.matrix(formula_obj, data=model_data)[,-1]  # 去掉截距项
    y <- model_data$vision_improvement_1m
    
    cat("特征矩阵大小:", nrow(x), "行 x", ncol(x), "列\n")
    
    # 检查特征矩阵是否有效
    if(ncol(x) == 0) {
      cat("特征矩阵为空，跳过LASSO，使用所有原始特征\n")
      selected_features <- key_features
    } else {
      # 使用交叉验证确定最佳lambda
      set.seed(123)
      # 使用tryCatch确保即使LASSO失败也能继续
      cv_fit <- tryCatch({
        cv.glmnet(x, y, alpha=1, nfolds=min(5, nrow(model_data)))
      }, error = function(e) {
        cat("LASSO特征选择失败:", e$message, "\n")
        return(NULL)
      })
      
      # 如果LASSO成功运行，提取非零系数的特征
      if(!is.null(cv_fit)) {
        # 获取最佳lambda对应的系数
        coefs <- coef(cv_fit, s="lambda.min")
        
        # 获取特征名称（从model.matrix）
        feature_names <- rownames(coefs)[-1]  # 去掉截距项
        
        # 获取非零系数的索引
        nonzero_idx <- which(coefs[-1] != 0)  # 忽略截距项
        
        if(length(nonzero_idx) > 0) {
          # 获取LASSO选择的特征名称
          selected_feature_names <- feature_names[nonzero_idx]
          
          # 将LASSO选择的特征映射回原始变量名
          # 这一步很关键，因为model.matrix可能创建了多个哑变量
          selected_features <- c()
          for(feat in key_features) {
            # 检查原始特征名是否被保留（对于连续变量）
            if(feat %in% selected_feature_names) {
              selected_features <- c(selected_features, feat)
            } 
            # 对于分类变量，检查任何相关的哑变量是否被保留
            else if(feat %in% categorical_vars && 
                    any(grepl(paste0("^", feat), selected_feature_names))) {
              selected_features <- c(selected_features, feat)
            }
          }
          
          cat("LASSO选择了", length(selected_features), "个特征:", 
              paste(selected_features, collapse=", "), "\n")
        } else {
          cat("LASSO未选择任何特征，将使用所有特征\n")
          selected_features <- key_features
        }
      } else {
        cat("LASSO执行失败，将使用所有特征\n")
        selected_features <- key_features
      }
    }
  }, error = function(e) {
    cat("特征选择过程出错:", e$message, "\n")
    cat("将使用所有可用特征继续\n")
    selected_features <- key_features
  })
  
  # 如果LASSO未选择任何特征，使用所有特征
  if(length(selected_features) == 0) {
    selected_features <- key_features
    cat("使用所有可用特征:", length(selected_features), "个\n")
  }
  
  # 使用选定的特征进行模型训练
  final_model_data <- model_data[, c(selected_features, "vision_improvement_1m")]
  
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
  formula_obj <- as.formula(paste("vision_improvement_1m ~", paste(selected_features, collapse = " + ")))
  
  # 线性回归 - 添加错误处理以应对奇异矩阵
  lm_result <- tryCatch({
    # 检查特征数量是否过多
    if(length(selected_features) > nrow(final_model_data) * 0.7) {
      # 使用岭回归代替普通线性回归
      cat("特征数量过多，使用岭回归代替线性回归\n")
      ridge_grid <- expand.grid(alpha = 0, lambda = 10^seq(-3, 2, length.out = 10))
      
      model <- train(
        formula_obj,
        data = final_model_data,
        method = "glmnet",
        tuneGrid = ridge_grid,
        trControl = ctrl,
        metric = "RMSE"
      )
      perf <- calculate_performance(model)
      perf$Model <- "Ridge Regression"
      perf
    } else {
      # 标准线性回归
      model <- train(
        formula_obj,
        data = final_model_data,
        method = "lm",
        trControl = ctrl,
        metric = "RMSE"
      )
      perf <- calculate_performance(model)
      perf$Model <- "Linear Regression"
      perf
    }
  }, error = function(e) {
    cat("线性回归模型训练失败:", e$message, "\n")
    cat("尝试使用岭回归替代\n")
    
    # 如果线性回归失败，尝试岭回归
    tryCatch({
      ridge_grid <- expand.grid(alpha = 0, lambda = 10^seq(-3, 2, length.out = 10))
      
      model <- train(
        formula_obj,
        data = final_model_data,
        method = "glmnet",
        tuneGrid = ridge_grid,
        trControl = ctrl,
        metric = "RMSE"
      )
      perf <- calculate_performance(model)
      perf$Model <- "Ridge Regression"
      perf
    }, error = function(e2) {
      cat("岭回归也失败了:", e2$message, "\n")
      return(NULL)
    })
  })
  
  if(!is.null(lm_result)) {
    results_list$lm <- lm_result
    cat(lm_result$Model, "模型训练成功, R² =", round(lm_result$R2_mean, 3), "\n")
  }
  
  # LASSO
  lasso_result <- tryCatch({
    lasso_grid <- expand.grid(
      alpha = 1,
      lambda = 10^seq(-4, 0, length.out = 10)  # 使用更宽的对数尺度范围
    )
    
    model <- train(
      formula_obj,
      data = final_model_data,
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
      alpha = c(0.2, 0.5, 0.8),  # 更均匀的混合参数
      lambda = 10^seq(-4, 0, length.out = 5)  # 更宽的正则化参数范围
    )
    
    model <- train(
      formula_obj,
      data = final_model_data,
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
      mtry = c(2, max(2, floor(sqrt(length(selected_features)))))
    )
    
    model <- train(
      formula_obj,
      data = final_model_data,
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
      data = final_model_data,
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
    combined_results$SampleSize <- nrow(final_model_data)
    combined_results$Features <- paste(selected_features, collapse="|")
    combined_results$FeatureCount <- length(selected_features)
    
    cat("分析完成，成功训练了", nrow(combined_results), "个模型\n")
    return(combined_results)
  } else {
    cat("没有产生任何有效结果\n")
    return(NULL)
  }
}

# 寻找所有CSV文件 - 仅匹配指定的时间窗口
all_csv_files <- list.files(data_dir, pattern = "^(pre|post)_.+_(NoD_Surg0|D_Surg1)\\.csv$", full.names = TRUE)
csv_files <- character(0)

# 遍历所有文件，只保留特定时间窗口的文件
for(file in all_csv_files) {
  file_name <- basename(file)
  file_info <- parse_filename(file_name)
  
  if(file_info$success && file_info$day_label %in% time_windows) {
    csv_files <- c(csv_files, file)
  }
}

if(length(csv_files) == 0) {
  stop("未找到任何符合指定时间窗口的CSV文件")
}

cat("找到", length(csv_files), "个符合时间窗口的CSV文件\n")

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
  
  return(heatmap_plot)
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
  
  # 创建带有标准时间窗口顺序的线图
  # 确保Day变量使用标准时间窗口顺序（当存在于结果中时）
  # 这将帮助线图显示按照正确的时间顺序
  if(nrow(all_results) > 0) {
    # 检查哪些标准时间窗口存在于我们的结果中
    present_windows <- unique(all_results$Day)
    valid_window_order <- intersect(time_windows, present_windows)
    
    if(length(valid_window_order) > 1) {
      # 如果有多个标准窗口，创建有序因子
      all_results$Day <- factor(all_results$Day, 
                                levels = valid_window_order, 
                                ordered = TRUE)
      cat("使用标准时间窗口顺序:", paste(valid_window_order, collapse=", "), "\n")
    }
  }
  
  # 1. 为每个组创建线图 (R²随时间段变化)
  message("创建组别模型表现线图...")
  
  # 使用map函数代替for循环
  purrr::walk(groups, function(group_name) {
    group_data <- all_results %>% filter(Group == group_name)
    
    if(nrow(group_data) == 0) return(NULL)
    
    # 线图 - 不依赖于Day的有序性
    line_plot <- ggplot(group_data, aes(x = Day, y = R2_mean, color = Model, group = Model)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      geom_text(aes(label = SizeLbl), vjust = -1.5, size = 2.5, show.legend = FALSE) +
      labs(
        title = paste("Vision Improvement Prediction -", group_name),
        subtitle = "Comparing Different Machine Learning Models with LOOCV",
        x = "Time Period Relative to Surgery",
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
      
      # 线图比较 - 修改x轴标签
      comparison_plot <- ggplot(model_data, aes(x = Day, y = R2_mean, color = Group, group = Group)) +
        geom_line(size = 1) +
        geom_point(size = 3) +
        geom_text(aes(label = SizeLbl), vjust = -1.5, size = 2.5, show.legend = FALSE) +
        labs(
          title = paste(model_name, "- Group Comparison"),
          subtitle = "Comparing prediction performance between groups using LOOCV",
          x = "Time Period Relative to Surgery",
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

# 检查是否有结果
if(nrow(all_results) > 0) {
  # 保存结果到CSV文件
  results_file <- file.path(output_dir, "prediction_results.csv")
  write.csv(all_results, results_file, row.names = FALSE)
  cat("结果已保存到:", results_file, "\n")
  
  # 保存特征选择结果
  features_file <- file.path(output_dir, "selected_features.csv")
  features_data <- all_results %>% 
    distinct(Day, Group, Features, FeatureCount) %>%
    arrange(Day, Group)
  write.csv(features_data, features_file, row.names = FALSE)
  cat("特征选择结果已保存到:", features_file, "\n")
  
  # 调用绘图函数创建所有图表
  create_all_plots(all_results, output_dir)
  
  cat("分析完成，结果和图表已保存到:", output_dir, "\n")
} else {
  cat("没有生成任何结果，请检查数据和模型训练过程\n")
}
