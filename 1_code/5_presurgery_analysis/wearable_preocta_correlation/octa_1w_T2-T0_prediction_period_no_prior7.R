# 术后OCTA指标改善预测模型 - 按时间窗口分析
# 使用不同时间窗口的可穿戴设备数据预测术后1个月OCTA指标改善值

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
# 1. 定义标准时间窗口 (从第一段代码中提取)
#-------------------------------

# 定义标准时间窗口
time_windows <- c(
  "pre_7d_all",   # 术前7天及以前
  "pre_3d_7d",    # 术前3-7天
  "pre_3d",       # 术前3天
  "post_3d",      # 术后3天
  "post_4d_6d",  # 术后4-6天
  "post_6d_all" # 术后7天
)

# 创建输出目录
output_dir <- "3_data_analysis/5_presurgery_analysis/octa_prediction/results_by_timewindow/1w/no_prior7"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------
# 2. 加载数据
#-------------------------------

# 加载基础数据 - 使用日期信息
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1w_prediction/daily_data/all_days_combined_data.csv")
octa_bloodflow <- read.csv("2_data/analysis_data/octa_data_bloodflow.csv")
octa_thickness <- read.csv("2_data/analysis_data/octa_data_thickness.csv")
baseline_info <- read.csv("2_data/analysis_data/baseline_info.csv")

#-------------------------------
# 3. 处理OCTA数据 - 区分血流和厚度
#-------------------------------

# 处理OCTA数据函数
process_octa_data <- function(baseline_data, octa_data, id_column = "id") {
  features <- baseline_data %>%
    filter(!is.na(surgery_eye_1)) %>%
    distinct(ID, surgery_eye_1, .keep_all = TRUE) %>%
    left_join(octa_data, by = c("ID" = id_column))
  
  return(features)
}

octa_bloodflow_features <- process_octa_data(baseline_info, octa_bloodflow)
octa_thickness_features <- process_octa_data(baseline_info, octa_thickness)

# 处理每个患者数据的函数
process_patient_data <- function(patient_data, time_points = c("T0", "T1")) {
  current_eye <- patient_data$surgery_eye_1[1]
  # 左眼手术用左眼数据，右眼和双眼手术用右眼数据
  pattern <- if(current_eye == 1) "_OS_" else "_OD_"
  
  result <- patient_data %>% dplyr::select(ID)
  
  # 处理每个时间点
  for(suffix in time_points) {
    # 选择当前时间点的列
    cols_to_keep <- grep(pattern, names(patient_data), value = TRUE)
    cols_to_keep <- cols_to_keep[grep(paste0(suffix, "$"), cols_to_keep)]
    
    if(length(cols_to_keep) > 0) {
      # 选择数据并重命名列
      time_data <- patient_data %>%
        dplyr::select("ID", all_of(cols_to_keep)) %>%
        rename_with(~ gsub("_(OD|OS)_", "_", .), -ID)
      
      result <- result %>% left_join(time_data, by = "ID")
    }
  }
  
  return(result)
}

# 处理所有患者的数据
process_all_patients <- function(features_data) {
  patient_list <- split(features_data, features_data$ID)
  processed_data <- purrr::map(patient_list, process_patient_data)
  return(bind_rows(processed_data))
}

# 分别处理血流和厚度数据
octa_bloodflow_features <- process_all_patients(octa_bloodflow_features)
octa_thickness_features <- process_all_patients(octa_thickness_features)

#-------------------------------
# 4. 计算OCTA指标的改善值(T1-T0)
#-------------------------------

# 计算改善值函数
calculate_improvement <- function(data, data_type = "unknown") {
  # 初始化结果数据框
  result <- data %>% dplyr::select(ID)
  
  # 对每个T0变量找到对应的T1变量并计算差值
  vars_T0 <- names(data)[grep("_T0$", names(data))]
  improvement_count <- 0
  
  for(t0_var in vars_T0) {
    t1_var <- gsub("_T0$", "_T1", t0_var)
    
    if(t1_var %in% names(data)) {
      # 构建改善值变量名
      imp_var <- gsub("_T0$", "_improvement", t0_var)
      
      # 计算改善值(T1-T0)
      result[[imp_var]] <- data[[t1_var]] - data[[t0_var]]
      improvement_count <- improvement_count + 1
    }
  }
  
  cat("计算了", improvement_count, "个", data_type, "改善值变量\n")
  return(result)
}

# 计算血流和厚度改善值
bloodflow_improvement <- calculate_improvement(octa_bloodflow_features, "血流")
thickness_improvement <- calculate_improvement(octa_thickness_features, "厚度")

#-------------------------------
# 5. 按时间窗口准备可穿戴设备预测变量
#-------------------------------

# 创建一个函数，根据时间窗口筛选数据并计算均值
prepare_timewindow_data <- function(data, time_window) {
  # 根据时间窗口定义日期范围
  date_range <- switch(time_window,
                       "pre_7d_all" = c(-6, -1),      # 术前7天以前(假设最多考虑术前30天)
                       "pre_4d_7d" = c(-6, -4),        # 术前3-7天
                       "pre_3d" = c(-3, -1),           # 术前3天
                       "post_3d" = c(1, 3),            # 术后3天
                       "post_4d_6d" = c(4, 6),       # 术后4-6天
                       "post_6d_all" = c(1, 6),    # 术后1-6天
                       c(-30, 30)                       # 默认范围，应该不会用到
  )
  
  # 筛选指定时间窗口的数据
  filtered_data <- data %>%
    filter(as.numeric(day) >= date_range[1] & as.numeric(day) <= date_range[2])
  
  # 按受试者ID分组并计算均值
  aggregated_data <- filtered_data %>%
    group_by(subject_id) %>%
    summarise(
      # 选择关键均值指标
      mean_rhr_1 = mean(mean_rhr_1, na.rm = TRUE),     # 平均静息心率
      mean_bo = mean(mean_bo, na.rm = TRUE),           # 平均血氧饱和度
      steps_total = mean(steps_total, na.rm = TRUE),   # 平均总步数
      total_sleep = mean(total_sleep, na.rm = TRUE),   # 平均总睡眠时间
      
      # 保留其他人口统计变量
      dm_2 = dplyr::first(dm_2),                    # 糖尿病状态
      age = dplyr::first(age),                      # 年龄
      gender = dplyr::first(gender),                # 性别
      bmi = dplyr::first(bmi),                      # BMI
      hypertension_2 = dplyr::first(hypertension_2) # 高血压状态
    )
  
  # 中位数填补缺失值
  for(col in c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep")) {
    if(sum(is.na(aggregated_data[[col]])) > 0) {
      aggregated_data[[col]][is.na(aggregated_data[[col]])] <- median(aggregated_data[[col]], na.rm = TRUE)
    }
  }
  
  return(aggregated_data)
}

#-------------------------------
# 6. 定义要分析的变量
#-------------------------------

# 血流变量 - PA_前缀为血流灌注区域(Perfusion Area)
significant_bloodflow_pairs <- list(
  list(improvement_var = "PA_Vitreous_0_6_improvement")
)

# 厚度变量
significant_thickness_pairs <- list(
  list(improvement_var = "Thickness_OuterRetina_0_6_improvement"),
  list(improvement_var = "Thickness_Retina_0_6_improvement"),
  list(improvement_var = "Thickness_RNFL_0_6_improvement")
)

# 定义可穿戴设备的预测变量
wearable_vars <- list(
  "mean_rhr_1",    # 静息心率
  "mean_bo",       # 血氧
  "steps_total",   # 步数
  "total_sleep"    # 睡眠时间
)

#-------------------------------
# 7. 创建一个函数来存储每个时间窗口的模型性能
#-------------------------------

# 预测模型构建函数 - 适用于血流和厚度改善
build_prediction_models <- function(model_data, target_pairs, model_type = "unknown", time_window = "unknown") {
  # 创建数据框来存储每个目标变量和预测器组合的模型性能指标
  performance_df <- data.frame(
    time_window = character(),
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
  
  # 对每个目标变量进行建模
  for(pair in target_pairs) {
    improvement_var <- pair$improvement_var
    
    # 检查改善值变量是否存在
    if(!(improvement_var %in% names(model_data))) {
      cat(model_type, "改善值变量不存在:", improvement_var, "\n")
      next
    }
    
    cat("\n===== 构建", model_type, "改善预测模型:", improvement_var, "- 时间窗口:", time_window, "=====\n")
    
    # 准备基本模型数据 - 移除含有NA的行
    base_vars <- c(improvement_var, unlist(wearable_vars), "age", "gender", "bmi", "dm_2")
    base_data <- model_data %>%
      dplyr::select(all_of(base_vars)) %>%
      na.omit()
    
    if(nrow(base_data) < 10) {  # 确保有足够的样本
      cat("样本量不足:", improvement_var, "\n")
      next
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
      improvement_var, "~", 
      wearable_formula_str, " + ",
      "age + gender + bmi + dm_2"
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
        time_window = time_window,
        target_var = improvement_var,
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
        improvement_var, "~", 
        wearable_var, "+", 
        "age + gender + bmi + dm_2"
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
          time_window = time_window,
          target_var = improvement_var,
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
  }
  
  return(performance_df)
}

#-------------------------------
# 8. 创建存储所有模型性能的数据框
#-------------------------------

# 为血流和厚度数据分别创建存储性能的数据框
all_bloodflow_performance <- data.frame()
all_thickness_performance <- data.frame()

#-------------------------------
# 9. 对每个时间窗口运行预测模型
#-------------------------------

cat("\n开始为每个时间窗口构建OCTA改善预测模型...\n")

for(time_window in time_windows) {
  cat("\n=== 处理时间窗口:", time_window, "===\n")
  
  # 获取当前时间窗口的可穿戴设备数据
  current_window_data <- prepare_timewindow_data(combined_data, time_window)
  
  # 合并血流和厚度数据
  model_data_bloodflow <- current_window_data %>%
    left_join(bloodflow_improvement, by = c("subject_id" = "ID"))
  
  model_data_thickness <- current_window_data %>%
    left_join(thickness_improvement, by = c("subject_id" = "ID"))
  
  # 构建血流改善预测模型
  cat("\n开始构建血流改善预测模型 - 时间窗口:", time_window, "...\n")
  bloodflow_performance <- tryCatch({
    build_prediction_models(model_data_bloodflow, significant_bloodflow_pairs, "bloodflow", time_window)
  }, error = function(e) {
    cat("Error in bloodflow model:", e$message, "\n")
    return(data.frame())
  })
  
  # 构建厚度改善预测模型
  cat("\n开始构建厚度改善预测模型 - 时间窗口:", time_window, "...\n")
  thickness_performance <- tryCatch({
    build_prediction_models(model_data_thickness, significant_thickness_pairs, "thickness", time_window)
  }, error = function(e) {
    cat("Error in thickness model:", e$message, "\n")
    return(data.frame())
  })
  
  # 将当前窗口的性能结果追加到总结果
  all_bloodflow_performance <- rbind(all_bloodflow_performance, bloodflow_performance)
  all_thickness_performance <- rbind(all_thickness_performance, thickness_performance)
}

#-------------------------------
# 10. 保存模型性能结果
#-------------------------------

# 合并所有性能结果
all_performance <- rbind(
  all_bloodflow_performance,
  all_thickness_performance
)

# 保存到CSV文件
write.csv(all_performance, file.path(output_dir, "octa_improvement_prediction_by_timewindow.csv"), row.names = FALSE)
cat("\n所有模型性能结果已保存到:", file.path(output_dir, "octa_improvement_prediction_by_timewindow.csv"), "\n")

#-------------------------------
# 11. 创建性能变化趋势图表 (折线图)
#-------------------------------

# 按指标和预测变量绘制不同时间窗口的模型性能折线图
create_trend_plots <- function(performance_data, data_type) {
  # 获取所有唯一的目标变量
  unique_targets <- unique(performance_data$target_var)
  
  # 对每个目标变量创建单独的图表
  for(target in unique_targets) {
    # 筛选当前目标变量的数据
    target_data <- performance_data %>%
      filter(target_var == target)
    
    # 准备绘图数据 - 按预测变量分组
    plot_data <- target_data %>%
      select(time_window, predictor, cv_r_squared) %>%
      # 确保时间窗口按照定义的顺序排序
      mutate(time_window = factor(time_window, levels = time_windows))
    
    # 简化指标名称
    simple_target <- gsub("_improvement$", "", target)
    simple_target <- gsub("^(PA|VD|SVD|Thickness)_", "", simple_target)
    
    # 为不同预测变量定义颜色
    color_values <- c(
      "All Predictors" = "#513a52", 
      "mean_rhr_1" = "#82a1bf", 
      "mean_bo" = "#faaa93", 
      "steps_total" = "#feefc4",
      "total_sleep" = "#7fb77e"  # 添加睡眠变量的颜色
    )
    
    # 创建折线图
    p <- ggplot(plot_data, aes(x = time_window, y = cv_r_squared, color = predictor, group = predictor)) +
      geom_line(size = 1) +
      geom_point(size = 3) +
      scale_color_manual(
        values = color_values,
        labels = c(
          "All Predictors" = "All Predictors", 
          "mean_rhr_1" = "Heart Rate", 
          "mean_bo" = "Blood Oxygen", 
          "steps_total" = "Total Steps",
          "total_sleep" = "Sleep Duration"
        )
      ) +
      labs(
        title = paste("Prediction of", ifelse(data_type == "bloodflow", "Blood Flow", "Thickness"), "Improvement:", simple_target),
        x = "Time Window",
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
    plot_filename <- paste0(data_type, "_", simple_target, "_performance_trend.pdf")
    ggsave(file.path(output_dir, plot_filename), p, width = 10, height = 6)
    
    # 同时保存PNG格式
    plot_filename_png <- paste0(data_type, "_", simple_target, "_performance_trend.png")
    ggsave(file.path(output_dir, plot_filename_png), p, width = 10, height = 6, dpi = 300)
    
    cat("已保存", data_type, "改善预测趋势图:", simple_target, "\n")
  }
}

# 绘制血流改善预测模型性能趋势图
create_trend_plots(all_bloodflow_performance, "bloodflow")

# 绘制厚度改善预测模型性能趋势图
create_trend_plots(all_thickness_performance, "thickness")

#-------------------------------
# 12. 创建热图显示不同时间窗口和预测变量的性能
#-------------------------------

create_heatmaps <- function(performance_data, data_type) {
  # 获取所有唯一的目标变量
  unique_targets <- unique(performance_data$target_var)
  
  # 对每个目标变量创建单独的热图
  for(target in unique_targets) {
    # 筛选当前目标变量的数据
    target_data <- performance_data %>%
      filter(target_var == target)
    
    # 简化指标名称
    simple_target <- gsub("_improvement$", "", target)
    simple_target <- gsub("^(PA|VD|SVD|Thickness)_", "", simple_target)
    
    # 创建热图数据 - 排除"All Predictors"
    heatmap_data <- target_data %>%
      filter(predictor != "All Predictors") %>%
      select(time_window, predictor, cv_r_squared) %>%
      # 确保时间窗口按照定义的顺序排序
      mutate(
        time_window = factor(time_window, levels = time_windows),
        predictor = factor(predictor, levels = c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep"),
                           labels = c("Heart Rate", "Blood Oxygen", "Total Steps", "Sleep Duration"))
      )
    
    # 创建热图
    p_heatmap <- ggplot(heatmap_data, aes(x = time_window, y = predictor, fill = cv_r_squared)) +
      geom_tile() +
      scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                           midpoint = median(heatmap_data$cv_r_squared, na.rm = TRUE)) +
      labs(
        title = paste("Predictive Power for", ifelse(data_type == "bloodflow", "Blood Flow", "Thickness"), "Improvement:", simple_target),
        x = "Time Window",
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
    heatmap_filename <- paste0(data_type, "_", simple_target, "_heatmap.pdf")
    ggsave(file.path(output_dir, heatmap_filename), p_heatmap, width = 10, height = 6)
    
    # 同时保存PNG格式
    heatmap_filename_png <- paste0(data_type, "_", simple_target, "_heatmap.png")
    ggsave(file.path(output_dir, heatmap_filename_png), p_heatmap, width = 10, height = 6, dpi = 300)
    
    cat("已保存", data_type, "改善预测热图:", simple_target, "\n")
  }
}

# 绘制血流改善预测模型热图
create_heatmaps(all_bloodflow_performance, "bloodflow")

# 绘制厚度改善预测模型热图
create_heatmaps(all_thickness_performance, "thickness")

