# 每日数据OCTA指标改善预测模型 (1个月)
# 使用每天的可穿戴设备数据预测术后1个月OCTA指标改善值

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
# 1. 定义每日时间点
#-------------------------------

# 定义每日时间点（从术前7天到术后29天）
daily_timepoints <- c(
   "-6", "-5", "-4", "-3", "-2", "-1",  # 术前7天到术前1天
  "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",  # 术后1-10天
  "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",  # 术后11-20天
  "21", "22", "23", "24", "25", "26", "27", "28", "29"  # 术后21-29天
)

# 创建输出目录
output_dir <- "3_data_analysis/5_presurgery_analysis/octa_prediction/results_by_day/1m/no_prior7"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------
# 2. 加载数据
#-------------------------------

# 加载基础数据 - 使用日期信息
combined_data <- read_csv("3_data_analysis/3_prediction_modeling/1m_prediction/daily_data/all_days_combined_data.csv")
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
process_patient_data <- function(patient_data, time_points = c("T0", "T2")) {
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
# 4. 计算OCTA指标的改善值(T2-T0)
#-------------------------------

# 计算改善值函数
calculate_improvement <- function(data, data_type = "unknown") {
  # 初始化结果数据框
  result <- data %>% dplyr::select(ID)
  
  # 对每个T0变量找到对应的T2变量并计算差值
  vars_T0 <- names(data)[grep("_T0$", names(data))]
  improvement_count <- 0
  
  for(t0_var in vars_T0) {
    t2_var <- gsub("_T0$", "_T2", t0_var)
    
    if(t2_var %in% names(data)) {
      # 构建改善值变量名
      imp_var <- gsub("_T0$", "_improvement", t0_var)
      
      # 计算改善值(T2-T0)
      result[[imp_var]] <- data[[t2_var]] - data[[t0_var]]
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
# 5. 按每日时间点准备可穿戴设备预测变量
#-------------------------------

# 创建一个函数，根据特定日期筛选数据
prepare_daily_data <- function(data, day_point) {
  # 筛选指定日期的数据
  filtered_data <- data %>%
    filter(as.numeric(day) == as.numeric(day_point))
  
  # 如果筛选后数据量很少，不继续处理
  if(nrow(filtered_data) < 5) {
    cat("日期", day_point, "的样本数量不足:", nrow(filtered_data), "\n")
    return(NULL)
  }
  
  # 数据预处理
  processed_data <- filtered_data %>%
    select(
      subject_id,
      # 选择关键指标
      mean_rhr_1,     # 平均静息心率
      mean_bo,        # 平均血氧饱和度
      steps_total,    # 平均总步数
      total_sleep,    # 平均总睡眠时间
      
      # 保留其他人口统计变量
      dm_2,           # 糖尿病状态
      age,            # 年龄
      gender,         # 性别
      bmi,            # BMI
      hypertension_2  # 高血压状态
    )
  
  # 中位数填补缺失值
  for(col in c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep")) {
    if(sum(is.na(processed_data[[col]])) > 0) {
      processed_data[[col]][is.na(processed_data[[col]])] <- median(processed_data[[col]], na.rm = TRUE)
    }
  }
  
  return(processed_data)
}

#-------------------------------
# 6. 定义要分析的变量
#-------------------------------

# 血流变量 - PA_前缀为血流灌注区域(Perfusion Area)
significant_bloodflow_pairs <- list(
  list(improvement_var = "PA_Choriocapillaris_0_21_improvement"),
  list(improvement_var = "VD_DCP_0_21_improvement")
)

# 厚度变量
significant_thickness_pairs <- list(
  list(improvement_var = "Thickness_GCL.IPL_0_21_improvement"),
  list(improvement_var = "Thickness_RNFL_0_21_improvement")
)

# 定义可穿戴设备的预测变量
wearable_vars <- list(
  "mean_rhr_1",    # 静息心率
  "mean_bo",       # 血氧
  "steps_total",   # 步数
  "total_sleep"    # 睡眠时间
)

#-------------------------------
# 7. 创建一个函数来存储每个日期的模型性能
#-------------------------------

# 预测模型构建函数 - 适用于血流和厚度改善
build_prediction_models <- function(model_data, target_pairs, model_type = "unknown", day_point = "unknown") {
  # 创建数据框来存储每个目标变量和预测器组合的模型性能指标
  performance_df <- data.frame(
    day_point = character(),
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
    
    cat("\n===== 构建", model_type, "改善预测模型:", improvement_var, "- 日期:", day_point, "=====\n")
    
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
        day_point = day_point,
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
          day_point = day_point,
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

# 创建存储每个日期样本量的数据框
sample_size_data <- data.frame(
  day_point = character(),
  sample_size = numeric(),
  stringsAsFactors = FALSE
)

#-------------------------------
# 9. 对每个日期运行预测模型
#-------------------------------

cat("\n开始为每个日期构建OCTA改善预测模型...\n")

for(day_point in daily_timepoints) {
  cat("\n=== 处理日期:", day_point, "===\n")
  
  # 获取当前日期的可穿戴设备数据
  current_day_data <- prepare_daily_data(combined_data, day_point)
  
  # 如果该日期的数据不足，跳过
  if(is.null(current_day_data)) {
    cat("跳过日期", day_point, "，因为数据不足\n")
    # 记录样本量为0
    sample_size_data <- rbind(sample_size_data, data.frame(
      day_point = day_point,
      sample_size = 0,
      stringsAsFactors = FALSE
    ))
    next
  }
  
  # 记录该日期的样本量
  sample_size_data <- rbind(sample_size_data, data.frame(
    day_point = day_point,
    sample_size = nrow(current_day_data),
    stringsAsFactors = FALSE
  ))
  
  # 合并血流和厚度数据
  model_data_bloodflow <- current_day_data %>%
    left_join(bloodflow_improvement, by = c("subject_id" = "ID"))
  
  model_data_thickness <- current_day_data %>%
    left_join(thickness_improvement, by = c("subject_id" = "ID"))
  
  # 构建血流改善预测模型
  cat("\n开始构建血流改善预测模型 - 日期:", day_point, "...\n")
  bloodflow_performance <- tryCatch({
    build_prediction_models(model_data_bloodflow, significant_bloodflow_pairs, "bloodflow", day_point)
  }, error = function(e) {
    cat("Error in bloodflow model:", e$message, "\n")
    return(data.frame())
  })
  
  # 构建厚度改善预测模型
  cat("\n开始构建厚度改善预测模型 - 日期:", day_point, "...\n")
  thickness_performance <- tryCatch({
    build_prediction_models(model_data_thickness, significant_thickness_pairs, "thickness", day_point)
  }, error = function(e) {
    cat("Error in thickness model:", e$message, "\n")
    return(data.frame())
  })
  
  # 将当前日期的性能结果追加到总结果
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
write.csv(all_performance, file.path(output_dir, "octa_improvement_prediction_by_day.csv"), row.names = FALSE)
cat("\n所有模型性能结果已保存到:", file.path(output_dir, "octa_improvement_prediction_by_day.csv"), "\n")

# 保存每日样本量数据
write.csv(sample_size_data, file.path(output_dir, "daily_sample_size.csv"), row.names = FALSE)
cat("\n每日样本量数据已保存到:", file.path(output_dir, "daily_sample_size.csv"), "\n")

#-------------------------------
# 11. 创建性能变化趋势图表 (折线图)
#-------------------------------

# 按指标和预测变量绘制不同日期的模型性能折线图
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
      select(day_point, predictor, cv_r_squared) %>%
      # 确保日期按照定义的顺序排序
      mutate(day_point = factor(day_point, levels = daily_timepoints))
    
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
    p <- ggplot(plot_data, aes(x = day_point, y = cv_r_squared, color = predictor, group = predictor)) +
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
        x = "Day Relative to Surgery",
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
    plot_filename <- paste0(data_type, "_", simple_target, "_daily_performance_trend.pdf")
    ggsave(file.path(output_dir, plot_filename), p, width = 10, height = 6)
    
    # 同时保存PNG格式
    plot_filename_png <- paste0(data_type, "_", simple_target, "_daily_performance_trend.png")
    ggsave(file.path(output_dir, plot_filename_png), p, width = 10, height = 6, dpi = 300)
    
    cat("已保存", data_type, "改善预测趋势图:", simple_target, "\n")
  }
}

# 绘制血流改善预测模型性能趋势图
create_trend_plots(all_bloodflow_performance, "bloodflow")

# 绘制厚度改善预测模型性能趋势图
create_trend_plots(all_thickness_performance, "thickness")

#-------------------------------
# 12. 创建热图显示不同日期和预测变量的性能
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
      select(day_point, predictor, cv_r_squared) %>%
      # 确保日期按照定义的顺序排序
      mutate(
        day_point = factor(day_point, levels = daily_timepoints),
        predictor = factor(predictor, levels = c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep"),
                           labels = c("Heart Rate", "Blood Oxygen", "Total Steps", "Sleep Duration"))
      )
    
    # 创建热图
    p_heatmap <- ggplot(heatmap_data, aes(x = day_point, y = predictor, fill = cv_r_squared)) +
      geom_tile() +
      scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                           midpoint = median(heatmap_data$cv_r_squared, na.rm = TRUE)) +
      labs(
        title = paste("Predictive Power for", ifelse(data_type == "bloodflow", "Blood Flow", "Thickness"), "Improvement:", simple_target),
        x = "Day Relative to Surgery",
        y = "Predictor",
        fill = "Cross-validated R²"
      ) + 
      scale_y_discrete(
        labels = c(
          "Heart Rate" = "Heart Rate", 
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
    heatmap_filename <- paste0(data_type, "_", simple_target, "_daily_heatmap.pdf")
    ggsave(file.path(output_dir, heatmap_filename), p_heatmap, width = 10, height = 6)
    
    # 同时保存PNG格式
    heatmap_filename_png <- paste0(data_type, "_", simple_target, "_daily_heatmap.png")
    ggsave(file.path(output_dir, heatmap_filename_png), p_heatmap, width = 10, height = 6, dpi = 300)
    
    cat("已保存", data_type, "改善预测热图:", simple_target, "\n")
  }
}

# 绘制血流改善预测模型热图
create_heatmaps(all_bloodflow_performance, "bloodflow")

# 绘制厚度改善预测模型热图
create_heatmaps(all_thickness_performance, "thickness")

#-------------------------------
# 13. 创建每日样本量图表
#-------------------------------

# 绘制每日样本量柱状图
plot_sample_size <- function(sample_data) {
  # 准备绘图数据
  plot_data <- sample_data %>%
    mutate(
      day_point = factor(day_point, levels = daily_timepoints),
      day_label = ifelse(as.numeric(day_point) < 0, 
                         paste0("Day ", day_point), 
                         paste0("Day ", day_point))
    )
  
  # 创建柱状图
  p <- ggplot(plot_data, aes(x = day_point, y = sample_size)) +
    geom_bar(stat = "identity", width = 0.7, fill = "#4292c6") +
    geom_text(aes(label = sample_size), vjust = -0.5, size = 3.5) +
    labs(
      title = "Daily sample size distribution",
      x = "Relative date to surgery",
      y = "Sample size"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) # 为文本标签留出空间
  
  return(p)
}

# 绘制样本量图表
p_sample_size <- plot_sample_size(sample_size_data)

# 保存样本量图表
sample_size_filename <- "daily_sample_size.pdf"
ggsave(file.path(output_dir, sample_size_filename), p_sample_size, width = 10, height = 6)





#-------------------------------
# 14. 二分类OCTA指标改善预测
#-------------------------------

cat("\n开始二分类OCTA指标改善预测分析...\n")

# 创建二分类改善变量的函数
create_binary_improvement <- function(data, data_type = "unknown") {
  # 初始化结果数据框
  result <- data %>% dplyr::select(ID)
  
  # 对每个连续改善值变量创建二分类变量
  improvement_vars <- names(data)[grep("_improvement$", names(data))]
  binary_count <- 0
  
  for(imp_var in improvement_vars) {
    # 构建二分类变量名
    bin_var <- gsub("_improvement$", "_binary_improvement", imp_var)
    
    # 创建二分类变量 - 血流增加为改善，厚度下降为改善
    if(grepl("^PA_|^VD_|^SVD_", imp_var)) {  # 血流变量
      result[[bin_var]] <- ifelse(data[[imp_var]] > 0, 1, 0)  # 血流增加(>0)为改善
    } else {  # 厚度变量
      result[[bin_var]] <- ifelse(data[[imp_var]] < 0, 1, 0)  # 厚度减少(<0)为改善
    }
    
    binary_count <- binary_count + 1
  }
  
  cat("创建了", binary_count, "个", data_type, "二分类改善变量\n")
  return(result)
}

# 为血流和厚度数据创建二分类改善变量
bloodflow_binary_improvement <- create_binary_improvement(bloodflow_improvement, "血流")
thickness_binary_improvement <- create_binary_improvement(thickness_improvement, "厚度")

# 定义要分析的二分类变量
significant_binary_bloodflow_pairs <- list(
  list(improvement_var = "PA_Choriocapillaris_0_21_binary_improvement"),
  list(improvement_var = "VD_DCP_0_21_binary_improvement")
)

significant_binary_thickness_pairs <- list(
  list(improvement_var = "Thickness_GCL.IPL_0_21_binary_improvement"),
  list(improvement_var = "Thickness_RNFL_0_21_binary_improvement")
)

#-------------------------------
# 15. 二分类预测模型构建函数
#-------------------------------

build_simple_binary_model <- function(model_data, target_pairs, model_type = "unknown", day_point = "unknown") {
  # 创建数据框来存储每个目标变量和预测器组合的模型性能指标
  performance_df <- data.frame(
    day_point = character(),
    target_var = character(),
    predictor = character(),
    auc = numeric(),
    sensitivity = numeric(),
    specificity = numeric(),
    accuracy = numeric(),
    sample_size = numeric(),
    predictor_coef = numeric(),
    predictor_p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 对每个目标变量进行建模
  for(pair in target_pairs) {
    improvement_var <- pair$improvement_var
    
    # 检查二分类改善变量是否存在
    if(!(improvement_var %in% names(model_data))) {
      cat(model_type, "二分类改善变量不存在:", improvement_var, "\n")
      next
    }
    
    cat("\n===== 构建", model_type, "二分类改善预测模型:", improvement_var, "- 日期:", day_point, "=====\n")
    
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
    
    # 计算目标变量的分布
    improvement_count <- sum(base_data[[improvement_var]] == 1)
    no_improvement_count <- sum(base_data[[improvement_var]] == 0)
    cat("改善数量:", improvement_count, "未改善数量:", no_improvement_count, "\n")
    
    # 验证是否有足够的类别分布
    if(improvement_count < 2 || no_improvement_count < 2) {
      cat("警告: 某一类别的样本数量太少，无法构建稳健的模型\n")
      next
    }
    
    # 合并可穿戴设备变量为公式字符串
    wearable_formula_str <- paste(unlist(wearable_vars), collapse = " + ")
    
    # 运行所有变量的模型
    run_model_and_record <- function(formula_str, predictor_name, base_data_input) {
      formula_obj <- as.formula(formula_str)
      
      # 构建简单逻辑回归模型
      model <- tryCatch({
        glm(formula_obj, data = base_data_input, family = binomial())
      }, error = function(e) {
        cat("构建模型错误:", e$message, "\n")
        return(NULL)
      })
      
      if(is.null(model)) return(NULL)
      
      # 获取预测概率
      prob_predictions <- tryCatch({
        predict(model, type = "response")
      }, error = function(e) {
        cat("预测错误:", e$message, "\n")
        return(NULL)
      })
      
      if(is.null(prob_predictions)) return(NULL)
      
      # 创建预测类别
      pred_class <- ifelse(prob_predictions > 0.5, 1, 0)
      actual_class <- base_data_input[[improvement_var]]
      
      # 创建混淆矩阵
      conf_matrix <- table(Predicted = pred_class, Actual = actual_class)
      cat("混淆矩阵:\n")
      print(conf_matrix)
      
      # 计算性能指标
      total <- sum(conf_matrix)
      accuracy <- sum(diag(conf_matrix)) / total
      
      # 只有在混淆矩阵是2x2时才计算敏感性和特异性
      sensitivity <- specificity <- NA
      if(nrow(conf_matrix) == 2 && ncol(conf_matrix) == 2) {
        true_pos <- conf_matrix[2, 2]
        false_neg <- conf_matrix[1, 2]
        true_neg <- conf_matrix[1, 1]
        false_pos <- conf_matrix[2, 1]
        
        sensitivity <- true_pos / (true_pos + false_neg)
        specificity <- true_neg / (true_neg + false_pos)
      }
      
      # 计算AUC
      auc_value <- tryCatch({
        roc_obj <- roc(actual_class, prob_predictions, quiet = TRUE)
        auc(roc_obj)
      }, error = function(e) {
        cat("计算AUC错误:", e$message, "\n")
        return(NA)
      })
      
      # 如果是单个变量预测，获取系数和p值
      predictor_coef <- predictor_p_value <- NA
      if(predictor_name != "All Predictors" && !is.null(model)) {
        tryCatch({
          if(predictor_name %in% names(coef(model))) {
            predictor_coef <- coef(model)[predictor_name]
            predictor_p_value <- summary(model)$coefficients[predictor_name, "Pr(>|z|)"]
          }
        }, error = function(e) {
          cat("提取系数错误:", e$message, "\n")
        })
      }
      
      # 打印模型性能
      cat("\n--- ", predictor_name, "模型性能 ---\n")
      cat("AUC:", auc_value, "\n")
      cat("敏感性:", sensitivity, "\n")
      cat("特异性:", specificity, "\n")
      cat("准确率:", accuracy, "\n")
      
      if(!is.na(predictor_coef) && !is.na(predictor_p_value)) {
        cat(predictor_name, "系数:", predictor_coef, "p值:", predictor_p_value, "\n")
      }
      
      # 尝试绘制ROC曲线
      tryCatch({
        if(!is.na(auc_value)) {
          roc_obj <- roc(actual_class, prob_predictions, quiet = TRUE)
          
          # 创建ROC曲线图
          filename_base <- paste0(model_type, "_", gsub("_binary_improvement", "", improvement_var), 
                                  "_day", day_point, "_", gsub(" ", "_", predictor_name))
          
          # 保存PDF格式
          pdf_path <- file.path(output_dir, paste0(filename_base, "_roc.pdf"))
          pdf(pdf_path, width = 8, height = 6)
          plot(roc_obj, main = paste("ROC Curve -", model_type, "-", 
                                     gsub("_binary_improvement", "", improvement_var), 
                                     "- Day", day_point, "-", predictor_name))
          abline(a = 0, b = 1, lty = 2, col = "gray")
          legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), cex = 1.2)
          dev.off()
          
          # 同时保存PNG格式
          png_path <- file.path(output_dir, paste0(filename_base, "_roc.png"))
          png(png_path, width = 480*2, height = 480*1.5, res = 144)
          plot(roc_obj, main = paste("ROC Curve -", model_type, "-", 
                                     gsub("_binary_improvement", "", improvement_var), 
                                     "- Day", day_point, "-", predictor_name))
          abline(a = 0, b = 1, lty = 2, col = "gray")
          legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), cex = 1.2)
          dev.off()
        }
      }, error = function(e) {
        cat("绘制ROC曲线错误:", e$message, "\n")
      })
      
      # 返回结果行
      return(data.frame(
        day_point = day_point,
        target_var = improvement_var,
        predictor = predictor_name,
        auc = auc_value,
        sensitivity = sensitivity,
        specificity = specificity,
        accuracy = accuracy,
        sample_size = nrow(base_data_input),
        predictor_coef = predictor_coef,
        predictor_p_value = predictor_p_value,
        stringsAsFactors = FALSE
      ))
    }
    
    # 1. 所有预测变量模型
    all_vars_formula <- paste(
      improvement_var, "~", 
      wearable_formula_str, " + ",
      "age + gender + bmi + dm_2"
    )
    
    all_model_results <- run_model_and_record(all_vars_formula, "All Predictors", base_data)
    if(!is.null(all_model_results)) {
      performance_df <- rbind(performance_df, all_model_results)
    }
    
    # 2. 对每个可穿戴设备变量构建单独的模型
    for(wearable_var in wearable_vars) {
      var_formula <- paste(
        improvement_var, "~", 
        wearable_var, "+", 
        "age + gender + bmi + dm_2"
      )
      
      var_model_results <- run_model_and_record(var_formula, wearable_var, base_data)
      if(!is.null(var_model_results)) {
        performance_df <- rbind(performance_df, var_model_results)
      }
    }
  }
  
  return(performance_df)
}

#-------------------------------
# 16. 创建存储所有二分类模型性能的数据框
#-------------------------------

# 为血流和厚度数据分别创建存储性能的数据框
all_binary_bloodflow_performance <- data.frame()
all_binary_thickness_performance <- data.frame()

#-------------------------------
# 17. 对每个日期运行二分类预测模型
#-------------------------------

cat("\n开始为每个日期构建OCTA二分类改善预测模型...\n")

for(day_point in daily_timepoints) {
  cat("\n=== 处理日期:", day_point, "===\n")
  
  # 获取当前日期的可穿戴设备数据
  current_day_data <- prepare_daily_data(combined_data, day_point)
  
  # 如果该日期的数据不足，跳过
  if(is.null(current_day_data)) {
    cat("跳过日期", day_point, "，因为数据不足\n")
    next
  }
  
  # 合并血流和厚度数据
  model_data_binary_bloodflow <- current_day_data %>%
    left_join(bloodflow_binary_improvement, by = c("subject_id" = "ID"))
  
  model_data_binary_thickness <- current_day_data %>%
    left_join(thickness_binary_improvement, by = c("subject_id" = "ID"))
  
  # 构建血流二分类改善预测模型
  cat("\n开始构建血流二分类改善预测模型 - 日期:", day_point, "...\n")
  binary_bloodflow_performance <- tryCatch({
    build_simple_binary_model(model_data_binary_bloodflow, significant_binary_bloodflow_pairs, "bloodflow", day_point)
  }, error = function(e) {
    cat("Error in binary bloodflow model:", e$message, "\n")
    return(data.frame())
  })
  
  # 构建厚度二分类改善预测模型
  cat("\n开始构建厚度二分类改善预测模型 - 日期:", day_point, "...\n")
  binary_thickness_performance <- tryCatch({
    build_simple_binary_model(model_data_binary_thickness, significant_binary_thickness_pairs, "thickness", day_point)
  }, error = function(e) {
    cat("Error in binary thickness model:", e$message, "\n")
    return(data.frame())
  })
  
  # 将当前日期的性能结果追加到总结果
  all_binary_bloodflow_performance <- rbind(all_binary_bloodflow_performance, binary_bloodflow_performance)
  all_binary_thickness_performance <- rbind(all_binary_thickness_performance, binary_thickness_performance)
}

#-------------------------------
# 18. 保存二分类模型性能结果
#-------------------------------

# 合并所有二分类性能结果
all_binary_performance <- rbind(
  all_binary_bloodflow_performance,
  all_binary_thickness_performance
)

# 保存到CSV文件
write.csv(all_binary_performance, file.path(output_dir, "octa_binary_improvement_prediction_by_day.csv"), row.names = FALSE)
cat("\n所有二分类模型性能结果已保存到:", file.path(output_dir, "octa_binary_improvement_prediction_by_day.csv"), "\n")

#-------------------------------
# 19. 创建二分类性能变化趋势图表 (折线图)
#-------------------------------

# 按指标和预测变量绘制不同日期的二分类模型性能折线图
create_binary_trend_plots <- function(performance_data, data_type) {
  # 检查数据框是否为空
  if(nrow(performance_data) == 0) {
    cat("警告：没有", data_type, "相关的性能数据可供绘图\n")
    return()
  }
  
  # 获取所有唯一的目标变量
  unique_targets <- unique(performance_data$target_var)
  
  # 对每个目标变量创建单独的图表
  for(target in unique_targets) {
    # 筛选当前目标变量的数据
    target_data <- performance_data %>%
      filter(target_var == target)
    
    # 检查是否有足够的数据点
    if(nrow(target_data) < 3) {
      cat("警告：", target, "没有足够的数据点来创建趋势图\n")
      next
    }
    
    # 准备绘图数据 - 按预测变量分组
    plot_data <- target_data %>%
      select(day_point, predictor, auc) %>%
      # 确保日期按照定义的顺序排序
      mutate(day_point = factor(day_point, levels = daily_timepoints))
    
    # 简化指标名称
    simple_target <- gsub("_binary_improvement$", "", target)
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
    p <- ggplot(plot_data, aes(x = day_point, y = auc, color = predictor, group = predictor)) +
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
        title = paste("Binary Prediction of", ifelse(data_type == "bloodflow", "Blood Flow", "Thickness"), "Improvement:", simple_target),
        x = "Day Relative to Surgery",
        y = "AUC",
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
    plot_filename <- paste0("binary_", data_type, "_", simple_target, "_daily_performance_trend.pdf")
    ggsave(file.path(output_dir, plot_filename), p, width = 10, height = 6)
    
    # 同时保存PNG格式
    plot_filename_png <- paste0("binary_", data_type, "_", simple_target, "_daily_performance_trend.png")
    ggsave(file.path(output_dir, plot_filename_png), p, width = 10, height = 6, dpi = 300)
    
    cat("已保存", data_type, "二分类改善预测趋势图:", simple_target, "\n")
  }
}

# 检查是否有足够的数据来创建图表
if(nrow(all_binary_bloodflow_performance) > 0) {
  # 绘制血流二分类改善预测模型性能趋势图
  create_binary_trend_plots(all_binary_bloodflow_performance, "bloodflow")
} else {
  cat("警告：没有足够的血流数据来创建趋势图\n")
}

if(nrow(all_binary_thickness_performance) > 0) {
  # 绘制厚度二分类改善预测模型性能趋势图
  create_binary_trend_plots(all_binary_thickness_performance, "thickness")
} else {
  cat("警告：没有足够的厚度数据来创建趋势图\n")
}

#-------------------------------
# 20. 创建二分类热图显示不同日期和预测变量的性能
#-------------------------------

create_binary_heatmaps <- function(performance_data, data_type) {
  # 检查数据框是否为空
  if(nrow(performance_data) == 0) {
    cat("警告：没有", data_type, "相关的性能数据可供绘制热图\n")
    return()
  }
  
  # 获取所有唯一的目标变量
  unique_targets <- unique(performance_data$target_var)
  
  # 对每个目标变量创建单独的热图
  for(target in unique_targets) {
    # 筛选当前目标变量的数据
    target_data <- performance_data %>%
      filter(target_var == target)
    
    # 检查是否有足够的数据点
    if(nrow(target_data) < 3) {
      cat("警告：", target, "没有足够的数据点来创建热图\n")
      next
    }
    
    # 简化指标名称
    simple_target <- gsub("_binary_improvement$", "", target)
    simple_target <- gsub("^(PA|VD|SVD|Thickness)_", "", simple_target)
    
    # 创建热图数据 - 排除"All Predictors"
    heatmap_data <- target_data %>%
      filter(predictor != "All Predictors") %>%
      select(day_point, predictor, auc) %>%
      # 确保日期按照定义的顺序排序
      mutate(
        day_point = factor(day_point, levels = daily_timepoints),
        predictor = factor(predictor, levels = c("mean_rhr_1", "mean_bo", "steps_total", "total_sleep"),
                           labels = c("Heart Rate", "Blood Oxygen", "Total Steps", "Sleep Duration"))
      )
    
    # 创建热图
    p_heatmap <- ggplot(heatmap_data, aes(x = day_point, y = predictor, fill = auc)) +
      geom_tile() +
      scale_fill_gradient2(low = "#ffffff", mid = "#bce4d8", high = "#02818a", 
                           midpoint = median(heatmap_data$auc, na.rm = TRUE),
                           na.value = "grey90") +
      labs(
        title = paste("Binary Predictive Power for", ifelse(data_type == "bloodflow", "Blood Flow", "Thickness"), "Improvement:", simple_target),
        x = "Day Relative to Surgery",
        y = "Predictor",
        fill = "AUC"
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
    heatmap_filename <- paste0("binary_", data_type, "_", simple_target, "_daily_heatmap.pdf")
    ggsave(file.path(output_dir, heatmap_filename), p_heatmap, width = 10, height = 6)
    
    # 同时保存PNG格式
    heatmap_filename_png <- paste0("binary_", data_type, "_", simple_target, "_daily_heatmap.png")
    ggsave(file.path(output_dir, heatmap_filename_png), p_heatmap, width = 10, height = 6, dpi = 300)
    
    cat("已保存", data_type, "二分类改善预测热图:", simple_target, "\n")
  }
}

# 检查是否有足够的数据来创建热图
if(nrow(all_binary_bloodflow_performance) > 0) {
  # 绘制血流二分类改善预测模型热图
  create_binary_heatmaps(all_binary_bloodflow_performance, "bloodflow")
} else {
  cat("警告：没有足够的血流数据来创建热图\n")
}

if(nrow(all_binary_thickness_performance) > 0) {
  # 绘制厚度二分类改善预测模型热图
  create_binary_heatmaps(all_binary_thickness_performance, "thickness")
} else {
  cat("警告：没有足够的厚度数据来创建热图\n")
}

# 打印完成信息
cat("\n===== 二分类OCTA指标改善预测分析完成 =====\n")
                       